# Commercial trawl CPUE index with a joint preferential-sampling model
#
# DRAFT where (cell size, fleet definition, depth trimming, ports, regions) are 
# placeholders to be checked against the data and the assessment's needs.
#
# Approach (docs/preferential-sampling.tex; Ward et al. 2026, docs/):
#   1. Catch model: delta-Gamma spatiotemporal GLMM of catch per tow with
#      log(hours) offset, year effects, depth smoother, month and vessel effects.
#   2. Sampling model: Bernoulli "was cell g fished by the fleet in year t" on a
#      gridded frame of eligible cells, sharing the catch model's latent fields
#      through the preference coefficient b.
#   3. Compare naive vs joint fits, test b = 0, and check time-varying b
#      (the case that biases index *trend*, not just scale, which catchability
#      in the assessment cannot absorb).
#
# Requires sdmTMB from the dev-preferential branch (RTMB backend).

# remotes::install_github("pbs-assess/sdmTMB@dev-preferential")
library(dplyr)
library(ggplot2)
library(sdmTMB)

dir.create("output", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

# Settings ---------------------------------------------------------------------

species <- "redbanded rockfish"
utm_crs <- 32609 # UTM zone 9N
cell_km <- 4 # sampling-frame cell size; b and gamma depend on this choice
min_yrs_cell <- 3 # eligible cell = fished by the fleet in >= this many years
min_positive_trips <- 5 # fleet definition (as in gfplot::tidy_cpue_index())
min_yrs_with_trips <- 5
max_hours <- 10 # drop implausible tow durations
depth_quantiles <- c(0.01, 0.99) # trim depth range using positive tows
majors <- formatC(3:9, width = 2, flag = "0") # 3C, 3D, 5A, 5B, 5C, 5D, 5E (GFFOS major stat area codes)
mesh_cutoff <- 15 # km; main mesh. Revisit after the mesh sensitivity section
mesh_cutoffs_test <- c(8, 10, 12, 15, 20) # km; naive fits compared in mesh sensitivity
run_mesh_test <- FALSE
refit <- FALSE # set TRUE to ignore cached fits

# Approximate trawl offload ports (lon, lat); replace with the ports that
# matter for this fleet, or with a steaming-distance layer.
ports <- data.frame(
  port = c("Prince Rupert", "Port Hardy", "Ucluelet", "Steveston"),
  longitude = c(-130.32, -127.49, -125.55, -123.18),
  latitude = c(54.31, 50.72, 48.94, 49.12)
)

# Data -------------------------------------------------------------------------

# gfdata::get_cpue_index() returns one row per fishing event x species caught
# (all species), so the target species is pulled out and zeros filled here.
raw <- readRDS("data/cpue-index-data.rds")

events <- raw |> mutate(
    fishing_event_id_unique = paste0(
    vessel_registration_number, "-",
    trip_id, "-", fishing_event_id)
  ) |>
  filter(
    toupper(fishery_sector) == "GROUNDFISH TRAWL", # check sector labels
    major_stat_area_code %in% majors
  ) |>
  group_by(fishing_event_id_unique) |> ### need to check all are unique
  summarise(
    fishery_sector = first(fishery_sector),
    trip_id = first(trip_id),
    vessel = first(vessel_registration_number),
    best_date = first(best_date),
    fe_start_date = first(fe_start_date),
    fe_end_date = first(fe_end_date),
    depth = first(best_depth),
    major = first(major_stat_area_code),
    latitude = first(latitude),
    longitude = first(longitude),
    # %in% not ==: some species codes have NA common names, and NA indices
    # would turn the whole event's sum into NA
    catch_kg = sum(catch_kg[toupper(species_common_name) %in% toupper(species)]),
    .groups = "drop"
  ) |>
  mutate(
    year = as.integer(format(best_date, "%Y")),
    month = as.integer(format(best_date, "%m")),
    hours = as.numeric(difftime(fe_end_date, fe_start_date, units = "hours")),
    longitude = -abs(longitude) # GFFOS longitudes are sometimes positive
  ) |>
  filter(
    !is.na(latitude), !is.na(longitude), !is.na(depth), depth > 0,
    !is.na(hours), hours > 0, hours <= max_hours
  )

# Fleet: vessels with >= min_positive_trips positive trips in
# >= min_yrs_with_trips years.
fleet_vessels <- events |>
  filter(catch_kg > 0) |>
  distinct(vessel, year, trip_id) |>
  count(vessel, year, name = "n_pos_trips") |>
  filter(n_pos_trips >= min_positive_trips) |>
  count(vessel, name = "n_years") |>
  filter(n_years >= min_yrs_with_trips) |>
  pull(vessel)


depth_range <- events |>
  filter(catch_kg > 0) |>
  pull(depth) |>
  quantile(depth_quantiles)

d <- events |>
  filter(
    vessel %in% fleet_vessels,
    depth >= depth_range[[1]], depth <= depth_range[[2]]
  ) |>
  add_utm_columns(utm_crs = utm_crs, units = "km")

# Shift coordinates to near zero (UTM northings are ~5500 km, which triggers
# sdmTMB's large-coordinate warning). Offsets are rounded to 100 km so the
# 4 km frame cells and 10 km exploratory bins stay aligned with the UTM grid.
# Add xy_offset back to map X/Y to true UTM km. Anything else placed in this
# space (ports, prediction grids) must be shifted by the same amounts.
xy_offset <- c(X = round(mean(d$X), -2), Y = round(mean(d$Y), -2))
d <- d |>
  mutate(
    X = X - xy_offset[["X"]],
    Y = Y - xy_offset[["Y"]],
    log_depth = log(depth),
    log_hours = log(hours),
    cell_X = floor(X / cell_km) * cell_km + cell_km / 2,
    cell_Y = floor(Y / cell_km) * cell_km + cell_km / 2
  )

# unique(d$fishery_sector)

# Exploratory plots -----------------------------------------------------

# Raw (unstandardized) catch rates by year in space. Binned more coarsely than
# the sampling frame so coastwide facets stay readable.
explore_km <- 10

d_explore <- d |>
  mutate(
    cpue = catch_kg / hours,
    bin_X = floor(X / explore_km) * explore_km + explore_km / 2,
    bin_Y = floor(Y / explore_km) * explore_km + explore_km / 2
  )

bin_year <- d_explore |>
  group_by(year, bin_X, bin_Y) |>
  summarise(
    n_tows = n(),
    prop_pos = mean(catch_kg > 0),
    mean_cpue = mean(cpue),
    .groups = "drop"
  )

# Mean CPUE (kg/h, zeros included) by bin and year; grey = fished, none caught
ggplot(bin_year, aes(bin_X, bin_Y, fill = ifelse(mean_cpue > 0, mean_cpue, NA))) +
  geom_tile(width = explore_km, height = explore_km) +
  facet_wrap(~year, ncol = 6) +
  scale_fill_viridis_c(name = "Mean CPUE\n(kg/h)", trans = "log10",
    na.value = "grey80") +
  coord_fixed() +
  labs(x = "X (km)", y = "Y (km)") +
  theme_light()
ggsave("figs/explore-cpue-map.png", width = 12, height = 14)

# Encounter rate by bin and year
ggplot(bin_year, aes(bin_X, bin_Y, fill = prop_pos)) +
  geom_tile(width = explore_km, height = explore_km) +
  facet_wrap(~year, ncol = 6) +
  scale_fill_viridis_c(name = "Proportion\npositive", limits = c(0, 1)) +
  coord_fixed() +
  labs(x = "X (km)", y = "Y (km)") +
  theme_light()
ggsave("figs/explore-encounter-map.png", width = 12, height = 14)

# Spread of bin-level CPUE within each year (positive bins only)
bin_year |>
  filter(mean_cpue > 0) |>
  ggplot(aes(factor(year), mean_cpue)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_y_log10() +
  labs(x = "Year", y = "Mean CPUE per bin (kg/h, positive bins)") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("figs/explore-cpue-spread.png", width = 10, height = 5)

# Annual summaries: effort, footprint, encounter rate, and raw catch rates
d_explore |>
  group_by(year) |>
  summarise(
    `Tows` = n(),
    `Vessels` = n_distinct(vessel),
    `Frame-size cells fished` = n_distinct(paste(cell_X, cell_Y)),
    `Proportion positive` = mean(catch_kg > 0),
    `Mean CPUE (kg/h)` = mean(cpue),
    `Mean positive CPUE (kg/h)` = mean(cpue[catch_kg > 0]),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(-year) |>
  ggplot(aes(year, value)) +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(~name, scales = "free_y", ncol = 2) +
  labs(x = "Year", y = NULL) +
  theme_light()
ggsave("figs/explore-annual-summary.png", width = 9, height = 8)

# Centre of gravity of effort vs catch: if catch moves away from effort (or
# effort tracks catch), that hints at shifting preferential sampling.
d_explore |>
  group_by(year) |>
  summarise(
    effort_X = mean(X), effort_Y = mean(Y),
    catch_X = weighted.mean(X, catch_kg, na.rm = TRUE), catch_Y = weighted.mean(Y, catch_kg, na.rm = TRUE),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(-year, names_to = c("what", "axis"), names_sep = "_") |>
  ggplot(aes(year, value, colour = what)) +
  geom_line() +
  facet_wrap(~axis, scales = "free_y", ncol = 1) +
  labs(x = "Year", y = "Centre of gravity (km)", colour = NULL) +
  theme_light()
ggsave("figs/explore-centre-of-gravity.png", width = 7, height = 6)


# Sampling frame ---------------------------------------------------------------

# One row per eligible cell x year. "Eligible" here = within the fleet's
# footprint (fished in >= min_yrs_cell years). Cells that could never be fished
# (land, permanent closures such as sponge-reef or RCA closures, out of depth
# range) must be omitted, not coded 0. ??Time-varying closures should be dropped
# from the years they were closed.??? I think this is wrong. 

region_lu <- c(
  "03" = "WCVI", "04" = "WCVI", "05" = "QCS", "06" = "QCS",
  "07" = "HS", "08" = "HS", "09" = "WCHG"
)

cells <- d |>
  group_by(cell_X, cell_Y) |>
  summarise(
    n_years = n_distinct(year),
    # Cell depth from observed tows; swap in a bathymetry layer if available.
    log_depth = mean(log_depth),
    major = names(which.max(table(major))),
    .groups = "drop"
  ) |>
  filter(n_years >= min_yrs_cell) |>
  mutate(region = factor(region_lu[as.character(major)]))

ports <- add_utm_columns(ports, utm_crs = utm_crs, units = "km") |>
  mutate(X = X - xy_offset[["X"]], Y = Y - xy_offset[["Y"]])
port_dist <- sqrt(
  outer(cells$cell_X, ports$X, "-")^2 + outer(cells$cell_Y, ports$Y, "-")^2
)
# Straight-line distance ignores land; a proxy for access cost.
cells$log_dist_port <- log(pmax(apply(port_dist, 1, min), cell_km))

# Keep catch data inside the frame footprint so both likelihoods describe the
# same fishery.
d <- semi_join(d, cells, by = c("cell_X", "cell_Y"))

years <- sort(unique(d$year))
stopifnot(all(diff(years) == 1)) # else use extra_time for coefficient = "rw"

ref_month <- as.integer(names(which.max(table(d$month))))
d$month <- factor(d$month)
d$vessel <- factor(d$vessel)

sampled_cy <- d |>
  distinct(cell_X, cell_Y, year) |>
  mutate(sampled = 1L)

frame <- cells |>
  select(cell_X, cell_Y, log_depth, log_dist_port, region) |>
  cross_join(data.frame(year = years)) |>
  left_join(sampled_cy, by = c("cell_X", "cell_Y", "year")) |>
  mutate(
    sampled = coalesce(sampled, 0L),
    # Catchability covariates held at reference values on the frame:
    month = factor(ref_month, levels = levels(d$month))
  ) |>
  rename(X = cell_X, Y = cell_Y) |>
  as.data.frame()

frame |>
  group_by(year) |>
  summarise(prop_sampled = mean(sampled), n_cells = n()) |>
  print(n = Inf)

# Index grid: frame cells, reference month, unit effort (offset 0). Used for
# the mesh sensitivity and the final indices.
grid <- frame |>
  select(X, Y, year, log_depth, month)

# Mesh -------------------------------------------------------------------------

# Every frame location must fall inside the mesh, so triangulate over both the
# tows and the frame cell centres, then attach the mesh to the catch data.
mesh_pts <- bind_rows(
  select(as.data.frame(d), X, Y),
  distinct(frame, X, Y)
)
build_mesh <- function(cutoff) {
  mesh0 <- make_mesh(mesh_pts, c("X", "Y"), cutoff = cutoff)
  make_mesh(d, c("X", "Y"), mesh = mesh0$mesh)
}


# Models -----------------------------------------------------------------------

# Year effects are fixed (needed for the preference-deviation centering).
# Depth is treated as biological (kept at cell values on the frame); month and
# vessel as catchability (month at reference level; vessel IID effects excluded
# from the preference surface by re_form_iid = NA, the default).
# Note: no scale()/mean()-type calls in the formula; they would be recomputed
# on the frame.
catch_formula <- catch_kg ~ 0 + factor(year) + s(log_depth, k = 5) +
  month + (1 | vessel)

# Free year intercepts absorb changes in total effort; distance to port and
# depth absorb access-driven (density-independent) sampling structure.
samp_formula <- sampled ~ 0 + factor(year) + log_dist_port + poly(log_depth, 2)

# Sensitivity: allow the regional distribution of effort to shift by year
# (unmodelled shifts get attributed to b; see docs). Watch for separation
# errors if a region-year has all or no cells fished.
samp_formula_region <- update(samp_formula, . ~ . + factor(year):region)

fit_cpue <- function(preferential = NULL, fit_mesh = mesh) {
  sdmTMB(
    catch_formula,
    data = d,
    mesh = fit_mesh,
    time = "year",
    family = delta_gamma(),
    offset = d$log_hours,
    spatial = "on",
    spatiotemporal = "iid",
    preferential = preferential,
    # Optional PC priors to stabilize range/SD estimates:
    # priors = sdmTMBpriors(matern_s = pc_matern(range_gt = 10, sigma_lt = 3),
    #                       matern_st = pc_matern(range_gt = 10, sigma_lt = 3)),
    control = sdmTMBcontrol(backend = "rtmb", multiphase = TRUE),
    silent = FALSE
  )
}

fit_file <- function(name, cutoff) {
  file.path("output", paste0("fit-", name, "-cutoff", cutoff, ".rds"))
}

# Mesh sensitivity -------------------------------------------------------------

# Fit the naive model at several cutoffs and compare cost, estimated ranges,
# and the index. Rule of thumb: the cutoff should be well below the estimated
# range (range/cutoff >= ~5). If the index trend is stable across cutoffs, use
# the coarsest adequate mesh for the (much slower) joint fits.
if (run_mesh_test) {
  mesh_test <- lapply(mesh_cutoffs_test, function(cutoff) {
    m <- build_mesh(cutoff)
    f <- fit_file("naive", cutoff)
    secs <- NA_real_
    if (file.exists(f) && !refit) {
      fit <- readRDS(f)
    } else {
      message("Fitting naive model, cutoff = ", cutoff, " km (",
        m$mesh$n, " vertices)")
      secs <- system.time(fit <- fit_cpue(fit_mesh = m))[["elapsed"]]
      saveRDS(fit, f)
    }
    ranges <- bind_rows(
      tidy(fit, "ran_pars", model = 1) |> mutate(component = "encounter"),
      tidy(fit, "ran_pars", model = 2) |> mutate(component = "positive")
    ) |>
      filter(term == "range") |>
      mutate(cutoff = cutoff)
    idx <- get_index(fit, newdata = grid, area = cell_km^2,
      bias_correct = TRUE, predict_args = list(re_form_iid = NA)) |>
      mutate(cutoff = cutoff)
    list(
      summary = data.frame(
        cutoff = cutoff,
        n_vertices = m$mesh$n,
        fit_seconds = secs, # NA when loaded from cache
        converged = isTRUE(sanity(fit, silent = TRUE)$all_ok),
        min_range = min(ranges$estimate),
        min_range_over_cutoff = min(ranges$estimate) / cutoff,
        slope = coef(lm(log_est ~ year, data = idx))[[2]]
      ),
      ranges = ranges,
      index = idx
    )
  })

  mesh_summary <- bind_rows(lapply(mesh_test, `[[`, "summary"))
  print(mesh_summary)
  bind_rows(lapply(mesh_test, `[[`, "ranges")) |>
    select(cutoff, component, estimate, conf.low, conf.high) |>
    print()

  bind_rows(lapply(mesh_test, `[[`, "index")) |>
    group_by(cutoff) |>
    mutate(rel_est = est / exp(mean(log_est))) |>
    ungroup() |>
    ggplot(aes(year, rel_est, colour = factor(cutoff))) +
    geom_line() +
    labs(x = "Year", y = "Relative index (geometric mean = 1)",
      colour = "Mesh cutoff (km)") +
    theme_light()
  ggsave("figs/mesh-sensitivity-index.png", width = 8, height = 5)
  write.csv(mesh_summary, "output/mesh-sensitivity.csv", row.names = FALSE)
}

# Main mesh for the model comparison below
mesh <- build_mesh(mesh_cutoff)
mesh$mesh$n
plot(mesh)

# Joint and naive fits ---------------------------------------------------------

pref <- function(coefficient = "constant", spatial = "off",
                 formula = samp_formula) {
  preferential_sampling(
    formula,
    data = frame,
    target = "fields",
    coefficient = coefficient,
    spatial = spatial
  )
}

model_specs <- list(
  naive = NULL, # catch only
  zero = pref("zero"), # joint, b = 0 (null for the LRT)
  constant = pref("constant"), # constant preference
  rw = pref("rw"), # time-varying preference (random walk)
  rw_xi = pref("rw", spatial = "on"), # + sampling-only spatial field
  rw_region = pref("rw", formula = samp_formula_region)
)

fits <- lapply(names(model_specs), function(nm) {
  f <- fit_file(nm, mesh_cutoff) # naive reuses the mesh-sensitivity fit
  if (file.exists(f) && !refit) return(readRDS(f))
  message("Fitting: ", nm)
  fit <- tryCatch(fit_cpue(model_specs[[nm]]), error = function(e) {
    message("  failed: ", conditionMessage(e))
    NULL
  })
  if (!is.null(fit)) saveRDS(fit, f)
  fit
})
names(fits) <- names(model_specs)
fits <- Filter(Negate(is.null), fits)

# Diagnostics ------------------------------------------------------------------

lapply(fits, sanity)
fits$rw
tidy(fits$rw, model = "sampling", conf.int = TRUE) # gamma and b_pref
tidy(fits$rw, "ran_pars", model = "sampling") # sigma_b_pref, etc.

# Test for preferential sampling: zero vs constant (same frame and sampling
# formula), LRT on 1 df. Joint fits are not comparable with the naive fit,
# whose likelihood lacks the sampling data.
lrt <- 2 * (as.numeric(logLik(fits$constant)) - as.numeric(logLik(fits$zero)))
data.frame(
  lrt = lrt,
  p = pchisq(lrt, df = 1, lower.tail = FALSE)
)
AIC(fits$zero, fits$constant, fits$rw)

# Preference coefficient by year
b_t <- tidy(fits$rw, "ran_vals", model = "sampling") |>
  filter(grepl("b_pref", term))
ggplot(b_t, aes(seq_along(estimate) + min(years) - 1, estimate)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Year", y = "Preference coefficient b_t")
ggsave("figs/b-pref-by-year.png", width = 6, height = 4)

# Fitted vs observed sampling for a few years
ps <- predict_sampling(fits$rw, type = "response")
ps |>
  filter(year %in% years[round(seq(1, length(years), length.out = 4))]) |>
  ggplot(aes(X, Y, fill = est)) +
  geom_tile(width = cell_km, height = cell_km) +
  geom_point(data = \(x) filter(x, sampled == 1), size = 0.2, colour = "white") +
  facet_wrap(~year) +
  scale_fill_viridis_c(name = "P(fished)") +
  coord_fixed()
ggsave("figs/sampling-probability.png", width = 8, height = 8)

# Catch-model residuals (conditional on fitted fields)
d$resid1 <- residuals(fits$rw, type = "mle-mvn", model = 1)
d$resid2 <- residuals(fits$rw, type = "mle-mvn", model = 2)
qqnorm(d$resid2[d$catch_kg > 0])
qqline(d$resid2[d$catch_kg > 0])

# Index ------------------------------------------------------------------------

# Standardized index over the fishery footprint (`grid`, built with the frame):
# reference month, unit effort, vessel effects excluded. Consider also an index
# over the full stock area (Ward et al. used both).
indices <- lapply(fits, function(fit) {
  get_index(
    fit,
    newdata = grid,
    area = cell_km^2,
    bias_correct = TRUE, # keep TRUE for joint delta models (see docs)
    predict_args = list(re_form_iid = NA)
  )
}) |>
  bind_rows(.id = "model") |>
  group_by(model) |>
  mutate(
    geo_mean = exp(mean(log_est)),
    rel_est = est / geo_mean,
    rel_lwr = lwr / geo_mean,
    rel_upr = upr / geo_mean
  ) |>
  ungroup()

ggplot(indices, aes(year, rel_est, colour = model, fill = model)) +
  geom_ribbon(aes(ymin = rel_lwr, ymax = rel_upr), alpha = 0.15, colour = NA) +
  geom_line() +
  labs(x = "Year", y = "Relative index (geometric mean = 1)") +
  theme_light()
ggsave("figs/index-comparison.png", width = 8, height = 5)

# Log-linear trend by model: the joint-vs-naive difference here is the
# quantity that matters for the assessment (scale is absorbed by q).
indices |>
  group_by(model) |>
  summarise(slope = coef(lm(log_est ~ year))[[2]])

saveRDS(indices, "output/indices.rds")

# Assessment input -------------------------------------------------------------

# Base-case choice is a judgement call: e.g., "rw" if sigma_b is supported and
# the LRT rejects b = 0, else "constant"/"naive"; report the others as
# sensitivities. SE is the log-scale SE from get_index().
base_model <- "rw"
ss_index <- indices |>
  filter(model == base_model) |>
  transmute(
    year,
    month = 7,
    fleet = NA_integer_, # set to the CPUE fleet number in the SS3 model
    obs = rel_est,
    se_log = se
  )
write.csv(ss_index, "output/cpue-index-ss3.csv", row.names = FALSE)
write.csv(indices, "output/cpue-index-all-models.csv", row.names = FALSE)
