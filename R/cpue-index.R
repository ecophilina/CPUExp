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
# Data prep, exploratory plots, the sampling frame, and the index grid are in
# R/data-prep.R, which saves data/cpue-index-clean.rds. Run that first.
#
# Requires sdmTMB from the dev-preferential branch (RTMB backend).

# remotes::install_github("pbs-assess/sdmTMB@dev-preferential")
library(dplyr)
library(ggplot2)
library(sdmTMB)

dir.create("output", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

# Settings ---------------------------------------------------------------------

mesh_cutoff <- 15 # km; main mesh. Revisit after the mesh sensitivity section
mesh_cutoffs_test <- c(8, 10, 12, 15, 20) # km; naive fits compared in mesh sensitivity
run_mesh_test <- FALSE
refit <- FALSE # set TRUE to ignore cached fits

# Data -------------------------------------------------------------------------

prepped <- readRDS("data/cpue-index-clean.rds") # from R/data-prep.R
d <- prepped$d # catch data: one row per tow
frame <- prepped$frame # sampling frame: one row per eligible cell x year
grid <- prepped$grid # index grid: frame cells, reference month, unit effort
years <- prepped$years
cell_km <- prepped$settings$cell_km

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


# multiphase = FALSE skips the first optimization phase (fields off). Note: it
# does not lower peak memory. Joint models with time-varying b (iid/rw) need
# ~16-18 GB to build the full model either way (vs ~9 GB for constant b).
fit_cpue <- function(preferential = NULL, fit_mesh = mesh,
                     family = delta_gamma(), multiphase = TRUE) {
  sdmTMB(
    catch_formula,
    data = d,
    mesh = fit_mesh,
    time = "year",
    family = family,
    offset = d$log_hours,
    spatial = "on",
    spatiotemporal = "iid",
    preferential = preferential,
    # Optional PC priors to stabilize range/SD estimates:
    # priors = sdmTMBpriors(matern_s = pc_matern(range_gt = 10, sigma_lt = 3),
    #                       matern_st = pc_matern(range_gt = 10, sigma_lt = 3)),
    control = sdmTMBcontrol(backend = "rtmb", multiphase = multiphase),
    silent = FALSE
  )
}


# Time-varying (random-walk) intercept in place of factor(year). Uses "rw0"
# with an intercept in the formula: with "rw" the first value is free, and
# `0 + ... + month` gives month all 12 levels (an implicit intercept), which
# made the model unidentifiable (non-positive-definite Hessians).
catch_formula_tv <- catch_kg ~ 1 + s(log_depth, k = 5) +
  month + (1 | vessel)

# The gamma prior on the random-walk SD (shape 1 / cv^2 = 4) has zero density
# at sigma_V = 0, so it keeps the time-varying SD from collapsing.
fit_cpue_tv <- function(preferential = NULL, fit_mesh = mesh,
                     family = delta_gamma(),
                     priors = sdmTMBpriors(sigma_V = gamma_cv(0.1, 0.5))) {
  sdmTMB(
    catch_formula_tv,
    data = d,
    mesh = fit_mesh,
    time = "year",
    time_varying = ~1,
    time_varying_type = "rw0",
    priors = priors,
    family = family,
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

# Free year intercepts absorb changes in total effort; distance to port and
# depth absorb access-driven (density-independent) sampling structure.
samp_formula <- sampled ~ 0 + factor(year) + log_dist_port + poly(log_depth, 2)

samp_formula2 <- sampled ~ 1 + log_dist_port + poly(log_depth, 2)

# Sensitivity: allow the regional distribution of effort to shift by year
# (unmodelled shifts get attributed to b; see docs). Watch for separation
# errors if a region-year has all or no cells fished.
samp_formula_region <- update(samp_formula, . ~ . + factor(year):region)


pref <- function(coefficient = "constant", spatial = "off",
                 target = "expected", # changed to expected after Oct 20
                 baseline = "off", # off, iid, or rw
                 formula = samp_formula) {
  # A baseline ("iid"/"rw") needs a formula with an intercept and no
  # factor(year) (e.g., samp_formula2); "off" pairs with samp_formula.
  preferential_sampling(
    formula,
    data = frame,
    re_form_iid = NA, # IID REs from the surface (default)
    offset = 0, #surface offset, 0 = unit exposure
    target = target, # fields or expected
    coefficient = coefficient, # constant, iid , rw, zero
    baseline = baseline,
    spatial = spatial
  )
}

# Joint models (target = "expected", pref()'s default) with the sampling field.
# Each zero_* model is the null for the matching constant_* model (same family,
# frame, formula, and field) for a likelihood-ratio test of b. Fit these in
# parallel with R/run-joint-fits.sh (via R/fit-joint.R); this script then loads
# the cached fits.
model_specs <- list(
  constant_xi_expected_dgamma = pref("constant", spatial = "on"),
  zero_xi_expected_dgamma = pref("zero", spatial = "on"),
  constant_xi_expected_dlnorm = pref("constant", spatial = "on"),
  zero_xi_expected_dlnorm = pref("zero", spatial = "on"),
  constant_rwbase_xi_expected_dlnorm = pref("constant", spatial = "on",
    baseline = "rw", formula = samp_formula2),
  # time-varying preference (random walk in b_t), rw baseline. With the
  # expected target, three attempts were killed out of memory (peak ~16 GB
  # resident, vs ~9 GB for constant b) when the full model was built after the
  # first optimization phase; the fields-target version is used instead.
  rw_rwbase_xi_expected_dlnorm = pref("rw", spatial = "on",
    baseline = "rw", formula = samp_formula2),
  rw_rwbase_xi_fields_dlnorm = pref("rw", spatial = "on", target = "fields",
    baseline = "rw", formula = samp_formula2),
  # time-varying preference (IID b_t deviations), rw baseline
  iid_rwbase_xi_expected_dlnorm = pref("iid", spatial = "on",
    baseline = "rw", formula = samp_formula2),
  # constant preference, generalized gamma
  constant_xi_expected_dgengamma = pref("constant", spatial = "on"),
  # same with the time-varying (rw0) catch intercept, and its b = 0 null.
  # Needs an sdmTMB dev-preferential build newer than 1.2.0.9001, which errors
  # "Unsupported: `time_varying`".
  constant_xi_expected_tv_dgengamma = pref("constant", spatial = "on"),
  zero_xi_expected_tv_dgengamma = pref("zero", spatial = "on")
  # # constant preference, fields target
  # constant_xi = pref("constant", spatial = "on", target = "fields"),
  # # constant preference, rw baseline, delta-Gamma
  # constant_rwbase_xi = pref("constant", spatial = "on", baseline = "rw",
  #   formula = samp_formula2)
  )

# Models whose catch model is fit_cpue_tv() (time-varying intercept) rather
# than fit_cpue() (factor(year) effects)
model_tv <- c("constant_xi_expected_tv_dgengamma", "zero_xi_expected_tv_dgengamma")

# Observation family per model (delta_gamma() if not listed)
model_families <- list(
  constant_xi_expected_dgamma = delta_gamma(),
  zero_xi_expected_dgamma = delta_gamma(),
  constant_xi_expected_dlnorm = delta_lognormal(),
  zero_xi_expected_dlnorm = delta_lognormal(),
  constant_rwbase_xi_expected_dlnorm = delta_lognormal(),
  rw_rwbase_xi_expected_dlnorm = delta_lognormal(),
  rw_rwbase_xi_fields_dlnorm = delta_lognormal(),
  iid_rwbase_xi_expected_dlnorm = delta_lognormal(),
  constant_xi_expected_dgengamma = delta_gengamma(),
  constant_xi_expected_tv_dgengamma = delta_gengamma(),
  zero_xi_expected_tv_dgengamma = delta_gengamma()
)


# model_specs <- list(
#   naive = NULL, # catch only
#   zero = pref("zero"), # joint, b = 0 (null for the LRT)
#   constant_xi = pref("constant", spatial = "on"), # constant preference
#   rw_xi = pref("rw", spatial = "on"), # + sampling-only spatial field
#   constant_expected = pref("constant", target = "expected"), # time-varying preference (random walk)
#   rw_region = pref("rw", formula = samp_formula_region)
# )

fits <- lapply(names(model_specs), function(nm) {
  f <- fit_file(nm, mesh_cutoff) # naive reuses the mesh-sensitivity fit
  if (file.exists(f) && !refit) return(readRDS(f))
  message("Fitting: ", nm)
  family <- model_families[[nm]]
  if (is.null(family)) family <- delta_gamma()
  fit_fun <- if (nm %in% model_tv) fit_cpue_tv else fit_cpue
  fit <- tryCatch(fit_fun(model_specs[[nm]], family = family), error = function(e) {
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
fits$constant_xi
tidy(fits$constant_xi, model = "sampling", conf.int = TRUE) # gamma and b_pref
tidy(fits$constant_xi, "ran_pars", model = "sampling") # sigma_b_pref, etc.

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
