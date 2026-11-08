# Estimated depth effects across the joint fits
#
# Why: the constant-b joint fits estimate b ~ -0.5 (the fleet favours cells
# with *lower* expected catch), and the sampling model's depth terms flip sign
# between those fits and their b = 0 nulls. Under target = "expected", b
# multiplies the catch model's log expected catch, which already contains the
# catch model's depth smoother. So depth reaches the sampling model by two
# routes, its own poly(log_depth, 2) terms and b x the catch depth effect, and
# the two can trade off. This script plots both routes to see whether the
# negative b is a depth artefact (R/08-check-negative-b.R motivates this).
#
# What: five views of depth for each fit, one per model component:
#   1. Catch model, encounter: the s(log_depth) smoother on the logit scale.
#   2. Catch model, positive catch: the s(log_depth) smoother on the log scale.
#   3. Catch model, expected catch: log encounter probability + log positive
#      catch, fixed effects only (the quantity the "expected" target uses).
#   4. Sampling model, direct: the poly(log_depth, 2) terms alone.
#   5. Sampling model, total: direct + b x catch depth effect, i.e. the whole
#      depth signal in the sampling linear predictor under target = "expected".
# Naive fits have views 1-3 only.
#
# How to read it: if view 3 is the same across b = 0 and constant-b fits while
# view 2 differs, the depth terms are just cancelling b's depth contribution;
# depth then neither drives nor rules out b, and b is identified from the
# catch model's random fields instead.
#
# Each curve is centred on its mean over the frame's cell-years, so only the
# shape is compared (intercepts, year effects, and the time-varying intercept
# shift curves up or down without changing their shape). Point estimates only.
#
# Usage (from the project root):
#   CPUE_SPECIES="pacific ocean perch" Rscript R/08-depth-effects.R
#   CPUE_SPECIES="pacific ocean perch" Rscript R/08-depth-effects.R <fit name> ...
# (named fits only; figure and CSV get a suffix)
# Writes <species>/figs/depth-effects-<cell_km>km.png and
# <species>/output/diagnostics/depth-effects.csv

library(dplyr)
library(ggplot2)
library(sdmTMB)

source("R/species.R")

mesh_cutoff <- 15

# Data -------------------------------------------------------------------------

# The frame (one row per eligible cell x year) sets the plotting range. Each
# joint fit's own frame (stored with the fit) sets its centring and polynomial
# basis, as frames differ between fits.
prepped <- read_prepped() # from R/01-data-prep.R
frame <- prepped$frame
cell_km <- prepped$settings$cell_km

# Fits to compare ----------------------------------------------------------------

# Per species. Redbanded: all 6 km joint fits in models/, each constant-b fit
# with its b = 0 null (same sampling formula and baseline), plus the footprint
# fit, which has complete separation in its footprint terms but a usable depth
# curve. POP: the naive fits and the core-footprint (>= 80%) joint fits.
fit_sets <- list(
  "redbanded rockfish" = tribble(
    ~label, ~fit,
    "zero (year effects)", "zero_xi_expected_tv_dgengamma",
    "zero, rw baseline", "zero_rwbase_xi_expected_tv_dgengamma",
    "constant, rw baseline", "constant_rwbase_xi_expected_tv_dgengamma",
    "zero, rw baseline + region", "zero_rwbase_region_xi_expected_tv_dgengamma",
    "constant, rw baseline + region", "constant_rwbase_region_xi_expected_tv_dgengamma",
    "constant, rw baseline + footprint", "ifootprint_rwbase_xi_expected_tv_dgengamma"
  ),
  "pacific ocean perch" = tribble(
    ~label, ~fit,
    "naive, IID", "naive_tv_iid-delta_gengamma",
    "naive, AR(1)", "naive_tv_ar1-delta_gengamma",
    "rw b, core, IID", "rw_iidbase_fp80frame_xi_expected_tv_dgengamma",
    "rw b, core, AR(1)", "rw_iidbase_fp80frame_xi_expected_tv_ar1_dgengamma",
    "constant b, core, IID", "constant_iidbase_fp80frame_xi_expected_tv_dgengamma",
    "constant b, core, AR(1)", "constant_iidbase_fp80frame_xi_expected_tv_ar1_dgengamma",
    "constant b, core, AR(1), target vessels", "constant_iidbase_fp80frame_target_xi_expected_tv_ar1_dgengamma"
  ),
  # Bocaccio: naive AR(1) fits without target_vessel, by fleet definition
  # (vessels with >= 5 positive trips in >= 1, 3 or 5 years)
  "bocaccio" = tribble(
    ~label, ~fit,
    "naive, AR(1), min 1 yr fleet", "naive_tv_ar1_current-delta_gengamma",
    "naive, AR(1), min 3 yr fleet", "naive_tv_ar1_min3yr-delta_gengamma",
    "naive, AR(1), min 5 yr fleet", "naive_tv_ar1_min5yr-delta_gengamma"
  )
)
fits <- fit_sets[[species]]
if (is.null(fits)) stop("No depth-effect fits listed for ", species)
# Optionally only some of them: Rscript R/08-depth-effects.R <fit name> ...
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  stopifnot(all(args %in% fits$fit))
  fits <- filter(fits, fit %in% args)
}

# Prediction inputs ----------------------------------------------------------------

# Depth grid over the frame's range. The other catch-model covariates are held
# at one value (median location, last year, first month level): with
# re_form = NA below, location drops out, and year and month only shift the
# curve, which centring removes.
depth_grid <- seq(min(frame$log_depth), max(frame$log_depth), length.out = 100)
# target_vessel (catchability, where in the model) at its reference value.
nd <- data.frame(
  X = median(frame$X), Y = median(frame$Y), year = max(frame$year),
  log_depth = depth_grid, month = frame$month[1], target_vessel = 0L
)

# Depth curves per fit -------------------------------------------------------------

# Returns the three curves for one fit, in long format.
depth_curves <- function(fit_name) {
  message("Depth effects: ", fit_name)
  fit <- readRDS(find_fit(fit_name, mesh_cutoff))

  # Centre on the fit's own frame: the sampling data for joint fits, the
  # current frame for naive fits.
  fit_frame <- fit$preferential$spec$data
  if (is.null(fit_frame)) fit_frame <- frame

  # 1-3. Catch model. re_form = NA drops the spatial and spatiotemporal fields
  # and re_form_iid = NA the vessel effects, leaving fixed effects only; the
  # offset is 0 (unit effort). For a delta model, on the link scale est1 is
  # logit(encounter probability), est2 log(positive catch), and est their
  # combination log(expected catch), the quantity the "expected" target uses.
  # Evaluated at the grid to plot and at every frame cell-year to centre on
  # the frame mean.
  catch_at <- function(x) {
    newdata <- nd[rep(1, length(x)), ]
    newdata$log_depth <- x
    predict(fit, newdata = newdata, re_form = NA, re_form_iid = NA,
      offset = rep(0, length(x)))[c("est1", "est2", "est")]
  }
  catch_grid <- catch_at(depth_grid)
  catch_mean <- colMeans(catch_at(fit_frame$log_depth))
  catch <- sweep(as.matrix(catch_grid), 2, catch_mean)
  out <- tibble(
    effect = rep(c("Catch model: encounter (logit)",
      "Catch model: positive catch (log)",
      "Catch model: log expected catch"), each = length(depth_grid)),
    value = c(catch),
    depth = rep(exp(depth_grid), 3)
  )
  b <- NA_real_

  if (!is.null(fit$preferential)) {
    # 4. Sampling model's own depth terms: basis x coefficients, on the logit
    # scale. poly() builds an orthogonal basis from the data it is given, so
    # rebuild it from the fit's frame depths (the data the sampling model was
    # fit to). b is 0 in the zero models (fixed, not estimated).
    depth_poly <- poly(fit_frame$log_depth, 2)
    s <- tidy(fit, model = "sampling", silent = TRUE)
    beta <- s$estimate[match(c("poly(log_depth, 2)1", "poly(log_depth, 2)2"), s$term)]
    b <- s$estimate[s$term == "b_pref"]
    direct <- drop(predict(depth_poly, depth_grid) %*% beta) -
      mean(predict(depth_poly, fit_frame$log_depth) %*% beta)

    # 5. Total depth signal in the sampling model: its own terms plus b times
    # the catch depth effect (both already centred, so their sum is too).
    out <- bind_rows(out, tibble(
      effect = rep(c("Sampling model: direct depth terms",
        "Sampling model: direct + b x catch depth"), each = length(depth_grid)),
      value = c(direct, direct + b * catch[, "est"]),
      depth = rep(exp(depth_grid), 2)
    ))
  }
  mutate(out, fit = fit_name, b = b)
}

# Each fit is ~100 MB on disk, so they are read one at a time.
curves <- fits |>
  rowwise() |>
  mutate(curve = list(depth_curves(fit))) |>
  ungroup() |>
  select(-fit) |>
  tidyr::unnest(curve) |>
  mutate(
    # Put each fit's b in its legend label; keep the table's order
    label = ifelse(is.na(b), label, paste0(label, " (b = ", round(b, 2), ")")),
    label = factor(label, levels = unique(label)),
    effect = factor(effect, levels = unique(effect))
  )

# Plot -------------------------------------------------------------------------

# One panel per view, free y scales (the catch curve spans ~7 log units, the
# sampling curves ~2 logit units). Dashed lines are the b = 0 nulls (redbanded)
# or the naive fits (POP; catch panels only), so a
# dashed/solid pair that differs in the middle panel but matches in the bottom
# one is the cancellation described above. The rug shows where frame cells
# are, so ignore curve ends with little data.
ggplot(curves, aes(depth, value, colour = label,
    linetype = is.na(b) | b == 0)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_line(linewidth = 0.7) +
  geom_rug(aes(x = exp(log_depth)), data = distinct(frame, X, Y, log_depth),
    inherit.aes = FALSE, alpha = 0.05, sides = "b") +
  facet_wrap(~effect, ncol = 1, scales = "free_y") +
  scale_x_log10() + # the models use log depth
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  scale_linetype_manual(values = c("TRUE" = "22", "FALSE" = "solid"),
    guide = "none") +
  labs(x = "Depth (m)", y = "Effect (log or logit scale, centred)",
    title = tools::toTitleCase(species),
    caption = "Dashed: b = 0 or naive models. Rug: frame cell depths.") +
  theme_light() +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(ncol = 2))
suffix <- if (length(args)) paste0("-", paste(args, collapse = "-")) else ""
ggsave(paste0(fig_dir, "/depth-effects-", cell_km, "km", suffix, ".png"),
  width = 8, height = if (all(is.na(curves$b))) 8 else 14)

# Save -------------------------------------------------------------------------

# The curves behind the figure (one row per fit x view x depth grid point)
dir.create(file.path(out_dir, "diagnostics"), showWarnings = FALSE)
write.csv(curves, file.path(out_dir, "diagnostics", paste0("depth-effects", suffix, ".csv")),
  row.names = FALSE)
