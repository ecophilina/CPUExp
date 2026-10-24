# Estimated depth effects across the joint fits
#
# Why: the constant-b joint fits estimate b ~ -0.5 (the fleet favours cells
# with *lower* expected catch), and the sampling model's depth terms flip sign
# between those fits and their b = 0 nulls. Under target = "expected", b
# multiplies the catch model's log expected catch, which already contains the
# catch model's depth smoother. So depth reaches the sampling model by two
# routes, its own poly(log_depth, 2) terms and b x the catch depth effect, and
# the two can trade off. This script plots both routes to see whether the
# negative b is a depth artefact (R/check-negative-b.R motivates this).
#
# What: three views of depth for each fit:
#   1. Catch model: the s(log_depth) smoother, as log expected catch
#      (log encounter probability + log positive catch), fixed effects only.
#   2. Sampling model, direct: the poly(log_depth, 2) terms alone.
#   3. Sampling model, total: direct + b x catch depth effect, i.e. the whole
#      depth signal in the sampling linear predictor under target = "expected".
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
# Usage (from the project root): Rscript R/depth-effects.R
# Writes figs/depth-effects-<cell_km>km.png and
# output/diagnostics/depth-effects.csv

library(dplyr)
library(ggplot2)
library(sdmTMB)

mesh_cutoff <- 15

# Data -------------------------------------------------------------------------

# The frame is the sampling model's data (one row per eligible cell x year);
# its depths set the plotting range, the centring, and the polynomial basis.
prepped <- readRDS("data/cpue-index-clean.rds") # from R/data-prep.R
frame <- prepped$frame
cell_km <- prepped$settings$cell_km

# Fits to compare ----------------------------------------------------------------

# All 6 km joint fits in models/: each constant-b fit with its b = 0 null
# (same sampling formula and baseline), plus the footprint fit, which has
# complete separation in its footprint terms but a usable depth curve.
fits <- tribble(
  ~label, ~fit,
  "zero (year effects)", "zero_xi_expected_tv_dgengamma",
  "zero, rw baseline", "zero_rwbase_xi_expected_tv_dgengamma",
  "constant, rw baseline", "constant_rwbase_xi_expected_tv_dgengamma",
  "zero, rw baseline + region", "zero_rwbase_region_xi_expected_tv_dgengamma",
  "constant, rw baseline + region", "constant_rwbase_region_xi_expected_tv_dgengamma",
  "constant, rw baseline + footprint", "ifootprint_rwbase_xi_expected_tv_dgengamma"
)

# Prediction inputs ----------------------------------------------------------------

# Depth grid over the frame's range. The other catch-model covariates are held
# at one value (median location, last year, first month level): with
# re_form = NA below, location drops out, and year and month only shift the
# curve, which centring removes.
depth_grid <- seq(min(frame$log_depth), max(frame$log_depth), length.out = 100)
nd <- data.frame(
  X = median(frame$X), Y = median(frame$Y), year = max(frame$year),
  log_depth = depth_grid, month = frame$month[1]
)
# poly() builds an orthogonal basis from the data it is given, so its
# coefficients only mean something on that basis. Rebuild the basis from the
# frame's depths (the data the sampling model was fit to) and evaluate it at
# the grid (for plotting) and at the frame (for centring).
depth_poly <- poly(frame$log_depth, 2)
poly_grid <- predict(depth_poly, depth_grid)
poly_frame <- predict(depth_poly, frame$log_depth)

# Depth curves per fit -------------------------------------------------------------

# Returns the three curves for one fit, in long format.
depth_curves <- function(fit_name) {
  message("Depth effects: ", fit_name)
  fit <- readRDS(file.path("models",
    paste0("fit-", fit_name, "-cutoff", mesh_cutoff, ".rds")))

  # 1. Catch model. re_form = NA drops the spatial and spatiotemporal fields
  # and re_form_iid = NA the vessel effects, leaving fixed effects only; the
  # offset is 0 (unit effort). For a delta model, `est` on the link scale is
  # log(encounter probability) + log(positive catch), i.e. log expected catch,
  # the quantity the "expected" target uses. Evaluated at the grid to plot and
  # at every frame cell-year to centre on the frame mean.
  catch_at <- function(x) {
    newdata <- nd[rep(1, length(x)), ]
    newdata$log_depth <- x
    predict(fit, newdata = newdata, re_form = NA, re_form_iid = NA,
      offset = rep(0, length(x)))$est
  }
  catch_grid <- catch_at(depth_grid)
  catch_mean <- mean(catch_at(frame$log_depth))

  # 2. Sampling model's own depth terms: basis x coefficients, on the logit
  # scale. b is 0 in the zero models (fixed, not estimated).
  s <- tidy(fit, model = "sampling", silent = TRUE)
  beta <- s$estimate[match(c("poly(log_depth, 2)1", "poly(log_depth, 2)2"), s$term)]
  b <- s$estimate[s$term == "b_pref"]
  direct_grid <- drop(poly_grid %*% beta)
  direct_mean <- mean(poly_frame %*% beta)

  # 3. Total depth signal in the sampling model: its own terms plus b times
  # the catch depth effect (both already centred, so their sum is too).
  catch <- catch_grid - catch_mean
  direct <- direct_grid - direct_mean
  bind_rows(
    tibble(effect = "Catch model: log expected catch", value = catch),
    tibble(effect = "Sampling model: direct depth terms", value = direct),
    tibble(effect = "Sampling model: direct + b x catch depth", value = direct + b * catch)
  ) |>
    mutate(depth = rep(exp(depth_grid), 3), fit = fit_name, b = b)
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
    label = factor(paste0(label, " (b = ", round(b, 2), ")"),
      levels = unique(paste0(label, " (b = ", round(b, 2), ")"))),
    effect = factor(effect, levels = unique(effect))
  )

# Plot -------------------------------------------------------------------------

# One panel per view, free y scales (the catch curve spans ~7 log units, the
# sampling curves ~2 logit units). Dashed lines are the b = 0 nulls, so a
# dashed/solid pair that differs in the middle panel but matches in the bottom
# one is the cancellation described above. The rug shows where frame cells
# are, so ignore curve ends with little data.
ggplot(curves, aes(depth, value, colour = label, linetype = b == 0)) +
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
    caption = "Dashed: b = 0 models. Rug: frame cell depths.") +
  theme_light() +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(ncol = 2))
ggsave(paste0("figs/depth-effects-", cell_km, "km.png"), width = 8, height = 10)

# Save -------------------------------------------------------------------------

# The curves behind the figure (one row per fit x view x depth grid point)
write.csv(curves, file.path("output", "diagnostics", "depth-effects.csv"),
  row.names = FALSE)
