# Is the negative preference coefficient (b) real or a depth artefact?
#
# constant_rwbase_xi_expected_tv_dgengamma estimates b ~ -0.5: cells with
# higher expected catch are fished *less* than depth and distance to port
# predict. Between the b = 0 null and that fit, the sampling model's depth
# terms flip sign, so b may be trading off against depth (the expected-catch
# target already contains the catch model's depth smoother). Checks:
#   1. How strongly expected catch on the frame tracks depth.
#   2. Sampling rate vs expected catch overall and *within depth bands*. If the
#      decline only appears across bands, b is absorbing depth (Simpson's
#      paradox); if it holds within bands, avoidance is in the data.
#   3. Outside-model logistic regressions of sampled on log expected catch,
#      adding year, depth, and port distance terms in turn.
#
# Expected catch comes from the b = 0 fit, so the surface is not shaped by the
# preference term itself; the constant-b fit's surface is used as a check.
# Predictions are plug-in (no bias correction).
#
# Usage (from the project root): Rscript R/check-negative-b.R

library(dplyr)
library(ggplot2)
library(sdmTMB)

source("R/species.R")

mesh_cutoff <- 15
n_depth_bands <- 5
n_mu_bins <- 10

prepped <- read_prepped() # from R/data-prep.R
frame <- prepped$frame
cell_km <- prepped$settings$cell_km

pred_dir <- file.path(out_dir, "diagnostics")
dir.create(pred_dir, showWarnings = FALSE)

fit_names <- c(
  zero = "zero_rwbase_xi_expected_tv_dgengamma",
  constant = "constant_rwbase_xi_expected_tv_dgengamma"
)

# Log expected catch per unit effort on the frame (reference month, vessel
# effects excluded), cached per fit
# fixed = TRUE drops the random fields (fixed effects and time-varying
# intercept only)
predict_log_mu <- function(fit_name, fixed = FALSE) {
  tag <- if (fixed) "logmu-fixed-" else "logmu-"
  f <- file.path(pred_dir, paste0(tag, fit_name, "-cutoff", mesh_cutoff, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  message("Predicting: ", fit_name, if (fixed) " (fixed effects only)")
  fit <- readRDS(find_fit(fit_name, mesh_cutoff))
  p <- predict(fit, newdata = select(frame, X, Y, year, log_depth, month),
    type = "response", re_form = if (fixed) NA else NULL, re_form_iid = NA,
    offset = rep(0, nrow(frame)))
  # Delta models: combined expected catch (encounter prob x positive mean)
  est <- if ("est" %in% names(p)) p$est else p$est1 * p$est2
  out <- log(est)
  saveRDS(out, f)
  out
}

fr <- frame |>
  mutate(
    log_mu = predict_log_mu(fit_names[["zero"]]),
    log_mu_constant = predict_log_mu(fit_names[["constant"]]),
    depth_band = cut(log_depth,
      quantile(log_depth, seq(0, 1, length.out = n_depth_bands + 1)),
      include.lowest = TRUE,
      labels = paste0("Q", seq_len(n_depth_bands)))
  ) |>
  group_by(year) |>
  mutate(log_mu_dev = log_mu - mean(log_mu)) |> # within-year deviation
  ungroup()

# 1. Expected catch vs depth -----------------------------------------------------

cat("\nCorrelation of log expected catch with log depth (cell-years):",
  round(cor(fr$log_mu, fr$log_depth), 2), "\n")
cat("Share of log expected catch variance explained by poly(log_depth, 2):",
  round(summary(lm(log_mu ~ poly(log_depth, 2), fr))$r.squared, 2), "\n")
cat("Correlation of zero vs constant fit surfaces:",
  round(cor(fr$log_mu, fr$log_mu_constant), 3), "\n")

ggplot(fr, aes(exp(log_depth), log_mu)) +
  geom_point(aes(colour = factor(sampled)), alpha = 0.15, size = 0.4) +
  geom_smooth(colour = "black", se = FALSE) +
  scale_x_log10() +
  scale_colour_manual(values = c("0" = "grey60", "1" = "#2c7fb8"),
    labels = c("not fished", "fished"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  labs(x = "Depth (m)", y = "Log expected catch (b = 0 fit)") +
  theme_light()
ggsave(paste0(fig_dir, "/negative-b-mu-vs-depth-", cell_km, "km.png"),
  width = 7, height = 5)

# 2. Sampling rate vs expected catch, overall and within depth bands ----------------

bin_rates <- function(x) {
  x |>
    mutate(mu_bin = ntile(log_mu_dev, n_mu_bins)) |>
    group_by(mu_bin) |>
    summarise(log_mu_dev = mean(log_mu_dev), prop_sampled = mean(sampled),
      n = n(), .groups = "drop")
}
rates <- bind_rows(
  bin_rates(fr) |> mutate(depth_band = "All depths"),
  fr |> group_by(depth_band) |> group_modify(~ bin_rates(.x)) |>
    ungroup() |> mutate(depth_band = paste("Depth", depth_band))
) |>
  mutate(depth_band = factor(depth_band, levels = unique(depth_band)))

ggplot(rates, aes(log_mu_dev, prop_sampled)) +
  geom_line() +
  geom_point(aes(size = n)) +
  facet_wrap(~depth_band, nrow = 1) +
  scale_size_area(max_size = 2.5, guide = "none") +
  labs(x = "Log expected catch, deviation from year mean (b = 0 fit)",
    y = "Proportion of cell-years fished",
    subtitle = paste("Depth bands are quintiles of cell depth (Q1 shallowest);",
      "points are deciles of expected catch within each panel")) +
  theme_light()
ggsave(paste0(fig_dir, "/negative-b-sampling-by-depth-", cell_km, "km.png"),
  width = 13, height = 4)

# Slope of logit(sampled) on log expected catch within each depth band
band_slopes <- fr |>
  group_by(depth_band) |>
  group_modify(~ {
    m <- glm(sampled ~ factor(year) + log_mu, family = binomial, data = .x)
    s <- summary(m)$coefficients["log_mu", ]
    tibble(slope = s[["Estimate"]], se = s[["Std. Error"]],
      prop_sampled = mean(.x$sampled), n = nrow(.x))
  }) |>
  ungroup()
cat("\nWithin-depth-band slopes, logit(sampled) ~ factor(year) + log_mu:\n")
print(band_slopes)

# 3. Outside-model logistic regressions -------------------------------------------

glm_slope <- function(form, mu = "log_mu") {
  m <- glm(update(form, paste(". ~ . +", mu)), family = binomial, data = fr)
  s <- summary(m)$coefficients[mu, ]
  tibble(terms = deparse(form[[3]]), surface = mu,
    slope = s[["Estimate"]], se = s[["Std. Error"]], AIC = AIC(m))
}
glm_forms <- list(
  sampled ~ 1,
  sampled ~ factor(year),
  sampled ~ factor(year) + log_dist_port,
  sampled ~ factor(year) + poly(log_depth, 2),
  sampled ~ factor(year) + log_dist_port + poly(log_depth, 2),
  sampled ~ factor(year) + log_dist_port + poly(log_depth, 2) + region
)
glm_slopes <- bind_rows(
  lapply(glm_forms, glm_slope),
  lapply(glm_forms[5], glm_slope, mu = "log_mu_constant")
)
cat("\nSlope of logit(sampled) on log expected catch, by sampling covariates:\n")
print(glm_slopes, width = Inf)

write.csv(band_slopes, file.path(pred_dir, "negative-b-band-slopes.csv"),
  row.names = FALSE)
write.csv(glm_slopes, file.path(pred_dir, "negative-b-glm-slopes.csv"),
  row.names = FALSE)

# 4. Field-level test -------------------------------------------------------------

# Depth's share of b cancels against the sampling model's depth terms
# (R/depth-effects.R), so b must be estimated from the catch model's random
# fields. Field contribution = full minus fixed-effects-only log expected
# catch, split into a cell mean over years (spatial) and the deviation from it
# by year (spatiotemporal). Every frame cell was fished in >= min_yrs_cell
# years, so every cell has catch data; cells fished in more years have
# better-constrained fields.
fields <- fr |>
  mutate(
    cell = paste(X, Y),
    field_zero = log_mu - predict_log_mu(fit_names[["zero"]], fixed = TRUE),
    field_constant = log_mu_constant -
      predict_log_mu(fit_names[["constant"]], fixed = TRUE)
  ) |>
  group_by(cell) |>
  mutate(
    n_fished = sum(sampled),
    across(c(field_zero, field_constant),
      list(cell = mean, dev = ~ .x - mean(.x)))
  ) |>
  ungroup()

# 4a. Between cells: share of years fished vs the cell's mean field, adjusted
# for the sampling model's covariates and region (quasibinomial for
# overdispersion)
cells <- fields |>
  group_by(cell, log_depth, log_dist_port, region) |>
  summarise(n_fished = sum(sampled), n_years = n(),
    zero = first(field_zero_cell), constant = first(field_constant_cell),
    .groups = "drop")
n_cells <- nrow(cells)
between <- bind_rows(lapply(c("zero", "constant"), function(surface) {
  m <- glm(reformulate(c(surface, "poly(log_depth, 2)", "log_dist_port", "region"),
    "cbind(n_fished, n_years - n_fished)"), family = quasibinomial, data = cells)
  s <- summary(m)$coefficients[surface, ]
  tibble(test = "between cells", surface, cells = "all",
    slope = s[["Estimate"]], se = s[["Std. Error"]], n = n_cells)
}))

# 4b. Within cells: fished in a given year vs that year's deviation from the
# cell's mean field, with cell and year fixed effects (cells fished in every
# or no year carry no information and are dropped). Split by how many years
# the cell was fished: conditioning on that count is fine here because the
# cell effects already condition on it.
fields <- fields |>
  mutate(data_group = cut(n_fished, c(0, 7, 15, 29),
    labels = c("fished 3-7 yrs", "fished 8-15 yrs", "fished 16-29 yrs")))
within_slope <- function(x, surface, group) {
  m <- glm(reformulate(c(paste0("field_", surface, "_dev"), "factor(year)", "factor(cell)"),
    "sampled"), family = binomial, data = x)
  s <- summary(m)$coefficients[paste0("field_", surface, "_dev"), ]
  tibble(test = "within cells", surface, cells = group,
    slope = s[["Estimate"]], se = s[["Std. Error"]], n = n_distinct(x$cell))
}
within_data <- filter(fields, n_fished > 0, n_fished < n_distinct(year))
within <- bind_rows(lapply(c("zero", "constant"), function(surface) {
  bind_rows(
    within_slope(within_data, surface, "all"),
    within_data |> group_by(data_group) |> group_split() |>
      lapply(function(x) within_slope(x, surface, as.character(x$data_group[1])))
  )
}))

field_slopes <- bind_rows(between, within) |>
  mutate(z = slope / se)
cat("\nSlope of logit(fished) on the catch model's field contribution:\n")
print(field_slopes, width = Inf)
write.csv(field_slopes, file.path(pred_dir, "negative-b-field-slopes.csv"),
  row.names = FALSE)

# How the constant-b fit reshapes the fields relative to b = 0, by how often
# the cell was fished
cells |>
  mutate(prop_fished = n_fished / n_years) |>
  tidyr::pivot_longer(c(zero, constant), names_to = "surface",
    values_to = "field") |>
  mutate(surface = factor(surface, levels = c("zero", "constant"),
    labels = c("b = 0 fit", "constant-b fit"))) |>
  ggplot(aes(prop_fished, field, colour = surface)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_point(alpha = 0.25, size = 0.7) +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
  scale_colour_manual(values = c("#7570b3", "#d95f02"), name = NULL) +
  labs(x = "Share of years the cell was fished",
    y = "Cell-mean field contribution\n(log expected catch)") +
  theme_light() +
  theme(legend.position = "bottom")
ggsave(paste0(fig_dir, "/negative-b-fields-by-effort-", cell_km, "km.png"),
  width = 7, height = 5)
