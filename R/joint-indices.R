# Index comparison: naive vs joint preferential-sampling fits
#
# Computes the bias-corrected index for each cached fit over the footprint grid
# (reference month, unit effort, vessel effects excluded), then compares
# relative indices and log-linear trends. Indices are cached per fit in
# output/indices/, so adding a fit to `models` only computes the new one.
#
# Usage (from the project root): Rscript R/joint-indices.R

library(dplyr)
library(ggplot2)
library(sdmTMB)

mesh_cutoff <- 15
prepped <- readRDS("data/cpue-index-clean.rds") # from R/data-prep.R
grid <- prepped$grid
cell_km <- prepped$settings$cell_km

index_dir <- file.path("output", "indices")
dir.create(index_dir, showWarnings = FALSE)

# label = fit name (output/fit-<name>-cutoff<mesh_cutoff>.rds). Joint fits use
# their converged (_newton) versions where a Newton step was needed.
models <- tribble(
  ~label, ~fit, ~family, ~model_type,
  "naive", "naive-delta_gamma", "delta-Gamma", "naive",
  "zero", "zero_xi_expected_dgamma_newton", "delta-Gamma", "joint, b = 0",
  "constant (fields)", "constant_xi_newton", "delta-Gamma", "joint",
  "constant (expected)", "constant_xi_expected_dgamma_newton", "delta-Gamma", "joint",
  "constant, rw baseline", "constant_rwbase_xi", "delta-Gamma", "joint",
  "naive", "naive-delta_lognormal", "delta-lognormal", "naive",
  "zero", "zero_xi_expected_dlnorm_newton", "delta-lognormal", "joint, b = 0",
  "constant (expected)", "constant_xi_expected_dlnorm_newton", "delta-lognormal", "joint",
  "constant, rw baseline", "constant_rwbase_xi_expected_dlnorm", "delta-lognormal", "joint",
  # The expected-target version ran out of memory (~16 GB+); fields instead
  "rw b_t, rw baseline (fields)", "rw_rwbase_xi_fields_dlnorm", "delta-lognormal", "joint",
  "naive", "naive-delta_gengamma", "delta-gengamma", "naive",
  "constant (expected)", "constant_xi_expected_dgengamma", "delta-gengamma", "joint"
)

# Naive indices were already computed by R/family-comparison.R
family_check_file <- function(fit_name) {
  fam <- sub("^naive-", "", fit_name)
  file.path("output", "family-checks",
    paste0("family-checks-", fam, "-cutoff", mesh_cutoff, ".rds"))
}

fit_path <- function(fit_name) {
  file.path("output", paste0("fit-", fit_name, "-cutoff", mesh_cutoff, ".rds"))
}

get_cached_index <- function(fit_name) {
  # Use the converged version from R/extra-optimization.R if there is one
  if (!grepl("_newton$", fit_name) && file.exists(fit_path(paste0(fit_name, "_newton")))) {
    fit_name <- paste0(fit_name, "_newton")
  }
  f <- file.path(index_dir, paste0("index-", fit_name, "-cutoff", mesh_cutoff, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  fc <- family_check_file(fit_name)
  if (grepl("^naive-", fit_name) && file.exists(fc)) {
    idx <- readRDS(fc)$index
  } else {
    message("Computing index: ", fit_name)
    fit <- readRDS(fit_path(fit_name))
    if (!isTRUE(sanity(fit, silent = TRUE)$all_ok)) {
      warning(fit_name, " does not pass sanity()")
    }
    idx <- get_index(fit, newdata = grid, area = cell_km^2,
      bias_correct = TRUE, # keep TRUE for joint delta models (see docs)
      predict_args = list(re_form_iid = NA))
  }
  saveRDS(idx, f)
  idx
}

indices <- models |>
  rowwise() |>
  mutate(index = list(get_cached_index(fit))) |>
  ungroup() |>
  tidyr::unnest(index) |>
  group_by(fit) |>
  mutate(
    geo_mean = exp(mean(log_est)),
    rel_est = est / geo_mean,
    rel_lwr = lwr / geo_mean,
    rel_upr = upr / geo_mean
  ) |>
  ungroup()

ggplot(indices, aes(year, rel_est, colour = label, fill = label)) +
  geom_ribbon(aes(ymin = rel_lwr, ymax = rel_upr), alpha = 0.1, colour = NA) +
  geom_line() +
  facet_wrap(~family, ncol = 1) +
  labs(x = "Year", y = "Relative index (geometric mean = 1)",
    colour = NULL, fill = NULL) +
  theme_light()
ggsave("figs/index-joint-comparison.png", width = 9, height = 8)

# Log-linear trend: the quantity preferential sampling can bias (the scale is
# absorbed by catchability in the assessment)
trends <- indices |>
  group_by(family, label, model_type, fit) |>
  summarise(
    slope = coef(lm(log_est ~ year))[[2]],
    change_30y = exp(slope * (max(year) - min(year))) - 1,
    mean_cv = mean(sqrt(exp(se^2) - 1)),
    .groups = "drop"
  ) |>
  arrange(family, model_type)
print(trends, width = Inf)

saveRDS(indices, file.path("output", "index-joint-comparison.rds"))
write.csv(trends, file.path("output", "index-joint-trends.csv"), row.names = FALSE)
