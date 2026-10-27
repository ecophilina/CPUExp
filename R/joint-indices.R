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

# label = fit name (output/ or models/fit-<name>-cutoff<mesh_cutoff>.rds; see
# fit_path()). Joint fits use their converged (_newton) versions where a
# Newton step was needed.
#
# 6 km frame: all fits are delta-gengamma with the time-varying (rw0) catch
# intercept; joint fits have the rw sampling baseline and sampling field.
# `sampling` = the sampling-model covariates beyond depth and port distance,
# one panel each (NA for the naive fit, which is drawn in every panel).
# Not included: ifootprint_rwbase_xi_expected_tv_dgengamma (complete
# separation in the data-driven footprint interaction; fails sanity()).
# 4 km frame fits are in models/4km frame/ and output/4km frame/.
models <- tribble(
  ~label, ~fit, ~sampling, ~model_type,
  "naive", "naive_tv-delta_gengamma", NA, "naive",
  "b = 0", "zero_rwbase_xi_expected_tv_dgengamma", "depth + port", "joint, b = 0",
  "constant b", "constant_rwbase_xi_expected_tv_dgengamma", "depth + port", "joint",
  "b = 0", "zero_rwbase_region_xi_expected_tv_dgengamma", "+ region", "joint, b = 0",
  "constant b", "constant_rwbase_region_xi_expected_tv_dgengamma", "+ region", "joint",
  "constant b, data-driven footprint", "footprint_rwbase_xi_expected_tv_dgengamma",
  "+ footprint", "joint",
  "b = 0, official footprint area x pre", "zero_iareafootprint_rwbase_xi_expected_tv_dgengamma",
  "+ footprint", "joint, b = 0",
  "constant b, official footprint area x pre", "iareafootprint_rwbase_xi_expected_tv_dgengamma",
  "+ footprint", "joint",
  # Depth + port, sampling frame from footprint_year (2013) on only
  "b = 0, frame 2013+", "zero_rwbase_post2013_xi_expected_tv_dgengamma",
  "frame 2013+ only", "joint, b = 0",
  "constant b, frame 2013+", "constant_rwbase_post2013_xi_expected_tv_dgengamma",
  "frame 2013+ only", "joint"
)

# Naive indices were already computed by R/family-comparison.R
# (naive-<family> -> family-checks-<family>; naive_tv-<family> -> family-checks_tv-<family>)
family_check_file <- function(fit_name) {
  tag <- sub("^naive(_[a-z0-9]+)?-.*$", "\\1", fit_name)
  fam <- sub("^naive(_[a-z0-9]+)?-", "", fit_name)
  file.path("output", "family-checks",
    paste0("family-checks", tag, "-", fam, "-cutoff", mesh_cutoff, ".rds"))
}

fit_path <- function(fit_name) {
  dir <- if (grepl("^naive", fit_name)) "output" else "models"
  file.path(dir, paste0("fit-", fit_name, "-cutoff", mesh_cutoff, ".rds"))
}

get_cached_index <- function(fit_name) {
  # Use the converged version from R/extra-optimization.R if there is one
  if (!grepl("_newton$", fit_name) && file.exists(fit_path(paste0(fit_name, "_newton")))) {
    fit_name <- paste0(fit_name, "_newton")
  }
  f <- file.path(index_dir, paste0("index-", fit_name, "-cutoff", mesh_cutoff, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  fc <- family_check_file(fit_name)
  if (grepl("^naive(_[a-z0-9]+)?-", fit_name) && file.exists(fc)) {
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

# The naive index has no `sampling` value, so dropping that column repeats it
# in every panel as the reference
joint <- filter(indices, !is.na(sampling)) |>
  mutate(sampling = factor(sampling, levels = unique(models$sampling[!is.na(models$sampling)])))
naive <- filter(indices, is.na(sampling)) |> select(-sampling)
ggplot(joint, aes(year, rel_est, colour = label, fill = label)) +
  geom_ribbon(aes(year, ymin = rel_lwr, ymax = rel_upr), data = naive,
    inherit.aes = FALSE, fill = "grey50", alpha = 0.15) +
  geom_line(aes(year, rel_est), data = naive, inherit.aes = FALSE,
    colour = "grey30", linewidth = 0.8) +
  geom_ribbon(aes(ymin = rel_lwr, ymax = rel_upr), alpha = 0.1, colour = NA) +
  geom_line(aes(linetype = model_type)) +
  facet_wrap(~sampling, ncol = 1) +
  scale_linetype_manual(values = c("joint" = "solid", "joint, b = 0" = "22"),
    guide = "none") +
  labs(x = "Year", y = "Relative index (geometric mean = 1)",
    colour = NULL, fill = NULL,
    caption = "Grey: naive (catch-only) index. Dashed: b = 0 joint fits.") +
  theme_light() +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(ncol = 2), fill = guide_legend(ncol = 2))
ggsave(paste0("figs/index-joint-comparison-", cell_km, "km.png"), width = 9, height = 10)

# Log-linear trend: the quantity preferential sampling can bias (the scale is
# absorbed by catchability in the assessment)
trends <- indices |>
  group_by(sampling, label, model_type, fit) |>
  summarise(
    slope = coef(lm(log_est ~ year))[[2]],
    change_30y = exp(slope * (max(year) - min(year))) - 1,
    mean_cv = mean(sqrt(exp(se^2) - 1)),
    .groups = "drop"
  ) |>
  arrange(sampling, model_type)
print(trends, width = Inf)

saveRDS(indices, file.path("output", "index-joint-comparison.rds"))
write.csv(trends, file.path("output", "index-joint-trends.csv"), row.names = FALSE)
