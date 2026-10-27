# Index comparison: IID vs AR(1) spatiotemporal catch fields
#
# Bias-corrected indices for the naive (catch-only) and joint (rw b, IID
# sampling baseline) fits, each with IID and AR(1) spatiotemporal fields, over
# the index grid (reference month, unit effort, vessel effects excluded), as in
# R/joint-indices.R. Indices are cached in output/indices/ under the same names
# as R/joint-indices.R uses.
#
# Joint-fit indices rebuild the full model (~20 GB), so compute one per process:
#   Rscript R/ar1-indices.R <fit name>   # compute and cache one index
#   Rscript R/ar1-indices.R              # compute any missing, then compare

library(dplyr)
library(ggplot2)
library(sdmTMB)

mesh_cutoff <- 15
prepped <- readRDS("data/cpue-index-clean.rds") # from R/data-prep.R
grid <- prepped$grid
cell_km <- prepped$settings$cell_km

index_dir <- file.path("output", "indices")
dir.create(index_dir, showWarnings = FALSE)

# naive_tv_iid is naive_tv refit with the current fit_cpue_tv() priors, so
# both naive fits share priors (naive_tv predates them)
models <- tribble(
  ~label, ~fit, ~fields,
  "naive", "naive_tv_iid-delta_gengamma", "IID",
  "naive", "naive_tv_ar1-delta_gengamma", "AR(1)",
  "joint, rw b, IID baseline", "rw_iidbase_xi_expected_tv_dgengamma", "IID",
  "joint, rw b, IID baseline", "rw_iidbase_xi_expected_tv_ar1_dgengamma", "AR(1)",
  "joint, rw b, rw baseline, area x footprint", "rw_iareafootprint_rwbase_xi_expected_tv_dgengamma", "IID",
  "joint, rw b, rw baseline, area x footprint", "rw_iareafootprint_rwbase_xi_expected_tv_ar1_dgengamma", "AR(1)"
)

# Joint fits may sit in a subfolder of models/ (e.g. "6km footprint models")
fit_path <- function(fit_name) {
  file_name <- paste0("fit-", fit_name, "-cutoff", mesh_cutoff, ".rds")
  if (grepl("^naive", fit_name)) return(file.path("output", file_name))
  f <- list.files("models", recursive = TRUE, full.names = TRUE)
  f <- f[basename(f) == file_name]
  stopifnot(length(f) == 1)
  f
}
index_path <- function(fit_name) {
  file.path(index_dir, paste0("index-", fit_name, "-cutoff", mesh_cutoff, ".rds"))
}

compute_index <- function(fit_name) {
  message("Computing index: ", fit_name)
  fit <- readRDS(fit_path(fit_name))
  if (!isTRUE(sanity(fit, silent = TRUE)$all_ok)) {
    warning(fit_name, " does not pass sanity()")
  }
  t0 <- Sys.time()
  idx <- get_index(fit, newdata = grid, area = cell_km^2,
    bias_correct = TRUE, # keep TRUE for joint delta models (see docs)
    predict_args = list(re_form_iid = NA))
  message("Index done: ", fit_name, " in ",
    round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
  saveRDS(idx, index_path(fit_name))
  invisible(idx)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  stopifnot(args[1] %in% models$fit)
  compute_index(args[1])
  quit(save = "no")
}

for (f in models$fit) if (!file.exists(index_path(f))) compute_index(f)

indices <- models |>
  rowwise() |>
  mutate(index = list(readRDS(index_path(fit)))) |>
  ungroup() |>
  tidyr::unnest(index) |>
  group_by(fit) |>
  mutate(
    geo_mean = exp(mean(log_est)),
    rel_est = est / geo_mean,
    rel_lwr = lwr / geo_mean,
    rel_upr = upr / geo_mean
  ) |>
  ungroup() |>
  mutate(fields = factor(fields, c("IID", "AR(1)")),
    label = factor(label, unique(models$label)))

ggplot(indices, aes(year, rel_est, colour = fields, fill = fields)) +
  geom_ribbon(aes(ymin = rel_lwr, ymax = rel_upr), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~label, ncol = 1) +
  scale_colour_manual(values = c("IID" = "#0072B2", "AR(1)" = "#D55E00")) +
  scale_fill_manual(values = c("IID" = "#0072B2", "AR(1)" = "#D55E00")) +
  labs(x = "Year", y = "Relative index (geometric mean = 1)",
    colour = "Spatiotemporal catch fields", fill = "Spatiotemporal catch fields",
    caption = "Ribbons: 95% CI. Delta-gengamma, rw0 catch intercept; 6 km index grid.") +
  theme_light() +
  theme(legend.position = "top")
ggsave(paste0("figs/index-ar1-comparison-", cell_km, "km.png"), width = 9, height = 10)

# Log-linear trend: the quantity preferential sampling can bias
trends <- indices |>
  group_by(label, fields, fit) |>
  summarise(
    slope = coef(lm(log_est ~ year))[[2]],
    change_30y = exp(slope * (max(year) - min(year))) - 1,
    mean_cv = mean(sqrt(exp(se^2) - 1)),
    .groups = "drop"
  ) |>
  arrange(label, fields)
print(trends, width = Inf)

saveRDS(indices, file.path("output", "index-ar1-comparison.rds"))
write.csv(trends, file.path("output", "index-ar1-trends.csv"), row.names = FALSE)
