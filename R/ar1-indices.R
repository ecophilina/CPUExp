# Index comparison: IID vs AR(1) spatiotemporal catch fields
#
# Bias-corrected indices for the naive (catch-only) and joint (rw b, IID
# sampling baseline) fits, each with IID and AR(1) spatiotemporal fields, over
# the index grid (reference month, unit effort, vessel effects excluded), as in
# R/joint-indices.R. Indices are cached in <species>/output/indices/ under
# the same names as R/joint-indices.R uses. Each index covers only the grid
# cells in its fit's sampling frame (index_grid() in R/species.R).
#
# Joint-fit indices rebuild the full model (~20 GB), so compute one per process:
#   Rscript R/ar1-indices.R <index key>  # compute and cache one index
#   Rscript R/ar1-indices.R              # compute any missing, then compare

library(dplyr)
library(ggplot2)
library(sdmTMB)

source("R/species.R")
mesh_cutoff <- 15
prepped <- read_prepped() # from R/data-prep.R
grid <- prepped$grid
cell_km <- prepped$settings$cell_km

index_dir <- file.path(out_dir, "indices")
dir.create(index_dir, showWarnings = FALSE)

# naive_tv_iid is naive_tv refit with the current fit_cpue_tv() priors, so
# both naive fits share priors (naive_tv predates them)
models <- tribble(
  ~label, ~fit, ~fields,
  "naive", "naive_tv_iid-delta_gengamma", "IID",
  "naive", "naive_tv_ar1-delta_gengamma", "AR(1)",
  "joint, rw b, IID baseline", "rw_iidbase_xi_expected_tv_dgengamma", "IID",
  "joint, rw b, IID baseline", "rw_iidbase_xi_expected_tv_ar1_dgengamma", "AR(1)",
  "joint, constant b, IID baseline", "constant_iidbase_xi_expected_tv_dgengamma", "IID",
  "joint, rw b, IID baseline, footprint frame + area + pre", "rw_areafootprint_iidbase_fpframe_xi_expected_tv_dgengamma", "IID",
  "joint, rw b, IID baseline, footprint frame + area + pre", "rw_areafootprint_iidbase_fpframe_xi_expected_tv_ar1_dgengamma", "AR(1)",
  "joint, constant b, IID baseline, footprint frame + area x pre", "constant_iareafootprint_iidbase_fpframe_xi_expected_tv_dgengamma", "IID",
  "joint, rw b, IID baseline, core footprint frame", "rw_iidbase_fp80frame_xi_expected_tv_dgengamma", "IID",
  "joint, rw b, IID baseline, core footprint frame", "rw_iidbase_fp80frame_xi_expected_tv_ar1_dgengamma", "AR(1)",
  "joint, constant b, IID baseline, core footprint frame", "constant_iidbase_fp80frame_xi_expected_tv_dgengamma", "IID",
  "joint, constant b, IID baseline, core footprint frame", "constant_iidbase_fp80frame_xi_expected_tv_ar1_dgengamma", "AR(1)",
  "joint, constant b, IID baseline, core footprint frame, b on fields", "constant_iidbase_fp80frame_xi_fields_tv_ar1_dgengamma", "AR(1)",
  "joint, constant b, IID baseline, core footprint frame, b on fields", "constant_iidbase_fp80frame_xi_fields_tv_dgengamma", "IID",
  "joint, rw b, IID baseline, core footprint frame, b on fields", "rw_iidbase_fp80frame_xi_fields_tv_ar1_dgengamma", "AR(1)",
  "joint, constant b, IID baseline, core footprint frame, target vessels", "constant_iidbase_fp80frame_target_xi_expected_tv_ar1_dgengamma", "AR(1)",
  "joint, constant b, IID baseline, core footprint frame, target vessels, b on fields", "constant_iidbase_fp80frame_target_xi_fields_tv_ar1_dgengamma", "AR(1)",
  "joint, rw b, rw baseline, area x footprint", "rw_iareafootprint_rwbase_xi_expected_tv_dgengamma", "IID",
  "joint, rw b, rw baseline, area x footprint", "rw_iareafootprint_rwbase_xi_expected_tv_ar1_dgengamma", "AR(1)"
)

# Core-footprint joint fits on the old (5-yr fleet) data: indexed over the
# current core frame (index_grids$fp80), like the fits on the current data, so
# all core indices cover the same cells (old frame: 647 cells, current: 655)
old_data_core_fits <- c(
  "rw_iidbase_fp80frame_xi_expected_tv_dgengamma",
  "rw_iidbase_fp80frame_xi_expected_tv_ar1_dgengamma",
  "constant_iidbase_fp80frame_xi_expected_tv_dgengamma",
  "constant_iidbase_fp80frame_xi_expected_tv_ar1_dgengamma",
  "constant_iidbase_fp80frame_xi_fields_tv_ar1_dgengamma",
  "rw_iidbase_fp80frame_xi_fields_tv_ar1_dgengamma",
  "constant_iidbase_fp80frame_xi_fields_tv_dgengamma"
)
models$index_frame <- ifelse(models$fit %in% old_data_core_fits, "fp80", NA)

# Naive fits indexed over the core footprint cells only (same fits, trimmed
# grid), to compare with the core-footprint joint fits. index_frame names an
# entry of index_grids; NA = the fit's own frame (index_grid()). key names the
# cached index: the fit name, plus "-<index_frame>" when there is one.
models <- bind_rows(models, tibble(
  label = "naive, core footprint grid",
  fit = c("naive_tv_iid-delta_gengamma", "naive_tv_ar1-delta_gengamma"),
  fields = c("IID", "AR(1)"), index_frame = "fp80"
), tibble(
  # catch model of the target-vessel joint fits (target_vessel effect), on the
  # current data; compare with those fits
  label = "naive, core footprint grid, target vessels",
  fit = "naive_tv_ar1_target-delta_gengamma", fields = "AR(1)", index_frame = "fp80"
)) |>
  mutate(key = ifelse(is.na(index_frame), fit, paste0(fit, "-", index_frame)))

# Cells at least 80% inside the official footprint (frame_fp80 in R/cpue-index.R),
# always from the current prepped data, even when CPUE_DATA selects older data
current <- readRDS(file.path(data_dir, "cpue-index-clean.rds"))
index_grids <- list(
  fp80 = semi_join(current$grid,
    distinct(filter(current$frame, footprint_area >= 0.8), X, Y), by = c("X", "Y"))
)
rm(current)

# Joint fits may sit in a subfolder of models/ (e.g. "6km footprint models")
fit_path <- function(fit_name) find_fit(fit_name, mesh_cutoff)
index_path <- function(key) {
  file.path(index_dir, paste0("index-", key, "-cutoff", mesh_cutoff, ".rds"))
}

compute_index <- function(key) {
  message("Computing index: ", key)
  m <- models[models$key == key, ]
  fit_name <- m$fit
  fit <- readRDS(fit_path(fit_name))
  if (!isTRUE(sanity(fit, silent = TRUE)$all_ok)) {
    warning(fit_name, " does not pass sanity()")
  }
  t0 <- Sys.time()
  newdata <- if (is.na(m$index_frame)) index_grid(fit, grid) else index_grids[[m$index_frame]]
  idx <- get_index(fit, newdata = newdata, area = cell_km^2,
    bias_correct = TRUE, # keep TRUE for joint delta models (see docs)
    predict_args = list(re_form_iid = NA))
  message("Index done: ", key, " in ",
    round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
  saveRDS(idx, index_path(key))
  invisible(idx)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  stopifnot(args[1] %in% models$key)
  compute_index(args[1])
  quit(save = "no")
}

# Compare only the models this species has fits (or cached indices) for
has_fit <- function(fit_name) {
  !inherits(try(fit_path(fit_name), silent = TRUE), "try-error")
}
available <- file.exists(index_path(models$key)) |
  vapply(models$fit, has_fit, logical(1))
if (any(!available)) message("No fit, skipping: ", paste(models$key[!available], collapse = ", "))
models <- models[available, ]

for (k in models$key) if (!file.exists(index_path(k))) compute_index(k)

indices <- models |>
  rowwise() |>
  mutate(index = list(readRDS(index_path(key)))) |>
  ungroup() |>
  tidyr::unnest(index) |>
  group_by(key) |>
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
    title = tools::toTitleCase(species),
    caption = paste0("Ribbons: 95% CI. Delta-gengamma, rw0 catch intercept; ",
      cell_km, " km index grid.")) +
  theme_light() +
  theme(legend.position = "top")
ggsave(paste0(fig_dir, "/index-ar1-comparison-", cell_km, "km.png"), width = 9,
  height = 2 + 2.7 * n_distinct(models$label))

# Same indices overlaid by catch-field structure: one panel for IID fields,
# one for AR(1). Free y scales, as AR(1) joint indices can run much higher.
model_cols <- c("#000000", "#0072B2", "#D55E00", "#009E73", "#CC79A7",
  "#E69F00", "#56B4E9", "#999999", "#882255", "#44AA99", "#AA4499",
  "#332288", "#DDCC77", "#661100")
stopifnot(length(levels(indices$label)) <= length(model_cols))
model_cols <- setNames(model_cols[seq_along(levels(indices$label))],
  levels(indices$label))
ggplot(indices, aes(year, rel_est, colour = label, fill = label)) +
  geom_ribbon(aes(ymin = rel_lwr, ymax = rel_upr), alpha = 0.08, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~fields, ncol = 1, scales = "free_y",
    labeller = labeller(fields = function(x) paste(x, "catch fields"))) +
  scale_colour_manual(values = model_cols) +
  scale_fill_manual(values = model_cols) +
  guides(colour = guide_legend(ncol = 1), fill = guide_legend(ncol = 1)) +
  labs(x = "Year", y = "Relative index (geometric mean = 1)",
    colour = NULL, fill = NULL,
    title = tools::toTitleCase(species),
    caption = paste0("Ribbons: 95% CI. Delta-gengamma, rw0 catch intercept; ",
      cell_km, " km index grid.")) +
  theme_light() +
  theme(legend.position = "top", legend.justification = "left")
ggsave(paste0(fig_dir, "/index-by-fields-", cell_km, "km.png"), width = 9,
  height = 8 + 0.25 * n_distinct(models$label))

# Core footprint only: naive and joint indices over the same cells (>= 80%
# inside the official footprint), same colours as above
core <- filter(indices, grepl("core footprint", label)) |> droplevels()
if (nrow(core)) {
  ggplot(core, aes(year, rel_est, colour = label, fill = label)) +
    geom_ribbon(aes(ymin = rel_lwr, ymax = rel_upr), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~fields, ncol = 1, scales = "free_y",
      labeller = labeller(fields = function(x) paste(x, "catch fields"))) +
    scale_colour_manual(values = model_cols) +
    scale_fill_manual(values = model_cols) +
    guides(colour = guide_legend(ncol = 1), fill = guide_legend(ncol = 1)) +
    labs(x = "Year", y = "Relative index (geometric mean = 1)",
      colour = NULL, fill = NULL,
      title = paste0(tools::toTitleCase(species), ": core footprint"),
      caption = paste0("Index over cells >= 80% inside the 2012 trawl footprint. ",
        "Ribbons: 95% CI. Delta-gengamma, rw0 catch intercept; ",
        cell_km, " km grid.")) +
    theme_light() +
    theme(legend.position = "top", legend.justification = "left")
  ggsave(paste0(fig_dir, "/index-core-footprint-", cell_km, "km.png"),
    width = 9, height = 8)
}

# Core footprint: b on expected catch vs b on the catch fields, against the
# matching naive fit. Panels: old (5-yr fleet) data with IID and with AR(1)
# catch fields; current data with target vessels (target_vessel in the catch
# model, sampled_target response), AR(1) catch fields. All quadratic sampling
# depth.
target_keys <- tribble(
  ~panel, ~model, ~key,
  "Old data (5-yr fleet), IID catch fields", "naive", "naive_tv_iid-delta_gengamma-fp80",
  "Old data (5-yr fleet), IID catch fields", "joint, b on expected catch", "constant_iidbase_fp80frame_xi_expected_tv_dgengamma-fp80",
  "Old data (5-yr fleet), IID catch fields", "joint, b on catch fields", "constant_iidbase_fp80frame_xi_fields_tv_dgengamma-fp80",
  "Old data (5-yr fleet), AR(1) catch fields", "naive", "naive_tv_ar1-delta_gengamma-fp80",
  "Old data (5-yr fleet), AR(1) catch fields", "joint, b on expected catch", "constant_iidbase_fp80frame_xi_expected_tv_ar1_dgengamma-fp80",
  "Old data (5-yr fleet), AR(1) catch fields", "joint, b on catch fields", "constant_iidbase_fp80frame_xi_fields_tv_ar1_dgengamma-fp80",
  "Current data, target vessels, AR(1) catch fields", "naive", "naive_tv_ar1_target-delta_gengamma-fp80",
  "Current data, target vessels, AR(1) catch fields", "joint, b on expected catch", "constant_iidbase_fp80frame_target_xi_expected_tv_ar1_dgengamma",
  "Current data, target vessels, AR(1) catch fields", "joint, b on catch fields", "constant_iidbase_fp80frame_target_xi_fields_tv_ar1_dgengamma"
)
target_cols <- c("naive" = "#4D4D4D", "joint, b on expected catch" = "#D55E00",
  "joint, b on catch fields" = "#0072B2")
by_target <- inner_join(target_keys, select(indices, -label, -fields), by = "key") |>
  mutate(panel = factor(panel, unique(target_keys$panel)),
    model = factor(model, names(target_cols)))
if (nrow(by_target)) {
  ggplot(by_target, aes(year, rel_est, colour = model, fill = model)) +
    geom_ribbon(aes(ymin = rel_lwr, ymax = rel_upr), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~panel, ncol = 1) +
    scale_colour_manual(values = target_cols) +
    scale_fill_manual(values = target_cols) +
    labs(x = "Year", y = "Relative index (geometric mean = 1)",
      colour = NULL, fill = NULL,
      title = paste0(tools::toTitleCase(species),
        ": preference for expected catch vs catch fields"),
      caption = paste0("Constant b, IID sampling baseline. ",
        "Index over cells >= 80% inside the 2012 trawl footprint (current frame).\n",
        "Ribbons: 95% CI. Delta-gengamma, rw0 catch intercept; ", cell_km, " km grid.")) +
    theme_light() +
    theme(legend.position = "top", legend.justification = "left")
  ggsave(paste0(fig_dir, "/index-target-fields-vs-expected-", cell_km, "km.png"),
    width = 9, height = 2 + 3 * n_distinct(by_target$panel))

  # Same, unscaled: est is predicted catch (kg) per hour of effort summed over
  # the cells times cell area (reference month, vessel effects excluded).
  # Shared y scale, so levels compare across panels.
  ggplot(by_target, aes(year, est, colour = model, fill = model)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~panel, ncol = 1) +
    scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    scale_colour_manual(values = target_cols) +
    scale_fill_manual(values = target_cols) +
    labs(x = "Year", y = "Index (predicted kg/h x cell area in km², summed over cells)",
      colour = NULL, fill = NULL,
      title = paste0(tools::toTitleCase(species),
        ": preference for expected catch vs catch fields (unscaled)"),
      caption = paste0("Constant b, IID sampling baseline. ",
        "Index over cells >= 80% inside the 2012 trawl footprint (current frame).\n",
        "Ribbons: 95% CI. Delta-gengamma, rw0 catch intercept; ", cell_km, " km grid.")) +
    theme_light() +
    theme(legend.position = "top", legend.justification = "left")
  ggsave(paste0(fig_dir, "/index-target-fields-vs-expected-absolute-", cell_km, "km.png"),
    width = 9, height = 2 + 3 * n_distinct(by_target$panel))
}

# Log-linear trend: the quantity preferential sampling can bias
trends <- indices |>
  group_by(label, fields, key) |>
  summarise(
    slope = coef(lm(log_est ~ year))[[2]],
    change_30y = exp(slope * (max(year) - min(year))) - 1,
    mean_cv = mean(sqrt(exp(se^2) - 1)),
    .groups = "drop"
  ) |>
  arrange(label, fields)
print(trends, width = Inf)

saveRDS(indices, file.path(out_dir, "index-ar1-comparison.rds"))
write.csv(trends, file.path(out_dir, "index-ar1-trends.csv"), row.names = FALSE)
