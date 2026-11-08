# Index comparison: IID vs AR(1) spatiotemporal catch fields
#
# Bias-corrected indices for the naive (catch-only) and joint (rw b, IID
# sampling baseline) fits, each with IID and AR(1) spatiotemporal fields, over
# the index grid (reference month, unit effort, vessel effects excluded), as in
# R/06-joint-indices.R. Indices are cached in <species>/output/indices/ under
# the same names as R/06-joint-indices.R uses. Each index covers only the grid
# cells in its fit's sampling frame (index_grid() in R/species.R).
#
# Joint-fit indices rebuild the full model (~20 GB), so compute one per process:
#   Rscript R/06-ar1-indices.R <index key>  # compute and cache one index
#   Rscript R/06-ar1-indices.R              # compute any missing, then compare

library(dplyr)
library(ggplot2)
library(sdmTMB)

source("R/species.R")
mesh_cutoff <- 15
prepped <- read_prepped() # from R/01-data-prep.R
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
  "joint, constant b, IID baseline, core footprint frame, target-vessel sampling", "constant_iidbase_fp80frame_target_xi_expected_tv_ar1_dgengamma", "AR(1)",
  "joint, constant b, IID baseline, core footprint frame, target-vessel sampling, b on fields", "constant_iidbase_fp80frame_target_xi_fields_tv_ar1_dgengamma", "AR(1)",
  "joint, constant b, IID baseline, core footprint frame, all-vessel sampling", "constant_iidbase_fp80frame_allves_xi_expected_tv_ar1_dgengamma", "AR(1)",
  "joint, constant b, IID baseline, core footprint frame, target-vessel sampling", "constant_iidbase_fp80frame_target_xi_expected_tv_dgengamma", "IID",
  "joint, constant b, IID baseline, core footprint frame, target-vessel sampling, b on fields", "constant_iidbase_fp80frame_target_xi_fields_tv_dgengamma", "IID",
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
  fit = c("naive_tv_iid_target-delta_gengamma", "naive_tv_ar1_target-delta_gengamma"),
  fields = c("IID", "AR(1)"), index_frame = "fp80"
), tibble(
  # old-data catch model (no target_vessel) on the current data: separates the
  # data change from the target_vessel term
  label = "naive, core footprint grid, current data, no target_vessel",
  fit = "naive_tv_ar1_current-delta_gengamma", fields = "AR(1)", index_frame = "fp80"
), tibble(
  # target-vessel joint fit on the min-3-year fleet data
  # (CPUE_DATA=cpue-index-clean-min3yr.rds, saved with --out under this name)
  label = "joint, constant b, IID baseline, core footprint frame, target-vessel sampling, min 3 yr fleet",
  fit = "constant_iidbase_fp80frame_target_xi_expected_tv_ar1_dgengamma-min3yr",
  fields = "AR(1)", index_frame = "fp80"
), tibble(
  label = "naive, core footprint grid, target vessels, min 3 yr fleet",
  fit = "naive_tv_ar1_target_min3yr-delta_gengamma", fields = "AR(1)", index_frame = "fp80"
), tibble(
  # old-data naive model refitted with the current code on the min-5-year data
  # (same tows as the old data)
  label = "naive, core footprint grid, min 5 yr fleet, current code",
  fit = "naive_tv_ar1_min5yr-delta_gengamma", fields = "AR(1)", index_frame = "fp80"
), tibble(
  label = "naive, core footprint grid, min 3 yr fleet, no target_vessel",
  fit = "naive_tv_ar1_min3yr-delta_gengamma", fields = "AR(1)", index_frame = "fp80"
), tibble(
  label = "naive, core footprint grid, target vessels, min 5 yr fleet",
  fit = "naive_tv_ar1_target_min5yr-delta_gengamma", fields = "AR(1)", index_frame = "fp80"
), tibble(
  label = "joint, constant b, IID baseline, core footprint frame, target-vessel sampling, min 5 yr fleet",
  fit = "constant_iidbase_fp80frame_target_xi_expected_tv_ar1_dgengamma-min5yr",
  fields = "AR(1)", index_frame = "fp80"
)) |>
  mutate(key = ifelse(is.na(index_frame), fit, paste0(fit, "-", index_frame)))

# Cells at least 80% inside the official footprint (frame_fp80 in R/02-cpue-index.R),
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
# more models than colours: extend with extra hues (these overview plots get
# crowded; the fleet comparisons have their own plot below)
n_extra <- length(levels(indices$label)) - length(model_cols)
if (n_extra > 0) model_cols <- c(model_cols, hcl.colors(n_extra, "Dark 3"))
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
# model, sampled_target response), IID and AR(1) catch fields. All quadratic sampling
# depth.
target_keys <- tribble(
  ~panel, ~model, ~key,
  "Old data (5-yr fleet), IID catch fields", "naive", "naive_tv_iid-delta_gengamma-fp80",
  "Old data (5-yr fleet), IID catch fields", "joint, b on expected catch", "constant_iidbase_fp80frame_xi_expected_tv_dgengamma-fp80",
  "Old data (5-yr fleet), IID catch fields", "joint, b on catch fields", "constant_iidbase_fp80frame_xi_fields_tv_dgengamma-fp80",
  "Old data (5-yr fleet), AR(1) catch fields", "naive", "naive_tv_ar1-delta_gengamma-fp80",
  "Old data (5-yr fleet), AR(1) catch fields", "joint, b on expected catch", "constant_iidbase_fp80frame_xi_expected_tv_ar1_dgengamma-fp80",
  "Old data (5-yr fleet), AR(1) catch fields", "joint, b on catch fields", "constant_iidbase_fp80frame_xi_fields_tv_ar1_dgengamma-fp80",
  "Current data, target-vessel sampling, IID catch fields", "naive", "naive_tv_iid_target-delta_gengamma-fp80",
  "Current data, target-vessel sampling, IID catch fields", "joint, b on expected catch", "constant_iidbase_fp80frame_target_xi_expected_tv_dgengamma",
  "Current data, target-vessel sampling, IID catch fields", "joint, b on catch fields", "constant_iidbase_fp80frame_target_xi_fields_tv_dgengamma",
  "Current data, target-vessel sampling, AR(1) catch fields", "naive", "naive_tv_ar1_target-delta_gengamma-fp80",
  "Current data, target-vessel sampling, AR(1) catch fields", "joint, b on expected catch", "constant_iidbase_fp80frame_target_xi_expected_tv_ar1_dgengamma",
  "Current data, target-vessel sampling, AR(1) catch fields", "joint, b on catch fields", "constant_iidbase_fp80frame_target_xi_fields_tv_ar1_dgengamma",
  "Current data, all-vessel sampling, AR(1) catch fields", "joint, b on expected catch", "constant_iidbase_fp80frame_allves_xi_expected_tv_ar1_dgengamma"
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
        "Catch models use all tows; current-data ones add a target_vessel catchability effect (index at target_vessel = 0).\n",
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
        "Catch models use all tows; current-data ones add a target_vessel catchability effect (index at target_vessel = 0).\n",
        "Ribbons: 95% CI. Delta-gengamma, rw0 catch intercept; ", cell_km, " km grid.")) +
    theme_light() +
    theme(legend.position = "top", legend.justification = "left")
  ggsave(paste0(fig_dir, "/index-target-fields-vs-expected-absolute-", cell_km, "km.png"),
    width = 9, height = 2 + 3 * n_distinct(by_target$panel))

  # Same, one panel per plot: colour = data/fields (fixed per data set across
  # plots), alpha = model (opaque naive, b on fields fainter, b on expected
  # faintest). No ribbons (overlapping CIs are unreadable; see the faceted
  # version). Plots: old vs current data, AR(1) catch fields; current data,
  # IID vs AR(1) catch fields. Current-data naive fits have no sampling model
  # (shared by target- and all-vessel sampling), so they get their own colour
  panel_cols <- c(
    "Old data (5-yr fleet), IID catch fields" = "#009E73",
    "Old data (5-yr fleet), AR(1) catch fields" = "#D55E00",
    "Current data, IID catch fields" = "#999999",
    "Current data, AR(1) catch fields" = "#4D4D4D",
    "Current data, target-vessel sampling, IID catch fields" = "#CC79A7",
    "Current data, target-vessel sampling, AR(1) catch fields" = "#0072B2",
    "Current data, all-vessel sampling, AR(1) catch fields" = "#E69F00"
  )
  model_alpha <- c("naive" = 1, "joint, b on expected catch" = 0.3,
    "joint, b on catch fields" = 0.6)
  one_panel_plot <- function(panels, file) {
    d <- filter(by_target, panel %in% panels) |>
      mutate(panel = as.character(panel),
        panel = if_else(model == "naive" & startsWith(panel, "Current data"),
          sub("target-vessel sampling, ", "", panel), panel),
        panel = factor(panel, intersect(names(panel_cols), panel)))
    if (!nrow(d)) return(invisible())
    ggplot(d, aes(year, est, colour = panel, alpha = model)) +
      geom_line(linewidth = 0.9) +
      scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
      scale_colour_manual(values = panel_cols) +
      scale_alpha_manual(values = model_alpha) +
      labs(x = "Year", y = "Index (predicted kg/h x cell area in km², summed over cells)",
        colour = NULL, alpha = NULL,
        title = paste0(tools::toTitleCase(species),
          ": preference for expected catch vs catch fields (unscaled)"),
        caption = paste0("Constant b, IID sampling baseline. ",
          "Index over cells >= 80% inside the 2012 trawl footprint (current frame).\n",
          "Catch models use all tows; current-data ones add a target_vessel catchability effect (index at target_vessel = 0).\n",
          "Delta-gengamma, rw0 catch intercept; ", cell_km, " km grid.")) +
      guides(alpha = guide_legend(ncol = 1, order = 1,
          override.aes = list(colour = "grey20")),
        colour = guide_legend(ncol = 1, order = 2)) +
      theme_light() +
      theme(legend.position = "top", legend.justification = "left",
        legend.key.width = unit(2.5, "lines"))
    ggsave(paste0(fig_dir, "/", file, "-", cell_km, "km.png"), width = 9, height = 6.5)
  }
  one_panel_plot(c("Current data, target-vessel sampling, IID catch fields",
    "Current data, target-vessel sampling, AR(1) catch fields"),
    "index-target-iid-vs-ar1-absolute-1panel")

  # Old vs current data, AR(1) catch fields. Redbanded: replaced by the
  # colour-only naive vs joint plot below (adds all-vessel sampling); bocaccio
  # (naive fits only so far) gets its naive panel
  if (!species %in% c("redbanded rockfish", "bocaccio")) {
    one_panel_plot(c("Old data (5-yr fleet), AR(1) catch fields",
      "Current data, target-vessel sampling, AR(1) catch fields"),
      "index-target-fields-vs-expected-absolute-1panel")
  } else {
    # Naive vs joint (b on expected catch), target- vs all-vessel sampling; AR(1)
    # catch fields. One colour per line, no alpha
    # Labelled by fleet: vessels with >= 5 positive trips in >= 1, 3 or 5 years
    # (min_yrs_with_trips; the 5-yr fleet is the old data). The naive 5-yr fit
    # on the old data equals its refit on the current code
    # (naive_tv_ar1_min5yr), so only one is drawn (other species: the refit)
    min5_notv <- if (species == "redbanded rockfish") "naive_tv_ar1-delta_gengamma-fp80" else
      "naive_tv_ar1_min5yr-delta_gengamma-fp80"
    sampling_lines <- tribble(
      ~line, ~key, ~col,
      "Min 1 yr fleet, naive", "naive_tv_ar1_target-delta_gengamma-fp80", "#4D4D4D",
      "Min 1 yr fleet, naive, no target_vessel", "naive_tv_ar1_current-delta_gengamma-fp80", "#999999",
      "Min 3 yr fleet, naive", "naive_tv_ar1_target_min3yr-delta_gengamma-fp80", "#009E73",
      "Min 3 yr fleet, naive, no target_vessel", "naive_tv_ar1_min3yr-delta_gengamma-fp80", "#7FD1B9",
      "Min 5 yr fleet, naive", "naive_tv_ar1_target_min5yr-delta_gengamma-fp80", "#B07800",
      "Min 5 yr fleet, naive, no target_vessel", min5_notv, "#E69F00",
      "Min 1 yr fleet, joint, target-vessel sampling", "constant_iidbase_fp80frame_target_xi_expected_tv_ar1_dgengamma", "#0072B2",
      "Min 3 yr fleet, joint, target-vessel sampling", "constant_iidbase_fp80frame_target_xi_expected_tv_ar1_dgengamma-min3yr-fp80", "#56B4E9",
      "Min 5 yr fleet, joint, target-vessel sampling", "constant_iidbase_fp80frame_target_xi_expected_tv_ar1_dgengamma-min5yr-fp80", "#882255",
      "Min 1 yr fleet, joint, all-vessel sampling", "constant_iidbase_fp80frame_allves_xi_expected_tv_ar1_dgengamma", "#D55E00"
    )
    d <- inner_join(select(sampling_lines, -col), select(indices, key, year, est, rel_est), by = "key") |>
      mutate(line = factor(line, sampling_lines$line), naive = grepl("naive", line))
    # Top: naive indices only, solid. Bottom: the same naive indices dotted, with
    # the joint fits solid
    panels <- c("Naive", "Naive (dotted) and joint (solid)")
    d2 <- bind_rows(
      mutate(filter(d, naive), panel = panels[1], lty = "solid"),
      mutate(d, panel = panels[2], lty = ifelse(naive, "dotted", "solid"))
    ) |> mutate(panel = factor(panel, panels))
    # No joint fits yet (bocaccio): the naive panel only
    naive_only <- all(d$naive)
    if (naive_only) d2 <- filter(d2, panel == panels[1])
    # Unscaled (est) and scaled (rel_est) versions
    for (scaled in c(FALSE, TRUE)) {
      v <- if (scaled) "rel_est" else "est"
      ggplot(d2, aes(year, .data[[v]], colour = line, linetype = lty)) +
        geom_line(linewidth = 0.9) +
        facet_wrap(~panel, ncol = 1) +
        scale_y_continuous(labels = if (scaled) waiver() else
          scales::label_number(scale_cut = scales::cut_short_scale())) +
        scale_colour_manual(values = setNames(sampling_lines$col, sampling_lines$line)) +
        scale_linetype_identity() +
        labs(x = "Year", y = if (scaled) "Relative index (geometric mean = 1)" else
          "Index (predicted kg/h x cell area in km², summed over cells)",
          colour = NULL,
          title = paste0(tools::toTitleCase(species),
            if (naive_only) ": naive indices by fleet definition" else
              ": naive vs joint, target- vs all-vessel sampling",
            if (scaled) "" else " (unscaled)"),
          caption = paste0("AR(1) catch fields.",
            if (!naive_only) " Joint: constant b on expected catch, IID sampling baseline.", " ",
            "Index over cells >= 80% inside the 2012 trawl footprint (current frame).\n",
            "Fleet: vessels with >= 5 positive trips in >= 1, 3 or 5 years",
            if (species == "redbanded rockfish") " (5 yr = old data)", ". Catch models use all fleet tows and a target_vessel\ncatchability effect (index at target_vessel = 0) unless marked. ",
            "Delta-gengamma, rw0 catch intercept; ", cell_km, " km grid.")) +
        guides(colour = guide_legend(ncol = 2)) +
        theme_light() +
        theme(legend.position = "top", legend.justification = "left",
          legend.key.width = unit(2.5, "lines"),
          strip.text = element_text(colour = "grey15", hjust = 0),
          strip.background = element_rect(fill = "grey92"))
      ggsave(paste0(fig_dir, "/index-sampling-comparison-",
        if (scaled) "" else "absolute-", if (naive_only) "1panel-" else "2panel-",
        cell_km, "km.png"), width = 10, height = if (naive_only) 6.5 else 10)
    }
  }
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
