# Diagnostics for the preferential-sampling correction
#
# The joint models (b ~ -0.6 to -0.8) turn a nearly flat naive index into a
# roughly doubling one while the fleet's footprint contracted (~60% -> ~23% of
# frame cells fished). These checks ask where the change comes from:
#   1. Fished vs unfished cells: predicted density on the frame, split by
#      whether each cell was fished that year (naive vs joint, both families).
#      How much of the joint increase sits in cells with no data?
#   2. Contraction: which cells the fleet dropped and when, and whether dropped
#      cells had higher *observed* redbanded catch rates while still fished
#      (model-free evidence for or against avoidance).
#
# Predictions are plug-in (no bias correction): fine for shares and trends.
# Usage (from the project root): Rscript R/index-diagnostics.R

library(dplyr)
library(ggplot2)
library(sdmTMB)

mesh_cutoff <- 15
early_years <- 1996:2005
late_years <- 2016:2025
min_years_period <- 3 # "regularly fished" in a period = fished in >= this many years

prepped <- readRDS("data/cpue-index-clean.rds") # from R/data-prep.R
d <- prepped$d
frame <- prepped$frame
cell_km <- prepped$settings$cell_km

pred_dir <- file.path("output", "diagnostics")
dir.create(pred_dir, showWarnings = FALSE)

fits <- tribble(
  ~label, ~fit, ~family, ~model,
  "naive, delta-Gamma", "naive-delta_gamma", "delta-Gamma", "naive",
  "joint, delta-Gamma", "constant_xi_expected_dgamma_newton", "delta-Gamma", "joint",
  "naive, delta-lognormal", "naive-delta_lognormal", "delta-lognormal", "naive",
  "joint, delta-lognormal", "constant_xi_expected_dlnorm_newton", "delta-lognormal", "joint"
)

# 1. Fished vs unfished cells ---------------------------------------------------

# Frame rows, in frame order, with the year's sampled indicator
newdata <- frame |> select(X, Y, year, log_depth, month)

predict_frame <- function(fit_name) {
  f <- file.path(pred_dir, paste0("pred-", fit_name, "-cutoff", mesh_cutoff, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  message("Predicting: ", fit_name)
  fit_dir <- if (grepl("^naive", fit_name)) "output" else "models"
  fit <- readRDS(file.path(fit_dir,
    paste0("fit-", fit_name, "-cutoff", mesh_cutoff, ".rds")))
  p <- predict(fit, newdata = newdata, type = "response", re_form_iid = NA)
  # Delta models: combined expected catch (encounter prob x positive mean)
  est <- if ("est" %in% names(p)) p$est else p$est1 * p$est2
  out <- frame |>
    select(X, Y, year, sampled, region) |>
    mutate(est = est, fit = fit_name)
  saveRDS(out, f)
  out
}

preds <- fits |>
  rowwise() |>
  mutate(pred = list(predict_frame(fit))) |>
  ungroup() |>
  select(-fit) |>
  tidyr::unnest(pred) |>
  mutate(status = if_else(sampled == 1, "fished that year", "not fished"))

# Mean predicted density (kg/h at reference conditions) by year and status
by_status <- preds |>
  group_by(label, family, model, year, status) |>
  summarise(mean_est = mean(est), n_cells = n(), total = sum(est),
    .groups = "drop")

ggplot(by_status, aes(year, mean_est, colour = model, linetype = status)) +
  geom_line() +
  facet_wrap(~family, ncol = 1, scales = "free_y") +
  scale_y_log10() +
  labs(x = "Year", y = "Mean predicted density per cell (log scale)",
    colour = NULL, linetype = NULL) +
  theme_light()
ggsave("figs/diag-fished-unfished-density.png", width = 9, height = 7)

# Share of each year's index total from unfished cells
unfished_share <- by_status |>
  group_by(label, family, model, year) |>
  summarise(share_unfished = sum(total[status == "not fished"]) / sum(total),
    prop_cells_unfished = sum(n_cells[status == "not fished"]) / sum(n_cells),
    .groups = "drop")

ggplot(unfished_share, aes(year, share_unfished, colour = model)) +
  geom_line() +
  geom_line(aes(y = prop_cells_unfished), colour = "grey50", linetype = 2) +
  facet_wrap(~family, ncol = 1) +
  labs(x = "Year", y = "Share of index from unfished cells\n(dashed: share of cells unfished)",
    colour = NULL) +
  theme_light()
ggsave("figs/diag-unfished-share.png", width = 8, height = 6)

# Relative indices: all cells vs fished cells only (mean density, since the
# fished area changes each year)
rel_index <- by_status |>
  group_by(label, family, model, year) |>
  summarise(
    all_cells = sum(total) / sum(n_cells),
    fished_only = mean_est[status == "fished that year"],
    unfished_only = mean_est[status == "not fished"],
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(c(all_cells, fished_only, unfished_only),
    names_to = "cells", values_to = "density") |>
  group_by(label, cells) |>
  mutate(rel = density / exp(mean(log(density)))) |>
  ungroup()

ggplot(rel_index, aes(year, rel, colour = model)) +
  geom_line() +
  facet_grid(family ~ cells) +
  labs(x = "Year", y = "Relative mean density (geometric mean = 1)",
    colour = NULL) +
  theme_light()
ggsave("figs/diag-index-by-cells.png", width = 11, height = 6)

trend_by_cells <- rel_index |>
  group_by(family, model, cells) |>
  summarise(
    slope = coef(lm(log(density) ~ year))[[2]],
    change_30y = exp(slope * (max(year) - min(year))) - 1,
    .groups = "drop"
  ) |>
  arrange(family, cells, model)
print(trend_by_cells, n = Inf)

# 2. Contraction ----------------------------------------------------------------

cell_use <- frame |>
  group_by(X, Y, region) |>
  summarise(
    n_early = sum(sampled[year %in% early_years]),
    n_late = sum(sampled[year %in% late_years]),
    first_year = min(year[sampled == 1]),
    last_year = max(year[sampled == 1]),
    log_depth = first(log_depth),
    log_dist_port = first(log_dist_port),
    .groups = "drop"
  ) |>
  mutate(class = case_when(
    n_early >= min_years_period & n_late >= min_years_period ~ "retained",
    n_early >= min_years_period ~ "dropped",
    n_late >= min_years_period ~ "added",
    TRUE ~ "occasional"
  ))

print(count(cell_use, region, class) |>
  tidyr::pivot_wider(names_from = class, values_from = n, values_fill = 0))

ggplot(cell_use, aes(X, Y, fill = class)) +
  geom_tile(width = cell_km, height = cell_km) +
  scale_fill_manual(values = c(retained = "grey60", dropped = "#d95f02",
    added = "#1b9e77", occasional = "#e6e6e6")) +
  coord_fixed() +
  labs(x = "X (km)", y = "Y (km)", fill = NULL,
    title = paste0("Cells fished in >= ", min_years_period, " years of ",
      min(early_years), "-", max(early_years), " vs ",
      min(late_years), "-", max(late_years))) +
  theme_light()
ggsave("figs/diag-contraction-map.png", width = 8, height = 9)

# Timing: cells fished per year by region
frame |>
  group_by(year, region) |>
  summarise(cells_fished = sum(sampled), .groups = "drop") |>
  ggplot(aes(year, cells_fished, colour = region)) +
  geom_line() +
  labs(x = "Year", y = "Frame cells fished", colour = "Region") +
  theme_light()
ggsave("figs/diag-cells-fished-by-region.png", width = 8, height = 5)

# Observed catch rates in dropped vs retained cells during the early period,
# when both were fished: if the fleet avoided redbanded, dropped cells should
# have had higher catch rates.
obs_early <- d |>
  filter(year %in% early_years) |>
  inner_join(select(cell_use, X, Y, class),
    by = c("cell_X" = "X", "cell_Y" = "Y")) |>
  group_by(class) |>
  summarise(
    tows = n(),
    cells = n_distinct(paste(cell_X, cell_Y)),
    prop_positive = mean(catch_kg > 0),
    mean_cpue = mean(catch_kg / hours),
    median_pos_cpue = median((catch_kg / hours)[catch_kg > 0]),
    .groups = "drop"
  )
print(obs_early)

# Same by year (dropped vs retained), to see whether a gap opens before cells
# are dropped
obs_by_year <- d |>
  inner_join(select(cell_use, X, Y, class),
    by = c("cell_X" = "X", "cell_Y" = "Y")) |>
  filter(class %in% c("retained", "dropped")) |>
  group_by(year, class) |>
  summarise(prop_positive = mean(catch_kg > 0),
    mean_cpue = mean(catch_kg / hours), tows = n(), .groups = "drop")

obs_by_year |>
  tidyr::pivot_longer(c(prop_positive, mean_cpue)) |>
  ggplot(aes(year, value, colour = class)) +
  geom_line() +
  facet_wrap(~name, scales = "free_y", ncol = 1) +
  labs(x = "Year", y = NULL, colour = NULL,
    title = "Observed redbanded catch in dropped vs retained cells") +
  theme_light()
ggsave("figs/diag-dropped-vs-retained-observed.png", width = 8, height = 6)

write.csv(trend_by_cells, file.path(pred_dir, "trend-by-cells.csv"), row.names = FALSE)
write.csv(obs_early, file.path(pred_dir, "observed-early-by-class.csv"), row.names = FALSE)
saveRDS(list(by_status = by_status, unfished_share = unfished_share,
  rel_index = rel_index, cell_use = cell_use, obs_by_year = obs_by_year),
  file.path(pred_dir, "index-diagnostics.rds"))
