# Depth mismatch between where the fleet fishes and where catch rates peak
#
# The joint fits with b on expected catch put the catch depth peak (~280 m)
# deeper than the sampling depth peak (~190 m). Cell depth in the frame is the
# mean log depth of the tows in the cell (R/data-prep.R), while the catch model
# uses each tow's depth, so part of the mismatch could come from aggregating to
# cells. Compares, over the core footprint frame (cells >= 80% inside the 2012
# trawl footprint):
#   1. Within-cell spread of tow depths, by cell depth
#   2. Effort and catch rate by depth at tow level vs at cell level
#   3. Within cells: catch rate vs the tow's depth relative to the cell mean
#      (removes between-cell differences in abundance)
# Raw catch rates only (no month, vessel or year adjustment).
#
# Usage (from the project root):
#   CPUE_SPECIES="pacific ocean perch" Rscript R/depth-mismatch.R

library(dplyr)
library(ggplot2)

source("R/species.R")
prepped <- read_prepped()
cell_km <- prepped$settings$cell_km

core <- prepped$frame |>
  filter(footprint_area >= 0.8) |>
  distinct(X, Y, cell_depth = exp(log_depth))
tows <- prepped$d |>
  inner_join(core, by = c(cell_X = "X", cell_Y = "Y")) |>
  mutate(cpue = catch_kg / hours, pos = catch_kg > 0)
frame <- prepped$frame |>
  filter(footprint_area >= 0.8) |>
  mutate(cell_depth = exp(log_depth))
message(nrow(tows), " tows in ", nrow(core), " core cells")

depth_bins <- seq(50, 600, by = 25)
bin_mid <- function(x) {
  b <- cut(x, depth_bins, include.lowest = TRUE)
  (depth_bins[-length(depth_bins)] + 12.5)[as.integer(b)]
}

# 1. Within-cell spread --------------------------------------------------------

spread <- tows |>
  group_by(cell_X, cell_Y, cell_depth) |>
  summarise(n_tows = n(), sd_log_depth = sd(log(depth)),
    q10 = quantile(depth, 0.1), q90 = quantile(depth, 0.9), .groups = "drop") |>
  mutate(range_80 = q90 - q10, depth_bin = bin_mid(cell_depth))
spread_summary <- spread |>
  group_by(depth_bin) |>
  summarise(cells = n(), median_tows = median(n_tows),
    median_range_80_m = median(range_80), median_sd_log = median(sd_log_depth),
    .groups = "drop")
cat("\nWithin-cell spread of tow depth (10-90% range), by cell depth:\n")
print(as.data.frame(spread_summary), digits = 3)

# 2. Effort and catch rate by depth, tow level vs cell level -------------------

tow_level <- tows |>
  mutate(depth_bin = bin_mid(depth)) |>
  group_by(depth_bin) |>
  summarise(effort = sum(hours), n = n(), encounter = mean(pos),
    cpue = sum(catch_kg) / sum(hours), .groups = "drop") |>
  mutate(effort_share = effort / sum(effort), level = "Tow depth")
cell_level <- tows |>
  mutate(depth_bin = bin_mid(cell_depth)) |>
  group_by(depth_bin) |>
  summarise(effort = sum(hours), n = n(), encounter = mean(pos),
    cpue = sum(catch_kg) / sum(hours), .groups = "drop") |>
  mutate(effort_share = effort / sum(effort), level = "Cell depth")
# Share of available cell-years fished, by cell depth (what the sampling
# model sees)
p_sampled <- frame |>
  mutate(depth_bin = bin_mid(cell_depth)) |>
  group_by(depth_bin) |>
  summarise(cell_years = n(), p_sampled = mean(sampled), .groups = "drop")

by_depth <- bind_rows(tow_level, cell_level) |>
  filter(n >= 100) |>
  mutate(level = factor(level, c("Tow depth", "Cell depth")))
peaks <- by_depth |>
  group_by(level) |>
  summarise(
    effort_peak = depth_bin[which.max(effort_share)],
    effort_median = NA_real_,
    cpue_peak = depth_bin[which.max(cpue)],
    encounter_peak = depth_bin[which.max(encounter)], .groups = "drop")
# Effort- and catch-weighted median depths at each level
wmed <- function(x, w) { o <- order(x); x[o][which(cumsum(w[o]) >= sum(w) / 2)[1]] }
peaks$effort_median <- c(wmed(tows$depth, tows$hours), wmed(tows$cell_depth, tows$hours))
peaks$catch_median <- c(wmed(tows$depth, tows$catch_kg), wmed(tows$cell_depth, tows$catch_kg))
cat("\nPeak and median depths (m; 25 m bins with >= 100 tows):\n")
print(as.data.frame(peaks))
cat("\nShare of available cell-years fished, by cell depth:\n")
print(as.data.frame(filter(p_sampled, cell_years >= 50)), digits = 3)

# 3. Within cells: catch rate vs relative tow depth ----------------------------

# Tow depth relative to its cell's mean (log ratio), binned; catch rate
# relative to the cell's overall catch rate. Cells with >= 20 tows and some
# POP catch.
within <- tows |>
  group_by(cell_X, cell_Y) |>
  filter(n() >= 20, sum(catch_kg) > 0) |>
  mutate(rel_depth = log(depth) - mean(log(depth)),
    cell_cpue = sum(catch_kg) / sum(hours),
    cell_enc = mean(pos)) |>
  ungroup() |>
  mutate(rel_bin = cut(rel_depth, c(-Inf, -0.3, -0.2, -0.1, -0.03, 0.03, 0.1, 0.2, 0.3, Inf)))
within_summary <- within |>
  group_by(rel_bin) |>
  summarise(tows = n(), effort_share = sum(hours),
    rel_cpue = (sum(catch_kg) / sum(hours)) / weighted.mean(cell_cpue, hours),
    rel_encounter = mean(pos) / mean(cell_enc), .groups = "drop") |>
  mutate(effort_share = effort_share / sum(effort_share))
cat("\nWithin cells: catch rate by tow depth relative to the cell mean",
  "(log ratio; 0.1 ~ 10% deeper):\n")
print(as.data.frame(within_summary), digits = 3)
# Same, split by cell depth: shelf-edge vs slope cells
within_by_cell <- within |>
  mutate(cell_zone = cut(cell_depth, c(0, 150, 250, 350, Inf),
    labels = c("< 150 m", "150-250 m", "250-350 m", "> 350 m")),
    deeper = rel_depth > 0) |>
  group_by(cell_zone, deeper) |>
  summarise(tows = n(), cpue = sum(catch_kg) / sum(hours), .groups = "drop") |>
  tidyr::pivot_wider(names_from = deeper, values_from = c(tows, cpue)) |>
  mutate(deeper_vs_shallower_cpue = cpue_TRUE / cpue_FALSE)
cat("\nWithin cells, deeper vs shallower half of the cell's tows:\n")
print(as.data.frame(within_by_cell), digits = 3)

# Figure -----------------------------------------------------------------------

level_cols <- c("Tow depth" = "#0072B2", "Cell depth" = "#D55E00")
fig_data <- by_depth |>
  mutate(effort_share = effort_share / max(effort_share), cpue = cpue / max(cpue)) |>
  tidyr::pivot_longer(c(effort_share, cpue), names_to = "measure") |>
  mutate(measure = factor(measure, c("effort_share", "cpue"),
    c("Fishing effort (hours, scaled to max = 1)", "POP catch rate (kg/h, scaled to max = 1)")))
p1 <- ggplot(fig_data, aes(depth_bin, value, colour = level)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
  facet_wrap(~measure, ncol = 1) +
  scale_colour_manual(values = level_cols) +
  labs(x = "Depth (m, 25 m bins)", y = NULL, colour = "Depth assigned by",
    title = paste0(tools::toTitleCase(species), ": effort and catch rate by depth"),
    caption = paste0("Core footprint frame (cells >= 80% inside the 2012 footprint), ",
      cell_km, " km cells. Raw catch rates; bins with >= 100 tows.")) +
  theme_light() + theme(legend.position = "top", legend.justification = "left")
ggsave(file.path(fig_dir, "depth-mismatch-effort-cpue.png"), p1, width = 8, height = 7)

p2 <- ggplot(spread, aes(cell_depth, range_80)) +
  geom_point(alpha = 0.4, size = 1.2, colour = "#4D4D4D") +
  geom_smooth(se = FALSE, colour = "#0072B2", method = "loess", formula = y ~ x) +
  labs(x = "Cell depth (m; mean log depth of its tows)",
    y = "10-90% range of tow depths in the cell (m)",
    title = paste0(tools::toTitleCase(species), ": spread of tow depths within cells"),
    caption = paste0("Core footprint frame, ", cell_km, " km cells.")) +
  theme_light()
ggsave(file.path(fig_dir, "depth-mismatch-within-cell.png"), p2, width = 7, height = 5)

write.csv(bind_rows(by_depth), file.path(out_dir, "depth-mismatch-by-depth.csv"), row.names = FALSE)
write.csv(within_summary, file.path(out_dir, "depth-mismatch-within-cell.csv"), row.names = FALSE)
