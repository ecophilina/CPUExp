# Catch random fields in fished vs unfished cell-years, by year
#
# A joint fit's index moves away from the naive index by roughly b times the
# difference in the catch random fields between fished and unfished cell-years
# of the sampling frame. With constant b, the joint and naive AR(1) indices
# agree until ~2013 and diverge after, which suggests this gap opens up late.
# This takes the naive AR(1) fits' fields (catch data only), predicts them over
# the core footprint frame, and compares fished and unfished cell-years each
# year:
#   - old data: naive_tv_ar1 with the old frame's sampled response
#   - current data: naive_tv_ar1_target with sampled_target
# "fields" is the random-field part of log expected catch: log expected catch
# minus the same with omega and epsilon set to zero (both delta components).
#
# Usage (from the project root; loads two naive fits, not the joint ones):
#   CPUE_SPECIES="pacific ocean perch" Rscript R/field-gap.R

library(dplyr)
library(ggplot2)
library(sdmTMB)

source("R/species.R")
mesh_cutoff <- 15

field_gap <- function(fit_name, data_file, response, label) {
  prepped <- read_prepped(file.path(data_dir, data_file))
  frame <- filter(prepped$frame, footprint_area >= 0.8)
  frame$vessel <- factor(levels(prepped$d$vessel)[1], levels(prepped$d$vessel))
  fit <- readRDS(find_fit(fit_name, mesh_cutoff))
  p <- predict(fit, newdata = frame, re_form_iid = NA)
  log_expected <- function(eta1, eta2) plogis(eta1, log.p = TRUE) + eta2
  p |>
    mutate(
      fished = .data[[response]] == 1,
      fields = log_expected(est1, est2) - log_expected(est_non_rf1, est_non_rf2),
      omega = omega_s1 + omega_s2,
      epsilon = epsilon_st1 + epsilon_st2
    ) |>
    group_by(year, fished) |>
    summarise(n = n(), fields = mean(fields), omega = mean(omega),
      epsilon = mean(epsilon), .groups = "drop") |>
    tidyr::pivot_wider(names_from = fished, values_from = c(n, fields, omega, epsilon)) |>
    transmute(year, data = label, n_fished = n_TRUE, n_unfished = n_FALSE,
      gap_fields = fields_TRUE - fields_FALSE,
      gap_omega = omega_TRUE - omega_FALSE,
      gap_epsilon = epsilon_TRUE - epsilon_FALSE)
}

gaps <- bind_rows(
  field_gap("naive_tv_ar1-delta_gengamma", "cpue-index-clean-5yr-fleet.rds",
    "sampled", "Old data (5-yr fleet)"),
  field_gap("naive_tv_ar1_target-delta_gengamma", "cpue-index-clean.rds",
    "sampled_target", "Current data, target vessels")
)
print(as.data.frame(gaps), digits = 3)
write.csv(gaps, file.path(out_dir, "field-gap-by-year.csv"), row.names = FALSE)

gap_long <- gaps |>
  tidyr::pivot_longer(c(gap_fields, gap_omega, gap_epsilon), names_to = "part") |>
  mutate(part = factor(part, c("gap_fields", "gap_omega", "gap_epsilon"),
    c("All fields (log expected catch)", "Spatial (omega, both components)",
      "Spatiotemporal (epsilon, both components)")))
part_cols <- c("#4D4D4D", "#0072B2", "#D55E00")
ggplot(gap_long, aes(year, value, colour = part)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~data, ncol = 1) +
  scale_colour_manual(values = setNames(part_cols, levels(gap_long$part))) +
  labs(x = "Year", y = "Mean in fished minus unfished cell-years (log scale)",
    colour = NULL,
    title = paste0(tools::toTitleCase(species),
      ": catch random fields where the fleet fished vs did not"),
    caption = paste0("Naive AR(1) fits (catch data only), predicted over the core ",
      "footprint frame (cells >= 80% inside the 2012 footprint).")) +
  guides(colour = guide_legend(ncol = 1)) +
  theme_light() +
  theme(legend.position = "top", legend.justification = "left")
ggsave(file.path(fig_dir, "field-gap-by-year.png"), width = 8, height = 7)
