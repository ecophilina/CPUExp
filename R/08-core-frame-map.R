# Map of the core footprint prediction frame used for the core indices
# (index_grids$fp80 in R/06-ar1-indices.R): frame cells at least 80% inside the
# official 2012 trawl footprint, from the current prepped data. Shows cells
# left out of the core (< 80% inside) and the footprint outline.
#
# Usage (from the project root):
#   CPUE_SPECIES="pacific ocean perch" Rscript R/08-core-frame-map.R

library(dplyr)
library(ggplot2)

source("R/species.R")

prepped <- read_prepped()
s <- prepped$settings
cell_km <- s$cell_km

cells <- prepped$frame |>
  distinct(X, Y, footprint_area) |>
  mutate(core = footprint_area >= 0.8)
n_core <- sum(cells$core)

footprint_poly <- sf::st_read(s$footprint_file, quiet = TRUE) |>
  sf::st_transform(s$utm_crs) |>
  sf::st_geometry() |>
  sf::st_union()
footprint_poly <- sf::st_set_crs(
  footprint_poly / 1000 - c(prepped$xy_offset[["X"]], prepped$xy_offset[["Y"]]), NA)

ggplot() +
  geom_tile(aes(X, Y, fill = ifelse(core, "Core (>= 80% inside footprint)",
    "Frame only (< 80% inside)")), data = cells, width = cell_km, height = cell_km) +
  geom_sf(data = footprint_poly, fill = NA, colour = "grey15", linewidth = 0.3) +
  geom_point(aes(X, Y), data = s$ports, shape = 17, size = 2) +
  geom_text(aes(X, Y, label = port), data = s$ports, hjust = -0.15, size = 3.2) +
  scale_fill_manual(values = c("Core (>= 80% inside footprint)" = "#0072B2",
    "Frame only (< 80% inside)" = "#BBBBBB")) +
  # room for the port labels
  coord_sf(xlim = range(c(cells$X, s$ports$X)) + c(-cell_km, 90)) +
  labs(x = "X (km)", y = "Y (km)", fill = NULL,
    title = paste0(tools::toTitleCase(species), ": core footprint prediction frame (",
      cell_km, " km cells)"),
    subtitle = paste0(n_core, " core cells of ", nrow(cells),
      " frame cells. Outline: official 2012 trawl footprint")) +
  theme_light() +
  theme(legend.position = "bottom")
ggsave(paste0(fig_dir, "/frame-core-footprint-", cell_km, "km.png"), width = 7, height = 8.5)
