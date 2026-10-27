# Data prep for the commercial trawl CPUE index
#
# Builds the catch data (one row per tow), runs exploratory plots, builds the
# preferential-sampling frame and the index grid, and saves them to
# <species>/data/cpue-index-clean.rds for R/cpue-index.R and later scripts.
#
# DRAFT where (cell size, fleet definition, depth trimming, ports, regions) are
# placeholders to be checked against the data and the assessment's needs.

library(dplyr)
library(ggplot2)
library(sdmTMB) # add_utm_columns()

source("R/species.R") # species (CPUE_SPECIES) and its folders

out_file <- file.path(data_dir, "cpue-index-clean.rds")

# Settings ---------------------------------------------------------------------

final_year <- 2025
utm_crs <- 32609 # UTM zone 9N
cell_km <- 6 # sampling-frame cell size; b and gamma depend on this choice
min_yrs_cell <- 3 # eligible cell = fished by the fleet in >= this many years
min_positive_trips <- 5 # fleet definition (as in gfplot::tidy_cpue_index())
min_yrs_with_trips <- 5
max_hours <- 10 # drop implausible tow durations
depth_quantiles <- c(0.01, 0.99) # trim depth range using positive tows
majors <- formatC(3:9, width = 2, flag = "0") # 3C, 3D, 5A, 5B, 5C, 5D, 5E (GFFOS major stat area codes)

# Official 2012 groundfish trawl footprint polygon (WGS 84)
footprint_file <- "data/TrawlFootprint/Trawl_Footprint_Final_May_14_2012_LL.shp"

# Approximate trawl offload ports (lon, lat); replace with the ports that
# matter for this fleet, or with a steaming-distance layer.
ports <- data.frame(
  port = c("Prince Rupert", "Port Hardy", "Ucluelet", "Steveston"),
  longitude = c(-130.32, -127.49, -125.55, -123.18),
  latitude = c(54.31, 50.72, 48.94, 49.12)
)

# Data -------------------------------------------------------------------------

# gfdata::get_cpue_index() returns one row per fishing event x species caught
# (all species), so the target species is pulled out and zeros filled here.
raw <- readRDS("data/cpue-index-data.rds")

events <- raw |> mutate(
    fishing_event_id_unique = paste0(
    vessel_registration_number, "-",
    trip_id, "-", fishing_event_id)
  ) |>
  filter(
    toupper(fishery_sector) == "GROUNDFISH TRAWL", # check sector labels
    major_stat_area_code %in% majors
  ) |>
  group_by(fishing_event_id_unique) |> ### need to check all are unique
  summarise(
    fishery_sector = first(fishery_sector),
    trip_id = first(trip_id),
    vessel = first(vessel_registration_number),
    best_date = first(best_date),
    fe_start_date = first(fe_start_date),
    fe_end_date = first(fe_end_date),
    depth = first(best_depth),
    major = first(major_stat_area_code),
    latitude = first(latitude),
    longitude = first(longitude),
    # %in% not ==: some species codes have NA common names, and NA indices
    # would turn the whole event's sum into NA
    catch_kg = sum(catch_kg[toupper(species_common_name) %in% toupper(species)]),
    .groups = "drop"
  ) |>
  mutate(
    year = as.integer(format(best_date, "%Y")),
    month = as.integer(format(best_date, "%m")),
    hours = as.numeric(difftime(fe_end_date, fe_start_date, units = "hours")),
    longitude = -abs(longitude) # GFFOS longitudes are sometimes positive
  ) |>
  filter(
    !is.na(latitude), !is.na(longitude), !is.na(depth), depth > 0,
    !is.na(hours), hours > 0, hours <= max_hours,
    year <= final_year
  )

# Fleet: vessels with >= min_positive_trips positive trips in
# >= min_yrs_with_trips years.
fleet_vessels <- events |>
  filter(catch_kg > 0) |>
  distinct(vessel, year, trip_id) |>
  count(vessel, year, name = "n_pos_trips") |>
  filter(n_pos_trips >= min_positive_trips) |>
  count(vessel, name = "n_years") |>
  filter(n_years >= min_yrs_with_trips) |>
  pull(vessel)


depth_range <- events |>
  filter(catch_kg > 0) |>
  pull(depth) |>
  quantile(depth_quantiles)

d <- events |>
  filter(
    vessel %in% fleet_vessels,
    depth >= depth_range[[1]], depth <= depth_range[[2]]
  ) |>
  add_utm_columns(utm_crs = utm_crs, units = "km")

# Shift coordinates to near zero (UTM northings are ~5500 km, which triggers
# sdmTMB's large-coordinate warning). Offsets are rounded to 100 km so the
# 4 km frame cells and 10 km exploratory bins stay aligned with the UTM grid.
# Add xy_offset back to map X/Y to true UTM km. Anything else placed in this
# space (ports, prediction grids) must be shifted by the same amounts.
xy_offset <- c(X = round(mean(d$X), -2), Y = round(mean(d$Y), -2))
d <- d |>
  mutate(
    X = X - xy_offset[["X"]],
    Y = Y - xy_offset[["Y"]],
    log_depth = log(depth),
    log_hours = log(hours),
    cell_X = floor(X / cell_km) * cell_km + cell_km / 2,
    cell_Y = floor(Y / cell_km) * cell_km + cell_km / 2
  )

# unique(d$fishery_sector)

# Exploratory plots -----------------------------------------------------

# Raw (unstandardized) catch rates by year in space. Binned more coarsely than
# the sampling frame so coastwide facets stay readable.
explore_km <- 10

d_explore <- d |>
  mutate(
    cpue = catch_kg / hours,
    bin_X = floor(X / explore_km) * explore_km + explore_km / 2,
    bin_Y = floor(Y / explore_km) * explore_km + explore_km / 2
  )

bin_year <- d_explore |>
  group_by(year, bin_X, bin_Y) |>
  summarise(
    n_tows = n(),
    prop_pos = mean(catch_kg > 0),
    mean_cpue = mean(cpue),
    .groups = "drop"
  )

# Mean CPUE (kg/h, zeros included) by bin and year; grey = fished, none caught
ggplot(bin_year, aes(bin_X, bin_Y, fill = ifelse(mean_cpue > 0, mean_cpue, NA))) +
  geom_tile(width = explore_km, height = explore_km) +
  facet_wrap(~year, ncol = 6) +
  scale_fill_viridis_c(name = "Mean CPUE\n(kg/h)", trans = "log10",
    na.value = "grey80") +
  coord_fixed() +
  labs(x = "X (km)", y = "Y (km)") +
  theme_light()
ggsave(file.path(fig_dir, "explore-cpue-map.png"), width = 12, height = 14)

# Encounter rate by bin and year
ggplot(bin_year, aes(bin_X, bin_Y, fill = prop_pos)) +
  geom_tile(width = explore_km, height = explore_km) +
  facet_wrap(~year, ncol = 6) +
  scale_fill_viridis_c(name = "Proportion\npositive", limits = c(0, 1)) +
  coord_fixed() +
  labs(x = "X (km)", y = "Y (km)") +
  theme_light()
ggsave(file.path(fig_dir, "explore-encounter-map.png"), width = 12, height = 14)

# Spread of bin-level CPUE within each year (positive bins only)
bin_year |>
  filter(mean_cpue > 0) |>
  ggplot(aes(factor(year), mean_cpue)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_y_log10() +
  labs(x = "Year", y = "Mean CPUE per bin (kg/h, positive bins)") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(file.path(fig_dir, "explore-cpue-spread.png"), width = 10, height = 5)

# Annual summaries: effort, footprint, encounter rate, and raw catch rates
d_explore |>
  group_by(year) |>
  summarise(
    `Tows` = n(),
    `Vessels` = n_distinct(vessel),
    `Frame-size cells fished` = n_distinct(paste(cell_X, cell_Y)),
    `Proportion positive` = mean(catch_kg > 0),
    `Mean CPUE (kg/h)` = mean(cpue),
    `Mean positive CPUE (kg/h)` = mean(cpue[catch_kg > 0]),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(-year) |>
  ggplot(aes(year, value)) +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(~name, scales = "free_y", ncol = 2) +
  labs(x = "Year", y = NULL) +
  theme_light()
ggsave(file.path(fig_dir, "explore-annual-summary.png"), width = 9, height = 8)

# Centre of gravity of effort vs catch: if catch moves away from effort (or
# effort tracks catch), that hints at shifting preferential sampling.
d_explore |>
  group_by(year) |>
  summarise(
    effort_X = mean(X), effort_Y = mean(Y),
    catch_X = weighted.mean(X, catch_kg, na.rm = TRUE), catch_Y = weighted.mean(Y, catch_kg, na.rm = TRUE),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(-year, names_to = c("what", "axis"), names_sep = "_") |>
  ggplot(aes(year, value, colour = what)) +
  geom_line() +
  facet_wrap(~axis, scales = "free_y", ncol = 1) +
  labs(x = "Year", y = "Centre of gravity (km)", colour = NULL) +
  theme_light()
ggsave(file.path(fig_dir, "explore-centre-of-gravity.png"), width = 7, height = 6)


# Sampling frame ---------------------------------------------------------------

# One row per eligible cell x year. "Eligible" here = within the fleet's
# footprint (fished in >= min_yrs_cell years). Cells that could never be fished
# (land, permanent closures such as sponge-reef or RCA closures, out of depth
# range) must be omitted, not coded 0. ??Time-varying closures should be dropped
# from the years they were closed.??? I think this is wrong.

# If not done previously
# d <- filter(d, year <= final_year)

region_lu <- c(
  "03" = "WCVI", "04" = "WCVI", "05" = "QCS", "06" = "QCS",
  "07" = "HS", "08" = "HS", "09" = "WCHG"
)

cells <- d |>
  group_by(cell_X, cell_Y) |>
  summarise(
    n_years = n_distinct(year),
    # Cell depth from observed tows; swap in a bathymetry layer if available.
    log_depth = mean(log_depth),
    major = names(which.max(table(major))),
    .groups = "drop"
  ) |>
  filter(n_years >= min_yrs_cell) |>
  mutate(region = factor(region_lu[as.character(major)]))

ports <- add_utm_columns(ports, utm_crs = utm_crs, units = "km") |>
  mutate(X = X - xy_offset[["X"]], Y = Y - xy_offset[["Y"]])
port_dist <- sqrt(
  outer(cells$cell_X, ports$X, "-")^2 + outer(cells$cell_Y, ports$Y, "-")^2
)
# Straight-line distance ignores land; a proxy for access cost.
cells$log_dist_port <- log(pmax(apply(port_dist, 1, min), cell_km))

# Keep catch data inside the frame footprint so both likelihoods describe the
# same fishery.
d <- semi_join(d, cells, by = c("cell_X", "cell_Y"))

years <- sort(unique(d$year))
stopifnot(all(diff(years) == 1)) # else use extra_time for coefficient = "rw"

ref_month <- as.integer(names(which.max(table(d$month))))
d$month <- factor(d$month)
d$vessel <- factor(d$vessel)

sampled_cy <- d |>
  distinct(cell_X, cell_Y, year) |>
  mutate(sampled = 1L)

frame <- cells |>
  select(cell_X, cell_Y, log_depth, log_dist_port, region) |>
  cross_join(data.frame(year = years)) |>
  left_join(sampled_cy, by = c("cell_X", "cell_Y", "year")) |>
  mutate(
    sampled = coalesce(sampled, 0L),
    # Catchability covariates held at reference values on the frame:
    month = factor(ref_month, levels = levels(d$month))
  ) |>
  rename(X = cell_X, Y = cell_Y) |>
  as.data.frame()

# Data-driven proxy for the 2012 trawl footprint freeze: footprint = 1 for
# cells fished from footprint_year on, 0 for cells only fished before it;
# pre_footprint = 1 for years before footprint_year. footprint = 0 cells are
# never sampled from footprint_year on (by construction), so a sampling model
# using footprint with those cell-years in the frame has complete separation.
footprint_year <- 2013
footprint_cells <- frame |>
  group_by(X, Y) |>
  summarise(footprint = as.integer(any(sampled == 1 & year >= footprint_year)),
    .groups = "drop")
frame <- frame |>
  left_join(footprint_cells, by = c("X", "Y")) |>
  mutate(pre_footprint = as.integer(year < footprint_year))

# Frame cells and the footprint proxy: cells fished from footprint_year on vs
# cells only fished before it
footprint_cells |>
  mutate(footprint = factor(footprint, levels = c(1, 0),
    labels = c(paste0("Fished ", footprint_year, "+"),
      paste0("Only fished before ", footprint_year)))) |>
  ggplot(aes(X, Y)) +
  geom_tile(aes(fill = footprint), width = cell_km, height = cell_km) +
  geom_point(data = ports, shape = 17, size = 2) +
  geom_text(aes(label = port), data = ports, hjust = -0.15, size = 3) +
  scale_fill_manual(values = c("#2c7fb8", "#f03b20"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) + # port labels
  coord_fixed() +
  labs(x = "X (km)", y = "Y (km)",
    title = paste0("Sampling frame (", nrow(footprint_cells), " cells, ",
      cell_km, " km)"),
    subtitle = paste0(sum(footprint_cells$footprint), " in footprint, ",
      sum(footprint_cells$footprint == 0), " only fished before ",
      footprint_year)) +
  theme_light() +
  theme(legend.position = "bottom")
ggsave(paste0(fig_dir, "/frame-footprint-map-", cell_km, "km.png"), width = 7, height = 8)

# Official footprint: the 2012 trawl footprint polygon (WGS 84), moved into the
# frame's coordinates (UTM km, shifted by xy_offset). footprint_area = share of
# each frame cell's area inside the polygon (6 km cells often straddle its
# edge, so this is more robust than testing the cell centre); footprint_2012 =
# 1 if at least half the cell is inside. An alternative to the data-driven
# footprint above: it does not depend on where the fleet happened to fish.
footprint_poly <- sf::st_read(footprint_file, quiet = TRUE) |>
  sf::st_transform(utm_crs) |>
  sf::st_geometry() |>
  sf::st_union()
footprint_poly <- sf::st_set_crs(
  footprint_poly / 1000 - c(xy_offset[["X"]], xy_offset[["Y"]]), NA)
cell_polys <- sf::st_sfc(lapply(seq_len(nrow(footprint_cells)), function(i) {
  x <- footprint_cells$X[i] + c(-1, 1, 1, -1, -1) * cell_km / 2
  y <- footprint_cells$Y[i] + c(-1, -1, 1, 1, -1) * cell_km / 2
  sf::st_polygon(list(cbind(x, y)))
}))
footprint_cells$footprint_area <- vapply(
  sf::st_intersects(cell_polys, footprint_poly),
  function(hit) length(hit) > 0, logical(1))
inside <- which(footprint_cells$footprint_area)
footprint_cells$footprint_area <- 0
footprint_cells$footprint_area[inside] <- as.numeric(sf::st_area(
  sf::st_intersection(cell_polys[inside], footprint_poly))) / cell_km^2
footprint_cells <- footprint_cells |>
  mutate(footprint_2012 = as.integer(footprint_area >= 0.5))
frame <- frame |>
  left_join(select(footprint_cells, X, Y, footprint_area, footprint_2012),
    by = c("X", "Y"))

# How the two footprint definitions agree. Cells outside the official
# footprint that were still fished from footprint_year on suggest edge
# mismatch (cell size, positional error) or the polygon not being the whole
# story; if there are none, footprint_2012 has the same separation problem as
# footprint.
print(count(footprint_cells, footprint_2012, footprint))
cat("Frame cell-years outside the official footprint from", footprint_year,
  "on that were fished:",
  sum(frame$sampled[frame$footprint_2012 == 0 & frame$year >= footprint_year]),
  "of", sum(frame$footprint_2012 == 0 & frame$year >= footprint_year), "\n")

footprint_cells |>
  mutate(agreement = factor(paste(footprint_2012, footprint),
    levels = c("1 1", "1 0", "0 1", "0 0"),
    labels = c("In both",
      paste0("Official only (not fished ", footprint_year, "+)"),
      paste0("Fished ", footprint_year, "+ only (outside official)"),
      "Neither"))) |>
  ggplot() +
  geom_tile(aes(X, Y, fill = agreement), width = cell_km, height = cell_km) +
  geom_sf(data = footprint_poly, fill = NA, colour = "grey20", linewidth = 0.2) +
  geom_point(aes(X, Y), data = ports, shape = 17, size = 2) +
  geom_text(aes(X, Y, label = port), data = ports, hjust = -0.15, size = 3) +
  scale_fill_manual(values = c("#2c7fb8", "#a6cee3", "#f03b20", "#fdae61"),
    name = NULL, drop = FALSE) +
  coord_sf(xlim = range(footprint_cells$X) + c(-20, 190),
    ylim = range(footprint_cells$Y) + c(-20, 20), expand = FALSE) +
  labs(x = "X (km)", y = "Y (km)",
    title = paste0("Official 2012 trawl footprint vs fishing from ",
      footprint_year, " (", cell_km, " km cells)"),
    subtitle = "Outline: official footprint polygon") +
  theme_light() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 2))
ggsave(paste0(fig_dir, "/frame-footprint-official-", cell_km, "km.png"),
  width = 7, height = 8.5)

frame |>
  group_by(year) |>
  summarise(prop_sampled = mean(sampled), n_cells = n()) |>
  print(n = Inf)

# Index grid: frame cells, reference month, unit effort (offset 0). Used for
# the mesh sensitivity and the final indices.
grid <- frame |>
  select(X, Y, year, log_depth, month)

# Save -------------------------------------------------------------------------

saveRDS(
  list(
    d = as.data.frame(d), # catch data: one row per tow
    frame = frame, # sampling frame: one row per eligible cell x year
    grid = grid, # index prediction grid
    years = years,
    ref_month = ref_month,
    xy_offset = xy_offset, # add back to X/Y for true UTM km
    depth_range = depth_range,
    fleet_vessels = fleet_vessels,
    settings = list(
      species = species, final_year = final_year, utm_crs = utm_crs,
      cell_km = cell_km, min_yrs_cell = min_yrs_cell,
      min_positive_trips = min_positive_trips,
      min_yrs_with_trips = min_yrs_with_trips, max_hours = max_hours,
      depth_quantiles = depth_quantiles, majors = majors, ports = ports,
      footprint_year = footprint_year, footprint_file = footprint_file
    ),
    created = Sys.time()
  ),
  out_file
)
message("Saved ", out_file, ": ", nrow(d), " tows, ", nrow(frame),
  " frame rows, years ", min(years), "-", max(years))
