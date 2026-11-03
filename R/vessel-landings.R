# Vessels with substantial landings of the species
#
# Ranks groundfish trawl vessels by landed catch over the CPUE years from a
# gfdata get_catch() extract in <species>/data/, and flags those in the CPUE
# data (R/data-prep.R). Writes <species>/output/vessel-landings.csv.
#
# Usage (from the project root):
#   CPUE_SPECIES="pacific ocean perch" Rscript R/vessel-landings.R pop-catch-data.rds

library(dplyr)

source("R/species.R")
catch_file <- file.path(data_dir, commandArgs(trailingOnly = TRUE)[1])
prepped <- read_prepped()
cpue_vessels <- levels(prepped$d$vessel) # vessel registration numbers
years <- prepped$years

catch <- readRDS(catch_file) |>
  filter(year %in% years, fishery_sector == "GROUNDFISH TRAWL",
    !is.na(vessel_registration_number))

vessels <- catch |>
  group_by(vessel_registration_number) |>
  summarise(
    vessel_name = names(which.max(table(vessel_name))),
    landed_t = sum(landed_kg) / 1e3,
    pct_midwater = 100 * sum(landed_kg[gear == "MIDWATER TRAWL"]) / sum(landed_kg),
    n_years = n_distinct(year),
    first_year = min(year),
    last_year = max(year),
    years_10t = sum(tapply(landed_kg, year, sum) >= 1e4), # years landing >= 10 t
    .groups = "drop"
  ) |>
  arrange(desc(landed_t)) |>
  mutate(
    pct_landed = 100 * landed_t / sum(landed_t),
    cum_pct = cumsum(pct_landed),
    in_cpue_data = vessel_registration_number %in% cpue_vessels
  )

write.csv(vessels, file.path(out_dir, "vessel-landings.csv"), row.names = FALSE)
message(nrow(vessels), " trawl vessels, ", min(years), "-", max(years), "; ",
  sum(vessels$pct_landed >= 1), " land >= 1% of the total; ",
  which(vessels$cum_pct >= 90)[1], " make up 90%")
