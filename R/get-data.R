
remotes::install_github("pbs-assess/gfdata")
library(gfdata)

# All trawl fishing events with every species caught, shared by all species:
# R/data-prep.R pulls out the species set in R/species.R
d <- get_cpue_index(min_cpue_year = 1996)
saveRDS(d, "data/cpue-index-data.rds")
