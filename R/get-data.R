
remotes::install_github("pbs-assess/gfdata")
library(gfdata)

d <- get_cpue_index(species = "redbanded rockfish", min_cpue_year = 1996)
saveRDS(d, "data/cpue-index-data.rds")



