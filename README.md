# CPUExp
Test sdmTMB preferential sampling models on commercial trawl CPUE data
(redbanded rockfish, Pacific ocean perch).

## Species
Each species has its own folder (`redbanded rockfish/`, `pacific ocean perch/`)
with `data/`, `output/`, `models/`, and `figs/`. Scripts pick the species from
the `CPUE_SPECIES` environment variable (see `R/species.R`; default redbanded
rockfish), run from the project root:

    CPUE_SPECIES="pacific ocean perch" Rscript R/data-prep.R

Shared inputs stay in `data/`: the all-species gfdata extract
(`R/get-data.R`) and the 2012 trawl footprint shapefile.
