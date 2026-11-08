# CPUExp
Test sdmTMB preferential sampling models on commercial trawl CPUE data
(redbanded rockfish, Pacific ocean perch).

## Species
Each species has its own folder (`redbanded rockfish/`, `pacific ocean perch/`)
with `data/`, `output/`, `models/`, and `figs/`. Scripts pick the species from
the `CPUE_SPECIES` environment variable (see `R/species.R`; default redbanded
rockfish), run from the project root:

    CPUE_SPECIES="pacific ocean perch" Rscript R/01-data-prep.R

Shared inputs stay in `data/`: the all-species gfdata extract
(`R/00-get-data.R`) and the 2012 trawl footprint shapefile.

## Scripts
Numbered in the order they are needed; scripts with the same number are
alternatives or can run in any order. `R/species.R` (species setup and helper
functions) is sourced by the others and not run on its own.

- `00-get-data.R`: all-species gfdata extract
- `01-data-prep.R`: cleaned data, sampling frame, and index grid per species
- `02-cpue-index.R`: catch and sampling models, mesh, and model specs; loads fits
- `03-family-comparison.R`: naive (catch-only) fits
- `04-fit-joint.R`, `04-run-joint-fits.sh`: joint fits, one or several at a time
- `05-extra-optimization.R`: extra Newton steps on a fit that needs them
- `06-ar1-indices.R`, `06-joint-indices.R`: indices and comparison plots
- `07-param-table.R`: parameter table across fits
- `08-*.R`: diagnostics and exploration (depth, fields, vessels, maps)
