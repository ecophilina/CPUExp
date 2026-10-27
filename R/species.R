# Species and folder setup shared by all scripts
#
# The species comes from the CPUE_SPECIES environment variable (default below),
# so both species can run side by side without editing scripts:
#   CPUE_SPECIES="pacific ocean perch" Rscript R/fit-joint.R <model name>
#
# Species files live in "<species>/" (the common name, lower case):
#   data/    cleaned data from R/data-prep.R (cpue-index-clean.rds)
#   output/  naive fits, indices, tables, logs
#   models/  joint fits; subfolders hold archived fits (see find_fit())
#   figs/    figures
# data/ at the project root holds inputs shared by all species: the raw gfdata
# extract (all species) and the trawl footprint shapefile.

species <- tolower(Sys.getenv("CPUE_SPECIES", "redbanded rockfish"))
species_dir <- species
if (!dir.exists(species_dir)) message("Creating folders for new species: ", species)
data_dir <- file.path(species_dir, "data")
out_dir <- file.path(species_dir, "output")
model_dir <- file.path(species_dir, "models")
fig_dir <- file.path(species_dir, "figs")
for (dir in c(data_dir, out_dir, model_dir, fig_dir)) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}
rm(dir)
message("Species: ", species)

# Cleaned data from R/data-prep.R, checked against the species setting
read_prepped <- function(file = file.path(data_dir, "cpue-index-clean.rds")) {
  if (!file.exists(file)) stop(file, " not found: run R/data-prep.R for ", species)
  prepped <- readRDS(file)
  if (!identical(tolower(prepped$settings$species), species)) {
    stop(file, " is for ", prepped$settings$species, ", not ", species)
  }
  prepped
}

# Where a fit is cached: joint fits (with a preference/sampling frame) in
# models/, naive fits in output/
fit_file <- function(name, cutoff) {
  dir <- if (grepl("^naive", name)) out_dir else model_dir
  file.path(dir, paste0("fit-", name, "-cutoff", cutoff, ".rds"))
}

# An existing fit: the cached file if there is one, otherwise the single match
# in a subfolder (e.g. models/6km footprint models/)
find_fit <- function(name, cutoff) {
  f <- fit_file(name, cutoff)
  if (file.exists(f)) return(f)
  hits <- list.files(dirname(f), recursive = TRUE, full.names = TRUE)
  hits <- hits[basename(hits) == basename(f)]
  if (length(hits) != 1) {
    stop(if (length(hits)) "Several fits named " else "No fit named ",
      basename(f), " in ", dirname(f), if (length(hits)) {
        paste0(":\n  ", paste(hits, collapse = "\n  "))
      })
  }
  hits
}
