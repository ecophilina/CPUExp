# Parameter table: naive vs joint fits
#
# Fixed effects and random-effect parameters of the catch model (encounter and
# positive components) and, for joint fits, the sampling model, side by side.
# Fits that don't exist yet for the species are skipped.
#
# Usage (from the project root):
#   CPUE_SPECIES="pacific ocean perch" Rscript R/07-param-table.R
# Writes <species>/output/param-table.csv (all terms, long format)

library(dplyr)
library(tidyr)
library(sdmTMB)

source("R/species.R")
mesh_cutoff <- 15

fits <- c(
  "naive, IID" = "naive_tv_iid-delta_gengamma",
  "naive, AR(1)" = "naive_tv_ar1-delta_gengamma",
  "joint, IID" = "rw_iidbase_xi_expected_tv_dgengamma",
  "joint, AR(1)" = "rw_iidbase_xi_expected_tv_ar1_dgengamma"
)
has_fit <- function(fit_name) {
  !inherits(try(find_fit(fit_name, mesh_cutoff), silent = TRUE), "try-error")
}
fits <- fits[vapply(fits, has_fit, logical(1))]

fit_params <- function(fit_name) {
  fit <- readRDS(find_fit(fit_name, mesh_cutoff))
  part <- c("encounter", "positive")
  catch <- lapply(1:2, function(k) {
    bind_rows(
      tidy(fit, "fixed", model = k, silent = TRUE),
      tidy(fit, "ran_pars", model = k, silent = TRUE) |>
        mutate(term = ifelse(is.na(group_name), term, paste(group_name, term)))
    ) |>
      mutate(part = part[k])
  })
  sampling <- if (!is.null(fit$preferential)) {
    bind_rows(
      tidy(fit, "fixed", model = "sampling", silent = TRUE),
      tidy(fit, "ran_pars", model = "sampling", silent = TRUE)
    ) |>
      mutate(part = "sampling")
  }
  bind_rows(catch, sampling) |>
    select(part, term, estimate, std.error, conf.low, conf.high) |>
    mutate(converged = isTRUE(sanity(fit, silent = TRUE)$all_ok))
}

params <- lapply(fits, fit_params) |>
  bind_rows(.id = "model") |>
  mutate(model = factor(model, names(fits)))
write.csv(params, file.path(out_dir, "param-table.csv"), row.names = FALSE)

# Display: estimate (95% CI), month and sampling-year effects left out
fmt <- function(x) formatC(x, digits = 3, format = "fg", flag = "#")
params |>
  filter(!grepl("^month|year", term)) |>
  mutate(value = paste0(fmt(estimate), " (", fmt(conf.low), ", ", fmt(conf.high), ")")) |>
  select(part, term, model, value) |>
  pivot_wider(names_from = model, values_from = value, values_fill = "") |>
  as.data.frame() |>
  print(right = FALSE)
