# Fit one model from `model_specs` in R/cpue-index.R and cache it.
#
# Usage (from the project root):
#   Rscript R/fit-joint.R <model name>
# e.g. Rscript R/fit-joint.R constant_xi_expected_dgamma
#
# Lets several models run as separate processes (see R/run-joint-fits.sh);
# R/cpue-index.R then loads the cached fits. Uses the same data, mesh,
# fit_cpue(), and cache file (fit_file()) as the main script.

args <- commandArgs(trailingOnly = TRUE)
nm <- args[1]
if (is.na(nm)) stop("Usage: Rscript R/fit-joint.R <model name> [--single-phase]")
# --single-phase: fit_cpue(multiphase = FALSE). Does not lower peak memory for
# time-varying-b models (~17.8 GB single-phase vs ~16 GB multiphase).
multiphase <- !"--single-phase" %in% args

# Evaluate R/cpue-index.R up to and including `model_families`, skipping the
# mesh sensitivity block: loads the clean data, builds the main mesh, and
# defines the formulas, fit_cpue(), fit_file(), pref(), and the model lists.
ex <- parse("R/cpue-index.R", keep.source = FALSE)
is_assign_to <- function(e, name) {
  is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
    identical(e[[2]], as.name(name))
}
is_mesh_test <- function(e) {
  is.call(e) && identical(e[[1]], as.name("if")) &&
    identical(e[[2]], as.name("run_mesh_test"))
}
stop_at <- which(vapply(ex, is_assign_to, logical(1), name = "model_families"))[1]
stopifnot(!is.na(stop_at))
pdf(NULL)
for (i in seq_len(stop_at)) {
  if (!is_mesh_test(ex[[i]])) eval(ex[[i]], globalenv())
}

if (!nm %in% names(model_specs)) {
  stop(nm, " is not in model_specs: ", paste(names(model_specs), collapse = ", "))
}
f <- fit_file(nm, mesh_cutoff)
if (file.exists(f) && !refit) stop(f, " already exists (set refit = TRUE to refit)")

family <- model_families[[nm]]
if (is.null(family)) family <- delta_gengamma()
spec <- model_specs[[nm]]
cat("Fitting", nm, "| family", paste(family$family, collapse = "/"),
  "| target", spec$target, "| coefficient", spec$coefficient,
  "| baseline", spec$baseline, "| spatial", spec$spatial,
  "| multiphase", multiphase, "\n  ->", f, "\n")

t0 <- Sys.time()
# Time-varying catch models use fit_cpue_tv() (no multiphase argument)
if (exists("model_tv") && nm %in% model_tv) {
  cat("  catch model: fit_cpue_tv() (time-varying intercept)\n")
  fit <- fit_cpue_tv(spec, family = family)
} else {
  fit <- fit_cpue(spec, family = family, multiphase = multiphase)
}
cat("FIT DONE:", nm, "in",
  round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
saveRDS(fit, f)

s <- sanity(fit)
g <- abs(fit$gradients)
cat(sprintf("all_ok: %s | optimizer: %s | max|grad| %.3g (%s)\n", s$all_ok,
  fit$model$message, max(g), names(fit$model$par)[which.max(g)]))
print(tidy(fit, model = "sampling", conf.int = TRUE, silent = TRUE) |>
  subset(!grepl("year", term)))
print(tidy(fit, "ran_pars", model = "sampling", silent = TRUE))
