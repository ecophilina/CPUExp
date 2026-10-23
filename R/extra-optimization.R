# Extra Newton optimization on a cached fit, to clear small gradients that fail
# sanity() when the optimizer has otherwise converged.
#
# Usage (from the project root):
#   Rscript R/extra-optimization.R <model name> [mesh cutoff, default 15]
# Reads  <dir>/fit-<model>-cutoff<cutoff>.rds
# Writes <dir>/fit-<model>_newton-cutoff<cutoff>.rds (the original is kept)
# where <dir> is output/ for naive fits and models/ for joint fits
#
# Uses ~9 GB of memory for a joint fit; run at most two at a time.

args <- commandArgs(trailingOnly = TRUE)
nm <- args[1]
if (is.na(nm)) stop("Usage: Rscript R/extra-optimization.R <model name> [cutoff]")
cutoff <- if (length(args) > 1) args[2] else 15

suppressMessages(library(sdmTMB))

fit_dir <- if (grepl("^naive", nm)) "output" else "models"
f_in <- file.path(fit_dir, paste0("fit-", nm, "-cutoff", cutoff, ".rds"))
f_out <- file.path(fit_dir, paste0("fit-", nm, "_newton-cutoff", cutoff, ".rds"))
if (!file.exists(f_in)) stop(f_in, " not found")
if (file.exists(f_out)) stop(f_out, " already exists")

fit <- readRDS(f_in)
g0 <- abs(fit$gradients)
t0 <- Sys.time()
fit2 <- run_extra_optimization(fit, nlminb_loops = 0, newton_loops = 1)
cat(nm, ": extra optimization took",
  round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
saveRDS(fit2, f_out)

g1 <- abs(fit2$gradients)
cat(sprintf("max|grad| %.3g -> %.3g; n > 0.001: %d -> %d; objective %.4f -> %.4f\n",
  max(g0), max(g1), sum(g0 > 1e-3), sum(g1 > 1e-3),
  fit$model$objective, fit2$model$objective))
s <- sanity(fit2)
cat("all_ok:", s$all_ok, "\n")
