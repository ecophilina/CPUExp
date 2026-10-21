# Observation-family comparison for the naive (catch-only) CPUE model
#
# DRAFT. Fits the naive model from R/cpue-index.R with each family that the
# joint preferential-sampling model also supports, so the winner can be carried
# into the joint fits. Candidates are limited to families that can handle
# catch_kg (continuous, with zeros):
#   - delta_gamma(), delta_lognormal()                  (standard delta)
#   - delta_gamma/lognormal(type = "poisson-link")      (Poisson-link delta)
#   - tweedie()
# Excluded: poisson()/nbinom2() (counts; kg isn't integer) and Gamma()/
# lognormal() (no zeros).
#
# Comparisons: convergence (sanity), AIC (valid here: all families model the
# same response with a point mass at zero), simulation-based residual checks,
# and sensitivity of the index and its trend.
#
# Runtime: roughly 15 min per fit at a 15 km mesh, plus simulation and
# bias-corrected indices.

# Setup from the main script ---------------------------------------------------

# Runs R/cpue-index.R up to (not including) the mesh sensitivity section: loads
# data/cpue-index-clean.rds (run R/data-prep.R first) and defines
# build_mesh(), the formulas, fit_cpue(), and fit_file().
ex <- parse("R/cpue-index.R", keep.source = FALSE)
stop_at <- which(vapply(ex, function(e) {
  is.call(e) && identical(e[[1]], as.name("if")) &&
    identical(e[[2]], as.name("run_mesh_test"))
}, logical(1)))
stopifnot(length(stop_at) == 1)
for (i in seq_len(stop_at - 1)) eval(ex[[i]], globalenv())

mesh <- build_mesh(mesh_cutoff)
n_sim <- 100 # simulated datasets per family for residual checks
# FALSE: compare only families with cached fits (e.g., while others are still
# fitting elsewhere); rerun later to add them. Per-family checks are cached, so
# a rerun only processes the new families.
fit_missing <- FALSE

families <- list(
  delta_gamma = delta_gamma(),
  delta_lognormal = delta_lognormal(),
  delta_gamma_plink = delta_gamma(type = "poisson-link"),
  delta_lognormal_plink = delta_lognormal(type = "poisson-link"),
  tweedie = tweedie()
)

# Fits -------------------------------------------------------------------------

fam_fits <- list()
fam_secs <- c()
for (nm in names(families)) {
  # Every family is fitted and cached under its own name (delta_gamma
  # included), so all fits are guaranteed to use the same data for AIC.
  f <- fit_file(paste0("naive-", nm), mesh_cutoff)
  if (file.exists(f) && !refit) {
    fam_fits[[nm]] <- readRDS(f)
    fam_secs[[nm]] <- NA_real_
    next
  }
  if (!fit_missing) {
    message("Skipping ", nm, ": no cached fit (set fit_missing = TRUE to fit)")
    next
  }
  message("Fitting naive model with ", nm)
  secs <- system.time(
    fit <- tryCatch(fit_cpue(family = families[[nm]]), error = function(e) {
      message("  failed: ", conditionMessage(e))
      NULL
    })
  )[["elapsed"]]
  if (!is.null(fit)) {
    saveRDS(fit, f)
    fam_fits[[nm]] <- fit
    fam_secs[[nm]] <- secs
  }
}

# Convergence and AIC ----------------------------------------------------------

fam_summary <- lapply(names(fam_fits), function(nm) {
  fit <- fam_fits[[nm]]
  ll <- logLik(fit)
  data.frame(
    family = nm,
    converged = isTRUE(sanity(fit, silent = TRUE)$all_ok),
    fit_minutes = round(fam_secs[[nm]] / 60, 1), # NA when loaded from cache
    n_par = attr(ll, "df"),
    logLik = as.numeric(ll),
    AIC = AIC(fit)
  )
}) |>
  bind_rows() |>
  mutate(dAIC = AIC - min(AIC[converged]))
print(fam_summary)

# Only compare families whose fits pass sanity()
ok <- fam_summary$family[fam_summary$converged]

# Simulation-based checks ------------------------------------------------------

# Simulate from each fit (fields drawn from their joint precision, "mle-mvn")
# and compare features that families differ on: the share of zeros and the
# upper tail of positive catches. Then DHARMa-style randomized quantile
# residuals. With ~170k tows the formal tests will reject almost anything, so
# judge the QQ plots and the size of the discrepancies, not the p-values.
obs <- d$catch_kg

checks_dir <- file.path("output", "family-checks")
dir.create(checks_dir, showWarnings = FALSE)

# Simulation checks, QQ plot, and index for one family; cached per family
family_checks <- function(nm) {
  f <- file.path(checks_dir,
    paste0("family-checks-", nm, "-cutoff", mesh_cutoff, ".rds"))
  if (file.exists(f) && !refit) return(readRDS(f))
  message("Checking ", nm)
  fit <- fam_fits[[nm]]
  s <- simulate(fit, nsim = n_sim, type = "mle-mvn", seed = 2026, silent = TRUE)

  q_pos <- function(x, p) quantile(x[x > 0], p, names = FALSE)
  res <- dharma_residuals(s, fit, plot = FALSE, return_DHARMa = TRUE)
  png(file.path("figs", paste0("family-qq-", nm, ".png")),
    width = 5, height = 5, units = "in", res = 150)
  DHARMa::plotQQunif(res, testUniformity = FALSE, testOutliers = FALSE,
    testDispersion = FALSE, main = nm)
  dev.off()

  sim <- data.frame(
    family = nm,
    obs_prop_zero = mean(obs == 0),
    sim_prop_zero = mean(s == 0),
    obs_q90_pos = q_pos(obs, 0.90),
    sim_q90_pos = median(apply(s, 2, q_pos, p = 0.90)),
    obs_q99_pos = q_pos(obs, 0.99),
    sim_q99_pos = median(apply(s, 2, q_pos, p = 0.99)),
    obs_max = max(obs),
    sim_max = median(apply(s, 2, max)),
    zero_inflation_ratio = DHARMa::testZeroInflation(res, plot = FALSE)$statistic[[1]],
    dispersion_ratio = DHARMa::testDispersion(res, plot = FALSE)$statistic[[1]]
  )
  index <- get_index(fit, newdata = grid, area = cell_km^2,
    bias_correct = TRUE, predict_args = list(re_form_iid = NA))

  out <- list(sim = sim, index = index)
  saveRDS(out, f)
  out
}

checks <- lapply(setNames(ok, ok), family_checks)
sim_checks <- bind_rows(lapply(checks, `[[`, "sim"))
print(sim_checks)

# Index sensitivity ------------------------------------------------------------

fam_index <- lapply(checks, `[[`, "index") |>
  bind_rows(.id = "family") |>
  group_by(family) |>
  mutate(
    rel_est = est / exp(mean(log_est)),
    rel_lwr = lwr / exp(mean(log_est)),
    rel_upr = upr / exp(mean(log_est))
  ) |>
  ungroup()

ggplot(fam_index, aes(year, rel_est, colour = family, fill = family)) +
  geom_ribbon(aes(ymin = rel_lwr, ymax = rel_upr), alpha = 0.1, colour = NA) +
  geom_line() +
  labs(x = "Year", y = "Relative index (geometric mean = 1)") +
  theme_light()
ggsave("figs/family-index.png", width = 8, height = 5)

fam_trend <- fam_index |>
  group_by(family) |>
  summarise(
    slope = coef(lm(log_est ~ year))[[2]],
    mean_cv = mean(sqrt(exp(se^2) - 1))
  )

# Summary ----------------------------------------------------------------------

# Prefer the lowest-AIC converged family unless the simulation checks show it
# misses the zeros or the upper tail badly; then check whether the index
# trend actually depends on the choice.
fam_results <- fam_summary |>
  left_join(sim_checks, by = "family") |>
  left_join(fam_trend, by = "family") |>
  arrange(dAIC)
print(fam_results)
write.csv(fam_results, "output/family-comparison.csv", row.names = FALSE)
