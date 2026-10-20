# Run from the repository root:
# Rscript CameraSpatialEncounterSurfaces2023_2024/tests/test_statistical_validation.R
# Base R only. No private data, INLA, or geospatial packages are required.
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
project <- if (length(script)) dirname(dirname(normalizePath(script[1L]))) else getwd()
source(file.path(project, "scripts", "statistical_validation_helpers.R"))
source(file.path(project, "scripts", "validation_models.R"))
for (file in list.files(file.path(project, "tests"), "\\.R$", full.names = TRUE))
  invisible(parse(file))
# Parse all code, then load function definitions without executing the pipeline.
for (file in list.files(file.path(project, "scripts"), "\\.R$", full.names = TRUE)) {
  expr <- parse(file)
  if (basename(file) == "wolf_encounter_surface_lib.R") {
    for (e in expr) {
      if (is.call(e) && identical(e[[1L]], as.name("<-")) &&
          ((is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name("function"))) ||
           startsWith(as.character(e[[2L]]), "PAT_"))) eval(e)
    }
  }
}

checks <- 0L
check <- function(x, label) {
  if (!isTRUE(x)) stop("FAIL: ", label)
  checks <<- checks + 1L
}
errors <- function(code) inherits(tryCatch({force(code); NULL}, error = identity), "error")

# Joint nonlinear transformations must preserve uncertainty and dependence.
ns <- 4000L
set.seed(2026)
gamma <- rnorm(ns, 0, 0.7)
intercept <- -gamma + rnorm(ns, -2, 0.1)
pi <- plogis(0.8 * gamma - 2)
spatial <- rbind(rep(0, ns), rep(log(2), ns))
factor <- 0.3 + 0.7 * exp(gamma)
got <- summarise_rate_draws(spatial, intercept, factor, pi)
manual <- 100 * exp(intercept) * factor * (1 - pi)
check(isTRUE(all.equal(got$mean[1], mean(manual))), "joint mean matches explicit draw transformation")
check(isTRUE(all.equal(got$sd[1], sd(manual))), "joint SD matches explicit draw transformation")
check(abs(got$mean[2] / got$mean[1] - 2) < 1e-12, "spatial projection remains on rate scale")
check(abs(got$mean[1] - 100 * mean(exp(intercept)) * mean(factor) * mean(1 - pi)) > 0.1,
      "joint result differs from multiplication of marginal means")
pi_only <- summarise_rate_draws(matrix(0, 1, ns), rep(0, ns), rep(1, ns), pi)
month_only <- summarise_rate_draws(matrix(0, 1, ns), rep(0, ns), factor, rep(0, ns))
check(pi_only$sd > 0 && month_only$sd > 0, "ZI and month uncertainty propagate with fixed spatial predictor")
check(errors(summarise_rate_draws(spatial, intercept, factor, rep(1, ns))), "invalid ZI draws rejected")

# Component extraction handles zero-padded INLA labels and fails on missing data.
samples <- lapply(c(-0.5, 0.5), function(g) list(
  latent = matrix(c(g, 1), ncol = 1, dimnames = list(c("month_2024_03:001", "spatial:001"), NULL)),
  hyperpar = c("zero-probability parameter" = 0.2)))
mix <- joint_month_factors(samples, c("2024-01", "2024-03"), c(0.25, 0.75), "2024-01")
check(isTRUE(all.equal(mix$factor, 0.25 + 0.75 * exp(c(-0.5, 0.5)))), "month weights transformed within draws")
check(errors(latent_component_draws(samples, "missing", 1L)), "missing posterior component fails")
check(errors(user_hyperparameter_draws(samples, "unavailable", "test")), "no point-estimate hyperparameter fallback")
check(all(user_hyperparameter_draws(samples, "zero", "test") == 0.2), "user-scale hyperparameter extraction")
check(all(extract_pi_samples(NULL, samples, "zeroinflatednbinomial1") == 0.2),
      "production ZI extractor uses shared hyperparameter pattern")
nb_samples <- lapply(samples, function(s) {s$hyperpar <- c("size for nbinomial" = 4); s})
check(all(extract_size_samples(NULL, nb_samples, "nbinomial") == 4),
      "production NB extractor uses shared hyperparameter pattern")

# Folds use coordinates only and do not depend on unrelated RNG use.
coords <- cbind(rep(0:5, each = 2), rep(c(0, 1), 6))
set.seed(8); before <- .Random.seed
folds <- strict_spatial_folds(coords, 3, seed = 101)
check(identical(before, .Random.seed), "fold construction preserves surrounding RNG state")
invisible(runif(30))
check(identical(folds, strict_spatial_folds(coords, 3, seed = 101)), "spatial folds reproducible")
check(length(unique(strict_spatial_folds(cbind(1:10, 0), 2))) == 2, "constant coordinate is supported")
check(errors(strict_spatial_folds(matrix(0, 10, 2), 2)), "coincident locations fail clearly")
check(errors(strict_spatial_folds(coords, 1)), "one-fold CV rejected")
original_kmeans <- kmeans
kmeans <- function(...) stop("deliberate clustering failure")
check(errors(strict_spatial_folds(coords, 3)), "clustering failure never produces random fallback")
kmeans <- original_kmeans

y <- c(2, 4, 1000); effort <- c(10, 20, 1)
center <- training_intercept_mean(y[1:2], effort[1:2])
y[3] <- 100000
check(identical(center, training_intercept_mean(y[1:2], effort[1:2])), "held-out outcomes cannot affect fold prior")
check(is.finite(training_intercept_mean(c(0, 0), c(10, 20))), "zero-event training fold has finite broad-prior center")
check(errors(training_intercept_mean(1, 0)), "invalid training exposure rejected")
PRIOR_INTERCEPT_MEAN <- 99
PRIOR_INTERCEPT_PREC <- 0.16
PRIOR_MONTH_LOG_RATE_RATIO_PREC <- 1
check(make_control_fixed("intercept", intercept_mean = center)$mean$intercept == center,
      "fold-specific prior overrides full-data center without mutation")
check(PRIOR_INTERCEPT_MEAN == 99, "CV prior does not mutate full-data setting")

# Dependence indexing preserves calendar gaps and holds unseen cameras separate.
dat <- data.frame(plotID = c("a", "a", "b", "b"),
                  month = c("2024-01", "2024-03", "2024-01", "2024-03"))
dep <- dependence_design(dat, "camera_ar1")
check(identical(dep$indices$camera_time, c(1L, 3L, 1L, 3L)), "AR1 retains missing-month gaps")
check(dep$n_camera == 2 && dep$n_time == 3, "AR1 replicates are camera-specific")
ar1_call <- str2lang(dep$term)
check(identical(eval(ar1_call$values), 1:3),
      "AR1 formula declares the complete calendar including unobserved months")
check(errors(dependence_design(transform(dat, month = "2024-01"), "camera_ar1")), "single-month AR1 rejected")

heldout <- data.frame(fold = rep(1:2, each = 5), y = c(rep(1, 8), 5, 5),
                      lo90 = 0, hi90 = 2, covered_90 = c(rep(TRUE, 8), FALSE, FALSE))
cal <- cv_calibration_summary(heldout, heldout)
check(all(cal$review_required) && all(cal$observed_coverage == 0.8), "undercoverage triggers explicit review")
check(all(cal$min_fold_coverage == 0.6) && all(cal$max_fold_coverage == 1), "fold heterogeneity retained")
check(interval_score_90(5, 0, 2) == 62, "interval score penalizes missed outcomes")
diag <- list(ppc_total_pass = TRUE, ppc_zero_pass = TRUE, ppc_max_pass = TRUE, moran_pass = TRUE)
status <- validation_review_status(diag, list(calibration = cal), NULL)
check(status$status == "review_required", "legacy diagnostic pass cannot hide undercoverage")
check(validation_review_status(diag, NULL, NULL)$status == "review_required", "skipped CV is not a validation pass")
good_cv <- list(calibration = transform(cal, review_required = FALSE))
check(validation_review_status(diag, good_cv, NULL, data.frame(status = "failed"))$status == "review_required",
      "failed model comparison cannot yield an overall pass")
fitted <- rbind(1:4, 1:4, 2 * (1:4))
mom <- aggregate_joint_predictive_moments(fitted, matrix(1, 3, 4), c("a", "a", "b"))
check(abs(mom$fit_var[1] - (2 + var(2 * (1:4)))) < 1e-12,
      "camera variance retains positive covariance between repeated observations")

cat(sprintf("PASS: %d statistical regression checks; all analysis scripts parsed.\n", checks))
