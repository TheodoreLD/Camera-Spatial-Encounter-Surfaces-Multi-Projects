# Native integration test: synthetic camera-month data, no private survey files.
# Requires the same packages as the analysis library.
script <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))])
project <- dirname(dirname(normalizePath(script[1L])))
source(file.path(project, "scripts", "statistical_validation_helpers.R"))
source(file.path(project, "scripts", "validation_models.R"))
for (e in parse(file.path(project, "scripts", "wolf_encounter_surface_lib.R"))) {
  if (is.call(e) && identical(e[[1L]], as.name("<-")) &&
      is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.name("function"))) eval(e)
}
suppressPackageStartupMessages({
  library(INLA); library(dplyr); library(tidyr); library(readr)
  library(sf); library(terra); library(ggplot2); library(viridis); library(scales)
})
sf::sf_use_s2(FALSE)
INLA::inla.setOption(num.threads = "1:1")
set.seed(500)
SURVEY_ID <- "synthetic"
SURVEY_LABEL <- "Synthetic integration survey"
SURVEY_PREFIX <- "synthetic"
TARGET_LABEL <- "synthetic"
TARGET_TITLE <- "Synthetic"
TARGET_SLUG <- "synthetic"
SURVEY_YEAR <- "2024"
FINAL_FAMILY <- "nbinomial"
FINAL_MODEL_NAME <- "synthetic_nb"
EPSG_UTM <- 32634L
OUTPUT_DIR <- tempfile("synthetic-inla-")
dir.create(OUTPUT_DIR)
PPC_NSIM <- CV_NSIM <- MAP_NSIM <- 64L
MAP_CHUNK_SIZE <- 7L
CV_K <- 2L
INLA_SAMPLE_SEED <- 201L
SPATIAL_FOLD_SEED <- 101L
RUN_PROFILE <- "integration"
RUN_DEPENDENCE_COMPARISON <- TRUE
MORAN_NPERM <- 99L
MORAN_ALPHA <- 0.05
PRIOR_INTERCEPT_PREC <- 0.16
PRIOR_MONTH_LOG_RATE_RATIO_PREC <- 1
PRIOR_ZI_LOGIT_MEAN <- qlogis(0.05)
PRIOR_ZI_LOGIT_PREC <- 1 / 1.5^2
PRIOR_NB_LOGSIZE_MEAN <- log(2)
PRIOR_NB_LOGSIZE_PREC <- 0.25
settings <- list(cell_size_m = 2000, pred_buffer_m = 500, mesh_cutoff_m = 100,
  mesh_max_edge = c(2000, 6000), mesh_offset = c(2000, 6000), fix_range_m = NULL,
  prior_range_m = c(3000, 0.5), prior_sigma = c(1, 0.05), use_month_effect = TRUE,
  month_reference = "2024-01", month_prediction = "2024-03")
MONTH_REFERENCE <- settings$month_reference
MONTH_PREDICTION <- settings$month_prediction
locations <- expand.grid(longitude = 21 + (0:3) * 0.025, latitude = 60 + (0:2) * 0.025)
dat <- expand.grid(plotID = sprintf("c%02d", 1:12), month = c("2024-01", "2024-03", "2024-04"),
                   stringsAsFactors = FALSE)
idx <- match(dat$plotID, sprintf("c%02d", 1:12))
dat$longitude <- locations$longitude[idx]
dat$latitude <- locations$latitude[idx]
dat$deploymentID <- paste(dat$plotID, dat$month)
dat$start <- as.POSIXct(paste0(dat$month, "-01"), tz = "UTC")
dat$end <- dat$start + 20 * 86400
dat$total_effort_days <- 20
camera_deviation <- rnorm(12, 0, 0.4)
dat$wolf_events <- rnbinom(nrow(dat), mu = 20 * exp(-1 + camera_deviation[idx]), size = 4)
dat$wolf_events_per_100_days <- 100 * dat$wolf_events / dat$total_effort_days
dat <- add_month_design(dat, settings)
PRIOR_INTERCEPT_MEAN <- training_intercept_mean(dat$wolf_events, dat$total_effort_days)

base <- fit_final_model(dat, settings, FINAL_FAMILY)
samples <- posterior_samples_safe(base$fit, 64L)
diag <- compute_diagnostics(base$fit, samples, base$model_dat, base$obs_index,
                            base$camera_sf, FINAL_FAMILY, write_files = FALSE)
map <- make_prediction_outputs(base, diag, settings, FINAL_FAMILY)
stopifnot(nrow(map$pred_sf) == nrow(base$pred_sf), all(is.finite(map$pred_sf$sd)),
          all(map$pred_sf$q95 >= map$pred_sf$q05), all(map$pred_sf$mean > 0),
          map$method == "joint_posterior_mesh_projection")
# Both likelihood hyperparameters must actually be present in user-scale samples.
stopifnot(all(extract_size_samples(base$fit, samples, FINAL_FAMILY) > 0))
cv <- spatial_block_cv(dat, settings, FINAL_FAMILY, K = 2)
stopifnot(nrow(cv$row) == nrow(dat), nrow(cv$camera) == 12L,
          nrow(cv$calibration) == 2L, nrow(cv$fold_metrics) == 4L)
prior <- read_csv(path_out("synthetic_final_spatial_cv_training_priors.csv"), show_col_types = FALSE)
for (f in prior$fold) {
  train_ids <- cv$folds$plotID[cv$folds$fold != f]
  train <- dat[dat$plotID %in% train_ids, ]
  stopifnot(abs(prior$intercept_prior_mean[prior$fold == f] -
    training_intercept_mean(train$wolf_events, train$total_effort_days)) < 1e-10)
}
td <- temporal_autocorrelation_diagnostics(diag$model_dat)
comparison <- run_dependence_comparison(dat, settings, FINAL_FAMILY, base, diag, cv, map, td)
stopifnot(all(comparison$status == "completed"), nrow(comparison) == 3L)
write_validation_report(dat, diag, cv, map, td, comparison)
write_manifest(dat, diag, cv, map, td, comparison)

# Refit a small ZINB model to verify strict user-scale probability extraction.
zi <- fit_dependence_model(dat, settings, "zeroinflatednbinomial1", "none", base$pred_sf)
zi_samples <- posterior_samples_safe(zi$fit, 32L)
p <- extract_pi_samples(zi$fit, zi_samples, "zeroinflatednbinomial1")
stopifnot(length(p) == 32L, all(p >= 0 & p < 1))
cat("PASS: synthetic native INLA mapping, strict CV, camera IID/AR1 comparisons, reports, and ZINB samples.\n")
cat("R/INLA versions:", as.character(getRversion()), as.character(packageVersion("INLA")), "\n")
