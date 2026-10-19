# Joint map summaries and explicit dependence-model comparisons.
# Functions are called by the shared workflow; sourcing has no side effects.

joint_surface_summary <- function(fit_obj, model_dat, settings, family,
                                   nsim = MAP_NSIM, seed = INLA_SAMPLE_SEED + 10L) {
  fit <- fit_obj$fit
  terms <- fixed_effect_terms(model_dat)
  selection <- c(list(spatial = 0L), setNames(rep(list(0L), length(terms)), terms))
  samples <- posterior_samples_safe(fit, nsim, selection = selection, seed = seed)
  spatial <- latent_component_draws(samples, "spatial", seq_len(fit_obj$spde_obj$spde$n.spde))
  intercept <- latent_component_draws(samples, "intercept", 1L)[1L, ]
  pi <- extract_pi_samples(fit, samples, family)

  weights <- aggregate(model_dat$total_effort_days, list(month = model_dat$month), sum)
  names(weights)[2L] <- "effort_days"
  weights <- weights[order(weights$month), , drop = FALSE]
  weights$weight <- weights$effort_days / sum(weights$effort_days)
  if (isTRUE(settings$use_month_effect)) {
    mix <- joint_month_factors(samples, weights$month, weights$weight,
                               settings$month_reference)
    prediction_gamma <- if (settings$month_prediction == settings$month_reference) {
      rep(0, nsim)
    } else latent_component_draws(samples, month_term_name(settings$month_prediction), 1L)[1L, ]
    relative_factor <- mix$factor / exp(prediction_gamma)
    weights$mean_rate_ratio_to_reference <- rowMeans(mix$rate_ratios)
    factor <- mix$factor
  } else {
    factor <- relative_factor <- rep(1, nsim)
    weights$mean_rate_ratio_to_reference <- 1
  }
  annualization <- list(
    factor = mean(relative_factor), factor_sd = sd(relative_factor),
    factor_q05 = unname(quantile(relative_factor, 0.05)),
    factor_q95 = unname(quantile(relative_factor, 0.95)), weights = weights,
    label = if (isTRUE(settings$use_month_effect)) {
      paste("effort-weighted sampled-month", SURVEY_YEAR, "surface")
    } else "single-period encounter-frequency surface"
  )
  A <- INLA::inla.spde.make.A(fit_obj$spde_obj$mesh,
                              loc = sf::st_coordinates(fit_obj$pred_sf))
  chunks <- split(seq_len(nrow(A)), ceiling(seq_len(nrow(A)) / MAP_CHUNK_SIZE))
  summaries <- lapply(chunks, function(index) {
    eta <- as.matrix(A[index, , drop = FALSE] %*% spatial)
    summarise_rate_draws(eta, intercept, factor, pi)
  })
  list(summary = do.call(rbind, summaries), annualization = annualization,
       method = "joint_posterior_mesh_projection", nsim = nsim, seed = seed)
}

fit_dependence_model <- function(model_dat, settings, family, kind, pred_sf) {
  dat <- model_dat %>% mutate(y = as.integer(wolf_events), intercept = 1)
  dep <- dependence_design(dat, kind)
  fixed_terms <- fixed_effect_terms(dat)
  camera_sf <- camera_to_utm(camera_summary_from_model(dat))
  coords <- st_coordinates(camera_sf)
  spde_obj <- build_spatial(coords, settings)
  obs_sf <- st_transform(st_as_sf(dat, coords = c("longitude", "latitude"),
                                  crs = 4326, remove = FALSE), EPSG_UTM)
  A <- INLA::inla.spde.make.A(spde_obj$mesh, loc = st_coordinates(obs_sf))
  effects <- cbind(as.data.frame(dat[, fixed_terms, drop = FALSE]), dep$indices)
  stack <- INLA::inla.stack(
    tag = "obs", data = list(y = dat$y, e = dat$total_effort_days),
    A = list(A, 1), effects = list(spatial = spde_obj$s_index, fixed = effects))
  stack_data <- INLA::inla.stack.data(stack)
  formula <- as.formula(paste("y ~ 0 +", paste(
    c(fixed_terms, "f(spatial, model=spde_obj$spde)", dep$term), collapse = " + ")))
  fit <- INLA::inla(
    formula, family = family, data = stack_data, E = stack_data$e,
    control.predictor = list(A = INLA::inla.stack.A(stack), compute = TRUE, link = 1),
    control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE),
    control.fixed = make_control_fixed(fixed_terms),
    control.family = make_control_family(family), verbose = FALSE)
  list(fit = fit, model_dat = dat, camera_sf = camera_sf, spde_obj = spde_obj,
       pred_sf = pred_sf, obs_index = INLA::inla.stack.index(stack, "obs")$data)
}

run_dependence_comparison <- function(model_dat, settings, family, baseline,
                                       baseline_diag, baseline_cv, prediction,
                                       temporal_diag) {
  if (!RUN_DEPENDENCE_COMPARISON || is.null(baseline_cv)) return(NULL)
  kinds <- "camera_iid"
  if (length(unique(model_dat$month)) > 1L) kinds <- c(kinds, "camera_ar1")
  summaries <- list()
  lag1_values <- function(td) {
    if (is.null(td) || is.null(td$lag_summary)) return(c(r = NA_real_, p = NA_real_))
    lag <- td$lag_summary[td$lag_summary$lag == 1L, , drop = FALSE]
    if (!nrow(lag)) return(c(r = NA_real_, p = NA_real_))
    c(r = lag$correlation[1L], p = lag$p_value[1L])
  }
  comparison_row <- function(kind, fit, diag, cv, td, map) {
    camera <- cv$camera
    lag <- lag1_values(td)
    data.frame(model = kind, status = "completed", waic = safe_value(fit$waic$waic),
      row_mean_lpd = mean(cv$row$lpd), camera_coverage_90 = mean(camera$covered_90),
      camera_interval_score_90 = mean(interval_score_90(camera$y, camera$lo90, camera$hi90)),
      camera_rmse_count = sqrt(mean((camera$y - camera$Ey)^2)),
      moran_p = diag$moran_p, residual_lag1_r = lag[["r"]], residual_lag1_p = lag[["p"]],
      map_rmse_difference = sqrt(mean((map$mean - prediction$pred_sf$mean)^2)),
      map_mean_absolute_difference = mean(abs(map$mean - prediction$pred_sf$mean)),
      mean_map_sd = mean(map$sd),
      conditional_camera_deviation = 0, failure = NA_character_)
  }
  summaries[[1L]] <- comparison_row("configured", baseline$fit, baseline_diag,
                                    baseline_cv, temporal_diag, prediction$pred_sf)
  failures <- character()
  for (kind in kinds) {
    result <- tryCatch({
      cat(sprintf("\n[%s] dependence sensitivity: %s\n", SURVEY_PREFIX, kind))
      obj <- fit_dependence_model(model_dat, settings, family, kind, baseline$pred_sf)
      samples <- posterior_samples_safe(obj$fit, PPC_NSIM, seed = INLA_SAMPLE_SEED + 20L)
      diag <- compute_diagnostics(obj$fit, samples, obj$model_dat, obj$obs_index,
                                   obj$camera_sf, family, write_files = FALSE)
      td <- temporal_autocorrelation_diagnostics(
        diag$model_dat, final_model = paste0(FINAL_MODEL_NAME, "_", kind))
      cv <- spatial_block_cv(model_dat, settings, family, K = CV_K,
                             folds = baseline_cv$folds, dependence = kind,
                             output_label = kind)
      map <- joint_surface_summary(obj, model_dat, settings, family)$summary
      delta <- data.frame(grid_id = baseline$pred_sf$grid_id,
        x = st_coordinates(baseline$pred_sf)[, 1L], y = st_coordinates(baseline$pred_sf)[, 2L],
        configured_mean = prediction$pred_sf$mean, candidate_mean = map$mean,
        mean_difference = map$mean - prediction$pred_sf$mean,
        configured_sd = prediction$pred_sf$sd, candidate_sd = map$sd,
        configured_q05 = prediction$pred_sf$q05, configured_q95 = prediction$pred_sf$q95,
        candidate_q05 = map$q05, candidate_q95 = map$q95)
      readr::write_csv(delta, path_out(paste0(SURVEY_PREFIX, "_", kind, "_map_sensitivity.csv")))
      row <- comparison_row(kind, obj$fit, diag, cv, td, map)
      rm(samples, obj, diag, cv, map)
      gc(verbose = FALSE)
      row
    }, error = function(e) {
      failures <<- c(failures, paste(kind, conditionMessage(e), sep = ": "))
      row <- summaries[[1L]]
      row[1L, ] <- NA
      row$model <- kind
      row$status <- "failed"
      row$failure <- conditionMessage(e)
      row
    })
    summaries[[length(summaries) + 1L]] <- result
  }
  out <- do.call(rbind, summaries)
  readr::write_csv(out, path_out(paste0(SURVEY_PREFIX, "_dependence_model_comparison.csv")))
  report <- c(
    "Dependence-model sensitivity and calibration review:",
    "Configured spatial model versus camera IID and (for multi-month data) replicated camera AR1 alternatives.",
    "Likelihood, spatial priors, training-only intercept centering, and saved spatial folds are held fixed.",
    "IID/AR1 precision uses PC P(SD > 1)=0.05; AR1 correlation uses pc.cor0 with u=0.7, alpha=0.5.",
    "AR1 steps are actual calendar months; missing months retain their gaps.",
    "New held-out cameras have unobserved random effects integrated through joint posterior simulation.",
    "Map comparisons set additional camera deviations to zero: standardized spatial surfaces, not population-marginal new-camera means.",
    "The map SD/interval summaries are conditional on those additional deviations, integrating spatial, fixed, and likelihood hyperparameter uncertainty.",
    "Lower interval score and higher held-out log predictive density are useful alongside coverage, RMSE, and residual diagnostics.",
    "These folds guide model refinement; selection on these results needs a further independent/nested assessment for a final performance claim.",
    "No alternative is automatically promoted, and no stability or absence-of-distortion conclusion is assumed.",
    capture.output(print(out, row.names = FALSE)),
    if (length(failures)) c("Incomplete comparisons (must be resolved):", failures) else NULL)
  writeLines(report, path_out(paste0(SURVEY_PREFIX, "_dependence_model_comparison_report.txt")))
  invisible(out)
}

validation_review_status <- function(diag, cv, temporal_diag,
                                      dependence_comparison = NULL, comparison_required = FALSE) {
  reasons <- diagnostic_failures(diag)
  if (is.null(cv)) reasons <- c(reasons, "spatial CV not run")
  else {
    low <- cv$calibration$level[cv$calibration$review_required]
    if (length(low)) reasons <- c(reasons, paste("held-out 90% undercoverage:", paste(low, collapse = ", ")))
  }
  if (!is.null(temporal_diag) && !is.null(temporal_diag$lag_summary)) {
    lag <- temporal_diag$lag_summary[temporal_diag$lag_summary$lag == 1L, , drop = FALSE]
    if (nrow(lag) && is.finite(lag$p_value[1L]) && lag$p_value[1L] < 0.05)
      reasons <- c(reasons, "residual within-camera temporal association needs review")
  }
  if ((comparison_required && is.null(dependence_comparison)) ||
      (!is.null(dependence_comparison) && any(dependence_comparison$status != "completed")))
    reasons <- c(reasons, "dependence-model sensitivity comparison incomplete")
  list(status = if (length(reasons)) "review_required" else "no_flag_detected",
       reasons = reasons)
}
