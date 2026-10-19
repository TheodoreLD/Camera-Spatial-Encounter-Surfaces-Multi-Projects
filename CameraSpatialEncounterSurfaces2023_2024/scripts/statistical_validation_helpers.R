# Pure R helpers shared by the analysis and dependency-free regression tests.
# Sourcing this file does not fit models, read survey inputs, or write outputs.

with_validation_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed)
  force(code)
}

strict_spatial_folds <- function(coords, K, seed = 101L) {
  coords <- as.matrix(coords)
  if (ncol(coords) != 2L || any(!is.finite(coords)) ||
      length(K) != 1L || !is.finite(K) || K != as.integer(K) || K < 2L ||
      nrow(coords) <= K || nrow(unique(coords)) < K) {
    stop("Spatial CV needs finite 2D coordinates, at least K distinct locations, and more than K cameras.")
  }
  # A single scaling constant preserves Euclidean geometry, including a survey
  # with one constant coordinate. Never substitute random folds on failure.
  centered <- sweep(coords, 2L, colMeans(coords), "-")
  scale_m <- max(abs(centered))
  if (scale_m <= 0) stop("Spatial CV cannot partition coincident camera locations.")
  fit <- tryCatch(
    with_validation_seed(seed, kmeans(centered / scale_m, centers = K, nstart = 20)),
    error = function(e) stop("Spatial fold construction failed: ", conditionMessage(e))
  )
  if (fit$ifault != 0L || length(unique(fit$cluster)) != K)
    stop("Spatial k-means did not converge to K nonempty folds.")
  # Co-located cameras must never appear on opposite sides of a held-out fold.
  location_key <- paste(coords[, 1], coords[, 2], sep = "|")
  if (any(vapply(split(fit$cluster, location_key), function(z) length(unique(z)) > 1L,
                 logical(1)))) stop("Co-located cameras were split between folds.")
  as.integer(fit$cluster)
}

training_intercept_mean <- function(y, effort) {
  if (!length(y) || length(y) != length(effort) ||
      any(!is.finite(y) | y < 0) || any(!is.finite(effort) | effort <= 0))
    stop("Training counts and exposure must be finite, aligned, and nonnegative/positive.")
  # For an all-zero training fold use half an event over its own total exposure.
  # This only centers a broad prior; no held-out outcome enters the calculation.
  log(max(sum(y), 0.5) / sum(effort))
}

latent_component_draws <- function(samples, component, ids) {
  if (!length(samples) || !length(ids)) stop("No posterior samples or component indices.")
  out <- vapply(samples, function(s) {
    rn <- trimws(rownames(s$latent))
    prefix <- paste0(component, ":")
    hit <- which(startsWith(rn, prefix))
    idx <- suppressWarnings(as.integer(substring(rn[hit], nchar(prefix) + 1L)))
    selected <- match(ids, idx)
    if (anyNA(selected) || anyDuplicated(idx))
      stop("Missing or ambiguous posterior component: ", component)
    val <- as.numeric(s$latent[hit[selected], 1L])
    if (any(!is.finite(val))) stop("Nonfinite posterior component: ", component)
    val
  }, numeric(length(ids)))
  matrix(out, nrow = length(ids), ncol = length(samples))
}

user_hyperparameter_draws <- function(samples, pattern, label) {
  vapply(samples, function(s) {
    hit <- grep(pattern, names(s$hyperpar), ignore.case = TRUE)
    if (length(hit) != 1L) stop("Missing or ambiguous posterior hyperparameter: ", label)
    value <- as.numeric(s$hyperpar[hit])
    if (!is.finite(value)) stop("Nonfinite posterior hyperparameter: ", label)
    value
  }, numeric(1))
}

joint_month_factors <- function(samples, months, weights, reference) {
  if (length(months) != length(weights) || any(!is.finite(weights) | weights < 0) ||
      abs(sum(weights) - 1) > 1e-8) stop("Invalid month aggregation weights.")
  gamma <- matrix(0, nrow = length(months), ncol = length(samples))
  for (i in seq_along(months)) {
    if (months[i] != reference)
      gamma[i, ] <- latent_component_draws(
        samples, month_term_name(months[i]), 1L)[1L, ]
  }
  factors <- as.numeric(crossprod(weights, exp(gamma)))
  if (any(!is.finite(factors) | factors <= 0)) stop("Invalid posterior month factors.")
  list(factor = factors, rate_ratios = exp(gamma))
}

summarise_rate_draws <- function(spatial_eta, intercept, month_factor, pi) {
  spatial_eta <- as.matrix(spatial_eta)
  nsim <- ncol(spatial_eta)
  if (nsim < 2L || any(lengths(list(intercept, month_factor, pi)) != nsim) ||
      any(!is.finite(c(spatial_eta, intercept, month_factor, pi))) ||
      any(month_factor <= 0 | pi < 0 | pi >= 1)) stop("Invalid joint map draws.")
  rate <- 100 * exp(sweep(spatial_eta, 2L, intercept, "+"))
  rate <- sweep(rate, 2L, month_factor * (1 - pi), "*")
  if (any(!is.finite(rate))) stop("Nonfinite rates from posterior map draws.")
  m <- rowMeans(rate)
  s <- apply(rate, 1L, sd)
  data.frame(mean = m, sd = s, cv = s / m,
             q05 = apply(rate, 1L, quantile, 0.05, names = FALSE),
             q95 = apply(rate, 1L, quantile, 0.95, names = FALSE))
}

interval_score_90 <- function(y, lo, hi) {
  if (length(y) != length(lo) || length(y) != length(hi) ||
      any(!is.finite(c(y, lo, hi))) || any(lo > hi)) stop("Invalid prediction interval.")
  hi - lo + 20 * pmax(lo - y, 0) + 20 * pmax(y - hi, 0)
}

aggregate_joint_predictive_moments <- function(fitted, conditional_variance, group) {
  fitted <- as.matrix(fitted)
  conditional_variance <- as.matrix(conditional_variance)
  if (!identical(dim(fitted), dim(conditional_variance)) || nrow(fitted) != length(group) ||
      anyNA(group) || ncol(fitted) < 2L) stop("Invalid joint predictive moments.")
  totals <- rowsum(fitted, group, reorder = TRUE)
  conditional_totals <- rowsum(conditional_variance, group, reorder = TRUE)
  data.frame(plotID = rownames(totals), fitted_count = rowMeans(totals),
             fit_var = rowMeans(conditional_totals) + apply(totals, 1L, var),
             row.names = NULL)
}

cv_calibration_summary <- function(rows, cameras, target = 0.90) {
  summarise_level <- function(dat, level) {
    if (!nrow(dat) || anyNA(dat$covered_90)) stop("Incomplete held-out predictions.")
    fold_coverage <- tapply(dat$covered_90, dat$fold, mean)
    data.frame(level = level, n = nrow(dat), folds = length(fold_coverage),
               nominal_coverage = target, observed_coverage = mean(dat$covered_90),
               min_fold_coverage = min(fold_coverage), max_fold_coverage = max(fold_coverage),
               mean_interval_width = mean(dat$hi90 - dat$lo90),
               mean_interval_score = mean(interval_score_90(dat$y, dat$lo90, dat$hi90)),
               review_required = mean(dat$covered_90) < target)
  }
  rbind(summarise_level(rows, "model_row"), summarise_level(cameras, "camera"))
}

dependence_design <- function(model_dat, kind = c("none", "camera_iid", "camera_ar1")) {
  kind <- match.arg(kind)
  if (kind == "none") return(list(indices = data.frame(row.names = seq_len(nrow(model_dat))),
                                  term = character(), n_camera = 0L, n_time = 0L))
  camera <- match(as.character(model_dat$plotID), sort(unique(as.character(model_dat$plotID))))
  if (anyNA(camera)) stop("Missing camera identifier in dependence model.")
  if (kind == "camera_iid") return(list(
    indices = data.frame(camera_effect = camera), n_camera = max(camera), n_time = 0L,
    term = paste0("f(camera_effect, model='iid', n=", max(camera),
                  ", hyper=list(prec=list(prior='pc.prec', param=c(1,0.05))))")))
  date <- as.Date(paste0(model_dat$month, "-01"))
  if (anyNA(date)) stop("AR1 comparison requires valid calendar months.")
  serial <- as.integer(format(date, "%Y")) * 12L + as.integer(format(date, "%m"))
  time <- serial - min(serial) + 1L
  if (max(time) < 2L) stop("AR1 comparison requires multiple calendar months.")
  list(indices = data.frame(camera_time = time, camera_rep = camera),
       n_camera = max(camera), n_time = max(time),
       term = paste0("f(camera_time, model='ar1', replicate=camera_rep, n=", max(time),
                     ", hyper=list(prec=list(prior='pc.prec', param=c(1,0.05)),",
                     "rho=list(prior='pc.cor0', param=c(0.7,0.5))))"))
}
