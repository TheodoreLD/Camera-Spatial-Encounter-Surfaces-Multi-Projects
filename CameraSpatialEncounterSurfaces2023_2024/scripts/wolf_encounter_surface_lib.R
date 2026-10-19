###############################################################################
# Wolf relative encounter frequency: shared INLA-SPDE analysis library
# -----------------------------------------------------------------------------
# Purpose
#   Single source of truth for the ordered workflow used by every surface in
#   this repository. All runs -- the three primary wolf surveys, the two
#   human-activity companions, and the two single-month March 2024 surfaces --
#   perform the same analyses, diagnostics, model comparison, sensitivity checks,
#   and map products from this one file. They differ only in the configuration a
#   thin runner sets before sourcing this library.
#
# How to run
#   Do not run this file directly. Run one of the survey runners, which set the
#   survey-specific configuration and then source() this library:
#     scripts/run_road_2023.R             wolf, road 2023        (NB, months)
#     scripts/run_road_2024.R             wolf, road 2024        (ZINB, months)
#     scripts/run_forest_2024.R           wolf, forest 2024      (NB, months)
#     scripts/run_road_2023_human.R       human activity, 2023   (NB, months)
#     scripts/run_road_2024_human.R       human activity, 2024   (NB, months)
#     scripts/run_road_march2024_wolf.R   wolf, March 2024       (Poisson, spatial-only)
#     scripts/run_road_march2024_human.R  human activity, Mar24  (Poisson, spatial-only)
#
# Runner contract (globals each runner must set before source()):
#   SURVEY_ID, SURVEY_LABEL, SURVEY_PREFIX  -- identity / output naming
#   FINAL_FAMILY, FINAL_MODEL_NAME          -- likelihood ("poisson" | "nbinomial"
#                                              | "zeroinflatednbinomial1") and name
#   SURVEY_DATA_SHAPE                        -- "road_deployments" | "forest_flat"
#   input_files_required                     -- input CSV name(s) for detection
#   settings                                 -- mesh, spatial PC priors, month
#                                              reference, and use_month_effect
#   (optional) TARGET_LABEL  -- report/map label for the modelled detections
#                              (default "wolf"); WOLF_NAMES -- scientificName
#                              strings to count (default Canis lupus);
#                              OUTPUT_SUBPATH, FOREST_INPUT_FILES, EPSG_UTM
#
# Model
#   Response: independent detection-event count per camera-month row -- wolf, or
#   a human-activity index (people + vehicles), selected via WOLF_NAMES.
#   Exposure: camera effort split across calendar months, passed to INLA via E.
#   Likelihood: Poisson, negative-binomial, or zero-inflated negative-binomial,
#   chosen per runner after a WAIC/DIC model comparison, with a spatial SPDE
#   field. Multi-month surveys add calendar-month fixed effects; single-month
#   surveys set use_month_effect = FALSE (intercept + spatial field only).
#
# Interpretation
#   Predictions are relative encounter frequency: expected events per 100
#   camera-days. They are not abundance, density, occupancy, or population size.
#
# Analyses run for every survey
#   * Exploratory checks before model selection.
#   * Candidate model comparison among spatial Poisson/NB/ZINB variants (WAIC/DIC).
#   * Prior-influence screen and targeted prior-sensitivity refits.
#   * Final mapping fit + observed-data diagnostic fit with joint posterior samples.
#   * Camera- and row-level posterior predictive checks, PIT, Pearson dispersion.
#   * Two-sided residual Moran's I permutation test (run-profile scaled).
#   * Residual temporal-autocorrelation diagnostics.
#   * Spatial block cross-validation with per-fold mesh rebuild.
#   * Mesh sensitivity checks.
#   * Posterior mean / SD / CV encounter-frequency maps (effort-weighted and
#     annualized for multi-month surveys; single-period for single-month ones).
###############################################################################


## 01. Runtime setup (survey configuration comes from the runner) --------------

# This file is a shared library. A thin runner (scripts/run_*.R) sets the
# survey-specific globals, then sources this file. Fail early and clearly if it
# was sourced without that configuration.
.required_runner_globals <- c(
  "SURVEY_ID", "SURVEY_LABEL", "SURVEY_PREFIX",
  "FINAL_FAMILY", "FINAL_MODEL_NAME",
  "SURVEY_DATA_SHAPE", "input_files_required", "settings"
)
.missing_runner_globals <- .required_runner_globals[
  !vapply(.required_runner_globals, exists, logical(1))
]
if (length(.missing_runner_globals)) {
  stop("wolf_encounter_surface_lib.R is a shared library, not an entry point. ",
       "It must be sourced from a runner that first sets: ",
       paste(.missing_runner_globals, collapse = ", "), ". ",
       "Run scripts/run_road_2023.R, scripts/run_road_2024.R, or ",
       "scripts/run_forest_2024.R instead.")
}

script_file <- tryCatch({
  ofile <- sys.frames()[[1]]$ofile
  if (is.null(ofile)) NA_character_ else normalizePath(ofile, winslash = "/", mustWork = FALSE)
}, error = function(e) NA_character_)

cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (is.na(script_file) && length(cmd_file)) {
  script_file <- normalizePath(sub("^--file=", "", cmd_file[[1]]),
                               winslash = "/", mustWork = FALSE)
}

script_dir <- if (is.na(script_file)) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  dirname(script_file)
}

source(file.path(script_dir, "statistical_validation_helpers.R"))
source(file.path(script_dir, "validation_models.R"))

# Check whether a directory (or its data/ subfolder) holds the required input files.
has_input_layout <- function(dir) {
  all(file.exists(file.path(dir, input_files_required))) ||
    all(file.exists(file.path(dir, "data", input_files_required)))
}

project_candidates <- unique(normalizePath(c(script_dir, getwd()),
                                           winslash = "/", mustWork = FALSE))
if (basename(script_dir) %in% c("scripts", "R")) {
  project_candidates <- unique(c(
    normalizePath(dirname(script_dir), winslash = "/", mustWork = FALSE),
    project_candidates
  ))
}
detected_project <- project_candidates[
  vapply(project_candidates, has_input_layout, logical(1))
]
default_project_dir <- if (length(detected_project)) detected_project[[1]] else getwd()

PROJECT_DIR <- Sys.getenv("WOLF_PROJECT_DIR", unset = default_project_dir)
PROJECT_DIR <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = TRUE)

DATA_DIR <- Sys.getenv("WOLF_DATA_DIR", unset = "")
if (!nzchar(DATA_DIR)) {
  data_subdir <- file.path(PROJECT_DIR, "data")
  DATA_DIR <- if (dir.exists(data_subdir)) data_subdir else PROJECT_DIR
}
DATA_DIR <- normalizePath(DATA_DIR, winslash = "/", mustWork = FALSE)

if (!exists("OUTPUT_SUBPATH")) OUTPUT_SUBPATH <- SURVEY_PREFIX
OUTPUT_DIR <- Sys.getenv(
  "WOLF_OUTPUT_DIR",
  unset = file.path(PROJECT_DIR, "outputs", OUTPUT_SUBPATH)
)
OUTPUT_DIR <- normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Spatial CRS for modelling and maps.
EPSG_UTM <- 32634L

# Detection labels to model, matched exactly against scientificName in the
# source files. A runner may override this before sourcing (to model a
# different species or detection category) by setting WOLF_NAMES first.
if (!exists("WOLF_NAMES")) WOLF_NAMES <- c("Canis_lupus", "Canis lupus")

# Human-readable label for the modelled detections, used in report prose, map
# titles, and raster band names. A runner may override it (defaults to "wolf").
if (!exists("TARGET_LABEL")) TARGET_LABEL <- "wolf"
TARGET_TITLE <- paste0(toupper(substr(TARGET_LABEL, 1, 1)),
                       substr(TARGET_LABEL, 2, nchar(TARGET_LABEL)))
TARGET_SLUG <- gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", tolower(TARGET_LABEL)))

# Set TRUE only if you explicitly want this script to install missing packages.
INSTALL_MISSING_PACKAGES <- FALSE

# Runtime profile:
#   quick    : fast testing; fewer posterior samples and no spatial CV
#   balanced : recommended default
#   final    : heavier publication rerun
RUN_PROFILE <- tolower(Sys.getenv("WOLF_RUN_PROFILE", unset = "balanced"))
if (!RUN_PROFILE %in% c("quick", "balanced", "final")) {
  stop("WOLF_RUN_PROFILE must be one of: quick, balanced, final.")
}

# Number of simulations and folds                        
PPC_NSIM <- switch(RUN_PROFILE, quick = 200L, balanced = 750L, final = 1500L)
CV_NSIM <- switch(RUN_PROFILE, quick = 150L, balanced = 300L, final = 600L)
RUN_SPATIAL_CV <- RUN_PROFILE != "quick"
CV_K <- switch(RUN_PROFILE, quick = 3L, balanced = 4L, final = 5L)
MAP_NSIM <- switch(RUN_PROFILE, quick = 200L, balanced = 750L, final = 1500L)
MAP_CHUNK_SIZE <- 1000L
INLA_SAMPLE_SEED <- 201L
SPATIAL_FOLD_SEED <- 101L
# Dependence alternatives are diagnostic comparisons, never automatic replacements.
RUN_DEPENDENCE_COMPARISON <- RUN_PROFILE != "quick"

# Diagnostics.
MORAN_ALPHA <- 0.05
MORAN_NPERM <- switch(RUN_PROFILE, quick = 199L, balanced = 499L, final = 999L)

# Survey identity (SURVEY_ID / SURVEY_LABEL / SURVEY_PREFIX), likelihood
# (FINAL_FAMILY / FINAL_MODEL_NAME), data shape (SURVEY_DATA_SHAPE), and the
# mesh / spatial-PC-prior `settings` list are all provided by the runner (see
# the runner contract in the header). Everything below is shared by all surveys.

MONTH_REFERENCE <- settings$month_reference
MONTH_PREDICTION <- settings$month_prediction
SURVEY_YEAR <- substr(MONTH_REFERENCE, 1, 4)

# Fixed-effect and likelihood priors (shared defaults for every survey; a runner
# may override any of them, and the human runners tune the spatial PC priors).
# The intercept prior mean is filled after loading the data, from the crude
# observed daily event rate; its SD is deliberately broad but finite. A runner
# may override any prior below by setting it before sourcing this library.
PRIOR_INTERCEPT_MEAN <- NA_real_
if (!exists("PRIOR_INTERCEPT_PREC")) PRIOR_INTERCEPT_PREC <- 1 / 2.5^2
# Month effects: weakly informative, SD = 1 on the log-rate-ratio scale.
if (!exists("PRIOR_MONTH_LOG_RATE_RATIO_PREC")) PRIOR_MONTH_LOG_RATE_RATIO_PREC <- 1
# Zero inflation: skeptical center at 5% structural-zero probability, SD = 1.5 on the logit scale.
if (!exists("PRIOR_ZI_LOGIT_MEAN")) PRIOR_ZI_LOGIT_MEAN <- qlogis(0.05)
if (!exists("PRIOR_ZI_LOGIT_PREC")) PRIOR_ZI_LOGIT_PREC <- 1 / 1.5^2
# Negative-binomial size prior on log(size): Gaussian(log(2), SD 2).
if (!exists("PRIOR_NB_LOGSIZE_MEAN")) PRIOR_NB_LOGSIZE_MEAN <- log(2)
if (!exists("PRIOR_NB_LOGSIZE_PREC")) PRIOR_NB_LOGSIZE_PREC <- 1 / 2^2

set.seed(1)


## 02. Packages ---------------------------------------------------------------

cat("Project directory: ", PROJECT_DIR, "\n", sep = "")
cat("Data directory:    ", DATA_DIR, "\n", sep = "")
cat("Output directory:  ", OUTPUT_DIR, "\n", sep = "")
cat(sprintf(
  "Run profile:       %s | PPC_NSIM=%d | spatial_CV=%s | CV_K=%d | CV_NSIM=%d\n",
  RUN_PROFILE, PPC_NSIM, RUN_SPATIAL_CV, CV_K, CV_NSIM
))

required_packages <- c(
  "readr", "dplyr", "tidyr", "sf", "terra",
  "ggplot2", "viridis", "scales", "INLA"
)

for (p in required_packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (!INSTALL_MISSING_PACKAGES) {
      stop("Package '", p, "' is missing. Install it or set INSTALL_MISSING_PACKAGES <- TRUE.")
    }
    if (p == "INLA") {
      install.packages(
        "INLA",
        dep = TRUE,
        repos = c(getOption("repos"),
                  INLA = "https://inla.r-inla-download.org/R/stable")
      )
    } else {
      install.packages(p)
    }
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(sf)
  library(terra)
  library(ggplot2)
  library(viridis)
  library(scales)
  library(INLA)
})

sf::sf_use_s2(FALSE)
try(INLA::inla.setOption(fmesher.evolution.warn = FALSE), silent = TRUE)


## 03. General helpers --------------------------------------------------------

# Build a file path inside the data directory.
path_in <- function(...) file.path(DATA_DIR, ...)
# Build a file path inside the output directory.
path_out <- function(...) file.path(OUTPUT_DIR, ...)

# Stop with a clear error if a data frame lacks any required column.
stop_missing_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(label, " missing required column(s): ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

# Parse a timestamp column into POSIXct (UTC), tolerating several input formats.
parse_time <- function(x) {
  if (inherits(x, "POSIXt")) return(x)
  if (inherits(x, "Date")) return(as.POSIXct(x, tz = "UTC"))

  s <- as.character(x)
  s <- gsub("T", " ", s, fixed = TRUE)
  s <- sub("([+-][0-9]{2}:?[0-9]{2}|Z)$", "", s)

  out <- as.POSIXct(substr(s, 1, 19),
                    format = "%Y-%m-%d %H:%M:%S",
                    tz = "UTC")
  bad <- is.na(out) & !is.na(s) & nzchar(s)
  if (any(bad)) {
    out[bad] <- as.POSIXct(substr(s[bad], 1, 10),
                           format = "%Y-%m-%d",
                           tz = "UTC")
  }
  out
}

# Numerically stable log of the mean of exponentials.
log_mean_exp <- function(x) {
  m <- max(x, na.rm = TRUE)
  if (!is.finite(m)) return(NA_real_)
  m + log(mean(exp(x - m), na.rm = TRUE))
}

# Pearson correlation that returns NA on too few pairs or zero variance instead of failing.
safe_cor <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3 || sd(a[ok]) == 0 || sd(b[ok]) == 0) return(NA_real_)
  cor(a[ok], b[ok])
}

# Correlation-test p-value, returning NA on too few pairs, zero variance, or error.
safe_cor_p_value <- function(a, b, min_pairs = 5L) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < min_pairs || sd(a[ok]) == 0 || sd(b[ok]) == 0) return(NA_real_)
  suppressWarnings(tryCatch(
    stats::cor.test(a[ok], b[ok])$p.value,
    error = function(e) NA_real_
  ))
}

# Coerce a value to a UTC Date.
as_utc_date <- function(x) {
  if (inherits(x, "POSIXt")) return(as.Date(x, tz = "UTC"))
  as.Date(x)
}


## 04. Family, prior, and hyperparameter helpers ------------------------------

# Normalise a family alias (case, punctuation, short forms) to INLA canonical family strings.
fit_family <- function(family) {
  key <- gsub("[^a-z0-9]", "", tolower(family))
  if (key %in% c("poisson")) {
    "poisson"
  } else if (key %in% c("nbinomial", "negativebinomial", "nb")) {
    "nbinomial"
  } else if (key %in% c("zeroinflatedpoisson1", "zip", "zip1")) {
    "zeroinflatedpoisson1"
  } else if (key %in% c("zeroinflatednbinomial1", "zinb", "zinb1")) {
    "zeroinflatednbinomial1"
  } else {
    family
  }
}

# TRUE if the family is zero-inflated.
is_zi <- function(family) grepl("zeroinflated", family)
# TRUE if the family is negative-binomial.
is_nb <- function(family) grepl("nbinomial", family)

# Row name this family takes in the model-comparison table. The comparison set
# always fits the three spatial-month candidates under those fixed names, so a
# survey configured as poisson_spatial (single month, no month term) still maps
# onto the poisson_spatial_month comparison row.
comparison_row_for_family <- function(family) {
  switch(fit_family(family),
         poisson = "poisson_spatial_month",
         nbinomial = "nb_spatial_month",
         zeroinflatednbinomial1 = "zinb_spatial_month",
         zeroinflatedpoisson1 = "zip_spatial_month",
         NA_character_)
}

# Short display tag for a comparison row name.
comparison_tag <- function(row_name) {
  switch(row_name,
         poisson_spatial_month = "Poisson",
         nb_spatial_month = "NB",
         zinb_spatial_month = "ZINB",
         zip_spatial_month = "ZIP",
         row_name)
}

# Return a usable NB size (finite, positive), else a large value approximating Poisson.
fam_nb_size <- function(size) {
  ifelse(is.finite(size) & size > 0, size, 1e6)
}

# Mean of the (optionally zero-inflated) response given mu, pi, and family.
fam_mean <- function(mu, pi, family, size = NA_real_) {
  p <- ifelse(is_zi(family) & is.finite(pi), pi, 0)
  (1 - p) * mu
}

# Variance of the (optionally zero-inflated) response given mu, pi, size, and family.
fam_var <- function(mu, pi, family, size = NA_real_) {
  p <- ifelse(is_zi(family) & is.finite(pi), pi, 0)
  base_var <- if (is_nb(family)) {
    s <- fam_nb_size(size)
    mu + mu^2 / s
  } else {
    mu
  }
  (1 - p) * base_var + p * (1 - p) * mu^2
}

# Log probability mass of an observed count under the (optionally zero-inflated) family.
fam_logpmf <- function(y, mu, pi, family, size = NA_real_) {
  p <- ifelse(is_zi(family) & is.finite(pi), pi, 0)
  base_prob <- if (is_nb(family)) {
    dnbinom(y, mu = mu, size = fam_nb_size(size))
  } else {
    dpois(y, mu)
  }
  log(pmax(p * (y == 0) + (1 - p) * base_prob,
           .Machine$double.xmin))
}

# Simulate counts from the (optionally zero-inflated) family.
fam_sim <- function(mu, pi, family, size = NA_real_) {
  y <- if (is_nb(family)) {
    rnbinom(length(mu), mu = mu, size = fam_nb_size(size))
  } else {
    rpois(length(mu), mu)
  }
  if (is_zi(family) && is.finite(pi) && pi > 0) {
    y * (1 - rbinom(length(mu), 1, pi))
  } else {
    y
  }
}

PAT_ZPROB <- "zero.*prob|probability.*zero|zero-probability"
PAT_NB_SIZE <- "size.*nbinomial|size for|size parameter"
PAT_RANGE <- "range.*spatial|range for spatial|spatial.*range"
PAT_SIGMA <- "stdev.*spatial|stdev for spatial|standard deviation.*spatial|sigma.*spatial"

# Posterior mean of the hyperparameter whose name matches a pattern.
hyp_point <- function(fit, pattern) {
  if (is.null(fit$summary.hyperpar)) return(NA_real_)
  i <- grep(pattern, rownames(fit$summary.hyperpar), ignore.case = TRUE)
  if (length(i)) fit$summary.hyperpar[i[[1]], "mean"] else NA_real_
}

# Posterior marginal of the hyperparameter whose name matches a pattern.
hyp_marg <- function(fit, pattern) {
  if (is.null(fit$marginals.hyperpar)) return(NULL)
  i <- grep(pattern, names(fit$marginals.hyperpar), ignore.case = TRUE)
  if (length(i)) fit$marginals.hyperpar[[i[[1]]]] else NULL
}

# Posterior mean of the NB size hyperparameter (NA if absent or invalid).
nb_size_point <- function(fit) {
  x <- hyp_point(fit, PAT_NB_SIZE)
  if (is.finite(x) && x > 0) x else NA_real_
}

# Build the INLA control.family list (NB size and/or zero-inflation priors) for a family.
make_control_family <- function(family) {
  hyper <- list()

  if (is_nb(family)) {
    hyper$size <- list(
      prior = "gaussian",
      param = c(PRIOR_NB_LOGSIZE_MEAN, PRIOR_NB_LOGSIZE_PREC)
    )
  }

  if (is_zi(family)) {
    hyper$prob <- list(
      prior = "gaussian",
      param = c(PRIOR_ZI_LOGIT_MEAN, PRIOR_ZI_LOGIT_PREC)
    )
  }

  if (length(hyper)) list(hyper = hyper) else list()
}

# Build the INLA fixed-effect prior means and precisions (intercept plus month terms).
make_control_fixed <- function(fixed_terms = "intercept",
                               intercept_mean = PRIOR_INTERCEPT_MEAN) {
  fixed_terms <- unique(fixed_terms)

  mean_values <- rep(0, length(fixed_terms))
  names(mean_values) <- fixed_terms
  if ("intercept" %in% names(mean_values)) {
    if (!is.finite(intercept_mean)) {
      stop("PRIOR_INTERCEPT_MEAN must be set before fitting.")
    }
    mean_values[["intercept"]] <- intercept_mean
  }
  mean_prior <- as.list(mean_values)
  mean_prior$default <- 0

  prec_values <- rep(PRIOR_MONTH_LOG_RATE_RATIO_PREC, length(fixed_terms))
  names(prec_values) <- fixed_terms
  if ("intercept" %in% names(prec_values)) {
    prec_values[["intercept"]] <- PRIOR_INTERCEPT_PREC
  }

  prec_prior <- as.list(prec_values)
  prec_prior$default <- PRIOR_MONTH_LOG_RATE_RATIO_PREC

  list(mean = mean_prior, prec = prec_prior)
}

# Turn a month string into its fixed-effect column name.
month_term_name <- function(month) {
  paste0("month_", gsub("[^A-Za-z0-9]", "_", month))
}

# Recover the month string from a month fixed-effect column name.
month_from_term <- function(term) {
  gsub("_", "-", sub("^month_", "", term), fixed = TRUE)
}

# List the month fixed-effect columns present in the data.
temporal_month_terms <- function(data) {
  grep("^month_[0-9]{4}_[0-9]{2}$", names(data), value = TRUE)
}

# List all fixed-effect terms (intercept plus month columns).
fixed_effect_terms <- function(data) {
  c("intercept", temporal_month_terms(data))
}


## 05. Data preparation -------------------------------------------------------

# Check that the required input file(s) exist before loading.
validate_inputs <- function() {
  if (identical(SURVEY_DATA_SHAPE, "forest_flat")) {
    # resolve_input_file() errors clearly if the flat file cannot be found.
    resolve_input_file(FOREST_INPUT_FILES,
                       paste0("[", SURVEY_PREFIX, "] forest camera-trap input"))
    return(invisible(TRUE))
  }
  missing <- input_files_required[!file.exists(path_in(input_files_required))]
  if (length(missing)) {
    stop("Missing input file(s): ", paste(missing, collapse = ", "),
         " in ", DATA_DIR)
  }
  invisible(TRUE)
}

# Resolve and validate the reference month from settings.
month_reference_from_settings <- function(months, settings) {
  ref <- settings$month_reference
  if (is.null(ref) || !nzchar(ref)) ref <- months[[1]]
  if (!ref %in% months) {
    stop("month_reference '", ref, "' is not present in deployment months: ",
         paste(months, collapse = ", "))
  }
  ref
}

# Resolve and validate the prediction/baseline month from settings.
month_prediction_from_settings <- function(months, settings) {
  pred <- settings$month_prediction
  if (is.null(pred) || !nzchar(pred)) pred <- month_reference_from_settings(months, settings)
  if (!pred %in% months) {
    stop("month_prediction '", pred, "' is not present in deployment months: ",
         paste(months, collapse = ", "))
  }
  pred
}

# Add month fixed-effect indicator columns and record the reference/prediction month.
add_month_design <- function(model_dat, settings) {
  months <- sort(unique(model_dat$month))
  if (!isTRUE(settings$use_month_effect)) {
    # Spatial-only run: no month fixed effects (e.g. a single-month surface).
    # Record the period but add no month indicator columns.
    model_dat$month_reference <- months[[1]]
    model_dat$month_prediction <- months[[1]]
    model_dat$model_row_type <- "deployment_month_no_month_effect"
    return(model_dat)
  }
  if (length(months) < 2) {
    stop("Month effect requires at least two deployment months.")
  }

  reference_month <- month_reference_from_settings(months, settings)
  prediction_month <- month_prediction_from_settings(months, settings)

  for (m in setdiff(months, reference_month)) {
    term <- month_term_name(m)
    model_dat[[term]] <- as.integer(model_dat$month == m)
  }

  model_dat$month_reference <- reference_month
  model_dat$month_prediction <- prediction_month
  model_dat$model_row_type <- "deployment_month"
  model_dat
}

# Aggregate camera-month rows to one row per camera (effort, events, rate, coords).
camera_summary_from_model <- function(model_dat) {
  has_deployment_id <- "deploymentID" %in% names(model_dat)
  has_n_deployments <- "n_deployments" %in% names(model_dat)

  model_dat %>%
    group_by(plotID) %>%
    summarise(
      longitude = mean(longitude, na.rm = TRUE),
      latitude = mean(latitude, na.rm = TRUE),
      total_effort_days = sum(total_effort_days, na.rm = TRUE),
      wolf_events = sum(wolf_events, na.rm = TRUE),
      n_model_rows = n(),
      n_deployments = if (has_deployment_id) {
        n_distinct(deploymentID)
      } else if (has_n_deployments) {
        sum(n_deployments, na.rm = TRUE)
      } else {
        n()
      },
      wolf_events_per_100_days = 100 * wolf_events / total_effort_days,
      .groups = "drop"
    ) %>%
    filter(is.finite(longitude),
           is.finite(latitude),
           is.finite(total_effort_days),
           total_effort_days > 0) %>%
    arrange(plotID)
}

# First instant (UTC) of the calendar month containing a timestamp.
month_period_start <- function(x) {
  as.POSIXct(
    paste0(format(x, "%Y-%m", tz = "UTC"), "-01 00:00:00"),
    tz = "UTC"
  )
}

# First instant (UTC) of the following calendar month.
next_month_start <- function(x) {
  lt <- as.POSIXlt(month_period_start(x), tz = "UTC")
  lt$mon <- lt$mon + 1L
  as.POSIXct(lt, tz = "UTC")
}

# Split one deployment into camera-month segments at calendar-month boundaries.
split_one_deployment_by_month <- function(deployment_row) {
  start <- deployment_row$start[[1]]
  end <- deployment_row$end[[1]]
  if (is.na(start) || is.na(end) || end <= start) return(NULL)

  cuts <- start
  boundary <- next_month_start(start)
  while (!is.na(boundary) && boundary < end) {
    cuts <- c(cuts, boundary)
    boundary <- next_month_start(boundary)
  }
  cuts <- c(cuts, end)

  if (length(cuts) < 2L) return(NULL)

  out <- deployment_row[rep(1L, length(cuts) - 1L), , drop = FALSE]
  out$original_deployment_start <- deployment_row$start[[1]]
  out$original_deployment_end <- deployment_row$end[[1]]
  out$start <- cuts[-length(cuts)]
  out$end <- cuts[-1L]
  out$deploymentEffort <- as.numeric(difftime(out$end, out$start, units = "days"))
  out$month <- format(out$start, "%Y-%m", tz = "UTC")
  out$model_row_id <- paste0(out$deploymentID, "__", out$month)
  out
}

# Split all deployments into camera-month segments with per-segment effort.
split_deployments_by_month <- function(deployments) {
  segments <- dplyr::bind_rows(lapply(seq_len(nrow(deployments)), function(i) {
    split_one_deployment_by_month(deployments[i, , drop = FALSE])
  }))

  segments %>%
    filter(is.finite(deploymentEffort), deploymentEffort > 0) %>%
    arrange(plotID, start, end)
}

# Count target detection events falling within each camera-month segment.
count_segment_wolf_events <- function(segments, wolf_events) {
  if (!nrow(segments) || !nrow(wolf_events)) return(rep(0L, nrow(segments)))

  vapply(seq_len(nrow(segments)), function(i) {
    ev <- wolf_events[wolf_events$deploymentID == segments$deploymentID[[i]], , drop = FALSE]
    if (!nrow(ev)) return(0L)

    is_last_segment <- abs(as.numeric(difftime(
      segments$end[[i]],
      segments$original_deployment_end[[i]],
      units = "secs"
    ))) < 1e-6
    in_segment <- ev$event_start >= segments$start[[i]] &
      (ev$event_start < segments$end[[i]] |
         (is_last_segment & ev$event_start <= segments$end[[i]]))
    as.integer(sum(in_segment, na.rm = TRUE))
  }, integer(1))
}

# Load and assemble road-camera camera-month rows from deployment and observation tables.
load_road_survey <- function(settings) {
  dep <- readr::read_csv(path_in(input_files_required[[1]]), show_col_types = FALSE)
  obs <- readr::read_csv(path_in(input_files_required[[2]]), show_col_types = FALSE)

  required_dep <- c("deploymentID", "locationID", "latitude", "longitude",
                    "deploymentStart", "deploymentEnd")
  required_obs <- c("deploymentID", "eventID", "scientificName", "eventStart")
  stop_missing_columns(dep, required_dep, "[wolf] deployments")
  stop_missing_columns(obs, required_obs, "[wolf] observations")

  deployments <- dep %>%
    transmute(
      deploymentID = na_if(as.character(deploymentID), ""),
      plotID = na_if(as.character(locationID), ""),
      latitude = as.numeric(latitude),
      longitude = as.numeric(longitude),
      start = parse_time(deploymentStart),
      end = parse_time(deploymentEnd),
      deploymentEffort = as.numeric(difftime(end, start, units = "days")),
      month = format(start, "%Y-%m", tz = "UTC")
    ) %>%
    filter(!is.na(deploymentID),
           !is.na(plotID),
           is.finite(latitude),
           is.finite(longitude),
           !is.na(start),
           !is.na(end),
           is.finite(deploymentEffort),
           deploymentEffort > 0,
           !is.na(month),
           nzchar(month))

  if (!nrow(deployments)) stop("[wolf] no valid dated deployments.")

  wolf_events <- obs %>%
    transmute(
      deploymentID = na_if(as.character(deploymentID), ""),
      eventID = na_if(as.character(eventID), ""),
      scientificName = as.character(scientificName),
      event_start = parse_time(eventStart)
    ) %>%
    filter(scientificName %in% WOLF_NAMES,
           !is.na(deploymentID),
           !is.na(eventID)) %>%
    distinct(deploymentID, eventID, .keep_all = TRUE)

  missing_event_time <- sum(is.na(wolf_events$event_start))
  if (missing_event_time > 0) {
    stop("[wolf] cannot split events by month: ",
         missing_event_time, " wolf event(s) have missing eventStart.")
  }

  model_dat <- split_deployments_by_month(deployments)
  model_dat$wolf_events <- count_segment_wolf_events(model_dat, wolf_events)
  model_dat <- model_dat %>%
    mutate(
      total_effort_days = deploymentEffort,
      wolf_events_per_100_days = 100 * wolf_events / total_effort_days
    ) %>%
    arrange(plotID, start)

  model_dat <- add_month_design(model_dat, settings)

  camera_rate <- camera_summary_from_model(model_dat)

  month_summary <- model_dat %>%
    group_by(month) %>%
    summarise(
      model_rows = n(),
      source_deployments = n_distinct(deploymentID),
      cameras = n_distinct(plotID),
      events = sum(wolf_events),
      effort_days = sum(total_effort_days),
      rate_per_100 = 100 * events / effort_days,
      .groups = "drop"
    )

  readr::write_csv(model_dat,
                   path_out(paste0(SURVEY_PREFIX, "_deployment_month_effort_rates.csv")))
  readr::write_csv(camera_rate,
                   path_out(paste0(SURVEY_PREFIX, "_camera_effort_rates.csv")))
  readr::write_csv(month_summary,
                   path_out(paste0(SURVEY_PREFIX, "_month_observed_summary.csv")))

  cat(sprintf(
    "[wolf] cameras %d | model rows %d | positive rows %d | events %d | effort %.1f camera-days | observed %.2f /100\n",
    nrow(camera_rate),
    nrow(model_dat),
    sum(model_dat$wolf_events > 0),
    sum(model_dat$wolf_events),
    sum(model_dat$total_effort_days),
    100 * sum(model_dat$wolf_events) / sum(model_dat$total_effort_days)
  ))
  cat(sprintf(
    "[wolf] month effect: reference=%s | prediction=%s | months=%s\n",
    unique(model_dat$month_reference),
    unique(model_dat$month_prediction),
    paste(sort(unique(model_dat$month)), collapse = ", ")
  ))

  model_dat
}


# Forest-camera flat-file loader --------------------------------------------
# The forest survey arrives as a single flat file (one row per detection event
# with its deployment's dates/coordinates) rather than separate deployment and
# observation tables. This loader rebuilds the same camera-month `model_dat`
# contract the road loader produces -- including per-row `start`/`end`
# timestamps -- so every downstream analysis in this library is shared.

# Find the first existing candidate input file (direct path or in the data dir).
resolve_input_file <- function(candidates, label) {
  candidates <- candidates[nzchar(candidates)]
  for (candidate in candidates) {
    if (file.exists(candidate)) return(candidate)
    data_path <- path_in(candidate)
    if (file.exists(data_path)) return(data_path)
  }
  stop("Could not find ", label, ". Checked: ",
       paste(candidates, collapse = ", "),
       ". Put the file in WOLF_DATA_DIR or set WOLF_FOREST_FILE.")
}

# Split flat-file deployments into camera-month rows with per-month effort and timestamps.
split_deployment_month_effort <- function(deployments) {
  rows <- vector("list", nrow(deployments))

  for (i in seq_len(nrow(deployments))) {
    start <- deployments$start_date[[i]]
    end <- deployments$end_date[[i]]
    month_starts <- seq(as.Date(format(start, "%Y-%m-01")),
                        as.Date(format(end - 1, "%Y-%m-01")),
                        by = "1 month")

    rows[[i]] <- do.call(rbind, lapply(month_starts, function(month_start) {
      next_month <- seq(month_start, by = "1 month", length.out = 2)[[2]]
      overlap_start <- max(start, month_start)
      overlap_end <- min(end, next_month)
      effort <- as.numeric(overlap_end - overlap_start)

      data.frame(
        deploymentID = deployments$deploymentID[[i]],
        plotID = deployments$plotID[[i]],
        longitude = deployments$longitude[[i]],
        latitude = deployments$latitude[[i]],
        month = format(month_start, "%Y-%m"),
        # Per-row timestamps for temporal ordering / diagnostics, so forest rows
        # are fully compatible with the road pipeline's start-based checks.
        start = as.POSIXct(paste0(format(overlap_start), " 00:00:00"), tz = "UTC"),
        end = as.POSIXct(paste0(format(overlap_end), " 00:00:00"), tz = "UTC"),
        total_effort_days = effort,
        stringsAsFactors = FALSE
      )
    }))
  }

  dplyr::bind_rows(rows) %>%
    dplyr::filter(is.finite(total_effort_days), total_effort_days > 0)
}

# Load and assemble forest-camera camera-month rows from the single flat file.
load_forest_flat_survey <- function(settings) {
  input_file <- resolve_input_file(
    FOREST_INPUT_FILES,
    paste0("[", SURVEY_PREFIX, "] forest camera-trap input"))
  dat <- readr::read_csv(input_file, show_col_types = FALSE)
  required <- c("deploymentID", "eventID", "eventStart", "scientificName",
                "plotID", "deploymentEffort", "latitude", "longitude",
                "startDate", "endDate")
  stop_missing_columns(dat, required, paste0("[", SURVEY_PREFIX, "] flat input"))

  deployments <- dat %>%
    dplyr::transmute(
      deploymentID = dplyr::na_if(as.character(deploymentID), ""),
      plotID = dplyr::na_if(as.character(plotID), ""),
      latitude = as.numeric(latitude),
      longitude = as.numeric(longitude),
      start_date = as.Date(startDate),
      end_date = as.Date(endDate)
    ) %>%
    dplyr::filter(!is.na(deploymentID),
                  !is.na(plotID),
                  is.finite(latitude),
                  is.finite(longitude),
                  !is.na(start_date),
                  !is.na(end_date),
                  end_date > start_date) %>%
    dplyr::distinct(deploymentID, .keep_all = TRUE)

  if (!nrow(deployments)) stop("[", SURVEY_PREFIX, "] no valid dated deployments.")

  month_effort <- split_deployment_month_effort(deployments)

  wolf_events <- dat %>%
    dplyr::transmute(
      deploymentID = dplyr::na_if(as.character(deploymentID), ""),
      eventID = dplyr::na_if(as.character(eventID), ""),
      scientificName = as.character(scientificName),
      event_time = parse_time(eventStart),
      event_month = format(event_time, "%Y-%m", tz = "UTC")
    ) %>%
    dplyr::filter(scientificName %in% WOLF_NAMES,
                  !is.na(deploymentID),
                  !is.na(eventID),
                  !is.na(event_month),
                  nzchar(event_month)) %>%
    dplyr::distinct(deploymentID, eventID, .keep_all = TRUE) %>%
    dplyr::count(deploymentID, month = event_month, name = "wolf_events")

  model_dat <- month_effort %>%
    dplyr::left_join(wolf_events, by = c("deploymentID", "month")) %>%
    dplyr::mutate(
      wolf_events = tidyr::replace_na(wolf_events, 0L),
      wolf_events_per_100_days = 100 * wolf_events / total_effort_days
    ) %>%
    add_month_design(settings) %>%
    dplyr::arrange(plotID, deploymentID, month)

  camera_rate <- camera_summary_from_model(model_dat)

  month_summary <- model_dat %>%
    dplyr::group_by(month) %>%
    dplyr::summarise(
      model_rows = dplyr::n(),
      source_deployments = dplyr::n_distinct(deploymentID),
      cameras = dplyr::n_distinct(plotID),
      events = sum(wolf_events),
      effort_days = sum(total_effort_days),
      rate_per_100 = 100 * events / effort_days,
      .groups = "drop"
    )

  readr::write_csv(model_dat,
                   path_out(paste0(SURVEY_PREFIX, "_deployment_month_effort_rates.csv")))
  readr::write_csv(camera_rate,
                   path_out(paste0(SURVEY_PREFIX, "_camera_effort_rates.csv")))
  readr::write_csv(month_summary,
                   path_out(paste0(SURVEY_PREFIX, "_month_observed_summary.csv")))

  cat(sprintf(
    "[wolf] cameras %d | model rows %d | positive rows %d | events %d | effort %.1f camera-days | observed %.2f /100\n",
    nrow(camera_rate),
    nrow(model_dat),
    sum(model_dat$wolf_events > 0),
    sum(model_dat$wolf_events),
    sum(model_dat$total_effort_days),
    100 * sum(model_dat$wolf_events) / sum(model_dat$total_effort_days)
  ))
  cat(sprintf(
    "[wolf] month effect: reference=%s | prediction=%s | months=%s\n",
    unique(model_dat$month_reference),
    unique(model_dat$month_prediction),
    paste(sort(unique(model_dat$month)), collapse = ", ")
  ))

  model_dat
}

# Dispatch on the runner-declared data shape. Both loaders return the same
# camera-month `model_dat` contract, so all analysis downstream is identical.
load_survey_data <- function(settings) {
  if (identical(SURVEY_DATA_SHAPE, "forest_flat")) {
    load_forest_flat_survey(settings)
  } else if (identical(SURVEY_DATA_SHAPE, "road_deployments")) {
    load_road_survey(settings)
  } else {
    stop("Unknown SURVEY_DATA_SHAPE: ", SURVEY_DATA_SHAPE,
         " (expected 'road_deployments' or 'forest_flat').")
  }
}


## 06. Spatial domain, mesh, and grid ----------------------------------------

# Convert camera lon/lat to an sf point layer in the study UTM CRS.
camera_to_utm <- function(camera_rate) {
  camera_rate %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(EPSG_UTM)
}

# Build the SPDE mesh and Matern spatial model (with PC priors) over the cameras.
build_spatial <- function(coords, settings) {
  range_arg <- if (!is.null(settings$fix_range_m)) {
    c(settings$fix_range_m, NA)
  } else {
    settings$prior_range_m
  }

  mesh <- INLA::inla.mesh.2d(
    loc = coords,
    cutoff = settings$mesh_cutoff_m,
    max.edge = settings$mesh_max_edge,
    offset = settings$mesh_offset
  )

  spde <- INLA::inla.spde2.pcmatern(
    mesh,
    alpha = 2,
    prior.range = range_arg,
    prior.sigma = settings$prior_sigma
  )

  list(
    mesh = mesh,
    spde = spde,
    s_index = INLA::inla.spde.make.index("spatial", n.spde = spde$n.spde),
    range_fixed = !is.null(settings$fix_range_m)
  )
}

# Build the prediction grid over the buffered convex hull of the cameras.
prediction_grid <- function(camera_sf, settings) {
  pts <- st_union(st_geometry(camera_sf))

  # Full prediction map: buffered convex hull around all cameras.
  area <- st_buffer(st_convex_hull(pts), settings$pred_buffer_m)

  grid <- st_make_grid(area, cellsize = settings$cell_size_m, what = "centers")
  pred_sf <- st_sf(grid_id = seq_along(grid), geometry = grid, crs = st_crs(camera_sf))
  pred_sf <- pred_sf[lengths(st_intersects(pred_sf, area)) > 0, ]


  if (!nrow(pred_sf)) stop("Prediction grid is empty.")
  pred_sf
}

# Assemble the fixed-effect design row used for map prediction.
prediction_fixed_effects <- function(model_dat, fixed_terms, settings, n_pred) {
  fixed_pred <- as.data.frame(matrix(0, nrow = n_pred, ncol = length(fixed_terms)))
  names(fixed_pred) <- fixed_terms

  if ("intercept" %in% fixed_terms) fixed_pred$intercept <- 1

  month_terms <- intersect(temporal_month_terms(model_dat), fixed_terms)
  if (length(month_terms)) {
    prediction_month <- settings$month_prediction
    prediction_term <- month_term_name(prediction_month)
    if (prediction_term %in% month_terms) fixed_pred[[prediction_term]] <- 1
  }

  fixed_pred
}


## 07. Posterior sample extraction -------------------------------------------

# Both R's RNG and INLA/GMRFLib's RNG are controlled. Never silently switch to
# internal-scale hyperparameters or replace failed joint draws with point values.
posterior_samples_safe <- function(fit, nsim, selection = list(),
                                   seed = INLA_SAMPLE_SEED) {
  nsim <- as.integer(nsim)
  samples <- tryCatch(
    with_validation_seed(seed, INLA::inla.posterior.sample(
      nsim, fit, intern = FALSE, selection = selection,
      seed = as.integer(seed), num.threads = 1L)),
    error = function(e) stop("Joint posterior sampling failed: ", conditionMessage(e))
  )
  if (is.null(samples) || !length(samples)) {
    stop("Posterior sampling failed. Make sure control.compute$config = TRUE.")
  }
  samples
}

# Extract the linear-predictor rows from one posterior sample.
predictor_rows_in_sample <- function(sample) {
  rn <- rownames(sample$latent)
  rows <- grep("^APredictor", rn)
  if (!length(rows)) rows <- grep("^Predictor", rn)
  rows
}

# Assemble the linear-predictor (eta) matrix across posterior samples for the stack rows.
extract_eta_matrix <- function(samples, stack_index, expected_n_stack = NULL) {
  first_rows <- predictor_rows_in_sample(samples[[1]])
  if (!length(first_rows)) stop("Could not find APredictor/Predictor rows in posterior samples.")

  if (!is.null(expected_n_stack) && length(first_rows) < max(stack_index)) {
    stop("Posterior sample predictor length is shorter than requested stack index.")
  }

  out <- vapply(samples, function(s) {
    rows <- predictor_rows_in_sample(s)
    as.numeric(s$latent[rows[stack_index], 1])
  }, numeric(length(stack_index)))

  if (is.null(dim(out))) out <- matrix(out, ncol = 1)
  out
}

# Extract a hyperparameter vector across posterior samples, with a fallback value.
sample_hyper_vector <- function(samples, pattern, fallback, transform = identity) {
  vapply(samples, function(s) {
    hp <- s$hyperpar
    if (!is.null(hp) && length(hp)) {
      i <- grep(pattern, names(hp), ignore.case = TRUE)
      if (length(i)) {
        val <- as.numeric(hp[i[[1]]])
        val <- transform(val)
        if (is.finite(val)) return(val)
      }
    }
    fallback
  }, numeric(1))
}

# Extract zero-inflation probability draws across posterior samples.
extract_pi_samples <- function(fit, samples, family) {
  if (!is_zi(family)) return(rep(0, length(samples)))
  p <- user_hyperparameter_draws(samples, PAT_ZPROB, "zero-inflation probability")
  if (any(p < 0 | p >= 1)) stop("Zero-inflation draws must be on user scale in [0, 1).")
  p
}

# Extract NB size draws across posterior samples.
extract_size_samples <- function(fit, samples, family) {
  if (!is_nb(family)) return(rep(NA_real_, length(samples)))
  s <- user_hyperparameter_draws(samples, PAT_NB_SIZE, "negative-binomial size")
  if (any(s <= 0)) stop("Negative-binomial size draws must be positive user-scale values.")
  fam_nb_size(s)
}

# Build per-row posterior draws of mean, variance, and fitted counts from samples.
build_posterior_draws <- function(fit, samples, index, effort, family,
                                  expected_n_stack = NULL) {
  eta <- extract_eta_matrix(samples, index, expected_n_stack)
  pi <- extract_pi_samples(fit, samples, family)
  size <- extract_size_samples(fit, samples, family)

  nsim <- ncol(eta)
  mu <- eta
  for (j in seq_len(nsim)) {
    mu[, j] <- effort * exp(eta[, j])
  }

  fitted <- mu
  fit_var <- mu
  for (j in seq_len(nsim)) {
    fitted[, j] <- fam_mean(mu[, j], pi[j], family, size[j])
    fit_var[, j] <- fam_var(mu[, j], pi[j], family, size[j])
  }

  list(eta = eta, mu = mu, fitted = fitted, fit_var = fit_var,
       pi = pi, size = size)
}

# Simulate replicate counts from posterior draws.
simulate_from_draws <- function(draws, family) {
  nsim <- ncol(draws$mu)
  sim <- draws$mu
  for (j in seq_len(nsim)) {
    sim[, j] <- fam_sim(draws$mu[, j], draws$pi[j], family, draws$size[j])
  }
  storage.mode(sim) <- "integer"
  sim
}


## 08. PPC, PIT, residual diagnostics ----------------------------------------

# Sum the rows of a draw matrix within groups (e.g. by camera).
aggregate_matrix_by_group <- function(mat, group) {
  group <- as.factor(group)
  apply(mat, 2, function(x) as.numeric(rowsum(x, group)[, 1]))
}

# Summarise posterior-predictive simulations against observed totals, zeros, max, and PIT.
summarise_ppc_simulations <- function(sim, model_dat, method) {
  yobs <- model_dat$y
  camera_group <- as.factor(model_dat$plotID)
  yobs_camera <- as.numeric(rowsum(yobs, camera_group)[, 1])

  if (is.null(dim(sim))) sim <- matrix(sim, ncol = 1)
  nsim <- ncol(sim)

  sim_camera <- aggregate_matrix_by_group(sim, camera_group)
  if (is.null(dim(sim_camera))) sim_camera <- matrix(sim_camera, ncol = 1)

  row_total <- colSums(sim)
  row_zero_fraction <- colMeans(sim == 0)
  row_max <- apply(sim, 2, max)

  camera_total <- colSums(sim_camera)
  camera_zero_fraction <- colMeans(sim_camera == 0)
  camera_max <- apply(sim_camera, 2, max)

  ppc_stat <- function(stat, observed, values, level) {
    data.frame(
      level = level,
      stat = stat,
      observed = observed,
      sim_median = median(values),
      sim_q025 = unname(quantile(values, 0.025)),
      sim_q975 = unname(quantile(values, 0.975)),
      pass = observed >= unname(quantile(values, 0.025)) &
        observed <= unname(quantile(values, 0.975)),
      method = method
    )
  }

  summary <- bind_rows(
    ppc_stat("total_events", sum(yobs), row_total, "model_row"),
    ppc_stat("zero_fraction", mean(yobs == 0), row_zero_fraction, "model_row"),
    ppc_stat("max_count", max(yobs), row_max, "model_row"),
    ppc_stat("total_events", sum(yobs_camera), camera_total, "camera"),
    ppc_stat("zero_fraction", mean(yobs_camera == 0), camera_zero_fraction, "camera"),
    ppc_stat("max_count", max(yobs_camera), camera_max, "camera")
  )

  row_pit <- (rowSums(sim < yobs) + runif(length(yobs)) * rowSums(sim == yobs)) / nsim
  cam_pit <- (rowSums(sim_camera < yobs_camera) +
                runif(length(yobs_camera)) * rowSums(sim_camera == yobs_camera)) / nsim

  list(summary = summary, row_pit = row_pit, camera_pit = cam_pit,
       sim = sim, sim_camera = sim_camera, nsim = nsim)
}

# Kolmogorov-Smirnov p-value testing PIT values against uniformity.
ks_uniform_p_value <- function(pit) {
  pit <- pit[is.finite(pit)]
  if (length(pit) < 5) return(NA_real_)
  suppressWarnings(tryCatch(
    ks.test(pit, "punif")$p.value,
    error = function(e) NA_real_
  ))
}

# Two-sided permutation test of residual spatial autocorrelation (Moran I).
moran_perm <- function(coords, x, nperm = MORAN_NPERM, two_sided = TRUE) {
  ok <- is.finite(x)
  coords <- coords[ok, , drop = FALSE]
  x <- x[ok]
  n <- length(x)

  if (n < 5) {
    return(list(I = NA_real_, expected = NA_real_, p_value = NA_real_,
                alternative = if (two_sided) "two_sided" else "greater"))
  }

  D <- as.matrix(dist(coords))
  W <- 1 / D
  diag(W) <- 0
  W[!is.finite(W)] <- 0

  rs <- rowSums(W)
  rs[rs == 0] <- 1
  W <- W / rs

  z <- x - mean(x)
  z2 <- sum(z^2)
  S0 <- sum(W)
  expected <- -1 / (n - 1)
  if (!is.finite(z2) || z2 <= 0 || !is.finite(S0) || S0 <= 0) {
    return(list(I = NA_real_, expected = expected, p_value = NA_real_,
                alternative = if (two_sided) "two_sided" else "greater"))
  }

  I_stat <- function(v) {
    (n / S0) * sum(W * outer(v, v)) / sum(v^2)
  }

  I0 <- I_stat(z)
  perm <- replicate(nperm, I_stat(sample(z)))

  p_value <- if (two_sided) {
    obs_dev <- abs(I0 - expected)
    perm_dev <- abs(perm - expected)
    (1 + sum(perm_dev >= obs_dev)) / (nperm + 1)
  } else {
    (1 + sum(perm >= I0)) / (nperm + 1)
  }

  list(I = I0, expected = expected, p_value = p_value,
       alternative = if (two_sided) "two_sided" else "greater")
}

# Empirical residual semivariogram by distance bin.
resid_variogram <- function(coords, residual, nbins = 12L) {
  ok <- is.finite(residual) & is.finite(coords[, 1]) & is.finite(coords[, 2])
  coords <- coords[ok, , drop = FALSE]
  residual <- residual[ok]
  if (length(residual) < 3) {
    return(data.frame(dist = numeric(), gamma = numeric(), n = integer()))
  }

  D <- as.matrix(dist(coords))
  G <- as.matrix(dist(residual))^2 / 2

  d <- D[lower.tri(D)]
  g <- G[lower.tri(G)]
  if (!length(d)) return(data.frame(dist = numeric(), gamma = numeric(), n = integer()))

  br <- seq(0, quantile(d, 0.9, na.rm = TRUE), length.out = nbins + 1)
  bin <- cut(d, br, include.lowest = TRUE)

  data.frame(
    dist = as.numeric(tapply(d, bin, mean)),
    gamma = as.numeric(tapply(g, bin, mean)),
    n = as.integer(table(bin))
  ) %>% filter(is.finite(dist), is.finite(gamma), n > 0)
}

# Compute camera-level residual diagnostics.
camera_residual_diagnostics <- function(model_dat) {
  model_dat %>%
    group_by(plotID) %>%
    summarise(
      longitude = mean(longitude, na.rm = TRUE),
      latitude = mean(latitude, na.rm = TRUE),
      y = sum(y, na.rm = TRUE),
      fitted_count = sum(fitted_count, na.rm = TRUE),
      fit_var = sum(fit_var, na.rm = TRUE),
      total_effort_days = sum(total_effort_days, na.rm = TRUE),
      wolf_events_per_100_days = 100 * y / total_effort_days,
      fitted_rate_per_100 = 100 * fitted_count / total_effort_days,
      n_model_rows = n(),
      month_first = min(month, na.rm = TRUE),
      month_last = max(month, na.rm = TRUE),
      pearson = (y - fitted_count) / sqrt(pmax(fit_var, 1e-9)),
      .groups = "drop"
    ) %>%
    arrange(plotID)
}

# Decide pass/fail for a posterior-predictive check from its tail probability.
ppc_pass_lookup <- function(ppc_summary, level, stat) {
  x <- ppc_summary$pass[ppc_summary$level == level & ppc_summary$stat == stat]
  if (length(x)) isTRUE(x[[1]]) else FALSE
}

# Run the full diagnostic suite (PPC, PIT, dispersion, Moran) and write its files.
compute_diagnostics <- function(fit, samples, model_dat, obs_index, camera_sf,
                                family, write_files = TRUE) {
  obs_draws <- build_posterior_draws(
    fit = fit,
    samples = samples,
    index = obs_index,
    effort = model_dat$total_effort_days,
    family = family,
    expected_n_stack = length(fit$summary.linear.predictor$mean)
  )

  sim <- simulate_from_draws(obs_draws, family)
  ppc <- summarise_ppc_simulations(sim, model_dat, "joint_posterior")

  model_dat$fitted_count <- rowMeans(obs_draws$fitted)
  model_dat$fit_var <- rowMeans(obs_draws$fit_var) +
    apply(obs_draws$fitted, 1, var, na.rm = TRUE)
  model_dat$eta_mean <- rowMeans(obs_draws$eta)
  model_dat$eta_sd <- apply(obs_draws$eta, 1, sd, na.rm = TRUE)
  model_dat$pearson <- (model_dat$y - model_dat$fitted_count) /
    sqrt(pmax(model_dat$fit_var, 1e-9))

  if (!is.null(fit$summary.fitted.values) &&
      nrow(fit$summary.fitted.values) >= max(obs_index)) {
    model_dat$inla_fitted_mean <- fit$summary.fitted.values$mean[obs_index]
  } else {
    model_dat$inla_fitted_mean <- NA_real_
  }

  camera_diag <- camera_residual_diagnostics(model_dat)
  # Aggregate conditional moments within each JOINT draw, then apply total
  # variance. Summing marginal row variances loses posterior covariance across
  # months, especially with camera IID or AR1 effects.
  camera_moments <- aggregate_joint_predictive_moments(
    obs_draws$fitted, obs_draws$fit_var, as.character(model_dat$plotID))
  camera_match <- match(as.character(camera_diag$plotID), camera_moments$plotID)
  if (anyNA(camera_match)) stop("Camera predictive moments are not aligned.")
  camera_diag$fitted_count <- camera_moments$fitted_count[camera_match]
  camera_diag$fit_var <- camera_moments$fit_var[camera_match]
  camera_diag$pearson <- (camera_diag$y - camera_diag$fitted_count) /
    sqrt(pmax(camera_diag$fit_var, 1e-9))
  coords_camera <- st_coordinates(camera_sf)
  moran <- moran_perm(coords_camera, camera_diag$pearson,
                      nperm = MORAN_NPERM, two_sided = TRUE)

  fitted_check <- data.frame(
    check = c("observed_total", "posterior_sample_fitted_total", "inla_fitted_total",
              "cor_posterior_fitted_vs_inla_fitted"),
    value = c(
      sum(model_dat$y),
      sum(model_dat$fitted_count),
      sum(model_dat$inla_fitted_mean, na.rm = TRUE),
      safe_cor(model_dat$fitted_count, model_dat$inla_fitted_mean)
    )
  )

  row_disp <- mean(model_dat$pearson^2, na.rm = TRUE)
  cam_disp <- mean(camera_diag$pearson^2, na.rm = TRUE)

  # Required-check gate: posterior predictive checks (total events, zero
  # fraction, max count) plus residual spatial autocorrelation (Moran's I).
  # PIT KS (ppc_pit_ks_row / ppc_pit_ks_camera, printed below) is computed
  # and reported as supporting evidence of calibration but deliberately does
  # NOT gate diagnostics_ok: PIT KS is sensitive to camera-level clustering
  # and small-count discreteness in a way that does not track whether the
  # mapped spatial surface is distorted, so a low PIT KS p-value alone is
  # not treated as disqualifying. It is never silently dropped -- it is
  # printed and logged alongside the gated checks.
  ppc_total_pass <- ppc_pass_lookup(ppc$summary, "camera", "total_events")
  ppc_zero_pass <- ppc_pass_lookup(ppc$summary, "camera", "zero_fraction")
  ppc_max_pass <- ppc_pass_lookup(ppc$summary, "camera", "max_count")
  moran_pass <- is.finite(moran$p_value) && moran$p_value >= MORAN_ALPHA

  diagnostics <- list(
    model_dat = model_dat,
    camera_diag = camera_diag,
    pi_hat = mean(obs_draws$pi, na.rm = TRUE),
    size_hat = mean(obs_draws$size, na.rm = TRUE),
    pearson_disp = row_disp,
    pearson_disp_camera = cam_disp,
    moran_I = moran$I,
    moran_expected = moran$expected,
    moran_p = moran$p_value,
    moran_alternative = moran$alternative,
    ppc_pit_ks_row = ks_uniform_p_value(ppc$row_pit),
    ppc_pit_ks_camera = ks_uniform_p_value(ppc$camera_pit),
    pit_mean_row = mean(ppc$row_pit, na.rm = TRUE),
    pit_mean_camera = mean(ppc$camera_pit, na.rm = TRUE),
    ppc_method = "joint_posterior",
    ppc_nsim = ppc$nsim,
    ppc = ppc$summary,
    row_pit = ppc$row_pit,
    camera_pit = ppc$camera_pit,
    ppc_total_pass = ppc_total_pass,
    ppc_zero_pass = ppc_zero_pass,
    ppc_max_pass = ppc_max_pass,
    moran_pass = moran_pass,
    diagnostics_ok = ppc_total_pass && ppc_zero_pass && ppc_max_pass && moran_pass,
    fitted_check = fitted_check,
    obs_draws = obs_draws
  )

  if (write_files) {
    readr::write_csv(model_dat,
                     path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                                     "_model_row_diagnostics.csv")))
    readr::write_csv(camera_diag,
                     path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                                     "_camera_residual_diagnostics.csv")))
    readr::write_csv(ppc$summary,
                     path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                                     "_posterior_predictive_check.csv")))
    readr::write_csv(fitted_check,
                     path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                                     "_fitted_scale_sanity_check.csv")))
  }

  diagnostics
}


## 09. Temporal autocorrelation diagnostics ----------------------------------

# These diagnostics test whether residuals remain temporally structured after
# the model has accounted for effort, spatial structure, and calendar-month
# effects. They are diagnostic checks only; the fitted model remains the same.
#
# Outputs:
#   * residuals versus deployment start date
#   * pooled within-camera residual lag correlations, including lag 1
#   * ACF of mean residuals by deployment date
#   * a text note explaining limitations when there are too few repeated rows

# Compute and write residual temporal-autocorrelation diagnostics.
temporal_autocorrelation_diagnostics <- function(model_dat,
                                                 prefix = SURVEY_PREFIX,
                                                 final_model = FINAL_MODEL_NAME,
                                                 max_lag = 5L,
                                                 min_pairs = 5L) {
  required <- c("plotID", "start", "month", "pearson")
  missing <- setdiff(required, names(model_dat))
  note_file <- path_out(paste0(prefix, "_", final_model,
                               "_temporal_autocorrelation_NOTE.txt"))

  if (length(missing)) {
    writeLines(
      c(
        "Temporal autocorrelation diagnostics skipped.",
        sprintf("Missing required column(s): %s", paste(missing, collapse = ", "))
      ),
      note_file
    )
    return(NULL)
  }

  dat <- model_dat %>%
    mutate(
      start_date = as_utc_date(start),
      start_doy = as.integer(format(start_date, "%j"))
    ) %>%
    filter(!is.na(plotID), !is.na(start_date), is.finite(pearson)) %>%
    arrange(plotID, start)

  if (!nrow(dat)) {
    writeLines(
      c("Temporal autocorrelation diagnostics skipped.",
        "No rows with finite Pearson residuals and valid deployment dates."),
      note_file
    )
    return(NULL)
  }

  readr::write_csv(
    dat,
    path_out(paste0(prefix, "_", final_model,
                    "_temporal_residual_input.csv"))
  )

  # --------------------------------------------------------------------------
  # 1. Residuals through calendar time.
  # --------------------------------------------------------------------------

  p_time <- ggplot(dat, aes(start_date, pearson)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_point(alpha = 0.75) +
    geom_smooth(se = TRUE, method = "loess", formula = y ~ x) +
    labs(
      title = paste0("Temporal residual trend: ", prefix),
      subtitle = paste0(final_model, " | residuals by deployment start date"),
      x = "deployment start date",
      y = "row Pearson residual"
    ) +
    theme_minimal(base_size = 12)

  ggsave(
    path_out(paste0(prefix, "_", final_model,
                    "_diag_temporal_residuals_vs_date.png")),
    p_time,
    width = 7,
    height = 5,
    dpi = 220
  )

  # --------------------------------------------------------------------------
  # 2. Pooled within-camera lag residual correlations.
  #    Lag is deployment sequence at the same camera, not a fixed time interval.
  # --------------------------------------------------------------------------

  make_lag_pairs <- function(lag_k) {
    pieces <- lapply(split(dat, dat$plotID), function(d) {
      d <- d[order(d$start), , drop = FALSE]
      n <- nrow(d)
      if (n <= lag_k) return(NULL)

      previous <- seq_len(n - lag_k)
      current <- previous + lag_k

      data.frame(
        plotID = d$plotID[current],
        lag = lag_k,
        start_previous = d$start[previous],
        start_current = d$start[current],
        residual_previous = d$pearson[previous],
        residual_current = d$pearson[current],
        days_between = as.numeric(difftime(d$start[current],
                                           d$start[previous],
                                           units = "days")),
        stringsAsFactors = FALSE
      )
    })
    dplyr::bind_rows(pieces)
  }

  lag_pairs <- dplyr::bind_rows(lapply(seq_len(max_lag), make_lag_pairs))

  if (nrow(lag_pairs)) {
    readr::write_csv(
      lag_pairs,
      path_out(paste0(prefix, "_", final_model,
                      "_temporal_within_camera_lag_pairs.csv"))
    )
  }

  safe_cor_p <- function(a, b) {
    safe_cor_p_value(a, b, min_pairs)
  }

  lag_summary <- if (nrow(lag_pairs)) {
    lag_pairs %>%
      group_by(lag) %>%
      summarise(
        n_pairs = n(),
        n_cameras = n_distinct(plotID),
        correlation = safe_cor(residual_previous, residual_current),
        p_value = safe_cor_p(residual_previous, residual_current),
        mean_days_between = mean(days_between, na.rm = TRUE),
        median_days_between = median(days_between, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    data.frame(
      lag = integer(),
      n_pairs = integer(),
      n_cameras = integer(),
      correlation = numeric(),
      p_value = numeric(),
      mean_days_between = numeric(),
      median_days_between = numeric()
    )
  }

  readr::write_csv(
    lag_summary,
    path_out(paste0(prefix, "_", final_model,
                    "_temporal_within_camera_lag_correlation.csv"))
  )

  lag1 <- lag_pairs %>% filter(lag == 1L)
  if (nrow(lag1) >= min_pairs &&
      is.finite(safe_cor(lag1$residual_previous, lag1$residual_current))) {
    lag1_row <- lag_summary %>% filter(lag == 1L)

    p_lag1 <- ggplot(lag1, aes(residual_previous, residual_current)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey70") +
      geom_point(alpha = 0.75) +
      geom_smooth(se = TRUE, method = "lm", formula = y ~ x) +
      labs(
        title = paste0("Within-camera lag-1 residual correlation: ", prefix),
        subtitle = sprintf(
          "%s | r = %.3f, p = %.3g, n = %d pairs, median gap = %.1f days",
          final_model,
          lag1_row$correlation[[1]],
          lag1_row$p_value[[1]],
          lag1_row$n_pairs[[1]],
          lag1_row$median_days_between[[1]]
        ),
        x = "previous residual at same camera",
        y = "current residual"
      ) +
      theme_minimal(base_size = 12)

    ggsave(
      path_out(paste0(prefix, "_", final_model,
                      "_diag_temporal_lag1_residual_correlation.png")),
      p_lag1,
      width = 6.5,
      height = 5.5,
      dpi = 220
    )
  } else {
    writeLines(
      c(
        "Within-camera lag-1 residual plot skipped.",
        sprintf("Usable lag-1 pairs: %d", nrow(lag1)),
        sprintf("Minimum required pairs: %d", min_pairs),
        "This usually means there are too few repeated camera-month rows per camera."
      ),
      path_out(paste0(prefix, "_", final_model,
                      "_temporal_lag1_residual_correlation_NOTE.txt"))
    )
  }

  if (nrow(lag_summary)) {
    p_lag_summary <- ggplot(lag_summary, aes(lag, correlation)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
      geom_point(aes(size = n_pairs)) +
      geom_line() +
      scale_x_continuous(breaks = seq_len(max_lag)) +
      labs(
        title = paste0("Within-camera residual autocorrelation by lag: ", prefix),
        subtitle = paste0(final_model, " | lag is deployment sequence within camera"),
        x = "within-camera deployment lag",
        y = "pooled residual correlation",
        size = "pairs"
      ) +
      theme_minimal(base_size = 12)

    ggsave(
      path_out(paste0(prefix, "_", final_model,
                      "_diag_temporal_lag_correlation_summary.png")),
      p_lag_summary,
      width = 6.5,
      height = 5,
      dpi = 220
    )
  }

  # --------------------------------------------------------------------------
  # 3. ACF of mean residuals by deployment start date.
  #    This is a broad date-ordered diagnostic. It mixes cameras, so it should
  #    be interpreted as evidence of survey-wide temporal residual structure,
  #    not camera-specific serial dependence.
  # --------------------------------------------------------------------------

  daily_resid <- dat %>%
    group_by(start_date) %>%
    summarise(
      mean_pearson = mean(pearson, na.rm = TRUE),
      median_pearson = median(pearson, na.rm = TRUE),
      n_rows = n(),
      .groups = "drop"
    ) %>%
    arrange(start_date)

  readr::write_csv(
    daily_resid,
    path_out(paste0(prefix, "_", final_model,
                    "_temporal_mean_residual_by_date.csv"))
  )

  acf_summary <- data.frame(lag = integer(), acf = numeric())
  if (nrow(daily_resid) >= 6 && sd(daily_resid$mean_pearson, na.rm = TRUE) > 0) {
    acf_obj <- stats::acf(daily_resid$mean_pearson,
                          plot = FALSE,
                          na.action = stats::na.pass)
    acf_summary <- data.frame(
      lag = as.integer(round(as.numeric(acf_obj$lag))),
      acf = as.numeric(acf_obj$acf)
    )
    readr::write_csv(
      acf_summary,
      path_out(paste0(prefix, "_", final_model,
                      "_temporal_acf_mean_residual_by_date.csv"))
    )

    png(
      filename = path_out(paste0(prefix, "_", final_model,
                                 "_diag_temporal_acf_mean_residual_by_date.png")),
      width = 1400,
      height = 1000,
      res = 180
    )
    stats::acf(
      daily_resid$mean_pearson,
      main = paste0("ACF of mean residuals by deployment date: ", prefix),
      xlab = "lag in ordered deployment dates"
    )
    dev.off()
  } else {
    writeLines(
      c(
        "Date-ordered residual ACF skipped.",
        sprintf("Unique deployment dates: %d", nrow(daily_resid)),
        "At least 6 unique dates and non-constant mean residuals are required."
      ),
      path_out(paste0(prefix, "_", final_model,
                      "_temporal_acf_mean_residual_by_date_NOTE.txt"))
    )
  }

  # --------------------------------------------------------------------------
  # 4. Compact text report.
  # --------------------------------------------------------------------------

  lag1_summary <- lag_summary %>% filter(lag == 1L)
  lag1_line <- if (nrow(lag1_summary)) {
    sprintf(
      "Within-camera lag-1 residual correlation: r = %.3f, p = %.4g, n pairs = %d, cameras = %d, median gap = %.1f days.",
      lag1_summary$correlation[[1]],
      lag1_summary$p_value[[1]],
      lag1_summary$n_pairs[[1]],
      lag1_summary$n_cameras[[1]],
      lag1_summary$median_days_between[[1]]
    )
  } else {
    "Within-camera lag-1 residual correlation was not evaluable."
  }

  acf1 <- if (nrow(acf_summary) >= 2 && any(acf_summary$lag == 1L)) {
    acf_summary$acf[acf_summary$lag == 1L][[1]]
  } else {
    NA_real_
  }
  acf_line <- if (is.finite(acf1)) {
    sprintf("Date-ordered ACF of mean residuals: lag-1 ACF = %.3f.", acf1)
  } else {
    "Date-ordered ACF of mean residuals was not evaluable."
  }

  report <- c(
    sprintf("Temporal autocorrelation diagnostics: %s", prefix),
    sprintf("Model: %s", final_model),
    sprintf("Rows used: %d", nrow(dat)),
    sprintf("Cameras used: %d", dplyr::n_distinct(dat$plotID)),
    sprintf("Unique deployment start dates: %d", dplyr::n_distinct(dat$start_date)),
    "",
    "Diagnostics are based on Pearson residuals from the fitted model.",
    "The month fixed effect already adjusts the mean by calendar camera-month; these diagnostics check residual temporal structure remaining after that adjustment.",
    "",
    lag1_line,
    acf_line
  )

  writeLines(
    report,
    path_out(paste0(prefix, "_", final_model,
                    "_temporal_autocorrelation_report.txt"))
  )

  list(
    data = dat,
    lag_pairs = lag_pairs,
    lag_summary = lag_summary,
    daily_resid = daily_resid,
    acf_summary = acf_summary,
    report = report
  )
}


## 10. Diagnostic plots -------------------------------------------------------

# Write the diagnostic figures (obs vs fitted, residual maps, PIT, variogram, etc.).
write_diagnostic_plots <- function(diag, camera_sf) {
  model_dat <- diag$model_dat
  camera_diag <- diag$camera_diag
  coords_camera <- st_coordinates(camera_sf)

  obs_fit_plot <- ggplot(camera_diag, aes(fitted_count, y)) +
    geom_point(alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    scale_x_continuous(trans = "sqrt") +
    scale_y_continuous(trans = "sqrt") +
    labs(
      title = "Observed vs fitted camera counts",
      subtitle = FINAL_MODEL_NAME,
      x = paste0("posterior mean fitted ", TARGET_LABEL, " events"),
      y = paste0("observed ", TARGET_LABEL, " events")
    ) +
    theme_minimal(base_size = 12)
  ggsave(path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                         "_diag_obs_vs_fitted.png")),
         obs_fit_plot, width = 6, height = 5.5, dpi = 220)

  residual_sf <- camera_sf
  residual_sf$pearson <- camera_diag$pearson
  residual_plot <- ggplot(residual_sf) +
    geom_sf(aes(colour = pearson, size = abs(pearson)), alpha = 0.9) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_size_continuous(range = c(2, 8)) +
    coord_sf(datum = NA) +
    labs(title = "Spatial Pearson residuals",
         subtitle = FINAL_MODEL_NAME,
         colour = "Pearson", size = "|Pearson|") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())
  ggsave(path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                         "_diag_spatial_residuals.png")),
         residual_plot, width = 7, height = 7, dpi = 220)

  pit_row_plot <- ggplot(data.frame(pit = diag$row_pit), aes(pit)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20,
                   fill = "grey55", colour = "white", boundary = 0) +
    geom_hline(yintercept = 1, linetype = 2, colour = "red") +
    labs(
      title = "Posterior predictive PIT, model rows",
      subtitle = sprintf("mean %.3f | KS p %.3g",
                         diag$pit_mean_row, diag$ppc_pit_ks_row),
      x = "PIT", y = "density"
    ) +
    xlim(0, 1) +
    theme_minimal(base_size = 12)
  ggsave(path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                         "_diag_pit_hist_model_rows.png")),
         pit_row_plot, width = 6, height = 5, dpi = 220)

  pit_camera_plot <- ggplot(data.frame(pit = diag$camera_pit), aes(pit)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20,
                   fill = "grey55", colour = "white", boundary = 0) +
    geom_hline(yintercept = 1, linetype = 2, colour = "red") +
    labs(
      title = "Posterior predictive PIT, camera aggregates",
      subtitle = sprintf("mean %.3f | KS p %.3g",
                         diag$pit_mean_camera, diag$ppc_pit_ks_camera),
      x = "PIT", y = "density"
    ) +
    xlim(0, 1) +
    theme_minimal(base_size = 12)
  ggsave(path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                         "_diag_pit_hist_camera.png")),
         pit_camera_plot, width = 6, height = 5, dpi = 220)

  variogram <- resid_variogram(coords_camera, camera_diag$pearson)
  readr::write_csv(variogram,
                   path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                                   "_diag_resid_variogram.csv")))
  if (nrow(variogram)) {
    variogram_plot <- ggplot(variogram, aes(dist, gamma)) +
      geom_point(aes(size = n)) +
      geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
      labs(
        title = "Residual semivariogram",
        subtitle = FINAL_MODEL_NAME,
        x = "distance (m)",
        y = "semivariance"
      ) +
      theme_minimal(base_size = 12)
    ggsave(path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                           "_diag_resid_variogram.png")),
           variogram_plot, width = 6.5, height = 5, dpi = 220)
  }

  residual_covariate <- model_dat %>%
    mutate(
      fitted_rate_per_100 = 100 * fitted_count / total_effort_days,
      northing_proxy = latitude
    )

  effort_plot <- ggplot(residual_covariate,
                        aes(total_effort_days, pearson)) +
    geom_point(alpha = 0.75) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
    labs(title = "Residuals vs effort",
         x = "deployment effort days", y = "row Pearson residual") +
    theme_minimal(base_size = 12)
  ggsave(path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                         "_diag_resid_vs_effort.png")),
         effort_plot, width = 6, height = 5, dpi = 220)

  month_plot <- ggplot(residual_covariate,
                       aes(month, pearson)) +
    geom_boxplot(outlier.alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(title = "Residuals by calendar month",
         x = "month", y = "row Pearson residual") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                         "_diag_resid_by_month.png")),
         month_plot, width = 7, height = 5, dpi = 220)

  lat_plot <- ggplot(residual_covariate,
                     aes(latitude, pearson)) +
    geom_point(alpha = 0.75) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
    labs(title = "Residuals vs latitude",
         x = "latitude", y = "row Pearson residual") +
    theme_minimal(base_size = 12)
  ggsave(path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                         "_diag_resid_vs_latitude.png")),
         lat_plot, width = 6, height = 5, dpi = 220)

  invisible(TRUE)
}


## 11. Prior-posterior plots and summaries -----------------------------------

# PC-prior density for the spatial range.
pc_range_density <- function(x, range0, prob_below_range0) {
  lambda <- -range0 * log(prob_below_range0)
  ifelse(x > 0, lambda * x^(-2) * exp(-lambda / x), 0)
}

# PC-prior quantile for the spatial range.
pc_range_quantile <- function(p, range0, prob_below_range0) {
  lambda <- -range0 * log(prob_below_range0)
  lambda / (-log(p))
}

# PC-prior density for the spatial marginal SD.
pc_sigma_density <- function(x, sigma0, prob_above_sigma0) {
  lambda <- -log(prob_above_sigma0) / sigma0
  ifelse(x >= 0, lambda * exp(-lambda * x), 0)
}

# PC-prior quantile for the spatial marginal SD.
pc_sigma_quantile <- function(p, sigma0, prob_above_sigma0) {
  lambda <- -log(prob_above_sigma0) / sigma0
  -log(1 - p) / lambda
}

# Prior density of the zero-inflation probability from the logit-Gaussian prior.
zip_prob_prior_density <- function(p) {
  sd_logit <- 1 / sqrt(PRIOR_ZI_LOGIT_PREC)
  dnorm(qlogis(p), mean = PRIOR_ZI_LOGIT_MEAN, sd = sd_logit) / (p * (1 - p))
}

# Prior density of the NB size from the log-Gaussian prior.
nb_size_prior_density <- function(x) {
  sd_logsize <- 1 / sqrt(PRIOR_NB_LOGSIZE_PREC)
  dnorm(log(x), mean = PRIOR_NB_LOGSIZE_MEAN, sd = sd_logsize) / x
}

# Prior quantile of the NB size.
nb_size_prior_quantile <- function(p) {
  sd_logsize <- 1 / sqrt(PRIOR_NB_LOGSIZE_PREC)
  exp(qnorm(p, mean = PRIOR_NB_LOGSIZE_MEAN, sd = sd_logsize))
}

# Plot one parameter prior against its posterior.
plot_prior_posterior_density <- function(parameter, prior_df,
                                         posterior_marginal, file_suffix,
                                         x_label, log_x = FALSE,
                                         reference_x = NULL,
                                         reference_label = NULL) {
  if (is.null(posterior_marginal) || nrow(posterior_marginal) < 2) {
    writeLines(
      paste("Posterior marginal not available for", parameter),
      path_out(paste0(SURVEY_PREFIX, "_prior_posterior_", file_suffix, "_NOTE.txt"))
    )
    return(invisible(NULL))
  }

  posterior_df <- data.frame(
    value = posterior_marginal[, 1],
    density = posterior_marginal[, 2],
    source = "posterior"
  ) %>% filter(is.finite(value), is.finite(density), density >= 0)

  prior_df <- prior_df %>%
    mutate(source = "prior") %>%
    filter(is.finite(value), is.finite(density), density >= 0)

  if (log_x) {
    posterior_df <- posterior_df %>% filter(value > 0)
    prior_df <- prior_df %>% filter(value > 0)
  }

  plot_df <- bind_rows(prior_df, posterior_df)
  if (!nrow(plot_df)) return(invisible(NULL))

  p <- ggplot(plot_df, aes(value, density, colour = source, linetype = source)) +
    geom_line(linewidth = 0.85) +
    labs(title = paste0("Prior vs posterior: ", parameter),
         subtitle = SURVEY_PREFIX,
         x = x_label, y = "density", colour = NULL, linetype = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")

  if (log_x) p <- p + scale_x_log10(labels = label_number())

  if (!is.null(reference_x) && is.finite(reference_x)) {
    p <- p + geom_vline(xintercept = reference_x, linetype = 2, colour = "grey35")
    if (!is.null(reference_label)) p <- p + labs(caption = reference_label)
  }

  ggsave(path_out(paste0(SURVEY_PREFIX, "_prior_posterior_", file_suffix, ".png")),
         p, width = 6.8, height = 4.8, dpi = 250)
  readr::write_csv(plot_df,
                   path_out(paste0(SURVEY_PREFIX, "_prior_posterior_", file_suffix, ".csv")))
  invisible(plot_df)
}

# Write all prior-versus-posterior overlay plots and tables.
write_prior_posterior_plots <- function(fit, settings, family) {
  intercept_post <- fit$marginals.fixed[["intercept"]]
  intercept_sd_prior <- 1 / sqrt(PRIOR_INTERCEPT_PREC)
  intercept_grid <- seq(qnorm(0.001, PRIOR_INTERCEPT_MEAN, intercept_sd_prior),
                        qnorm(0.999, PRIOR_INTERCEPT_MEAN, intercept_sd_prior),
                        length.out = 700)
  plot_prior_posterior_density(
    parameter = "intercept",
    prior_df = data.frame(value = intercept_grid,
                          density = dnorm(intercept_grid, PRIOR_INTERCEPT_MEAN, intercept_sd_prior)),
    posterior_marginal = intercept_post,
    file_suffix = "intercept",
    x_label = "log-rate intercept",
    reference_x = PRIOR_INTERCEPT_MEAN,
    reference_label = sprintf("Prior center: crude observed log daily rate = %.2f",
                              PRIOR_INTERCEPT_MEAN)
  )

  month_terms <- grep("^month_[0-9]{4}_[0-9]{2}$", names(fit$marginals.fixed), value = TRUE)
  if (length(month_terms)) {
    month_sd_prior <- 1 / sqrt(PRIOR_MONTH_LOG_RATE_RATIO_PREC)
    month_grid <- seq(qnorm(0.001, 0, month_sd_prior),
                      qnorm(0.999, 0, month_sd_prior),
                      length.out = 700)
    month_prior <- data.frame(value = month_grid,
                              density = dnorm(month_grid, 0, month_sd_prior))
    for (term in month_terms) {
      plot_prior_posterior_density(
        parameter = paste0("month log-rate ratio: ", month_from_term(term),
                           " vs ", settings$month_reference),
        prior_df = month_prior,
        posterior_marginal = fit$marginals.fixed[[term]],
        file_suffix = paste0("fixed_", term),
        x_label = "log-rate ratio"
      )
    }
  }

  if (!is.null(settings$fix_range_m)) {
    writeLines(c("Spatial range was fixed.",
                 sprintf("fix_range_m = %s", settings$fix_range_m)),
               path_out(paste0(SURVEY_PREFIX, "_prior_posterior_range_NOTE.txt")))
  } else {
    range0 <- settings$prior_range_m[1]
    prob_range <- settings$prior_range_m[2]
    range_grid <- exp(seq(log(pc_range_quantile(0.001, range0, prob_range)),
                          log(pc_range_quantile(0.995, range0, prob_range)),
                          length.out = 700))
    range_prior <- data.frame(value = range_grid,
                              density = pc_range_density(range_grid, range0, prob_range))
    plot_prior_posterior_density(
      parameter = "spatial range",
      prior_df = range_prior,
      posterior_marginal = hyp_marg(fit, PAT_RANGE),
      file_suffix = "range",
      x_label = "range (m)",
      log_x = TRUE,
      reference_x = range0,
      reference_label = sprintf("Prior statement: P(range < %.0f m) = %.2f",
                                range0, prob_range)
    )
  }

  sigma0 <- settings$prior_sigma[1]
  prob_sigma <- settings$prior_sigma[2]
  sigma_grid <- seq(0, pc_sigma_quantile(0.995, sigma0, prob_sigma), length.out = 700)
  sigma_prior <- data.frame(value = sigma_grid,
                            density = pc_sigma_density(sigma_grid, sigma0, prob_sigma))
  plot_prior_posterior_density(
    parameter = "spatial marginal SD",
    prior_df = sigma_prior,
    posterior_marginal = hyp_marg(fit, PAT_SIGMA),
    file_suffix = "spatial_sd",
    x_label = "spatial marginal SD",
    reference_x = sigma0,
    reference_label = sprintf("Prior statement: P(SD > %.2f) = %.2f",
                              sigma0, prob_sigma)
  )

  if (is_zi(family)) {
    prob_grid <- seq(0.001, 0.999, length.out = 700)
    plot_prior_posterior_density(
      parameter = "zero-inflation probability",
      prior_df = data.frame(value = prob_grid,
                            density = zip_prob_prior_density(prob_grid)),
      posterior_marginal = hyp_marg(fit, PAT_ZPROB),
      file_suffix = "zero_inflation_probability",
      x_label = "zero-inflation probability",
      reference_x = plogis(PRIOR_ZI_LOGIT_MEAN),
      reference_label = sprintf("Prior center on probability scale: %.2f",
                                plogis(PRIOR_ZI_LOGIT_MEAN))
    )
  }

  if (is_nb(family)) {
    size_grid <- exp(seq(log(nb_size_prior_quantile(0.001)),
                         log(nb_size_prior_quantile(0.995)),
                         length.out = 700))
    prior_median_size <- exp(PRIOR_NB_LOGSIZE_MEAN)
    plot_prior_posterior_density(
      parameter = "negative-binomial size",
      prior_df = data.frame(value = size_grid,
                            density = nb_size_prior_density(size_grid)),
      posterior_marginal = hyp_marg(fit, PAT_NB_SIZE),
      file_suffix = "nb_size",
      x_label = "negative-binomial size",
      log_x = TRUE,
      reference_x = prior_median_size,
      reference_label = sprintf(
        "Prior on log(size): Normal(%.2f, SD %.2f); median size %.2f",
        PRIOR_NB_LOGSIZE_MEAN,
        1 / sqrt(PRIOR_NB_LOGSIZE_PREC),
        prior_median_size
      )
    )
  }

  invisible(TRUE)
}



# Quantitative prior-influence screen ---------------------------------------
# These diagnostics do not prove causality by themselves. They summarize how
# much the posterior narrowed relative to the prior and whether the posterior
# lies in a low-prior-density region. Parameters flagged here are the ones that
# should be prioritized in the prior-sensitivity reruns.

# Trapezoidal numerical integration.
trapz_numeric <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2) return(NA_real_)
  o <- order(x)
  x <- x[o]
  y <- y[o]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

# Rescale a density so it integrates to one.
standardise_density <- function(df) {
  df <- df %>%
    filter(is.finite(value), is.finite(density), density >= 0) %>%
    arrange(value)
  if (nrow(df) < 2) return(df[FALSE, ])
  area <- trapz_numeric(df$value, df$density)
  if (!is.finite(area) || area <= 0) return(df[FALSE, ])
  df$density <- df$density / area
  df
}

# Transform a density onto the scale used for a prior-influence metric.
transform_density_for_metric <- function(df, transform = "identity") {
  df <- df %>% filter(is.finite(value), is.finite(density), density >= 0)
  if (identical(transform, "log")) {
    df <- df %>% filter(value > 0)
    if (!nrow(df)) return(data.frame(value = numeric(), density = numeric()))
    z <- log(df$value)
    # density for z = log(x): f_Z(z) = f_X(exp(z)) * exp(z)
    out <- data.frame(value = z, density = df$density * exp(z))
  } else {
    out <- data.frame(value = df$value, density = df$density)
  }
  standardise_density(out)
}

# Summary statistics (mean, SD, quantiles) of a density.
density_summary_stats <- function(df) {
  df <- standardise_density(df)
  if (nrow(df) < 2) {
    return(list(mean = NA_real_, sd = NA_real_, q025 = NA_real_,
                median = NA_real_, q975 = NA_real_))
  }
  x <- df$value
  d <- df$density
  dx <- diff(x)
  increments <- dx * (head(d, -1) + tail(d, -1)) / 2
  cdf <- c(0, cumsum(increments))
  if (max(cdf, na.rm = TRUE) > 0) cdf <- cdf / max(cdf, na.rm = TRUE)
  qfun <- function(p) {
    keep <- !duplicated(cdf) & is.finite(cdf) & is.finite(x)
    if (sum(keep) < 2) return(NA_real_)
    approx(cdf[keep], x[keep], xout = p, rule = 2)$y
  }
  m <- trapz_numeric(x, x * d)
  v <- trapz_numeric(x, (x - m)^2 * d)
  list(
    mean = m,
    sd = sqrt(pmax(v, 0)),
    q025 = qfun(0.025),
    median = qfun(0.5),
    q975 = qfun(0.975)
  )
}

# Overlap area between two densities (prior versus posterior).
density_overlap <- function(prior_df, post_df) {
  prior_df <- standardise_density(prior_df)
  post_df <- standardise_density(post_df)
  if (nrow(prior_df) < 2 || nrow(post_df) < 2) return(NA_real_)

  lo <- max(min(prior_df$value), min(post_df$value))
  hi <- min(max(prior_df$value), max(post_df$value))
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) return(0)
  grid <- seq(lo, hi, length.out = 1500)
  p <- approx(prior_df$value, prior_df$density, xout = grid, rule = 1, yleft = 0, yright = 0)$y
  q <- approx(post_df$value, post_df$density, xout = grid, rule = 1, yleft = 0, yright = 0)$y
  ov <- trapz_numeric(grid, pmin(p, q))
  if (is.finite(ov)) pmin(pmax(ov, 0), 1) else NA_real_
}

# Probability mass of a density between two bounds.
density_mass_between <- function(df, lower, upper) {
  df <- standardise_density(df)
  if (nrow(df) < 2 || !is.finite(lower) || !is.finite(upper) || lower >= upper) return(NA_real_)
  lo <- max(min(df$value), lower)
  hi <- min(max(df$value), upper)
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) return(0)
  grid <- seq(lo, hi, length.out = 1000)
  dens <- approx(df$value, df$density, xout = grid, rule = 1, yleft = 0, yright = 0)$y
  mass <- trapz_numeric(grid, dens)
  if (is.finite(mass)) pmin(pmax(mass, 0), 1) else NA_real_
}

# Compute one parameter prior-influence summary row.
prior_influence_row <- function(parameter, prior_df, posterior_marginal,
                                metric_scale = "identity",
                                class = "model_parameter") {
  if (is.null(posterior_marginal) || nrow(posterior_marginal) < 2) {
    return(data.frame(
      parameter = parameter,
      class = class,
      metric_scale = metric_scale,
      prior_mean = NA_real_, prior_sd = NA_real_, prior_q025 = NA_real_,
      prior_median = NA_real_, prior_q975 = NA_real_,
      posterior_mean = NA_real_, posterior_sd = NA_real_, posterior_q025 = NA_real_,
      posterior_median = NA_real_, posterior_q975 = NA_real_,
      posterior_prior_sd_ratio = NA_real_,
      posterior_shift_prior_sd = NA_real_,
      prior_mass_in_posterior_95 = NA_real_,
      density_overlap = NA_real_,
      sensitivity_priority_score = NA_real_,
      sensitivity_priority = NA,
      reason = "posterior marginal not available",
      row.names = NULL
    ))
  }

  post_df <- data.frame(value = posterior_marginal[, 1], density = posterior_marginal[, 2])
  p0 <- transform_density_for_metric(prior_df, metric_scale)
  p1 <- transform_density_for_metric(post_df, metric_scale)
  s0 <- density_summary_stats(p0)
  s1 <- density_summary_stats(p1)

  sd_ratio <- if (is.finite(s0$sd) && s0$sd > 0 && is.finite(s1$sd)) s1$sd / s0$sd else NA_real_
  shift <- if (is.finite(s0$sd) && s0$sd > 0 && is.finite(s0$median) && is.finite(s1$median)) {
    abs(s1$median - s0$median) / s0$sd
  } else NA_real_
  mass <- density_mass_between(p0, s1$q025, s1$q975)
  ov <- density_overlap(p0, p1)

  # Heuristic priority score for sensitivity. A high value can mean either:
  #   1) posterior remains wide relative to prior -> possible prior influence;
  #   2) posterior is in low-prior-density space -> prior-data tension.
  comp_width <- if (is.finite(sd_ratio)) pmin(1, sd_ratio / 0.60) else 0
  comp_mass <- if (is.finite(mass)) pmin(1, pmax(0, (0.10 - mass) / 0.10)) else 0
  comp_overlap <- if (is.finite(ov)) pmin(1, pmax(0, (0.20 - ov) / 0.20)) else 0
  score <- max(comp_width, comp_mass, comp_overlap, na.rm = TRUE)

  reasons <- character()
  if (is.finite(sd_ratio) && sd_ratio > 0.60) {
    reasons <- c(reasons, "posterior SD is >60% of prior SD; data may not dominate this prior")
  }
  if (is.finite(mass) && mass < 0.10) {
    reasons <- c(reasons, "posterior 95% interval had <10% prior mass; prior-data tension")
  }
  if (is.finite(ov) && ov < 0.20) {
    reasons <- c(reasons, "prior-posterior density overlap <20%; strong update/tension")
  }
  if (!length(reasons)) reasons <- "no obvious prior-sensitivity signal"

  data.frame(
    parameter = parameter,
    class = class,
    metric_scale = metric_scale,
    prior_mean = s0$mean,
    prior_sd = s0$sd,
    prior_q025 = s0$q025,
    prior_median = s0$median,
    prior_q975 = s0$q975,
    posterior_mean = s1$mean,
    posterior_sd = s1$sd,
    posterior_q025 = s1$q025,
    posterior_median = s1$median,
    posterior_q975 = s1$q975,
    posterior_prior_sd_ratio = sd_ratio,
    posterior_shift_prior_sd = shift,
    prior_mass_in_posterior_95 = mass,
    density_overlap = ov,
    sensitivity_priority_score = score,
    sensitivity_priority = is.finite(score) && score >= 0.5,
    reason = paste(reasons, collapse = "; "),
    row.names = NULL
  )
}

# Compute and write the prior-influence screen across parameters.
write_prior_influence_diagnostics <- function(fit, settings, family,
                                              label = "pre_sensitivity") {
  rows <- list()

  # Intercept on the log-rate scale.
  intercept_sd_prior <- 1 / sqrt(PRIOR_INTERCEPT_PREC)
  intercept_grid <- seq(qnorm(0.001, PRIOR_INTERCEPT_MEAN, intercept_sd_prior),
                        qnorm(0.999, PRIOR_INTERCEPT_MEAN, intercept_sd_prior),
                        length.out = 1200)
  rows[[length(rows) + 1L]] <- prior_influence_row(
    parameter = "intercept",
    prior_df = data.frame(value = intercept_grid,
                          density = dnorm(intercept_grid, PRIOR_INTERCEPT_MEAN, intercept_sd_prior)),
    posterior_marginal = fit$marginals.fixed[["intercept"]],
    metric_scale = "identity",
    class = "fixed_effect"
  )

  # Month fixed effects on the log-rate-ratio scale.
  month_terms <- grep("^month_[0-9]{4}_[0-9]{2}$", names(fit$marginals.fixed), value = TRUE)
  if (length(month_terms)) {
    month_sd_prior <- 1 / sqrt(PRIOR_MONTH_LOG_RATE_RATIO_PREC)
    month_grid <- seq(qnorm(0.001, 0, month_sd_prior),
                      qnorm(0.999, 0, month_sd_prior),
                      length.out = 1200)
    month_prior <- data.frame(value = month_grid,
                              density = dnorm(month_grid, 0, month_sd_prior))
    for (term in month_terms) {
      rows[[length(rows) + 1L]] <- prior_influence_row(
        parameter = paste0("month_", month_from_term(term), "_vs_", settings$month_reference),
        prior_df = month_prior,
        posterior_marginal = fit$marginals.fixed[[term]],
        metric_scale = "identity",
        class = "fixed_effect"
      )
    }
  }

  # Spatial range, evaluated on log(range) scale for interpretability.
  if (is.null(settings$fix_range_m)) {
    range0 <- settings$prior_range_m[1]
    prob_range <- settings$prior_range_m[2]
    range_grid <- exp(seq(log(pc_range_quantile(0.001, range0, prob_range)),
                          log(pc_range_quantile(0.995, range0, prob_range)),
                          length.out = 1200))
    rows[[length(rows) + 1L]] <- prior_influence_row(
      parameter = "spatial_range",
      prior_df = data.frame(value = range_grid,
                            density = pc_range_density(range_grid, range0, prob_range)),
      posterior_marginal = hyp_marg(fit, PAT_RANGE),
      metric_scale = "log",
      class = "spatial_hyperparameter"
    )
  }

  # Spatial marginal SD, native scale.
  sigma0 <- settings$prior_sigma[1]
  prob_sigma <- settings$prior_sigma[2]
  sigma_grid <- seq(0, pc_sigma_quantile(0.995, sigma0, prob_sigma), length.out = 1200)
  rows[[length(rows) + 1L]] <- prior_influence_row(
    parameter = "spatial_sd",
    prior_df = data.frame(value = sigma_grid,
                          density = pc_sigma_density(sigma_grid, sigma0, prob_sigma)),
    posterior_marginal = hyp_marg(fit, PAT_SIGMA),
    metric_scale = "identity",
    class = "spatial_hyperparameter"
  )

  # Zero inflation, only if relevant to the family.
  if (is_zi(family)) {
    prob_grid <- seq(0.001, 0.999, length.out = 1200)
    rows[[length(rows) + 1L]] <- prior_influence_row(
      parameter = "zero_inflation_probability",
      prior_df = data.frame(value = prob_grid,
                            density = zip_prob_prior_density(prob_grid)),
      posterior_marginal = hyp_marg(fit, PAT_ZPROB),
      metric_scale = "identity",
      class = "likelihood_hyperparameter"
    )
  }

  # NB size, evaluated on log(size) scale.
  if (is_nb(family)) {
    size_grid <- exp(seq(log(nb_size_prior_quantile(0.001)),
                         log(nb_size_prior_quantile(0.995)),
                         length.out = 1200))
    rows[[length(rows) + 1L]] <- prior_influence_row(
      parameter = "negative_binomial_size",
      prior_df = data.frame(value = size_grid,
                            density = nb_size_prior_density(size_grid)),
      posterior_marginal = hyp_marg(fit, PAT_NB_SIZE),
      metric_scale = "log",
      class = "likelihood_hyperparameter"
    )
  }

  out <- do.call(rbind, rows) %>%
    arrange(desc(sensitivity_priority_score), parameter)

  readr::write_csv(out, path_out(paste0(SURVEY_PREFIX, "_prior_influence_screen_", label, ".csv")))

  plot_df <- out %>%
    mutate(parameter = factor(parameter, levels = rev(parameter)))

  p <- ggplot(plot_df, aes(parameter, sensitivity_priority_score)) +
    geom_col() +
    coord_flip() +
    geom_hline(yintercept = 0.5, linetype = 2, colour = "grey35") +
    scale_y_continuous(limits = c(0, 1), labels = label_number(accuracy = 0.01)) +
    labs(
      title = paste0("Prior-sensitivity priority screen: ", SURVEY_PREFIX),
      subtitle = "Higher values flag parameters where sensitivity runs are most important",
      x = NULL,
      y = "sensitivity priority score"
    ) +
    theme_minimal(base_size = 12)

  ggsave(path_out(paste0(SURVEY_PREFIX, "_prior_influence_priority_", label, ".png")),
         p, width = 8, height = 5.5, dpi = 250)

  detail_plot <- ggplot(out,
                        aes(posterior_prior_sd_ratio,
                            prior_mass_in_posterior_95,
                            label = parameter)) +
    geom_vline(xintercept = 0.60, linetype = 2, colour = "grey55") +
    geom_hline(yintercept = 0.10, linetype = 2, colour = "grey55") +
    geom_point(size = 2.5) +
    geom_text(size = 3, check_overlap = TRUE, nudge_y = 0.03) +
    scale_x_continuous(labels = label_number(accuracy = 0.01)) +
    scale_y_continuous(limits = c(0, 1), labels = label_number(accuracy = 0.01)) +
    labs(
      title = paste0("Prior influence components: ", SURVEY_PREFIX),
      subtitle = "Upper-right means posterior remains prior-width; lower-left means posterior lies in low-prior region",
      x = "posterior SD / prior SD",
      y = "prior mass inside posterior 95% interval"
    ) +
    theme_minimal(base_size = 12)

  ggsave(path_out(paste0(SURVEY_PREFIX, "_prior_influence_components_", label, ".png")),
         detail_plot, width = 8, height = 5.5, dpi = 250)

  flagged <- out %>% filter(isTRUE(sensitivity_priority) | sensitivity_priority == TRUE)
  report <- c(
    "Prior-influence screen:",
    "  Purpose: decide which priors deserve explicit sensitivity reruns before final interpretation.",
    "  This is a heuristic screen, not a formal proof that a prior caused the posterior.",
    "",
    "How to read the metrics:",
    "  posterior_prior_sd_ratio close to 1 means the posterior remains almost as wide as the prior, so data may be weak for that parameter.",
    "  prior_mass_in_posterior_95 below 0.10 means the posterior sits in a low-prior-probability region, suggesting prior-data tension.",
    "  density_overlap below 0.20 means the prior and posterior distributions are very different.",
    "",
    if (nrow(flagged)) {
      c("Parameters prioritized for sensitivity:",
        paste0("  - ", flagged$parameter, ": ", flagged$reason))
    } else {
      "No parameter crossed the automatic sensitivity-priority thresholds."
    },
    "",
    "Full numeric table written to:",
    paste0("  ", path_out(paste0(SURVEY_PREFIX, "_prior_influence_screen_", label, ".csv")))
  )
  writeLines(report, path_out(paste0(SURVEY_PREFIX, "_PRIOR_INFLUENCE_SCREEN_", label, ".txt")))

  invisible(out)
}

# Write month fixed-effect rate ratios with credible intervals.
write_month_coefficients <- function(fit, model_dat, settings) {
  month_terms <- temporal_month_terms(model_dat)
  if (!length(month_terms)) return(invisible(NULL))
  month_summary <- fit$summary.fixed[month_terms, , drop = FALSE]
  out <- data.frame(
    term = rownames(month_summary),
    month = month_from_term(rownames(month_summary)),
    reference_month = settings$month_reference,
    prediction_month = settings$month_prediction,
    mean_log_rate_ratio = month_summary[, "mean"],
    q025_log_rate_ratio = month_summary[, "0.025quant"],
    q975_log_rate_ratio = month_summary[, "0.975quant"],
    mean_rate_ratio = exp(month_summary[, "mean"]),
    q025_rate_ratio = exp(month_summary[, "0.025quant"]),
    q975_rate_ratio = exp(month_summary[, "0.975quant"]),
    row.names = NULL
  )
  readr::write_csv(out, path_out(paste0(SURVEY_PREFIX, "_month_coefficients.csv")))
  invisible(out)
}

# Posterior mean log-rate-ratio of a month versus the reference month.
month_log_rate_ratio_mean <- function(fit, month, settings) {
  if (is.null(month) || is.na(month) || !nzchar(month)) return(0)
  if (!is.null(settings$month_reference) &&
      !is.na(settings$month_reference) &&
      month == settings$month_reference) {
    return(0)
  }

  term <- month_term_name(month)
  fixed_names <- rownames(fit$summary.fixed)
  if (term %in% fixed_names) {
    return(as.numeric(fit$summary.fixed[term, "mean"]))
  }

  0
}

# Save summaries of the draw-wise aggregation factor. The map computation uses
# each joint draw's own factor, never this summary as a constant multiplier.
write_annualization_weights <- function(annualization) {
  weights <- annualization$weights
  weights$factor_posterior_mean <- annualization$factor
  weights$factor_posterior_sd <- annualization$factor_sd
  weights$factor_q05 <- annualization$factor_q05
  weights$factor_q95 <- annualization$factor_q95
  readr::write_csv(weights, path_out(paste0(SURVEY_PREFIX, "_annualization_weights.csv")))
  invisible(annualization)
}


## 12. Fit model and prediction maps -----------------------------------------

# Fit the final INLA-SPDE mapping model used to produce the map surfaces.
fit_final_model <- function(model_dat, settings, family) {
  cat(sprintf("\n[wolf] fitting %s: family=%s\n", FINAL_MODEL_NAME, family))

  model_dat <- model_dat %>% mutate(y = as.integer(wolf_events), intercept = 1)

  camera_summary <- camera_summary_from_model(model_dat)
  camera_sf <- camera_to_utm(camera_summary)
  coords_camera <- st_coordinates(camera_sf)
  colnames(coords_camera) <- c("x", "y")

  obs_sf <- model_dat %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(EPSG_UTM)
  coords_obs <- st_coordinates(obs_sf)
  colnames(coords_obs) <- c("x", "y")

  pred_sf <- prediction_grid(camera_sf, settings)
  coords_pred <- st_coordinates(pred_sf)
  colnames(coords_pred) <- c("x", "y")
  cat(sprintf("[wolf] prediction cells: %d at %.0f m\n",
              nrow(pred_sf), settings$cell_size_m))

  # The SPDE mesh is built on the camera locations; prediction cells are mapped
  # onto that mesh by the projector matrix (A_pred) rather than becoming nodes.
  spde_obj <- build_spatial(coords_camera, settings)

  fixed_terms <- fixed_effect_terms(model_dat)
  fixed_obs <- as.data.frame(model_dat[, fixed_terms, drop = FALSE])
  fixed_pred <- prediction_fixed_effects(model_dat, fixed_terms, settings, nrow(pred_sf))

  A_obs <- INLA::inla.spde.make.A(spde_obj$mesh, loc = coords_obs)
  A_pred <- INLA::inla.spde.make.A(spde_obj$mesh, loc = coords_pred)

  stack_obs <- INLA::inla.stack(
    tag = "obs",
    data = list(y = model_dat$y, e = model_dat$total_effort_days),
    A = list(A_obs, 1),
    effects = list(spatial = spde_obj$s_index, fixed = fixed_obs)
  )

  stack_pred <- INLA::inla.stack(
    tag = "pred",
    data = list(y = rep(NA_real_, nrow(pred_sf)), e = rep(100, nrow(pred_sf))),
    A = list(A_pred, 1),
    effects = list(spatial = spde_obj$s_index, fixed = fixed_pred)
  )

  stack_all <- INLA::inla.stack(stack_obs, stack_pred)
  stack_data <- INLA::inla.stack.data(stack_all)

  obs_index <- INLA::inla.stack.index(stack_all, tag = "obs")$data
  pred_index <- INLA::inla.stack.index(stack_all, tag = "pred")$data

  formula <- as.formula(
    paste("y ~ 0 +",
          paste(c(fixed_terms, "f(spatial, model = spde_obj$spde)"), collapse = " + "))
  )

  fit <- INLA::inla(
    formula,
    family = family,
    data = stack_data,
    E = stack_data$e,
    control.predictor = list(
      A = INLA::inla.stack.A(stack_all),
      compute = TRUE,
      link = 1
    ),
    control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE),
    control.fixed = make_control_fixed(fixed_terms),
    control.family = make_control_family(family),
    verbose = FALSE
  )

  list(fit = fit, stack_all = stack_all, obs_index = obs_index,
       pred_index = pred_index, model_dat = model_dat,
       camera_sf = camera_sf, pred_sf = pred_sf,
       spde_obj = spde_obj, fixed_terms = fixed_terms)
}


# Fit the observed-data diagnostic model used for joint posterior sampling.
fit_diagnostic_model <- function(model_dat, settings, family) {
  cat(sprintf("\n[wolf] fitting observed-data diagnostic model for posterior sampling\n"))

  model_dat <- model_dat %>% mutate(y = as.integer(wolf_events), intercept = 1)

  camera_summary <- camera_summary_from_model(model_dat)
  camera_sf <- camera_to_utm(camera_summary)
  coords_camera <- st_coordinates(camera_sf)
  colnames(coords_camera) <- c("x", "y")

  obs_sf <- model_dat %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(EPSG_UTM)
  coords_obs <- st_coordinates(obs_sf)
  colnames(coords_obs) <- c("x", "y")

  spde_obj <- build_spatial(coords_camera, settings)
  fixed_terms <- fixed_effect_terms(model_dat)
  fixed_obs <- as.data.frame(model_dat[, fixed_terms, drop = FALSE])
  A_obs <- INLA::inla.spde.make.A(spde_obj$mesh, loc = coords_obs)

  stack_obs <- INLA::inla.stack(
    tag = "obs",
    data = list(y = model_dat$y, e = model_dat$total_effort_days),
    A = list(A_obs, 1),
    effects = list(spatial = spde_obj$s_index, fixed = fixed_obs)
  )

  stack_data <- INLA::inla.stack.data(stack_obs)
  obs_index <- INLA::inla.stack.index(stack_obs, tag = "obs")$data

  formula <- as.formula(
    paste("y ~ 0 +",
          paste(c(fixed_terms, "f(spatial, model = spde_obj$spde)"), collapse = " + "))
  )

  fit <- INLA::inla(
    formula,
    family = family,
    data = stack_data,
    E = stack_data$e,
    control.predictor = list(
      A = INLA::inla.stack.A(stack_obs),
      compute = TRUE,
      link = 1
    ),
    control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE),
    control.fixed = make_control_fixed(fixed_terms),
    control.family = make_control_family(family),
    verbose = FALSE
  )

  list(fit = fit, stack_obs = stack_obs, obs_index = obs_index,
       model_dat = model_dat, camera_sf = camera_sf,
       spde_obj = spde_obj, fixed_terms = fixed_terms)
}

# Build the annualized posterior mean/SD/CV prediction surfaces and rasters.
make_prediction_outputs <- function(fit_obj, diag, settings, family) {
  fit <- fit_obj$fit
  pred_sf <- fit_obj$pred_sf
  pred_index <- fit_obj$pred_index
  joint <- joint_surface_summary(fit_obj, diag$model_dat, settings, family)
  annualization <- joint$annualization
  write_annualization_weights(annualization)
  pred_sf <- cbind(pred_sf, joint$summary)
  pred_sf$posterior_nsim <- joint$nsim
  pred_sf$summary_method <- joint$method


  model_dat <- diag$model_dat
  overall_rate <- 100 * sum(model_dat$y) / sum(model_dat$total_effort_days)

  coords_pred <- st_coordinates(pred_sf)
  pred_sf$x <- coords_pred[, 1]
  pred_sf$y <- coords_pred[, 2]

  wkt <- st_crs(pred_sf)$wkt
  pred_table <- st_drop_geometry(pred_sf)

  readr::write_csv(pred_table,
                   path_out(paste0(SURVEY_PREFIX, "_final_prediction_grid.csv")))

  make_raster <- function(col) {
    terra::rast(pred_table[, c("x", "y", col)], type = "xyz", crs = wkt)
  }

  r_mean <- make_raster("mean")
  r_sd <- make_raster("sd")
  r_cv <- make_raster("cv")
  r_q05 <- make_raster("q05")
  r_q95 <- make_raster("q95")
  names(r_mean) <- paste0(TARGET_SLUG, "_events_per_100_camera_days")
  names(r_sd) <- "posterior_sd"
  names(r_cv) <- "posterior_cv"
  names(r_q05) <- "posterior_q05"
  names(r_q95) <- "posterior_q95"

  terra::writeRaster(r_mean,
                     path_out(paste0(SURVEY_PREFIX,
                                     "_final_predicted_events_per_100_days_mean.tif")),
                     overwrite = TRUE)
  terra::writeRaster(r_sd,
                     path_out(paste0(SURVEY_PREFIX,
                                     "_final_predicted_events_per_100_days_sd.tif")),
                     overwrite = TRUE)
  terra::writeRaster(r_cv,
                     path_out(paste0(SURVEY_PREFIX,
                                     "_final_predicted_events_per_100_days_cv.tif")),
                     overwrite = TRUE)

  terra::writeRaster(r_q05, path_out(paste0(SURVEY_PREFIX,
    "_final_predicted_events_per_100_days_q05.tif")), overwrite = TRUE)
  terra::writeRaster(r_q95, path_out(paste0(SURVEY_PREFIX,
    "_final_predicted_events_per_100_days_q95.tif")), overwrite = TRUE)
  rasters <- list(mean = r_mean, sd = r_sd, cv = r_cv, q05 = r_q05, q95 = r_q95)

  plot_map_outputs(fit_obj$camera_sf, diag$model_dat, rasters, overall_rate,
                   annualization)
  invisible(list(pred_sf = pred_sf, rasters = rasters,
                 overall_rate = overall_rate,
                 annualization = annualization, method = joint$method,
                 nsim = joint$nsim, seed = joint$seed))
}

# Render and save the mean, SD, and CV map figures.
plot_map_outputs <- function(camera_sf, model_dat, rasters, overall_rate,
                             annualization = NULL) {
  plot_label <- sub(" survey$", "", SURVEY_LABEL)
  agg_word <- if (!is.null(annualization) &&
                  grepl("single-period", annualization$label, fixed = TRUE)) {
    ""
  } else {
    "annualized "
  }

  raster_to_df <- function(r, name) {
    d <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
    names(d) <- c("x", "y", name)
    d
  }

  mean_df <- raster_to_df(rasters$mean, "rate")
  cap <- quantile(mean_df$rate, 0.98, na.rm = TRUE)

  camera_obs <- camera_summary_from_model(model_dat)
  positive_sf <- camera_sf
  positive_sf$wolf_events_per_100_days <- camera_obs$wolf_events_per_100_days
  positive_sf$wolf_events <- camera_obs$wolf_events
  positive_sf <- positive_sf[camera_obs$wolf_events > 0, ]

  mean_plot <- ggplot() +
    geom_raster(data = mean_df,
                aes(x, y, fill = pmin(rate, cap)),
                interpolate = TRUE) +
    geom_sf(data = camera_sf,
            shape = 21, size = 1.4, fill = "white",
            colour = "grey35", stroke = 0.25) +
    geom_sf(data = positive_sf,
            aes(size = wolf_events_per_100_days),
            shape = 21, fill = "black", colour = "white",
            stroke = 0.25, alpha = 0.9) +
    # Guide order is set explicitly: with two legends and no `order`, ggplot's
    # box packing decides the stacking, which is not stable across versions and
    # made this figure the only non-reproducible output in the pipeline.
    scale_fill_viridis_c(option = "magma", na.value = NA,
                         name = "predicted events\n/100 camera-days",
                         labels = label_number(accuracy = 0.01),
                         guide = guide_colourbar(order = 1)) +
    scale_size_continuous(range = c(2, 7),
                          name = "observed events\n/100 camera-days",
                          labels = label_number(accuracy = 0.01),
                          guide = guide_legend(order = 2)) +
    coord_sf(datum = NA) +
    labs(title = paste0(TARGET_TITLE, " encounter-frequency surface: ", plot_label),
         subtitle = paste0("posterior mean of ", agg_word, "expected encounter frequency"),
         x = "Easting, UTM 34N", y = "Northing, UTM 34N") +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(), legend.position = "right")

  ggsave(path_out(paste0(SURVEY_PREFIX, "_final_event_frequency_mean.png")),
         mean_plot, width = 9.5, height = 9, dpi = 150)

  sd_df <- raster_to_df(rasters$sd, "sd")
  sd_cap <- quantile(sd_df$sd, 0.98, na.rm = TRUE)
  sd_plot <- ggplot() +
    geom_raster(data = sd_df, aes(x, y, fill = pmin(sd, sd_cap)),
                interpolate = TRUE) +
    geom_sf(data = camera_sf,
            shape = 21, size = 1.4, fill = "white",
            colour = "grey35", stroke = 0.25) +
    scale_fill_viridis_c(option = "viridis", na.value = NA,
                         name = "posterior SD\n(events /100 camera-days)",
                         labels = label_number(accuracy = 0.01)) +
    coord_sf(datum = NA) +
    labs(title = paste0("Uncertainty surface: ", plot_label),
         subtitle = paste0("posterior standard deviation of ", agg_word, "expected encounter frequency"),
         x = "Easting, UTM 34N", y = "Northing, UTM 34N") +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(), legend.position = "right")

  ggsave(path_out(paste0(SURVEY_PREFIX, "_final_event_frequency_sd.png")),
         sd_plot, width = 9.5, height = 9, dpi = 150)

  cv_df <- raster_to_df(rasters$cv, "cv")
  cv_cap <- quantile(cv_df$cv, 0.98, na.rm = TRUE)
  cv_plot <- ggplot() +
    geom_raster(data = cv_df, aes(x, y, fill = pmin(cv, cv_cap)),
                interpolate = TRUE) +
    geom_sf(data = camera_sf,
            shape = 21, size = 1.4, fill = "white",
            colour = "grey35", stroke = 0.25) +
    scale_fill_viridis_c(option = "cividis", na.value = NA,
                         name = "posterior CV\n(SD / mean, unitless)",
                         labels = label_number(accuracy = 0.01)) +
    coord_sf(datum = NA) +
    labs(title = paste0("Relative uncertainty surface: ", plot_label),
         subtitle = paste0("posterior coefficient of variation (SD / mean) of ", agg_word, "expected encounter frequency"),
         x = "Easting, UTM 34N", y = "Northing, UTM 34N") +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(), legend.position = "right")

  ggsave(path_out(paste0(SURVEY_PREFIX, "_final_event_frequency_cv.png")),
         cv_plot, width = 9.5, height = 9, dpi = 150)

  invisible(TRUE)
}


## 13. Spatial block cross-validation ----------------------------------------

# Run spatial block cross-validation with a per-fold mesh rebuild.
spatial_block_cv <- function(model_dat, settings, family, K = CV_K,
                              folds = NULL, dependence = "none", output_label = "final") {
  cat(sprintf("\n[wolf] spatial block CV for %s (K=%d)\n", FINAL_MODEL_NAME, K))

  model_dat <- model_dat %>% mutate(y = as.integer(wolf_events), intercept = 1)
  fixed_terms <- fixed_effect_terms(model_dat)

  camera_summary <- camera_summary_from_model(model_dat)
  camera_sf <- camera_to_utm(camera_summary)
  coords_camera <- st_coordinates(camera_sf)
  colnames(coords_camera) <- c("x", "y")

  row_sf <- model_dat %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(EPSG_UTM)
  coords_row <- st_coordinates(row_sf)
  colnames(coords_row) <- c("x", "y")

  y <- model_dat$y
  effort <- model_dat$total_effort_days

  if (is.null(folds)) {
    K_final <- min(K, nrow(coords_camera) - 1L)
    fold <- strict_spatial_folds(coords_camera, K_final, seed = SPATIAL_FOLD_SEED)
    folds <- data.frame(plotID = camera_summary$plotID,
                        x = coords_camera[, 1L], y = coords_camera[, 2L], fold = fold)
    readr::write_csv(folds, path_out(paste0(SURVEY_PREFIX, "_spatial_cv_folds.csv")))
  } else {
    if (anyDuplicated(folds$plotID) || !setequal(folds$plotID, camera_summary$plotID))
      stop("Saved spatial folds do not match the complete camera set.")
    matched <- match(camera_summary$plotID, folds$plotID)
    if (any(abs(folds$x[matched] - coords_camera[, 1L]) > 1e-6) ||
        any(abs(folds$y[matched] - coords_camera[, 2L]) > 1e-6))
      stop("Camera coordinates do not match the saved folds.")
    fold <- folds$fold[matched]
    if (anyNA(fold) || length(unique(fold)) < 2L)
      stop("Saved folds must assign every camera to at least two nonempty folds.")
  }
  row_fold <- fold[match(model_dat$plotID, camera_summary$plotID)]
  if (anyNA(row_fold)) stop("Some model rows lack a spatial fold.")
  dep <- dependence_design(model_dat, dependence)
  prior_rows <- list()

  rows <- list()
  cam_rows <- list()
  failed_folds <- character()

  for (f in sort(unique(fold))) {
    cat(sprintf("[wolf]   CV fold %s\n", f))
    test <- which(row_fold == f)
    train <- which(row_fold != f)

    result <- tryCatch({
      # Strict CV: train mesh only. Offset is large enough for edge-fold prediction.
      train_camera_ids <- unique(model_dat$plotID[train])
      train_camera_coords <- coords_camera[camera_summary$plotID %in% train_camera_ids, , drop = FALSE]
      spde_obj <- build_spatial(train_camera_coords, settings)

      A_train <- INLA::inla.spde.make.A(spde_obj$mesh,
                                        loc = coords_row[train, , drop = FALSE])
      A_test <- INLA::inla.spde.make.A(spde_obj$mesh,
                                       loc = coords_row[test, , drop = FALSE])

      fixed_train <- cbind(as.data.frame(model_dat[train, fixed_terms, drop = FALSE]),
                           dep$indices[train, , drop = FALSE])
      fixed_test <- cbind(as.data.frame(model_dat[test, fixed_terms, drop = FALSE]),
                          dep$indices[test, , drop = FALSE])
      fold_prior_mean <- training_intercept_mean(y[train], effort[train])

      stack_train <- INLA::inla.stack(
        tag = "train",
        data = list(y = y[train], e = effort[train]),
        A = list(A_train, 1),
        effects = list(spatial = spde_obj$s_index, fixed = fixed_train)
      )

      stack_test <- INLA::inla.stack(
        tag = "test",
        data = list(y = rep(NA_real_, length(test)), e = effort[test]),
        A = list(A_test, 1),
        effects = list(spatial = spde_obj$s_index, fixed = fixed_test)
      )

      stack_all <- INLA::inla.stack(stack_train, stack_test)
      stack_data <- INLA::inla.stack.data(stack_all)

      formula <- as.formula(
        paste("y ~ 0 +",
              paste(c(fixed_terms, "f(spatial, model = spde_obj$spde)", dep$term),
                    collapse = " + "))
      )

      fit_fold <- INLA::inla(
        formula,
        family = family,
        data = stack_data,
        E = stack_data$e,
        control.predictor = list(
          A = INLA::inla.stack.A(stack_all),
          compute = TRUE,
          link = 1
        ),
        control.compute = list(config = TRUE),
        control.fixed = make_control_fixed(fixed_terms, intercept_mean = fold_prior_mean),
        control.family = make_control_family(family),
        verbose = FALSE
      )

      test_index <- INLA::inla.stack.index(stack_all, tag = "test")$data
      samples <- posterior_samples_safe(fit_fold, CV_NSIM,
                                        seed = INLA_SAMPLE_SEED + 100L + as.integer(f))
      test_draws <- build_posterior_draws(fit_fold, samples, test_index,
                                          effort[test], family,
                                          expected_n_stack = length(fit_fold$summary.linear.predictor$mean))
      sim <- simulate_from_draws(test_draws, family)
      Ey <- rowMeans(test_draws$fitted)
      lo <- apply(sim, 1, quantile, 0.05, na.rm = TRUE)
      hi <- apply(sim, 1, quantile, 0.95, na.rm = TRUE)

      lpd <- vapply(seq_along(test), function(j) {
        log_mean_exp(fam_logpmf(y[test[j]], test_draws$mu[j, ],
                                test_draws$pi, family, test_draws$size))
      }, numeric(1))

      row_out <- data.frame(
        fold = f,
        plotID = model_dat$plotID[test],
        deploymentID = model_dat$deploymentID[test],
        month = model_dat$month[test],
        y = y[test],
        Ey = Ey,
        effort_days = effort[test],
        rate_obs = 100 * y[test] / effort[test],
        rate_pred = 100 * Ey / effort[test],
        lpd = lpd,
        lo90 = lo,
        hi90 = hi,
        covered_90 = y[test] >= lo & y[test] <= hi
      )

      group <- as.factor(model_dat$plotID[test])
      sim_cam <- aggregate_matrix_by_group(sim, group)
      if (is.null(dim(sim_cam))) sim_cam <- matrix(sim_cam, nrow = 1)
      y_cam <- as.numeric(rowsum(y[test], group)[, 1])
      effort_cam <- as.numeric(rowsum(effort[test], group)[, 1])
      Ey_cam <- rowMeans(sim_cam)
      lo_cam <- apply(sim_cam, 1, quantile, 0.05, na.rm = TRUE)
      hi_cam <- apply(sim_cam, 1, quantile, 0.95, na.rm = TRUE)
      cam_ids <- rownames(rowsum(y[test], group))

      cam_out <- data.frame(
        fold = f,
        plotID = cam_ids,
        y = y_cam,
        Ey = Ey_cam,
        effort_days = effort_cam,
        rate_obs = 100 * y_cam / effort_cam,
        rate_pred = 100 * Ey_cam / effort_cam,
        lo90 = lo_cam,
        hi90 = hi_cam,
        covered_90 = y_cam >= lo_cam & y_cam <= hi_cam
      )

      list(row = row_out, camera = cam_out,
           prior = data.frame(fold = f, dependence = dependence,
             training_rows = length(train), held_out_rows = length(test),
             training_events = sum(y[train]), training_effort = sum(effort[train]),
             intercept_prior_mean = fold_prior_mean))
    }, error = function(e) {
      failed_folds <<- c(failed_folds, paste0(f, ": ", conditionMessage(e)))
      NULL
    })

    if (!is.null(result)) {
      rows[[length(rows) + 1L]] <- result$row
      cam_rows[[length(cam_rows) + 1L]] <- result$camera
      prior_rows[[length(prior_rows) + 1L]] <- result$prior
    }
  }

  if (length(failed_folds)) {
    stop("[wolf] spatial CV failed for fold(s): ",
         paste(failed_folds, collapse = "; "))
  }
  if (!length(rows)) stop("[wolf] spatial CV produced no successful folds.")

  cv_row <- do.call(rbind, rows)
  cv_cam <- do.call(rbind, cam_rows)

  readr::write_csv(cv_row, path_out(paste0(SURVEY_PREFIX, "_", output_label, "_spatial_block_cv_rows.csv")))
  readr::write_csv(cv_cam, path_out(paste0(SURVEY_PREFIX, "_", output_label, "_spatial_block_cv_camera.csv")))

  summary <- bind_rows(
    data.frame(
      level = "model_row",
      metric = c("mean_log_predictive_density", "rmse_count",
                 "rmse_rate_per100", "coverage_90"),
      value = c(mean(cv_row$lpd),
                sqrt(mean((cv_row$y - cv_row$Ey)^2)),
                sqrt(mean((cv_row$rate_obs - cv_row$rate_pred)^2)),
                mean(cv_row$covered_90))
    ),
    data.frame(
      level = "camera",
      metric = c("rmse_count", "rmse_rate_per100", "coverage_90"),
      value = c(sqrt(mean((cv_cam$y - cv_cam$Ey)^2)),
                sqrt(mean((cv_cam$rate_obs - cv_cam$rate_pred)^2)),
                mean(cv_cam$covered_90))
    )
  )

  readr::write_csv(summary,
                   path_out(paste0(SURVEY_PREFIX, "_", output_label, "_spatial_block_cv_summary.csv")))

  cat(sprintf(
    "[wolf] spatial CV rows: mean LPD %.3f | RMSE count %.2f | RMSE rate %.2f | 90%% coverage %.2f\n",
    summary$value[summary$level == "model_row" & summary$metric == "mean_log_predictive_density"],
    summary$value[summary$level == "model_row" & summary$metric == "rmse_count"],
    summary$value[summary$level == "model_row" & summary$metric == "rmse_rate_per100"],
    summary$value[summary$level == "model_row" & summary$metric == "coverage_90"]
  ))
  cat(sprintf(
    "[wolf] spatial CV cameras: RMSE count %.2f | RMSE rate %.2f | 90%% coverage %.2f\n",
    summary$value[summary$level == "camera" & summary$metric == "rmse_count"],
    summary$value[summary$level == "camera" & summary$metric == "rmse_rate_per100"],
    summary$value[summary$level == "camera" & summary$metric == "coverage_90"]
  ))

  calibration <- cv_calibration_summary(cv_row, cv_cam)
  fold_metrics <- do.call(rbind, lapply(sort(unique(cv_row$fold)), function(f) {
    z <- cv_calibration_summary(cv_row[cv_row$fold == f, , drop = FALSE],
                                cv_cam[cv_cam$fold == f, , drop = FALSE])
    z$fold <- f
    z
  }))
  readr::write_csv(calibration, path_out(paste0(SURVEY_PREFIX, "_", output_label,
                                               "_spatial_cv_calibration.csv")))
  readr::write_csv(fold_metrics, path_out(paste0(SURVEY_PREFIX, "_", output_label,
                                                "_spatial_cv_fold_metrics.csv")))
  readr::write_csv(do.call(rbind, prior_rows), path_out(paste0(SURVEY_PREFIX, "_", output_label,
                                                             "_spatial_cv_training_priors.csv")))
  list(row = cv_row, camera = cv_cam, summ = summary, folds = folds,
       calibration = calibration, fold_metrics = fold_metrics)
}


## 14. Reporting --------------------------------------------------------------

# List which required diagnostics failed (empty if all pass).
diagnostic_failures <- function(diag) {
  c(
    if (!isTRUE(diag$ppc_total_pass)) "camera-level PPC total events" else NULL,
    if (!isTRUE(diag$ppc_zero_pass)) "camera-level PPC zero fraction" else NULL,
    if (!isTRUE(diag$ppc_max_pass)) "camera-level PPC max count" else NULL,
    if (!isTRUE(diag$moran_pass)) {
      if (is.finite(diag$moran_p)) {
        sprintf("residual spatial autocorrelation (Moran's I = %.3f, p = %.3f)",
                diag$moran_I, diag$moran_p)
      } else {
        "residual Moran's I was not evaluable"
      }
    } else NULL
  )
}

# Format the prior specification as report text lines.
prior_lines_for_report <- function(settings, family) {
  c(
    "",
    "Priors:",
    sprintf("  Field range: PC prior, P(range < %d m) = %.2f",
            as.integer(settings$prior_range_m[1]), settings$prior_range_m[2]),
    sprintf("  Field marginal SD: PC prior, P(SD > %.2f) = %.2f",
            settings$prior_sigma[1], settings$prior_sigma[2]),
    sprintf("  Intercept: Gaussian(mean = %.3f, prec = %.3f), SD %.1f on log scale",
            PRIOR_INTERCEPT_MEAN, PRIOR_INTERCEPT_PREC, 1 / sqrt(PRIOR_INTERCEPT_PREC)),
    sprintf("  Month log-rate ratios: Gaussian(0, prec = %g), SD %.1f",
            PRIOR_MONTH_LOG_RATE_RATIO_PREC,
            1 / sqrt(PRIOR_MONTH_LOG_RATE_RATIO_PREC)),
    if (is_zi(family)) {
      sprintf("  Zero-inflation logit(p): Gaussian(%.2f, prec = %.2f)",
              PRIOR_ZI_LOGIT_MEAN, PRIOR_ZI_LOGIT_PREC)
    } else NULL,
    if (is_nb(family)) {
      sprintf("  Negative-binomial log(size): Gaussian(mean = %.2f, prec = %.3f), SD %.1f",
              PRIOR_NB_LOGSIZE_MEAN,
              PRIOR_NB_LOGSIZE_PREC,
              1 / sqrt(PRIOR_NB_LOGSIZE_PREC))
    } else NULL
  )
}

# Write the consolidated validation report.
write_validation_report <- function(model_dat, diag, cv, prediction,
                                    temporal_diag = NULL, dependence_comparison = NULL) {
  failures <- diagnostic_failures(diag)
  passes_required <- isTRUE(diag$diagnostics_ok)
  n_cameras <- dplyr::n_distinct(model_dat$plotID)
  annualization <- prediction$annualization
  review <- validation_review_status(diag, cv, temporal_diag, dependence_comparison,
                                     comparison_required = RUN_DEPENDENCE_COMPARISON)

  cv_lines <- if (!is.null(cv)) {
    c(
      "",
      "Spatial block cross-validation:",
      sprintf("  Row mean LPD: %.3f",
              cv$summ$value[cv$summ$level == "model_row" &
                              cv$summ$metric == "mean_log_predictive_density"]),
      sprintf("  Row RMSE count: %.2f",
              cv$summ$value[cv$summ$level == "model_row" &
                              cv$summ$metric == "rmse_count"]),
      sprintf("  Row RMSE rate /100: %.2f",
              cv$summ$value[cv$summ$level == "model_row" &
                              cv$summ$metric == "rmse_rate_per100"]),
      sprintf("  Row 90%% coverage: %.2f",
              cv$summ$value[cv$summ$level == "model_row" &
                              cv$summ$metric == "coverage_90"]),
      sprintf("  Camera RMSE count: %.2f",
              cv$summ$value[cv$summ$level == "camera" &
                              cv$summ$metric == "rmse_count"]),
      sprintf("  Camera RMSE rate /100: %.2f",
              cv$summ$value[cv$summ$level == "camera" &
                              cv$summ$metric == "rmse_rate_per100"]),
      sprintf("  Camera 90%% coverage: %.2f",
              cv$summ$value[cv$summ$level == "camera" &
                              cv$summ$metric == "coverage_90"])
    )
  } else {
    c("", "Spatial block cross-validation: not run")
  }



  temporal_lines <- if (!is.null(temporal_diag) &&
                        !is.null(temporal_diag$lag_summary)) {
    lag1 <- temporal_diag$lag_summary %>% filter(lag == 1L)
    acf1 <- if (!is.null(temporal_diag$acf_summary) &&
                nrow(temporal_diag$acf_summary) >= 2 &&
                any(temporal_diag$acf_summary$lag == 1L)) {
      temporal_diag$acf_summary$acf[temporal_diag$acf_summary$lag == 1L][[1]]
    } else {
      NA_real_
    }

    c(
      "",
      "Temporal residual autocorrelation diagnostics:",
      "  Month is included as a fixed effect; these diagnostics check residual temporal structure after that adjustment.",
      if (nrow(lag1)) {
        sprintf(
          "  Within-camera lag-1 residual correlation: r = %.3f, p = %.4g, n pairs = %d, median gap = %.1f days",
          lag1$correlation[[1]],
          lag1$p_value[[1]],
          lag1$n_pairs[[1]],
          lag1$median_days_between[[1]]
        )
      } else {
        "  Within-camera lag-1 residual correlation: not evaluable"
      },
      if (is.finite(acf1)) {
        sprintf("  Date-ordered mean residual ACF, lag 1: %.3f", acf1)
      } else {
        "  Date-ordered mean residual ACF: not evaluable"
      }
    )
  } else {
    c(
      "",
      "Temporal residual autocorrelation diagnostics: not run"
    )
  }

  report <- c(
    sprintf("Survey: %s", SURVEY_LABEL),
    sprintf("Final model: %s (family = %s)", FINAL_MODEL_NAME, FINAL_FAMILY),
    sprintf("Run profile: %s", RUN_PROFILE),
    sprintf("Joint posterior PPC simulations: %d", diag$ppc_nsim),
    sprintf(
      "Cameras: %d | model rows: %d | positive rows: %d | events: %d | effort: %.1f camera-days | observed mean %.3f /100",
      n_cameras,
      nrow(model_dat),
      sum(model_dat$wolf_events > 0),
      sum(model_dat$wolf_events),
      sum(model_dat$total_effort_days),
      100 * sum(model_dat$wolf_events) / sum(model_dat$total_effort_days)
    ),
    "",
    "Validation:",
    sprintf("  PPC method: %s", diag$ppc_method),
    sprintf("  Pearson dispersion, model rows: %.3f", diag$pearson_disp),
    sprintf("  Pearson dispersion, camera aggregates: %.3f", diag$pearson_disp_camera),
    sprintf("  Camera-level PPC total events pass: %s", isTRUE(diag$ppc_total_pass)),
    sprintf("  Camera-level PPC zero fraction pass: %s", isTRUE(diag$ppc_zero_pass)),
    sprintf("  Camera-level PPC max count pass: %s", isTRUE(diag$ppc_max_pass)),
    sprintf("  Residual Moran's I: %.3f; expected %.3f; two-sided p = %.3f",
            diag$moran_I, diag$moran_expected, diag$moran_p),
    sprintf("  Residual Moran pass: %s", isTRUE(diag$moran_pass)),
    sprintf("  Row PIT KS p: %.4g", diag$ppc_pit_ks_row),
    sprintf("  Camera PIT KS p: %.4g", diag$ppc_pit_ks_camera),
    if (is_zi(FINAL_FAMILY)) sprintf("  Zero-inflation probability posterior mean: %.3f", diag$pi_hat) else NULL,
    if (is_nb(FINAL_FAMILY)) sprintf("  Negative-binomial size posterior mean: %.3f", diag$size_hat) else NULL,
    sprintf("  Passes required checks: %s", passes_required),
    if (!passes_required) {
      sprintf("  Stated limitation: %s", paste(failures, collapse = "; "))
    } else NULL,
    sprintf("  Overall review status: %s", review$status),
    if (length(review$reasons)) paste0("  Review reason: ", review$reasons) else NULL,
    "  A legacy PPC/Moran pass does not establish predictive calibration or temporal robustness.",
    cv_lines,
    "  Coverage below 0.90 is a descriptive review flag, not an independence-based significance test.",
    "  See spatial_cv_calibration.csv and spatial_cv_fold_metrics.csv for counts, fold ranges, widths, and interval scores.",
    temporal_lines,
    "  Residual p-values are exploratory screens; nonsignificance does not rule out a mechanism.",
    "  Spatial CV and mesh stability alone cannot establish absence of temporal distortion.",
    sprintf("  Dependence comparison: %s", if (is.null(dependence_comparison)) "not run" else
      if (any(dependence_comparison$status != "completed")) "incomplete; review required" else "completed; inspect comparison report"),
    prior_lines_for_report(settings, FINAL_FAMILY),
    "",
    "Temporal structure:",
    if (isTRUE(settings$use_month_effect)) {
      c("  Calendar camera-month is included as a fixed effect after aligning effort and events by eventStart month.",
        sprintf("  Coefficient-coding baseline month: %s", MONTH_REFERENCE),
        sprintf("  Prediction-stack baseline month used internally: %s", MONTH_PREDICTION),
        sprintf("  Map aggregation: effort-weighted sampled-month %s surface; draw-wise factor mean %.3f (SD %.3f).",
                SURVEY_YEAR, annualization$factor, annualization$factor_sd))
    } else {
      c(sprintf("  Single period only (%s); no month fixed effect fitted.", MONTH_REFERENCE),
        "  Map aggregation: single-period prediction surface (no annualization; scale factor 1.000).")
    },
    "",
    "Prediction:",
    sprintf("  Map units: expected %s events per 100 camera-days.", TARGET_LABEL),
    sprintf("  Map summary method: %s; joint posterior draws: %d.", prediction$method, prediction$nsim),
    "  Spatial field, intercept, month effects, and zero-inflation probability are transformed within each joint draw.",
    "",
    "Interpretation:",
    sprintf("  Relative %s encounter frequency only.", TARGET_LABEL),
    "  Not abundance, density, occupancy, or population size.",
    if (isTRUE(settings$use_month_effect)) {
      sprintf("  Spatial surface is month-adjusted and annualized over the sampled %s months.",
              SURVEY_YEAR)
    } else {
      sprintf("  Spatial surface is for a single period (%s); no month adjustment.",
              MONTH_REFERENCE)
    }
  )

  writeLines(report, path_out(paste0(SURVEY_PREFIX, "_validation_report.txt")))
  invisible(report)
}

# Write the one-row run manifest.
write_manifest <- function(model_dat, diag, cv, prediction, temporal_diag, dependence_comparison) {
  review <- validation_review_status(diag, cv, temporal_diag, dependence_comparison,
                                     comparison_required = RUN_DEPENDENCE_COMPARISON)
  manifest <- data.frame(
    survey = SURVEY_ID,
    prefix = SURVEY_PREFIX,
    model = FINAL_MODEL_NAME,
    family = FINAL_FAMILY,
    cameras = dplyr::n_distinct(model_dat$plotID),
    model_rows = nrow(model_dat),
    events = sum(model_dat$wolf_events),
    effort_days = sum(model_dat$total_effort_days),
    month_effect = isTRUE(settings$use_month_effect),
    diagnostics_ok = isTRUE(diag$diagnostics_ok),
    spatial_cv_run = !is.null(cv),
    review_status = review$status,
    map_summary_method = prediction$method,
    map_nsim = prediction$nsim,
    map_seed = prediction$seed,
    spatial_fold_seed = SPATIAL_FOLD_SEED,
    training_only_cv_prior = TRUE,
    dependence_comparison_complete = !is.null(dependence_comparison) &&
      all(dependence_comparison$status == "completed"),
    run_profile = RUN_PROFILE
  )
  readr::write_csv(manifest, path_out(paste0(SURVEY_PREFIX, "_run_manifest.csv")))
  invisible(manifest)
}



## 15. Additional science checks: model comparison and mesh ----------------

# These checks are deliberately separated from the main fit. They are not used
# to draw the final map directly; they document whether the selected model is
# actually needed and whether the spatial result is sensitive to common modelling
# choices that reviewers usually ask about.

RUN_MODEL_COMPARISON <- tolower(Sys.getenv("WOLF_RUN_MODEL_COMPARISON", unset = ifelse(RUN_PROFILE == "quick", "false", "true"))) %in%
  c("true", "1", "yes", "y")
RUN_MESH_SENSITIVITY <- tolower(Sys.getenv("WOLF_RUN_MESH_SENSITIVITY", unset = ifelse(RUN_PROFILE == "quick", "false", "true"))) %in%
  c("true", "1", "yes", "y")

# Return a scalar value, or NA if it is missing or empty.
safe_value <- function(x, default = NA_real_) {
  if (is.null(x) || !length(x)) default else as.numeric(x[[1]])
}

# Extract WAIC, DIC, and marginal likelihood from a fit.
inla_criteria <- function(fit) {
  c(
    dic = safe_value(fit$dic$dic),
    p_dic = safe_value(fit$dic$p.eff),
    waic = safe_value(fit$waic$waic),
    p_waic = safe_value(fit$waic$p.eff),
    marginal_loglik = safe_value(fit$mlik[1])
  )
}

# Fraction of CPO/PIT values flagged as failures.
cpo_failure_rate <- function(fit) {
  if (is.null(fit$cpo)) return(NA_real_)
  vals <- fit$cpo$cpo
  fail <- fit$cpo$failure
  bad_vals <- !is.finite(vals) | vals <= 0
  bad_fail <- if (is.null(fail)) rep(FALSE, length(vals)) else fail > 0
  mean(bad_vals | bad_fail, na.rm = TRUE)
}

# Number of mesh vertices used by a fit.
mesh_n_vertices <- function(fit_result) {
  if (!is.null(fit_result$spde_obj) && !is.null(fit_result$spde_obj$mesh$n)) {
    return(as.integer(fit_result$spde_obj$mesh$n))
  }
  NA_integer_
}

# Fit an observed-data model for a given family (used by comparison and screening).
fit_observed_model_generic <- function(model_dat, settings, family,
                                       spatial = TRUE,
                                       model_label = "model") {
  family <- fit_family(family)
  dat <- model_dat %>% mutate(y = as.integer(wolf_events), intercept = 1)
  fixed_terms <- fixed_effect_terms(dat)
  formula_fixed <- paste(fixed_terms, collapse = " + ")

  if (!isTRUE(spatial)) {
    cat(sprintf("[wolf] model comparison fit: %s | family=%s | spatial=false\n",
                model_label, family))
    fit <- INLA::inla(
      as.formula(paste("y ~ 0 +", formula_fixed)),
      family = family,
      data = dat,
      E = dat$total_effort_days,
      control.predictor = list(compute = TRUE, link = 1),
      control.compute = list(config = FALSE, dic = TRUE, waic = TRUE, cpo = TRUE),
      control.fixed = make_control_fixed(fixed_terms),
      control.family = make_control_family(family),
      verbose = FALSE
    )
    return(list(fit = fit, model_label = model_label, family = family,
                spatial = FALSE, spde_obj = NULL))
  }

  cat(sprintf("[wolf] model comparison fit: %s | family=%s | spatial=true\n",
              model_label, family))
  camera_summary <- camera_summary_from_model(dat)
  camera_sf <- camera_to_utm(camera_summary)
  coords_camera <- st_coordinates(camera_sf)
  colnames(coords_camera) <- c("x", "y")

  obs_sf <- dat %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(EPSG_UTM)
  coords_obs <- st_coordinates(obs_sf)
  colnames(coords_obs) <- c("x", "y")

  spde_obj <- build_spatial(coords_camera, settings)
  A_obs <- INLA::inla.spde.make.A(spde_obj$mesh, loc = coords_obs)
  fixed_obs <- as.data.frame(dat[, fixed_terms, drop = FALSE])

  stack_obs <- INLA::inla.stack(
    tag = "obs",
    data = list(y = dat$y, e = dat$total_effort_days),
    A = list(A_obs, 1),
    effects = list(spatial = spde_obj$s_index, fixed = fixed_obs)
  )
  stack_data <- INLA::inla.stack.data(stack_obs)

  formula <- as.formula(
    paste("y ~ 0 +",
          paste(c(fixed_terms, "f(spatial, model = spde_obj$spde)"), collapse = " + "))
  )

  fit <- INLA::inla(
    formula,
    family = family,
    data = stack_data,
    E = stack_data$e,
    control.predictor = list(
      A = INLA::inla.stack.A(stack_obs),
      compute = TRUE,
      link = 1
    ),
    control.compute = list(config = FALSE, dic = TRUE, waic = TRUE, cpo = TRUE),
    control.fixed = make_control_fixed(fixed_terms),
    control.family = make_control_family(family),
    verbose = FALSE
  )

  list(fit = fit, model_label = model_label, family = family,
       spatial = TRUE, spde_obj = spde_obj)
}

# Summarise a fitted model (criteria and hyperparameters) into one row.
summarise_fitted_model <- function(fit_result, note = NA_character_) {
  fit <- fit_result$fit
  crit <- inla_criteria(fit)
  data.frame(
    model = fit_result$model_label,
    family = fit_result$family,
    spatial = isTRUE(fit_result$spatial),
    note = note,
    dic = crit[["dic"]],
    p_dic = crit[["p_dic"]],
    waic = crit[["waic"]],
    p_waic = crit[["p_waic"]],
    marginal_loglik = crit[["marginal_loglik"]],
    cpo_failure_rate = cpo_failure_rate(fit),
    zi_prob_mean = if (is_zi(fit_result$family)) hyp_point(fit, PAT_ZPROB) else NA_real_,
    nb_size_mean = if (is_nb(fit_result$family)) nb_size_point(fit) else NA_real_,
    spatial_range_mean_m = if (isTRUE(fit_result$spatial)) hyp_point(fit, PAT_RANGE) else NA_real_,
    spatial_sd_mean = if (isTRUE(fit_result$spatial)) hyp_point(fit, PAT_SIGMA) else NA_real_,
    mesh_vertices = mesh_n_vertices(fit_result),
    row.names = NULL
  )
}

# Fit and compare Poisson/NB/ZINB candidates by WAIC and DIC.
run_model_comparison <- function(model_dat, settings) {
  if (!isTRUE(RUN_MODEL_COMPARISON)) {
    writeLines("Model comparison skipped by WOLF_RUN_MODEL_COMPARISON.",
               path_out(paste0(SURVEY_PREFIX, "_MODEL_COMPARISON_SKIPPED.txt")))
    return(NULL)
  }

  cat("\n[wolf] running model comparison\n")
  specs <- data.frame(
    model = c(
      "poisson_spatial_month",
      "nb_spatial_month",
      "zinb_spatial_month"
    ),
    family = c(
      "poisson",
      "nbinomial",
      "zeroinflatednbinomial1"
    ),
    spatial = c(TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  rows <- list()
  failures <- character()
  for (i in seq_len(nrow(specs))) {
    sp <- specs[i, ]
    res <- tryCatch(
      fit_observed_model_generic(model_dat, settings, sp$family,
                                 spatial = sp$spatial,
                                 model_label = sp$model),
      error = function(e) {
        failures <<- c(failures, paste0(sp$model, ": ", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(res)) rows[[length(rows) + 1L]] <- summarise_fitted_model(res)
  }

  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(out)) {
    out <- out %>%
      mutate(delta_waic = waic - min(waic, na.rm = TRUE),
             delta_dic = dic - min(dic, na.rm = TRUE)) %>%
      arrange(waic)
  }
  readr::write_csv(out, path_out(paste0(SURVEY_PREFIX, "_model_comparison.csv")))

  report <- c(
    "Model comparison:",
    "  Purpose: compare spatial Poisson, NB, and ZINB month models to check whether overdispersion or zero inflation is needed.",
    "  Lower WAIC/DIC is better, but spatial block CV and PPC diagnostics should be given more weight than information criteria alone.",
    if (nrow(out)) capture.output(print(out, row.names = FALSE)) else "  No comparison models completed.",
    if (length(failures)) c("", "Failures:", paste0("  ", failures)) else NULL
  )
  writeLines(report, path_out(paste0(SURVEY_PREFIX, "_model_comparison_report.txt")))
  invisible(out)
}

# Build finer and coarser mesh-settings variants for the sensitivity check.
mesh_settings_variants <- function(settings) {
  list(
    current = settings,
    finer = modifyList(settings, list(
      mesh_cutoff_m = max(100, settings$mesh_cutoff_m * 0.70),
      mesh_max_edge = pmax(100, settings$mesh_max_edge * 0.70),
      mesh_offset = pmax(500, settings$mesh_offset * 0.80)
    )),
    coarser = modifyList(settings, list(
      mesh_cutoff_m = settings$mesh_cutoff_m * 1.45,
      mesh_max_edge = settings$mesh_max_edge * 1.45,
      mesh_offset = settings$mesh_offset * 1.20
    ))
  )
}

# Refit under mesh variants and report WAIC/hyperparameter stability.
run_mesh_sensitivity <- function(model_dat, settings, family) {
  if (!isTRUE(RUN_MESH_SENSITIVITY)) {
    writeLines("Mesh sensitivity skipped by profile or WOLF_RUN_MESH_SENSITIVITY.",
               path_out(paste0(SURVEY_PREFIX, "_MESH_SENSITIVITY_SKIPPED.txt")))
    return(NULL)
  }

  cat("\n[wolf] running mesh sensitivity\n")
  variants <- mesh_settings_variants(settings)
  rows <- list()
  failures <- character()

  for (nm in names(variants)) {
    stg <- variants[[nm]]
    label <- paste0("final_", nm, "_mesh")
    res <- tryCatch(
      fit_observed_model_generic(model_dat, stg, family,
                                 spatial = TRUE,
                                 model_label = label),
      error = function(e) {
        failures <<- c(failures, paste0(label, ": ", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(res)) {
      sm <- summarise_fitted_model(res, note = nm)
      sm$mesh_cutoff_m <- stg$mesh_cutoff_m
      sm$mesh_max_edge_inner_m <- stg$mesh_max_edge[1]
      sm$mesh_max_edge_outer_m <- stg$mesh_max_edge[2]
      rows[[length(rows) + 1L]] <- sm
    }
  }

  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(out)) {
    out <- out %>% mutate(delta_waic = waic - min(waic, na.rm = TRUE))
  }
  readr::write_csv(out, path_out(paste0(SURVEY_PREFIX, "_mesh_sensitivity.csv")))

  report <- c(
    "Mesh sensitivity:",
    "  Purpose: check whether spatial hyperparameters and information criteria are stable under finer/coarser SPDE meshes.",
    "  This is a diagnostic, not an automatic model-selection rule.",
    if (nrow(out)) capture.output(print(out, row.names = FALSE)) else "  No mesh variants completed.",
    if (length(failures)) c("", "Failures:", paste0("  ", failures)) else NULL
  )
  writeLines(report, path_out(paste0(SURVEY_PREFIX, "_mesh_sensitivity_report.txt")))
  invisible(out)
}

## 16. Ordered workflow helpers: exploratory checks, prior sensitivity --------

RUN_PRIOR_SENSITIVITY <- tolower(Sys.getenv("WOLF_RUN_PRIOR_SENSITIVITY", unset = ifelse(RUN_PROFILE == "quick", "false", "true"))) %in%
  c("true", "1", "yes", "y")

# Write the ordered-workflow description report.
write_workflow_order_report <- function() {
  lines <- c(
    sprintf("Ordered %s spatial-modelling workflow:", SURVEY_YEAR),
    "",
    "1. Data preparation and quality control",
    sprintf("   Load deployments and observations, check coordinates, effort, months, and independent %s-event counts.", TARGET_LABEL),
    "",
    "2. Exploratory checks",
    "   Summarise observed encounter rates, month structure, effort, raw spatial pattern, and deployment timing versus northing.",
    "",
    "3. Candidate model comparison",
    "   Fit spatial Poisson/NB/ZINB month models to check whether overdispersion and zero inflation are needed.",
    "",
    "4. Final-model decision",
    "   Keep the configured final model only if it is supported by diagnostics, information criteria, parsimony, and ecological interpretation.",
    "",
    "5. Prior-influence screen",
    "   Fit the selected model with the chosen priors and quantify which priors may be influential or in tension with the data.",
    "",
    "6. Prior sensitivity",
    "   Rerun reasonable alternatives for the priors flagged by the screen, plus core spatial and likelihood priors.",
    "",
    "7. Final model fit and mapping",
    "   Fit the configured final model with the chosen priors and predict relative encounter frequency per 100 camera-days.",
    "",
    "8. Final model diagnostics",
    "   Run posterior predictive checks, PIT, Pearson residual checks, Moran's I, semivariogram, residual covariate plots, and temporal lag/ACF-style checks.",
    "",
    "9. Robustness checks",
    "   Run spatial block cross-validation and mesh sensitivity.",
    "",
    "10. Interpretation and limitations",
    "   Report that the output is relative encounter frequency, not abundance, density, occupancy, or population size."
  )
  writeLines(lines, path_out(paste0(SURVEY_PREFIX, "_ORDERED_WORKFLOW.txt")))
  invisible(lines)
}

# Run and write the pre-modelling exploratory checks (rates, timing).
write_exploratory_checks <- function(model_dat) {
  cat("\n[wolf] writing exploratory checks\n")

  overall <- data.frame(
    cameras = dplyr::n_distinct(model_dat$plotID),
    deployments = nrow(model_dat),
    positive_deployment_rows = sum(model_dat$wolf_events > 0, na.rm = TRUE),
    zero_deployment_rows = sum(model_dat$wolf_events == 0, na.rm = TRUE),
    events = sum(model_dat$wolf_events, na.rm = TRUE),
    effort_days = sum(model_dat$total_effort_days, na.rm = TRUE),
    observed_rate_per_100_days = 100 * sum(model_dat$wolf_events, na.rm = TRUE) /
      sum(model_dat$total_effort_days, na.rm = TRUE),
    first_start = as.character(min(model_dat$start, na.rm = TRUE)),
    last_start = as.character(max(model_dat$start, na.rm = TRUE)),
    months = paste(sort(unique(model_dat$month)), collapse = ", "),
    row.names = NULL
  )
  readr::write_csv(overall, path_out(paste0(SURVEY_PREFIX, "_exploratory_overall_summary.csv")))

  by_month <- model_dat %>%
    group_by(month) %>%
    summarise(
      deployment_rows = n(),
      cameras = n_distinct(plotID),
      positive_rows = sum(wolf_events > 0, na.rm = TRUE),
      events = sum(wolf_events, na.rm = TRUE),
      effort_days = sum(total_effort_days, na.rm = TRUE),
      rate_per_100 = 100 * events / effort_days,
      .groups = "drop"
    ) %>%
    arrange(month)
  readr::write_csv(by_month, path_out(paste0(SURVEY_PREFIX, "_exploratory_month_summary.csv")))

  camera_summary <- camera_summary_from_model(model_dat %>% mutate(y = wolf_events))
  camera_sf <- camera_to_utm(camera_summary)
  camera_sf$wolf_events <- camera_summary$wolf_events
  camera_sf$wolf_events_per_100_days <- camera_summary$wolf_events_per_100_days

  raw_map <- ggplot(camera_sf) +
    geom_sf(aes(size = wolf_events_per_100_days), shape = 21,
            fill = "black", colour = "white", alpha = 0.85, stroke = 0.25) +
    scale_size_continuous(range = c(1.5, 7), labels = label_number(accuracy = 0.01)) +
    coord_sf(datum = NA) +
    labs(title = "Observed wolf encounter frequency by camera",
         subtitle = "Raw observed events per 100 camera-days before model smoothing",
         size = "observed events\n/100 days",
         x = "Easting, UTM 34N", y = "Northing, UTM 34N") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())
  ggsave(path_out(paste0(SURVEY_PREFIX, "_exploratory_raw_camera_rates_map.png")),
         raw_map, width = 7.5, height = 7, dpi = 250)

  model_dat$start_date <- as.Date(model_dat$start)
  date_plot <- ggplot(model_dat, aes(start_date, wolf_events_per_100_days)) +
    geom_point(alpha = 0.65) +
    geom_smooth(se = TRUE, method = "loess", formula = y ~ x) +
    scale_y_continuous(trans = "sqrt") +
    labs(title = "Observed camera-month-row rates through time",
         x = "deployment start date",
         y = "observed wolf events / 100 camera-days, square-root scale") +
    theme_minimal(base_size = 12)
  ggsave(path_out(paste0(SURVEY_PREFIX, "_exploratory_observed_rate_vs_date.png")),
         date_plot, width = 7, height = 5, dpi = 220)

  month_plot <- ggplot(by_month, aes(month, rate_per_100)) +
    geom_col() +
    labs(title = "Observed monthly wolf encounter frequency",
         x = "calendar month", y = "observed events / 100 camera-days") +
    theme_minimal(base_size = 12)
  ggsave(path_out(paste0(SURVEY_PREFIX, "_exploratory_month_rates.png")),
         month_plot, width = 6.5, height = 4.8, dpi = 220)

  obs_sf <- model_dat %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(EPSG_UTM)
  xy <- st_coordinates(obs_sf)
  timing <- model_dat %>%
    mutate(
      utm_x = xy[, 1],
      utm_y = xy[, 2],
      start_date = as.Date(start),
      start_doy = as.numeric(format(start, "%j"))
    )

  timing_cor <- tryCatch(
    suppressWarnings(cor.test(timing$start_doy, timing$utm_y, method = "spearman")),
    error = function(e) NULL
  )
  timing_summary <- data.frame(
    diagnostic = "deployment_start_day_of_year_vs_northing",
    spearman_rho = if (!is.null(timing_cor)) unname(timing_cor$estimate) else NA_real_,
    p_value = if (!is.null(timing_cor)) timing_cor$p.value else NA_real_,
    n = nrow(timing),
    row.names = NULL
  )
  readr::write_csv(timing_summary, path_out(paste0(SURVEY_PREFIX, "_exploratory_timing_vs_northing.csv")))

  timing_plot <- ggplot(timing, aes(start_doy, utm_y)) +
    geom_point(alpha = 0.65) +
    geom_smooth(se = TRUE, method = "lm", formula = y ~ x) +
    labs(title = "Deployment timing versus northing",
         subtitle = sprintf("Spearman rho = %.3f, p = %.3g",
                            timing_summary$spearman_rho, timing_summary$p_value),
         x = "deployment start day of year", y = "Northing, UTM 34N") +
    theme_minimal(base_size = 12)
  ggsave(path_out(paste0(SURVEY_PREFIX, "_exploratory_deployment_timing_vs_northing.png")),
         timing_plot, width = 6.8, height = 5, dpi = 220)

  report <- c(
    "Exploratory checks:",
    "",
    capture.output(print(overall, row.names = FALSE)),
    "",
    "Month summary:",
    capture.output(print(by_month, row.names = FALSE)),
    "",
    sprintf("Deployment timing vs northing: Spearman rho = %.3f, p = %.3g, n = %d.",
            timing_summary$spearman_rho,
            timing_summary$p_value,
            timing_summary$n),
    "",
    "Interpretation:",
    "  These checks are descriptive. They justify including month and a spatial field but do not replace posterior diagnostics."
  )
  writeLines(report, path_out(paste0(SURVEY_PREFIX, "_EXPLORATORY_REPORT.txt")))
  invisible(list(overall = overall, by_month = by_month, timing = timing_summary))
}

# Write the narrative model-choice justification.
write_model_choice_report <- function(model_comparison) {
  lines <- c(
    "Final model-choice report:",
    "",
    sprintf("Configured final model: %s, family = %s.", FINAL_MODEL_NAME, fit_family(FINAL_FAMILY)),
    "The configured final model is used for mapping, diagnostics, mesh sensitivity, and spatial CV.",
    "Use model comparison together with posterior predictive checks, spatial CV, residual diagnostics, prior sensitivity, and parsimony."
  )

  if (!is.null(model_comparison) && nrow(model_comparison)) {
    best <- model_comparison[order(model_comparison$waic), ][1, ]
    lines <- c(lines, "", sprintf("Best WAIC model in the comparison table: %s, WAIC = %.2f.", best$model, best$waic))

    # Everything below is stated relative to the family this survey actually
    # configured, not to a fixed assumption that NB is the configured model.
    cfg_row <- comparison_row_for_family(FINAL_FAMILY)
    cfg <- model_comparison[model_comparison$model == cfg_row, ]
    cfg_tag <- comparison_tag(cfg_row)

    if (nrow(cfg)) {
      lines <- c(lines, sprintf("Configured %s spatial-month WAIC = %.2f, DIC = %.2f.",
                                cfg_tag, cfg$waic[[1]], cfg$dic[[1]]))
      if (is.finite(cfg$zi_prob_mean[[1]])) {
        lines <- c(lines, sprintf("%s estimated zero-inflation probability mean in comparison fit: %.3f.",
                                  cfg_tag, cfg$zi_prob_mean[[1]]))
      }
      if (is.finite(cfg$nb_size_mean[[1]])) {
        lines <- c(lines, sprintf("%s estimated negative-binomial size mean in comparison fit: %.3f.",
                                  cfg_tag, cfg$nb_size_mean[[1]]))
      }
    } else {
      lines <- c(lines, sprintf(
        "Configured model (%s) has no matching row in the comparison table; deltas below are omitted.",
        FINAL_MODEL_NAME))
    }

    # Every other candidate, reported as a delta against the configured model,
    # each followed by its own fitted overdispersion / zero-inflation terms so
    # the parsimony argument below can be checked against them.
    others <- model_comparison[model_comparison$model != cfg_row, ]
    others <- others[order(others$waic), ]
    for (i in seq_len(nrow(others))) {
      o <- others[i, ]
      o_tag <- comparison_tag(as.character(o$model[[1]]))
      if (nrow(cfg) && is.finite(o$waic[[1]]) && is.finite(cfg$waic[[1]])) {
        lines <- c(lines, sprintf("%s spatial-month WAIC = %.2f; delta(%s - %s) = %.2f.",
                                  o_tag, o$waic[[1]], o_tag, cfg_tag,
                                  o$waic[[1]] - cfg$waic[[1]]))
      } else {
        lines <- c(lines, sprintf("%s spatial-month WAIC = %.2f.", o_tag, o$waic[[1]]))
      }
      if (is.finite(o$zi_prob_mean[[1]])) {
        lines <- c(lines, sprintf("%s estimated zero-inflation probability mean in comparison fit: %.3f.",
                                  o_tag, o$zi_prob_mean[[1]]))
      }
      if (is.finite(o$nb_size_mean[[1]])) {
        lines <- c(lines, sprintf("%s estimated negative-binomial size mean in comparison fit: %.3f.",
                                  o_tag, o$nb_size_mean[[1]]))
      }
    }

    # CPO failure rate for every candidate that reports one. A high rate means
    # INLA's leave-one-out approximation was unreliable for many observations,
    # so WAIC/DIC for that row should be read with caution.
    for (i in seq_len(nrow(model_comparison))) {
      m <- model_comparison[i, ]
      if (is.finite(m$cpo_failure_rate[[1]]) && m$cpo_failure_rate[[1]] > 0) {
        lines <- c(lines, sprintf("%s spatial-month CPO failure rate = %.3f.",
                                  comparison_tag(as.character(m$model[[1]])),
                                  m$cpo_failure_rate[[1]]))
      }
    }

    if (nrow(cfg)) {
      best_tag <- comparison_tag(as.character(best$model[[1]]))
      gap <- best$waic[[1]] - cfg$waic[[1]]
      if (identical(as.character(best$model[[1]]), cfg_row)) {
        lines <- c(lines, sprintf(
          "%s spatial-month is the best WAIC model in the comparison set; retain it unless PPC/CV diagnostics fail.",
          cfg_tag))
      } else if (is.finite(gap) && abs(gap) < 2) {
        # When the marginally-better model is zero-inflated, name its fitted
        # zero-inflation probability: a low value is the substantive reason the
        # extra parameter is not worth carrying.
        best_zi <- best$zi_prob_mean[[1]]
        zi_clause <- if (is.finite(best_zi) && best_zi < 0.10) {
          sprintf(" and its estimated zero-inflation probability is low (%.3f)", best_zi)
        } else {
          ""
        }
        lines <- c(lines, sprintf(
          paste0("%s is marginally lower by WAIC, but the difference is < 2 WAIC units%s;",
                 " the configured %s spatial-month model is retained for parsimony."),
          best_tag, zi_clause, cfg_tag))
      } else {
        lines <- c(lines, sprintf(
          paste("The configured %s model is not the best WAIC model (%s is lower by %.2f WAIC units);",
                "decide final use after PPC, residual diagnostics, spatial CV, and parsimony."),
          cfg_tag, best_tag, abs(gap)))
      }
    }
  } else {
    lines <- c(lines, "", "Model comparison was skipped or failed, so the configured final model remains a prior decision.")
  }

  writeLines(lines, path_out(paste0(SURVEY_PREFIX, "_model_choice_report.txt")))
  invisible(lines)
}

# Snapshot the current prior global variables.
save_prior_state <- function(settings_current) {
  list(
    settings = settings_current,
    PRIOR_INTERCEPT_MEAN = PRIOR_INTERCEPT_MEAN,
    PRIOR_INTERCEPT_PREC = PRIOR_INTERCEPT_PREC,
    PRIOR_MONTH_LOG_RATE_RATIO_PREC = PRIOR_MONTH_LOG_RATE_RATIO_PREC,
    PRIOR_ZI_LOGIT_MEAN = PRIOR_ZI_LOGIT_MEAN,
    PRIOR_ZI_LOGIT_PREC = PRIOR_ZI_LOGIT_PREC,
    PRIOR_NB_LOGSIZE_MEAN = PRIOR_NB_LOGSIZE_MEAN,
    PRIOR_NB_LOGSIZE_PREC = PRIOR_NB_LOGSIZE_PREC
  )
}

# Restore the prior global variables from a snapshot.
restore_prior_state <- function(state) {
  PRIOR_INTERCEPT_MEAN <<- state$PRIOR_INTERCEPT_MEAN
  PRIOR_INTERCEPT_PREC <<- state$PRIOR_INTERCEPT_PREC
  PRIOR_MONTH_LOG_RATE_RATIO_PREC <<- state$PRIOR_MONTH_LOG_RATE_RATIO_PREC
  PRIOR_ZI_LOGIT_MEAN <<- state$PRIOR_ZI_LOGIT_MEAN
  PRIOR_ZI_LOGIT_PREC <<- state$PRIOR_ZI_LOGIT_PREC
  PRIOR_NB_LOGSIZE_MEAN <<- state$PRIOR_NB_LOGSIZE_MEAN
  PRIOR_NB_LOGSIZE_PREC <<- state$PRIOR_NB_LOGSIZE_PREC
  invisible(TRUE)
}

# Refit under perturbed priors and report WAIC/hyperparameter stability.
run_prior_sensitivity <- function(model_dat, settings, family, observed_daily_rate) {
  if (!isTRUE(RUN_PRIOR_SENSITIVITY)) {
    writeLines("Prior sensitivity skipped by profile or WOLF_RUN_PRIOR_SENSITIVITY.",
               path_out(paste0(SURVEY_PREFIX, "_PRIOR_SENSITIVITY_SKIPPED.txt")))
    return(NULL)
  }

  cat("\n[wolf] running prior sensitivity\n")
  old <- save_prior_state(settings)
  on.exit(restore_prior_state(old), add = TRUE)

  variants <- list(
    final_current = list(
      settings = settings,
      intercept_mean = log(observed_daily_rate),
      intercept_prec = 1 / 2.5^2,
      zi_mean = qlogis(0.05),
      zi_prec = 1 / 1.5^2,
      nb_mean = log(2),
      nb_prec = 1 / 2^2,
      note = "final selected NB priors: spatial SD prior, current range, broad log-size prior"
    ),
    wider_spatial_sd_1p5 = list(
      settings = modifyList(settings, list(prior_sigma = c(1.50, 0.05))),
      intercept_mean = log(observed_daily_rate),
      intercept_prec = 1 / 2.5^2,
      zi_mean = qlogis(0.05),
      zi_prec = 1 / 1.5^2,
      nb_mean = log(2),
      nb_prec = 1 / 2^2,
      note = "wider spatial SD prior for comparison"
    ),
    tighter_original_like = list(
      settings = modifyList(settings, list(prior_sigma = c(0.85, 0.05))),
      intercept_mean = 0,
      intercept_prec = 0.01,
      zi_mean = -1,
      zi_prec = 0.2,
      nb_mean = log(2),
      nb_prec = 1 / 2^2,
      note = "tighter original-like spatial SD and broad intercept priors"
    ),
    extra_wide_spatial_sd = list(
      settings = modifyList(settings, list(prior_sigma = c(2.50, 0.05))),
      intercept_mean = log(observed_daily_rate),
      intercept_prec = 1 / 2.5^2,
      zi_mean = qlogis(0.05),
      zi_prec = 1 / 1.5^2,
      nb_mean = log(2),
      nb_prec = 1 / 2^2,
      note = "extra-wide spatial SD prior, checks robustness beyond final prior"
    ),
    shorter_spatial_range = list(
      settings = modifyList(settings, list(prior_range_m = c(2500, 0.5))),
      intercept_mean = log(observed_daily_rate),
      intercept_prec = 1 / 2.5^2,
      zi_mean = qlogis(0.05),
      zi_prec = 1 / 1.5^2,
      nb_mean = log(2),
      nb_prec = 1 / 2^2,
      note = "shorter spatial-range prior median"
    ),
    wider_spatial_range = list(
      settings = modifyList(settings, list(prior_range_m = c(10000, 0.5))),
      intercept_mean = log(observed_daily_rate),
      intercept_prec = 1 / 2.5^2,
      zi_mean = qlogis(0.05),
      zi_prec = 1 / 1.5^2,
      nb_mean = log(2),
      nb_prec = 1 / 2^2,
      note = "wider spatial-range prior median"
    )
  )

  rows <- list()
  failures <- character()
  for (nm in names(variants)) {
    v <- variants[[nm]]
    PRIOR_INTERCEPT_MEAN <<- v$intercept_mean
    PRIOR_INTERCEPT_PREC <<- v$intercept_prec
    PRIOR_ZI_LOGIT_MEAN <<- v$zi_mean
    PRIOR_ZI_LOGIT_PREC <<- v$zi_prec
    PRIOR_NB_LOGSIZE_MEAN <<- v$nb_mean
    PRIOR_NB_LOGSIZE_PREC <<- v$nb_prec

    res <- tryCatch(
      fit_observed_model_generic(model_dat, v$settings, family,
                                 spatial = TRUE,
                                 model_label = paste0("prior_", nm)),
      error = function(e) {
        failures <<- c(failures, paste0(nm, ": ", conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(res)) {
      sm <- summarise_fitted_model(res, note = v$note)
      sm$prior_variant <- nm
      sm$intercept_prior_mean <- v$intercept_mean
      sm$intercept_prior_sd <- 1 / sqrt(v$intercept_prec)
      sm$spatial_sd_prior_sigma0 <- v$settings$prior_sigma[1]
      sm$spatial_sd_prior_prob_above <- v$settings$prior_sigma[2]
      sm$spatial_range_prior_range0 <- v$settings$prior_range_m[1]
      sm$spatial_range_prior_prob_below <- v$settings$prior_range_m[2]
      sm$zi_prior_center <- plogis(v$zi_mean)
      sm$zi_prior_logit_sd <- 1 / sqrt(v$zi_prec)
      sm$nb_logsize_prior_mean <- v$nb_mean
      sm$nb_logsize_prior_sd <- 1 / sqrt(v$nb_prec)
      rows[[length(rows) + 1L]] <- sm
    }
  }

  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(out)) {
    out <- out %>%
      mutate(delta_waic = waic - min(waic, na.rm = TRUE),
             delta_dic = dic - min(dic, na.rm = TRUE)) %>%
      arrange(waic)
  }
  readr::write_csv(out, path_out(paste0(SURVEY_PREFIX, "_prior_sensitivity.csv")))

  report <- c(
    "Prior sensitivity:",
    "  Purpose: check whether the selected spatial model changes under reasonable prior choices.",
    "  This observed-data sensitivity is not the final map; it screens spatial/likelihood hyperparameters and information criteria.",
    if (nrow(out)) capture.output(print(out, row.names = FALSE)) else "  No prior-sensitivity variants completed.",
    if (length(failures)) c("", "Failures:", paste0("  ", failures)) else NULL
  )
  writeLines(report, path_out(paste0(SURVEY_PREFIX, "_prior_sensitivity_report.txt")))
  invisible(out)
}

# Write the ordered science-checks summary.
write_science_checks_summary_ordered <- function(model_comparison, prior_influence, prior_sensitivity, mesh_sensitivity) {
  best_model_line <- if (!is.null(model_comparison) && nrow(model_comparison)) {
    best <- model_comparison[order(model_comparison$waic), ][1, ]
    sprintf("  Best WAIC model among comparison set: %s (WAIC %.2f).", best$model, best$waic)
  } else {
    "  Model comparison was not available."
  }

  prior_influence_line <- if (!is.null(prior_influence) && nrow(prior_influence)) {
    flagged <- prior_influence %>% filter(isTRUE(sensitivity_priority) | sensitivity_priority == TRUE)
    if (nrow(flagged)) {
      paste0("  Prior-influence screen flagged: ", paste(flagged$parameter, collapse = ", "), ".")
    } else {
      "  Prior-influence screen did not flag any parameter at the automatic thresholds."
    }
  } else {
    "  Prior-influence screen was not available."
  }

  prior_line <- if (!is.null(prior_sensitivity) && nrow(prior_sensitivity)) {
    rng <- range(prior_sensitivity$waic, na.rm = TRUE)
    sprintf("  Prior-sensitivity WAIC range across variants: %.2f to %.2f.", rng[1], rng[2])
  } else {
    "  Prior sensitivity was not available."
  }

  mesh_line <- if (!is.null(mesh_sensitivity) && nrow(mesh_sensitivity)) {
    rng <- range(mesh_sensitivity$waic, na.rm = TRUE)
    sprintf("  Mesh-sensitivity WAIC range across variants: %.2f to %.2f.", rng[1], rng[2])
  } else {
    "  Mesh sensitivity was not available."
  }


  lines <- c(
    "Science-check summary, ordered workflow:",
    "",
    "1. Model comparison:",
    best_model_line,
    "  Use this together with PPC and spatial block CV; do not select solely by WAIC/DIC.",
    "",
    "2. Prior-influence screen:",
    prior_influence_line,
    "  Use this screen to decide which priors need explicit sensitivity reruns.",
    "",
    "3. Prior sensitivity:",
    prior_line,
    "  Stable hyperparameters and similar criteria across priors support prior robustness.",
    "",
    "4. Mesh sensitivity:",
    mesh_line,
    "  Stable spatial range/SD and similar WAIC across meshes support mesh robustness.",
    "",
    "5. Full-map prediction domain:",
    "  Final maps use the buffered convex-hull domain only. Disk-based prediction maps have been removed by design.",
    "",
    "Interpretation reminder:",
    "  These models estimate relative encounter frequency only."
  )
  writeLines(lines, path_out(paste0(SURVEY_PREFIX, "_science_checks_summary.txt")))
  invisible(lines)
}


## 17. Main run ---------------------------------------------------------------

validate_inputs()
write_workflow_order_report()

# 1. Data preparation and quality control.
family <- fit_family(FINAL_FAMILY)
model_dat <- load_survey_data(settings)

observed_daily_rate <- sum(model_dat$wolf_events, na.rm = TRUE) /
  sum(model_dat$total_effort_days, na.rm = TRUE)
if (!is.finite(observed_daily_rate) || observed_daily_rate <= 0) {
  stop("Cannot set intercept prior: observed daily rate is not positive.")
}

# Intercept prior must be set before any model fit because all model-comparison
# and diagnostic fits use the same fixed-effect prior machinery.
PRIOR_INTERCEPT_MEAN <- log(observed_daily_rate)
cat(sprintf(
  "[%s] prior setup: intercept mean %.3f = log(observed daily rate %.4f); spatial SD prior P(SD > %.2f) = %.2f; ZI candidate prior center %.3f; NB candidate log-size prior mean %.3f SD %.2f\n",
  SURVEY_PREFIX,
  PRIOR_INTERCEPT_MEAN,
  observed_daily_rate,
  settings$prior_sigma[1],
  settings$prior_sigma[2],
  plogis(PRIOR_ZI_LOGIT_MEAN),
  PRIOR_NB_LOGSIZE_MEAN,
  1 / sqrt(PRIOR_NB_LOGSIZE_PREC)
))

# 2. Exploratory checks before model selection.
exploratory <- write_exploratory_checks(model_dat)

# 3. Candidate model comparison before deciding whether NB is supported relative to Poisson/NB alternatives.
model_comparison <- run_model_comparison(model_dat, settings)

# 4. Provisional final-model decision report. The configured family is fixed in this script; this report documents whether the configured NB model is supported.
write_model_choice_report(model_comparison)

# 5. Prior-influence screen after model comparison and before sensitivity runs.
# This observed-data fit is inexpensive relative to the final mapping fit and is
# used only to decide which priors deserve explicit sensitivity checks.
cat("\n[wolf] running prior-influence screen before sensitivity reruns\n")
prior_screen_fit <- fit_observed_model_generic(
  model_dat,
  settings,
  family,
  spatial = TRUE,
  model_label = "prior_influence_screen"
)
prior_influence <- write_prior_influence_diagnostics(
  prior_screen_fit$fit,
  settings,
  family,
  label = "pre_sensitivity"
)

# 6. Prior sensitivity after the influence screen and before the final map.
prior_sensitivity <- run_prior_sensitivity(model_dat, settings, family, observed_daily_rate)

# Restore the chosen final relaxed priors after the sensitivity loop.
PRIOR_INTERCEPT_MEAN <- log(observed_daily_rate)
PRIOR_INTERCEPT_PREC <- 1 / 2.5^2
PRIOR_ZI_LOGIT_MEAN <- qlogis(0.05)  # skeptical but flexible NB prior
PRIOR_ZI_LOGIT_PREC <- 1 / 1.5^2
PRIOR_NB_LOGSIZE_MEAN <- log(2)
PRIOR_NB_LOGSIZE_PREC <- 1 / 2^2

# 7. Main final model fit and mapping fit.
fit_obj <- fit_final_model(model_dat, settings, family)

# 8. Final diagnostics from observed-data diagnostic fit with posterior samples.
diag_fit_obj <- fit_diagnostic_model(model_dat, settings, family)

cat(sprintf("\n[wolf] drawing %d posterior samples for PPC and diagnostics from observed-data fit\n", PPC_NSIM))
ppc_samples <- posterior_samples_safe(diag_fit_obj$fit, PPC_NSIM)

diag <- compute_diagnostics(
  fit = diag_fit_obj$fit,
  samples = ppc_samples,
  model_dat = diag_fit_obj$model_dat,
  obs_index = diag_fit_obj$obs_index,
  camera_sf = diag_fit_obj$camera_sf,
  family = family,
  write_files = TRUE
)

write_diagnostic_plots(diag, diag_fit_obj$camera_sf)
temporal_diag <- temporal_autocorrelation_diagnostics(diag$model_dat)
write_prior_posterior_plots(diag_fit_obj$fit, settings, family)
write_month_coefficients(diag_fit_obj$fit, diag$model_dat, settings)
write_model_hyperparameters(diag_fit_obj$fit)

prediction <- make_prediction_outputs(fit_obj, diag, settings, family)

# 9. Robustness checks after the final model is fitted.
cv <- NULL
if (RUN_SPATIAL_CV) {
  cv <- spatial_block_cv(model_dat, settings, family, K = CV_K)
}
dependence_comparison <- run_dependence_comparison(
  model_dat, settings, family, fit_obj, diag, cv, prediction, temporal_diag)
mesh_sensitivity <- run_mesh_sensitivity(model_dat, settings, family)
# 10. Final summaries.
write_science_checks_summary_ordered(model_comparison, prior_influence, prior_sensitivity, mesh_sensitivity)
write_validation_report(model_dat, diag, cv, prediction, temporal_diag, dependence_comparison)
write_manifest(model_dat, diag, cv, prediction, temporal_diag, dependence_comparison)

cat(sprintf("\n[wolf] VALIDATION\n"))
cat(sprintf("[wolf]   PPC method: %s\n", diag$ppc_method))
cat(sprintf("[wolf]   Pearson dispersion, model rows: %.3f\n", diag$pearson_disp))
cat(sprintf("[wolf]   Pearson dispersion, camera aggregates: %.3f\n", diag$pearson_disp_camera))
cat(sprintf("[wolf]   PPC pass, camera total / zero / max: %s / %s / %s\n",
            isTRUE(diag$ppc_total_pass),
            isTRUE(diag$ppc_zero_pass),
            isTRUE(diag$ppc_max_pass)))
cat(sprintf("[wolf]   residual Moran's I: %.3f (expected %.3f, two-sided p = %.3f)\n",
            diag$moran_I, diag$moran_expected, diag$moran_p))
cat(sprintf("[wolf]   row PIT KS p: %.3g | camera PIT KS p: %.3g\n",
            diag$ppc_pit_ks_row, diag$ppc_pit_ks_camera))
if (!is.null(temporal_diag) && nrow(temporal_diag$lag_summary)) {
  lag1_console <- temporal_diag$lag_summary %>% filter(lag == 1L)
  if (nrow(lag1_console)) {
    cat(sprintf("[wolf]   temporal lag-1 residual correlation: r = %.3f (p = %.3g, n pairs = %d)\n",
                lag1_console$correlation[[1]],
                lag1_console$p_value[[1]],
                lag1_console$n_pairs[[1]]))
  }
}
if (is_zi(family)) {
  cat(sprintf("[wolf]   zero-inflation probability mean: %.3f\n", diag$pi_hat))
}
if (is_nb(family)) {
  cat(sprintf("[wolf]   negative-binomial size mean: %.3f\n", diag$size_hat))
}
cat(sprintf("[wolf]   legacy PPC/Moran checks: %s\n", isTRUE(diag$diagnostics_ok)))
review <- validation_review_status(diag, cv, temporal_diag, dependence_comparison,
                                   comparison_required = RUN_DEPENDENCE_COMPARISON)
cat(sprintf("[wolf]   overall review status: %s\n", review$status))
if (length(review$reasons)) cat(paste0("[wolf]   review: ", review$reasons, "\n"))

if (!isTRUE(diag$diagnostics_ok)) {
  cat(sprintf("[wolf]   NOTE: mapped despite failing %s.\n",
              paste(diagnostic_failures(diag), collapse = "; ")))
  cat("[wolf]         Treat flagged aspects as stated limitations.\n")
}

if (!is.null(model_comparison) && nrow(model_comparison)) {
  cat(sprintf("[wolf]   model comparison written: %s\n",
              path_out(paste0(SURVEY_PREFIX, "_model_comparison.csv"))))
}
if (!is.null(prior_sensitivity) && nrow(prior_sensitivity)) {
  cat(sprintf("[wolf]   prior sensitivity written: %s\n",
              path_out(paste0(SURVEY_PREFIX, "_prior_sensitivity.csv"))))
}
if (!is.null(mesh_sensitivity) && nrow(mesh_sensitivity)) {
  cat(sprintf("[wolf]   mesh sensitivity written: %s\n",
              path_out(paste0(SURVEY_PREFIX, "_mesh_sensitivity.csv"))))
}
cat(sprintf("\nAll %s workflow steps completed in the ordered analysis sequence.\n",
            SURVEY_YEAR))
cat("Final outputs are in:\n  ", OUTPUT_DIR, "\n", sep = "")
cat(sprintf("Key files are prefixed '%s_' (diagnostic files also carry the model tag '%s_'):\n",
            SURVEY_PREFIX, FINAL_MODEL_NAME))
cat("  _ORDERED_WORKFLOW.txt, _EXPLORATORY_REPORT.txt, _model_choice_report.txt\n")
cat("  _validation_report.txt, _science_checks_summary.txt\n")
cat("  _model_comparison.csv / .txt, _prior_sensitivity.csv / .txt, _mesh_sensitivity.csv / .txt\n")
cat("  _final_predicted_events_per_100_days_mean/sd/cv.tif and _final_event_frequency_*.png\n")
cat("  <model>_posterior_predictive_check.csv, <model>_temporal_autocorrelation_report.txt\n")
cat("  _hyperparameters.csv, _month_coefficients.csv, _month_observed_summary.csv\n")
cat("  _final_spatial_block_cv_summary.csv (if CV was run), _run_manifest.csv\n")

###############################################################################

