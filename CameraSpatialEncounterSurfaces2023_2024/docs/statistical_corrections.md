# Statistical corrections and rerun requirements

These changes address four issues identified in the October 2026 code review.
The code fixes and diagnostic comparisons do not establish that the existing
survey models now have adequate calibration. The raw survey inputs are private;
the committed `results/` remain a legacy snapshot until an authorized rerun.

## 1. Joint uncertainty in mapped rates

`joint_surface_summary()` draws spatial coefficients, the intercept, month
effects, and likelihood hyperparameters from the same fitted posterior. For
each draw it projects the spatial field onto the prediction grid and computes
the effort-weighted sampled-month expected rate, including that draw's
zero-inflation probability. It then calculates mean, SD, CV, and 5th/95th
percentiles. This removes constant plug-in month/zero-inflation scaling and the
Gaussian marginal moment approximation. Chunked grid projection avoids storing
a complete grid-by-month-by-simulation array.

The target remains expected events per 100 camera-days over the sampled-month
effort distribution. Unsampled months are not estimated. These latent-mean
credible intervals do not include future observation noise. Cross-validation
predictive intervals include count noise separately.

Camera-level diagnostic variance is also computed from joint aggregate moments:
the mean of the summed conditional variances plus the variance of the summed
conditional means. Summing marginal row variances omitted covariance across
repeated observations from the same camera.

## 2. Undercoverage and dependence alternatives

Each CV run now writes pooled and fold-specific coverage, camera/row sample
counts, interval widths, and the 90% interval score (lower is better). Coverage
below 0.90 is an explicit descriptive review flag, not a formal test based on
independent camera-month rows. Coverage, scores, residual checks, and application
needs must be considered together; increasing widths alone is not success.

Balanced/final profiles compare the configured model with a camera IID random
intercept and, for multiple months, independent camera-specific AR1 sequences
sharing hyperparameters. AR1 indices preserve missing calendar months. New
held-out camera effects are unobserved, and their uncertainty is integrated in
predictive simulations. Spatial meshes are rebuilt from training cameras only.

The additional-effect SD prior uses `pc.prec`, `P(SD > 1)=0.05` on the log-rate
scale. AR1 correlation uses `pc.cor0`, `u=0.7, alpha=0.5`, allowing positive or
negative correlation with base model zero. These are explicit sensitivity-model
priors, not a finding that such effects exist. If they materially affect model
choice, repeat with reasonable alternative priors and retain that evidence.

The comparison exports WAIC, held-out log predictive density, camera interval
scores/coverage/RMSE, residual Moran's I p-values, temporal diagnostic summaries,
and spatial-map differences. Additional-effect map comparisons set the camera
deviation to zero. They assess standardized spatial components; they do not
claim to estimate the population-marginal mean for a randomly chosen new camera.

## 3. Temporal robustness claims

Neither mesh stability nor spatial CV establishes that residual temporal
dependence has no effect on a surface. A nonsignificant timing/northing
correlation does not rule out timing mechanisms. Those claims have been removed.
The IID/AR1 comparison produces actual sensitivity evidence; no candidate is
automatically promoted and no map-stability conclusion is generated in advance.

These CV results are model-development evidence. Choosing a new model from the
same folds requires another independent or nested assessment before claiming
unbiased final predictive performance. Spatial camera holdouts assess new-camera
prediction, not forecasting future months.

## 4. Strict spatial folds and training-only prior centering

Fold construction validates coordinates, clustering convergence, and co-located
cameras, uses an explicit seed, saves assignments with coordinates, and stops on
failure. It never falls back to random folds. Every dependence alternative uses
the same saved assignments. Its intercept prior center is calculated solely from
training counts and exposure and saved for each fold. An all-zero training fold
uses half an event over its own total effort to center the broad prior.

Posterior sampling controls both R and INLA/GMRFLib seeds and requests one thread.
Missing hyperparameter draws now fail explicitly instead of falling back to
point estimates or guessing the parameter scale.

## Tests

Run the dependency-free mathematical and regression checks from the repository
root:

```sh
Rscript CameraSpatialEncounterSurfaces2023_2024/tests/test_statistical_validation.R
```

The synthetic native integration test additionally requires the analysis
packages, including INLA, sf, and terra:

```sh
Rscript CameraSpatialEncounterSurfaces2023_2024/tests/test_inla_integration.R
```

The GitHub Actions workflow runs both suites. It exercises joint grid projection,
NB/ZINB user-scale posterior extraction, two spatial holdout folds, the IID/AR1
alternatives, and report writing on simulated inputs. Passing these tests checks
implementation; it does not validate the real survey's scientific conclusions.

## Private-data rerun and release

From `CameraSpatialEncounterSurfaces2023_2024/`, with the original required input
files available, run each survey in `final` profile, writing to a fresh directory:

```sh
WOLF_RUN_PROFILE=final WOLF_OUTPUT_DIR=outputs/corrected/road_2024 Rscript scripts/run_road_2024.R
WOLF_RUN_PROFILE=final WOLF_OUTPUT_DIR=outputs/corrected/human_2023 Rscript scripts/run_road_2023_human.R
WOLF_RUN_PROFILE=final WOLF_OUTPUT_DIR=outputs/corrected/human_2024 Rscript scripts/run_road_2024_human.R
```

Repeat for the other four runners before updating their maps and reported
uncertainty. Preserve the prior archive. Inspect each new validation report,
`*_dependence_model_comparison.csv`, `*_map_sensitivity.csv`,
`*_spatial_cv_calibration.csv`, `*_spatial_cv_fold_metrics.csv`, and
`*_spatial_cv_training_priors.csv`. Resolve failed comparisons and calibration or
residual flags before promoting results. If any remain, report the unresolved
limitation explicitly and retain exploratory status.

Record R/INLA versions, input hashes, seeds, and the commit used. Update the
curated files and README numbers together after reviewing the rerun. Existing
legacy reports must not be edited to pretend the corrected analyses were run.

## Primary software references

- [INLA posterior sampling and evaluating joint draws](https://www.r-inla.org/learnmore/docs/reference/posterior.sample.html)
- [INLA replicated random effects](https://www.r-inla.org/models/tools)
- [INLA AR1 correlation PC prior](https://www.r-inla.org/learnmore/docs/reference/pc.cor0.html)
