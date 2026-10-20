# Bayesian Spatial Encounter-Surface Models For 2023-2024 Wolf Camera-Trap Detections

> **Methods update; survey rerun required.** The code now calculates map means,
> SDs, CVs, and 90% credible bounds from joint posterior draws; uses saved spatial
> folds and training-only prior centering; and compares camera IID/AR1 alternatives.
> The committed `results/` and the numerical summaries below are **legacy outputs
> from the earlier implementation**. They have not been regenerated with these
> corrections and must not be presented as validation of the updated code.
> In particular, legacy map uncertainty used plug-in month/zero-inflation factors.
> See [the correction and rerun guide](docs/statistical_corrections.md).


This project contains the 2023-2024 wolf relative encounter-frequency models
from camera-trap data. This file is the single reference for the ecological
question, the methodology, the input data structure, the three primary wolf models
(plus four companion surfaces) with their legacy diagnostic numbers,
and the outputs needed to reproduce or audit the analysis. Each `results/` subfolder has its own short README that only
lists and describes the files it contains; all explanation lives here.

The analyses model the number of independent wolf event IDs recorded in
camera-month rows. The statistical approach is Bayesian count modelling with
INLA-SPDE spatial random fields. Active camera-days are used as an exposure
term, calendar month is included as a fixed temporal control, and outputs are
relative encounter-frequency surfaces expressed as expected wolf events per 100
camera-days across the sampled survey-year period. The maps should not be
interpreted as abundance, density, occupancy, or population size.

The study's primary target is the wolf. Companion analyses use the same
pipeline: a human-activity index (people and vehicles) on the two road-camera
surveys, giving a relative human-disturbance surface
([Human-Activity Companion Surfaces](#human-activity-companion-surfaces)); and
single-month, spatial-only surfaces for March 2024 -- the only spring maps in
the repository -- for both
[wolf](#single-month-spatial-surface-wolf-march-2024) and
[human activity](#single-month-spatial-surface-human-activity-march-2024).

All survey numbers below describe earlier `WOLF_RUN_PROFILE=final` runs against
the private camera-trap data and match the legacy files under `results/`.

> **Refactor note (shared-library consolidation).** The three surveys now share
> one analysis library, `scripts/wolf_encounter_surface_lib.R`, run through thin
> per-survey runners (`run_road_2023.R`, `run_road_2024.R`, `run_forest_2024.R`),
> so every survey runs the identical analyses, diagnostics, and outputs. The
> original consolidation preserved the road-camera 2023 and 2024 behaviour.
> The subsequent statistical corrections intentionally change uncertainty and
> validation outputs, so the current library does not reproduce the legacy
> archive unchanged. The forest-camera
> 2024 survey now runs that same pipeline, which **added the candidate-likelihood
> model-comparison table it previously lacked** and standardized its
> prior-sensitivity set to the shared six variants. `results/forest_2024/` and
> the forest numbers in this document have been regenerated from a
> `WOLF_RUN_PROFILE=final` run of the unified pipeline; the mapping-fit
> hyperparameters and month effects are unchanged from before, while the
> diagnostic-refit figures (dispersion, PIT, temporal correlation) reflect the
> shared road-style diagnostic pipeline. See
> [Reproduction And Validation](#reproduction-and-validation) below.

## Final Models

Three wolf analyses are included, plus two human-activity companion surfaces
(people and vehicles on the road cameras; see
[Human-Activity Companion Surfaces](#human-activity-companion-surfaces)):

| Survey | Final model | Cameras | Events | Effort | Final output |
| --- | --- | ---: | ---: | ---: | --- |
| Road-camera 2023 | Negative-binomial spatial-month INLA-SPDE model | 60 | 586 | 5222.2 camera-days | `results/road_2023/` |
| Forest-camera 2024 | Negative-binomial spatial-month INLA-SPDE model | 53 | 46 | 4423.0 camera-days | `results/forest_2024/` |
| Road-camera 2024 | Zero-inflated negative-binomial spatial-month INLA-SPDE model | 60 | 479 | 3574.0 camera-days | `results/road_2024/` |
| Human-activity 2023 *(companion)* | Negative-binomial spatial-month INLA-SPDE model | 60 | 8284 | 5222.2 camera-days | `results/human_2023/` |
| Human-activity 2024 *(companion)* | Negative-binomial spatial-month INLA-SPDE model | 60 | 6781 | 3574.0 camera-days | `results/human_2024/` |
| Wolf March 2024 *(companion)* | Poisson spatial-only INLA-SPDE model (single month) | 60 | 303 | 1783.4 camera-days | `results/wolf_march_2024/` |
| Human-activity March 2024 *(companion)* | Poisson spatial-only INLA-SPDE model (single month) | 60 | 2805 | 1783.4 camera-days | `results/human_march_2024/` |

Legacy diagnostic status (all surveys require corrected reruns; full historical
numbers are in each survey's section below):

| Survey | Likelihood | Required diagnostics | Open caveat |
| --- | --- | --- | --- |
| Road-camera 2023 | Negative-binomial | Pass | None |
| Forest-camera 2024 | Negative-binomial | Pass | Only 46 independent events, so posterior uncertainty on month/spatial effects is wide |
| Road-camera 2024 | Zero-inflated negative-binomial | Legacy gate passes; further review required | Residual temporal association; its effect on the surface has not been established. Compare camera IID/AR1 alternatives and map differences |
| Human-activity 2023 *(companion)* | Negative-binomial | Legacy gate passes; calibration review required | Legacy 90% camera-total interval coverage is 77%; inspect fold results and dependence alternatives |
| Human-activity 2024 *(companion)* | Negative-binomial | Review required | Legacy camera-total 90% interval coverage is 80%, and Moran p = 0.017; model adequacy remains unresolved |
| Wolf March 2024 *(companion)* | Poisson | Pass | Single month, spatial-only (no month effect); the only spring wolf surface. Model comparison prefers Poisson (no overdispersion) |
| Human-activity March 2024 *(companion)* | Poisson | Pass | Single-month spatial-only human-activity surface; model comparison prefers Poisson |

For the five multi-month surveys, the legacy WAIC comparison favors NB or ZINB
over Poisson. The two single-month March 2024 comparisons favor Poisson.
This contrast does not establish which process causes overdispersion; a single
month can also contain overdispersed counts. For
road-camera 2023 specifically, a ZINB variant scored a fractionally lower
(better) WAIC than NB (delta WAIC 0.25), but its estimated zero-inflation
probability was small (3.2%) and not mechanistically justified, so the
simpler NB model was kept. See each survey's section below for the
comparison tables.

### Reading The `cpo_failure_rate` Column

Every `_model_comparison.csv` carries a `cpo_failure_rate` column. It is the
fraction of observations for which INLA flagged its leave-one-out
conditional predictive ordinate as unreliable, and it is a numerical
diagnostic of the approximation, not a measure of fit. A high rate means that
row's CPO-based leave-one-out quantities need numerical investigation. CPO
failure does not by itself establish that WAIC or DIC is unreliable; those
criteria summarize different quantities. INLA's
[`inla.cpo()`](https://www.r-inla.org/learnmore/docs/reference/cpo.html)
can recompute flagged CPO/PIT values by leaving out observations explicitly.

In the five multi-month surveys the selected NB or ZINB model has a rate of
0, and only the rejected Poisson candidate shows failures. The two
single-month March 2024 surveys are the exception: there Poisson is both the
selected model and the one with failures (wolf 0.690, human-activity 0.207).
Those two surfaces are retained on the strength of their WAIC margin
(5.6 and 20.5 units clear of NB), their posterior predictive checks, and
their spatial cross-validation, all of which are reported in the same
folder. Their CPO failures remain a separate unresolved numerical limitation;
WAIC comparisons also require the usual model and diagnostic assessment.

### Comparing Across Years

The 2023 and 2024 road-camera surfaces use **different likelihoods** --
negative-binomial for 2023, zero-inflated negative-binomial for 2024 -- each
selected by the same rule applied to that year's own data (2023 kept NB for
parsimony at delta WAIC 0.25; 2024 took ZINB at delta WAIC 3.67, beyond the
2-unit parsimony margin). The two surfaces use different observation models,
and their adequacy must be assessed separately after the corrected reruns.
For year-on-year comparison, prefer within-year standardized values or
relative spatial pattern over differences in raw events per 100 camera-days.
The same caution applies to the March 2024 single-month surfaces, which are
Poisson and cover one month rather than an annualized period.

## Result Maps

Each survey below is mapped as three matching GeoTIFF/PNG surfaces, all
effort-weighted and annualized over that survey's sampled months (see
[Annualized Map Surface](#annualized-map-surface) below) -- not abundance,
density, occupancy, or population size.

The updated code evaluates the expected encounter rate within every joint
posterior draw, projecting sampled SPDE mesh coefficients onto the grid:

```text
rate_draw(s) = 100 * (1 - pi_draw) * exp(intercept_draw + spatial_draw(s))
               * sum_month(weight_month * exp(month_effect_draw))
```

The reference month has effect zero. Weights are the survey's fixed observed
camera-effort proportions; unsampled months receive no inferred annual coverage.
Single-month analyses use a factor of one. Spatial, intercept, month, and
zero-inflation draws come from the **same fit and same joint draw**. Means, SDs,
CVs, and 5th/95th percentiles are then calculated from the transformed rates.
Grid projection is processed in chunks to limit memory use. These are credible
summaries of the latent mean, without future count noise. Predictive count
intervals are generated separately for cross-validation.

The legacy results below instead use Gaussian marginal moment formulas and
fixed month/zero-inflation multipliers. Their SD/CV values omit uncertainty in
those multipliers and require regeneration before interpretation with the new method.

**Display note.** The PNG figures are 9.5 x 9 inches at 150 dpi, sized for
on-screen reading and ordinary print: they are drawn with bilinear
interpolation between the prediction-grid cells (150 m for the six
road-camera surfaces, 60 m for forest-camera 2024;
`geom_raster(interpolate = TRUE)`), and each map's color scale is
capped at that surface's own 98th percentile, so a handful of extreme cells
do not dominate the palette. Every cell above the cap is drawn in the top
color, which means the PNGs deliberately understate their peaks: the
road-camera 2023 mean surface, for example, reaches 88.8 events per 100
camera-days while its color scale stops at 29.5. These are display-only
steps and never modify the underlying values. The GeoTIFF raster in each
results folder holds the raw, uncapped model output and is the file to use
for any quantitative or GIS work. Opening a GeoTIFF in a GIS with a default
min-max stretch will therefore look much darker and flatter than the PNG; a
2-98 percent stretch reproduces the published appearance.

The same mechanism governs both the mean map's behavior between cameras and
the SD/CV maps' uncertainty: the spatial random field `u(s)` is a
stationary Gaussian field with an estimated **range**, the distance over
which nearby cameras' data pull it away from zero. Within roughly one
range-length of a camera, `u(s)` is informed by that camera's data and the
mean estimate can sit meaningfully above or below the population baseline.
Beyond that distance, `u(s)` reverts toward its prior mean of zero, so the
mean estimate reverts toward the fixed-effect baseline (intercept plus
month effect) -- not because the model "gives up" at an arbitrary radius,
but because there is no data left nearby to inform it. The same reversion
shows up as rising, then saturating, SD and CV: once the field's
contribution to `eta_sd(s)` reaches its marginal (prior) variance, both are
at their local ceiling. That ceiling sits slightly above the field-only
value, because `eta_sd(s)` also carries the posterior uncertainty of the
intercept and month effects: forest-camera 2024 saturates near CV 1.28
against a field-only 0.97, road-camera 2024 near 1.36 against 1.23.

This range is an estimated hyperparameter, not a fixed radius we chose, and
it differs by survey:

| Survey | Spatial range posterior mean | 95% credible interval |
| --- | ---: | ---: |
| Road-camera 2023 | 2965 m | 1337 to 5394 m |
| Road-camera 2024 | 4175 m | 2007 to 7618 m |
| Forest-camera 2024 | 584 m | 150 to 1599 m |

Forest-camera 2024's range is roughly 5-7x shorter than the two road
surveys', so its mean map falls back to baseline much closer to each
camera, and its CV map saturates to its ceiling almost everywhere outside
the immediate camera clusters (see that survey's map above). The range
itself is also fairly uncertain in all three surveys -- forest's 95%
interval spans 150 to 1599 m -- which is a direct consequence of having
comparatively few independent wolf events to estimate it from.

### Road-Camera 2023

| Posterior mean | Posterior SD | Posterior CV |
| --- | --- | --- |
| ![Road-camera 2023 posterior mean encounter frequency](results/road_2023/wolf_2023_final_event_frequency_mean.png) | ![Road-camera 2023 posterior SD](results/road_2023/wolf_2023_final_event_frequency_sd.png) | ![Road-camera 2023 posterior CV](results/road_2023/wolf_2023_final_event_frequency_cv.png) |

### Forest-Camera 2024

| Posterior mean | Posterior SD | Posterior CV |
| --- | --- | --- |
| ![Forest-camera 2024 posterior mean encounter frequency](results/forest_2024/wolf_forest_2024_final_event_frequency_mean.png) | ![Forest-camera 2024 posterior SD](results/forest_2024/wolf_forest_2024_final_event_frequency_sd.png) | ![Forest-camera 2024 posterior CV](results/forest_2024/wolf_forest_2024_final_event_frequency_cv.png) |

### Road-Camera 2024

| Posterior mean | Posterior SD | Posterior CV |
| --- | --- | --- |
| ![Road-camera 2024 posterior mean encounter frequency](results/road_2024/wolf_2024_final_event_frequency_mean.png) | ![Road-camera 2024 posterior SD](results/road_2024/wolf_2024_final_event_frequency_sd.png) | ![Road-camera 2024 posterior CV](results/road_2024/wolf_2024_final_event_frequency_cv.png) |

## Common Modelling Target

The response is an independent wolf event count per camera-month row. The
exposure is camera effort in camera-days, passed to INLA through the `E`
argument, so the linear predictor describes the expected daily encounter
rate and map outputs are converted to expected events per 100 camera-days.

The models estimate relative encounter frequency, not abundance, density,
occupancy, or population size. Camera-trap encounter rates are a standard
relative-abundance index in ecology when detection probability cannot be
estimated directly (Rowcliffe et al. 2008; O'Brien 2011).

For camera-month row `i`, the shared model structure is:

```text
log(mu_i) = log(E_i) + beta_0 + gamma[m_i] + u(s_i)
```

where `mu_i` is the expected event count, `E_i` is active camera-days,
`beta_0` is the model intercept on the log encounter-rate scale, `gamma[m_i]`
is the fixed effect for calendar month, and `u(s_i)` is a spatial random
field estimated with the INLA-SPDE method: the field is represented on a
triangulated mesh over the study area and given a Matérn covariance via a
stochastic partial differential equation, fitted by integrated nested
Laplace approximation rather than MCMC (Rue, Martino & Chopin 2009; Lindgren,
Rue & Lindström 2011). Two of the three models use a negative-binomial
likelihood for `y_i` to absorb overdispersion beyond what a Poisson count
model allows (Hilbe 2011); the road-camera 2024 model additionally uses a
zero-inflated negative-binomial likelihood (see that survey's section below).

## Data Units: Camera-Month Rows

All final models are fit on camera-month rows, not raw per-day or
per-deployment records:

1. A camera deployment is split at every calendar-month boundary it crosses.
2. Each resulting row's exposure is the number of active camera-days inside
   that specific camera-month interval.
3. Wolf events are assigned to the row whose month contains the event's
   `eventStart` timestamp.

For example, a camera deployed 2023-07-20 to 2023-08-15 produces two rows:
one for July (12 active days, 2023-07-20 to 2023-07-31) and one for August
(15 active days, 2023-08-01 to 2023-08-15). Any wolf event recorded in that
window is assigned to whichever row's month contains it. This keeps effort
and events aligned and avoids mixing two different months' exposure into one
row.

The calendar month itself enters the model as a fixed effect (`gamma[m_i]`
above), estimated jointly with the spatial field, not fit separately.

## Annualized Map Surface

Month is a fixed effect in the model, but the published maps are not a
single-month prediction. Each survey's sampled months are combined into one
effort-weighted, annualized surface:

```text
lambda_year(s) = sum_m w_m * 100 * exp(beta_0 + gamma[m] + u(s))
```

`w_m` is the share of that survey's total camera-days that fell in month
`m`. Equivalently, the reference-month prediction is scaled by an
annualization factor: the effort-weighted average of `exp(gamma[m] -
gamma[m_ref])` over all sampled months. Concretely, for road-camera 2023
(reference month 2023-08):

| Month | Share of camera-days | Rate ratio vs. 2023-08 |
| --- | ---: | ---: |
| 2023-07 | 0.296 | 0.969 |
| 2023-08 (reference) | 0.341 | 1.000 |
| 2023-09 | 0.334 | 1.371 |
| 2023-10 | 0.030 | 1.078 |

The effort-weighted average of the rate-ratio column is 1.117 -- the
annualization factor reported for that survey below. It means the
sampled-year average rate runs about 11.7% above the August-only rate,
mostly because September (33% of camera-days) has a higher fitted rate. Each
survey's full weight table is in
`results/<survey>/wolf_<survey>_annualization_weights.csv`.

## Diagnostic Gate

A model is called final only if it clears a specific gate: the camera-level
posterior predictive checks (total events, zero fraction, max count; Gelman,
Meng & Stern 1996) and the residual Moran's I spatial-autocorrelation test
(a two-sided permutation test, scaled by `WOLF_RUN_PROFILE`, applied to model
residuals following Dormann et al. 2007; Moran 1950). Every survey runs this
identical gate from the shared analysis library.

PIT (probability integral transform) KS p-values (Czado, Gneiting & Held
2009, for discrete/count outcomes) are also computed and reported for every
model as a supporting calibration diagnostic, but are not part of the gate.
Camera-level PIT is sensitive to small-count discreteness and camera-level
clustering in ways that do not track whether the mapped spatial surface
itself is distorted, so a low PIT KS p-value is reported but does not on its
own fail a model. For example, the road-camera 2023 model reports "required
diagnostics pass: TRUE" alongside a camera PIT KS p-value near 0.0015: the
low value is shown, not dropped, but does not gate the pass/fail call.

The updated validation report distinguishes the legacy PPC/Moran gate from an
`overall review status`. Held-out coverage below the nominal 90% is an explicit
descriptive review flag; fold-specific coverage, sample counts, interval widths,
and interval scores accompany it. This flag is not a significance test treating
correlated camera-month observations as independent.

Spatial folds are constructed once, saved with camera coordinates, and reused by
all dependence alternatives. A clustering failure stops validation. Each fold's
broad intercept prior is centered using **training outcomes/exposure only**; the
saved training-prior table makes this auditable. Fixed month designs may use known
held-out calendar labels, but never their responses.

Model comparison across candidate likelihoods (Poisson / NB / ZINB) uses
WAIC (Watanabe 2010), cross-checked against DIC (Spiegelhalter et al. 2002).
Out-of-sample predictive performance is checked with spatial block
cross-validation: cameras are grouped into spatial folds by k-means
clustering, each fold is held out in turn, and the SPDE mesh for that fold
is rebuilt from the training cameras only, so no held-out location leaks
into the training mesh (following the general spatially-blocked
cross-validation approach of Roberts et al. 2017). Held-out counts are
simulated from full joint posterior draws of the fitted model.

For every survey, this posterior-predictive step refits a separate
diagnostic model with `config = TRUE` so that joint posterior samples can be
drawn. Two kinds of number come out of that refit and they are not
interchangeable: `hyperparameters.csv` and `month_coefficients.csv` hold its
marginal posterior summaries, while each `validation_report.txt` prints the
mean of the joint posterior *draws* (for example, road-camera 2024
negative-binomial size 3.62 from the draws vs. 3.30 in the marginal summary).
The mapping fit is a third INLA call, on the same data and model but with the
prediction stack attached; it supplies the linear predictor behind the
published surfaces and agrees with the diagnostic fit to roughly five
significant figures (compare `month_coefficients.csv`, written from the
diagnostic fit, with `annualization_weights.csv`, written from the mapping
fit). The values quoted throughout this document are the marginal summaries
in `hyperparameters.csv`.

**Sensitivity checks.** For every survey, the shared library refits the final
model under perturbed priors (six variants) and perturbed SPDE mesh resolution
(finer/coarser) and reports whether WAIC, DIC, and posterior hyperparameters
stay stable. Because all surveys run the same code, these checks are identical
in scope across surveys: they test WAIC/DIC/hyperparameter stability rather
than recomputing the full PPC/Moran's I gate at each variant.

**Prior interpretation.** The intercept and NB-size priors use SD 2.5 and 2
on the log scale, month log-rate ratios use SD 1, and spatial range/SD use
PC priors (Fuglstad et al. 2019; Simpson et al. 2017). Their influence depends
on the parameter scale and the information in each survey; a broad or PC prior
can still materially affect a sparse-data posterior. The six legacy
perturbations assess stability of selected criteria and hyperparameters, not
prior irrelevance or stability under every plausible prior. This matters
especially for forest-camera 2024, with only 46 events. Corrected reruns must
review the sensitivity evidence together with mapped uncertainty and predictive
checks before drawing a robustness conclusion.

## Repository Layout

```text
CameraSpatialEncounterSurfaces2023_2024/
  README.md
  data/
    README.md
  results/
    README.md
    road_2023/
    forest_2024/
    road_2024/
    human_2023/                    # companion: human-activity surface
    human_2024/                    # companion: human-activity surface
    wolf_march_2024/               # companion: single-month spatial-only wolf surface
    human_march_2024/              # companion: single-month spatial-only human surface
  scripts/
    wolf_encounter_surface_lib.R   # shared analysis library: all workflow logic
    run_road_2023.R                # runner: road-camera 2023 (negative-binomial)
    run_road_2024.R                # runner: road-camera 2024 (zero-inflated NB)
    run_forest_2024.R              # runner: forest-camera 2024 (negative-binomial)
    run_road_2023_human.R          # runner: 2023 human-activity companion surface
    run_road_2024_human.R          # runner: 2024 human-activity companion surface
    run_road_march2024_wolf.R      # runner: March 2024 wolf, single-month spatial-only
    run_road_march2024_human.R     # runner: March 2024 human-activity, single-month spatial-only
    capture_session_info.R         # records R/INLA versions for reproducibility
```

Every survey runs the **same** analyses, diagnostics, and outputs, because they
all share one analysis library (`wolf_encounter_surface_lib.R`). Each `run_*.R`
runner is a thin entry point that only sets its survey's likelihood family,
priors, mesh settings, and data paths, then sources the library. Surveys differ
only in that configuration and in the shape of their raw input data (the two
road surveys read separate deployment/observation tables; the forest survey
reads a single dated flat file), never in the analysis code itself. The same
library also produces the four companion surfaces -- two multi-month
[human-activity](#human-activity-companion-surfaces) surfaces and two
single-month March 2024 surfaces
([wolf](#single-month-spatial-surface-wolf-march-2024) and
[human activity](#single-month-spatial-surface-human-activity-march-2024)) --
by changing only the configuration (detection labels, likelihood, and whether a
month fixed effect is fitted).

Raw camera-trap CSV files are not committed here. The `data/README.md` file
lists the expected input files and where to place them for reproduction.
Each `results/<survey>/` folder has its own README listing exactly what each
committed file contains; this file is where all explanation of methodology
and numbers lives.

## Main Scripts

Each survey is run through its thin runner, which sources the shared
`wolf_encounter_surface_lib.R` and runs the identical ordered workflow.

### Road-Camera 2023 Final Model

```sh
Rscript scripts/run_road_2023.R
```

Final model: negative-binomial likelihood, INLA-SPDE spatial component,
calendar-month fixed effects, active camera-days as exposure, effort-weighted
annualized 2023 map target. Full priors and diagnostics are in
[Road-Camera 2023 Model](#road-camera-2023-model) below.

### Forest-Camera 2024 Final Model

```sh
Rscript scripts/run_forest_2024.R
```

Final model: negative-binomial likelihood, INLA-SPDE spatial component,
calendar-month fixed effects, active camera-days as exposure, effort-weighted
annualized 2024 map target. Full priors and diagnostics are in
[Forest-Camera 2024 Model](#forest-camera-2024-model) below.

### Road-Camera 2024 Final Model

```sh
Rscript scripts/run_road_2024.R
```

Final model: zero-inflated negative-binomial (type 1) likelihood, INLA-SPDE
spatial component, calendar-month fixed effects, active camera-days as
exposure, effort-weighted annualized 2024 map target. Full priors and
diagnostics are in [Road-Camera 2024 Model](#road-camera-2024-model) below.

## Runtime Profiles

All workflows use the `WOLF_RUN_PROFILE` environment variable:

```sh
# Fast development run
WOLF_RUN_PROFILE=quick

# Recommended reproducible analysis
WOLF_RUN_PROFILE=balanced

# Heavier final run with more posterior simulations
WOLF_RUN_PROFILE=final
```

On Windows PowerShell, for example:

```powershell
$env:WOLF_RUN_PROFILE = "balanced"
& "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" scripts\run_road_2024.R
```

Path overrides are available:

```powershell
$env:WOLF_PROJECT_DIR = "C:\path\to\CameraSpatialEncounterSurfaces2023_2024"
$env:WOLF_DATA_DIR = "C:\path\to\CameraSpatialEncounterSurfaces2023_2024\data"
$env:WOLF_OUTPUT_DIR = "C:\path\to\outputs"
```

## Road-Camera 2023 Model

Final script:

```text
scripts/run_road_2023.R
```

Final output folder:

```text
results/road_2023/
```

Model:

```text
y_i ~ NegativeBinomial(mu_i, size)
log(mu_i) = log(effort_i) + intercept + month_i + spatial(s_i)
```

Key settings:

- 490 camera-month rows;
- 60 cameras;
- 586 independent wolf events;
- 5222.2 camera-days;
- observed mean encounter frequency: 11.221 events per 100 camera-days;
- map target: effort-weighted annualized 2023 surface (annualization factor
  1.117, see [Annualized Map Surface](#annualized-map-surface) above).

Weakly informative priors:

- intercept: Gaussian(mean = -2.187, SD 2.5 on log scale), centered on the
  crude observed daily rate;
- month log-rate ratios: Gaussian(0, SD 1);
- negative-binomial log(size): Gaussian(log(2), SD 2);
- spatial range: PC prior, `P(range < 5000 m) = 0.5` (Fuglstad et al. 2019);
- spatial marginal SD: PC prior, `P(SD > 2.0) = 0.05` (Simpson et al. 2017).

Fitted hyperparameters (from the mapping fit,
`results/road_2023/wolf_2023_hyperparameters.csv`):

- negative-binomial size posterior mean: 1.733;
- spatial range posterior mean: 2965 m (95% CrI 1337 to 5394 m);
- spatial SD posterior mean: 1.184 (95% CrI 0.890 to 1.541).

Model comparison:

| Model | WAIC | Delta WAIC |
| --- | ---: | ---: |
| ZINB spatial-month | 1162.65 | 0.00 |
| NB spatial-month | 1162.90 | 0.25 |
| Poisson spatial-month | 1303.06 | 140.41 |

ZINB is only marginally lower by WAIC and has a low estimated zero-inflation
probability (mean 0.032), so the negative-binomial spatial-month model is
retained for parsimony.

Main diagnostics:

- posterior predictive row and camera total events / zero fraction /
  maximum count: all pass;
- row Pearson dispersion: 0.661; camera Pearson dispersion: 0.286;
- residual Moran's I: -0.004 (expected -0.017), two-sided p = 0.494;
- row PIT KS p-value: 0.2056; camera PIT KS p-value: 0.0006432
  (supporting diagnostic, not part of the gate -- see
  [Diagnostic Gate](#diagnostic-gate) above);
- required diagnostics pass: TRUE;
- temporal residual autocorrelation: within-camera lag-1 r = 0.016,
  p = 0.7345 (n = 430 pairs); no evidence of residual temporal
  autocorrelation; date-ordered mean-residual lag-1 ACF: -0.162;
- spatial block cross-validation: row 90 percent coverage = 0.96, camera 90
  percent coverage = 0.95;
- prior sensitivity: WAIC, DIC, and posterior hyperparameters (NB size,
  spatial range, spatial SD) are stable across 6 prior variants (WAIC 1162.80
  to 1163.56; delta WAIC 0.00 to 0.76; stability checked, gate not
  recomputed per variant -- see "Sensitivity checks" above);
- mesh sensitivity: WAIC and hyperparameters are stable across the final,
  finer, and coarser mesh variants (WAIC 1162.81 to 1163.05; delta WAIC 0.00
  to 0.24; same basis as above).

Exploratory check: a Spearman correlation of deployment start day-of-year
against UTM northing (`results/road_2023/wolf_2023_exploratory_timing_vs_northing.csv`)
gives rho = 0.042, p = 0.357 (n = 490) -- deployment timing is not
meaningfully correlated with camera location in this survey.

## Forest-Camera 2024 Model

Final script:

```text
scripts/run_forest_2024.R
```

Final output folder:

```text
results/forest_2024/
```

Model:

```text
y_i ~ NegativeBinomial(mu_i, size)
log(mu_i) = log(effort_i) + intercept + month_i + spatial(s_i)
```

Key settings:

- 356 camera-month rows;
- 53 cameras;
- 46 independent wolf events;
- 4423.0 camera-days;
- observed mean encounter frequency: 1.040 events per 100 camera-days;
- map target: effort-weighted annualized 2024 surface (annualization factor
  1.035).

Weakly informative priors:

- intercept: Gaussian(mean = -4.566, SD 2.5 on log scale), centered on the
  crude observed daily rate;
- month log-rate ratios: Gaussian(0, SD 1);
- negative-binomial log(size): Gaussian(log(2), SD 2);
- spatial range: PC prior, `P(range < 1000 m) = 0.5`;
- spatial marginal SD: PC prior, `P(SD > 1.5) = 0.05`.

Fitted hyperparameters:

- negative-binomial size posterior mean: 1.836;
- spatial range posterior mean: 584.1 m (95% CrI 150.1 to 1598.9 m);
- spatial SD posterior mean: 0.816 (95% CrI 0.389 to 1.454).

Month effects (rate ratio vs. reference month 2024-08):

- 2024-03: 0.34 (95% CrI 0.08 to 1.47);
- 2024-04: 0.65 (95% CrI 0.22 to 1.92);
- 2024-05: 0.91 (95% CrI 0.36 to 2.31);
- 2024-06: 1.69 (95% CrI 0.74 to 3.79);
- 2024-07: 1.28 (95% CrI 0.53 to 3.07);
- 2024-09: 0.67 (95% CrI 0.23 to 1.99).

Model comparison:

| Model | WAIC | Delta WAIC |
| --- | ---: | ---: |
| NB spatial-month | 270.22 | 0.00 |
| ZINB spatial-month | 270.51 | 0.29 |
| Poisson spatial-month | 274.02 | 3.81 |

The negative-binomial model has the best (lowest) WAIC. A ZINB variant is only
marginally behind (delta WAIC 0.29) with a small estimated zero-inflation
probability (mean 0.093), and Poisson is clearly rejected (delta WAIC 3.81), so
the negative-binomial spatial-month model is retained.

Main diagnostics:

- posterior predictive row and camera total events / zero fraction /
  maximum count: all pass;
- row Pearson dispersion: 0.420; camera Pearson dispersion: 0.410;
- residual Moran's I: -0.041 (expected -0.019), two-sided p = 0.590;
- row PIT KS p-value: 0.3918; camera PIT KS p-value: 0.08281;
- required diagnostics pass: TRUE;
- temporal residual autocorrelation: within-camera lag-1 r = 0.004,
  p = 0.9412 (n = 303 pairs); date-ordered mean-residual lag-1 ACF: -0.009;
  no evidence of residual temporal autocorrelation;
- spatial block cross-validation: row 90 percent coverage = 0.98, camera 90
  percent coverage = 0.92. Each fold's SPDE mesh is rebuilt from the
  training-fold cameras only, and held-out counts are simulated from full
  joint posterior samples;
- prior sensitivity: WAIC, DIC, and posterior hyperparameters stay stable
  across the shared six prior variants (WAIC 269.54 to 277.26);
- mesh sensitivity: final, finer, and coarser mesh variants; WAIC range =
  269.49 to 270.21.

Main limitation: the forest-camera dataset contains only 46 independent wolf
events. Weakly informative priors do not add information the data lack, so
with this few events the posterior for month and spatial effects stays wide
-- the wide credible intervals above reflect data sparsity, not the priors
being informative. The model remains valid for relative encounter-frequency
mapping, but month and spatial effects should be read with that wide
uncertainty in mind.

## Road-Camera 2024 Model

Final script:

```text
scripts/run_road_2024.R
```

Final output folder:

```text
results/road_2024/
```

Model:

```text
y_i ~ ZeroInflatedNegativeBinomial1(mu_i, size, pi)
log(mu_i) = log(effort_i) + intercept + month_i + spatial(s_i)
```

`ZeroInflatedNegativeBinomial1` is INLA's type-1 zero-inflated
negative-binomial parameterization (see Martin et al. 2005 for the ecological
rationale for modelling structural zeros separately from count-process
zeros). `pi` is the additional structural-zero probability, while `mu_i` and
`size` describe the negative-binomial count component.

Key settings:

- 344 camera-month rows;
- 60 cameras;
- 479 independent wolf events;
- 3574.0 camera-days;
- observed mean encounter frequency: 13.402 events per 100 camera-days;
- map target: effort-weighted annualized 2024 surface (annualization factor
  1.195).

Weakly informative priors:

- intercept: Gaussian(mean = -2.010, SD 2.5 on log scale), centered on the
  crude observed daily rate;
- month log-rate ratios: Gaussian(0, SD 1);
- zero-inflation logit probability: Gaussian(mean = -2.94, SD approximately
  1.5 on logit scale);
- negative-binomial log(size): Gaussian(log(2), SD 2);
- spatial range: PC prior, `P(range < 5000 m) = 0.5`;
- spatial marginal SD: PC prior, `P(SD > 2.5) = 0.05`.

Fitted hyperparameters (from the mapping fit,
`results/road_2024/wolf_2024_hyperparameters.csv`):

- zero-inflation probability posterior mean: 0.063;
- negative-binomial size posterior mean: 3.297;
- spatial range posterior mean: 4175 m (95% CrI 2007 to 7618 m);
- spatial SD posterior mean: 0.961 (95% CrI 0.703 to 1.281).

Model comparison:

| Model | WAIC | Delta WAIC |
| --- | ---: | ---: |
| ZINB spatial-month | 933.64 | 0.00 |
| NB spatial-month | 937.31 | 3.67 |
| Poisson spatial-month | 997.30 | 63.66 |

Main diagnostics:

- posterior predictive row and camera total events / zero fraction /
  maximum count: all pass;
- row Pearson dispersion: 0.584; camera Pearson dispersion: 0.238;
- residual Moran's I: -0.036 (expected -0.017), two-sided p = 0.335;
- row PIT KS p-value: 0.03376; camera PIT KS p-value: 0.0001554;
- required diagnostics pass: TRUE;
- temporal residual autocorrelation: within-camera lag-1 r = -0.181,
  p = 0.002179 (n = 284 pairs); date-ordered mean-residual lag-1 ACF: 0.256.
  Residual deployment-order temporal structure remains detectable here,
  unlike the other two surveys. An exploratory correlation between deployment
  timing and northing was near zero, but that result does not rule out timing
  mechanisms or establish absence of bias. The effect of residual temporal
  structure on the mapped surface remains unresolved. The updated code fits
  camera IID and replicated monthly AR1 alternatives, reuses identical spatial
  folds, and exports spatial-surface differences and uncertainty summaries;
- spatial block cross-validation: row 90 percent coverage = 0.96, camera 90
  percent coverage = 0.93;
- prior sensitivity: WAIC, DIC, and posterior hyperparameters are stable
  across the retained prior variants (WAIC 933.42 to 933.89; delta WAIC 0.00
  to 0.48; stability checked, gate not recomputed per variant -- see
  "Sensitivity checks" above);
- mesh sensitivity: WAIC and hyperparameters are stable across the final,
  finer, and coarser mesh variants (WAIC 933.33 to 933.64; delta WAIC 0.00
  to 0.31; same basis as above).

## Final Interpretation

All three wolf models are final for relative encounter-frequency mapping and
pass the required diagnostic gate (camera-level PPC plus residual Moran's I).

The road-camera 2023 model shows no evidence of residual temporal
autocorrelation and is retained as a parsimonious NB model over the
marginally-better-fitting ZINB alternative. The forest-camera 2024 model's
WAIC and posterior hyperparameters hold steady across all six prior variants
and both mesh variants (see "Sensitivity checks" above for what those refits
do and do not recompute); its main limitation is the small number of
independent events (46), which widens posterior uncertainty on month and
spatial effects without indicating a model problem. The road-camera 2024 model passes the required
posterior-predictive and spatial diagnostics and is supported over NB/Poisson
by WAIC; its one open issue is a residual within-camera temporal correlation
whose mechanism and effect on the mapped spatial surface remain unresolved.
Spatial CV and mesh sensitivity alone do not establish absence of temporal distortion.

## Human-Activity Companion Surfaces

Because the camera detections file also records people and vehicles, the
identical pipeline was run for a **human-activity index** on the two road-camera
surveys, as a companion to the wolf maps (a relative human-disturbance surface,
not the study's primary target). "Human activity" is the count of independent
events labelled *Homo sapiens*, `cars`, `bikes`, or `motorcycle`. Only the
detection labels, output prefix, mesh resolution, and a slightly re-centred
spatial-range prior differ from the wolf road runners; the model, diagnostics,
and outputs are otherwise the same.

```sh
Rscript scripts/run_road_2023_human.R
Rscript scripts/run_road_2024_human.R
```

Final outputs are in `results/human_2023/` and `results/human_2024/` (the same
file set as the wolf surveys, with the `human_2023_` / `human_2024_` prefix).

| Survey | Cameras | Events | Effort | Observed rate /100 | Final model | Required diagnostics |
| --- | ---: | ---: | --- | ---: | --- | --- |
| Human-activity 2023 | 60 | 8284 | 5222.2 camera-days | 158.6 | NB spatial-month | Pass |
| Human-activity 2024 | 60 | 6781 | 3574.0 camera-days | 189.7 | NB spatial-month | Fail (residual Moran; see caveat) |

The maps are the same three effort-weighted, annualized surfaces as the wolf
maps -- posterior mean, SD, and CV of expected human-activity events per 100
camera-days (see [Result Maps](#result-maps) for the derivation and display
notes).

**Human-activity 2023**

| Posterior mean | Posterior SD | Posterior CV |
| --- | --- | --- |
| ![Human-activity 2023 posterior mean](results/human_2023/human_2023_final_event_frequency_mean.png) | ![Human-activity 2023 posterior SD](results/human_2023/human_2023_final_event_frequency_sd.png) | ![Human-activity 2023 posterior CV](results/human_2023/human_2023_final_event_frequency_cv.png) |

**Human-activity 2024**

| Posterior mean | Posterior SD | Posterior CV |
| --- | --- | --- |
| ![Human-activity 2024 posterior mean](results/human_2024/human_2024_final_event_frequency_mean.png) | ![Human-activity 2024 posterior SD](results/human_2024/human_2024_final_event_frequency_sd.png) | ![Human-activity 2024 posterior CV](results/human_2024/human_2024_final_event_frequency_cv.png) |

Model comparison strongly prefers the negative-binomial likelihood for both
years (2023: NB best, ZINB delta WAIC 2.78, Poisson 372.6; 2024: NB best, ZINB
2.93, Poisson 806.2 -- Poisson is decisively rejected by the large counts).

Fitted hyperparameters:

- 2023: negative-binomial size 12.64; spatial range 1349 m (95% CrI 123 to
  4620 m); spatial SD 1.16.
- 2024: negative-binomial size 6.51; spatial range 1204 m (95% CrI 328 to
  2643 m); spatial SD 1.06.

Human activity is far less overdispersed than wolves (NB size ~6-13 vs ~1.7-3.3)
and varies on a shorter spatial scale (range ~1.2-1.3 km vs the ~3-4 km of the
road-camera wolf surveys) -- as expected for activity tied to roads and
settlements.

Diagnostics:

- 2023: row Pearson dispersion 0.916; residual Moran's I -0.010 (two-sided
  p = 0.712); spatial block CV row/camera 90% coverage 0.87 / 0.77; required
  diagnostics pass: TRUE.
- 2024: row Pearson dispersion 0.649; residual Moran's I 0.030 (two-sided
  p = 0.017); spatial block CV row/camera 90% coverage 0.89 / 0.80; required
  diagnostics pass: FALSE.

**Mesh and priors.** The human mesh is finer than the road-wolf mesh (inner
edge 300 m vs 400 m) to resolve the shorter human spatial range. The
mesh-sensitivity check confirms the fit is good where it matters: WAIC and the
mapped surface are essentially unchanged across coarser, final, and finer
meshes (17k to 62k vertices). The priors are weakly informative: across the
prior-sensitivity perturbations WAIC, the negative-binomial size, and the
spatial SD barely move. The spatial range is more prior-sensitive (as expected
for a weakly-identified SPDE hyperparameter), but the model choice and the
mapped surface are not.

**Caveat (human-activity 2024).** This survey shows a small but statistically
significant residual spatial autocorrelation (Moran's I 0.030, p = 0.017) that a
finer mesh does not remove -- genuine fine-scale human structure (roads,
settlements) that a smooth stationary field cannot fully absorb. The surface is
retained as a relative human-activity index with this caveat (analogous to the
road-camera 2024 wolf model's residual temporal-autocorrelation caveat); it
should not be read as a fully clean fit.

## Single-Month Spatial Surface: Wolf, March 2024

A companion surface built from the same pipeline with **no month fixed effect**
(`use_month_effect = FALSE`): an intercept + INLA-SPDE spatial field over the
March 2024 road-camera wolf data. It demonstrates the pipeline on a single
calendar month -- where month effects do not apply and there is nothing to
annualize -- and it is the repository's only **spring** wolf surface (the
road-camera multi-month surveys sample summer into autumn only). March 2024 is
drawn from the full multi-year camera dataset (303 distinct wolf events).

```sh
Rscript scripts/run_road_march2024_wolf.R
```

Final outputs are in `results/wolf_march_2024/`.

| Cameras | Wolf events | Effort | Observed rate /100 | Final model | Required diagnostics |
| ---: | ---: | --- | ---: | --- | --- |
| 60 | 303 | 1783.4 camera-days | 16.99 | Poisson spatial (single-period) | Pass |

**Why Poisson.** Unlike the multi-month wolf surveys (where NB/ZINB decisively
beat Poisson), the candidate model comparison here prefers the simpler **Poisson**
likelihood -- the March counts show no overdispersion once the spatial field is
included (fitted Pearson dispersion 0.61):

| Model | WAIC | Delta WAIC |
| --- | ---: | ---: |
| Poisson spatial | 432.38 | 0.00 |
| NB spatial | 438.02 | 5.64 |
| ZINB spatial | 440.20 | 7.82 |

Fitted hyperparameters: spatial range 4039 m (95% CrI 1899 to 7365 m); spatial
SD 0.90. (Poisson has no negative-binomial size parameter.)

Diagnostics: residual Moran's I -0.039 (two-sided p = 0.246); spatial block CV
row/camera 90% coverage 0.94 / 0.92; required diagnostics pass: TRUE. It is a
single-period surface -- expected wolf events per 100 camera-days in March 2024,
with no month adjustment or annualization (scale factor 1.000).

| Posterior mean | Posterior SD | Posterior CV |
| --- | --- | --- |
| ![Wolf March 2024 posterior mean](results/wolf_march_2024/wolf_march_2024_final_event_frequency_mean.png) | ![Wolf March 2024 posterior SD](results/wolf_march_2024/wolf_march_2024_final_event_frequency_sd.png) | ![Wolf March 2024 posterior CV](results/wolf_march_2024/wolf_march_2024_final_event_frequency_cv.png) |

## Single-Month Spatial Surface: Human Activity, March 2024

The human-activity companion to the wolf March 2024 surface above: the same
single-month, spatial-only pipeline (no month fixed effect) applied to the March
2024 human-activity index (Homo sapiens + cars + bikes + motorcycle) -- 2,805
distinct events over the same 60 road cameras.

```sh
Rscript scripts/run_road_march2024_human.R
```

Final outputs are in `results/human_march_2024/`.

| Cameras | Human events | Effort | Observed rate /100 | Final model | Required diagnostics |
| ---: | ---: | --- | ---: | --- | --- |
| 60 | 2805 | 1783.4 camera-days | 157.3 | Poisson spatial (single-period) | Pass |

Like the wolf March surface, the candidate model comparison prefers **Poisson**
(WAIC 20.5 below NB, 21.8 below ZINB): within a single month the spatial field
absorbs the variation and no overdispersion remains (fitted Pearson dispersion
0.58). This is the opposite of the multi-month human surveys, where NB
decisively beats Poisson -- that multi-month overdispersion came from pooling
months with different activity rates.

| Model | WAIC | Delta WAIC |
| --- | ---: | ---: |
| Poisson spatial | 697.66 | 0.00 |
| NB spatial | 718.13 | 20.47 |
| ZINB spatial | 719.45 | 21.78 |

Fitted hyperparameters: spatial range 2615 m (95% CrI 1327 to 4395 m); spatial
SD 1.23. The human range (~2.6 km) is shorter than the wolf March range
(~3.9 km), consistent with human activity being tied to roads and settlements.

Diagnostics: residual Moran's I -0.027 (two-sided p = 0.582); spatial block CV
row/camera 90% coverage 0.91 / 0.87; required diagnostics pass: TRUE. It is a
single-period surface (March 2024); no month adjustment or annualization.

| Posterior mean | Posterior SD | Posterior CV |
| --- | --- | --- |
| ![Human activity March 2024 posterior mean](results/human_march_2024/human_march_2024_final_event_frequency_mean.png) | ![Human activity March 2024 posterior SD](results/human_march_2024/human_march_2024_final_event_frequency_sd.png) | ![Human activity March 2024 posterior CV](results/human_march_2024/human_march_2024_final_event_frequency_cv.png) |

## Outputs Included Here

Across the final-results folders, the curated outputs include:

- final validation report;
- posterior predictive checks;
- hyperparameter summaries;
- month-effect summaries;
- prior sensitivity reports and tables;
- mesh sensitivity reports and tables for every survey;
- model-comparison report and table for every survey (the forest-camera table
  is added by the unified pipeline; see the refactor note near the top);
- spatial block cross-validation summaries;
- temporal residual diagnostics;
- posterior mean encounter-frequency map as PNG and GeoTIFF;
- posterior-SD (absolute uncertainty) map as PNG and GeoTIFF;
- posterior-CV (relative uncertainty) map as PNG and GeoTIFF.

The full generated output folders also contain exploratory plots, full
prediction grids, and additional intermediate diagnostics. Those full scratch
archives are not part of the curated GitHub result set. The
[`results/README.md`](results/README.md) index lists exactly which files are
included in each `results/<survey>/` folder and what each one contains.

## Required R Packages

- `readr`
- `dplyr`
- `tidyr`
- `sf`
- `terra`
- `ggplot2`
- `viridis`
- `scales`
- `sn` (required by INLA joint posterior sampling)
- `INLA`

The scripts check these dependencies before running. Install them once in R:

```r
install.packages(
  c("readr", "dplyr", "tidyr", "sf", "terra", "ggplot2", "viridis", "scales", "sn", "INLA"),
  repos = c(INLA = "https://inla.r-inla-download.org/R/stable", getOption("repos"))
)
```

Package installation remains opt-in; analysis runs do not modify the R environment.

## Reproduction And Validation

Every survey is produced by one command that sources the shared library:

```sh
Rscript scripts/run_road_2023.R           # negative-binomial, road-camera 2023
Rscript scripts/run_road_2024.R           # zero-inflated negative-binomial, road 2024
Rscript scripts/run_forest_2024.R         # negative-binomial, forest-camera 2024
Rscript scripts/run_road_2023_human.R     # companion: 2023 human-activity surface
Rscript scripts/run_road_2024_human.R     # companion: 2024 human-activity surface
Rscript scripts/run_road_march2024_wolf.R # companion: March 2024 wolf (single-month, Poisson)
Rscript scripts/run_road_march2024_human.R # companion: March 2024 human (single-month, Poisson)
```

with `WOLF_RUN_PROFILE` set to `final` for a publication run (see
[Runtime Profiles](#runtime-profiles)). Each writes its full output archive to
`outputs/…`; the curated subset is what lives under `results/`. The raw
camera-trap data are private and will be released with the associated
publication. The public integration example generates synthetic data without
using survey observations; see the [public checks](../README.md#run-the-public-checks).

**Reproducing and checking the committed results.** Because the raw survey data
are private, the committed `results/` are regenerated on the machine that holds
them, then checked:

1. Run each runner with `WOLF_RUN_PROFILE=final`.
2. The corrected method intentionally changes map uncertainty and validation outputs.
   Preserve the old result archive, inspect the new reports, and compare models
   before replacing curated results. After an approved corrected run, record
   numerical tolerances for repeat runs; do not require binary-identical GeoTIFFs
   or claim identity across different INLA builds.
3. Re-run `Rscript scripts/capture_session_info.R` and commit
   `results/session_info.txt` so the reported numbers are tied to a known
   R/INLA build.

## Reproducibility

INLA results can shift across package versions and INLA builds, and this
repository does not pin one. To make a specific run reproducible:

- Run `Rscript scripts/capture_session_info.R` in the same R environment used
  for a final run and commit the resulting `results/session_info.txt`. It
  records `sessionInfo()` plus the exact INLA package version and build.
- For a fuller environment pin, initialize [`renv`](https://rstudio.github.io/renv/)
  in this project (`renv::init()`) and commit the generated `renv.lock`.

The raw camera-trap data are private and are not distributed here; they will be
released with the associated publication (see `data/README.md`).

## References

- Czado, C., Gneiting, T. & Held, L. (2009). Predictive model assessment for
  count data. *Biometrics*, 65(4), 1254-1261.
- Dormann, C. F. et al. (2007). Methods to account for spatial
  autocorrelation in the analysis of species distributional data: a review.
  *Ecography*, 30(5), 609-628.
- Fuglstad, G.-A., Simpson, D., Lindgren, F. & Rue, H. (2019). Constructing
  priors that penalize the complexity of Gaussian random fields. *Journal of
  the American Statistical Association*, 114(525), 445-452.
- Gelman, A., Meng, X.-L. & Stern, H. (1996). Posterior predictive assessment
  of model fitness via realized discrepancies. *Statistica Sinica*, 6(4),
  733-760.
- Hilbe, J. M. (2011). *Negative Binomial Regression* (2nd ed.). Cambridge
  University Press.
- Lindgren, F., Rue, H. & Lindström, J. (2011). An explicit link between
  Gaussian fields and Gaussian Markov random fields: the stochastic partial
  differential equation approach. *Journal of the Royal Statistical Society:
  Series B*, 73(4), 423-498.
- Martin, T. G. et al. (2005). Zero tolerance ecology: improving ecological
  inference by modelling the source of zero observations. *Ecology Letters*,
  8(11), 1235-1246.
- Moran, P. A. P. (1950). Notes on continuous stochastic phenomena.
  *Biometrika*, 37(1/2), 17-23.
- O'Brien, T. G. (2011). Abundance, density and relative abundance: a
  conceptual framework. In *Camera Traps in Animal Ecology* (pp. 71-96).
  Springer.
- Roberts, D. R. et al. (2017). Cross-validation strategies for data with
  temporal, spatial, hierarchical, or phylogenetic structure. *Ecography*,
  40(8), 913-929.
- Rowcliffe, J. M., Field, J., Turvey, S. T. & Carbone, C. (2008). Estimating
  animal density using camera traps without the need for individual
  recognition. *Journal of Applied Ecology*, 45(4), 1228-1236.
- Rue, H., Martino, S. & Chopin, N. (2009). Approximate Bayesian inference
  for latent Gaussian models by using integrated nested Laplace
  approximations. *Journal of the Royal Statistical Society: Series B*,
  71(2), 319-392.
- Simpson, D., Rue, H., Riebler, A., Sørbye, S. H. & Fuglstad, G.-A. (2017).
  Penalising model component complexity: a principled, practical approach to
  constructing priors. *Statistical Science*, 32(1), 1-28.
- Spiegelhalter, D. J., Best, N. G., Carlin, B. P. & van der Linde, A.
  (2002). Bayesian measures of model complexity and fit. *Journal of the
  Royal Statistical Society: Series B*, 64(4), 583-639.
- Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
  widely applicable information criterion in singular learning theory.
  *Journal of Machine Learning Research*, 11, 3571-3594.

