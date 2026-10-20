# Bayesian Camera-Trap Encounter Surfaces

[![Statistical validation](https://github.com/TheodoreLD/Camera-Spatial-Encounter-Surfaces-Multi-Projects/actions/workflows/statistical-validation.yml/badge.svg)](https://github.com/TheodoreLD/Camera-Spatial-Encounter-Surfaces-Multi-Projects/actions/workflows/statistical-validation.yml)

A reusable R workflow for estimating spatial patterns of wolf and human activity
from camera-trap detections. Seven survey analyses share one INLA-SPDE modelling
library, with exposure adjustment, count-likelihood comparison, spatial
cross-validation, and uncertainty maps.

The target is **expected independent events per 100 camera-days** over the
sampled period. These surfaces describe relative encounter frequency; they do
not estimate abundance, density, or occupancy.

| Area | What the code demonstrates |
| --- | --- |
| Bayesian statistics | Poisson, negative-binomial and zero-inflated models; explicit priors; joint posterior simulation |
| Spatial analysis | SPDE meshes, coordinate transformations, raster prediction, spatial camera holdouts |
| Validation | Posterior predictive checks, calibration and interval scores, residual diagnostics, camera IID/AR1 sensitivity |
| Software engineering | Shared functions and thin survey runners, strict input checks, seeded simulations, regression tests and native INLA CI |

## Start here

- [Methods and survey results](CameraSpatialEncounterSurfaces2023_2024/README.md)
- [Input schemas and data access](CameraSpatialEncounterSurfaces2023_2024/data/README.md)
- [Curated maps and diagnostic files](CameraSpatialEncounterSurfaces2023_2024/results/README.md)
- [Statistical corrections and remaining survey validation](CameraSpatialEncounterSurfaces2023_2024/docs/statistical_corrections.md)
- [Shared analysis library](CameraSpatialEncounterSurfaces2023_2024/scripts/wolf_encounter_surface_lib.R)

## Run the public checks

From the repository root, run the base-R regression suite:

```sh
Rscript CameraSpatialEncounterSurfaces2023_2024/tests/test_statistical_validation.R
```

After installing the [analysis dependencies](CameraSpatialEncounterSurfaces2023_2024/README.md#required-r-packages),
run the native INLA example and integration checks:

```sh
CAMERA_SYNTHETIC_OUTPUT_DIR=synthetic-integration-outputs Rscript CameraSpatialEncounterSurfaces2023_2024/tests/test_inla_integration.R
```

This generates seeded synthetic camera-month data, fits models, builds maps,
runs spatial holdouts and IID/AR1 comparisons, and writes diagnostic reports.
It requires no private survey data. GitHub Actions runs both suites and retains
the synthetic outputs and session information as a downloadable artifact for
30 days. The synthetic example tests implementation, not ecological conclusions.

## Survey validation status

The original survey inputs are private. The committed survey maps and diagnostic
numbers are a **legacy snapshot** from before the joint-uncertainty and
spatial-validation corrections. They require regeneration with the corrected
code. Human-activity predictive undercoverage and residual temporal dependence
remain open scientific checks; the [rerun guide](CameraSpatialEncounterSurfaces2023_2024/docs/statistical_corrections.md)
states what evidence is still needed.

## License and citation

[CC BY 4.0](LICENSE). Citation metadata is in [CITATION.cff](CITATION.cff).
