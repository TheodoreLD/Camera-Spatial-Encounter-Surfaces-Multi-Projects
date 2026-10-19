# Bayesian Spatial Encounter-Surface Modelling using Camera-Trap Detections

This repository contains a 2023-2024 camera-trap modelling project focused on
the spatial pattern of wolf and human relative encounter frequency. The analyses model
independent detections recorded by camera traps using Bayesian count models
with INLA-SPDE spatial random fields. Active camera-days are used as exposure,
calendar month is included as a fixed temporal control, and outputs are
relative encounter-frequency surfaces expressed as expected independent
events per 100 camera-days across the sampled survey-year period.

## Project

- [CameraSpatialEncounterSurfaces2023_2024](CameraSpatialEncounterSurfaces2023_2024/):
  seven Bayesian spatial encounter surfaces from camera-trap data. Three are
  the primary wolf surfaces -- a 2023 road-camera negative-binomial
  spatial-month model, a 2024 forest-camera negative-binomial spatial-month
  model, and a 2024 road-camera zero-inflated negative-binomial spatial-month
  model. Four are companions: two multi-month **human-activity** surfaces
  (people and vehicles) as a relative human-disturbance index, and two
  single-month (March 2024) spatial-only surfaces for wolf and human activity.

All seven run the same analyses, diagnostics, and outputs through a single
shared analysis library
([`scripts/wolf_encounter_surface_lib.R`](CameraSpatialEncounterSurfaces2023_2024/scripts/wolf_encounter_surface_lib.R)),
driven by one thin runner per survey. See the
[project README](CameraSpatialEncounterSurfaces2023_2024/README.md) for the full
methodology, diagnostics, and outputs.

## License And Citation

This repository is licensed under [CC BY 4.0](LICENSE). If you use it, please
cite it using the metadata in [CITATION.cff](CITATION.cff).
