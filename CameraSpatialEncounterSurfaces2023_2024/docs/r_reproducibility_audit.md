# Scoped R reproducibility audit: four statistical corrections

The implementation defects in map uncertainty and spatial-fold handling are
addressed, and explicit calibration/dependence checks are implemented. Survey
adequacy remains unverified until corrected runs on the private inputs complete.

| ID | Severity | File/function | Evidence | Issue and consequence | Fix | Scientific status |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | Major | `make_prediction_outputs`, `joint_surface_summary` | Observed: legacy maps used fixed aggregation/ZI factors | Omitted parameter uncertainty and joint dependence | Joint sampled-rate summaries, credible bounds, and aggregate covariance | Implemented; legacy survey maps require regeneration |
| 2 | Major | `spatial_block_cv`, `run_dependence_comparison` | Observed in legacy reports: 77%/80% camera coverage for 90% intervals | Camera-total predictive undercoverage | Fold coverage and scores, IID/AR1 alternatives, overall review flags | Model adequacy and candidate choice pending private-data reruns |
| 3 | Major | Project README; dependence comparison | Observed: absence of temporal distortion was asserted without a temporal alternative | Unsupported robustness interpretation | Removed claims; added residual and standardized map sensitivity comparisons | Real-data temporal robustness unresolved |
| 4 | Major | `strict_spatial_folds`, CV prior construction | Observed: random fallback and full-data prior centering | Spatial-validation design could change and held-out responses could enter the prior | Fail explicitly; saved geometry-based folds; training-only prior center | Implementation regression checks available |

The approved scope is the user's request to solve these four review findings.
No separate versioned analysis plan was available. Alternative model promotion
requires scientific assessment of the resulting evidence; the code preserves the
configured model and does not manufacture a better-fitting result.

Local verification: all R analysis scripts parse and 31 mathematical/regression
checks pass in R 4.6.0 via webR 0.6.0. This environment has no native R/INLA
installation, so native integration execution is delegated to the committed
GitHub Actions test. See its actual run result before treating integration as
verified. Private-data reproduction has not been performed.

Smallest remaining verification: pass the native synthetic integration job,
rerun the seven final survey commands against the original inputs, inspect the
saved fold/prior evidence and model comparisons, and update numerical results
and figures together. The committed raw input contract and data remain unchanged.
