# Spiral Feature Specificity Report

Status: frozen result

Residual under audit

`SPIRAL_HOLDOUT_BASELINE_SUPERIOR`

Failure class

`FEATURE_BUNDLE_SPECIFICITY_INSUFFICIENT`

What was tested

The frozen spiral holdout was evaluated with the existing latent proxy features
already present in the intrinsic-geometry matrix:

- `shortest_path_distance` as the best baseline
- `haos_distance`
- `spectral_distance`
- `invariant_overlap`
- `persistence_score`
- `causal_depth`
- `temporal_order_stability`

The decisive comparison was:

- best baseline
- versus best baseline + one precommitted proxy invariant

The training freeze was respected:

- train on ring + grid
- hold out spiral
- no score tuning after inspection

Observed holdout ranking

- best baseline (`shortest_path_distance`): mean Spearman `0.829858`
- best augmented proxy (`shortest_path_distance + persistence_score`): mean Spearman `0.829858`

No candidate produced a stable positive increment over the best baseline.

Per-candidate mean holdout Spearman

- `shortest_path_distance + haos_distance`: `0.707593`
- `shortest_path_distance + spectral_distance`: `0.706821`
- `shortest_path_distance + invariant_overlap`: `0.676319`
- `shortest_path_distance + persistence_score`: `0.829858`
- `shortest_path_distance + causal_depth`: `0.143436`
- `shortest_path_distance + temporal_order_stability`: `0.447362`

Interpretation

1. The spiral baseline remains strong on its own.
2. The currently frozen HAOS proxy bundle does not add stable holdout value.
3. The best augmented proxy merely ties the baseline on Spearman and does not
   establish a distinct gain.
4. There is no evidence here that a latent invariant is already recoverable
   from the current bundle in a way that improves spiral transfer.

Frozen outcome

`SPIRAL_SPECIFIC_STRUCTURE_NOT_RECOVERED`

Boundary

This is a synthetic calibration result only. It does not authorize any
external geometry, physics, or mechanism claim.
