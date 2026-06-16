# Synthetic Intrinsic Geometry Recovery Report

- bridge id: GEO-01
- version: synthetic-intrinsic-geometry-recovery-v0.1
- verdict: CONTROL_INVALID, GEOMETRY_NOT_DEMONSTRATED, MIXED_OPEN
- calibration pass: False
- holdout pass: False

## Key results
- holdout HAOS spearman: 0.687323
- best baseline spearman: 0.829858
- holdout top-k precision: 0.648148
- best baseline top-k precision: 0.814815

## Interpretation
This is a synthetic blind geometry calibration. It tests whether operator-only transport features
recover latent intrinsic distances better than conventional baselines under frozen holdout.
It does not authorize Bell, quantum, or physical mechanism claims.
