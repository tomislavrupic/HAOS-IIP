# PB-01 Network Recovery Bridge

Status: `PREDICTION_NOT_DISTINCT_FROM_BASELINES`

This is a synthetic network-recovery calibration benchmark. It is not empirical physical validation.

## Verdict Labels

- `CONTROL_INVALID`
- `HOLDOUT_TRANSFER_PASS`
- `MIXED_OPEN`
- `PHYSICAL_MECHANISM_NOT_ESTABLISHED`
- `PREDICTION_NOT_DISTINCT_FROM_BASELINES`

## Holdout Comparison

- HAOS holdout Spearman: `0.912028`
- best baseline: `supervised_graph_spectral_model`
- best baseline holdout Spearman: `0.800812`
- margin: `0.111217`
- holdout degradation: `-0.021880`

## Boundary

- No empirical bridge claim.
- No physical mechanism claim.
- HAOS-specific scores are secondary diagnostics.
- A `PREDICTION_NOT_DISTINCT_FROM_BASELINES` verdict is a valid benchmark outcome.
