# EL-R4-OC-01 Result Classification

Classification: `NEGATIVE_RESULT`

## What Passed

- Validation aggregate internal advantage: `0.074074`.
- Holdout aggregate internal advantage: `0.067901`.
- Positive holdout strata fraction: `0.833333`.
- Transition-shuffled training degraded the internal advantage by `0.179012`.
- Context ablation degraded the internal advantage by `0.339506`.
- Source, split, control, and deterministic-replay checks passed.

## What Failed

The predeclared run-level uncertainty gate failed. The holdout 95% bootstrap
interval for internal advantage is:

```text
[-0.037037, 0.046296]
```

Because the interval includes zero, Rung 4 operational closure is not
supported by the compressed Phase XVI event-chain mechanism.

## Frozen Interpretation

The current mechanism contains a reproducible aggregate predictive signal, but
the advantage is not stable enough across stored run keys. No threshold,
predictor key, split, or control will be changed in this version.

The next legitimate mechanism must use intervention-native state trajectories
with separately measured internal state and external forcing.
