# GEO-01 Spiral Failure Localization

The open GEO-01 status is not a generic geometry failure. It is specifically
an unrecovered spiral-holdout transfer failure.

Observed comparison on the spiral family:

- HAOS holdout spearman: 0.687323
- best baseline spearman: 0.829858
- HAOS holdout top-k precision: 0.648148
- best baseline top-k precision: 0.814815

Control pattern on the spiral family:

- label permutation collapses the signal sharply
- degree-preserving rewiring degrades, but not enough to close the gap
- parameter-matched null stays above the HAOS holdout score on average
- seed-repeat matches the HAOS score, showing the current bundle is stable but
  not uniquely predictive
- topology-destroyed control does not uniformly degrade, which suggests the
  current feature bundle is partly insensitive to the spiral-specific structure

Interpretation:

1. The baseline is already strong on spiral geometry.
2. The current HAOS feature bundle is too coarse to add unique predictive value.
3. The open criterion is therefore a real transfer-specific gap, not a
   threshold artifact.

Failure class:

`FEATURE_BUNDLE_SPECIFICITY_INSUFFICIENT`

Declared residual criterion:

`holdout transfer does not outperform the best baseline on the spiral family`

Until that changes, GEO-01 stays synthetic OPEN.
