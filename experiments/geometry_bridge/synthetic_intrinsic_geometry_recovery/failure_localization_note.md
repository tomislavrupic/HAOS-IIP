# GEO-01 Failure Localization Note

The intrinsic-geometry benchmark remains open on the spiral holdout because
the best conventional baseline still outperforms the frozen HAOS-derived
geometry score on that family.

Observed failure pattern:

- holdout transfer on spiral does not beat the best baseline
- the HAOS score remains below the shortest-path baseline
- topology destruction does not uniformly collapse all control signals
- degree-preserving rewiring and parameter-matched nulls preserve enough
  structure that the HAOS feature bundle does not separate them cleanly from
  the target

Most likely contributors:

1. Baseline already captures much of the spiral geometry
2. The current HAOS feature bundle is still too coarse to uniquely encode the
   spiral manifold
3. The holdout transfer is therefore a genuine specificity failure, not a
   threshold artifact

What this does not show:

- that geometry is absent in general
- that the residual failure is Bell-related
- that the synthetic chain is invalid

Declared residual criterion:

`holdout transfer does not outperform the best baseline on the spiral family`

Until that criterion is closed, GEO-01 stays open.

