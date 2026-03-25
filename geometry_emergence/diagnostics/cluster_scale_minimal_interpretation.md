# Cluster-Scale Sweep: Minimal Interpretation Layer

This note is bounded to the deterministic outputs in:

- `geometry_emergence/diagnostics/cluster_scale_sweep_summary.json`
- `geometry_emergence/diagnostics/cluster_scale_sweep_full.json`

## Supported statements

1. Transition onset depends on cluster scale.
   The mean onset rises from `0.1118` at scale `0.50` to `0.1562` at scale `1.25`, with a small relaxation to `0.1490` at scale `1.50`.

2. Detection is stable across the full scale panel.
   The detected fraction is `1.0` for every tested scale.

3. Seed dependence becomes structured rather than arbitrary.
   The sweep resolves a low-variance regime at scale `0.50`, a broader mid-scale regime at `0.75` to `1.00`, and a tighter late regime at `1.25` to `1.50`.

4. Cluster scale acts as a meaningful control parameter for the coupled transition signal.
   The onset-vs-scale monotonicity score is `0.9`, and onset variance is reduced relative to the previous clustered canonical panel.

## Unsupported statements

This experiment does not support:

- spacetime claims
- curvature claims
- continuum claims
- physical geometry claims

It supports only the narrower diagnostic statement that coupled transport-geometry signatures in interaction graphs vary systematically with cluster structural scale.

