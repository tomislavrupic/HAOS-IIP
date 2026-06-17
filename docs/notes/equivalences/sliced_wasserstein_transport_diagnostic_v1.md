# Sliced Wasserstein Distance as a Scalable Transport Diagnostic

## Status

- Classification: structural equivalence / operational analogy
- Scope: frozen discrete transport comparisons, projected neighborhood measures, and synthetic approximation regimes
- Claim ceiling: no claim that sliced transport equals exact Wasserstein distance, no claim that projection averaging proves geometry, no claim that the approximation closes the geometry bridge

## Candidate equivalence statement

In the discrete HAOS-IIP setting, sliced Wasserstein is the fast approximation layer for comparing neighborhood measures when full transport is too expensive.

The correspondence is:

`sliced Wasserstein ≈ projected transport diagnostic`

This note does not claim that sliced transport recovers exact Wasserstein geometry. It says that random projection averaging is a practical way to get a cheaper, lower-fidelity transport signal when the support geometry is large or anisotropic.

## Why this comparison is justified

The repository already uses:

- Wasserstein distance as the transport-cost primitive;
- Sinkhorn as the practical approximation layer;
- Ollivier-Ricci curvature as the neighborhood transport diagnostic;
- graph and normalized Laplacian diagnostics;
- Cheeger conductance and Fiedler-style bottleneck estimates.

Sliced Wasserstein appears one rung further down the approximation stack: instead of solving the full transport problem, it projects the measures to 1D and averages the exact one-dimensional transport costs. That makes it the fastest practical approximation in the transport family recorded here.

## What is shared

- probability measures on a frozen graph or complex;
- transport-cost language with a declared projection rule;
- sensitivity to bottlenecks, bridges, and directional mismatch;
- compatibility with Ollivier-style curvature diagnostics;
- use as a fast probe in a frozen synthetic regime.

## What is not shared

- no claim that sliced transport gives exact optimal transport;
- no claim that projection averaging preserves every geometric feature;
- no claim that sliced transport replaces spectral or Hodge diagnostics;
- no claim that the approximation settles the GEO-HIDDEN-01 transformation bottleneck.

## Evidence gates

The comparison is only considered operationally valid if:

1. the measures are defined on the frozen discrete graph or complex;
2. the projection rule is declared before evaluation;
3. the approximation is compared against the same declared controls as the exact and Sinkhorn probes;
4. the approximation does not collapse into a label proxy or a degree proxy;
5. the note remains a diagnostic aid, not a verdict engine.

If the approximation is not control-separated, the note degrades to a descriptive analogy.

## Failure conditions

This note fails if any of the following happen:

- sliced transport is treated as exact Wasserstein distance;
- projection averaging is promoted to an ontological statement;
- the note is used to infer a physical mechanism;
- the approximation is tuned on holdout outcomes;
- the note is used to reopen a frozen physics bridge.

## Audit summary

- Implemented fact: exact transport is the most expensive local option in the Ollivier-style family.
- Design choice: use sliced Wasserstein as the fastest projection-based approximation layer.
- Heuristic: averaging many 1D transport costs is a good tradeoff when support geometry is large or noisy.
- Analogy: sliced transport behaves like a projection-based transport-cost estimator in the existing discrete operator story.
- Unverified hypothesis: a frozen sliced probe may make the transport diagnostics more scalable without changing the open verdict ceiling.

## Related frozen sources

- `docs/notes/equivalences/wasserstein_distance_transport_cost_diagnostic_v1.md`
- `docs/notes/equivalences/sinkhorn_approximation_transport_diagnostic_v1.md`
- `docs/notes/equivalences/ollivier_ricci_curvature_local_transport_diagnostic_v1.md`
- `docs/notes/equivalences/discrete_curvature_models_local_to_global_diagnostic_v1.md`
- `experiments/geometry_bridge/spectral_diagnostics_summary.md`
