# Sinkhorn Approximation as a Practical Transport Diagnostic

## Status

- Classification: structural equivalence / operational analogy
- Scope: frozen discrete transport comparisons on synthetic graphs and neighborhood measures
- Claim ceiling: no claim that entropic transport equals exact Wasserstein distance, no claim that the approximation improves ontology, no claim that the approximation closes the geometry bridge

## Candidate equivalence statement

In the discrete HAOS-IIP setting, Sinkhorn approximation is the practical version of transport cost computation when exact optimal transport is too expensive.

The correspondence is:

`Sinkhorn approximation ≈ entropically regularized transport diagnostic`

This note does not claim that Sinkhorn recovers exact Wasserstein geometry. It says that entropic regularization is a reasonable computational compromise for frozen neighborhood transport probes, especially when the exact edge-wise cost is too heavy to use directly.

## Why this comparison is justified

The repository already uses:

- Wasserstein distance as the transport-cost primitive;
- Ollivier-Ricci curvature as the neighborhood transport diagnostic;
- graph and normalized Laplacian diagnostics;
- Cheeger conductance and Fiedler-style bottleneck estimates;
- rewiring and bridge-removal controls.

Sinkhorn appears in that stack as the standard approximation layer: it trades exactness for speed while keeping the transport framing explicit. That makes it a practical companion to the exact transport note rather than a new claim.

## What is shared

- probability measures on a frozen graph or complex;
- transport-plan language with a declared cost matrix;
- sensitivity to bottlenecks, bridges, and dense overlap;
- compatibility with Ollivier-style curvature diagnostics;
- use as a computational probe in a frozen synthetic regime.

## What is not shared

- no claim that Sinkhorn gives exact optimal transport;
- no claim that the entropic regularizer is physically meaningful;
- no claim that Sinkhorn replaces spectral or Hodge diagnostics;
- no claim that the approximation settles the GEO-HIDDEN-01 transformation bottleneck.

## Evidence gates

The comparison is only considered operationally valid if:

1. the cost matrix and support sets are defined on the frozen discrete graph or complex;
2. the entropic parameter is declared before evaluation;
3. the approximation is compared against the same declared controls as the exact transport probe;
4. the approximation does not collapse into a label proxy or a degree proxy;
5. the note remains a diagnostic aid, not a verdict engine.

If the approximation is not control-separated, the note degrades to a descriptive analogy.

## Failure conditions

This note fails if any of the following happen:

- Sinkhorn is treated as exact Wasserstein distance;
- the entropic regularizer is promoted to an ontological statement;
- the note is used to infer a physical mechanism;
- the approximation is tuned on holdout outcomes;
- the note is used to reopen a frozen physics bridge.

## Audit summary

- Implemented fact: exact transport is expensive in edge-wise Ollivier-style probes.
- Design choice: use Sinkhorn as the practical approximation layer for transport diagnostics.
- Heuristic: entropic transport is the right speed/precision compromise when the support sets are large.
- Analogy: the approximation behaves like a smoother transport-cost estimator in the existing discrete operator story.
- Unverified hypothesis: a frozen Sinkhorn probe may make Ollivier-style diagnostics more scalable without changing the open verdict ceiling.

## Related frozen sources

- `docs/notes/equivalences/wasserstein_distance_transport_cost_diagnostic_v1.md`
- `docs/notes/equivalences/ollivier_ricci_curvature_local_transport_diagnostic_v1.md`
- `docs/notes/equivalences/sliced_wasserstein_transport_diagnostic_v1.md`
- `docs/notes/equivalences/discrete_curvature_models_local_to_global_diagnostic_v1.md`
- `experiments/geometry_bridge/spectral_diagnostics_summary.md`
