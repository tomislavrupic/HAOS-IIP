# Wasserstein Distance as a Transport-Cost Diagnostic

## Status

- Classification: structural equivalence / operational analogy
- Scope: frozen probability measures on discrete graphs, weighted neighborhoods, and synthetic transport comparisons
- Claim ceiling: no theorem about physical transport, no claim that Wasserstein geometry is derived from HAOS-IIP, no assertion that the metric alone explains recovery

## Candidate equivalence statement

In the discrete HAOS-IIP setting, Wasserstein distance is the natural cost language for comparing two neighborhood measures.

The correspondence is:

`Wasserstein distance ≈ minimal transport cost between frozen local distributions`

This note does not claim that HAOS-IIP derives optimal transport. It says that transport cost is the correct local primitive for the Ollivier-Ricci comparison class already recorded in the equivalences layer.

## Why this comparison is justified

The repository already uses:

- graph and normalized Laplacian diagnostics;
- Cheeger conductance and Fiedler-style bottleneck estimates;
- discrete curvature as a local feature;
- Ollivier-Ricci curvature as a neighborhood transport diagnostic;
- rewiring and bridge-removal controls.

Wasserstein distance is the explicit cost functional underneath that Ollivier-style neighborhood comparison. It gives a transparent way to ask how expensive it is to move one lazy neighborhood measure to another across an edge or local orbit.

## What is shared

- local probability measures on a frozen graph or complex;
- cost-sensitive comparison of neighborhoods;
- sensitivity to bottlenecks, bridges, and dense overlap;
- compatibility with curvature, Cheeger, and spectral-gap diagnostics;
- use as a feature, control, or diagnostic ingredient in a frozen synthetic regime.

## What is not shared

- no claim that the transport metric proves smooth geometry;
- no claim that the cost alone determines curvature sign;
- no claim that Wasserstein distance replaces spectral or Hodge diagnostics;
- no claim that the metric settles the GEO-HIDDEN-01 transformation bottleneck.

## Evidence gates

The comparison is only considered operationally valid if:

1. the measures are defined on the frozen discrete graph or complex;
2. the transport cost is computed from the declared metric and supports the declared support sets;
3. the cost moves for the intended rewiring / bridge-removal / topology-destruction controls;
4. the measure does not collapse into a trivial label proxy or degree proxy;
5. the note remains a diagnostic aid, not a verdict engine.

If the feature is not control-separated, the note degrades to a descriptive analogy.

## Failure conditions

This note fails if any of the following happen:

- Wasserstein distance is promoted to a smooth-geometry theorem;
- the transport cost is treated as an ontological statement by default;
- the note is used to infer a physical mechanism;
- the diagnostic is tuned on holdout outcomes;
- the note is used to reopen a frozen physics bridge.

## Audit summary

- Implemented fact: the repo already records spectral, Hodge, curvature, and control-separated bottleneck diagnostics.
- Design choice: add Wasserstein as the explicit transport-cost primitive that supports Ollivier-Ricci.
- Heuristic: minimal transport cost is the right way to quantify neighborhood overlap or mismatch.
- Analogy: transport cost behaves like a local cost-to-shape signal in the existing discrete operator story.
- Unverified hypothesis: a frozen Wasserstein probe may sharpen local transport diagnostics if the feature contract is explicit and holdout remains untouched.

## Related frozen sources

- `docs/notes/equivalences/ollivier_ricci_curvature_local_transport_diagnostic_v1.md`
- `docs/notes/equivalences/discrete_curvature_models_local_to_global_diagnostic_v1.md`
- `docs/notes/equivalences/spectral_graph_clustering_diagnostic_v1.md`
- `docs/notes/equivalences/hodge_laplacian_extensions_higher_order_sector_recovery_v1.md`
- `experiments/geometry_bridge/spectral_diagnostics_summary.md`
