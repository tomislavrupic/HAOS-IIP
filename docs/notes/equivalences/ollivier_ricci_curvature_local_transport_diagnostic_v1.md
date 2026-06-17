# Ollivier-Ricci Curvature as a Local Transport Diagnostic

## Status

- Classification: structural equivalence / operational analogy
- Scope: frozen discrete graphs and weighted neighborhood transport on synthetic operator regimes
- Claim ceiling: no smooth-manifold theorem, no physical-curvature claim, no assertion that transport curvature alone explains the geometry bridge

## Candidate equivalence statement

In the discrete HAOS-IIP setting, Ollivier-Ricci curvature is a strong local transport language for explaining why some edges or neighborhoods behave like bridges and others behave like cohesive interior structure.

The correspondence is:

`Ollivier-Ricci curvature ≈ neighborhood transport overlap diagnostic on frozen graphs`

This note does not claim that HAOS-IIP derives Ricci curvature or external geometry. It says that transport-based curvature is a useful companion language for diagnosing whether a frozen synthetic graph has bridge-like edges, community-like interiors, or fragile bottlenecks. The transport cost underneath that comparison is Wasserstein distance, which is the explicit primitive recorded in the companion note.

## Why this comparison is justified

The repository already uses:

- graph and normalized Laplacian diagnostics;
- Cheeger conductance and Fiedler-style bottleneck estimates;
- discrete curvature as a local feature;
- rewiring and bridge-removal controls;
- frozen synthetic holdout separation.

Ollivier-Ricci fits that stack because it asks a local question about how much neighborhood mass can be transported cheaply across an edge. That makes it a natural diagnostic companion to the repo’s spectral and Hodge-style probes.

## What is shared

- local graph or cell-complex structure;
- sensitivity to bottlenecks, bridges, and dense communities;
- compatibility with Laplacian, Cheeger, and spectral-gap diagnostics;
- use as a feature, control, or regularizer on frozen synthetic complexes;
- ability to distinguish edges that act like transport bridges from edges that sit inside coherent neighborhoods.

## What is not shared

- no claim that the curvature definition proves smooth Ricci geometry;
- no claim that positive or negative curvature is an ontological property by itself;
- no claim that transport curvature replaces spectral or Hodge diagnostics;
- no claim that one local transport model settles the GEO-HIDDEN-01 transformation recovery bottleneck.

## Evidence gates

The comparison is only considered operationally valid if:

1. the transport measures are defined on the frozen discrete graph or complex;
2. the transport feature is stable enough to be compared against the declared controls;
3. the feature moves for the intended rewiring / bridge-removal / topology-destruction controls;
4. the feature does not reduce to a trivial degree proxy or a label proxy;
5. the note remains a diagnostic aid, not a verdict engine.

If the feature is not control-separated, the note degrades to a descriptive analogy.

## Failure conditions

This note fails if any of the following happen:

- the curvature model is promoted to a smooth-manifold theorem;
- the curvature value is treated as physical curvature by default;
- the note is used to infer a physical mechanism;
- the diagnostic is tuned on holdout outcomes;
- the note is used to reopen a frozen physics bridge.

## Audit summary

- Implemented fact: the repo already uses spectral gap, Cheeger, Hodge, and local curvature diagnostics to interpret bottlenecks.
- Design choice: add Ollivier-Ricci as a transport-based local diagnostic rather than a claim about ontology.
- Heuristic: a neighborhood transport comparison is the right way to ask whether an edge is bridge-like or community-like.
- Analogy: curvature behaves like a local transport regularizer for the existing discrete operator story.
- Unverified hypothesis: Ollivier-Ricci features may sharpen the GEO-HIDDEN-01 transformation bottleneck diagnosis if the feature contract is frozen and holdout remains untouched.

## Related frozen sources

- `docs/notes/equivalences/wasserstein_distance_transport_cost_diagnostic_v1.md`
- `docs/notes/equivalences/discrete_curvature_models_local_to_global_diagnostic_v1.md`
- `docs/notes/equivalences/spectral_graph_clustering_diagnostic_v1.md`
- `docs/notes/equivalences/spectral_address_laplacian_eigenmodes_v1.md`
- `docs/notes/equivalences/hodge_laplacian_extensions_higher_order_sector_recovery_v1.md`
- `docs/notes/equivalences/current_closure_discrete_conservation_balance_v1.md`
- `experiments/geometry_bridge/spectral_diagnostics_summary.md`
