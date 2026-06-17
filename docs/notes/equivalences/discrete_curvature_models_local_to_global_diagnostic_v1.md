# Discrete Curvature Models and Local-to-Global Diagnostics

## Status

- Classification: structural equivalence / operational analogy
- Scope: frozen discrete graphs, cell complexes, and simplicial / clique-complex regimes
- Claim ceiling: no smooth-manifold claim, no ontological curvature claim, no assertion that curvature alone explains the geometry bridge

## Candidate equivalence statement

In the discrete HAOS-IIP setting, local curvature models are mathematically adjacent to the project’s existing operator diagnostics.

The correspondence is:

`discrete curvature features ≈ local-to-global diagnostics for graph and complex structure`

This note does not claim that curvature is derived from HAOS-IIP. It says that curvature is a useful companion language for understanding why some frozen graphs are easy or hard to recover with spectral and Hodge-style probes.

## Why this comparison is justified

The repository already uses:

- graph Laplacians and normalized Laplacians;
- Cheeger conductance and spectral gaps;
- Hodge-style exact / coexact / harmonic bookkeeping;
- higher-order cochain language;
- control families that compare label permutation, rewiring, and topology destruction.

Discrete curvature lives in the same neighborhood: it is a local descriptor that can explain global separation, diffusion, and recovery behavior without assuming a smooth manifold. In the current note set, Ollivier-Ricci curvature is the stronger transport-based comparison class, while Forman curvature stays the lightweight combinatorial baseline.

## What is shared

- local graph or cell-complex structure;
- sensitivity to bottlenecks, bridges, and hubs;
- compatibility with Laplacian and Cheeger-style diagnostics;
- use as a feature, control, or regularizer on frozen synthetic complexes;
- ability to distinguish regions that are structurally easy versus structurally hard to separate.

## What is not shared

- no claim that curvature proves continuum geometry;
- no claim that a curvature sign determines ontology;
- no claim that curvature replaces spectral or Hodge diagnostics;
- no claim that one local curvature model settles all hidden-geometry recovery questions.

## Evidence gates

The comparison is only considered operationally valid if:

1. the curvature quantity is defined on the frozen discrete graph or complex;
2. the curvature feature is stable enough to be compared against the declared controls;
3. the feature moves for the intended rewiring / bridge-removal controls;
4. the feature does not collapse into a trivial degree proxy;
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

- Implemented fact: the repo already uses spectral gap, Cheeger, Hodge, and control-separated recovery diagnostics.
- Design choice: add curvature as a local diagnostic layer rather than a claim about ontology.
- Heuristic: Forman curvature is the lightest useful probe; Ollivier-Ricci curvature is the stronger geometric comparison class.
- Analogy: curvature behaves like a local-to-global regularizer for the existing discrete operator story.
- Unverified hypothesis: curvature features may help diagnose the GEO-HIDDEN-01 transformation bottleneck if a frozen feature contract is defined.

## Related frozen sources

- `docs/notes/equivalences/wasserstein_distance_transport_cost_diagnostic_v1.md`
- `docs/notes/equivalences/ollivier_ricci_curvature_local_transport_diagnostic_v1.md`
- `docs/notes/equivalences/spectral_graph_clustering_diagnostic_v1.md`
- `docs/notes/equivalences/spectral_address_laplacian_eigenmodes_v1.md`
- `docs/notes/equivalences/hodge_laplacian_extensions_higher_order_sector_recovery_v1.md`
- `docs/notes/equivalences/hodge_decomposition_discrete_sector_split_v1.md`
- `docs/notes/equivalences/current_closure_discrete_conservation_balance_v1.md`
- `experiments/geometry_bridge/spectral_diagnostics_summary.md`
