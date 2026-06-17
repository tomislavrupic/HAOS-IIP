# Discrete Ricci Flow on Graphs as a Diagnostic Probe

## Status

- Classification: structural equivalence / operational analogy
- Scope: frozen discrete graphs and weighted cell-complex approximations
- Claim ceiling: no smooth-manifold theorem, no geometric-optimization claim, no physical-curvature claim

## Candidate equivalence statement

In the discrete HAOS-IIP setting, Ricci-flow-like updates can be read as a local-to-global diagnostic for graph geometry.

The correspondence is:

`discrete Ricci flow ≈ curvature-weighted diagnostic evolution on frozen graphs`

This note does not claim that the repository has implemented a true discrete Ricci flow theorem. It says that curvature-aware reweighting is a useful companion language for watching bridges shrink, communities sharpen, and bottlenecks reorganize in a controlled synthetic benchmark.

## Why this comparison is justified

The repository already uses:

- graph and normalized Laplacian diagnostics;
- Cheeger conductance and Fiedler-style bottleneck estimates;
- local curvature as a discrete feature;
- rewiring and bridge-removal controls;
- frozen synthetic holdout separation.

Those ingredients make a flow-like curvature update a natural diagnostic extension, especially if it is used only to compare pre/post structure on synthetic graphs.

## What is shared

- edge-weight evolution under a local diagnostic rule;
- sensitivity to bottlenecks and bridge edges;
- compatibility with curvature features such as Forman or Ollivier-style estimates;
- ability to emphasize communities while weakening sparse cuts;
- use as a synthetic probe, not a claim about external geometry.

## What is not shared

- no claim that the update solves Ricci flow in the smooth setting;
- no claim that the flow converges to a canonical geometry;
- no claim that weight evolution proves a physical metric flow;
- no claim that a short flow pass closes the GEO-HIDDEN-01 transformation gap.

## Evidence gates

The comparison is only considered operationally valid if:

1. the update rule is defined on the frozen discrete graph;
2. the update is local and transparent enough to audit;
3. the update moves the intended curvature or bottleneck statistics;
4. control graphs show the expected change or failure to change;
5. the flow is reported as a diagnostic layer, not a verdict engine.

If the update cannot be control-separated, the note degrades to a descriptive analogy.

## Failure conditions

This note fails if any of the following happen:

- the flow is promoted to a smooth Ricci flow theorem;
- the update is treated as ontic curvature evolution by default;
- the diagnostic is tuned on holdout outcomes;
- the flow is used to infer a physical mechanism;
- the note is used to reopen a frozen physics bridge.

## Audit summary

- Implemented fact: the repo already uses curvature, Cheeger, and spectral diagnostics to interpret bottlenecks.
- Design choice: describe a Ricci-flow-like pass as a diagnostic evolution rule rather than a claim about geometry.
- Heuristic: a short 5–10 step curvature-weighted update is enough to expose whether bottlenecks tighten or dissolve in the synthetic regime.
- Analogy: the update behaves like a discrete flow on graph geometry, not a physical metric evolution.
- Unverified hypothesis: a frozen Ricci-flow diagnostic may sharpen the transformation bottleneck if a contract is defined and holdout remains untouched.

## Related frozen sources

- `docs/notes/equivalences/discrete_curvature_models_local_to_global_diagnostic_v1.md`
- `docs/notes/equivalences/spectral_graph_clustering_diagnostic_v1.md`
- `docs/notes/equivalences/hodge_laplacian_extensions_higher_order_sector_recovery_v1.md`
- `experiments/geometry_bridge/spectral_diagnostics_summary.md`
