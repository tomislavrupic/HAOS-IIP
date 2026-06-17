# Maximum Mean Discrepancy as a Kernel Diagnostic

## Status

- Classification: structural equivalence / operational analogy
- Scope: frozen distributions, kernel embeddings, and synthetic comparison regimes
- Claim ceiling: no claim that MMD is a transport metric, no claim that kernel discrepancy proves geometry, no claim that the diagnostic closes the GEO-HIDDEN-01 gap

## Candidate equivalence statement

In the discrete HAOS-IIP setting, Maximum Mean Discrepancy is a fast kernel-side diagnostic for comparing two frozen distributions or feature clouds.

The correspondence is:

`MMD ≈ RKHS mean-embedding discrepancy`

This note does not claim that MMD recovers Wasserstein geometry. It says that kernel mean discrepancy is the appropriate lightweight comparison class when the question is whether two frozen embeddings, histograms, or curvature summaries look distributionally different.

## Why this comparison is justified

The repository already uses:

- spectral embeddings and Laplacian eigenmodes;
- curvature summaries and transport diagnostics;
- holdout-separated synthetic comparison regimes;
- control families that include shuffling, rewiring, and topology destruction.

MMD fits that stack because it compares mean embeddings in a reproducing kernel Hilbert space rather than solving transport. That makes it useful as a quick distribution test alongside the heavier transport family.

## What is shared

- distributions on frozen synthetic features;
- sensitivity to class separation in embedding space;
- compatibility with spectral or curvature-derived features;
- use as a control, sanity check, or fast screening probe;
- ability to compare pre/post or class-vs-class summaries cheaply.

## What is not shared

- no claim that MMD measures geometric mass transport;
- no claim that the RKHS embedding is a physical manifold;
- no claim that MMD replaces Wasserstein, Sinkhorn, or sliced transport;
- no claim that the kernel discrepancy settles the GEO-HIDDEN-01 transformation bottleneck.

## Evidence gates

The comparison is only considered operationally valid if:

1. the kernel is declared before evaluation;
2. the feature map is defined on the frozen synthetic inputs;
3. the comparison is run against the same declared controls as the transport probes;
4. the statistic does not collapse into a trivial label proxy or a training artifact;
5. the note remains a diagnostic aid, not a verdict engine.

If the statistic is not control-separated, the note degrades to a descriptive analogy.

## Failure conditions

This note fails if any of the following happen:

- MMD is promoted to a transport metric;
- the RKHS discrepancy is treated as an ontological statement by default;
- the note is used to infer a physical mechanism;
- the diagnostic is tuned on holdout outcomes;
- the note is used to reopen a frozen physics bridge.

## Audit summary

- Implemented fact: kernel mean discrepancy is fast and differentiable compared with full transport.
- Design choice: add MMD as the lightweight distribution-test layer beside the transport stack.
- Heuristic: kernel embeddings are a good fit when the question is simply whether two feature clouds separate.
- Analogy: MMD behaves like a mean-embedding separation score in the existing discrete operator story.
- Unverified hypothesis: a frozen MMD probe may provide a cheap screening test before heavier transport diagnostics.

## Related frozen sources

- `docs/notes/equivalences/mmd_kernel_selection_strategies_v1.md`
- `docs/notes/equivalences/spectral_address_laplacian_eigenmodes_v1.md`
- `docs/notes/equivalences/spectral_graph_clustering_diagnostic_v1.md`
- `docs/notes/equivalences/wasserstein_distance_transport_cost_diagnostic_v1.md`
- `docs/notes/equivalences/sinkhorn_approximation_transport_diagnostic_v1.md`
- `docs/notes/equivalences/sliced_wasserstein_transport_diagnostic_v1.md`
- `experiments/geometry_bridge/spectral_diagnostics_summary.md`
