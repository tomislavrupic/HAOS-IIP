# MMD Kernel Selection Strategies

## Status

- Classification: operational guidance / structural analogy
- Scope: frozen kernel-based distribution tests on synthetic feature clouds
- Claim ceiling: no claim that any kernel choice closes the GEO-HIDDEN-01 gap, no claim that kernel selection proves geometry, no claim that a chosen kernel is universally optimal

## Candidate equivalence statement

In the discrete HAOS-IIP setting, kernel choice is the main practical lever for Maximum Mean Discrepancy.

The correspondence is:

`kernel selection ≈ sensitivity control for MMD on frozen features`

This note does not claim that one kernel is objectively correct. It says that the kernel defines which discrepancies the MMD probe will notice most easily, so kernel choice should be treated as a frozen diagnostic parameter rather than an afterthought.

## Why this comparison is justified

The repository already uses:

- spectral embeddings and Laplacian eigenmodes;
- curvature histograms and transport diagnostics;
- holdout-separated synthetic comparisons;
- control families that include shuffling, rewiring, and topology destruction.

MMD is only as good as the kernel that underlies it. In practice, the kernel determines whether the probe is mostly sensitive to shifts in mean location, scale, tails, or multi-scale differences. That is exactly the sort of choice that should be frozen before comparison.

## Common kernel options

| Kernel | Strengths | Weaknesses | Best use here |
| --- | --- | --- | --- |
| Gaussian / RBF | Flexible, universal, smooth | Bandwidth-sensitive | Default for spectral embeddings |
| Laplace / exponential | More robust to outliers | Less smooth | Curvature histograms or noisy features |
| Polynomial | Captures higher-order moments | Not characteristic at finite degree | Small low-dimensional embeddings |
| Rational quadratic | Multi-scale behavior | More expensive | Mixed-scale features |
| Linear | Very fast, simple | Only first-moment differences | Quick sanity checks |
| Energy-distance kernel | No bandwidth tuning | Less flexible than RBF | Baseline comparison and cross-check |

## Selection strategies

1. **Median heuristic**
   - set the Gaussian bandwidth from the median pairwise distance
   - good default when no frozen bandwidth is already justified

2. **Cross-validation**
   - choose a kernel parameter on a development split only
   - useful when test power matters and the split discipline is frozen

3. **Multiple-kernel ensemble**
   - combine several kernels with fixed or learned weights
   - useful when no single scale is obviously correct

4. **Domain-specific choice**
   - Gaussian for spectral embeddings
   - Laplace or rational quadratic for curvature histograms
   - linear for very fast first-pass screening

## What is shared

- kernel-based comparison of frozen distributions;
- sensitivity to class separation in embedding space;
- compatibility with spectral or curvature-derived features;
- use as a screening probe before heavier transport diagnostics;
- ability to compare pre/post or class-vs-class summaries cheaply.

## What is not shared

- no claim that kernel tuning recovers Wasserstein geometry;
- no claim that bandwidth choice is a physical parameter;
- no claim that a selected kernel proves a mechanism;
- no claim that kernel selection should be optimized on holdout outcomes.

## Evidence gates

The comparison is only considered operationally valid if:

1. the kernel family is declared before evaluation;
2. the parameterization rule is declared before evaluation;
3. the comparison is run against the same declared controls as the other probes;
4. the statistic does not collapse into a trivial label proxy or a training artifact;
5. the note remains a diagnostic aid, not a verdict engine.

If the kernel choice is not control-separated, the note degrades to a descriptive analogy.

## Failure conditions

This note fails if any of the following happen:

- kernel selection is treated as a transport metric;
- bandwidth tuning is treated as an ontological statement;
- the note is used to infer a physical mechanism;
- the diagnostic is tuned on holdout outcomes;
- the note is used to reopen a frozen physics bridge.

## Audit summary

- Implemented fact: MMD’s sensitivity is dominated by the kernel.
- Design choice: freeze kernel selection as part of the diagnostic contract.
- Heuristic: the median heuristic is the safest default when no domain-specific kernel is already frozen.
- Analogy: kernel choice behaves like a sensitivity dial in the existing discrete operator story.
- Unverified hypothesis: a carefully frozen kernel may make MMD a reliable lightweight screening test before heavier transport diagnostics.

## Related frozen sources

- `docs/notes/equivalences/maximum_mean_discrepancy_kernel_diagnostic_v1.md`
- `docs/notes/equivalences/spectral_address_laplacian_eigenmodes_v1.md`
- `docs/notes/equivalences/spectral_graph_clustering_diagnostic_v1.md`
- `docs/notes/equivalences/wasserstein_distance_transport_cost_diagnostic_v1.md`
- `experiments/geometry_bridge/spectral_diagnostics_summary.md`
