# Spectral Diagnostics Summary

This note records the current hardening pass for `GEO-HIDDEN-01` without
changing its verdict ceiling.

## What changed

The hidden-geometry probe now reports three spectral views side by side:

| Variant | Role | Current reading |
| --- | --- | --- |
| Unnormalized Laplacian | Baseline spectral path | Useful, but sensitive to degree variation |
| Symmetric normalized Laplacian (`L_sym`) | Primary normalized path | More stable than the unnormalized baseline, but still below the frozen transformation threshold |
| Random-walk normalized Laplacian (`L_rw`) | Alternative normalized path | Slightly better than `L_sym` on the current holdout, but still not enough to close the gap |
| Cheeger sweep-cut estimate | Bottleneck diagnostic | Detects a sparse cut, but does not recover the transformation classes |

## Frozen values

The current frozen result remains:

- `transform_accuracy`: `0.500000`
- `fiedler_transform_accuracy`: `0.250000`
- `fiedler_sym_accuracy`: `0.250000`
- `fiedler_rw_accuracy`: `0.312500`
- `fiedler_sign_stability`: `0.513889`
- `cheeger_conductance`: `0.304962`

## Transport-layer companion

The curvature notes now separate the transport primitive from the curvature
comparison class:

| Layer | Role | Current reading |
| --- | --- | --- |
| Wasserstein distance | Transport cost between frozen neighborhood measures | Explicit primitive for local transport comparisons |
| Ollivier-Ricci curvature | Neighborhood overlap / bridge-vs-community diagnostic | Companion diagnostic built on Wasserstein transport cost |

### Wasserstein variants

| Metric | Order | Diagnostic reading | Practical note |
| --- | --- | --- | --- |
| `W1` / Earth Mover's distance | 1 | Primary transport cost used by Ollivier-style curvature | Most direct fit for local neighborhood overlap |
| `W2` | 2 | Smooth quadratic transport cost | Potentially useful for embedding-style comparisons, but heavier than `W1` |
| Sliced `W1` / `W2` | Approx. projected | Faster approximation by 1D projections | Good for scalable diagnostics, but less exact than full transport |
| Sinkhorn-regularized transport | Regularized | Entropic approximation to transport cost | Practical for larger graphs, but introduces regularization choice |

Interpretation:
- Wasserstein distance is the cost functional underneath the local transport comparison.
- Ollivier-Ricci turns that cost into a curvature-style bridge diagnostic.
- Sinkhorn is recorded as the practical approximation layer when exact transport
  is too expensive to use edge by edge.
- Neither note changes the GEO-HIDDEN-01 verdict ceiling.
- The metric family is recorded for diagnostic completeness; no code path or
  holdout verdict is changed here.

### Kernel companion

| Statistic | Role | Current reading |
| --- | --- | --- |
| Maximum Mean Discrepancy (`MMD`) | Kernel-side distribution test | Lightweight comparison for frozen embeddings and curvature summaries |

### MMD kernel selection

| Kernel family | Diagnostic reading | Practical note |
| --- | --- | --- |
| Gaussian / RBF | Default kernel for spectral embeddings | Median heuristic is the safest frozen starting point |
| Laplace | Robust to outliers | Good for noisy curvature summaries |
| Rational quadratic | Multi-scale discrepancy | Useful when the feature scale is mixed |
| Linear | Fast first-pass check | Only sees coarse mean shifts |

Interpretation:
- MMD is the kernel-side companion to the transport ladder, not part of the
  transport metric family itself.
- It is useful for quick distribution separation on spectral embeddings or
  curvature histograms.
- Kernel choice is the main sensitivity dial for MMD, so it should be frozen
  alongside the diagnostic contract rather than tuned informally.
- It does not change the GEO-HIDDEN-01 verdict ceiling.

## Interpretation

Implemented fact:
- The normalized spectrum and Cheeger sweep add real diagnostic structure.

Design choice:
- `L_sym` is the primary normalized path; `L_rw` is reported as a secondary
  comparison rather than a replacement verdict gate.

Heuristic:
- The Cheeger sweep is a local bottleneck estimate from the normalized Fiedler
  vector, not a proof of recoverability.
- Wasserstein and Ollivier-Ricci now sit beside the spectral probes as local
  transport diagnostics, not as an escalation in claim scope.

Analogy:
- `L_sym` / `L_rw` behave like alternative views on the same low-mode geometry,
  but they are not equivalent to a successful transformation recovery.
- Wasserstein distance behaves like the cost primitive beneath neighborhood
  transport comparison, while Ollivier-Ricci behaves like the curvature-style
  reading of that cost.
- `W1` is the direct local choice for neighborhood transport; `W2`, sliced
  variants, and Sinkhorn are recorded as practical approximations or adjacent
  comparisons, not as stronger claims.
- Sinkhorn is the speed-oriented approximation layer, not a replacement for
  exact Wasserstein or a stronger transport claim.
- MMD is the kernel-side separation score beside the transport family, useful
  when the question is distributional difference rather than mass transport.
- For kernel choice, Gaussian / RBF is the default; Laplace and rational
  quadratic are recorded as practical alternatives.

Unverified hypothesis:
- The remaining transformation gap may be caused by geometry that is too weak
  in the current synthetic family, but the benchmark does not establish that.

## Boundary

The diagnostic layer is stronger than before, but the open verdict remains
correct: the probe still does not recover transformation classes robustly on
holdout.
