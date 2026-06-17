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
- Neither note changes the GEO-HIDDEN-01 verdict ceiling.
- The metric family is recorded for diagnostic completeness; no code path or
  holdout verdict is changed here.

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

Unverified hypothesis:
- The remaining transformation gap may be caused by geometry that is too weak
  in the current synthetic family, but the benchmark does not establish that.

## Boundary

The diagnostic layer is stronger than before, but the open verdict remains
correct: the probe still does not recover transformation classes robustly on
holdout.
