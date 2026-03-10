# L1 Large-n Scaling Test

## Purpose

Stress-test the restricted transverse spectrum on larger periodic lattices without changing the operator.

## Setup

- branch: baseline periodic torus
- sizes: `[12, 16, 20, 24, 28, 32]`
- operator: `T = d1* d1 | ker(d0*) intersect (H1)^perp`
- restricted modes: `30`
- phase modes: `30`
- epsilon: `0.2`

## Lowest restricted mode

| `n` | lowest `lambda` | IPR | `||d0* a||` | `||d1 a||` |
| --- | ---: | ---: | ---: | ---: |
| 12 | 0.263337 | 0.000487 | 1.266e-16 | 0.513164 |
| 16 | 0.150761 | 0.000194 | 1.454e-16 | 0.388280 |
| 20 | 0.097277 | 0.000109 | 2.066e-16 | 0.311893 |
| 24 | 0.067853 | 0.000060 | 1.361e-16 | 0.260486 |
| 28 | 0.049985 | 0.000044 | 1.161e-16 | 0.223572 |
| 32 | 0.038336 | 0.000032 | 1.066e-16 | 0.195795 |

## Fit

- `lambda_min ~ 34.9777 / n^1.966` with `R^2 = 1.0000`
- spectral-collapse spread of first 10 scaled modes: `0.0118`

## Direct result

- observation: the lowest restricted transverse eigenvalue keeps decreasing on the baseline torus, the lowest-mode IPR continues to fall, and the first thirty scaled eigenvalues show partial n^2 collapse
- conclusion: the large-n baseline scan remains consistent with a continuum Laplacian-type transverse band in the tested range

## Artifacts

- results: `data/20260310_145017_L1_large_n_scaling_test.json`
- plots: `plots/20260310_145017_transverse_floor_vs_n_large.png`, `plots/20260310_145017_transverse_band_collapse_large.png`, `plots/20260310_145017_transverse_ipr_vs_n_large.png`, `plots/20260310_145017_divergence_curl_phase_large.png`
- timestamp: `20260310_145017`
