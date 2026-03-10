# Transverse Effective Operator Fit

## Purpose

Fit the low restricted transverse spectrum on the clean periodic branch to the minimal effective form `lambda_k ≈ c_eff |q_k|^2 + m_eff`.

## Setup

- branch: baseline periodic torus
- sizes: `[12, 16, 20, 24]`
- restricted modes: `24`
- fit modes: `2`

## Fitted parameters

| `n` | `c_eff` | `m_eff` | `R^2` |
| --- | ---: | ---: | ---: |
| 12 | 37.920592 | -8.275114e-17 | 1.00000 |
| 16 | 38.594929 | -2.482534e-16 | 1.00000 |
| 20 | 38.910833 | -2.068778e-17 | 1.00000 |
| 24 | 39.083446 | 0.000000e+00 | 1.00000 |

## Direct result

- observation: the low restricted spectrum on the clean periodic branch is well fit by a linear function of effective momentum squared across the tested lattice sizes
- conclusion: the fitted gap stays small, from -8.275e-17 at n=12 to 0.000e+00 at n=24, while the effective stiffness varies only mildly across the tested sizes, consistent with a nearly gapless continuum transverse branch in the tested range

## Artifacts

- results: `data/20260310_155713_transverse_effective_operator_fit.json`
- plots: `plots/20260310_155713_effective_operator_fit.png`, `plots/20260310_155713_effective_gap_vs_n.png`, `plots/20260310_155713_effective_stiffness_vs_n.png`, `plots/20260310_155713_divergence_curl_phase_effective_fit.png`
- timestamp: `20260310_155713`
