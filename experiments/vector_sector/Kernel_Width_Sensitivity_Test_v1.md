# Kernel Width Sensitivity Test

## Purpose

Measure how the restricted transverse band changes as the kernel width varies at fixed grid size.

## Setup

- branch: baseline periodic torus
- `n = 16`
- `epsilon_k = c_epsilon h^2`
- `c_epsilon = [0.25, 0.5, 1.0, 1.5]`

## Summary metrics

| `c_epsilon` | lowest `lambda` | `lambda_10 - lambda_1` | lowest-mode IPR |
| --- | ---: | ---: | ---: |
| 0.25 | 0.020604 | 0.020604 | 0.000183 |
| 0.5 | 0.056006 | 0.056006 | 0.000174 |
| 1.0 | 0.092339 | 0.092339 | 0.000177 |
| 1.5 | 0.109085 | 0.109085 | 0.000193 |

## Direct result

- observation: the restricted transverse spectrum remains present across the tested kernel-width range, while the lowest eigenvalue, first-band spread, and lowest-mode IPR shift smoothly with c_epsilon
- conclusion: within the tested local-kernel range, the restricted transverse band is robust to moderate kernel-width changes

## Artifacts

- results: `data/20260310_144513_kernel_width_sensitivity_test.json`
- plots: `plots/20260310_144513_transverse_band_vs_kernel_width.png`, `plots/20260310_144513_divergence_curl_phase_kernel_width.png`
- timestamp: `20260310_144513`
