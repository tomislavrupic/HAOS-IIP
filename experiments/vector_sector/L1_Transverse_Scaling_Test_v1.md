# L1 Transverse Scaling Test

## Purpose

Determine whether the descending restricted transverse floor is organizing like a continuum transverse band or remaining as a defect-localized low mode.

## Setup

- operator unchanged:
  - `L1 = d0 d0* + d1* d1`
  - `T = d1* d1` restricted to `ker(d0*) intersect (H1)^perp`
- substrate branches:
  - baseline periodic torus
  - single cubic puncture
  - line defect
- sizes: `[6, 8, 10, 12, 14, 16]`
- kernel parameter: `epsilon = 0.2`

## Lowest restricted eigenvalue scaling

- baseline periodic torus: `lambda_min ~ 26.7146 / n^1.862` (`R^2 = 0.9996`)
- single cubic puncture: `lambda_min ~ 11.7572 / n^1.578` (`R^2 = 0.9989`)
- line defect: `lambda_min ~ 16.8815 / n^1.698` (`R^2 = 0.9965`)

## Lowest restricted mode diagnostics

### Baseline periodic torus

| `n` | lowest restricted `lambda` | IPR | near-defect fraction | coexact fraction |
| --- | ---: | ---: | ---: | ---: |
| 6 | 0.932912 | 0.003196 | 0.189282 | 1.000000 |
| 8 | 0.563345 | 0.001462 | 0.061848 | 1.000000 |
| 10 | 0.372535 | 0.000808 | 0.050651 | 1.000000 |
| 12 | 0.263337 | 0.000560 | 0.027000 | 1.000000 |
| 14 | 0.195552 | 0.000263 | 0.007191 | 1.000000 |
| 16 | 0.150761 | 0.000191 | 0.008584 | 1.000000 |

### Single cubic puncture

| `n` | lowest restricted `lambda` | IPR | near-defect fraction | coexact fraction |
| --- | ---: | ---: | ---: | ---: |
| 6 | 0.684887 | 0.006034 | 0.547102 | 1.000000 |
| 8 | 0.441271 | 0.002257 | 0.306207 | 1.000000 |
| 10 | 0.317111 | 0.001024 | 0.157783 | 1.000000 |
| 12 | 0.237634 | 0.000548 | 0.077671 | 1.000000 |
| 14 | 0.182925 | 0.000297 | 0.038638 | 1.000000 |
| 16 | 0.144128 | 0.000183 | 0.020121 | 1.000000 |

### Line defect

| `n` | lowest restricted `lambda` | IPR | near-defect fraction | coexact fraction |
| --- | ---: | ---: | ---: | ---: |
| 6 | 0.766501 | 0.006203 | 0.513371 | 1.000000 |
| 8 | 0.514298 | 0.004109 | 0.168196 | 1.000000 |
| 10 | 0.350303 | 0.001593 | 0.078371 | 1.000000 |
| 12 | 0.252087 | 0.000945 | 0.040018 | 1.000000 |
| 14 | 0.189324 | 0.000793 | 0.022190 | 1.000000 |
| 16 | 0.147058 | 0.000533 | 0.013179 | 1.000000 |

## Direct result

- observation: the restricted floor keeps descending while the lowest-mode IPR drops with lattice size in every branch, and the defect branches do not retain a growing fraction of norm near the defect
- conclusion: the descending restricted floor is behaving more like a continuum transverse-band candidate than a defect-pinned low mode in the tested range

## Artifacts

- results: `data/20260310_132323_L1_transverse_scaling.json`
- plots: `plots/20260310_132323_transverse_floor_scaling_loglog.png`, `plots/20260310_132323_transverse_floor_ipr_vs_n.png`, `plots/20260310_132323_transverse_mode_profiles.png`, `plots/20260310_132323_divergence_curl_phase_scaling.png`
- timestamp: `20260310_132323`
