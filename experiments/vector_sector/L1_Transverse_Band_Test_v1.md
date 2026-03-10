# L1 Transverse Band Test

## Purpose

Determine whether the restricted transverse sector forms a stable low spectral band consistent with continuum scaling, or only an isolated descending mode.

## Setup

- operator unchanged:
  - `L1 = d0 d0* + d1* d1`
  - `T = d1* d1` restricted to `ker(d0*) intersect (H1)^perp`
- branches:
  - baseline periodic torus
  - puncture defect
  - line defect
  - flux-tube defect
- sizes: `[6, 8, 10, 12, 14, 16]`
- kernel parameter: `epsilon = 0.2`
- flux-tube phase: `1.5707963267948966`

## Lowest restricted eigenvalue fits

- baseline periodic torus: `lambda_min ~ 26.7146 / n^1.862` (`R^2 = 0.9996`)
- single cubic puncture: `lambda_min ~ 11.7572 / n^1.578` (`R^2 = 0.9989`)
- line defect: `lambda_min ~ 16.8815 / n^1.698` (`R^2 = 0.9965`)
- flux-tube defect: `lambda_min ~ 15.9007 / n^1.696` (`R^2 = 0.9984`)

## Lowest restricted mode diagnostics

### Baseline periodic torus

| `n` | lowest `lambda` | IPR | near-defect fraction | `||d0* a||` | `||d1 a||` |
| --- | ---: | ---: | ---: | ---: | ---: |
| 6 | 0.932912 | 0.002718 | 0.169853 | 1.568e-16 | 0.965874 |
| 8 | 0.563345 | 0.001371 | 0.059140 | 1.620e-16 | 0.750563 |
| 10 | 0.372535 | 0.000580 | 0.046261 | 1.547e-16 | 0.610357 |
| 12 | 0.263337 | 0.000399 | 0.030616 | 1.426e-16 | 0.513164 |
| 14 | 0.195552 | 0.000243 | 0.007284 | 1.392e-16 | 0.442213 |
| 16 | 0.150761 | 0.000235 | 0.003692 | 1.605e-16 | 0.388280 |

### Single cubic puncture

| `n` | lowest `lambda` | IPR | near-defect fraction | `||d0* a||` | `||d1 a||` |
| --- | ---: | ---: | ---: | ---: | ---: |
| 6 | 0.684887 | 0.004395 | 0.547102 | 1.697e-16 | 0.827579 |
| 8 | 0.441271 | 0.002323 | 0.306207 | 1.815e-16 | 0.664282 |
| 10 | 0.317111 | 0.000974 | 0.157783 | 2.098e-16 | 0.563126 |
| 12 | 0.237634 | 0.000358 | 0.077671 | 1.701e-16 | 0.487477 |
| 14 | 0.182925 | 0.000291 | 0.038638 | 1.837e-16 | 0.427698 |
| 16 | 0.144128 | 0.000133 | 0.020121 | 1.823e-16 | 0.379642 |

### Line defect

| `n` | lowest `lambda` | IPR | near-defect fraction | `||d0* a||` | `||d1 a||` |
| --- | ---: | ---: | ---: | ---: | ---: |
| 6 | 0.766501 | 0.005238 | 0.451382 | 1.488e-16 | 0.875500 |
| 8 | 0.514298 | 0.002688 | 0.058905 | 1.575e-16 | 0.717146 |
| 10 | 0.350303 | 0.002075 | 0.019069 | 1.353e-16 | 0.591864 |
| 12 | 0.252087 | 0.000998 | 0.011612 | 1.395e-16 | 0.502083 |
| 14 | 0.189324 | 0.000690 | 0.010591 | 1.450e-16 | 0.435114 |
| 16 | 0.147058 | 0.000329 | 0.005228 | 1.459e-16 | 0.383482 |

### Flux-tube defect

| `n` | lowest `lambda` | IPR | near-defect fraction | `||d0* a||` | `||d1 a||` |
| --- | ---: | ---: | ---: | ---: | ---: |
| 6 | 0.738660 | 0.004831 | 0.340171 | 2.402e-16 | 0.859453 |
| 8 | 0.477984 | 0.001771 | 0.227188 | 1.970e-16 | 0.691364 |
| 10 | 0.329125 | 0.000835 | 0.153581 | 2.092e-16 | 0.573694 |
| 12 | 0.238671 | 0.000462 | 0.101545 | 1.925e-16 | 0.488540 |
| 14 | 0.180326 | 0.000284 | 0.068613 | 1.708e-16 | 0.424648 |
| 16 | 0.140752 | 0.000187 | 0.058057 | 1.758e-16 | 0.375169 |

## Direct result

- observation: the first twenty restricted transverse eigenvalues show a partial n^2 spectral collapse while the lowest-mode IPR and defect concentration both fall with lattice size
- conclusion: the restricted transverse sector is organizing as a stable low spectral band consistent with continuum scaling in the tested range

## Artifacts

- results: `data/20260310_133322_L1_transverse_band_scan.json`
- plots: `plots/20260310_133322_transverse_band_collapse.png`, `plots/20260310_133322_transverse_floor_vs_n.png`, `plots/20260310_133322_transverse_ipr_vs_n.png`, `plots/20260310_133322_divergence_curl_phase_band.png`, `plots/20260310_133322_transverse_mode_profiles.png`
- timestamp: `20260310_133322`
