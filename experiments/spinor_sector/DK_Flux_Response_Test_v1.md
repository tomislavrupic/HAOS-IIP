# DK Flux Response Test

## Purpose

Measure how the 2D Dirac-Kaehler spectrum responds to flat torus-cycle holonomy while preserving the graded square identity.

## Setup

- lattice: 2D periodic cell complex with `n=12`
- epsilon: `0.2`
- cycle phases: `[0.0, 0.2617993877991494, 0.5235987755982988, 0.7853981633974483]`
- tracked positive branch size: `12`

## Flux diagnostics

| cycle phase | relative square residual | first positive mode | last tracked positive mode | mean `Omega1` fraction |
| ---: | ---: | ---: | ---: | ---: |
| 0.000000 | 0.000e+00 | 0.497265 | 0.725724 | 0.500 |
| 0.261799 | 0.000e+00 | 0.021628 | 0.711084 | 0.500 |
| 0.523599 | 0.000e+00 | 0.043253 | 0.696725 | 0.500 |
| 0.785398 | 0.000e+00 | 0.064873 | 0.682673 | 0.500 |

## Direct result

- observation: `the low Dirac-Kaehler branch shifts smoothly under flat torus-cycle holonomy while the graded square residual remains at numerical precision`
- conclusion: `the Stage 6 2D Dirac-Kaehler lift preserves its square identity under mild background cycle holonomy and shows stable spectral flow`

## Artifacts

- results: `data/20260310_175553_DK_flux_response_test.json`
- plots: `plots/20260310_175553_DK_flux_spectral_flow.png, plots/20260310_175553_DK_flux_square_residual.png`
- timestamp: `20260310_175553`
