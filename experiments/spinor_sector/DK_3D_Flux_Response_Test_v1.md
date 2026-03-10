# DK 3D Flux Response Test

## Purpose

Measure how the 3D Dirac-Kaehler spectrum responds to flat torus-cycle holonomy while preserving the graded square identity.

## Setup

- lattice: 3D periodic cell complex with `n=4`
- epsilon: `0.2`
- x-cycle phases: `[0.0, 0.2617993877991494, 0.5235987755982988, 0.7853981633974483]`
- tracked positive branch size: `12`

## Flux diagnostics

| x-cycle phase | relative square residual | first positive mode | last tracked positive mode | mean `Omega1` fraction |
| ---: | ---: | ---: | ---: | ---: |
| 0.000000 | 1.388e-16 | 1.307934 | 1.307934 | 0.386 |
| 0.261799 | 1.388e-16 | 0.060520 | 1.309333 | 0.386 |
| 0.523599 | 1.388e-16 | 0.120976 | 1.313517 | 0.374 |
| 0.785398 | 1.388e-16 | 0.181302 | 1.320440 | 0.367 |

## Direct result

- observation: `the low 3D Dirac-Kaehler branch shifts smoothly under flat x-cycle holonomy while the graded square residual remains at numerical precision`
- conclusion: `the Stage 7 3D Dirac-Kaehler lift preserves its square identity under mild torus-cycle holonomy and shows stable low spectral flow`

## Artifacts

- results: `data/20260310_180436_DK_3D_flux_response_test.json`
- plots: `plots/20260310_180436_DK_3D_flux_spectral_flow.png, plots/20260310_180436_DK_3D_flux_square_residual.png`
- timestamp: `20260310_180436`
