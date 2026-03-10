# DH Square-Root Test

## Purpose

Test whether the Stage 5 Dirac-type lift behaves as a square root of the scalar operator on the same kernel substrate.

## Setup

- prototype: full chiral 4-spinor lift
- 2D baseline: `n=6` periodic torus
- 3D baseline: `n=8` periodic lattice
- scalar comparison: repeated low spectrum of `L0`

## Square-root diagnostics

| system | mean relative eigenvalue error | max relative eigenvalue error | eigenvalue correlation | operator relative Frobenius error |
| --- | ---: | ---: | ---: | ---: |
| 2d_periodic_n6 | 3.065e-01 | 7.321e-01 | 0.943113 | 5.477e-01 |
| 3d_periodic_n8 | 7.612e-01 | 1.000e+00 | 0.934254 | 6.547e-01 |

## Direct result

- observation: `the chiral Dirac prototype tracks the low scalar spectrum only qualitatively: eigenvalue correlations remain high, but the relative spectral and operator errors stay large on both prototype grids`
- conclusion: `the current Stage 5 prototype does not satisfy the square-root criterion tightly enough; under the prompt rule, this lift is rejected in its present form`

## Artifacts

- results: `data/20260310_172136_DH_square_root_test.json`
- plots: `plots/20260310_172136_DH2_vs_L0_spectrum.png`
- timestamp: `20260310_172136`
