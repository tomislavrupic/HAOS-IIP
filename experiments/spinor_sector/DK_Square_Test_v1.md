# DK Square Test

## Purpose

Test whether the graded Dirac-Kaehler operator squares to the graded Hodge Laplacian on the validated 2D cochain complex.

## Setup

- epsilon: `0.2`
- lattice sizes: `[8, 12, 16]`
- compared low modes: `40`

## Square diagnostics

| system | relative square residual | mean relative eigenvalue error | `Omega0` block error | `Omega1` block error | `Omega2` block error |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2d_periodic_n8 | 0.000e+00 | 4.777e-03 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| 2d_periodic_n12 | 0.000e+00 | 9.922e-03 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| 2d_periodic_n16 | 0.000e+00 | 9.278e-03 | 0.000e+00 | 0.000e+00 | 0.000e+00 |

## Direct result

- observation: `the graded Dirac-Kaehler square matches the graded Hodge Laplacian to numerical precision across the tested 2D lattice sizes`
- conclusion: `the Stage 6 Dirac-Kaehler lift satisfies the defining square identity on the validated 2D cochain architecture`

## Artifacts

- results: `data/20260310_175001_DK_square_test.json`
- plots: `plots/20260310_175001_DK_square_vs_Hodge.png`
- timestamp: `20260310_175001`
