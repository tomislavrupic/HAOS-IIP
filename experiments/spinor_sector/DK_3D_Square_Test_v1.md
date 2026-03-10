# DK 3D Square Test

## Purpose

Test whether the graded 3D Dirac-Kaehler operator squares to the graded Hodge Laplacian on the validated 3D cochain complex.

## Setup

- epsilon: `0.2`
- lattice sizes: `[4, 6, 8]`
- compared low modes: `40`

## Grade-by-grade square diagnostics

| system | relative square residual | mean relative eigenvalue error | `Omega0` block error | `Omega1` block error | `Omega2` block error | `Omega3` block error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 3d_periodic_n4 | 1.388e-16 | 5.390e-02 | 0.000e+00 | 1.602e-16 | 1.602e-16 | 0.000e+00 |
| 3d_periodic_n6 | 0.000e+00 | 5.755e-02 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| 3d_periodic_n8 | 0.000e+00 | 5.378e-02 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 |

## Direct result

- observation: `the graded 3D Dirac-Kaehler square matches the graded Hodge Laplacian to numerical precision across the tested lattice sizes and on every form degree`
- conclusion: `the Stage 7 3D Dirac-Kaehler lift satisfies the grade-by-grade square identity on the validated 3D cochain architecture`

## Artifacts

- results: `data/20260310_180419_DK_3D_square_test.json`
- plots: `plots/20260310_180419_DK_3D_square_vs_Hodge.png`
- timestamp: `20260310_180419`
