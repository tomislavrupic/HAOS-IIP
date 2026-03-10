# DK 3D Cochain Infrastructure

## Purpose

Validate the weighted 3D periodic cochain complex used by the Stage 7 Dirac-Kaehler prototype.

## Setup

- prototype: weighted periodic cubic cell complex (`3`-torus)
- epsilon: `0.2`
- lattice sizes: `[4, 6, 8]`
- cycle phases: zero (untwisted infrastructure check)

## Cochain identities

| `n` | nodes | edges | faces | volumes | `||d1 d0||` | `||d2 d1||` | `||delta1 delta2||` | `||delta2 delta3||` |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 64 | 192 | 192 | 64 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| 6 | 216 | 648 | 648 | 216 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| 8 | 512 | 1536 | 1536 | 512 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 0.000e+00 |

## Direct result

- observation: `the weighted periodic 3D cochain complex satisfies the nilpotency and adjoint-nilpotency identities to numerical precision across the tested lattice sizes`
- conclusion: `the Stage 7 3D cochain infrastructure is internally consistent and admissible for the graded square test`

## Artifacts

- results: `data/20260310_180412_DK_3D_cochain_infrastructure.json`
- timestamp: `20260310_180412`
