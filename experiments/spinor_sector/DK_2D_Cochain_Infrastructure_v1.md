# DK 2D Cochain Infrastructure

## Purpose

Validate the weighted 2D periodic cochain complex used by the Stage 6 Dirac-Kaehler prototype.

## Setup

- prototype: weighted periodic square cell complex (`2`-torus)
- epsilon: `0.2`
- lattice sizes: `[8, 12, 16]`
- cycle phases: zero (untwisted infrastructure check)

## Cochain identities

| `n` | nodes | edges | faces | `||d1 d0||` | `||delta1 delta2||` |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 64 | 128 | 64 | 0.000e+00 | 0.000e+00 |
| 12 | 144 | 288 | 144 | 0.000e+00 | 0.000e+00 |
| 16 | 256 | 512 | 256 | 0.000e+00 | 0.000e+00 |

## Direct result

- observation: `the weighted periodic 2D cochain complex satisfies the coboundary and adjoint-coboundary identities to numerical precision across the tested lattice sizes`
- conclusion: `the Stage 6 2D cochain infrastructure is internally consistent and admissible for the Dirac-Kaehler square test`

## Artifacts

- results: `data/20260310_174955_DK_2D_cochain_infrastructure.json`
- timestamp: `20260310_174955`
