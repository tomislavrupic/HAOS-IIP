# DH Hermiticity Test

## Purpose

Verify Hermiticity of the Stage 5 Dirac-type lift on the prototype periodic grids.

## Setup

- prototype: full chiral 4-spinor lift
- epsilon: `0.2`
- 2D periodic torus: `n=6`
- 3D periodic lattice: `n=8`

## Hermiticity summary

| system | total dimension | `||D_H - D_H^dagger||` |
| --- | ---: | ---: |
| 2d_periodic_n6 | 144 | 0.000e+00 |
| 3d_periodic_n8 | 2048 | 0.000e+00 |

## Direct result

- observation: `the chiral Dirac-type lift is Hermitian to numerical precision on both the 2D and 3D prototype grids`
- conclusion: `the Stage 5 prototype passes the Hermiticity gate and is admissible for the square-root and spectral tests`

## Artifacts

- results: `data/20260310_172558_DH_hermiticity_test.json`
- timestamp: `20260310_172558`
