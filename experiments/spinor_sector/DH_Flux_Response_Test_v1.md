# DH Flux Response Test

## Purpose

Measure how the Stage 5 Dirac-type spectrum shifts under background torus flux in the 2D prototype.

## Setup

- prototype: full chiral 4-spinor lift
- lattice: 2D periodic torus with `n=6`
- flux quanta: `[0, 1, 2, 3]`
- tracked positive branch size: `12`

## Flux-flow summary

| flux quanta | plaquette angle | first positive mode | last tracked positive mode |
| ---: | ---: | ---: | ---: |
| 0 | 0.000000 | 0.499337 | 1.197020 |
| 1 | 0.174533 | 0.003144 | 0.894280 |
| 2 | 0.349066 | 0.036391 | 0.838798 |
| 3 | 0.523599 | 0.091980 | 1.084837 |

## Direct result

- observation: `the low positive DH branch shifts smoothly as torus flux is increased in the 2D prototype`
- conclusion: `the prototype has smooth 2D background-flux spectral flow, but it remains rejected overall because the square-root criterion against the scalar operator failed`

## Artifacts

- results: `data/20260310_172754_DH_flux_response_test.json`
- plots: `plots/20260310_172754_DH_flux_spectral_flow.png`
- timestamp: `20260310_172754`
