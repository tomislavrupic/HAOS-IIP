# DK Spectrum Structure

## Purpose

Measure the low graded spectral structure of the Stage 6 Dirac-Kaehler operator.

## Setup

- lattice: 2D periodic cell complex with `n=12`
- epsilon: `0.2`
- spectrum window: `60` modes nearest zero

## Low-spectrum diagnostics

- pairing error: `5.380e-03`
- lowest `|lambda|`: `1.089e-17`
- mean IPR: `3.703e-03`
- max IPR: `4.763e-03`

## Degeneracy summary

| `|lambda|` | degeneracy |
| ---: | ---: |
| 0.000000 | 4 |
| 0.497265 | 1 |
| 0.513164 | 14 |
| 0.541547 | 1 |
| 0.704585 | 1 |
| 0.725724 | 14 |
| 0.728703 | 1 |
| 0.985738 | 1 |

## Direct result

- observation: `the low graded Dirac-Kaehler spectrum shows clean sign symmetry, low near-zero modes, and stable weight distribution across 0-, 1-, and 2-form sectors`
- conclusion: `the Stage 6 2D prototype exhibits a stable graded spectral structure consistent with a Dirac-Kaehler lift on the validated cochain complex`

## Artifacts

- results: `data/20260310_175158_DK_spectrum_structure.json`
- plots: `plots/20260310_175158_DK_spectrum.png, plots/20260310_175158_DK_mode_grading.png, plots/20260310_175158_DK_ipr_vs_mode.png`
- timestamp: `20260310_175158`
