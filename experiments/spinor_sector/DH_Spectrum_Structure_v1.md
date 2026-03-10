# DH Spectrum Structure

## Purpose

Measure the low spectral structure of the Stage 5 Dirac-type lift.

## Setup

- prototype: full chiral 4-spinor lift
- spectrum window: `50` modes nearest zero
- systems: 2D periodic torus and 3D periodic lattice

## Spectral diagnostics

| system | pairing error | lowest | mean IPR | max IPR |
| --- | ---: | ---: | ---: | ---: |
| 2d_periodic_n6 | 1.281e-02 | 1.069e-16 | 1.172e-02 | 1.959e-02 |
| 3d_periodic_n8 | 3.844e-01 | 1.126e-17 | 1.037e-03 | 2.676e-03 |

## Direct result

- observation: `the low DH spectrum is sign-symmetric and extended in the 2D prototype, while the 3D prototype remains extended but shows only partial near-zero pairing`
- conclusion: `the prototype shows qualitative Dirac-like spectral structure, but this does not overcome the failed square-root test against the scalar operator`

## Artifacts

- results: `data/20260310_173100_DH_spectrum_structure.json`
- plots: `plots/20260310_173100_DH_spectrum.png, plots/20260310_173100_DH_ipr_vs_mode.png`
- timestamp: `20260310_173100`
