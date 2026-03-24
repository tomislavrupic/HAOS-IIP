# Transverse Gap/Stiffness Response

## Purpose

Track how the effective low-band gap and stiffness respond to scalar/background geometry deformations.

## Setup

- sizes: `[12, 16]`
- deformation families: `['radial', 'anisotropic', 'bump']`
- amplitudes: `[0.0, 0.02, 0.05, 0.1]`

## Fit summary

| family | n | eta | `c_eff` | `m_eff` |
| --- | ---: | ---: | ---: | ---: |
| radial | 12 | 0.00 | 37.920592 | -8.275114e-17 |
| radial | 12 | 0.02 | 37.904238 | 7.671609e-05 |
| radial | 12 | 0.05 | 37.858546 | 3.039171e-04 |
| radial | 12 | 0.10 | 37.796904 | 3.542848e-04 |
| radial | 16 | 0.00 | 38.594929 | -2.482534e-16 |
| radial | 16 | 0.02 | 38.583870 | 2.923865e-05 |
| radial | 16 | 0.05 | 38.566889 | 4.453159e-05 |
| radial | 16 | 0.10 | 38.525843 | 4.899001e-05 |
| anisotropic | 12 | 0.00 | 37.920592 | -8.275114e-17 |
| anisotropic | 12 | 0.02 | 37.913232 | 4.398324e-05 |
| anisotropic | 12 | 0.05 | 37.900382 | 9.597146e-05 |
| anisotropic | 12 | 0.10 | 37.859478 | 2.478771e-04 |
| anisotropic | 16 | 0.00 | 38.594929 | -2.482534e-16 |
| anisotropic | 16 | 0.02 | 38.589773 | 1.746538e-05 |
| anisotropic | 16 | 0.05 | 38.579559 | 4.065150e-05 |
| anisotropic | 16 | 0.10 | 38.561330 | 5.394913e-05 |
| bump | 12 | 0.00 | 37.920592 | -8.275114e-17 |
| bump | 12 | 0.02 | 37.918453 | 1.537721e-05 |
| bump | 12 | 0.05 | 37.915982 | 3.446754e-05 |
| bump | 12 | 0.10 | 37.910853 | 6.879248e-05 |
| bump | 16 | 0.00 | 38.594929 | -2.482534e-16 |
| bump | 16 | 0.02 | 38.593952 | 4.751863e-06 |
| bump | 16 | 0.05 | 38.592538 | 9.682135e-06 |
| bump | 16 | 0.10 | 38.591698 | 1.300837e-05 |

## Direct result

- observation: small scalar/background deformations shift the effective transverse fit mainly through smooth stiffness renormalization, with no large unstable gap opening in the tested range
- conclusion: for the anisotropic family at n=16 and eta=0.10, the fitted stiffness is 38.561330 and the fitted gap is 5.394913e-05, consistent with a structured operator-level response rather than chaotic spectral rearrangement

## Artifacts

- results: `data/20260310_165717_transverse_gap_response.json`
- plots: `plots/20260310_165717_transverse_gap_vs_deformation.png`, `plots/20260310_165717_transverse_stiffness_vs_deformation.png`
- timestamp: `20260310_165717`
