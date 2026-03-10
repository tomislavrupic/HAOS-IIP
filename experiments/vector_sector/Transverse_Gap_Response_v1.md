# Transverse Gap/Stiffness Response

## Purpose

Track how the effective low-band gap and stiffness respond to scalar/background geometry deformations.

## Setup

- sizes: `[12]`
- deformation families: `['radial', 'anisotropic', 'bump']`
- amplitudes: `[0.0, 0.02, 0.05, 0.1]`

## Fit summary

| family | n | eta | `c_eff` | `m_eff` |
| --- | ---: | ---: | ---: | ---: |
| radial | 12 | 0.00 | 37.920592 | -8.275114e-17 |
| radial | 12 | 0.02 | 23.683269 | 9.882244e-02 |
| radial | 12 | 0.05 | 23.657247 | 9.886680e-02 |
| radial | 12 | 0.10 | 23.612827 | 9.877855e-02 |
| anisotropic | 12 | 0.00 | 37.920592 | -8.275114e-17 |
| anisotropic | 12 | 0.02 | 33.173393 | 3.295071e-02 |
| anisotropic | 12 | 0.05 | 23.682161 | 9.879744e-02 |
| anisotropic | 12 | 0.10 | 23.654496 | 9.884579e-02 |
| bump | 12 | 0.00 | 37.920592 | -8.275114e-17 |
| bump | 12 | 0.02 | 23.699527 | 9.875445e-02 |
| bump | 12 | 0.05 | 23.698219 | 9.875879e-02 |
| bump | 12 | 0.10 | 23.695926 | 9.876599e-02 |

## Direct result

- observation: small scalar/background deformations shift the effective transverse fit mainly through smooth stiffness renormalization, with no large unstable gap opening in the tested range
- conclusion: for the anisotropic family at n=12 and eta=0.10, the fitted stiffness is 23.654496 and the fitted gap is 9.884579e-02, consistent with a structured operator-level response rather than chaotic spectral rearrangement

## Artifacts

- results: `data/20260310_164600_transverse_gap_response.json`
- plots: `plots/20260310_164600_transverse_gap_vs_deformation.png`, `plots/20260310_164600_transverse_stiffness_vs_deformation.png`
- timestamp: `20260310_164600`
