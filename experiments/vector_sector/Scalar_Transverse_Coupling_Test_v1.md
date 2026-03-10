# Scalar-Transverse Coupling Test

## Purpose

Measure how the restricted transverse spectrum deforms under small scalar/background geometry perturbations without changing operator definitions.

## Setup

- sizes: `[12]`
- variants: `['baseline', 'puncture', 'line_defect']`
- deformation families: `['radial', 'anisotropic', 'bump']`
- amplitudes: `[0.02, 0.05, 0.1]`

## Response summary

| branch | family | eta | deformed `lambda_1` | shift `Delta lambda_1` | mean overlap | lowest-mode IPR |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| baseline | radial | 0.02 | 0.263235 | -0.000103 | 0.4111 | 0.000499 |
| baseline | radial | 0.05 | 0.263039 | -0.000298 | 0.4079 | 0.000415 |
| baseline | radial | 0.10 | 0.262605 | -0.000733 | 0.4452 | 0.000409 |
| baseline | anisotropic | 0.02 | 0.263245 | -0.000092 | 0.4229 | 0.000672 |
| baseline | anisotropic | 0.05 | 0.263103 | -0.000234 | 0.4185 | 0.000436 |
| baseline | anisotropic | 0.10 | 0.262858 | -0.000479 | 0.4121 | 0.000466 |
| baseline | bump | 0.02 | 0.263323 | -0.000015 | 0.4467 | 0.000317 |
| baseline | bump | 0.05 | 0.263300 | -0.000037 | 0.4644 | 0.000412 |
| baseline | bump | 0.10 | 0.263262 | -0.000076 | 0.4469 | 0.000405 |
| puncture | radial | 0.02 | 0.237654 | 0.000019 | 0.7522 | 0.000462 |
| puncture | radial | 0.05 | 0.237598 | -0.000036 | 0.7397 | 0.000383 |
| puncture | radial | 0.10 | 0.237283 | -0.000351 | 0.7558 | 0.000446 |
| puncture | anisotropic | 0.02 | 0.237629 | -0.000006 | 0.7490 | 0.000587 |
| puncture | anisotropic | 0.05 | 0.237620 | -0.000014 | 0.7478 | 0.000587 |
| puncture | anisotropic | 0.10 | 0.237473 | -0.000161 | 0.8125 | 0.000382 |
| puncture | bump | 0.02 | 0.237614 | -0.000021 | 0.7483 | 0.000487 |
| puncture | bump | 0.05 | 0.237582 | -0.000052 | 0.7501 | 0.000347 |
| puncture | bump | 0.10 | 0.237529 | -0.000105 | 0.7289 | 0.000444 |
| line_defect | radial | 0.02 | 0.252133 | 0.000045 | 0.7489 | 0.001260 |
| line_defect | radial | 0.05 | 0.252101 | 0.000013 | 0.7512 | 0.001260 |
| line_defect | radial | 0.10 | 0.251784 | -0.000304 | 0.7571 | 0.001261 |
| line_defect | anisotropic | 0.02 | 0.252086 | -0.000001 | 0.7263 | 0.000867 |
| line_defect | anisotropic | 0.05 | 0.252084 | -0.000003 | 0.7380 | 0.000867 |
| line_defect | anisotropic | 0.10 | 0.252033 | -0.000054 | 0.7297 | 0.000867 |
| line_defect | bump | 0.02 | 0.252084 | -0.000003 | 0.7492 | 0.001260 |
| line_defect | bump | 0.05 | 0.252080 | -0.000008 | 0.7492 | 0.001260 |
| line_defect | bump | 0.10 | 0.252070 | -0.000017 | 0.7492 | 0.001259 |

## Direct result

- observation: small scalar/background deformations shift the restricted transverse spectrum smoothly, while the matched low-mode overlaps remain moderate to high depending on branch and deformation family
- conclusion: for the baseline anisotropic branch at eta=0.10, the lowest restricted eigenvalue shifts by -0.000479 with mean matched overlap 0.4121, indicating a measurable but bounded operator-level spectral response

## Artifacts

- results: `data/20260310_164041_scalar_transverse_coupling_test.json`
- plots: `plots/20260310_164041_scalar_transverse_eigenvalue_shifts.png`, `plots/20260310_164041_scalar_transverse_mode_overlap.png`, `plots/20260310_164041_scalar_transverse_deformation_response.png`, `plots/20260310_164041_divergence_curl_phase_scalar_transverse.png`
- timestamp: `20260310_164041`
