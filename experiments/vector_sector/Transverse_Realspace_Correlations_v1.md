# Transverse Real-Space Correlations

## Purpose

Measure whether the lowest restricted transverse modes develop extended, approximately isotropic correlations on the periodic branch.

## Setup

- branch: baseline periodic torus
- sizes: `[12, 16]`
- pair samples per size: `8000`
- bins: `24`

## Correlation summary

| `n` | correlation length | anisotropy ratio | lowest-mode IPR |
| --- | ---: | ---: | ---: |
| 12 | 0.094035 | 4.6129 | 0.000531 |
| 16 | 0.088414 | 1.1975 | 0.000116 |

## Direct result

- observation: the lowest restricted transverse mode shows extended two-point correlations while directional anisotropy remains modest across the tested sizes
- conclusion: the correlation length remains of the same order, from 0.094 at n=12 to 0.088 at n=16, while the anisotropy ratio stays moderate at 1.197 for the largest size, consistent with an extended continuum-like branch in the tested range

## Artifacts

- results: `data/20260310_155153_transverse_realspace_correlations.json`
- plots: `plots/20260310_155153_transverse_correlation_decay.png`, `plots/20260310_155153_transverse_anisotropy.png`, `plots/20260310_155153_divergence_curl_phase_correlations.png`
- timestamp: `20260310_155153`
