# Transverse Continuum Comparison

## Purpose

Compare the low restricted transverse spectrum against the continuum transverse mode counting and spacing on a periodic box, while explicitly checking where that projected branch sits relative to the full low `L1` sector split.

## Setup

- sizes: `[12, 16, 20]`
- variants: `['baseline', 'puncture', 'line_defect', 'flux_tube']`
- restricted modes: `20`
- coexact threshold for full-mode scan: `0.8`
- continuum reference built from nonzero integer wavevectors with transverse multiplicity two
- restricted band is built after exact/harmonic projection; the full-mode columns below check whether the actual low `L1` floor already sits in that same coexact sector

## Comparison metrics

| branch | `n` | best continuum scale | relative error (first 10) | spacing error (first 10) | lowest full harmonic frac | lowest full coexact frac | coexact candidate `j` | meets threshold? | restricted / candidate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline | 12 | 45.677077 | 0.2189 | 0.0000 | 1.000 | 0.000 | 10 | no | 1.000 |
| baseline | 16 | 46.489347 | 0.2189 | 0.0000 | 1.000 | 0.000 | 3 | yes | 1.000 |
| baseline | 20 | 46.869867 | 0.2189 | 0.0000 | 1.000 | 0.000 | 8 | no | 1.000 |
| flux_tube | 12 | 41.959400 | 0.1592 | 0.0000 | 1.000 | 0.000 | 8 | yes | 0.906 |
| flux_tube | 16 | 42.810237 | 0.1627 | 0.0000 | 1.000 | 0.000 | 8 | yes | 0.933 |
| flux_tube | 20 | 43.216282 | 0.1647 | 0.0000 | 1.000 | 0.000 | 8 | yes | 0.949 |
| line_defect | 12 | 36.611429 | 0.0180 | 0.0000 | 1.000 | 0.000 | 3 | yes | 1.000 |
| line_defect | 16 | 37.874858 | 0.0103 | 0.0000 | 1.000 | 0.000 | 3 | yes | 1.000 |
| line_defect | 20 | 38.460608 | 0.0066 | 0.0000 | 1.000 | 0.000 | 3 | yes | 1.000 |
| puncture | 12 | 36.047822 | 0.0526 | 0.0000 | 1.000 | 0.000 | 3 | yes | 1.000 |
| puncture | 16 | 37.600157 | 0.0251 | 0.0000 | 1.000 | 0.000 | 3 | yes | 1.000 |
| puncture | 20 | 38.372356 | 0.0131 | 0.0000 | 1.000 | 0.000 | 3 | yes | 1.000 |

## Direct result

- observation: after n^2 rescaling, the low restricted transverse spectrum still follows the continuum ordering across the tested branches, but the full low L1 floor remains sector-split from that branch rather than coinciding with it
- conclusion: for the largest baseline case, the first-ten restricted spectrum error is 0.219, the lowest full mode carries harmonic fraction 1.000, and none of the first 20 full modes reaches the coexact threshold 0.80, and the best low coexact candidate within that window sits at index 8; the restricted floor tracks that candidate with ratio 1.000, so CP2 now makes the harmonic/coexact split explicit without yet establishing a clean low continuum vector band

## Artifacts

- results: `data/20260412_141530_transverse_continuum_comparison.json`
- plots: `plots/20260412_141530_transverse_continuum_comparison.png`, `plots/20260412_141530_transverse_mode_spacing.png`, `plots/20260412_141530_divergence_curl_phase_continuum.png`, `plots/20260412_141530_transverse_sector_split.png`
- timestamp: `20260412_141530`
