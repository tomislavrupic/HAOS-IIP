# Transverse Continuum Comparison

## Purpose

Compare the low restricted transverse spectrum against the continuum transverse mode counting and spacing on a periodic box.

## Setup

- sizes: `[12, 16, 20]`
- variants: `['baseline', 'puncture', 'line_defect', 'flux_tube']`
- restricted modes: `20`
- continuum reference built from nonzero integer wavevectors with transverse multiplicity two

## Comparison metrics

| branch | `n` | best continuum scale | relative error (first 10) | spacing error (first 10) |
| --- | ---: | ---: | ---: | ---: |
| baseline | 12 | 45.677077 | 0.2189 | 0.0000 |
| baseline | 16 | 46.489347 | 0.2189 | 0.0000 |
| baseline | 20 | 46.869867 | 0.2189 | 0.0000 |
| flux_tube | 12 | 41.959400 | 0.1592 | 0.0000 |
| flux_tube | 16 | 42.810237 | 0.1627 | 0.0000 |
| flux_tube | 20 | 43.216282 | 0.1647 | 0.0000 |
| line_defect | 12 | 36.611429 | 0.0180 | 0.0000 |
| line_defect | 16 | 40.000677 | 0.0428 | 0.0000 |
| line_defect | 20 | 40.460753 | 0.0430 | 0.0000 |
| puncture | 12 | 36.047822 | 0.0526 | 0.0000 |
| puncture | 16 | 37.600157 | 0.0251 | 0.0000 |
| puncture | 20 | 38.372356 | 0.0131 | 0.0000 |

## Direct result

- observation: after n^2 rescaling, the low restricted transverse spectrum aligns with the continuum transverse mode ordering and the mode-spacing pattern remains stable across the tested branches
- conclusion: for the largest baseline case, the first-ten relative spectrum error is 0.219 and the first-ten spacing error is 0.000, consistent with a low continuum transverse operator in the tested range

## Artifacts

- results: `data/20260310_150006_transverse_continuum_comparison.json`
- plots: `plots/20260310_150006_transverse_continuum_comparison.png`, `plots/20260310_150006_transverse_mode_spacing.png`, `plots/20260310_150006_divergence_curl_phase_continuum.png`
- timestamp: `20260310_150006`
