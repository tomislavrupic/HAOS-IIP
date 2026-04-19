# Transverse Sector Window Scan

## Purpose

Push the CP2 sector split deeper into the full `L1` spectrum and test whether defect-localized backgrounds pull a genuinely coexact branch down or localize it in a controlled way.

## Setup

- sizes: `[12, 16, 20, 24, 28]`
- variants: `['baseline', 'puncture', 'line_defect', 'flux_tube']`
- restricted modes: `20`
- full-mode scan count: `40`
- coexact threshold: `0.8`

## Sector-window metrics

| branch | `n` | lowest full harmonic frac | candidate coexact frac | candidate `j` | meets threshold? | scaled gap | candidate near-defect frac | restricted / candidate | support |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| baseline | 12 | 1.000 | 0.864 | 34 | yes | 113.762 | 0.017 | 0.333 | coexact z-biased circulation |
| baseline | 16 | 1.000 | 0.861 | 11 | yes | 38.595 | 0.012 | 1.000 | coexact x-biased circulation |
| baseline | 20 | 1.000 | 0.812 | 10 | yes | 38.911 | 0.004 | 1.000 | coexact z-biased circulation |
| baseline | 24 | 1.000 | 0.842 | 6 | yes | 39.083 | 0.002 | 1.000 | coexact x-biased circulation |
| baseline | 28 | 1.000 | 0.806 | 19 | yes | 78.376 | 0.002 | 0.500 | coexact z-biased circulation |
| flux_tube | 12 | 1.000 | 0.883 | 6 | yes | 37.921 | 0.042 | 0.906 | coexact y-biased circulation |
| flux_tube | 16 | 1.000 | 0.807 | 4 | yes | 38.595 | 0.018 | 0.934 | coexact z-biased circulation |
| flux_tube | 20 | 1.000 | 0.907 | 6 | yes | 38.911 | 0.015 | 0.949 | coexact z-biased circulation |
| flux_tube | 24 | 1.000 | 0.812 | 6 | yes | 39.083 | 0.014 | 0.959 | coexact y-biased circulation |
| flux_tube | 28 | 1.000 | 1.000 | 11 | yes | 39.192 | 0.011 | 0.965 | coexact x-biased circulation |
| line_defect | 12 | 1.000 | 1.000 | 3 | yes | 36.301 | 0.011 | 1.000 | coexact z-biased circulation |
| line_defect | 16 | 1.000 | 1.000 | 3 | yes | 37.647 | 0.006 | 1.000 | coexact z-biased circulation |
| line_defect | 20 | 1.000 | 1.000 | 3 | yes | 38.296 | 0.001 | 1.000 | coexact z-biased circulation |
| line_defect | 24 | 1.000 | 1.000 | 3 | yes | 38.654 | 0.001 | 1.000 | coexact z-biased circulation |
| line_defect | 28 | 1.000 | 1.000 | 3 | yes | 38.871 | 0.000 | 1.000 | coexact z-biased circulation |
| puncture | 12 | 1.000 | 1.000 | 3 | yes | 34.219 | 0.078 | 1.000 | coexact y-biased circulation |
| puncture | 16 | 1.000 | 1.000 | 3 | yes | 36.897 | 0.020 | 1.000 | coexact y-biased circulation |
| puncture | 20 | 1.000 | 1.000 | 3 | yes | 38.032 | 0.006 | 1.000 | coexact x-biased circulation |
| puncture | 24 | 1.000 | 1.000 | 3 | yes | 38.576 | 0.003 | 1.000 | coexact z-biased circulation |
| puncture | 28 | 1.000 | 1.000 | 3 | yes | 38.870 | 0.001 | 1.000 | coexact z-biased circulation |

## Direct result

- observation: the true low full-spectrum floor remains harmonic across the scanned backgrounds, but the direct harmonic-to-coexact gap stays finite and defect backgrounds can reduce that separation while pulling the coexact candidate deeper into the low window
- conclusion: at the largest size, the clean baseline still has harmonic floor fraction 1.000, reaches its first coexact candidate at j=19, and keeps a scaled harmonic-to-coexact gap of 78.376; by contrast, puncture reaches its coexact candidate at j=3, puncture gives the smallest scaled gap at 38.870, and flux_tube shows the strongest candidate defect localization with near-defect fraction 0.011, so CP2c says defects compress the harmonic-to-coexact separation but still do not collapse it into a clean propagating vector floor

## Artifacts

- results: `data/20260412_144123_transverse_sector_window_scan.json`
- plots: `plots/20260412_144123_transverse_sector_candidate_index.png`, `plots/20260412_144123_transverse_sector_localization.png`, `plots/20260412_144123_transverse_sector_gap.png`, `plots/20260412_144123_transverse_sector_candidate_metrics.png`
- timestamp: `20260412_144123`
