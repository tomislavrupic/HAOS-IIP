# Transverse Clean Holonomy Symmetry Labels

## Purpose

Execute `CB3` in bounded form on the already-resolved clean holonomy family. The goal is not a full point-group irrep package, but a stable momentum-family / circulation-axis label on the rescued active coexact window.

## Setup

- branch: `baseline`
- sizes: `[12, 16, 20]`
- holonomy phase: `0.05`
- restricted modes: `8`
- resolved active window: start `2`, size `4`

## Window summary

| `n` | first-shell fraction | `A_x` fraction | `k_y` fraction | `k_z` fraction | `curl_y` fraction | `curl_z` fraction | family counts |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 12 | 1.000 | 1.000 | 0.465 | 0.535 | 0.535 | 0.465 | `{'kz-Ax-curly': 3, 'ky-Ax-curlz': 1}` |
| 16 | 1.000 | 1.000 | 0.508 | 0.492 | 0.492 | 0.508 | `{'ky-Ax-curlz': 2, 'kz-Ax-curly': 2}` |
| 20 | 1.000 | 1.000 | 0.457 | 0.543 | 0.543 | 0.457 | `{'kz-Ax-curly': 2, 'ky-Ax-curlz': 2}` |

## Family sequences

- `n = 12`: `['kz-Ax-curly', 'kz-Ax-curly', 'ky-Ax-curlz', 'kz-Ax-curly']`
- `n = 16`: `['ky-Ax-curlz', 'kz-Ax-curly', 'ky-Ax-curlz', 'kz-Ax-curly']`
- `n = 20`: `['kz-Ax-curly', 'ky-Ax-curlz', 'kz-Ax-curly', 'ky-Ax-curlz']`

## Direct result

- observation: the already-resolved holonomy window can be labeled at the family level by combining first-shell momentum support with coarse field-axis and curl-axis readout, rather than demanding a rigid per-mode identity at exact clean symmetry
- conclusion: the holonomy-resolved clean window carries the same bounded first-shell family across n=[12, 16, 20]: the first-shell support stays at least 1.000, the coarse field stays fully x-polarized at 1.000, and the window-level ky/kz split only ranges from 0.457-0.508 versus 0.492-0.543; the raw mode ordering still shuffles=True, but the family label set is stable=True, so the clean-baseline identity problem closes in the bounded sense as a symmetry-degenerate first-shell family rather than a missing branch

## Artifacts

- results: `data/20260419_221932_transverse_clean_holonomy_symmetry_labels.json`
- plots: `plots/20260419_221932_transverse_clean_holonomy_symmetry_labels_fractions.png`, `plots/20260419_221932_transverse_clean_holonomy_symmetry_labels_scaled_eigenvalues.png`
- timestamp: `20260419_221932`
