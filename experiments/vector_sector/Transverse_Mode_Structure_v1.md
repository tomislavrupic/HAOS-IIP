# Transverse Mode Structure

## Purpose

Inspect the lowest restricted transverse modes directly to check whether they appear as extended circulation fields.

## Setup

- branch: `baseline`
- `n = 16`
- `epsilon = 0.2`
- restricted modes analyzed: `10`

## Lowest 10 restricted transverse modes

| mode | eigenvalue | `||d0* a||` | `||d1 a||` | IPR |
| --- | ---: | ---: | ---: | ---: |
| 1 | 0.150761 | 1.382e-16 | 0.388280 | 0.000181 |
| 2 | 0.150761 | 1.167e-16 | 0.388280 | 0.000227 |
| 3 | 0.150761 | 1.813e-16 | 0.388280 | 0.000195 |
| 4 | 0.150761 | 1.429e-16 | 0.388280 | 0.000235 |
| 5 | 0.150761 | 1.245e-16 | 0.388280 | 0.000182 |
| 6 | 0.301523 | 1.664e-16 | 0.549111 | 0.000206 |
| 7 | 0.301523 | 2.435e-16 | 0.549111 | 0.000197 |
| 8 | 0.301523 | 1.439e-16 | 0.549111 | 0.000238 |
| 9 | 0.301523 | 1.557e-16 | 0.549111 | 0.000223 |
| 10 | 0.452284 | 1.622e-16 | 0.672521 | 0.000195 |

## Direct result

- observation: the lowest restricted transverse modes remain numerically divergence-free while their spatial support extends over the periodic lattice rather than concentrating at a single location
- conclusion: the inspected low restricted modes behave like extended circulation fields in the tested baseline branch

## Artifacts

- results: `data/20260310_145035_transverse_mode_structure.json`
- plots: `plots/20260310_145035_transverse_vector_field_profiles.png`, `plots/20260310_145035_divergence_curl_phase_mode_structure.png`
- timestamp: `20260310_145035`
