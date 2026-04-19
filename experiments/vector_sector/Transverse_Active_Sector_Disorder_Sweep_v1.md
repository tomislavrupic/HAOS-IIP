# Transverse Active Sector Disorder Sweep

## Purpose

Execute `V2` in bounded form by keeping the same active-sector transport map `J_n` fixed and asking whether the clean-torus branch-mixing weakness decays continuously once bounded smooth disorder is introduced.

## Setup

- sizes: `[12, 16, 20]`
- disorder strengths: `[0.0, 0.015, 0.03, 0.06, 0.09, 0.12]`
- restricted modes: `6`
- transported low-window modes: `4`
- probe lattice side: `6`

## Strength summary

| disorder strength | mean pair overlap | worst pair mean overlap | worst matched overlap | `16->20` mean overlap | `16->20` mean principal cosine | `16->20` exact-match fraction | `16->20` mean assignment margin | `16->20` max drift |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.000 | 0.327 | 0.317 | 0.222 | 0.337 | 0.518 | 0.500 | 0.017 | 0.008 |
| 0.015 | 0.976 | 0.971 | 0.967 | 0.981 | 0.981 | 1.000 | 0.972 | 0.008 |
| 0.030 | 0.976 | 0.971 | 0.967 | 0.981 | 0.981 | 1.000 | 0.972 | 0.008 |
| 0.060 | 0.976 | 0.970 | 0.967 | 0.981 | 0.981 | 1.000 | 0.971 | 0.008 |
| 0.090 | 0.976 | 0.970 | 0.967 | 0.981 | 0.981 | 1.000 | 0.971 | 0.008 |
| 0.120 | 0.976 | 0.970 | 0.967 | 0.981 | 0.981 | 1.000 | 0.970 | 0.008 |

## Direct result

- observation: with the active-sector comparison map held fixed, the disorder sweep isolates whether branch identity improves continuously under bounded symmetry breaking instead of only under topological defects
- conclusion: the largest refinement pair rises from mean overlap 0.337 at disorder strength 0.000 to 0.981 at strength 0.090, and it first crosses the present bounded identity threshold at strength 0.015; because the same `J_n` map is used throughout, this points to a clean-background degeneracy / branch-mixing problem rather than a broken transport construction

## Artifacts

- results: `data/20260419_101910_transverse_active_sector_disorder_sweep.json`
- plots: `plots/20260419_101910_transverse_active_sector_disorder_sweep_overlap.png`, `plots/20260419_101910_transverse_active_sector_disorder_sweep_drift.png`
- timestamp: `20260419_101910`
