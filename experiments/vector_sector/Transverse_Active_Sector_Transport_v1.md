# Transverse Active Sector Transport

## Purpose

Execute `N1` of the continuum-closure program by making the `T8` comparison maps `J_n` explicit on the restricted coexact sector, then use those maps for a first low-window branch-identity test.

## Transport map

For each refinement level, the active restricted edge modes are transported into one common probe space:

- the unit torus is partitioned into a fixed `6^3` probe lattice
- each edge amplitude is assigned to the probe cell containing its midpoint
- amplitudes are separated into `x`, `y`, and `z` channels
- each channel is normalized by the square root of the number of contributing edges

This defines the explicit comparison map `J_n` used in the present test.

## Setup

- sizes: `[12, 16, 20, 24, 28]`
- variants: `['baseline', 'mild_disorder', 'line_defect']`
- restricted modes: `6`
- transported low-window modes: `4`
- probe lattice side: `6`
- mild-disorder strength: `0.12`

## Transport metrics

| branch | refinement | mean overlap | min overlap | mean principal cosine | max scaled-eigen drift | mode matching |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| baseline | 12 -> 16 | 0.351 | 0.307 | 0.510 | 0.017 | 0->1, 1->3, 2->2, 3->0 |
| baseline | 16 -> 20 | 0.468 | 0.306 | 0.544 | 0.008 | 0->1, 1->0, 2->2, 3->3 |
| baseline | 20 -> 24 | 0.343 | 0.276 | 0.441 | 0.004 | 0->2, 1->3, 2->0, 3->1 |
| baseline | 24 -> 28 | 0.414 | 0.262 | 0.532 | 0.003 | 0->1, 1->2, 2->0, 3->3 |
| line_defect | 12 -> 16 | 0.626 | 0.223 | 0.665 | 0.036 | 0->1, 1->0, 2->3, 3->2 |
| line_defect | 16 -> 20 | 0.659 | 0.283 | 0.762 | 0.017 | 0->1, 1->0, 2->3, 3->2 |
| line_defect | 20 -> 24 | 0.687 | 0.418 | 0.816 | 0.009 | 0->0, 1->1, 2->3, 3->2 |
| line_defect | 24 -> 28 | 0.605 | 0.284 | 0.746 | 0.006 | 0->1, 1->0, 2->2, 3->3 |
| mild_disorder | 12 -> 16 | 0.970 | 0.967 | 0.971 | 0.018 | 0->0, 1->1, 2->2, 3->3 |
| mild_disorder | 16 -> 20 | 0.981 | 0.979 | 0.981 | 0.008 | 0->0, 1->1, 2->2, 3->3 |
| mild_disorder | 20 -> 24 | 0.989 | 0.988 | 0.989 | 0.005 | 0->0, 1->1, 2->2, 3->3 |
| mild_disorder | 24 -> 28 | 0.992 | 0.991 | 0.992 | 0.003 | 0->0, 1->1, 2->2, 3->3 |

## Direct result

- observation: explicit probe-space comparison maps now carry the restricted coexact branch across refinement, so active-sector identity can be measured directly by transported overlaps, principal-angle alignment, and scaled-eigenvalue drift rather than inferred only from separate spectra
- conclusion: within the tested low transported window, the strongest refinement pair is mild_disorder 24->28 with mean overlap 0.992, the weakest pair is baseline 20->24 with mean overlap 0.343, and the smallest max scaled-eigen drift is baseline 24->28 at 0.003; this gives `T8` a real `J_n` mechanism and an initial branch-identity test, but the continuum question still depends on whether these transported overlaps stay strong as the window and background family are enlarged

## Artifacts

- results: `data/20260416_232531_transverse_active_sector_transport.json`
- plots: `plots/20260416_232531_transverse_active_sector_transport_overlap.png`, `plots/20260416_232531_transverse_active_sector_transport_stability.png`
- timestamp: `20260416_232531`
