# Transverse Clean Holonomy Continuation

## Purpose

Execute `CB4` in bounded form by using tiny cycle holonomy as a controlled symmetry splitter on the clean torus while keeping the same active coexact subspace metric used in `CB1-lite`.

This scan preserves periodic topology and does not replace the clean branch with a defect background. It only changes the boundary-cycle phase.

## Setup

- sizes: `[16, 20]`
- holonomy phases: `[0.0, 0.05]`
- restricted modes analyzed: `8`
- window size: `4`
- scan depth: `8`
- probe lattice side: `6`

## Strength summary

| holonomy phase | resolved pairs | worst fixed affinity | worst best affinity | worst best cosine | worst best drift | largest affinity gain | CB4 status |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0.000 | 0/1 | 0.330 | 0.342 | 0.542 | 0.008 | 0.011 | OPEN |
| 0.050 | 1/1 | 0.743 | 0.961 | 0.980 | 0.008 | 0.218 | RESOLVED |

## Direct result

- observation: the clean-baseline identity question can be sharpened without changing topology by introducing a tiny cycle holonomy and measuring the same low-window active coexact subspace under refinement
- conclusion: the clean family resolves under cycle holonomy: the zero-holonomy case starts at worst best affinity 0.342, and all adjacent refinement pairs first resolve at holonomy phase 0.050; this supports the reading that exact clean symmetry is what obscures branch identity

## Artifacts

- results: `data/20260419_220322_transverse_clean_holonomy_continuation.json`
- plots: `plots/20260419_220322_transverse_clean_holonomy_continuation_affinity.png`, `plots/20260419_220322_transverse_clean_holonomy_continuation_resolution.png`
- timestamp: `20260419_220322`
