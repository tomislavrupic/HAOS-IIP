# Transverse Clean Pinning Continuation

## Purpose

Execute `CB2` in bounded form by testing whether infinitesimal symmetry pinning resolves the clean-baseline identity question when measured with the same windowed subspace metric introduced in `CB1-lite`.

This scan keeps the same active coexact sector, the same probe-space transport rule, and the same low-window scan. Only the smooth pinning strength changes.

## Setup

- sizes: `[12, 16]`
- pinning strengths: `[0.0, 0.015]`
- restricted modes analyzed: `8`
- window size: `4`
- scan depth: `8`
- probe lattice side: `6`

## Strength summary

| pinning strength | resolved pairs | worst fixed affinity | worst best affinity | worst best cosine | worst best drift | largest affinity gain | CB2 status |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0.000 | 0/1 | 0.333 | 0.366 | 0.530 | 0.017 | 0.032 | OPEN |
| 0.015 | 1/1 | 0.942 | 0.956 | 0.978 | 0.017 | 0.014 | RESOLVED |

## Direct result

- observation: the clean-baseline identity question can be sharpened by introducing an infinitesimal symmetry pinning and measuring the same low-window subspace metric as the pinning is removed
- conclusion: the clean family resolves under infinitesimal pinning: the baseline starts at worst best affinity 0.366, and all adjacent refinement pairs first resolve at pinning strength 0.015; this strongly supports a degenerate clean-family interpretation rather than branch failure

## Artifacts

- results: `data/20260419_215411_transverse_clean_pinning_continuation.json`
- plots: `plots/20260419_215411_transverse_clean_pinning_continuation_affinity.png`, `plots/20260419_215411_transverse_clean_pinning_continuation_resolution.png`
- timestamp: `20260419_215411`
