# Transverse Clean Subspace Transport

## Purpose

Fast-track the clean-baseline identity question by treating the low active coexact branch as a potentially rotating low-dimensional subspace rather than as a fixed mode-by-mode labeling problem.

This bounded `CB1-lite` scan keeps the same probe-space transport map, the same active coexact sector, and the same refinement pairs. The only added freedom is that the low transported window may slide inside the first `8` restricted coexact modes.

## Setup

- sizes: `[12, 16, 20]`
- variants: `['baseline', 'mild_disorder']`
- restricted modes analyzed: `8`
- window size: `4`
- scan depth: `8`
- probe lattice side: `6`

## CB1-lite criteria

- best projector affinity >= `0.80`
- best mean principal cosine >= `0.80`
- best scaled-window drift <= `0.05`
- best minus fixed affinity gain >= `0.20`

If the fixed window already meets them, the branch is treated as already stable. If only the shifted best window meets them with a substantial affinity gain, the pair is treated as rescued by windowed subspace transport. If neither happens, the clean-baseline ambiguity remains open.

## Pair results

| branch | refinement | fixed affinity | best affinity | best mean principal cosine | best scaled-window drift | best window shift | status |
| --- | --- | ---: | ---: | ---: | ---: | --- | --- |
| baseline | 12 -> 16 | 0.311 | 0.320 | 0.527 | 0.017 | 2->2 | open |
| baseline | 16 -> 20 | 0.316 | 0.323 | 0.499 | 0.008 | 2->2 | open |
| mild_disorder | 12 -> 16 | 0.942 | 0.957 | 0.978 | 0.018 | 3->3 | stable_fixed |
| mild_disorder | 16 -> 20 | 0.963 | 0.963 | 0.981 | 0.008 | 0->0 | stable_fixed |

## Variant summary

| branch | resolved pairs | worst fixed affinity | worst best affinity | largest affinity gain | worst best mean principal cosine | worst best drift | CB1-lite status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| baseline | 0/2 | 0.311 | 0.320 | 0.009 | 0.499 | 0.017 | OPEN |
| mild_disorder | 2/2 | 0.942 | 0.957 | 0.015 | 0.978 | 0.018 | RESOLVED |

## Direct result

- observation: the clean-baseline identity question can be fast-tracked by scanning low contiguous coexact windows as transported subspaces; this keeps the same active-sector transport map while allowing the clean basis and the low window to rotate or slide inside a bounded scan depth
- conclusion: the fast-track CB1-lite scan improves low-window identification but does not yet close the clean baseline: the weakest best-window case is baseline 12->16 with best projector affinity 0.320; this means the next shortest move is still CB2-style infinitesimal symmetry pinning, but we now know whether simple window sliding already resolves the issue

## Artifacts

- results: `data/20260419_213842_transverse_clean_subspace_transport.json`
- plots: `plots/20260419_213842_transverse_clean_subspace_transport_affinity.png`, `plots/20260419_213842_transverse_clean_subspace_transport_shift.png`
- timestamp: `20260419_213842`
