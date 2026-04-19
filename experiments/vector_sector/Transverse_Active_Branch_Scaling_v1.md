# Transverse Active Branch Scaling

## Purpose

Execute a first bounded `V4` read by asking whether one coherent scaling family `a_n = n^alpha` keeps the same tracked active window finite and comparable across refinement.

This note reuses the same transported active window as `V1` and `V2`. It does not retune the window by branch or by pair. Formal `V4` closure is evaluated only on adjacent pairs that already pass `V1` identity.

## Setup

- sizes: `[12, 16, 20]`
- variants: `['baseline', 'line_defect', 'mild_disorder']`
- restricted modes: `6`
- transported active window: `4`
- candidate exponents: `[1.5, 1.75, 2.0, 2.25, 2.5]`
- chosen global exponent: `2.00`

## V4 criteria

- evaluate formal `V4` pass only on adjacent pairs that already pass `V1` identity
- global active-window scaling family must satisfy max scaled-eigen drift <= `0.05` on those `V1`-passing pairs
- pairwise best exponent may differ from the chosen global exponent by at most `0.25`
- at least `2.0` identity-passing pairs are required for a bounded closure call

## Pair results

| branch | refinement | V1 identity | pair-best alpha | pair-best mean drift | global-alpha mean drift | global-alpha max drift | |pair-best - global| | V4 status |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| baseline | 12 -> 16 | OPEN | 2.00 | 0.017 | 0.017 | 0.017 | 0.00 | OPEN |
| baseline | 16 -> 20 | OPEN | 2.00 | 0.008 | 0.008 | 0.008 | 0.00 | OPEN |
| line_defect | 12 -> 16 | OPEN | 2.00 | 0.034 | 0.034 | 0.036 | 0.00 | OPEN |
| line_defect | 16 -> 20 | OPEN | 2.00 | 0.015 | 0.015 | 0.017 | 0.00 | OPEN |
| mild_disorder | 12 -> 16 | PASS | 2.00 | 0.017 | 0.017 | 0.018 | 0.00 | PASS |
| mild_disorder | 16 -> 20 | PASS | 2.00 | 0.008 | 0.008 | 0.008 | 0.00 | PASS |

## Direct result

- observation: with the transported active window fixed, the same exponent scan can now be applied across every adjacent refinement pair, so V4 asks directly whether one scaling family stabilizes the active branch instead of rescuing it with pair-specific rescaling
- conclusion: the first bounded V4 read passes: every tested pair has the same pair-best exponent alpha=2.00, and on the V1-passing subset the global choice alpha=2.00 keeps the max scaled-eigen drift bounded at 0.018; this means the tracked active window already closes onto one coherent scaling family in the current bounded regime

## Artifacts

- results: `data/20260419_165739_transverse_active_branch_scaling.json`
- plots: `plots/20260419_165739_transverse_active_branch_scaling.png`
- timestamp: `20260419_165739`
