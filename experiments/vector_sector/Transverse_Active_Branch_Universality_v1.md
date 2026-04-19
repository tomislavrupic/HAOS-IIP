# Transverse Active Branch Universality

## Purpose

Execute a first bounded `V5` read by asking whether the current active-branch contract survives a small operator-family and substrate-family variation without changing the contract itself.

This first freeze stays on the already-validated `mild_disorder` branch. It does not reopen the clean baseline or enlarge the physical-sector definition.

## Setup

- sizes: `[12, 16, 20]`
- variant family: `mild_disorder`
- disorder strengths: `[0.015, 0.06, 0.12]`
- penalties: `[8.0, 10.0, 12.0]`
- restricted modes: `6`
- transported active window: `4`
- fixed scaling exponent: `2.00`

## V5 criteria

- keep the same active-window contract fixed: `V1` identity thresholds, `V3` purity thresholds, and `V4` scaling exponent `alpha = 2.00`
- test bounded operator-family variation by changing penalty across `[8.0, 10.0, 12.0]`
- test bounded substrate-family variation by changing mild-disorder strength across `[0.015, 0.06, 0.12]`
- a combo counts as `PASS` only if every adjacent refinement pair passes `V1`, every tested size passes `V3`, and the fixed-`alpha` scaling drift stays <= `0.05`

## Combo results

| penalty | disorder | identity passes | worst mean overlap | worst scaling drift | worst harmonic fraction | worst coexact fraction | V5 status |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 8.0 | 0.015 | 2/2 | 0.971 | 0.018 | 1.238e-31 | 1.000000000 | PASS |
| 8.0 | 0.060 | 2/2 | 0.970 | 0.018 | 1.818e-31 | 1.000000000 | PASS |
| 8.0 | 0.120 | 2/2 | 0.970 | 0.018 | 8.162e-32 | 1.000000000 | PASS |
| 10.0 | 0.015 | 2/2 | 0.971 | 0.018 | 9.415e-32 | 1.000000000 | PASS |
| 10.0 | 0.060 | 2/2 | 0.970 | 0.018 | 1.530e-31 | 1.000000000 | PASS |
| 10.0 | 0.120 | 2/2 | 0.970 | 0.018 | 2.104e-31 | 1.000000000 | PASS |
| 12.0 | 0.015 | 2/2 | 0.971 | 0.018 | 4.800e-32 | 1.000000000 | PASS |
| 12.0 | 0.060 | 2/2 | 0.970 | 0.018 | 1.364e-31 | 1.000000000 | PASS |
| 12.0 | 0.120 | 2/2 | 0.970 | 0.018 | 6.067e-32 | 1.000000000 | PASS |

## Direct result

- observation: the first bounded V5 scan keeps the active-branch contract fixed and asks whether it survives a small family of penalty choices and mild-disorder substrates instead of being re-earned after each parameter change
- conclusion: all 9 tested operator/substrate combos pass the first bounded V5 universality check: the weakest combo is penalty=8.0, disorder=0.120 with worst mean overlap 0.970, while the largest fixed-alpha drift over the whole family is 0.018; this means the current active-branch contract survives a real bounded mild-disorder and penalty family without changing the physical-sector rule, purity rule, or scaling exponent

## Artifacts

- results: `data/20260419_172934_transverse_active_branch_universality.json`
- plots: `plots/20260419_172934_transverse_active_branch_universality.png`
- timestamp: `20260419_172934`
