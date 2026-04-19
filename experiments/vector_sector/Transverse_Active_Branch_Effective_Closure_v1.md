# Transverse Active Branch Effective Closure

## Purpose

Execute a first bounded `V6` read by asking whether the validated active branch closes onto one local coarse-grained law family and supports a smooth bounded response read on the same family.

This first freeze stays on the `V5`-validated `mild_disorder` branch and does not reopen the clean baseline.

## Setup

- sizes: `[12, 16, 20]`
- variant family: `mild_disorder`
- disorder strengths: `[0.015, 0.06, 0.12]`
- penalties: `[8.0, 10.0, 12.0]`
- effective mode: `0`

## V6 criteria

- every tested combo must keep the coarse-grained active mode inside the local law window `curl curl A ≈ lambda A` with relative residual <= `0.08`
- every tested combo must keep relative divergence norm <= `1.0e-08`
- at the largest tested size, adjacent mild-disorder steps at fixed penalty must keep the scaled-eigen response jump <= `0.12`

## Combo results

| penalty | disorder | worst curl-curl residual | worst divergence norm | largest-size scaled eigenvalue | max adjacent response jump | V6 status |
| ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 8.0 | 0.015 | 0.050 | 1.005e-03 | 38.864 | 0.006 | OPEN |
| 8.0 | 0.060 | 0.050 | 4.026e-03 | 38.711 | 0.006 | OPEN |
| 8.0 | 0.120 | 0.056 | 8.115e-03 | 38.482 | 0.006 | OPEN |
| 10.0 | 0.015 | 0.050 | 1.005e-03 | 38.864 | 0.006 | OPEN |
| 10.0 | 0.060 | 0.050 | 4.026e-03 | 38.711 | 0.006 | OPEN |
| 10.0 | 0.120 | 0.056 | 8.114e-03 | 38.482 | 0.006 | OPEN |
| 12.0 | 0.015 | 0.050 | 1.005e-03 | 38.864 | 0.006 | OPEN |
| 12.0 | 0.060 | 0.050 | 4.026e-03 | 38.711 | 0.006 | OPEN |
| 12.0 | 0.120 | 0.056 | 8.114e-03 | 38.482 | 0.006 | OPEN |

## Direct result

- observation: the validated mild-disorder active branch admits the same coarse-grained local law test across the full bounded V5 family, and the corresponding scaled law parameter moves smoothly as the disorder strength changes
- conclusion: the first bounded V6 closure remains open: failing combos = ['(penalty=8.0, disorder=0.015)', '(penalty=8.0, disorder=0.060)', '(penalty=8.0, disorder=0.120)', '(penalty=10.0, disorder=0.015)', '(penalty=10.0, disorder=0.060)', '(penalty=10.0, disorder=0.120)', '(penalty=12.0, disorder=0.015)', '(penalty=12.0, disorder=0.060)', '(penalty=12.0, disorder=0.120)']; max residual = 0.056, max divergence = 8.115e-03, max adjacent response jump = 0.006

## Artifacts

- results: `data/20260419_201944_transverse_active_branch_effective_closure.json`
- plots: `plots/20260419_201944_transverse_active_branch_effective_closure_residuals.png`, `plots/20260419_201944_transverse_active_branch_effective_closure_response.png`
- timestamp: `20260419_201944`
