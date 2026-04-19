# Transverse Active Branch Effective Closure V6b

## Purpose

Execute `V6b` by rerunning the bounded `V6` family with a divergence-aware coarse-graining / projection step, without changing the validated branch, the response family, or the law target.

## Setup

- sizes: `[12, 16, 20]`
- variant family: `mild_disorder`
- disorder strengths: `[0.015, 0.06, 0.12]`
- penalties: `[8.0, 10.0, 12.0]`
- effective mode: `0`

## V6b criteria

- apply a centered-difference transverse projection after coarse-graining, then judge the same local-law window `curl curl A ≈ lambda A` with residual <= `0.08`
- require projected divergence norm <= `1.0e-08`
- keep the same response criterion as `V6`: largest-size adjacent scaled-eigen response jump <= `0.12`

## Combo results

| penalty | disorder | worst raw divergence norm | worst projected divergence norm | worst projected residual | max adjacent response jump | V6b status |
| ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 8.0 | 0.015 | 1.005e-03 | 1.115e-16 | 0.050 | 0.006 | PASS |
| 8.0 | 0.060 | 4.026e-03 | 1.190e-16 | 0.049 | 0.006 | PASS |
| 8.0 | 0.120 | 8.115e-03 | 1.238e-16 | 0.055 | 0.006 | PASS |
| 10.0 | 0.015 | 1.005e-03 | 1.150e-16 | 0.050 | 0.006 | PASS |
| 10.0 | 0.060 | 4.026e-03 | 1.178e-16 | 0.049 | 0.006 | PASS |
| 10.0 | 0.120 | 8.114e-03 | 1.194e-16 | 0.055 | 0.006 | PASS |
| 12.0 | 0.015 | 1.005e-03 | 1.333e-16 | 0.050 | 0.006 | PASS |
| 12.0 | 0.060 | 4.026e-03 | 1.221e-16 | 0.049 | 0.006 | PASS |
| 12.0 | 0.120 | 8.114e-03 | 1.168e-16 | 0.055 | 0.006 | PASS |

## Direct result

- observation: once the coarse-grained active field is projected with the same centered-difference transverse symbol used by the reconstruction, the divergence defect collapses while the bounded effective-law response remains smooth across the validated mild-disorder family
- conclusion: all 9 tested combos pass `V6b`: the raw divergence defect reaches 8.115e-03, but the projected divergence defect drops to 1.333e-16, the worst projected local-law residual remains 0.055, and the largest adjacent scaled-response jump remains 0.006; this means the original V6 failure was a coarse-graining artifact rather than a breakdown of bounded effective-law closure on the validated active branch

## Artifacts

- results: `data/20260419_210029_transverse_active_branch_effective_closure_v6b.json`
- plots: `plots/20260419_210029_transverse_active_branch_effective_closure_v6b_divergence.png`, `plots/20260419_210029_transverse_active_branch_effective_closure_v6b_response.png`
- timestamp: `20260419_210029`
