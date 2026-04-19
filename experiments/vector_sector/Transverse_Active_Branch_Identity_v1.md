# Transverse Active Branch Identity

## Purpose

Execute the first bounded `V1` validation read by asking whether the same low active coexact family survives as the same physical object across refinement.

This first freeze uses the executed bounded window:

- sizes: `[12, 16, 20]`
- variants: `['baseline', 'line_defect', 'mild_disorder']`
- restricted modes: `6`
- transported low-window modes: `4`

## Source artifact

- results: `data/20260419_095948_transverse_active_branch_identity.json`
- plots: `plots/20260419_095948_transverse_active_branch_identity_overlap.png`, `plots/20260419_095948_transverse_active_branch_identity_margin.png`

## V1 criteria

For this first bounded `V1` freeze, an adjacent refinement pair counts as `PASS` only if all of the following hold:

- mean matched overlap >= `0.75`
- minimum matched overlap >= `0.50`
- mean principal cosine >= `0.80`
- max scaled-eigen drift <= `0.05`
- mean assignment margin >= `0.10`

Assignment margin is the chosen-match overlap minus the strongest competing overlap on the same row / column pair, clipped below at zero. In this first `V1` read it is the bounded operational proxy for "no arbitrary relabeling."

## Pair results

| branch | refinement | mean overlap | min overlap | mean principal cosine | max scaled-eigen drift | mean assignment margin | exact-match fraction | status |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| baseline | 12 -> 16 | 0.482 | 0.192 | 0.552 | 0.017 | 0.222 | 0.250 | FAIL |
| baseline | 16 -> 20 | 0.328 | 0.197 | 0.495 | 0.008 | 0.000 | 0.500 | FAIL |
| line_defect | 12 -> 16 | 0.791 | 0.627 | 0.815 | 0.036 | 0.366 | 0.500 | PASS |
| line_defect | 16 -> 20 | 0.719 | 0.506 | 0.876 | 0.017 | 0.208 | 0.500 | FAIL |
| mild_disorder | 12 -> 16 | 0.970 | 0.967 | 0.971 | 0.018 | 0.948 | 1.000 | PASS |
| mild_disorder | 16 -> 20 | 0.981 | 0.979 | 0.981 | 0.008 | 0.970 | 1.000 | PASS |

## Variant summary

| branch | passing pairs | worst mean overlap | worst mean principal cosine | worst max drift | worst mean margin | mean exact-match fraction | V1 status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| baseline | 0 / 2 | 0.328 | 0.495 | 0.017 | 0.000 | 0.375 | OPEN |
| line_defect | 1 / 2 | 0.719 | 0.815 | 0.036 | 0.208 | 0.500 | OPEN |
| mild_disorder | 2 / 2 | 0.970 | 0.971 | 0.018 | 0.948 | 1.000 | PASS |

## Direct result

- observation: with the active coexact restriction and fixed probe-space transport map held fixed, branch identity can now be tested directly against explicit `V1` thresholds instead of being inferred from qualitative transport language
- conclusion: `V1` is not yet closed on the bounded trio; `mild_disorder` passes every adjacent refinement pair, `line_defect` passes one pair and remains open on the next, and the clean baseline fails both tested pairs. In plain language: active-branch identity is already robust under mild symmetry breaking, partly stabilized on the line-defect background, and still unstable on the clean-background side

## What this means for the next step

This first bounded `V1` read points directly to `V2`.

The present failure pattern is not:

```text
the active branch does not exist
```

It is:

```text
the active branch exists, but clean-background identity remains unstable,
while bounded symmetry breaking stabilizes it strongly
```

So the next clean move is:

- run `V2` explicitly as a disorder-threshold / symmetry-breaking validation step
- measure where the branch crosses from unstable relabeling to robust transported identity under the same fixed `J_n` map
