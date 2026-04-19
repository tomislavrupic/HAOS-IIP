# Transverse Active Branch Purity

## Purpose

Execute a first bounded `V3` read by testing whether the active restricted branch remains physically coexact under refinement rather than drifting back into harmonic or mixed sectors.

This first freeze uses the restricted active window itself:

- sizes: `[12, 16, 20]`
- variants: `['baseline', 'line_defect', 'mild_disorder']`
- restricted modes: `6`
- purity window: first `4` restricted modes

## V3 criteria

For this first bounded `V3` freeze, a size slice counts as `PASS` only if all of the following hold on the tested low restricted window:

- worst exact fraction <= `1.0e-08`
- worst harmonic fraction <= `1.0e-08`
- worst coexact fraction >= `0.999999`

This is an internal active-branch purity test. It does not yet test transported purity against a separately transported harmonic subspace.

## Case results

| branch | size | worst exact fraction | worst harmonic fraction | worst coexact fraction | status |
| --- | ---: | ---: | ---: | ---: | --- |
| baseline | 12 | 1.478e-32 | 1.376e-31 | 1.000000000 | PASS |
| line_defect | 12 | 1.409e-32 | 2.689e-31 | 1.000000000 | PASS |
| mild_disorder | 12 | 7.222e-33 | 1.801e-31 | 1.000000000 | PASS |
| baseline | 16 | 4.430e-33 | 1.056e-31 | 1.000000000 | PASS |
| line_defect | 16 | 1.228e-32 | 8.718e-33 | 1.000000000 | PASS |
| mild_disorder | 16 | 4.705e-33 | 1.037e-31 | 1.000000000 | PASS |
| baseline | 20 | 5.886e-33 | 1.267e-31 | 1.000000000 | PASS |
| line_defect | 20 | 1.086e-32 | 5.366e-32 | 1.000000000 | PASS |
| mild_disorder | 20 | 1.543e-32 | 5.527e-32 | 1.000000000 | PASS |

## Variant summary

| branch | passing sizes | worst exact fraction | worst harmonic fraction | worst coexact fraction | V3 status |
| --- | ---: | ---: | ---: | ---: | --- |
| baseline | 3/3 | 1.478e-32 | 1.376e-31 | 1.000000000 | PASS |
| line_defect | 3/3 | 1.409e-32 | 2.689e-31 | 1.000000000 | PASS |
| mild_disorder | 3/3 | 1.543e-32 | 1.801e-31 | 1.000000000 | PASS |

## Direct result

- observation: the restricted active branch now carries explicit exact / harmonic / coexact fractions, so harmonic recontamination can be tested directly on the low active window instead of being assumed away
- conclusion: all tested variants pass the first bounded V3 purity thresholds, and the worst harmonic leakage remains 2.689e-31; in this first bounded read, the low restricted active window remains physically coexact to leading order under refinement

## Artifacts

- results: `data/20260419_164039_transverse_active_branch_purity.json`
- plots: `plots/20260419_164039_transverse_active_branch_purity.png`
- timestamp: `20260419_164039`
