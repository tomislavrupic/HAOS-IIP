# Phase VIII Trace Summary

Timestamp: `2026-03-21T15:16:37Z`.

## Objective

Test whether the frozen branch-local cochain-Laplacian hierarchy supports a reproducible short-time trace scaling regime across refinement, without making any continuum or geometric claim.

## Frozen References

- Phase V authority manifest: `phase5-readout/phase5_authoritative_manifest.json`.
- Phase VI operator manifest: `phase6-operator/phase6_operator_manifest.json`.
- Phase VII spectral manifest: `phase7-spectral/phase7_spectral_manifest.json`.

## Trace Proxy Method

- Trace proxy: `T_h(t) = Tr(exp(-t * delta_h))`.
- Evaluation method: exact torus-mode spectrum induced by the frozen periodic DK2D cochain-Laplacian block structure.
- Master probe grid: `[0.01, 0.015, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0]`.
- Randomness used: `False`.
- Direct dense validation on the coarsest level: max eigenvalue diff `0.0`, max trace diff `0.0`.

## Operational Short-Time Window

- Operational short-time window: `[0.02, 0.03, 0.05, 0.075, 0.1]`.
- Max slope dispersion in window: `0.005211400534`.
- Max normalized trace relative range in window: `0.005820971051`.
- Max kernel share in window: `0.010092438116`.
- Max relative trace error for `lambda <= 7.5`: `0.039036855531`.
- Max relative trace error for `lambda <= 7.8`: `0.014977531175`.

The operational window is selected automatically as the longest contiguous probe band that satisfies the slope-dispersion, normalized-trace-spacing, kernel-share, truncation-error, and coefficient-fit robustness thresholds.

## Trace Curve and Slope Behaviour

- t = `0.02`: slope dispersion `0.001260704553`, normalized trace rel. range `0.001264880157`, kernel share `0.007506673711`.
- t = `0.03`: slope dispersion `0.001872181624`, normalized trace rel. range `0.001878373242`, kernel share `0.007800114643`.
- t = `0.05`: slope dispersion `0.002972961777`, normalized trace rel. range `0.003067515519`, kernel share `0.008412117856`.
- t = `0.075`: slope dispersion `0.004161996207`, normalized trace rel. range `0.004483237187`, kernel share `0.009225092725`.
- t = `0.1`: slope dispersion `0.005211400534`, normalized trace rel. range `0.005820971051`, kernel share `0.010092438116`.
- The full trace is monotone decreasing in `t` at every refinement level.

## Coefficient-Like Surrogate

- Surrogate fit: `T_h(t) ≈ h^(-2) * (b0(h) + b1(h) * t + b2(h) * t^2)` on the operational short-time window.

| level | n_side | h | b0(h) | b1(h) | b2(h) | R^2 |
|---|---:|---:|---:|---:|---:|---:|
| R1 | 12 | 0.083333333333 | 3.992541196696 | -15.196408636509 | 27.978374055231 | 0.999996635281 |
| R2 | 24 | 0.041666666667 | 3.992270826129 | -15.383158273830 | 28.598645942437 | 0.999996461490 |
| R3 | 36 | 0.027777777778 | 3.992219728715 | -15.417962259342 | 28.714851153143 | 0.999996428345 |
| R4 | 48 | 0.020833333333 | 3.992201767437 | -15.430159990250 | 28.755622542545 | 0.999996416672 |

- Relative coefficient spans across refinement: `{'b0': 8.5020802e-05, 'b1': 0.015221236998, 'b2': 0.027260519706}`.
- Coefficient trends: `{'b0': 'nonincreasing', 'b1': 'nonincreasing', 'b2': 'nondecreasing'}`.
- Minimum fit R^2: `0.999996416672`.
- Max relative coefficient span: `0.027260519706`.
- Max probe-subgrid coefficient drift in the selected window: `0.019376305038`.

## Truncation Robustness

- `lambda_le_7_5`: max coefficient rel. diff. vs full `0.110877321888`, max slope diff. `0.010940479377`.
- `lambda_le_7_8`: max coefficient rel. diff. vs full `0.044047599547`, max slope diff. `0.004232595036`.

## Kernel Sensitivity

- Max kernel share in the selected window: `0.010092438116`.
- Max log-slope difference after kernel subtraction in the selected window: `0.003044349031`.

## Bounded Interpretation

- The frozen hierarchy supports a reproducible operational short-time window on the declared probe grid.
- Within that window, the log-trace slope is stable across refinement, the size-factored coefficient-like terms drift smoothly, and deterministic spectral cutoffs do not overturn the observed regime.
- This is a feasibility result only. It does not identify geometric coefficients, universal asymptotic coefficients, or continuum structure.

## Non-Claims

- No continuum correspondence is established.
- No geometric or universal coefficient interpretation is made.
- No Seeley-DeWitt claim is made.

Phase VIII establishes short-time trace asymptotic feasibility for the frozen operator hierarchy.
