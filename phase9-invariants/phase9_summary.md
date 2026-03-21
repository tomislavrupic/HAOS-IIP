# Phase IX Summary

Timestamp: `2026-03-21T15:55:56Z`.

## Objective

Test whether the coefficient-like structures extracted from the frozen Phase VIII short-time trace regime stabilize across refinement, whether normalized combinations behave invariant-like, and whether this behavior differs from a deterministic generic-graph control hierarchy.

## Frozen References

- Phase VI operator manifest: `phase6-operator/phase6_operator_manifest.json`.
- Phase VII spectral manifest: `phase7-spectral/phase7_spectral_manifest.json`.
- Phase VIII trace manifest: `phase8-trace/phase8_trace_manifest.json`.
- Phase VIII contract note: `PHASE_VIII_TRACE_CONTRACT_NOTE.md`.

## Frozen Trace Contract

- Operational short-time window: `[0.02, 0.03, 0.05, 0.075, 0.1]`.
- Trace method continuity: exact deterministic spectrum with the frozen `full_exact`, `lambda_le_7_5`, and `lambda_le_7_8` truncation policies.
- Dominant Phase VII low-mode exponent used for rescaled descriptors: `1.972182228379`.

## Frozen Branch Coefficient Stabilization

| coefficient | values across refinement | relative span | best extrapolation | limit estimate | fit R^2 |
|---|---|---:|---|---:|---:|
| b0 | [3.992541196696, 3.992270826129, 3.992219728715, 3.992201767437] | 8.5020802e-05 | linear_in_h^2 | 3.992179654625 | 0.999988793788 |
| b1 | [-15.196408636509, -15.38315827383, -15.417962259342, -15.43015999025] | 0.015221236998 | linear_in_h^2 | -15.445632893201 | 0.999998915236 |
| b2 | [27.978374055231, 28.598645942437, 28.714851153144, 28.755622542545] | 0.027260519706 | linear_in_h^2 | 28.806767863436 | 0.999996380389 |

- Successive coefficient-vector distances: `[0.647775203445, 0.121305280239, 0.04255691666]`.
- Distances shrink across refinement: `True`.
- Largest raw coefficient relative span: `2.007728209099`.

## Ratio Stabilization

| ratio | values across refinement | relative span | best extrapolation | limit estimate |
|---|---|---:|---|---:|
| b1_over_b0 | [-3.806199582633, -3.853235149567, -3.862002421471, -3.865075186357] | 0.015305769735 | linear_in_h^2 | -3.868971692687 |
| b2_over_b0 | [7.007660704512, 7.163503476583, 7.192703083602, 7.202948201941] | 0.027344664673 | linear_in_h^2 | 7.215797137512 |

- Maximum branch ratio span: `0.027344664673`.
- Ratio spans are smaller than raw coefficient spans: `True`.
- Maximum subset-induced ratio drift: `0.069683572506`.
- Maximum truncation-induced ratio drift: `0.07202504059`.

## Rescaled Invariant-Like Descriptors

| descriptor | values across refinement | relative span | best extrapolation | limit estimate |
|---|---|---:|---|---:|
| I0 | [4.278286099098, 4.361284281729, 4.410697810555, 4.446116798825] | 0.038369228803 | linear_in_h^1 | 4.48488158993 |
| I1 | [-16.284010764769, -16.805053891612, -17.034125624739, -17.184575714784] | 0.053519229865 | linear_in_h^1 | -17.435256281568 |
| I2 | [29.980777379304, 31.242075114536, 31.724839742812, 32.025149001718] | 0.065434108792 | linear_in_h^1 | 32.63855294095 |

- Maximum branch rescaled-invariant span: `0.065434108792`.
- Maximum subset-induced invariant drift: `0.070571240424`.
- Maximum truncation-induced invariant drift: `0.110877321888`.

## Generic-Graph Control Comparison

- Control hierarchy: `generic_open_grid_scalar_block_control`.
- Control description: Deterministic control hierarchy with the same total node counts as the frozen branch, constructed as four repeated scalar open-grid Laplacian blocks with open boundary connectivity.

| control coefficient | relative span | best extrapolation | limit estimate |
|---|---:|---|---:|
| b0 | 0.000302231342 | linear_in_h^1 | 3.992101755014 |
| b1 | 0.078116509612 | linear_in_h^1 | -15.517159116109 |
| b2 | 0.121073848083 | linear_in_h^1 | 29.012971564412 |

- Branch max ratio span: `0.027344664673`.
- Control max ratio span: `0.121364602609`.
- Ratio-span multiplier (control / branch): `4.438328429342`.
- Final coefficient-vector distance multiplier (control / branch): `8.109704279502`.
- Branch best extrapolation models: `{'b0': 'linear_in_h^2', 'b1': 'linear_in_h^2', 'b2': 'linear_in_h^2'}`.
- Control best extrapolation models: `{'b0': 'linear_in_h^1', 'b1': 'linear_in_h^1', 'b2': 'linear_in_h^1'}`.

## Bounded Interpretation

- The frozen branch exhibits reproducible coefficient stabilization in the Phase VIII window, with the strongest plateau in `b0` and shrinking successive drift across all three scaled coefficients.
- Dimensionless ratios are materially more stable than the raw trace coefficients, and the rescaled descriptors form compact refinement bands on the frozen branch.
- The deterministic open-grid scalar control does not behave identically to the frozen branch: its ratio spans are larger, its successive coefficient-vector distances remain larger, and its simple extrapolation model family differs from the branch.
- This is a feasibility result only. It does not assign geometric, continuum, or universal meaning to any coefficient or descriptor.

## Non-Claims

- No continuum limit is established.
- No geometric coefficient interpretation is made.
- No universal invariant is claimed.
- No continuum or physical correspondence is asserted.

Phase IX establishes coefficient stabilization and invariant-tracking feasibility for the frozen operator hierarchy.
