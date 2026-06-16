# PB-01 Control Failure Localization

Status: frozen failure map

Current PB-01 verdict ceiling

- status: `PREDICTION_NOT_DISTINCT_FROM_BASELINES`
- labels: `CONTROL_INVALID`, `HOLDOUT_TRANSFER_PASS`, `MIXED_OPEN`, `PHYSICAL_MECHANISM_NOT_ESTABLISHED`

What the controls say

- `shuffled_target_labels`: PASS
- `intentional_leakage_positive_control`: `TARGET_LEAKAGE_DETECTED`
- `seed_repeat`: PASS
- `topology_destroyed_graph`: `CONTROL_INVALID`
- `degree_preserving_rewire`: `CONTROL_INVALID`
- `weight_shuffled_graph`: `CONTROL_INVALID`
- `parameter_matched_null`: `CONTROL_INVALID`
- `perturbation_free_baseline`: PASS

Failure pattern

The leakage detector works, but the destructive controls do not move the
recoverability signal enough to satisfy the frozen contract.

Observed control weakness

- topology-destroyed graph:
  - recovery shift: `-0.010769`
  - degree delta: `2.537352`
  - density delta: `0.198068`
  - spectral gap delta: `0.212957`
  - shortest-path delta: `2.736111`
- degree-preserving rewire:
  - recovery shift: `-0.000063`
  - degree delta: `0.000000`
  - density delta: `0.000000`
  - spectral gap delta: `0.574012`
  - shortest-path delta: `0.273148`
- weight-shuffled graph:
  - recovery shift: `-0.002548`
  - degree delta: `0.316328`
  - density delta: `0.000000`
  - spectral gap delta: `0.149992`
  - shortest-path delta: `0.000000`
- parameter-matched null:
  - recovery shift: `-0.013710`
  - degree delta: `1.863243`
  - density delta: `0.141304`
  - spectral gap delta: `0.197116`
  - shortest-path delta: `1.023148`

Declared control threshold

- destructive controls are expected to move by at least `0.02` in absolute
  recovery-quality shift

Interpretation

1. PB-01 is still usable as a synthetic predictive calibration benchmark.
2. PB-01 is not yet control-valid.
3. The topology-destroyed control now breaks degree, density, spectral gap, and
   shortest-path structure, but the recoverability signal still does not move
   enough.
4. The degree-preserving rewire keeps degree and density fixed, so it is
   specifically checking whether the signal is merely degree-driven. It is not
   yet moving the score.
5. The weight-shuffled control is changing weighting structure more than
   topology, but not enough to break the recovery ranking.
6. The parameter-matched null still leaves too much of the ranking intact.
7. The leakage positive control is already doing its job, so the failure is not
   a blind detector.

Repair direction

- isolate whether HAOS is duplicating degree/spectral structure or whether the
  synthetic dynamics are simply too recoverable
- keep the leakage control unchanged so the guard remains honest
- freeze PB-01 as control-invalid unless a future control contract can move the
  recovery signal by the declared threshold without changing the target score

Boundary

This note localizes the control problem only. It does not upgrade PB-01 to an
external validation claim.
