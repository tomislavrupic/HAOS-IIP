# Post-67.1 Bridge Audit

Status: continuity note, not a new phase and not a new mechanism.

This note ties the hardened 66.5 scale-bridge baseline to the 67.1 HBP
recovery companion release. It exists to keep the validation language aligned
across the paper spine and the experiment registry.

## Frozen Baselines

- 66.5 scale-bridge hardening is the continuity baseline.
- 67.1 HBP recovery benchmarks inherit that boundary discipline.
- quick reproduction remains the public smoke test for frozen artifacts.

## Bridge Rows

### PB-01 Network Recovery

- status: `PREDICTION_NOT_DISTINCT_FROM_BASELINES`
- labels: `CONTROL_INVALID`, `HOLDOUT_TRANSFER_PASS`, `MIXED_OPEN`, `PHYSICAL_MECHANISM_NOT_ESTABLISHED`, `PREDICTION_NOT_DISTINCT_FROM_BASELINES`
- reading: the benchmark remains a valid synthetic calibration case, but the current control contract is still not strong enough to support a predictive bridge claim.

### PB-03 Temporal Recovery Under Damage

- status: `PREDICTION_NOT_DISTINCT_FROM_BASELINES`
- labels: `CONTROL_INVALID`, `MIXED_OPEN`, `PREDICTION_NOT_DISTINCT_FROM_BASELINES`
- reading: the frozen holdout still does not distinguish HAOS-style descriptors from the best conventional baselines.

### PB-04 Cross-Domain Structure Transfer

- status: `PREDICTION_NOT_DISTINCT_FROM_BASELINES`
- labels: `HAOS_BELL_STATUS_OPEN`, `MIXED_OPEN`, `PREDICTION_NOT_DISTINCT_FROM_BASELINES`
- reading: the transfer sidecar remains open at the boundary, but the frozen target still does not outperform the reference baselines.

## Continuity Claim

The 66.5 hardening pass clarified where the scale-bridge line stops. The 67.1
HBP companion keeps that discipline: the registry is live, the benchmarks are
frozen, and the negative results stay negative.

