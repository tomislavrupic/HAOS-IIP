# HBP Status Snapshot

This page gives a compact view of the current HBP branches and their terminal labels.
It is a snapshot, not a claim upgrade.

## PB-01

- Bridge ID: `PB-01-NETWORK-RECOVERY`
- Current status: `PREDICTION_NOT_DISTINCT_FROM_BASELINES`
- Terminal labels:
  - `CONTROL_INVALID`
  - `HOLDOUT_TRANSFER_PASS`
  - `MIXED_OPEN`
  - `PHYSICAL_MECHANISM_NOT_ESTABLISHED`
  - `PREDICTION_NOT_DISTINCT_FROM_BASELINES`
- Reading: PB-01 remains the synthetic network-recovery calibration benchmark; the target does not separate cleanly from baselines.

## PB-02

- Bridge ID: `PB-02-EXTERNAL-NETWORK-RECOVERY`
- Current status: `PRECOMMITMENT_DRAFT`
- Terminal labels:
  - `OPERATIONAL_MAPPING_PARTIAL`
  - `OPERATIONAL_MAPPING_VALID`
  - `PHYSICAL_MECHANISM_NOT_ESTABLISHED`
- Reading: PB-02 is still a frozen precommitment draft for external PowerGraph holdout work; it has not yet been run as a scored holdout bridge.

## PB-03

- Bridge ID: `PB-03-TEMPORAL-RECOVERY`
- Current status: `PREDICTION_NOT_DISTINCT_FROM_BASELINES`
- Terminal labels:
  - `DIMENSIONAL_ANALYSIS_FAIL`
  - `FORMAL_CORRESPONDENCE_ONLY`
  - `OPERATIONAL_MAPPING_PARTIAL`
  - `PHYSICAL_MECHANISM_NOT_ESTABLISHED`
- Reading: PB-03 remains a temporal recovery boundary case; the target does not beat the baseline.

## PB-04

- Bridge ID: `PB-04-CROSS-DOMAIN-TRANSFER`
- Current status: `PREDICTION_NOT_DISTINCT_FROM_BASELINES`
- Terminal labels:
  - `CROSS_DOMAIN_TRANSFER_BOUNDARY_OPEN`
  - `DIMENSIONAL_ANALYSIS_FAIL`
  - `FORMAL_CORRESPONDENCE_ONLY`
  - `OPERATIONAL_MAPPING_PARTIAL`
  - `PHYSICAL_MECHANISM_NOT_ESTABLISHED`
- Reading: PB-04 remains a cross-domain transfer boundary case; it does not upgrade the bridge.

## One-line view

| Bridge | Status | Terminal reading |
| --- | --- | --- |
| PB-01 | `PREDICTION_NOT_DISTINCT_FROM_BASELINES` | control-invalid synthetic calibration |
| PB-02 | `PRECOMMITMENT_DRAFT` | frozen draft for external PowerGraph holdout |
| PB-03 | `PREDICTION_NOT_DISTINCT_FROM_BASELINES` | temporal recovery boundary open |
| PB-04 | `PREDICTION_NOT_DISTINCT_FROM_BASELINES` | cross-domain transfer boundary open |

## Boundary

This snapshot does not change the registry, contracts, or holdout rules.
It simply records the current terminal labels in one place.
