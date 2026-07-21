# HBP Status Snapshot

This page gives a compact view of the current HBP branches and their terminal labels.
It is a snapshot, not a claim upgrade.

Lifecycle decision: PB-01 through PB-04 are `QUARANTINED_INVALID`. Their
artifacts remain preserved, but none is an active development target. The only
eligible instrument successor was the separately versioned `HBP-IR-01`; its
minimal synthetic integrity run now passes and is retained without predictive
promotion. See the
[branch lifecycle registry](../../docs/branch_governance/branch_lifecycle_summary.md).

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
- Contract status: `PRECOMMITMENT_DRAFT`
- Stored result status: `CONTROL_INVALID`
- Audit status: quarantined pending integrity repair
- Terminal labels:
  - `CONTROL_INVALID`
  - `PREDICTION_NOT_DISTINCT_FROM_BASELINES`
  - `MIXED_OPEN`
- Reading: scored PowerGraph artifacts exist, but the baseline reconstruction and shuffled-label control are internally inconsistent. Hash reproducibility does not validate this result.

## PB-03

- Bridge ID: `PB-03-TEMPORAL-RECOVERY`
- Current status: `PREDICTION_NOT_DISTINCT_FROM_BASELINES`
- Terminal labels:
  - `DIMENSIONAL_ANALYSIS_FAIL`
  - `FORMAL_CORRESPONDENCE_ONLY`
  - `OPERATIONAL_MAPPING_PARTIAL`
  - `PHYSICAL_MECHANISM_NOT_ESTABLISHED`
- Reading: PB-03 remains a temporal recovery boundary case. Its stored result is reproducible, but the runner does not implement the full frozen baseline manifest and its seed-repeat control is tautological, so no predictive promotion is supported.

## PB-04

- Bridge ID: `PB-04-CROSS-DOMAIN-TRANSFER`
- Current status: `PREDICTION_NOT_DISTINCT_FROM_BASELINES`
- Terminal labels:
  - `CROSS_DOMAIN_TRANSFER_BOUNDARY_OPEN`
  - `DIMENSIONAL_ANALYSIS_FAIL`
  - `FORMAL_CORRESPONDENCE_ONLY`
  - `OPERATIONAL_MAPPING_PARTIAL`
  - `PHYSICAL_MECHANISM_NOT_ESTABLISHED`
- Reading: PB-04 remains a cross-domain transfer boundary case. Its baseline and HAOS models currently use the same feature path, and its seed-repeat and topology controls do not execute their declared mechanisms.

## One-line view

| Bridge | Status | Terminal reading |
| --- | --- | --- |
| PB-01 | `PREDICTION_NOT_DISTINCT_FROM_BASELINES` | control-invalid synthetic calibration |
| PB-02 | `CONTROL_INVALID` | scored artifact quarantined; contract remains draft |
| PB-03 | `PREDICTION_NOT_DISTINCT_FROM_BASELINES` | reproducible result, incomplete baseline/control implementation |
| PB-04 | `PREDICTION_NOT_DISTINCT_FROM_BASELINES` | non-discriminative candidate path and invalid control mechanisms |

## Boundary

This snapshot does not change the registry, contracts, or holdout rules.
It simply records the current terminal labels in one place.
