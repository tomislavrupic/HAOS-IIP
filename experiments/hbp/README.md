# HAOS-IIP Physical Bridge Program

Lifecycle status: PB-01 through PB-04 are `QUARANTINED_INVALID` and may not be
repaired in place. The versioned `HBP-IR-01` synthetic integrity calibration is
complete and retained as a supporting instrument; it does not promote any
historical result. See the
[branch lifecycle registry](../../docs/branch_governance/branch_lifecycle_summary.md).

HBP is an experimental bridge registry and verifier. It classifies candidate
bridges from HAOS-IIP telemetry to externally measurable or real-world systems.

It does not claim HAOS-IIP is fundamental physics. It does not reopen Bell,
hydrogen, semiconductor, CST, or frozen phase branches.

## Run

```bash
uv run python experiments/hbp/run_hbp_registry.py
uv run python experiments/hbp/check_hbp_bundle.py
```

PB-01 can be run directly:

```bash
uv run python experiments/hbp/pb01_network_recovery/run_pb01_network_recovery.py
```

PB-02 is frozen as a precommitment draft:

```bash
cat experiments/hbp/pb02_external_network_recovery/README.md
cat experiments/hbp/pb02_external_network_recovery/precommitment_contract.json
uv run python experiments/hbp/pb02_external_network_recovery/run_pb02_external_network_recovery.py
```

PB-03 and PB-04 are the next contract-first branches:

```bash
cat experiments/hbp/pb03_temporal_recovery/README.md
cat experiments/hbp/pb03_temporal_recovery/precommitment_contract.json
cat experiments/hbp/pb03_temporal_recovery/dataset_selection.md
cat experiments/hbp/pb03_temporal_recovery/source_manifest.json
cat experiments/hbp/pb03_temporal_recovery/execution_readiness.md
cat experiments/hbp/pb03_temporal_recovery/split_manifest.json
cat experiments/hbp/pb03_temporal_recovery/execution_contract.json
cat experiments/hbp/pb03_temporal_recovery/metrics_manifest.json
cat experiments/hbp/pb03_temporal_recovery/baselines_manifest.json
cat experiments/hbp/pb03_temporal_recovery/data_schema.json
cat experiments/hbp/pb03_temporal_recovery/control_manifest.json
cat experiments/hbp/pb04_cross_domain_transfer/README.md
cat experiments/hbp/pb04_cross_domain_transfer/precommitment_contract.json
cat experiments/hbp/pb04_cross_domain_transfer/dataset_selection.md
cat experiments/hbp/pb04_cross_domain_transfer/source_manifest.json
cat experiments/hbp/pb04_cross_domain_transfer/execution_readiness.md
cat experiments/hbp/pb04_cross_domain_transfer/split_manifest.json
cat experiments/hbp/pb04_cross_domain_transfer/execution_contract.json
cat experiments/hbp/pb04_cross_domain_transfer/metrics_manifest.json
cat experiments/hbp/pb04_cross_domain_transfer/baselines_manifest.json
cat experiments/hbp/pb04_cross_domain_transfer/data_schema.json
cat experiments/hbp/pb04_cross_domain_transfer/control_manifest.json
```

## Outputs

- [HBP Status Snapshot](hbp_status_snapshot.md)
- `results/hbp_bridge_registry.json`
- `results/hbp_bridge_registry.csv`
- `results/hbp_result.json`
- `results/hbp_report.md`
- `pb01_network_recovery/results/precommitment_contract.json`
- `pb01_network_recovery/results/pb01_result.json`
- `pb02_external_network_recovery/README.md`
- `pb02_external_network_recovery/precommitment_contract.json`
- `pb02_external_network_recovery/run_pb02_external_network_recovery.py`
- `pb03_temporal_recovery/README.md`
- `pb03_temporal_recovery/precommitment_contract.json`
- `pb03_temporal_recovery/dataset_selection.md`
- `pb03_temporal_recovery/source_manifest.json`
- `pb03_temporal_recovery/execution_readiness.md`
- `pb03_temporal_recovery/split_manifest.json`
- `pb03_temporal_recovery/execution_contract.json`
- `pb03_temporal_recovery/metrics_manifest.json`
- `pb03_temporal_recovery/baselines_manifest.json`
- `pb03_temporal_recovery/data_schema.json`
- `pb03_temporal_recovery/control_manifest.json`
- `pb04_cross_domain_transfer/README.md`
- `pb04_cross_domain_transfer/precommitment_contract.json`
- `pb04_cross_domain_transfer/dataset_selection.md`
- `pb04_cross_domain_transfer/source_manifest.json`
- `pb04_cross_domain_transfer/execution_readiness.md`
- `pb04_cross_domain_transfer/split_manifest.json`
- `pb04_cross_domain_transfer/execution_contract.json`
- `pb04_cross_domain_transfer/metrics_manifest.json`
- `pb04_cross_domain_transfer/baselines_manifest.json`
- `pb04_cross_domain_transfer/data_schema.json`
- `pb04_cross_domain_transfer/control_manifest.json`

## Bridge Levels

- `FORMAL_CORRESPONDENCE`
- `ANALOGY_BRIDGE`
- `OPERATIONAL_MAPPING`
- `PREDICTIVE_BRIDGE`
- `EMPIRICALLY_SUPPORTED_BRIDGE`
- `PHYSICAL_MECHANISM_CANDIDATE`

No level may be skipped. Missing required bridge fields, weak mapping status,
missing holdout, missing external data, missing controls, or missing replication
automatically lower the effective classification.

## PB-01 Boundary

PB-01 is a synthetic network-recovery calibration benchmark. It asks whether a
precommitted HAOS-style calibrated score predicts post-perturbation recovery
ranking better than graph, spectral, and domain-diffusion baselines.

A clean `PREDICTION_NOT_DISTINCT_FROM_BASELINES` result is valid. The benchmark
is not empirical physical validation and not a physical mechanism claim.

## PB-02 Boundary

PB-02 is a quarantined historical artifact. Its scored output remains preserved,
but control routing, holdout selection, and feature reconstruction are not valid
enough for interpretation. Any replacement must use the new `HBP-IR-01`
precommitment and a new result namespace.

## PB-03 Boundary

PB-03 is a quarantined temporal-recovery artifact. Its present controls and
baseline coverage do not satisfy their contracts.

## PB-04 Boundary

PB-04 is a quarantined cross-domain transfer artifact. Its baseline and HAOS
candidate currently alias the same model path, so no transfer interpretation is
authorized.
