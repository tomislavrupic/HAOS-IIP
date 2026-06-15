# HAOS-IIP Physical Bridge Program

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

## Outputs

- `results/hbp_bridge_registry.json`
- `results/hbp_bridge_registry.csv`
- `results/hbp_result.json`
- `results/hbp_report.md`
- `pb01_network_recovery/results/precommitment_contract.json`
- `pb01_network_recovery/results/pb01_result.json`

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
