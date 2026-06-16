# PB-02 External Network Recovery

Status: precommitment draft

Purpose

Test whether frozen HAOS-IIP features add predictive value on an independently
defined real-world network-recovery system.

Frozen dataset manifest

- `dataset_manifest.json`

Frozen split manifests

- `development_split.json`
- `calibration_split.json`
- `holdout_split.json`

Bundle check

```bash
uv run python experiments/hbp/pb02_external_network_recovery/check_pb02_bundle.py
```

Execution contract

- `execution_contract.json`

Primary candidate

- dataset: PowerGraph
- task family: cascading failure analysis on power grids
- public source: <https://arxiv.org/abs/2402.02827>
- dataset landing page: <https://figshare.com/articles/dataset/PowerGraph/22820534>
- code: <https://github.com/PowerGraph-Datasets>

Why this target

- direct semantics: power-grid recovery and cascading failure behavior
- public benchmark framing
- not designed around HAOS-IIP
- supports baseline-vs-HAOS comparisons on a real-world system

Precommitted comparison

- best conventional baseline
- versus best conventional baseline + frozen HAOS features

Frozen rule

- no new features after holdout inspection
- no post hoc threshold movement
- no feature selection using holdout outcomes

Claim ceiling

- operational mapping only until holdout transfer is demonstrated
- no physical mechanism claim
- no empirical bridge claim
- no claim of universal recoverability

Recommended benchmark shape

- frozen development split: tune baselines and frozen HAOS feature usage
- frozen calibration split: choose one HAOS feature bundle and stop
- untouched holdout split: score once, report baseline-only, HAOS-only, baseline+HAOS, ablated HAOS, matched null

Next repair dependency

PB-01 control validity should be improved before PB-02 is executed so that the frozen HAOS features are not carried forward through a known control leak.
