# PB-02 Holdout Report

Status: `CONTROL_INVALID`

This report executes the frozen PB-02 comparison on the local PowerGraph dataset.
It does not upgrade the bridge claim ceiling.

## Dataset Validation

- family: `ieee24`
- sample count: `21500`
- dataset validation hash: `pb02_dataset_validation_e49fce0968836133f7148809`

## Split Manifest

- development: `12938`
- calibration: `4241`
- holdout: `4321`

## Best Baseline

- model: `closeness_centrality`
- holdout spearman: `0.991454`

## Baseline + HAOS
- model: `baseline_plus_haos`
- holdout spearman: `0.679913`

## Controls
- `shuffled_target_labels`: `FAIL` (1.000000)
- `topology_destroyed_graph`: `FAIL` (0.991454)
- `degree_preserving_rewire`: `FAIL` (0.991454)
- `weight_shuffled_graph`: `FAIL` (0.991454)
- `parameter_matched_null`: `PASS` (-0.006747)
- `perturbation_free_baseline`: `FAIL` (0.991454)
- `seed_repeat`: `PASS` (0.679913)
- `intentional_leakage_positive_control`: `TARGET_LEAKAGE_DETECTED` (1.000000)

## Uncertainty
- bootstrap replicates: `120`
- incremental prediction CI: `-0.000355` to `0.001704`

## Non-Claims
- operational mapping only until holdout transfer is demonstrated
- no physical mechanism claim
- no empirical bridge claim
- no universal recoverability claim
- no claim that HAOS explains power-grid physics

Result hash: `pb02_result_c84f151554dc98485a16bee4`
