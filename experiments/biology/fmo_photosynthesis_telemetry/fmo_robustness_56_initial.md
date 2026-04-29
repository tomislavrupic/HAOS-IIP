# Phase 56 FMO Photosynthesis Telemetry Robustness

Toy FMO-like 7-site weighted-network telemetry. This is not a molecular simulation or a quantum-biology claim.

## Configuration

- seed_count: 10
- thermal_noise: 0.00,0.04,0.08,0.12,0.16
- disorder_scale: 0.18
- damage_scale: 0.42
- null_level: 5
- permutation_trials: 32
- null_candidates: 8

## Overall Summary

- runs: 50
- pass_rate: 0.040000
- recoverability_score: 0.803724 +/- 0.039498
- site_identity_retention: 0.984934 +/- 0.005867
- pathway_identity_retention: 0.546597 +/- 0.094994
- active_null_z: 2.472024 +/- 0.574870
- control_pass_max: 2

## By Noise Level

| thermal_noise | runs | pass_rate | recoverability | site_identity | pathway_identity | active_null_z | max_control_passes |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.000 | 10 | 0.000 | 0.810337 +/- 0.035822 | 0.985998 +/- 0.005710 | 0.562496 +/- 0.087049 | 2.543880 +/- 0.528820 | 2 |
| 0.040 | 10 | 0.100 | 0.802751 +/- 0.039829 | 0.984781 +/- 0.005772 | 0.544286 +/- 0.095509 | 2.462474 +/- 0.586447 | 2 |
| 0.080 | 10 | 0.100 | 0.802331 +/- 0.040002 | 0.984713 +/- 0.005829 | 0.543256 +/- 0.095936 | 2.462263 +/- 0.593233 | 2 |
| 0.120 | 10 | 0.000 | 0.801861 +/- 0.040276 | 0.984635 +/- 0.005903 | 0.542107 +/- 0.096672 | 2.460967 +/- 0.600445 | 2 |
| 0.160 | 10 | 0.000 | 0.801341 +/- 0.040652 | 0.984546 +/- 0.005993 | 0.540842 +/- 0.097716 | 2.430537 +/- 0.556018 | 2 |
