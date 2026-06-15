# CST Toy Benchmark Report

Status: FAIL

Implemented fact: this report is generated from frozen HAOS-IIP ledgers only.
Design choice: branch identity is an operational hash of the compact address summary.
Heuristic: non-front event times are rank-interpolated from Phase 18 recovery-time fields.
Analogy: generalized CST selectivity is an engineering signal/noise ratio, not physical resonance Q.
Unverified hypothesis: repeated recovery of this address might indicate closure stability in this toy slice.

## Reasons
- negative control family periodic_diagonal_augmented_control performs equally or better than target mean
- negative control family randomized_edge_control performs equally or better than target mean
- negative control family shuffled_structure_control performs equally or better than target mean
- available seed count is below configured minimum for PASS

## Warnings
- one or more requested controls are unavailable in the frozen artifact slice

## Recoverability Vector
- recovery_rate: 1
- branch_identity_persistence: 1
- closure_fidelity: 0.790262870038
- control_margin: 0.0766236067677
- variance_bound: 0.803229236982
- ablation_balance: 1
- optional_scalar: 0.7783526189646687

## Distance Summary
- target_scalar_mean: 0.209737129962
- target_scalar_std: 0.0157416610414
- control_scalar_mean: 0.286360736729
- control_scalar_std: 0.1487835779

## Unavailable Controls
- perturbation_free_control: not present in frozen Phase 15-18 artifact slice

## Output Files
- benchmark_report: experiments/cst/runs/benchmark_report.md
- benchmark_result: experiments/cst/runs/benchmark_result.json
- closure_distance_components: experiments/cst/runs/closure_distance_components.csv
- closure_distance_matrix: experiments/cst/runs/closure_distance_matrix.csv
- closure_signatures: experiments/cst/runs/closure_signatures.json
- control_distributions: experiments/cst/runs/control_distributions.csv
- cst_runs: experiments/cst/runs/cst_runs.json
- recoverability_vectors: experiments/cst/runs/recoverability_vectors.csv
- seed_manifest: experiments/cst/runs/seed_manifest.json

## Limitations
- This benchmark does not test Bell correlations, spectra, semiconductors, biology, consciousness, or cosmology.
- A PASS/OPEN/FAIL verdict is scoped only to this frozen toy slice and configured thresholds.
- Scalar compression is optional and heuristic; component metrics and the vector are the auditable result.
