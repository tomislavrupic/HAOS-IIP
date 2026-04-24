# Microtubule Lattice Demo Validation

## System
The system is a deterministic coarse-grained microtubule-inspired interaction lattice with 13 protofilaments and 24 dimers per protofilament, for 312 total nodes.

Each node represents one tubulin dimer position. Edge classes are lateral, longitudinal, seam_or_diagonal, weak_support. Longitudinal edges are the strongest local supports, lateral edges provide medium circumferential support, seam_or_diagonal edges provide weak deterministic geometric offsets, and weak_support edges provide very weak longer-range support.

The primary perturbation mode is lateral weakening: lateral edge strengths are multiplied by `1 - p` for p in [0, 1]. A deterministic defect patch mode is implemented in `microtubule_lattice_model.py` around protofilaments 4-6 and z indices 10-14, but it is not the primary validation sweep.

## Metrics
These are lightweight external proxy diagnostics aligned with HAOS-IIP stability logic. They are not frozen HAOS core metrics.

- recoverability: `0.05 * largest_component_fraction + 0.90 * weighted_degree_retention + 0.05 * propagation_retention`, clipped to [0, 1].
- delta_persistence: change in recoverability between consecutive perturbation levels.
- k_star: first perturbation index where recoverability remains below 0.70 for the current sample plus 2 following samples.
- safety_margin: distance in perturbation parameter from current p to p at k_star. At and beyond k_star this value is zero or negative.
- confidence: magnitude of the strongest negative delta_persistence observed so far.
- visible_failure: True when `largest_component_fraction < 0.85` or `propagation_retention < 0.50`.

The edge score is `interaction_strength * locality_stability * perturbation_factor`, with `locality_stability = exp(-distance / lambda_locality)` and cylindrical distance used only as a coarse locality measure.

The filtration order sorts edge supports by this score, so strongest and most local supports enter before weaker or less local supports.

## Result
Pass status: PASS

- baseline_stable: PASS
- recoverability_declines_gradually: PASS
- early_detection: PASS
- deterministic_repeated_runs: PASS
- robustness_early_detection: PASS

- baseline stable: True
- k_star perturbation level: 0.8333333333333334
- first visible failure level: 0.8666666666666667
- early_detection: True
- robustness variation: longitudinal_strength_x_1_02
- robustness early_detection: True
- minimum recoverability: 0.578094

## Limitations
This is a coarse-grained toy model.

It is not molecular dynamics.

It is not real biological validation.

It makes no consciousness claims.

It makes no Orch-OR claims.

It makes no quantum computation claims.

It does not modify or prove HAOS-IIP core.
