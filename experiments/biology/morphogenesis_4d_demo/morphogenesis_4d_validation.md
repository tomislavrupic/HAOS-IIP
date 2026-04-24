# 4D Morphogenesis Demo Validation

## System
The system is a deterministic toy 4D developmental structure: a 3D grid with size (8, 8, 8), 512 cells, and 24 developmental time steps.

Each cell has a scalar state in [0, 1]. The target morphology is deterministic: early development is a compact central form, mid development introduces an axis polarity gradient, and late development transitions toward a core-shell/anterior-posterior pattern.

Edge classes are developmental_signal, diagonal_neighbor, local_neighbor. The local_neighbor class gives 6-connected nearest-neighbor support, diagonal_neighbor gives weak local support, and developmental_signal gives deterministic longer coordination links.

The primary perturbation is developmental drift. It delays and weakens target pull in a deterministic region centered at (5.1, 3.5, 3.5) during mid-to-late development, while weakening developmental_signal links in the same affected region. A local developmental lesion mode is implemented but is not the primary validation sweep.

## Metrics
These are lightweight external proxy diagnostics aligned with HAOS-IIP stability logic. They are not frozen HAOS core metrics.

- recoverability: `0.35 * final_morphology_match + 0.32 * trajectory_coherence + 0.25 * interaction_support_retention + 0.08 * spatial_continuity_retention`, clipped to [0, 1].
- delta_persistence: change in recoverability between consecutive perturbation levels.
- k_star: first perturbation index where recoverability remains below 0.70 for the current sample plus 2 following samples.
- safety_margin: distance in perturbation parameter from current p to p at k_star. At and beyond k_star this value is zero or negative.
- confidence: magnitude of the strongest negative delta_persistence observed so far.
- visible_failure: True when `final_morphology_match < 0.50` or `spatial_continuity_retention < 0.50`.

## Result
Pass status: PASS

- baseline_stable: PASS
- recoverability_declines_gradually: PASS
- early_detection: PASS
- deterministic_repeated_runs: PASS
- robustness_early_detection: PASS

- baseline stable: True
- k_star perturbation level: 0.8666666666666667
- first visible failure level: 0.9666666666666667
- early_detection: True
- robustness variation: target_weight_x_1_02
- robustness early_detection: True
- minimum recoverability: 0.654913

## Limitations
This is a toy 4D developmental model.

It is not real morphogenesis.

It is not a biological validation claim.

It makes no consciousness claims.

It makes no quantum biology claims.

It does not modify or prove HAOS-IIP core.
