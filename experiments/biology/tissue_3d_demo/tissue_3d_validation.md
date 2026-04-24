# 3D Tissue Demo Validation

## System
The system is a deterministic toy 3D tissue-like lattice with grid size (8, 8, 8), for 512 cells.

Each node is a cell with a scalar propagation state in [0, 1]. Edge classes are diagonal_neighbor, local_neighbor, long_range_signal. The local_neighbor class is 6-connected nearest-neighbor support, diagonal_neighbor provides weak local support, and long_range_signal provides sparse deterministic coordination links.

The primary perturbation is a local lesion centered at (3.5, 3.5, 3.5) with radius 4.1. Perturbation level p weakens edge support in that spherical region. A secondary signaling_weakening mode is implemented for long_range_signal edges but is not the primary validation sweep.

## Metrics
These are lightweight external proxy diagnostics aligned with HAOS-IIP stability logic. They are not frozen HAOS core metrics.

- recoverability: `0.15 * largest_component_fraction + 0.45 * weighted_degree_retention + 0.25 * spatial_coherence_retention + 0.15 * propagation_retention`, clipped to [0, 1].
- delta_persistence: change in recoverability between consecutive perturbation levels.
- k_star: first perturbation index where recoverability remains below 0.70 for the current sample plus 2 following samples.
- safety_margin: distance in perturbation parameter from current p to p at k_star. At and beyond k_star this value is zero or negative.
- confidence: magnitude of the strongest negative delta_persistence observed so far.
- visible_failure: True when `largest_component_fraction < 0.85` or `spatial_coherence_retention < 0.50`.

## Result
Pass status: PASS

- baseline_stable: PASS
- recoverability_declines_gradually: PASS
- early_detection: PASS
- deterministic_repeated_runs: PASS
- robustness_early_detection: PASS

- baseline_stable: True
- k_star perturbation level: 0.8
- first visible failure level: 1.0
- early_detection: True
- robustness variation: local_neighbor_strength_x_1_02
- robustness early_detection: True
- minimum recoverability: 0.563310

## Limitations
This is a toy 3D model.

It is not real tissue simulation.

It is not morphogenesis.

It makes no biological claims.

It does not modify or prove HAOS-IIP core.
