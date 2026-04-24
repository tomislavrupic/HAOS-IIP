# Gene Network Demo Validation

## System
The system is a deterministic toy gene regulatory network with 12 genes: G0, G1, G2, G3, G4, G5, G6, G7, G8, G9, G10, G11.

G0 is the hub regulator. The network includes a negative feedback loop where G1 activates G2 and G2 inhibits G1 and G0. It includes a positive feedback motif where G3 and G4 mutually activate. The fragile downstream branch is G8, G9, G10, G11; its support depends on weakened inputs from G0 to G8, G4 to G8, and G5 to G9.

## Perturbation
The primary v0.1 experiment uses edge weakening. The perturbation parameter p is swept from 0.0 to 1.0 in 31 deterministic samples. At p = 0.0 no fragile-branch support edges are weakened. At p = 1.0 those selected support edges are fully removed.

Hub knockdown is implemented in `gene_network_model.py` as an available perturbation mode, but it is not the primary validation path for this artifact.

## Metrics
These are lightweight external proxy diagnostics aligned with HAOS-IIP stability logic. They are not frozen HAOS core metrics.

- recoverability: `1 - normalized_distance(perturbed_final_state, baseline_final_state)`, clipped to [0, 1]. The final state is represented by the mean activity over the final 30 time steps.
- delta_persistence: change in recoverability between consecutive perturbation levels.
- k_star: first perturbation index where recoverability stays below 0.65 for the current sample plus 2 following samples.
- safety_margin: distance in perturbation parameter from current p to p at k_star. At and beyond k_star this value is zero or negative.
- confidence: magnitude of the strongest negative delta_persistence signal.

The effective graph uses `effective_weight(source, target) = abs(W[source, target]) * activity_factor`, where `activity_factor` is source-gene mean activity over the final trajectory window.

## Result
Pass status: PASS

- baseline_stable: PASS
- recoverability_declines_gradually: PASS
- early_detection: PASS
- deterministic_repeated_runs: PASS
- fragile_branch_degrades_earlier: PASS

- baseline_stable: True
- k_star_level: 0.7
- first_visible_failure_level: 0.9333333333333333
- early_detection: True
- fragile_branch_k_star_level: 0.4
- minimum_fragile_branch_recoverability: 0.311028
- minimum_whole_network_recoverability: 0.602222

The collapse origin for this toy run is the fragile downstream branch, because its branch-level recoverability crosses the sustained-collapse threshold before the whole-network recoverability does.

## Limitations
This is a toy model.

It is not real biological validation.

It is not a full systems biology model.

It does not modify or prove HAOS-IIP core.
