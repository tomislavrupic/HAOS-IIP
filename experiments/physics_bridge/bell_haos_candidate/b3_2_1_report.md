# B3.2.1 Null Failure Localization

Implemented fact: B3.2 remains frozen and CHSH is not computed.
Design choice: this harness compares whole G distributions, relation matrices, topology dependence, and operator variants.
Heuristic: matrix similarity and rank-order similarity localize whether nulls preserve geometry provenance.
Unverified hypothesis: HAOS-specific Bell geometry remains unestablished.

## Labels
- NULL_FAILURE_LOCALIZED
- GENERIC_BILINEAR_GEOMETRY_SUSPECTED
- TARGET_SPECIFIC_GEOMETRY_NOT_ESTABLISHED
- CHSH_SCORING_NOT_AUTHORIZED
- HAOS_BELL_DERIVATION_NOT_ESTABLISHED
- MIXED / OPEN
- TOPOLOGY_ZERO_CONTROL_REJECTED
- REFINEMENT_SPECIFICITY_NOT_ESTABLISHED
- SETTING_SEMANTICS_NOT_ESTABLISHED

## Target Distribution
- mean G: `0.0468316377276`
- variance G: `0.0631102135971`
- positive/negative count: `19` / `15`
- pair-order signature: `holdout:h3_h1>holdout:h0_h2>chsh_held_back:a1_b0>chsh_held_back:a1_b1>chsh_held_back:a0_b1>holdout:h3_h2>chsh_held_back:a0_b0>holdout:h0_h1`

## Failing Nulls
- label_permuted_settings
- refinement_broken_control

## Operator Dependence
- target_J: status `REFERENCE_OPERATOR`, pair-order rho `1`, aligned distance `0`
- degree_matched_J: status `OPERATOR_GEOMETRY_PERSISTS`, pair-order rho `0.833333333333`, aligned distance `0.577342599933`
- spectrum_matched_J: status `OPERATOR_GEOMETRY_DEGRADED`, pair-order rho `-0.309523809524`, aligned distance `0.989605900683`
- locality_preserving_randomized_J: status `OPERATOR_GEOMETRY_DEGRADED`, pair-order rho `0.047619047619`, aligned distance `0.925875340412`
- orientation_reversed_J: status `OPERATOR_GEOMETRY_PERSISTS`, pair-order rho `-1`, aligned distance `0`

## Boundary
- This is a failure-localization harness.
- No Bell/CHSH score is computed.
- No B3.2 field, threshold, setting, vector, or operator is changed.
- HAOS Bell derivation is not established.
