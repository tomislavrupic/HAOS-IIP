# B3.2 Relational Geometry Audit

Implemented fact: this harness computes a frozen HAOS-IIP relational invariant before any CHSH score.
Design choice: CHSH scoring is authorization-gated and not computed here.
Heuristic: `J_lambda` is an oriented chain-transport pairing from Phase 18 chain signatures.
Unverified hypothesis: HAOS-IIP Bell derivation remains unestablished.

## Gate Results
- G1_RELATIONAL_SENSITIVITY: `RELATIONAL_STRUCTURE_DETECTED`
- G2_SIGN_STRUCTURE: `SIGN_CHANGING_GEOMETRY_DETECTED`
- G3_COVARIANCE: `COVARIANCE_PASS`
- G4_HOLDOUT_TRANSFER: `HOLDOUT_TRANSFER_PASS`
- G5_NULL_REJECTION: `NULL_REJECTION_FAIL`
- G6_CHSH_SCOREBOARD: `CHSH_SCORING_NOT_AUTHORIZED`

## Target Geometry Summary
- all-pair range: `1.3159719224`
- all-pair positives/negatives: `19` / `15`
- holdout range: `1.28288438842`
- covariance max error: `0`

## Null Controls
- identity_pairing: score `0`, ratio `0`, status `NULL_REJECTION_PASS`
- shuffled_pairing: score `0.217664424953`, ratio `0.529286490078`, status `NULL_REJECTION_PASS`
- random_orthogonal_pairing: score `0.302888955932`, ratio `0.73652381368`, status `NULL_REJECTION_FAIL`
- topology_destroyed_pairing: score `0`, ratio `0`, status `NULL_REJECTION_PASS`
- label_permuted_settings: score `0.316852730851`, ratio `0.770479006021`, status `NULL_REJECTION_FAIL`
- closure_vector_permutation: score `0.167704136943`, ratio `0.407799915092`, status `NULL_REJECTION_PASS`
- refinement_broken_control: score `0.498655720208`, ratio `1.21256257638`, status `NULL_REJECTION_FAIL`

## Boundary
- This is a geometry audit harness, not a Bell-score optimizer.
- No CHSH S value is computed here.
- A failed earlier gate blocks later promotion.
- HAOS Bell derivation is not established.
