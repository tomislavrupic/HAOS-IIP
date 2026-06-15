# Synthetic Relational Geometry Calibration v0.1

Implemented fact: this suite calibrates semantic/refinement diagnostics on known synthetic geometry.
Design choice: the synthetic domain has known positives, destructive controls, and an ambiguous partial-preservation case.
Heuristic: success means the instrument distinguishes known structure from declared controls, not that any physics claim is supported.
Unverified hypothesis: no HAOS-IIP Bell derivation is tested here.

## Labels
- KNOWN_POSITIVE_PASS
- SEMANTIC_CONTROL_REJECTED
- TOPOLOGY_CONTROL_REJECTED
- REFINEMENT_CONTROL_REJECTED
- AMBIGUOUS_CASE_OPEN
- CALIBRATION_PASS
- NO_PHYSICAL_EXPERIMENT
- HAOS_BELL_DERIVATION_NOT_ESTABLISHED

## Conditions
- known_positive: expected `PASS`, observed `PASS`
- semantic_permutation_control: expected `FAIL_SEMANTIC`, observed `FAIL_SEMANTIC`
- topology_destroyed_control: expected `FAIL_SEMANTIC`, observed `FAIL_BOTH`
- refinement_broken_control: expected `FAIL_REFINEMENT`, observed `FAIL_BOTH`
- ambiguous_partial_preservation: expected `OPEN`, observed `OPEN`

## Boundary
- This is synthetic calibration, not a physics bridge.
- No Bell/CHSH score is computed.
- Passing calibration does not establish HAOS-IIP physics recovery.
