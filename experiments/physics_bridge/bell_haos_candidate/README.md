# Bell HAOS-IIP Candidate Bridge

This isolated bundle builds the first HAOS-IIP-to-Bell candidate bridge with the
frozen CHSH sidecar acting only as a scoreboard convention.

Implemented stages:

- B0: neutral adapter records and diagnostics.
- B1: deliberately local sanity candidate.
- B2: HAOS-IIP-local recoverability candidate using frozen Phase 18 source rows.
- B3: explicit joint-closure-cost candidate using frozen Phase 18 source rows.
- B3 precommitment: present, generated before scoring, with factorization and
  outcome independence explicitly relaxed.
- B3.2: relational geometry audit harness. It tests for a sign-changing
  invariant `G(a,b,lambda)` before CHSH scoring is allowed.
- B3.2.1: null failure localization. It keeps the B3.2 invariant frozen and
  compares target/null G distributions, relation matrices, topology dependence,
  operator variants, and pair ordering.
- B3.2.2: semantic and refinement specificity audit. It keeps B3.2 frozen and
  asks whether the relation matrix recovers setting semantics and refinement
  continuity.
- v0.3 freeze: terminal negative result for the current `G / J_lambda` branch.

It is not a physical Bell experiment, not a loophole-free Bell test, not a
derivation of quantum mechanics, and not evidence that CST recovers Bell
correlations.

## Run

```bash
uv run python experiments/physics_bridge/bell_haos_candidate/run_bell_haos_candidate.py
uv run python experiments/physics_bridge/bell_haos_candidate/run_b3_geometry_audit.py
uv run python experiments/physics_bridge/bell_haos_candidate/run_b3_null_failure_localization.py
uv run python experiments/physics_bridge/bell_haos_candidate/run_b3_semantic_refinement_audit.py
```

## Outputs

- `precommitment_contract.json`
- `assumption_ledger.json`
- `candidate_trials.csv`
- `candidate_correlations.json`
- `no_signalling_diagnostics.csv`
- `setting_independence_diagnostics.csv`
- `control_results.csv`
- `bell_candidate_report.md`
- `bell_candidate_result.json`
- `b3_precommitment_contract.json`
- `b3_precommitment_report.md`
- `b3_joint_cost_audit.csv`
- `b3_2_precommitment_contract.json`
- `b3_2_geometry_matrix.csv`
- `b3_2_covariance_diagnostics.csv`
- `b3_2_control_results.csv`
- `b3_2_geometry_report.md`
- `b3_2_result.json`
- `b3_2_1_precommitment_contract.json`
- `b3_2_1_distribution_comparison.csv`
- `b3_2_1_matrix_comparison.csv`
- `b3_2_1_topology_dependence.csv`
- `b3_2_1_operator_dependence.csv`
- `b3_2_1_pair_ordering.csv`
- `b3_2_1_null_invariant_ledger.json`
- `b3_2_1_report.md`
- `b3_2_1_result.json`
- `b3_2_2_precommitment_contract.json`
- `b3_2_2_semantic_ordering.csv`
- `b3_2_2_semantic_edges.csv`
- `b3_2_2_refinement_persistence.csv`
- `b3_2_2_refinement_pairs.csv`
- `b3_2_2_report.md`
- `b3_2_2_result.json`
- `bell_bridge_v0_3_freeze.json`
- `bell_bridge_v0_3_technical_note.md`

## Boundary

A numerical `S > 2` is not sufficient for a positive result. Any candidate must
also pass no-signalling, setting independence or explicitly declared relaxation,
no target leakage, no post-selection, seed reproducibility, and matched local
controls.

## B3 Gate

`b3_precommitment_contract.json` records an explicit model:

`P(A,B | a,b,lambda) = 0.5 * P(product=A*B | a,b,lambda)`

The product probability is generated from a precommitted Phase 18 chain-order
closure cost. This relaxes factorization and outcome independence by design.
Measurement independence and no-signalling remain diagnostics, not assumptions
granted for free.

## B3.2 Geometry Gate

`run_b3_geometry_audit.py` does not compute CHSH. It applies the gate ladder:

- G1 relational sensitivity
- G2 sign-changing structure
- G3 covariance
- G4 holdout transfer
- G5 null rejection
- G6 CHSH scoreboard authorization

A failed earlier gate returns `CHSH_SCORING_NOT_AUTHORIZED`. This is a geometry
audit harness, not a Bell-score optimizer.

## B3.2.1 Null Failure Localization

`run_b3_null_failure_localization.py` does not change B3.2. It localizes why
nulls preserve sign-changing geometry. Current diagnostics classify the failure
as target-specificity/provenance loss, with CHSH scoring still unauthorized.

## B3.2.2 Provenance Gate

`run_b3_semantic_refinement_audit.py` has two independent gates:

- S1 semantic ordering: target G must recover intended setting relations better
  than label-permuted controls.
- S2 refinement persistence: target matrices must persist across refinements
  better than refinement-broken controls.

Dual failure returns `GENERIC_RELATIONAL_GEOMETRY`. CHSH remains forbidden.

## v0.3 Freeze

The current `G / J_lambda` construction is terminal:

`GENERIC_RELATIONAL_GEOMETRY`

It produces genuine relational structure, but not structure specific enough to
support a HAOS-IIP Bell bridge. No B3.3 may reuse the same `G`, `J_lambda`,
setting map, signature representation, thresholds, or normalization.
