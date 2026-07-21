# HAOS-IIP Branch Lifecycle Summary

> Generated from `branch_lifecycle_registry.json`. Edit the JSON authority, then regenerate this view.

Audit date: `2026-07-11`

This registry closes exhausted mechanisms without deleting their evidence. A closed branch may only be followed by a new candidate with a new identifier, a new precommitment, independent motivation or evidence, and the prior negative result preserved.

## Active Queue

| Priority | Candidate | Lifecycle | Implementation | Next gate |
| ---: | --- | --- | --- | --- |
| 1 | `EL-R3-THEORY-REVISION-01` | `ACTIVE_CANDIDATE` | `blocked` | `REQUIRED`: `docs/roadmaps/EL_R3_THEORY_REVISION_01.md` |
| 2 | `EL-R4-OC-02-INTERVENTION-NATIVE` | `ACTIVE_CANDIDATE` | `blocked` | `REQUIRED`: `experiments/emergence_ladder/rung4_operational_closure/oc_02_precommitment_contract.json` |
| 3 | `EL-R5-CSR-01` | `ACTIVE_CANDIDATE` | `blocked` | `REQUIRED`: `experiments/emergence_ladder/rung5_cross_scale_recovery/precommitment_contract.json` |

The queue is ranked by directness, falsifiability, control quality, holdout readiness, and dependency readiness. It is not ranked by which target is easiest to make pass.

## Closed And Quarantined

| Branch | Lifecycle | Frozen reading | Why work stops |
| --- | --- | --- | --- |
| `BELL-HAOS-B0-B3.2.2` | `TERMINAL_NEGATIVE` | TERMINAL FAIL / OPEN | The current G/J_lambda branch yields generic relational geometry without setting semantics or refinement specificity. CHSH scoring is not authorized. |
| `CELESTIAL-BOUNDARY-V1` | `TERMINAL_NEGATIVE` | OPEN / REQUIREMENTS NOT CLOSED | The current proxy path does not close the declared celestial requirements. Further work requires a genuinely new operational mapping. |
| `CST-V0.2.2` | `TERMINAL_NEGATIVE` | OPEN / TARGET NOT DISTINCT / UNDERPOWERED | The repaired instrument is partially discriminative, the target is not distinct, and the evidence is underpowered. Target-adjacent repair is closed. |
| `EL-R3-RP-01-DISTRIBUTED-PARITY` | `TERMINAL_NEGATIVE` | VALIDATION GATE FAILED / RUNG 3 UNSUPPORTED | Within-radius local parity decoding restores the relation/operator region, but functional recovery remains absent, the full recovery basin is empty, and target performance does not separate from RT-02 or matched controls. The current alphabet/code/bridge path is closed. |
| `EL-R3-RT-01` | `TERMINAL_NEGATIVE` | RUNG 3 NEGATIVE RESULT | The tested grade-locked transport moves away from its unperturbed trajectory, loses the grade signature, and does not beat passive or operator-only controls. |
| `EL-R3-RT-02-RELATIONAL-FEEDBACK` | `TERMINAL_NEGATIVE` | PARTIAL RECOVERY ONLY / RUNG 3 UNSUPPORTED | One-bit relational feedback repairs constraints and preserves identity but has zero fully recovered runs, no functional recovery, no qualifying basin, and failed validation/control-ablation gates. |
| `EL-R4-OC-01` | `TERMINAL_NEGATIVE` | RUNG 4 NEGATIVE RESULT | Aggregate internal-context advantage survives validation and holdout, but the run-level confidence interval includes zero. The compressed event-chain mechanism is closed. |
| `GEO-HIDDEN-01-V1` | `TERMINAL_NEGATIVE` | BENCHMARK OPEN / CURRENT FEATURE PATH CLOSED | Distance and orientation signals survive, but transformation-class recovery is not robust on holdout under the current feature representation. |
| `GEO-INTRINSIC-SPIRAL-V1` | `TERMINAL_NEGATIVE` | FEATURE BUNDLE SPECIFICITY INSUFFICIENT | The current feature bundle does not establish spiral-specific transfer beyond the strongest baseline. |
| `HBP-PB01-V1` | `QUARANTINED_INVALID` | PREDICTION NOT DISTINCT / CONTROL INVALID | The stored result is not distinct from baselines and reports invalid controls. It remains reproducible historical evidence only. |
| `HBP-PB02-V1` | `QUARANTINED_INVALID` | CONTROL INVALID | Control routing, holdout model selection, feature-slice consistency, and draft-contract enforcement are invalid in the current artifact. |
| `HBP-PB03-V1` | `QUARANTINED_INVALID` | PREDICTION NOT DISTINCT / INSTRUMENT INVALID | The baseline manifest is incomplete and seed/topology controls do not execute their declared perturbations through the target path. |
| `HBP-PB04-V1` | `QUARANTINED_INVALID` | PREDICTION NOT DISTINCT / INSTRUMENT INVALID | The baseline and HAOS candidate share the same implementation path, and declared controls do not establish valid transfer discrimination. |
| `PHYSICS-RAW-GRADIENT-V1` | `TERMINAL_NEGATIVE` | OPEN / RAW PROXY NOT RECOVERED | Raw local-gradient observables fail the current shape thresholds; only the separately bounded shell-native proxies survive. |
| `SCALE-BRIDGE-66.5-CURRENT-MECHANISM` | `TERMINAL_NEGATIVE` | OPEN / CURRENT MECHANISM EXHAUSTED | CP2, CP3, comparative, and CP5 gates remain open under the current mechanism. More diagnostics around the same object are not authorized. |

## Retained But Inactive

| Branch | Lifecycle | Role |
| --- | --- | --- |
| `BELL-CHSH-REFERENCE` | `FROZEN_REFERENCE` | The sidecar is a valid frozen calibration target, not a HAOS mechanism branch. |
| `BELL-SYNTHETIC-RELATIONAL-CALIBRATION` | `SUPPORTING_CALIBRATION` | The instrument calibration remains reusable, but it cannot reopen the terminal Bell mechanism branch. |
| `CORE-PHASE-SPINE-03-19` | `FROZEN_REFERENCE` | The public phase spine is a reference dependency, not an expansion target. |
| `EOT-01-SYNTHETIC-CALIBRATION` | `SUPPORTING_CALIBRATION` | EOT-01 recovers terms generated by construction and remains a valid method calibration, not an artifact-derived result. |
| `EOT-02-ARTIFACT-DERIVED` | `SUPPORTING_CALIBRATION` | The candidate remains scientifically useful but is deferred until enough independent refinements exist and lower recovery/closure rungs are supported. |
| `GEO-MM-01-ARITHMETIC-CALIBRATION` | `SUPPORTING_CALIBRATION` | The shared-support diagnostic passes with a narrow claim ceiling; local formal graph invariance remains unproved. |
| `GEO-MM-02-FORMAL-INVARIANCE` | `SUPPORTING_CALIBRATION` | The finite formal target remains legitimate but is deferred because it does not advance the current emergence frontier. |
| `GEO-SYNTHETIC-CHAIN-CALIBRATION` | `SUPPORTING_CALIBRATION` | Transformation, probability, and observable sub-harnesses remain useful calibration instruments even though external semantics are open. |
| `HBP-IR-01` | `SUPPORTING_CALIBRATION` | The minimal synthetic integrity run passes and detects forced leakage and no-op controls. It is retained for future predictive branches without promoting any historical HBP result. |
| `HYDROGEN-REFERENCE` | `FROZEN_REFERENCE` | The reference law reproduces; the HAOS derivation slot remains empty. |
| `SEMICONDUCTOR-REFERENCE` | `FROZEN_REFERENCE` | The sidecar is a reference calibration and contains no independent HAOS derivation. |
| `SPECULATIVE-BRIDGE` | `SPECULATIVE_ONLY` | Speculative notes may inspire questions but cannot enter the active evidence queue or reopen closed branches. |

## Enforcement

- Terminal and quarantined branches cannot authorize expansion or implementation.
- Active candidates cannot authorize implementation until their new precommitment exists and is marked `FROZEN`.
- Reference and calibration bundles remain reusable but cannot promote a mechanism claim.
- Reopening requires a new candidate ID, new precommitment, new mechanism or independent evidence, preserved prior results, and no threshold or post-hoc field movement.
- The registry governs development attention only. It does not rewrite any scientific result artifact.

Validation command:

```bash
uv run python scripts/check_branch_lifecycle.py
```
