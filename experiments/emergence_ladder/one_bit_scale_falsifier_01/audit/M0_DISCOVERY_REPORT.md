# EL-R3-ONE-BIT-SCALE-FALSIFIER-01 M0 Discovery Report

Date: `2026-08-18`

Terminal M0 decision: `BLOCKED_BY_THEORY_REVISION`

## Starting Repository State

- Repository: `/Volumes/Samsung T5/2026/HAOS/HAOS DOCS/HAOS-IIP`
- Branch: `main`
- Starting commit: `f4742c49f4a50b18f7b92b98b18f176cd1637987`
- Upstream relation: `main...origin/main`
- Instructions read completely: `AGENTS.md`
- Pre-existing modified files, preserved and not edited by this audit:
  - `README.md`
  - `docs/notes/foundations/HAOS_IIP_Claim_Gating_Template_v1.md`
  - `docs/notes/foundations/HAOS_IIP_Falsification_Engine_v1.md`
  - `docs/notes/foundations/HAOS_IIP_Scale_Bridge_Legitimacy_Audit_v1.md`
- Pre-existing untracked files, preserved and not edited by this audit (the
  repository's configured short status initially hid untracked files; their
  presence was independently established before this bundle was written by
  reading and executing the existing claim-contract checker):
  - `docs/claim_governance/README.md`
  - `docs/claim_governance/claim_contract.schema.json`
  - `docs/claim_governance/scale_bridge_66_5.claim.json`
  - `scripts/check_claim_contracts.py`
  - `tests/test_claim_contracts.py`
- Starting combined SHA-256 of the tracked frozen Rung 3, telemetry, Phase
  V/VI/X/XI, and lifecycle slice:
  `7651740fbebe8d9c9a7e728f50617c225ec719a826c6be381a6e2a7e25aa0c66`

The fingerprint was computed over tracked files in the three frozen Rung 3
bundles, `telemetry/frozen_metrics.py`, Phase V, Phase VI, Phase X, Phase XI,
and `docs/branch_governance/`.

## Repository Map

| Requested item | Exact authority or closest artifact | M0 classification | Finding |
| --- | --- | --- | --- |
| One-bit relation representation | `experiments/emergence_ladder/rung3_recovery_trajectory_v2/mechanism.py` (`RelationalMemory.orientation`) | present and reusable only as frozen historical mechanism | One sign bit is stored for each valid oriented edge, but incidence, valid-edge masks, and addressing are also available to the mechanism. |
| RT-02 mechanism | `experiments/emergence_ladder/rung3_recovery_trajectory_v2/mechanism.py` | present and reusable only for an authorized audit | The correction is local cochain-constraint feedback. Re-executing it as the next Rung 3 candidate is forbidden by current lifecycle governance. |
| Decoder/update rule | `correction_vector` and `simulate_feedback` in the same file | present and reusable only for an authorized audit | The frozen rule uses `eta * D^T [availability * s * max(0, margin - s * D x)] / max_degree`. |
| Frozen parameters | `calibration/parameter_selection.json`, `calibration/derived_thresholds.json`, and `precommitment_contract.json` in `rung3_recovery_trajectory_v2/` | present and reusable | `eta=0.7`, `inequality_margin=0.05`, `steps=80`, `edge_dropout_fraction=0.2`, minimum encoded edge magnitude `0.08`. |
| Perturbation construction | `experiments/emergence_ladder/rung3_recovery_trajectory_v2/fixtures.py` | ambiguous across refinement | The three RT-02 families exist, but their support rules are not frozen as scale-equivalent semantics. A fixed `2x2` block and a node-count-proportional sparse corruption imply different normalized interventions. |
| RT-02 seeds | `rung3_recovery_trajectory_v2/precommitment_contract.json` and `final/seed_registry.json` | present but incompatible with scaling | Calibration `6101-6106`, validation `6201-6204`, final `6301-6308`; all apply only to the `8x8` system. No cross-level seed-generation policy is frozen. |
| Passive control | `rung3_recovery_trajectory_v2/controls.py` (`passive_relaxation`) | present and reusable | Same-size RT-02 path only. |
| Ordinary filtering control | same file (`operator_only_filtering`) | present and reusable | Same-size RT-02 path only. |
| Trivial-attractor control | same file (`trivial_attractor`) | present and reusable | Same-size RT-02 path only. |
| Altered-connectivity control | same file (`topology_altered`) | present and reusable | Same-size RT-02 path only. |
| Matched random-bit/local-correction control | `signal_blind` and RP-01 `memory_budget_random_bits` are the closest artifacts | present but incompatible | Neither is a preregistered, common refinement-level implementation of the exact requested control. |
| Shuffled-relation-bit control | RT-02 `identity_scrambled` and RP-01 parity permutation controls | present but incompatible | Available only at `8x8`; matching of addressing and compute budgets is not audited across scale. |
| Independent functional target | `experiments/emergence_ladder/rung3_distributed_parity/FUNCTIONAL_TARGET.md` | present but incompatible with a scale claim | Four low-frequency projection responses were frozen for RP-01. They are independent of parity, but no common normalization or acceptance rule is frozen across the Phase VI/X/XI hierarchy. H2 must not be scored from current artifacts. |
| Aggregate generation | `rung3_recovery_trajectory_v2/run_experiment.py` and `rung3_distributed_parity/run_experiment.py` | present but single-level | Both aggregate valid frozen experiments at `n_side=8`; neither aggregates refinement scaling. |
| Semantic checkers | `rung3_recovery_trajectory_v2/check_bundle.py`, `rung3_distributed_parity/check_bundle.py` | present but incompatible | They validate their own frozen bundles, not scale invariance, information overhead, common control routing, or three-level coverage. |
| Lifecycle registry | `docs/branch_governance/branch_lifecycle_registry.json` | present and authoritative | It is the controlling M0 blocker. |
| Lifecycle summary | `docs/branch_governance/branch_lifecycle_summary.md` | present and generated | It places `EL-R3-THEORY-REVISION-01` first and implementation-blocked. |
| Emergence ladder | `docs/roadmaps/HAOS_IIP_EMERGENCE_LADDER_2026-07-11.md` | present and authoritative for direction | It explicitly authorizes no next Rung 3 experiment until theory revision supplies independently justified mesoscopic information. |
| Required theory-revision precommitment | `docs/roadmaps/EL_R3_THEORY_REVISION_01.md` | missing | The registry marks this path `REQUIRED`; the file does not exist at M0. |
| Frozen telemetry | `telemetry/frozen_metrics.py` | present and reusable | It defines overlap, localization width, concentration retention, participation ratio, frozen recovery score, and persistence helpers. It does not define relation identity, operator repair, functional restoration, or intervention accounting. |
| Phase V evidence | `phase5-readout/phase5_authoritative_manifest.json` and `phase5-readout/runs/phase5_runs_ledger_latest.json` | present but incompatible | It supports recovery-distribution separation for stable, degraded, and shuffled classes, not RT-02 scale execution. |
| Frozen operator/refinement hierarchy | `phase6-operator/phase6_operator_manifest.json` and `phase6_refinement_ledger.csv` | present but incompatible | Four admissible cochain-Laplacian levels exist at `n_side=12,24,36,48`, but no RT-02 mechanism, matched RT-02 controls, or intervention budget was run on them. |
| Phase X scale bridge | `phase10-bridge/phase10_manifest.json` and `runs/phase10_effective_scaling_ledger.csv` | present but incompatible | Five operator/descriptor levels exist at `n_side=12,24,36,48,60`; they are not one-bit recovery trials. |
| Phase XI persistence | `phase11-protection/phase11_manifest.json`, `runs/phase11_perturbation_survival_ledger.csv`, and `runs/phase11_persistence_scaling.csv` | present but incompatible | Four levels exist at `n_side=48,60,72,84`; they measure persistence under Phase XI perturbations, not recovery by RT-02. |
| Same-surrogate audits | `continuum-sketch/same_surrogate_coarse_graining_recovery.csv` and `same_surrogate_control_integrity_report.md` | present but adverse/incompatible | The report records failed same-surrogate control integrity; it cannot supply missing one-bit scale runs. |
| Literal information-budget accounting | RT-02 contract, mechanism schema, RT-02 run results, and RP-01 `information_accounting.json` | missing for the requested scale audit | Representational orientation bits are bounded at one size, but addressing/index overhead, valid-edge masks, incidence access, control-policy information, hyperparameter information, correction-event counts, wall time, and peak memory are not jointly accounted by refinement level. |

## Frozen Parent Evidence

- `EL-R3-RT-01`: frozen `NEGATIVE_RESULT`; the tested transport did not beat
  passive or operator-only controls.
- `EL-R3-RT-02`: frozen `PARTIAL_RECOVERY_ONLY`; result hash
  `el_r3_rt_02_cade4e057a54d119e6af143c`. It has median relational identity
  `1.0` and operator recovery `0.919147`, but median functional recovery `0.0`,
  target recovery rate `0.0`, and no qualifying recovery basin.
- `EL-R3-RP-01`: frozen `VALIDATION_GATE_FAILED`; result hash
  `el_r3_rp_01_e7806ce1a296242b75ad3d3b`. Final seed count consumed is `0`.
- Highest supported emergence rung remains Rung 2. Rung 3 remains
  `NEGATIVE_RESULT` for the tested mechanism families.

## Artifact Sufficiency Matrix

Legend: `P` present; `I` present but incompatible with the proposed one-bit
scale audit; `M` missing; `A` ambiguous. “Bits/interventions” distinguishes
recorded representational bits from the complete information and intervention
account demanded by this audit.

| Artifact / level | System size | Hierarchy depth | State dimension | Seed coverage | Perturbation family | Mechanism parameters | Bits / interventions | Relational repair | Operator repair | Functional restoration | Passive | Filter | Trivial attractor | Altered connectivity | Uncertainty / repeats | Frozen artifact |
| --- | ---: | ---: | ---: | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RT-02 only level | `8x8` / 64 nodes | 1 | 64 scalars | P: `6101-6308` partitions | P: 3 families | P: frozen | A: orientation bits and corrective cost only; overhead/events M | P | P | P, result 0 | P | P | P | P | P: bootstrap | P |
| RP-01 only level | `8x8` / 64 nodes | 1 | 64 scalars | P: `7101-7204`; final `7301-7308` untouched | P: 8 families | P: frozen | P for parity/checkpoint comparison; scale overhead M | P | P | P, result 0 | P | P | P | P | P: seed-block bootstrap | P |
| Phase VI R1 | `12x12` | 1 of 4 | 576 cochain components / 144 zero forms | M for RT-02 | I | M | M | M | P: operator diagnostics only | M | M | M | M | M | deterministic artifact only | P |
| Phase VI R2 | `24x24` | 2 of 4 | 2304 / 576 zero forms | M for RT-02 | I | M | M | M | P: operator diagnostics only | M | M | M | M | M | deterministic artifact only | P |
| Phase VI R3 | `36x36` | 3 of 4 | 5184 / 1296 zero forms | M for RT-02 | I | M | M | M | P: operator diagnostics only | M | M | M | M | M | deterministic artifact only | P |
| Phase VI / X R4 | `48x48` | 4 | 9216 / 2304 zero forms | M for RT-02 | I | M | M | M | P: operator and Phase XI persistence diagnostics | M | M | M | M | P only in Phase XI, not same path | I | P |
| Phase X R5 / Phase XI 60 | `60x60` | 5 / 2 of XI | operator hierarchy artifact | M for RT-02 | I | M | M | M | P: scale/persistence diagnostics | M | M | M | M | P only in Phase XI, not same path | I | P |
| Phase XI 72 | `72x72` | 3 of XI | Phase XI state artifact | M for RT-02 | I | M | M | M | P: persistence diagnostics | M | M | M | M | P, not same path | I | P |
| Phase XI 84 | `84x84` | 4 of XI | Phase XI state artifact | M for RT-02 | I | M | M | M | P: persistence diagnostics | M | M | M | M | P, not same path | I | P |

No row set supplies three comparable levels containing the same one-bit
mechanism, same semantic perturbation, same controls through the same path,
the independent function, uncertainty, and complete information accounting.

## Sufficiency Decision

1. Artifact-only retrospective scale audit: **not supported**. The only frozen
   one-bit records are at `n_side=8`. The refinement ledgers measure different
   mechanisms and state spaces.
2. New preregistered execution: **scientifically designable but not currently
   admissible**. It would need a theory-revision closure, an authorized new
   lifecycle entry, at least three common levels, scale-invariant perturbation
   semantics, new seed partitions, complete bit/address/intervention
   accounting, and a scale-normalized functional contract.
3. Current terminal route: **neither**, because governance forbids
   implementation and the frozen record is insufficient for an artifact-only
   classification.

## Governance Decision

The authoritative registry entry for `EL-R3-THEORY-REVISION-01` has:

- lifecycle status `ACTIVE_CANDIDATE`;
- implementation authorization `false`;
- next precommitment path `docs/roadmaps/EL_R3_THEORY_REVISION_01.md` with
  status `REQUIRED`;
- a stop condition rejecting proposals that merely re-encode RT-02 orientation
  or RP-01 parity without a new information representation.

The required file is absent. The emergence ladder independently states that no
new Rung 3 experiment is authorized until theory revision explains and freezes
function-relevant mesoscopic restorative information. The proposed experiment
tests the unchanged RT-02 orientation mechanism, so it cannot satisfy that
dependency by construction.

M0 therefore stops with `BLOCKED_BY_THEORY_REVISION`. The authoritative
lifecycle registry was not edited. The accompanying registry JSON is a proposal
for review only. The accompanying contract is a non-frozen draft. No final
seeds were inspected or consumed.

## Baseline Validation

| Command | Outcome |
| --- | --- |
| `uv sync` | pass; 14 packages resolved, 13 audited |
| `uv run python examples/quick_reproduce.py` | pass; reproduction spine passed |
| `uv run python run_phase.py 10 --check` | pass; `success: true` |
| `uv run python run_phase.py 18 --check` | pass; all declared gates true |
| `uv run python run_phase.py 19 --check` | pass; all declared gates true |
| `uv run python scripts/check_branch_lifecycle.py` | pass; status `ok`; active priority begins with theory revision |
| RT-01 unit-test discovery | pass; 5 tests |
| RT-02 unit-test discovery | pass; 12 tests |
| RP-01 unit-test discovery | pass; 15 tests |
| `uv run python -m unittest tests.test_branch_lifecycle -v` | pass; 9 tests |
| `uv run python experiments/emergence_ladder/rung3_recovery_trajectory_v2/check_bundle.py` | pass; frozen `PARTIAL_RECOVERY_ONLY` hash verified |
| `uv run python experiments/emergence_ladder/rung3_distributed_parity/check_bundle.py` | pass; frozen `VALIDATION_GATE_FAILED`, 0 final seeds consumed |
| `uv run python scripts/check_claim_contracts.py` | pass; 1 claim contract validated |

## M0 Claim Ceiling

Demonstrated: repository governance and artifact coverage are sufficient to
reject execution of the proposed audit in the current lifecycle state.

Not demonstrated: scale stability, scale failure, intervention-cost growth,
functional recovery, or any finite-size scaling law of the one-bit mechanism.
