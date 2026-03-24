# Stage IV A1 Smeared Sector Observable Ledger v1

## Stage purpose
This note defines the minimal operational observable basis of the surviving smeared-transfer sector on the frozen clustered Dirac-Kahler branch. It is a reduction-only ledger extraction from the frozen Phase III closure outputs.

## Frozen inputs
- `23.11 / III-C2`: json=`data/20260316_104146_stage23_11_dk_minimal_effective_model_extraction.json`, csv=`data/20260316_104146_stage23_11_dk_minimal_effective_model_extraction.csv`
- `23.12 / III-C3`: json=`data/20260316_105448_stage23_12_dk_effective_phase_diagram_closure.json`, csv=`data/20260316_105448_stage23_12_dk_effective_phase_diagram_closure.csv`
- `23.13 / III-C4`: json=`data/20260316_110320_stage23_13_dk_mechanism_isolation_probe.json`, csv=`data/20260316_110320_stage23_13_dk_mechanism_isolation_probe.csv`
- `23.14 / III-C5`: json=`data/20260316_111201_stage23_14_phaseIII_final_consolidation_and_freeze.json`, csv=`data/20260316_111201_stage23_14_phaseIII_claims_ledger.csv`

## Stable object and contrast policy
- stable smeared-sector cells = `5`
- contrast transient-mixed cells = `10`
- stable phase span = `0.375 -> 0.575`
- derived phase-corridor width = `0.200`
- contrast controls are used only to test separation; they are not promoted into the Phase IV core.

## Validation checks
- `stage23_11_effective_model_valid_true` = `True`
- `stage23_12_phase_diagram_closed_true` = `True`
- `stage23_12_stable_region_is_smeared_only` = `True`
- `stage23_13_phase_corridor_primary` = `True`
- `stage23_14_freeze_valid_true` = `True`
- `primary_basis_smaller_than_candidate_list` = `True`
- `primary_basis_distinguishes_stable_region` = `True`
- `no_braid_quantity_promoted_to_primary` = `True`

## Primary observable basis
- `topology_survival_time`: sensitivity=`0.715`; stable=`0.0626 -> 0.2501`; contrast=`0.1253 -> 0.5011`; rationale=This is the minimal persistence-carrying observable on the closed branch. It is not sufficient alone, but it is required to carry smeared-sector lifetime structure into 24.2.
- `refinement_stability_flag`: sensitivity=`1.000`; stable=`1.0000 -> 1.0000`; contrast=`0.0000 -> 0.0000`; rationale=Refinement stability is a necessary closure-support observable and perfectly separates the stable smeared cells from the transient mixed controls on the frozen lattice.
- `weak_coupling_stability_flag`: sensitivity=`0.600`; stable=`1.0000 -> 1.0000`; contrast=`0.0000 -> 1.0000`; rationale=The stable smeared component is beta-stable by construction. This flag is not sufficient on its own, but it is indispensable because the closed regime is defined only after the weak-coupling consistency check is passed.
- `reverse_stability_flag`: sensitivity=`1.000`; stable=`1.0000 -> 1.0000`; contrast=`0.0000 -> 0.0000`; rationale=Reverse stability is required by the Phase III closure rule and cleanly separates the stable smeared cells from the transient mixed controls.
- `flow_concentration_index`: sensitivity=`1.000`; stable=`0.8691 -> 0.8824`; contrast=`0.8862 -> 0.8997`; rationale=This is the strongest continuous transport-like separator on the frozen branch: the stable smeared cells occupy a clean low-flow-concentration band disjoint from the transient mixed controls.

## Secondary and derived observables
- `phase_corridor_position`: status=`secondary_observable`; stable=`0.3750 -> 0.5750`; contrast=`0.3750 -> 0.5750`; rationale=Phase position organizes where the stable smeared row sits, but the same sampled phases also appear in transient mixed controls, so it is contextual rather than a primary separator.
- `grade_exchange_coherence_level`: status=`secondary_observable`; stable=`0.4348 -> 0.4430`; contrast=`0.4329 -> 0.4410`; rationale=Coherence remains well-defined across the smeared sector, but its stable and contrast ranges overlap too strongly to justify promotion into the minimal Phase IV core.
- `phase_corridor_width`: status=`derived`; stable=`0.200`; contrast=`n/a`; rationale=Corridor width is a boundary summary of the stable smeared component, not a per-cell operational observable.
- `smeared_regime_membership_flag`: status=`derived`; stable=`1`; contrast=`0`; rationale=This is the target label produced by the Phase III closure pass. It is necessary as an outcome variable, but it is not itself a measurable Phase IV observable.

## Explicit reject table

| Observable | Reason rejected |
| --- | --- |
| `transfer_asymmetry` | No branch-stable transfer-asymmetry quantity is present in the frozen Phase III closure outputs without introducing a new bookkeeping transport model. |
| `coherence_decay_rate` | The frozen closure outputs retain coherence levels, not a regime-stable coherence-decay-rate observable. Promoting a decay metric here would require a new analysis layer not licensed by Phase III. |
| `topology_return_error` | Return error exists only as a branch-specific ablation diagnostic in Stage 23.13 and does not define the stable smeared sector across the frozen closure lattice. |
| `braid_like_exchange_indicator` | Braid-centered quantities are explicitly excluded from promotion into the Phase IV core. They may appear only as contrast descriptors, because the family-wide braid claim and stable braid-phase claim failed in Phase III. |
| `mechanism_sensitivity_score` | Mechanism sensitivity is tied to the ablation design of Stage 23.13 and does not characterize the closed smeared sector itself. |

## Carry-forward set for 24.2
- core carry-forward observables: `topology_survival_time, refinement_stability_flag, weak_coupling_stability_flag, reverse_stability_flag, flow_concentration_index`
- context-only carry-forward observables: `phase_corridor_position, phase_corridor_width`
- excluded braid-centered quantities remain outside the Phase IV core and may appear only as contrast descriptors if explicitly needed.

## Final read
On the frozen clustered DK branch, the closed smeared-transfer regime is minimally characterized by refinement stability, weak-coupling stability, reverse stability, topology survival time, and flow concentration index, while phase-corridor position and grade-exchange coherence remain secondary or derived descriptors.

## Output artifact list
- stamped JSON summary: `data/20260316_114317_stage24_1_smeared_sector_observable_ledger.json`
- stamped CSV ledger: `data/20260316_114317_stage24_1_smeared_sector_observable_ledger.csv`
- plot: `plots/20260316_114317_stage24_1_smeared_sector_observable_table.png`
- plot: `plots/20260316_114317_stage24_1_smeared_sector_correlation_panel.png`
- plot: `plots/20260316_114317_stage24_1_smeared_sector_stability_panel.png`
