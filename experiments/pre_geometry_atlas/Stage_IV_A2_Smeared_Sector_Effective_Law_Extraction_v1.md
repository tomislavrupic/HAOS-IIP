# Stage IV A2 Smeared Sector Effective Law Extraction v1

## Stage purpose
This stage runs the first true effective-law extraction pass for the surviving smeared-transfer sector on the frozen clustered DK branch. It compares a minimal family of reduced laws under the Stage 24.1 observable ledger and keeps only the smallest branch-valid law that achieves genuine predictive compression.

## Authority note
This note is authoritative for the stamped run `20260316_115346`.

Superseded run:
- `20260316_115312` is superseded.
- reason: tie-break selection was sharpened so equal-fit threshold candidates prefer the transport-side continuous separator over a binary closure-support proxy.
- interpretation change: yes. The authoritative selector became `flow_concentration_index <= 0.884308` instead of the earlier refinement-flag proxy.

## Authoritative inputs
- `23.11 / III-C2`: json=`data/20260316_104146_stage23_11_dk_minimal_effective_model_extraction.json`, csv=`data/20260316_104146_stage23_11_dk_minimal_effective_model_extraction.csv`
- `23.12 / III-C3`: json=`data/20260316_105448_stage23_12_dk_effective_phase_diagram_closure.json`, csv=`data/20260316_105448_stage23_12_dk_effective_phase_diagram_closure.csv`
- `23.13 / III-C4`: json=`data/20260316_110320_stage23_13_dk_mechanism_isolation_probe.json`, csv=`data/20260316_110320_stage23_13_dk_mechanism_isolation_probe.csv`
- `23.14 / III-C5`: json=`data/20260316_111201_stage23_14_phaseIII_final_consolidation_and_freeze.json`, csv=`data/20260316_111201_stage23_14_phaseIII_claims_ledger.csv`
- `23.15 / III-P1`: paper=`papers/drafts/Phase_III_Clustered_DK_Effective_Sector_Draft_v1.md`
- `24.1 / IV-A1`: json=`data/20260316_114317_stage24_1_smeared_sector_observable_ledger.json`, csv=`data/20260316_114317_stage24_1_smeared_sector_observable_ledger.csv`

## Law target definition
stable_closed_smeared = 1 iff Stage 23.12 derived_phase_label == smeared_transfer_phase; otherwise 0 for transient mixed controls on the frozen closure lattice.

## Candidate model families
- `IVA2_M1` `threshold_rule`: core=`flow_concentration_index` context=`none`
- `IVA2_M2` `simple_reduced_surrogate`: core=`topology_survival_time, refinement_stability_flag, weak_coupling_stability_flag, reverse_stability_flag, flow_concentration_index` context=`none`
- `IVA2_M3` `context_conditioned_threshold_rule`: core=`flow_concentration_index` context=`phase_corridor_position, phase_corridor_width`
- `IVA2_M4` `closure_support_conjunction`: core=`refinement_stability_flag, weak_coupling_stability_flag, reverse_stability_flag` context=`none`

## Consistency checks
- `stage23_11_effective_model_valid_true` = `True`
- `stage23_13_phase_corridor_primary` = `True`
- `stage23_14_freeze_valid_true` = `True`
- `only_accepted_primary_observables_used` = `True`
- `context_variables_not_promoted` = `True`
- `rejected_variables_unused` = `True`
- `winning_law_simpler_than_lookup_table` = `True`
- `winning_law_beats_naive_baseline` = `True`
- `winning_law_smeared_side_consistent` = `True`

## Model comparison results
- `IVA2_M1`: balanced_accuracy=`1.0000`, baseline_delta=`0.5000`, complexity=`2`, branch_valid=`True`, selected=`True`
- `IVA2_M2`: balanced_accuracy=`1.0000`, baseline_delta=`0.5000`, complexity=`7`, branch_valid=`True`, selected=`False`
- `IVA2_M3`: balanced_accuracy=`1.0000`, baseline_delta=`0.5000`, complexity=`5`, branch_valid=`True`, selected=`False`
- `IVA2_M4`: balanced_accuracy=`1.0000`, baseline_delta=`0.5000`, complexity=`3`, branch_valid=`True`, selected=`False`

## Selected effective law
- selected model: `IVA2_M1`
- law: `stable_closed_smeared = 1 iff flow_concentration_index <= 0.884308`
- core inputs: `flow_concentration_index`
- fit: accuracy=`1.0000`, balanced_accuracy=`1.0000`, baseline_delta=`0.5000`

## Failed candidates and why they failed
- `IVA2_M2`: unselected because a simpler branch-valid model achieved equal or better compression; notes=Least-squares linear score model over the accepted core observables.
- `IVA2_M3`: unselected because a simpler branch-valid model achieved equal or better compression; notes=Two-bin context-conditioned threshold rule. Phase enters only as an external selector, not as a state variable.
- `IVA2_M4`: unselected because a simpler branch-valid model achieved equal or better compression; notes=Minimal conjunction over the binary closure-support observables.

## Minimal interpretation
On the frozen clustered DK branch, the smeared sector admits a compact effective law built from the accepted Phase IV observable core. The selected law predicts branch-valid smeared-sector persistence/stability without reintroducing rejected braid- or asymmetry-based quantities. Phase-corridor variables remain external conditioning coordinates rather than promoted state variables.

## Output artifact list
- stamped JSON summary: `data/20260316_115346_stage24_2_smeared_sector_effective_law_extraction.json`
- stamped CSV model ledger: `data/20260316_115346_stage24_2_smeared_sector_effective_law_model_ledger.csv`
- plot: `plots/20260316_115346_stage24_2_law_family_comparison_table.png`
- plot: `plots/20260316_115346_stage24_2_minimal_law_decision_panel.png`
- plot: `plots/20260316_115346_stage24_2_prediction_vs_truth_panel.png`
- plot: `plots/20260316_115346_stage24_2_context_conditioning_panel.png`

## Commit recommendation
- recommended status: commit
- rationale: Phase IV now has a valid reduced law candidate for the smeared sector using only the frozen observable core.
- branch effect: enables 24.3 quasi-invariant probing from a compact law basis.
- caution: rejected variables remain excluded unless a later explicit branch changes the measurement regime.
