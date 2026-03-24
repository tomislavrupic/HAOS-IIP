# Stage IV A3 Stable Smeared Quasi-Invariant Probe v1

## Stage purpose
This stage probes whether the law-defined stable smeared sector on the frozen clustered Dirac-Kahler branch carries any branch-stable quasi-invariant, bounded composite, or preserved ordering relation. It is a reduction-only pass anchored to the Stage 24.2 threshold law.

## Authoritative inputs
- `23.12 / III-C3`: json=`data/20260316_105448_stage23_12_dk_effective_phase_diagram_closure.json`, csv=`data/20260316_105448_stage23_12_dk_effective_phase_diagram_closure.csv`
- `23.14 / III-C5`: json=`data/20260316_111201_stage23_14_phaseIII_final_consolidation_and_freeze.json`, csv=`data/20260316_111201_stage23_14_phaseIII_claims_ledger.csv`
- `24.1 / IV-A1`: json=`data/20260316_114317_stage24_1_smeared_sector_observable_ledger.json`, csv=`data/20260316_114317_stage24_1_smeared_sector_observable_ledger.csv`
- `24.2 / IV-A2`: json=`data/20260316_115346_stage24_2_smeared_sector_effective_law_extraction.json`, csv=`data/20260316_115346_stage24_2_smeared_sector_effective_law_model_ledger.csv`

## Law-anchored positive and contrast sets
- selected threshold law: `stable_closed_smeared = 1 iff flow_concentration_index <= 0.884308`
- positive stable smeared cells = `5`
- full contrast cells = `10`
- near-threshold contrast cells = `5`
- positive rows are defined by the intersection of the Stage 24.2 threshold law and the frozen Stage 23.12 stable-smeared truth label.
- contrast rows are transient mixed controls above the threshold boundary; the near-threshold subset is retained as the main failure comparison.

## Candidate quasi-invariant families
- `IVA3_Q1` `support_signature`: Exact closure-support signature built from refinement, weak-coupling, and reverse stability flags.
- `IVA3_Q2` `restated_law_band`: Flow-threshold band reported only as a direct restatement of the Stage 24.2 law.
- `IVA3_Q3` `ordering_relation`: Inverse ordering relation between flow concentration and topology survival time.
- `IVA3_Q4` `affine_envelope`: Stable-sector affine envelope linking topology survival time to flow concentration.
- `IVA3_Q5` `bounded_quantity`: Single-observable survival-time envelope tested as a quasi-invariant candidate.

## Consistency checks
- `stage24_1_observable_ledger_valid_true` = `True`
- `stage24_2_effective_law_found_true` = `True`
- `stage24_2_selected_core_is_flow_only` = `True`
- `law_matches_truth_on_closure_lattice` = `True`
- `positive_set_is_smeared_only` = `True`
- `no_forbidden_observables_in_candidates` = `True`
- `no_braid_quantity_promoted` = `True`
- `selected_candidate_uses_core_only` = `True`
- `selected_candidate_not_law_restatement` = `True`
- `phaseIII_smeared_side_still_frozen` = `True`

## Candidate comparison results
- `IVA3_Q1`: stable_score=`1.0000`, contrast_score=`0.0000`, separation=`1.0000`, selected=`False`
- `IVA3_Q2`: stable_score=`0.9948`, contrast_score=`0.9954`, separation=`-0.0007`, selected=`False`
- `IVA3_Q3`: stable_score=`0.9747`, contrast_score=`0.0000`, separation=`0.9747`, selected=`True`
- `IVA3_Q4`: stable_score=`0.9613`, contrast_score=`0.4075`, separation=`0.5539`, selected=`False`
- `IVA3_Q5`: stable_score=`0.7240`, contrast_score=`0.7983`, separation=`-0.0742`, selected=`False`

## Selected quasi-invariant
- selected candidate: `IVA3_Q3` `inverse_flow_survival_ordering`
- expression: `lower flow_concentration_index orders with longer topology_survival_time inside the stable smeared sector`
- stable metric: `spearman_rho=-0.9747`
- contrast metric: `spearman_rho=0.5476`
- separation score: `0.9747`

## Rejected and secondary candidates
- `IVA3_Q1`: reported as a support signature, but not selected because it does not add structure beyond the existing law/closure support; notes=Exact branch-local support signature, but it rephrases the closure-support flags already available in Stage 24.2 rather than adding a new relation.
- `IVA3_Q2`: reported but rejected because it only restates the Stage 24.2 threshold law; notes=Reported for completeness only. This is the Stage 24.2 law itself and may not be selected as the quasi-invariant read.
- `IVA3_Q4`: left as a secondary candidate because a simpler or stronger branch-valid relation dominated; notes=Compact stable-sector envelope. The fit is tight on the positive set but is treated as secondary to the nonparametric ordering relation because it is fitted on only five stable cells.
- `IVA3_Q5`: rejected because it is not more stable inside the positive set than in the contrast controls; notes=Reported as a bounded quantity candidate, but the survival envelope alone is not more stable than the contrast controls.

## Minimal interpretation
On the frozen clustered DK branch, the law-defined stable smeared sector carries a branch-local preserved ordering relation: lower flow concentration co-occurs with longer topology survival across the stable set, while the contrast controls do not preserve that inverse order. This read uses only accepted core observables, does not restate the threshold law, and remains bounded to the stable smeared sector.

## Output artifact list
- stamped JSON summary: `data/20260316_120334_stage24_3_stable_smeared_quasi_invariant_probe.json`
- stamped CSV candidate ledger: `data/20260316_120334_stage24_3_stable_smeared_quasi_invariant_ledger.csv`
- plot: `plots/20260316_120334_stage24_3_quasi_invariant_candidate_table.png`
- plot: `plots/20260316_120334_stage24_3_selected_relation_panel.png`
- plot: `plots/20260316_120334_stage24_3_stable_vs_contrast_panel.png`

## Commit recommendation
- recommended status: commit
- rationale: Phase IV now has a branch-local quasi-invariant read for the law-defined stable smeared sector using only the accepted observable core.
- branch effect: enables 24.4 extension-boundary testing from a compact law-and-quasi-invariant basis.
- caution: the selected read is branch-local and does not reopen braid or enlarge the observable ledger.
