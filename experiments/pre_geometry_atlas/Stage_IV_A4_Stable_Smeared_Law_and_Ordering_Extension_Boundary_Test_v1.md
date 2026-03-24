# Stage IV A4 Stable Smeared Law and Ordering Extension Boundary Test v1

## Stage purpose
This stage tests how far the fixed Stage 24.2 threshold law and the fixed Stage 24.3 inverse flow-survival ordering survive under bounded structural extension on the frozen clustered DK branch. The threshold remains fixed at `0.884308` throughout; no retuning is allowed.

## Authority note
This note is authoritative for the rerun stamped `20260316_121603`.

Superseded run:
- `20260316_121048` is superseded.
- reason: the earlier read only covered a narrower width-retreat manifold, while the authoritative rerun executes the full four-class extension lattice.
- interpretation change: yes. The authoritative run replaces the narrower anchor/retreat read with the full bounded extension-boundary map.

## Authoritative inputs
- `23.10 / III-C1`: json=`data/20260316_102709_stage23_10_dk_braid_mechanism_isolation_probe.json`, csv=`data/20260316_102709_stage23_10_dk_braid_mechanism_isolation_probe.csv`
- `23.12 / III-C3`: json=`data/20260316_105448_stage23_12_dk_effective_phase_diagram_closure.json`, csv=`data/20260316_105448_stage23_12_dk_effective_phase_diagram_closure.csv`
- `23.14 / III-C5`: json=`data/20260316_111201_stage23_14_phaseIII_final_consolidation_and_freeze.json`, csv=`data/20260316_111201_stage23_14_phaseIII_claims_ledger.csv`
- `24.1 / IV-A1`: json=`data/20260316_114317_stage24_1_smeared_sector_observable_ledger.json`, csv=`data/20260316_114317_stage24_1_smeared_sector_observable_ledger.csv`
- `24.2 / IV-A2`: json=`data/20260316_115346_stage24_2_smeared_sector_effective_law_extraction.json`, csv=`data/20260316_115346_stage24_2_smeared_sector_effective_law_model_ledger.csv`
- `24.3 / IV-A3`: json=`data/20260316_120334_stage24_3_stable_smeared_quasi_invariant_probe.json`, csv=`data/20260316_120334_stage24_3_stable_smeared_quasi_invariant_ledger.csv`

## Extension lattice definition
- fixed threshold law: `stable_closed_smeared = 1 iff flow_concentration_index <= 0.884308`
- carried ordering relation: `lower flow_concentration_index orders with longer topology_survival_time inside the stable smeared sector`
- extension classes: clustered texture perturbation, support anisotropy / skew, motif structure perturbation, corridor geometry perturbation
- extension levels: baseline, mild, moderate
- law-validity rule: exact five-phase family accuracy >= `1.0` without threshold retuning
- ordering-validity rule: stable subset size >= `3`, stable Spearman rho <= `-0.9`, and ordering gap >= `0.8`

## Law survival test
- `IVA4_E1_baseline`: class=`clustered_texture_perturbation`, level=`baseline`, law_valid=`True`, stable_subset_size=`5`, law_fit_score=`1.0000`
- `IVA4_E1_mild`: class=`clustered_texture_perturbation`, level=`mild`, law_valid=`True`, stable_subset_size=`4`, law_fit_score=`1.0000`
- `IVA4_E1_moderate`: class=`clustered_texture_perturbation`, level=`moderate`, law_valid=`True`, stable_subset_size=`2`, law_fit_score=`1.0000`
- `IVA4_E4_baseline`: class=`corridor_geometry_perturbation`, level=`baseline`, law_valid=`True`, stable_subset_size=`5`, law_fit_score=`1.0000`
- `IVA4_E4_mild`: class=`corridor_geometry_perturbation`, level=`mild`, law_valid=`True`, stable_subset_size=`5`, law_fit_score=`1.0000`
- `IVA4_E4_moderate`: class=`corridor_geometry_perturbation`, level=`moderate`, law_valid=`False`, stable_subset_size=`0`, law_fit_score=`0.8000`
- `IVA4_E3_baseline`: class=`motif_structure_perturbation`, level=`baseline`, law_valid=`True`, stable_subset_size=`5`, law_fit_score=`1.0000`
- `IVA4_E3_mild`: class=`motif_structure_perturbation`, level=`mild`, law_valid=`True`, stable_subset_size=`5`, law_fit_score=`1.0000`
- `IVA4_E3_moderate`: class=`motif_structure_perturbation`, level=`moderate`, law_valid=`False`, stable_subset_size=`0`, law_fit_score=`0.8000`
- `IVA4_E2_baseline`: class=`support_anisotropy_skew`, level=`baseline`, law_valid=`True`, stable_subset_size=`5`, law_fit_score=`1.0000`
- `IVA4_E2_mild`: class=`support_anisotropy_skew`, level=`mild`, law_valid=`True`, stable_subset_size=`4`, law_fit_score=`1.0000`
- `IVA4_E2_moderate`: class=`support_anisotropy_skew`, level=`moderate`, law_valid=`True`, stable_subset_size=`4`, law_fit_score=`1.0000`

## Ordering survival test
- `IVA4_E1_baseline`: ordering_applicable=`True`, ordering_valid=`True`, stable_rho=`-0.9747`, contrast_rho=`undefined`, ordering_gap=`0.9747`
- `IVA4_E1_mild`: ordering_applicable=`True`, ordering_valid=`True`, stable_rho=`-1.0000`, contrast_rho=`undefined`, ordering_gap=`1.0000`
- `IVA4_E1_moderate`: ordering_applicable=`False`, ordering_valid=`False`, stable_rho=`undefined`, contrast_rho=`undefined`, ordering_gap=`undefined`
- `IVA4_E4_baseline`: ordering_applicable=`True`, ordering_valid=`True`, stable_rho=`-0.9747`, contrast_rho=`undefined`, ordering_gap=`0.9747`
- `IVA4_E4_mild`: ordering_applicable=`True`, ordering_valid=`True`, stable_rho=`-0.9487`, contrast_rho=`undefined`, ordering_gap=`0.9487`
- `IVA4_E4_moderate`: ordering_applicable=`False`, ordering_valid=`False`, stable_rho=`undefined`, contrast_rho=`undefined`, ordering_gap=`undefined`
- `IVA4_E3_baseline`: ordering_applicable=`True`, ordering_valid=`True`, stable_rho=`-0.9747`, contrast_rho=`undefined`, ordering_gap=`0.9747`
- `IVA4_E3_mild`: ordering_applicable=`True`, ordering_valid=`False`, stable_rho=`-0.8660`, contrast_rho=`undefined`, ordering_gap=`0.8660`
- `IVA4_E3_moderate`: ordering_applicable=`False`, ordering_valid=`False`, stable_rho=`undefined`, contrast_rho=`undefined`, ordering_gap=`undefined`
- `IVA4_E2_baseline`: ordering_applicable=`True`, ordering_valid=`True`, stable_rho=`-0.9747`, contrast_rho=`undefined`, ordering_gap=`0.9747`
- `IVA4_E2_mild`: ordering_applicable=`True`, ordering_valid=`True`, stable_rho=`-0.9487`, contrast_rho=`undefined`, ordering_gap=`0.9487`
- `IVA4_E2_moderate`: ordering_applicable=`True`, ordering_valid=`True`, stable_rho=`-0.9487`, contrast_rho=`undefined`, ordering_gap=`0.9487`

## Extension status map
- `IVA4_E1_baseline`: status=`law_survives_ordering_survives`; notes=Both the fixed threshold law and the preserved inverse ordering survive.
- `IVA4_E1_mild`: status=`law_survives_ordering_survives`; notes=Both the fixed threshold law and the preserved inverse ordering survive.
- `IVA4_E1_moderate`: status=`law_survives_ordering_degrades`; notes=The threshold law survives, but the preserved inverse ordering no longer satisfies the frozen 24.3 criteria.
- `IVA4_E4_baseline`: status=`law_survives_ordering_survives`; notes=Both the fixed threshold law and the preserved inverse ordering survive.
- `IVA4_E4_mild`: status=`law_survives_ordering_survives`; notes=Both the fixed threshold law and the preserved inverse ordering survive.
- `IVA4_E4_moderate`: status=`law_fails_ordering_not_applicable`; notes=The fixed Stage 24.2 threshold law no longer classifies the stable smeared subset exactly on this five-phase family.
- `IVA4_E3_baseline`: status=`law_survives_ordering_survives`; notes=Both the fixed threshold law and the preserved inverse ordering survive.
- `IVA4_E3_mild`: status=`law_survives_ordering_degrades`; notes=The threshold law survives, but the preserved inverse ordering no longer satisfies the frozen 24.3 criteria.
- `IVA4_E3_moderate`: status=`law_fails_ordering_not_applicable`; notes=The fixed Stage 24.2 threshold law no longer classifies the stable smeared subset exactly on this five-phase family.
- `IVA4_E2_baseline`: status=`law_survives_ordering_survives`; notes=Both the fixed threshold law and the preserved inverse ordering survive.
- `IVA4_E2_mild`: status=`law_survives_ordering_survives`; notes=Both the fixed threshold law and the preserved inverse ordering survive.
- `IVA4_E2_moderate`: status=`law_survives_ordering_survives`; notes=Both the fixed threshold law and the preserved inverse ordering survive.

## First extension boundary
- `clustered_texture_perturbation`: first boundary at `IVA4_E1_moderate` (`moderate`) with status `law_survives_ordering_degrades`
- `support_anisotropy_skew`: no failure boundary detected on the tested levels
- `motif_structure_perturbation`: first boundary at `IVA4_E3_mild` (`mild`) with status `law_survives_ordering_degrades`
- `corridor_geometry_perturbation`: first boundary at `IVA4_E4_moderate` (`moderate`) with status `law_fails_ordering_not_applicable`

## Minimal interpretation
On the frozen clustered DK branch, the Phase IV law-and-ordering package survives only within a bounded extension radius. The fixed 24.2 threshold law and the 24.3 inverse flow-survival ordering can be tracked separately under local structural extension, allowing a clean identification of where law survival persists, where ordering degrades, and where the branch breaks down.

## Output artifact list
- stamped JSON summary: `data/20260316_121603_stage24_4_stable_smeared_law_and_ordering_extension_boundary_test.json`
- stamped CSV extension ledger: `data/20260316_121603_stage24_4_stable_smeared_extension_boundary_ledger.csv`
- plot: `plots/20260316_121603_stage24_4_extension_status_matrix.png`
- plot: `plots/20260316_121603_stage24_4_law_survival_panel.png`
- plot: `plots/20260316_121603_stage24_4_ordering_survival_panel.png`
- plot: `plots/20260316_121603_stage24_4_boundary_transition_map.png`

## Commit recommendation
- recommended status: commit
- rationale: Phase IV now has a boundary-qualified law-and-ordering package for the stable smeared sector.
- branch effect: clarifies the local extension radius of the surviving transport structure.
- caution: survival remains branch-local and does not imply a family-wide transport law.
