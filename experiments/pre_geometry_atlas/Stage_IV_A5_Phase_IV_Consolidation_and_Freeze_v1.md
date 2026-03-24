# Stage IV A5 Phase IV Consolidation and Freeze v1

## Stage purpose
This stage freezes the current Phase IV law-and-ordering branch on the frozen clustered DK sector. It is a reduction-only consolidation pass over Stages 24.1 through 24.4, with the 24.4 rerun stamped `20260316_121603` treated as authoritative.

## Authority note
This note is authoritative for the rerun stamped `20260316_122247`.

Superseded run:
- `20260316_122224` is superseded.
- reason: the initial freeze attempt used an overly blunt portable-core contradiction check.
- interpretation change: no. The rerun corrects the freeze-validation logic but preserves the same scientific Phase IV read.

## Inputs consolidated
- `23.14 / III-C5`: json=`data/20260316_111201_stage23_14_phaseIII_final_consolidation_and_freeze.json`, csv=`data/20260316_111201_stage23_14_phaseIII_claims_ledger.csv`
- `24.1 / IV-A1`: json=`data/20260316_114317_stage24_1_smeared_sector_observable_ledger.json`, csv=`data/20260316_114317_stage24_1_smeared_sector_observable_ledger.csv`
- `24.2 / IV-A2`: json=`data/20260316_115346_stage24_2_smeared_sector_effective_law_extraction.json`, csv=`data/20260316_115346_stage24_2_smeared_sector_effective_law_model_ledger.csv`
- `24.3 / IV-A3`: json=`data/20260316_120334_stage24_3_stable_smeared_quasi_invariant_probe.json`, csv=`data/20260316_120334_stage24_3_stable_smeared_quasi_invariant_ledger.csv`
- `24.4 / IV-A4`: json=`data/20260316_121603_stage24_4_stable_smeared_law_and_ordering_extension_boundary_test.json`, csv=`data/20260316_121603_stage24_4_stable_smeared_extension_boundary_ledger.csv`

## Consistency checks
- `stage24_1_expected_primary_observables` = `True`
- `stage24_1_rejected_family_intact` = `True`
- `stage24_2_selector_law_exact` = `True`
- `stage24_3_selected_q3` = `True`
- `stage24_3_contrast_separated_spearman_signs` = `True`
- `stage24_4_status_counts_match` = `True`
- `stage24_4_boundaries_match` = `True`
- `portable_claims_noncontradictory` = `True`
- `phaseIII_negatives_preserved` = `True`

## Portable but bounded positive core
- `P4_POS_01`: The smeared sector is carried by the accepted 24.1 observable core and not by rejected braid- or asymmetry-based quantities. [survives_portably]
- `P4_POS_02`: The stable closed smeared sector is selected by the fixed threshold law: stable_closed_smeared = 1 iff flow_concentration_index <= 0.884308. [survives_portably]
- `P4_POS_03`: The fixed threshold selector remains valid across most tested extension conditions without retuning. [survives_portably]
- `P4_POS_04`: Ordering boundaries appear before the first selector-failure boundary on the tested extension lattice. [survives_portably]
- `P4_POS_05`: phase_corridor_position and phase_corridor_width remain external conditioning coordinates rather than Phase IV state variables. [survives_portably]

## Locally portable internal structure
- `P4_LOC_01`: Inside the law-defined stable smeared sector, lower flow_concentration_index orders with longer topology_survival_time. [survives_locally]
- `P4_LOC_02`: The accepted inverse flow-survival ordering is not a restatement of the selector threshold law. [survives_locally]
- `P4_LOC_03`: The stable-set inverse ordering differs sharply from contrast behavior on the frozen branch. [survives_locally]
- `P4_LOC_04`: The inverse flow-survival ordering has bounded local portability but does not survive uniformly across the tested extension lattice. [survives_locally]
- `P4_LOC_05`: Support anisotropy / skew stayed inside the tested survival radius on the current extension lattice. [bounded_support]

## Failed / excluded claims
- `P4_NEG_01`: The inverse flow-survival ordering does not survive uniformly across the tested extension lattice. [fails]
- `P4_NEG_02`: Motif structure perturbation reaches the earliest tested ordering boundary at IVA4_E3_mild. [fails]
- `P4_NEG_03`: Clustered texture perturbation reaches ordering degradation at IVA4_E1_moderate. [fails]
- `P4_NEG_04`: Moderate corridor geometry perturbation reaches the first tested selector-failure boundary at IVA4_E4_moderate. [fails]
- `P4_NEG_05`: No braid-labeled quantity re-enters the Phase IV core. [excluded]
- `P4_NEG_06`: The rejected observable family remains excluded from Phase IV core structure. [excluded]
- `P4_NEG_07`: The current law-and-ordering package does not justify a family-wide transport law claim. [excluded]
- `P4_NEG_08`: Nothing in Phase IV reopens family-wide braid mechanism, stable braid phase, or localized encounter claims from Phase III. [excluded]

## Final closure statement
On the frozen clustered DK branch, Phase IV identifies a compact threshold selector for the stable closed smeared sector, flow_concentration_index <= 0.884308, together with a non-trivial inverse flow-survival ordering relation that survives only within a bounded local extension radius. Under controlled structural extension, ordering degradation appears before selector failure in some classes, while moderate corridor-geometry perturbation provides the first tested boundary at which the fixed selector itself fails. The surviving Phase IV content is therefore a bounded local transport package with separable selector and ordering boundaries, not a family-wide transport structure.

## Frozen claims ledger summary
- portable core claims = `5`
- local structure claims = `5`
- failed / excluded claims = `8`
- extension status counts = `{'law_survives_ordering_survives': 8, 'law_survives_ordering_degrades': 2, 'law_fails_ordering_not_applicable': 2}`

## Output artifact list
- stamped JSON summary: `data/20260316_122247_stage24_5_phaseIV_consolidation_and_freeze.json`
- stamped CSV claims ledger: `data/20260316_122247_stage24_5_phaseIV_claims_ledger.csv`
- plot: `plots/20260316_122247_stage24_5_phaseIV_portable_core_table.png`
- plot: `plots/20260316_122247_stage24_5_phaseIV_local_structure_table.png`
- plot: `plots/20260316_122247_stage24_5_phaseIV_boundary_and_negative_table.png`
- plot: `plots/20260316_122247_stage24_5_phaseIV_status_map_panel.png`

## Commit recommendation
- recommended status: commit
- rationale: Phase IV selector, local ordering structure, and extension boundaries are now internally consolidated.
- branch effect: closes the current Phase IV law-and-ordering arc on the frozen clustered DK branch.
- caution: reopening portability or internal-structure claims requires a new branch with explicitly changed extension manifold, operator class, or structural prior.
