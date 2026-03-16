# Stage III-C3 Effective Phase Diagram Closure v1

This closure pass stays inside the frozen HAOS (Harmonic Address Operating System) / IIP (Interaction Invariance Physics) clustered Dirac-Kahler collision sector.

Timestamped JSON: `data/20260316_105448_stage23_12_dk_effective_phase_diagram_closure.json`
Timestamped CSV: `data/20260316_105448_stage23_12_dk_effective_phase_diagram_closure.csv`

phase_diagram_closed = `TRUE`
closure_classification = `smeared_dominant_closed_phase_diagram`

## Closure note
The source brief specified a three-axis manifold but a 15-run primary grid. This implementation therefore uses a `5 x 3` phase-width lattice and treats `beta in {0.01, 0.02}` as an internal weak-coupling stability probe inside each primary cell rather than as an external 45-point expansion.

## Summary metrics
- stable connected region present: `True`
- finite boundary flag: `True`
- bidirectional preservation flag: `True`
- boundary segment count: `5`
- stable phase counts: `{'transient_mixed_phase': 10, 'smeared_transfer_phase': 5}`
- dominant topology counts: `{'braid_like_exchange': 6, 'transfer_smeared': 9}`

## Per-cell summary
- `S23_12_p0375_w100`: phase=`0.375`, width=`1.00`, topology=`braid_like_exchange`, phase_label=`transient_mixed_phase`, R=`0`, B=`1`, D=`0`
- `S23_12_p0475_w100`: phase=`0.475`, width=`1.00`, topology=`braid_like_exchange`, phase_label=`transient_mixed_phase`, R=`0`, B=`1`, D=`0`
- `S23_12_p0500_w100`: phase=`0.500`, width=`1.00`, topology=`transfer_smeared`, phase_label=`transient_mixed_phase`, R=`0`, B=`0`, D=`0`
- `S23_12_p0550_w100`: phase=`0.550`, width=`1.00`, topology=`transfer_smeared`, phase_label=`transient_mixed_phase`, R=`0`, B=`0`, D=`0`
- `S23_12_p0575_w100`: phase=`0.575`, width=`1.00`, topology=`transfer_smeared`, phase_label=`transient_mixed_phase`, R=`0`, B=`0`, D=`0`
- `S23_12_p0375_w115`: phase=`0.375`, width=`1.15`, topology=`braid_like_exchange`, phase_label=`transient_mixed_phase`, R=`0`, B=`1`, D=`0`
- `S23_12_p0475_w115`: phase=`0.475`, width=`1.15`, topology=`braid_like_exchange`, phase_label=`transient_mixed_phase`, R=`0`, B=`0`, D=`0`
- `S23_12_p0500_w115`: phase=`0.500`, width=`1.15`, topology=`braid_like_exchange`, phase_label=`transient_mixed_phase`, R=`0`, B=`0`, D=`0`
- `S23_12_p0550_w115`: phase=`0.550`, width=`1.15`, topology=`braid_like_exchange`, phase_label=`transient_mixed_phase`, R=`0`, B=`0`, D=`0`
- `S23_12_p0575_w115`: phase=`0.575`, width=`1.15`, topology=`transfer_smeared`, phase_label=`transient_mixed_phase`, R=`0`, B=`1`, D=`0`
- `S23_12_p0375_w135`: phase=`0.375`, width=`1.35`, topology=`transfer_smeared`, phase_label=`smeared_transfer_phase`, R=`1`, B=`1`, D=`1`
- `S23_12_p0475_w135`: phase=`0.475`, width=`1.35`, topology=`transfer_smeared`, phase_label=`smeared_transfer_phase`, R=`1`, B=`1`, D=`1`
- `S23_12_p0500_w135`: phase=`0.500`, width=`1.35`, topology=`transfer_smeared`, phase_label=`smeared_transfer_phase`, R=`1`, B=`1`, D=`1`
- `S23_12_p0550_w135`: phase=`0.550`, width=`1.35`, topology=`transfer_smeared`, phase_label=`smeared_transfer_phase`, R=`1`, B=`1`, D=`1`
- `S23_12_p0575_w135`: phase=`0.575`, width=`1.35`, topology=`transfer_smeared`, phase_label=`smeared_transfer_phase`, R=`1`, B=`1`, D=`1`

## Minimal regime boundary summary
- smeared_transfer_phase: size=5, phase in [0.375, 0.575], width ratio in [1.35, 1.35]

## Interpretation
The clustered DK sector admits a closed effective phase description on this manifold, but the closure is one-sided: the stable connected region is a smeared-transfer phase, while the nominal braid corridor collapses into transient mixed cells once refinement, weak-coupling consistency, and reverse evolution are enforced together. Phase III therefore closes this branch as a sector-local phenomenology map rather than a general intrinsic braid mechanism.
No localized encounter phase closed on this lattice.
No stable braid phase survived the full closure criteria on this lattice.

## Plots
- `plots/20260316_105448_stage23_12_effective_phase_map.png`
- `plots/20260316_105448_stage23_12_topology_survival_heatmap.png`
- `plots/20260316_105448_stage23_12_refinement_stability_table.png`
