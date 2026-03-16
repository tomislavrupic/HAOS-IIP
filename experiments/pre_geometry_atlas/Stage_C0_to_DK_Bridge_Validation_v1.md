# Stage C0-to-DK Bridge Validation v1

Timestamped JSON: `data/20260316_084155_stage_c0_dk_bridge_validation.json`
Timestamped CSV: `data/20260316_084155_stage_c0_dk_bridge_validation.csv`

Purpose: test a genuinely executable C0-to-DK bridge candidate using only signed combinatorial incidence, deterministic harmonic-address dressing, and bridge-specific dynamics on the derived graph complex.

## Algebraic checks
- `clustered_composite_anchor`: `||d1 d0|| = 0.000e+00`, `||D^2 - Delta|| = 0.000e+00`, `||d0^T d0 - L0|| = 0.000e+00`, `triangles = 6`
- `counter_propagating_corridor`: `||d1 d0|| = 0.000e+00`, `||D^2 - Delta|| = 0.000e+00`, `||d0^T d0 - L0|| = 0.000e+00`, `triangles = 6`
- `phase_ordered_symmetric_triad`: `||d1 d0|| = 0.000e+00`, `||D^2 - Delta|| = 0.000e+00`, `||d0^T d0 - L0|| = 0.000e+00`, `triangles = 3`

## Detuning scan
- `clustered_composite_anchor`: first class change at `delta = 0.250`; trace = `0.000:braid_like_exchange, 0.125:braid_like_exchange, 0.250:transfer_smeared, 0.375:transfer_smeared, 0.500:transfer_smeared, 0.625:transfer_smeared, 0.750:transfer_smeared, 0.875:transfer_smeared, 1.000:transfer_smeared`
- `counter_propagating_corridor`: first class change at `delta = 0.125`; trace = `0.000:transfer_smeared, 0.125:unresolved_mixed, 0.250:unresolved_mixed, 0.375:unresolved_mixed, 0.500:unresolved_mixed, 0.625:unresolved_mixed, 0.750:transfer_smeared, 0.875:unresolved_mixed, 1.000:unresolved_mixed`
- `phase_ordered_symmetric_triad`: first class change at `delta = none`; trace = `0.000:transfer_smeared, 0.125:transfer_smeared, 0.250:transfer_smeared, 0.375:transfer_smeared, 0.500:transfer_smeared, 0.625:transfer_smeared, 0.750:transfer_smeared, 0.875:transfer_smeared, 1.000:transfer_smeared`

## Focused protection panel
- `balanced_baseline`: topology=`braid_like_exchange`, flow=`0.2908`, coherence=`0.4969`, selectivity=`1.0000`, loop score=`0.2185`
- `phase_edge_probe`: topology=`transfer_smeared`, flow=`0.4116`, coherence=`0.4850`, selectivity=`0.9235`, loop score=`0.1136`
- `width_asymmetry_probe`: topology=`transfer_smeared`, flow=`0.2258`, coherence=`0.5477`, selectivity=`1.0000`, loop score=`0.2652`
- `degree_skew_probe`: topology=`braid_like_exchange`, flow=`0.2834`, coherence=`0.5846`, selectivity=`1.0000`, loop score=`0.2120`
- `zero_form_enhanced_probe`: topology=`transfer_smeared`, flow=`0.3050`, coherence=`0.4527`, selectivity=`1.0000`, loop score=`0.2056`

## Obligation status
- `O_C0_ORIENT`: open; the tested signed complex closes exactly, but the orientation source is still a mismatch surrogate rather than the protected sequencing ledger itself.
- `O_C0_EVOLUTION`: open; the bridge uses an explicit derived `D_H` evolution, but equivalence to the protected sequencing rule remains unproved.
- `O_C0_GRADING`: open; harmonic address remains a dressing, not yet an intrinsic grade source.
- `O_C0_RECOVERY`: partially closed; the clustered family opens a real braid-to-smear bridge corridor and the focused clustered protection panel matches `1.00` of the Stage-23-style qualitative pattern, but the first class change occurs at `delta = 0.25` and the corridor/triad families do not recover a broad braid sector.

## Canonical statement
A combinatorial Dirac-Kahler candidate can now be constructed from C0 primitives, and on the tested representative graph complexes its signed chain identities close exactly. The present run therefore establishes a real bridge scaffold rather than a finished derivation: the clustered family shows a measurable harmonic detuning boundary and a focused protection-panel recovery, but the orientation source, evolution-law equivalence, intrinsic grading source, and family-wide protection recovery remain conditional.

## Plots
- `plots/20260316_084155_stage_c0_dk_bridge_validation_detuning_panel.png`
- `plots/20260316_084155_stage_c0_dk_bridge_validation_algebraic_panel.png`
- `plots/20260316_084155_stage_c0_dk_bridge_validation_protection_panel.png`
