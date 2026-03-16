# Stage III-C4 Mechanism Isolation Probe v1

This mechanism probe stays inside the frozen HAOS (Harmonic Address Operating System) / IIP (Interaction Invariance Physics) Dirac-Kahler clustered collision sector.

Timestamped JSON: `data/20260316_110320_stage23_13_dk_mechanism_isolation_probe.json`
Timestamped CSV: `data/20260316_110320_stage23_13_dk_mechanism_isolation_probe.csv`

overall_mechanism_label = `phase_corridor_primary`

## Assumption note
The grade-transfer branch is implemented as bounded blockwise grade-fraction anchoring, because the frozen DK branch does not expose a literal transfer-off switch without changing operator class. The operator-coupling branch scales only the `0 <-> 1` DK blocks and leaves the rest of the propagator intact.

## Branch outcomes
- `Grade-transfer suppression`: label=`no_clear_mechanism_signal`, control=`braid_like_exchange`, weak=`braid_like_exchange`, moderate=`braid_like_exchange`, mean_sensitivity=`0.000`
- `Phase-corridor flattening`: label=`phase_corridor_primary`, control=`braid_like_exchange`, weak=`transfer_smeared`, moderate=`transfer_smeared`, mean_sensitivity=`0.552`
- `Geometry-overlap weakening`: label=`coupled_mechanism_no_single_primary`, control=`braid_like_exchange`, weak=`braid_like_exchange`, moderate=`transfer_smeared`, mean_sensitivity=`0.325`
- `Operator cross-coupling suppression`: label=`coupled_mechanism_no_single_primary`, control=`braid_like_exchange`, weak=`braid_like_exchange`, moderate=`transfer_smeared`, mean_sensitivity=`0.399`

## Per-run summary
- `S23_13_A_control`: topology=`braid_like_exchange`, braid_survival=`0.376`, coherence=`0.483`, flow=`0.889`, sensitivity=`0.000`
- `S23_13_A_weak`: topology=`braid_like_exchange`, braid_survival=`0.626`, coherence=`0.564`, flow=`0.889`, sensitivity=`0.000`
- `S23_13_A_moderate`: topology=`braid_like_exchange`, braid_survival=`0.501`, coherence=`0.563`, flow=`0.886`, sensitivity=`0.000`
- `S23_13_B_control`: topology=`braid_like_exchange`, braid_survival=`0.376`, coherence=`0.483`, flow=`0.889`, sensitivity=`0.000`
- `S23_13_B_weak`: topology=`transfer_smeared`, braid_survival=`0.626`, coherence=`0.478`, flow=`0.891`, sensitivity=`0.552`
- `S23_13_B_moderate`: topology=`transfer_smeared`, braid_survival=`0.626`, coherence=`0.475`, flow=`0.893`, sensitivity=`0.553`
- `S23_13_C_control`: topology=`braid_like_exchange`, braid_survival=`0.376`, coherence=`0.483`, flow=`0.889`, sensitivity=`0.000`
- `S23_13_C_weak`: topology=`braid_like_exchange`, braid_survival=`0.376`, coherence=`0.484`, flow=`0.887`, sensitivity=`0.000`
- `S23_13_C_moderate`: topology=`transfer_smeared`, braid_survival=`0.251`, coherence=`0.485`, flow=`0.884`, sensitivity=`0.650`
- `S23_13_D_control`: topology=`braid_like_exchange`, braid_survival=`0.376`, coherence=`0.483`, flow=`0.889`, sensitivity=`0.000`
- `S23_13_D_weak`: topology=`braid_like_exchange`, braid_survival=`0.501`, coherence=`0.420`, flow=`0.892`, sensitivity=`0.020`
- `S23_13_D_moderate`: topology=`transfer_smeared`, braid_survival=`0.125`, coherence=`0.396`, flow=`0.882`, sensitivity=`0.777`

## Mechanism conclusion
The strongest discriminative collapse comes from phase-corridor flattening: both weak and moderate phase flattening move the clustered control out of the braid class while the control remains braid-like and the run stays bounded rather than dissolving into noise.
- Conclusion statement: `braid-like exchange in the DK sector is primarily controlled by phase corridor, with geometry overlap and operator cross-coupling acting as secondary support`

## Plots
- `plots/20260316_110320_S23_13_A_control_trace_panel.png`
- `plots/20260316_110320_S23_13_A_weak_trace_panel.png`
- `plots/20260316_110320_S23_13_A_moderate_trace_panel.png`
- `plots/20260316_110320_S23_13_B_control_trace_panel.png`
- `plots/20260316_110320_S23_13_B_weak_trace_panel.png`
- `plots/20260316_110320_S23_13_B_moderate_trace_panel.png`
- `plots/20260316_110320_S23_13_C_control_trace_panel.png`
- `plots/20260316_110320_S23_13_C_weak_trace_panel.png`
- `plots/20260316_110320_S23_13_C_moderate_trace_panel.png`
- `plots/20260316_110320_S23_13_D_control_trace_panel.png`
- `plots/20260316_110320_S23_13_D_weak_trace_panel.png`
- `plots/20260316_110320_S23_13_D_moderate_trace_panel.png`
- `plots/20260316_110320_stage23_13_mechanism_branch_comparison_matrix.png`
- `plots/20260316_110320_stage23_13_braid_survival_vs_ablation_panel.png`
- `plots/20260316_110320_stage23_13_coherence_collapse_panel.png`
