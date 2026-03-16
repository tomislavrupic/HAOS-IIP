# Stage 12 Intrinsic Grading and Structural Irreversibility v1

Timestamped JSON: `data/20260316_093434_stage_12_intrinsic_grading_and_structural_irreversibility_probe.json`
Timestamped CSV: `data/20260316_093434_stage_12_intrinsic_grading_and_structural_irreversibility_probe.csv`

Global read: `negative_closure`

## Obligation diagnostics
- `O_C0_ORIENT`: positive-family count = `0`; success = `False`
- `O_C0_EVOLUTION`: median refinement ratio = `1.0381`; success = `True`
- `Intrinsic grading`: entropy = `0.2570`, separation = `0.2019`, monotonicity = `-0.3616`; success = `False`
- `Family-wide recovery`: ordered settings fraction = `0.0000`; success = `False`
- `Structural irreversibility`: bidirectional positive fraction = `0.3333`; success = `True`

## Family medians
- `clustered_braid_seed`: alignment=`0.1161`, protected_fraction=`1.0000`, irreversibility=`0.0426`, entropy=`0.2458`
- `corridor_transport_seed`: alignment=`0.0339`, protected_fraction=`0.5541`, irreversibility=`1.0503`, entropy=`0.3501`
- `triad_competition_seed`: alignment=`0.1081`, protected_fraction=`0.9965`, irreversibility=`0.0269`, entropy=`0.2305`

## Per-run summary
- `S12_clustered_braid_seed_forward_only_base`: family=`Clustered protected seed`, protocol=`Forward only`, resolution=`Base resolution`, forward=`transfer_smeared`, return=`transfer_smeared`, final=`transfer_smeared`, protected=`1.0000`, irreversibility=`0.0000`, entropy=`0.0000`
- `S12_clustered_braid_seed_forward_only_refined`: family=`Clustered protected seed`, protocol=`Forward only`, resolution=`Refinement doubling`, forward=`transfer_smeared`, return=`transfer_smeared`, final=`transfer_smeared`, protected=`1.0000`, irreversibility=`0.0000`, entropy=`0.0507`
- `S12_clustered_braid_seed_forward_reverse_restore_base`: family=`Clustered protected seed`, protocol=`Forward, reverse, restore`, resolution=`Base resolution`, forward=`transfer_smeared`, return=`transfer_smeared`, final=`transfer_smeared`, protected=`1.0000`, irreversibility=`0.0426`, entropy=`0.2832`
- `S12_clustered_braid_seed_forward_reverse_restore_refined`: family=`Clustered protected seed`, protocol=`Forward, reverse, restore`, resolution=`Refinement doubling`, forward=`transfer_smeared`, return=`transfer_smeared`, final=`transfer_smeared`, protected=`0.9911`, irreversibility=`0.0456`, entropy=`0.2630`
- `S12_clustered_braid_seed_forward_reverse_reperturb_base`: family=`Clustered protected seed`, protocol=`Forward, reverse, re-perturb`, resolution=`Base resolution`, forward=`transfer_smeared`, return=`transfer_smeared`, final=`transfer_smeared`, protected=`1.0000`, irreversibility=`0.0426`, entropy=`0.2648`
- `S12_clustered_braid_seed_forward_reverse_reperturb_refined`: family=`Clustered protected seed`, protocol=`Forward, reverse, re-perturb`, resolution=`Refinement doubling`, forward=`transfer_smeared`, return=`transfer_smeared`, final=`transfer_smeared`, protected=`0.9861`, irreversibility=`0.0456`, entropy=`0.2287`
- `S12_corridor_transport_seed_forward_only_base`: family=`Corridor weak-protection seed`, protocol=`Forward only`, resolution=`Base resolution`, forward=`localized_capture`, return=`transfer_smeared`, final=`transfer_smeared`, protected=`0.4167`, irreversibility=`1.0000`, entropy=`0.3261`
- `S12_corridor_transport_seed_forward_only_refined`: family=`Corridor weak-protection seed`, protocol=`Forward only`, resolution=`Refinement doubling`, forward=`localized_capture`, return=`transfer_smeared`, final=`transfer_smeared`, protected=`0.4792`, irreversibility=`1.0000`, entropy=`0.3342`
- `S12_corridor_transport_seed_forward_reverse_restore_base`: family=`Corridor weak-protection seed`, protocol=`Forward, reverse, restore`, resolution=`Base resolution`, forward=`localized_capture`, return=`transfer_smeared`, final=`localized_capture`, protected=`0.5179`, irreversibility=`1.0503`, entropy=`0.4287`
- `S12_corridor_transport_seed_forward_reverse_restore_refined`: family=`Corridor weak-protection seed`, protocol=`Forward, reverse, restore`, resolution=`Refinement doubling`, forward=`localized_capture`, return=`transfer_smeared`, final=`localized_capture`, protected=`0.6339`, irreversibility=`1.0602`, entropy=`0.3436`
- `S12_corridor_transport_seed_forward_reverse_reperturb_base`: family=`Corridor weak-protection seed`, protocol=`Forward, reverse, re-perturb`, resolution=`Base resolution`, forward=`localized_capture`, return=`transfer_smeared`, final=`transfer_smeared`, protected=`0.6389`, irreversibility=`1.0503`, entropy=`0.3997`
- `S12_corridor_transport_seed_forward_reverse_reperturb_refined`: family=`Corridor weak-protection seed`, protocol=`Forward, reverse, re-perturb`, resolution=`Refinement doubling`, forward=`localized_capture`, return=`transfer_smeared`, final=`transfer_smeared`, protected=`0.5903`, irreversibility=`1.0602`, entropy=`0.3566`
- `S12_triad_competition_seed_forward_only_base`: family=`Triad competitive seed`, protocol=`Forward only`, resolution=`Base resolution`, forward=`transfer_smeared`, return=`transfer_smeared`, final=`transfer_smeared`, protected=`1.0000`, irreversibility=`0.0000`, entropy=`0.1725`
- `S12_triad_competition_seed_forward_only_refined`: family=`Triad competitive seed`, protocol=`Forward only`, resolution=`Refinement doubling`, forward=`transfer_smeared`, return=`transfer_smeared`, final=`transfer_smeared`, protected=`1.0000`, irreversibility=`0.0000`, entropy=`0.2505`
- `S12_triad_competition_seed_forward_reverse_restore_base`: family=`Triad competitive seed`, protocol=`Forward, reverse, restore`, resolution=`Base resolution`, forward=`transfer_smeared`, return=`transfer_smeared`, final=`transfer_smeared`, protected=`1.0000`, irreversibility=`0.0304`, entropy=`0.1892`
- `S12_triad_competition_seed_forward_reverse_restore_refined`: family=`Triad competitive seed`, protocol=`Forward, reverse, restore`, resolution=`Refinement doubling`, forward=`transfer_smeared`, return=`transfer_smeared`, final=`transfer_smeared`, protected=`0.9911`, irreversibility=`0.0269`, entropy=`0.2429`
- `S12_triad_competition_seed_forward_reverse_reperturb_base`: family=`Triad competitive seed`, protocol=`Forward, reverse, re-perturb`, resolution=`Base resolution`, forward=`transfer_smeared`, return=`transfer_smeared`, final=`transfer_smeared`, protected=`0.9722`, irreversibility=`0.0304`, entropy=`0.2732`
- `S12_triad_competition_seed_forward_reverse_reperturb_refined`: family=`Triad competitive seed`, protocol=`Forward, reverse, re-perturb`, resolution=`Refinement doubling`, forward=`transfer_smeared`, return=`transfer_smeared`, final=`transfer_smeared`, protected=`0.9931`, irreversibility=`0.0269`, entropy=`0.2181`

## Plots
- `plots/20260316_093434_S12_clustered_braid_seed_forward_only_base_trace.png`
- `plots/20260316_093434_S12_clustered_braid_seed_forward_only_refined_trace.png`
- `plots/20260316_093434_S12_clustered_braid_seed_forward_reverse_restore_base_trace.png`
- `plots/20260316_093434_S12_clustered_braid_seed_forward_reverse_restore_refined_trace.png`
- `plots/20260316_093434_S12_clustered_braid_seed_forward_reverse_reperturb_base_trace.png`
- `plots/20260316_093434_S12_clustered_braid_seed_forward_reverse_reperturb_refined_trace.png`
- `plots/20260316_093434_S12_corridor_transport_seed_forward_only_base_trace.png`
- `plots/20260316_093434_S12_corridor_transport_seed_forward_only_refined_trace.png`
- `plots/20260316_093434_S12_corridor_transport_seed_forward_reverse_restore_base_trace.png`
- `plots/20260316_093434_S12_corridor_transport_seed_forward_reverse_restore_refined_trace.png`
- `plots/20260316_093434_S12_corridor_transport_seed_forward_reverse_reperturb_base_trace.png`
- `plots/20260316_093434_S12_corridor_transport_seed_forward_reverse_reperturb_refined_trace.png`
- `plots/20260316_093434_S12_triad_competition_seed_forward_only_base_trace.png`
- `plots/20260316_093434_S12_triad_competition_seed_forward_only_refined_trace.png`
- `plots/20260316_093434_S12_triad_competition_seed_forward_reverse_restore_base_trace.png`
- `plots/20260316_093434_S12_triad_competition_seed_forward_reverse_restore_refined_trace.png`
- `plots/20260316_093434_S12_triad_competition_seed_forward_reverse_reperturb_base_trace.png`
- `plots/20260316_093434_S12_triad_competition_seed_forward_reverse_reperturb_refined_trace.png`
- `plots/20260316_093434_stage_12_entropy_panel.png`
- `plots/20260316_093434_stage_12_irreversibility_panel.png`
- `plots/20260316_093434_stage_12_refinement_panel.png`
