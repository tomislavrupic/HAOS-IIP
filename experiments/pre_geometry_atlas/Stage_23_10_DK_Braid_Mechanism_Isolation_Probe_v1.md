# Stage 23.10 DK Braid Mechanism Isolation Probe v1

Timestamped JSON: `data/20260316_102709_stage23_10_dk_braid_mechanism_isolation_probe.json`
Timestamped CSV: `data/20260316_102709_stage23_10_dk_braid_mechanism_isolation_probe.csv`

Mechanism classification: `clustered_texture_artefact`

## Per-run readout
- `S23_10_clustered_baseline`: base=`braid_like_exchange`, refined=`braid_like_exchange`, repeatability=`0.75`, stability=`1.33`
- `S23_10_clustered_degraded`: base=`transfer_smeared`, refined=`transfer_smeared`, repeatability=`1.00`, stability=`1.33`
- `S23_10_clustered_topology_broken`: base=`braid_like_exchange`, refined=`braid_like_exchange`, repeatability=`1.00`, stability=`1.00`
- `S23_10_corridor_baseline`: base=`transfer_smeared`, refined=`transfer_smeared`, repeatability=`1.00`, stability=`1.00`
- `S23_10_corridor_motif_injected`: base=`transfer_smeared`, refined=`transfer_smeared`, repeatability=`1.00`, stability=`1.00`
- `S23_10_corridor_anisotropic_kernel`: base=`transfer_smeared`, refined=`transfer_smeared`, repeatability=`1.00`, stability=`1.00`
- `S23_10_triad_baseline`: base=`transfer_smeared`, refined=`transfer_smeared`, repeatability=`1.00`, stability=`1.00`
- `S23_10_triad_degree_skewed`: base=`transfer_smeared`, refined=`transfer_smeared`, repeatability=`1.00`, stability=`1.00`
- `S23_10_distributed_random_seed`: base=`transfer_smeared`, refined=`transfer_smeared`, repeatability=`1.00`, stability=`1.00`

## Decision summary
- robust non-clustered braid hits: `0`
- non-clustered braid hits (raw): `0`
- clustered braid hits: `2`
- clustered degraded collapse: `True`

## Interpretation
The braid window remains clustered-texture dependent in this pass: non-clustered families do not recover a robust braid sector, and the clustered family either stays isolated or collapses under mild degradation. Phase III therefore closes this line as a seed-dependent texture rather than a general exchange law.

## Assumption note
Because the frozen DK lattice here is regular, the requested degree-spectrum perturbation was represented by a mild localized graded-support skew surrogate rather than a literal graph-degree modification. The kernel class, grading rule, timestep policy, and normalization were otherwise held fixed.

## Plots
- `plots/20260316_102709_S23_10_clustered_baseline_n12_trajectory.png`
- `plots/20260316_102709_S23_10_clustered_baseline_n24_trajectory.png`
- `plots/20260316_102709_S23_10_clustered_degraded_n12_trajectory.png`
- `plots/20260316_102709_S23_10_clustered_degraded_n24_trajectory.png`
- `plots/20260316_102709_S23_10_clustered_topology_broken_n12_trajectory.png`
- `plots/20260316_102709_S23_10_clustered_topology_broken_n24_trajectory.png`
- `plots/20260316_102709_S23_10_corridor_baseline_n12_trajectory.png`
- `plots/20260316_102709_S23_10_corridor_baseline_n24_trajectory.png`
- `plots/20260316_102709_S23_10_corridor_motif_injected_n12_trajectory.png`
- `plots/20260316_102709_S23_10_corridor_motif_injected_n24_trajectory.png`
- `plots/20260316_102709_S23_10_corridor_anisotropic_kernel_n12_trajectory.png`
- `plots/20260316_102709_S23_10_corridor_anisotropic_kernel_n24_trajectory.png`
- `plots/20260316_102709_S23_10_triad_baseline_n12_trajectory.png`
- `plots/20260316_102709_S23_10_triad_baseline_n24_trajectory.png`
- `plots/20260316_102709_S23_10_triad_degree_skewed_n12_trajectory.png`
- `plots/20260316_102709_S23_10_triad_degree_skewed_n24_trajectory.png`
- `plots/20260316_102709_S23_10_distributed_random_seed_n12_trajectory.png`
- `plots/20260316_102709_S23_10_distributed_random_seed_n24_trajectory.png`
- `plots/20260316_102709_stage23_10_topology_survival_matrix.png`
- `plots/20260316_102709_stage23_10_family_comparison_persistence_panel.png`
