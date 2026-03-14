# Atlas-1 Perturbation Resilience v1

Timestamped summary: `data/20260314_153736_stage10_atlas1_perturbation_resilience.json`

Timestamped run table: `data/20260314_153736_stage10_atlas1_perturbation_resilience.csv`

Pilot baseline seeds: A0_transverse_periodic_wide_midk_amp05, A0_transverse_open_wide_lowk_amp05, A0_scalar_periodic_wide_lowk_amp05, A0_scalar_periodic_medium_lowk_amp05, A0_scalar_open_medium_midk_amp05, A0_transverse_open_wide_highk_amp05, A0_scalar_periodic_medium_highk_amp05, A0_transverse_open_narrow_highk_amp05, A0_scalar_periodic_wide_highk_amp05, A0_scalar_open_narrow_highk_amp05

Perturbed regime counts: {'Ballistic coherent': 12, 'Ballistic dispersive': 14, 'Diffusive': 6, 'Chaotic or irregular': 4, 'Fragmenting': 4}

Transition counts: {'Ballistic coherent -> Ballistic coherent': 12, 'Ballistic dispersive -> Ballistic dispersive': 12, 'Diffusive -> Diffusive': 6, 'Diffusive -> Ballistic dispersive': 2, 'Chaotic or irregular -> Chaotic or irregular': 4, 'Fragmenting -> Fragmenting': 4}

Persistence by baseline regime: {'Ballistic coherent': 1.0, 'Ballistic dispersive': 1.0, 'Chaotic or irregular': 1.0, 'Diffusive': 0.75, 'Fragmenting': 1.0}

Persistence by perturbation axis: {'anisotropic_bias': 0.9, 'edge_weight_noise': 1.0, 'kernel_width_jitter': 1.0, 'sparse_graph_defects': 0.9}

Official Atlas-1 labels are evaluated with the frozen divergence constraint because the transverse projector remains frozen from Atlas-0. The perturbed divergence constraint is retained in the run table as an auxiliary sensitivity diagnostic.

Interpretation boundary: Atlas-1 maps regime persistence and transition under mild perturbation only. It does not assert continuum reconstruction, particle content, or geometry.
