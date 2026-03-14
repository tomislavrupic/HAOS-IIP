# Atlas-2 Transition Structure v1

Timestamped summary: `data/20260314_161318_stage10_atlas2_transition_structure.json`

Timestamped run table: `data/20260314_161318_stage10_atlas2_transition_structure.csv`

Pilot baseline seeds: A0_transverse_periodic_wide_midk_amp05, A0_transverse_open_wide_lowk_amp05, A0_scalar_periodic_wide_lowk_amp05, A0_scalar_periodic_medium_lowk_amp05, A0_scalar_open_medium_midk_amp05, A0_transverse_open_wide_highk_amp05, A0_scalar_periodic_medium_highk_amp05, A0_transverse_open_narrow_highk_amp05, A0_scalar_periodic_wide_highk_amp05, A0_scalar_open_narrow_highk_amp05

Perturbed regime counts: {'Ballistic coherent': 6, 'Ballistic dispersive': 8, 'Diffusive': 2, 'Chaotic or irregular': 2, 'Fragmenting': 2}

Transition counts: {'Ballistic coherent -> Ballistic coherent': 6, 'Ballistic dispersive -> Ballistic dispersive': 6, 'Diffusive -> Ballistic dispersive': 2, 'Diffusive -> Diffusive': 2, 'Chaotic or irregular -> Chaotic or irregular': 2, 'Fragmenting -> Fragmenting': 2}

Transition-type counts: {'stable': 18, 'regime_shift': 1, 'regime_softening': 1}

Persistence by perturbation pair: {'defects_plus_kernel_jitter': 0.9, 'edge_noise_plus_anisotropic_bias': 0.9}

Official Atlas-2 labels are evaluated with the frozen divergence constraint because the transverse projector remains frozen from Atlas-0. The paired perturbed divergence constraint is retained in the run table as an auxiliary sensitivity diagnostic.

Interpretation boundary: Atlas-2 maps transition structure in regime space under paired stress only. It does not assert continuum reconstruction, particle content, or geometry.
