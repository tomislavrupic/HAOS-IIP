# EXPERIMENT_LOG

This file records orchestrated experiment runs.

## Laplacian geometry test
- Date: 2026-03-09T22:31:38
- Config: substrate=random_geometric, nodes=500, epsilon=0.2, seed=42
- Results: `data/20260309_223137_laplacian_modes.json`
- Plots: `plots/20260309_223137_laplacian_modes_laplacian_modes_spectrum.png`, `plots/20260309_223137_laplacian_modes_laplacian_modes_mode1.png`
- Observation: low spectrum contains localized modes
- Conclusion: geometry sector mixed with substrate or boundary artifacts

## Gauge sector test
- Date: 2026-03-09T22:31:38
- Config: substrate=cubic_lattice, nodes=216, epsilon=0.2, seed=42
- Results: `data/20260309_223137_gauge_modes.json`
- Plots: `plots/20260309_223137_gauge_modes_gauge_modes_spectrum.png`, `plots/20260309_223137_gauge_modes_gauge_modes_phase_mode.png`
- Observation: phase dressing produces circulating low modes
- Conclusion: connection-sensitive gauge branch present, but still scalar in background transport

## Laplacian geometry test
- Date: 2026-03-09T22:33:52
- Config: substrate=random_geometric, nodes=500, epsilon=0.2, seed=42
- Results: `data/20260309_223352_laplacian_modes.json`
- Plots: `plots/20260309_223352_laplacian_modes_laplacian_modes_spectrum.png`, `plots/20260309_223352_laplacian_modes_laplacian_modes_mode1.png`
- Observation: low spectrum contains localized modes
- Conclusion: geometry sector mixed with substrate or boundary artifacts

## Gauge sector test
- Date: 2026-03-09T22:33:52
- Config: substrate=cubic_lattice, nodes=216, epsilon=0.2, seed=42
- Results: `data/20260309_223352_gauge_modes.json`
- Plots: `plots/20260309_223352_gauge_modes_gauge_modes_spectrum.png`, `plots/20260309_223352_gauge_modes_gauge_modes_phase_mode.png`
- Observation: phase dressing produces circulating low modes
- Conclusion: connection-sensitive gauge branch present, but still scalar in background transport

## Laplacian geometry test
- Date: 2026-03-09T22:44:48
- Config: substrate=random_geometric, nodes=500, epsilon=0.2, seed=42
- Results: `data/20260309_224448_laplacian_modes.json`
- Plots: `plots/20260309_224448_laplacian_modes_laplacian_modes_spectrum.png`, `plots/20260309_224448_laplacian_modes_laplacian_modes_mode1.png`
- Observation: low spectrum contains localized modes
- Conclusion: geometry sector mixed with substrate or boundary artifacts

## Gauge sector test
- Date: 2026-03-09T22:44:48
- Config: substrate=cubic_lattice, nodes=216, epsilon=0.2, seed=42
- Results: `data/20260309_224448_gauge_modes.json`
- Plots: `plots/20260309_224448_gauge_modes_gauge_modes_spectrum.png`, `plots/20260309_224448_gauge_modes_gauge_modes_phase_mode.png`
- Observation: phase dressing produces circulating low modes
- Conclusion: connection-sensitive gauge branch present, but still scalar in background transport

## Hodge L1 test
- Date: 2026-03-09T22:44:49
- Config: substrate=cubic_lattice, nodes=216, epsilon=0.2, seed=42
- Results: `data/20260309_224448_hodge_modes.json`
- Plots: `plots/20260309_224448_hodge_modes_hodge_modes_spectrum.png`, `plots/20260309_224448_hodge_modes_hodge_modes_sector_split.png`, `plots/20260309_224448_hodge_modes_hodge_modes_mode.png`
- Observation: low edge spectrum is dominated by exact gradient-like modes
- Conclusion: edge branch is present but not yet separated into a clear vector sector

## Parameter sweep
- Date: 2026-03-09T22:44:49
- Config: substrate=None, nodes=None, epsilon=None, seed=None
- Results: `data/20260309_224448_parameter_sweep.json`
- Plots: `plots/20260309_224448_parameter_sweep_parameter_sweep_laplacian_gap.png`, `plots/20260309_224448_parameter_sweep_parameter_sweep_laplacian_ipr.png`, `plots/20260309_224448_parameter_sweep_parameter_sweep_hodge_coexact.png`
- Observation: cubic substrates stay smoother than random graphs across the sweep
- Conclusion: substrate regularity matters strongly for clean low-mode geometry, while the vector branch remains only partially separated

## Laplacian geometry test
- Date: 2026-03-09T22:45:36
- Config: substrate=random_geometric, nodes=500, epsilon=0.2, seed=42
- Results: `data/20260309_224536_laplacian_modes.json`
- Plots: `plots/20260309_224536_laplacian_modes_laplacian_modes_spectrum.png`, `plots/20260309_224536_laplacian_modes_laplacian_modes_mode1.png`
- Observation: low spectrum contains localized modes
- Conclusion: geometry sector mixed with substrate or boundary artifacts

## Gauge sector test
- Date: 2026-03-09T22:45:36
- Config: substrate=cubic_lattice, nodes=216, epsilon=0.2, seed=42
- Results: `data/20260309_224536_gauge_modes.json`
- Plots: `plots/20260309_224536_gauge_modes_gauge_modes_spectrum.png`, `plots/20260309_224536_gauge_modes_gauge_modes_phase_mode.png`
- Observation: phase dressing produces circulating low modes
- Conclusion: connection-sensitive gauge branch present, but still scalar in background transport

## Hodge L1 test
- Date: 2026-03-09T22:45:36
- Config: substrate=cubic_lattice, nodes=216, epsilon=0.2, seed=42
- Results: `data/20260309_224536_hodge_modes.json`
- Plots: `plots/20260309_224536_hodge_modes_hodge_modes_spectrum.png`, `plots/20260309_224536_hodge_modes_hodge_modes_sector_split.png`, `plots/20260309_224536_hodge_modes_hodge_modes_mode.png`
- Observation: low edge spectrum is dominated by exact gradient-like modes
- Conclusion: edge branch is present but not yet separated into a clear vector sector

## Parameter sweep
- Date: 2026-03-09T22:45:37
- Config: substrates=['random_geometric', 'cubic_lattice'], nodes=[216, 343], epsilons=[0.12, 0.2, 0.28]
- Results: `data/20260309_224536_parameter_sweep.json`
- Plots: `plots/20260309_224536_parameter_sweep_parameter_sweep_laplacian_gap.png`, `plots/20260309_224536_parameter_sweep_parameter_sweep_laplacian_ipr.png`, `plots/20260309_224536_parameter_sweep_parameter_sweep_hodge_coexact.png`
- Observation: cubic substrates stay smoother than random graphs across the sweep
- Conclusion: substrate regularity matters strongly for clean low-mode geometry, while the vector branch remains only partially separated

## Periodic / twisted L1 experiment
- Date: 2026-03-09T23:02:39
- Config: epsilon=0.2, sizes=[4, 5], flux_quanta=[0, 1], open_compare_side=5
- Results: `data/20260309_230015_periodic_twisted_l1.json`
- Plots: `plots/20260309_230015_periodic_twisted_l1_periodic_twisted_l1_spectral_flow.png`, `plots/20260309_230015_periodic_twisted_l1_periodic_twisted_l1_open_vs_periodic.png`, `plots/20260309_230015_periodic_twisted_l1_periodic_twisted_l1_fractions.png`, `plots/20260309_230015_periodic_twisted_l1_periodic_twisted_l1_divergence.png`, `plots/20260309_230015_periodic_twisted_l1_periodic_twisted_l1_curl.png`, `plots/20260309_230015_periodic_twisted_l1_periodic_twisted_l1_modes.png`
- Observation: periodic topology produces a low divergence-free family absent in the open box, and nontrivial flux shifts it cleanly
- Conclusion: the periodic/twisted L1 branch shows a robust non-scalar flux-sensitive edge sector, but it is not yet a Maxwell-like propagating family

## Periodic / twisted L1 flux scan
- Date: 2026-03-09T23:21:41
- Config: epsilon=0.2, sizes=[4, 5], flux_quanta=[0, 1, 2, 3, 4], open_compare_sizes=[4, 5], low_modes=4
- Results: `data/20260309_232007_periodic_twisted_l1_flux_scan.json`
- Plots: `plots/20260309_232007_periodic_twisted_l1_flux_scan_periodic_twisted_l1_flux_scan_spectral_flow.png`, `plots/20260309_232007_periodic_twisted_l1_flux_scan_periodic_twisted_l1_flux_scan_fractions.png`, `plots/20260309_232007_periodic_twisted_l1_flux_scan_periodic_twisted_l1_flux_scan_divergence.png`, `plots/20260309_232007_periodic_twisted_l1_flux_scan_periodic_twisted_l1_flux_scan_curl.png`, `plots/20260309_232007_periodic_twisted_l1_flux_scan_periodic_twisted_l1_flux_scan_persistence.png`, `plots/20260309_232007_periodic_twisted_l1_flux_scan_periodic_twisted_l1_flux_scan_open_vs_periodic.png`, `plots/20260309_232007_periodic_twisted_l1_flux_scan_periodic_twisted_l1_flux_scan_modes.png`
- Observation: the low L1 branch stays mostly coexact across the scan, but its ordering under flux is not yet a single clean smooth band
- Conclusion: the scan shows some departure from pure torus-cycle behavior, but no isolated propagating Maxwell-like band is established yet

## Periodic / twisted L1 Hodge projection
- Date: 2026-03-09T23:58:31
- Config: epsilon=0.2, sizes=[4, 5], flux_quanta=[0, 1, 2, 3, 4], low_modes=6
- Results: `data/20260309_235831_periodic_twisted_l1_hodge_projection.json`
- Plots: `plots/20260309_235831_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_full_vs_coexact.png`, `plots/20260309_235831_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_harmonic_fraction.png`, `plots/20260309_235831_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_coexact_fraction.png`, `plots/20260309_235831_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_coexact_flow.png`, `plots/20260309_235831_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_persistence.png`, `plots/20260309_235831_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_modes.png`
- Observation: after harmonic removal, a low coexact branch persists but its ordering remains irregular or finite-size sensitive
- Conclusion: a low coexact band exists, but it still looks dominated by topological remnants or finite-size effects rather than a clean vector band

## Periodic / twisted L1 Hodge projection
- Date: 2026-03-10T00:00:47
- Config: epsilon=0.2, sizes=[4, 5], flux_quanta=[0, 1, 2, 3, 4], low_modes=6
- Results: `data/20260310_000047_periodic_twisted_l1_hodge_projection.json`
- Plots: `plots/20260310_000047_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_full_vs_coexact.png`, `plots/20260310_000047_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_harmonic_fraction.png`, `plots/20260310_000047_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_coexact_fraction.png`, `plots/20260310_000047_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_coexact_flow.png`, `plots/20260310_000047_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_persistence.png`, `plots/20260310_000047_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_modes.png`
- Observation: after harmonic removal, a low coexact branch persists but its ordering remains irregular or finite-size sensitive
- Conclusion: a low coexact band exists, but it still looks dominated by topological remnants or finite-size effects rather than a clean vector band

## Periodic / twisted L1 Hodge projection
- Date: 2026-03-10T00:02:45
- Config: epsilon=0.2, sizes=[4, 5], flux_quanta=[0, 1, 2, 3, 4], low_modes=6
- Results: `data/20260310_000245_periodic_twisted_l1_hodge_projection.json`
- Plots: `plots/20260310_000245_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_full_vs_coexact.png`, `plots/20260310_000245_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_harmonic_fraction.png`, `plots/20260310_000245_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_coexact_fraction.png`, `plots/20260310_000245_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_coexact_flow.png`, `plots/20260310_000245_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_persistence.png`, `plots/20260310_000245_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_modes.png`
- Observation: after harmonic removal, the coexact floor stays well above the mixed full low branch across the flux scan
- Conclusion: no clear coexact low band yet; the very low modes remain dominated by harmonic or mixed topological remnants

## Periodic / twisted L1 Hodge projection
- Date: 2026-03-10T00:03:08
- Config: epsilon=0.2, sizes=[4, 5], flux_quanta=[0, 1, 2, 3, 4], low_modes=6
- Results: `data/20260310_000308_periodic_twisted_l1_hodge_projection.json`
- Plots: `plots/20260310_000308_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_full_vs_coexact.png`, `plots/20260310_000308_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_harmonic_fraction.png`, `plots/20260310_000308_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_coexact_fraction.png`, `plots/20260310_000308_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_coexact_flow.png`, `plots/20260310_000308_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_persistence.png`, `plots/20260310_000308_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_modes.png`
- Observation: after harmonic removal, the coexact floor stays well above the mixed full low branch across the flux scan
- Conclusion: no clear coexact low band yet; the very low modes remain dominated by harmonic or mixed topological remnants

## Recoverable hydrogenic spectrum
- Date: 2026-03-10T10:46:08
- Config: D=1.0, alpha=1.0, r_min=0.001, r_max=80.0, n_grid=1400, l_values=[0, 1, 2], states_per_l=4
- Results: `data/20260310_104607_recoverable_hydrogenic_spectrum.json`
- Plots: `plots/20260310_104607_recoverable_hydrogenic_energy_scaling.png`, `plots/20260310_104607_recoverable_hydrogenic_n2_energy.png`, `plots/20260310_104607_recoverable_hydrogenic_modes.png`
- Observation: bound radial modes form a discrete ladder with near-constant n^2 |E_n|
- Conclusion: discreteness emerges from quadratic dispersion, inverse-distance geometry, and normalizability without a separate quantization postulate

## Interaction attractor spectrum
- Date: 2026-03-10T10:50:45
- Config: D=1.0, alpha=1.0, r_min=0.001, r_max=80.0, n_grid=1200, states=4, modal_basis_size=12, dt=6.0, max_steps=400, tol=1e-10, recovery_noise=0.12, seed=42
- Results: `data/20260310_105044_interaction_attractor_spectrum.json`
- Plots: `plots/20260310_105044_interaction_attractor_energy_descent.png`, `plots/20260310_105044_interaction_attractor_scaling.png`, `plots/20260310_105044_interaction_attractor_recovery.png`, `plots/20260310_105044_interaction_attractor_modes.png`
- Observation: normalized interaction flow converges to a discrete ladder of recoverable bound attractors
- Conclusion: discrete states arise from the inverse-distance interaction geometry and boundary-constrained radial operator; the attractor flow stabilizes and recovers these modes, but unconstrained flow selects only the lowest sector

## Emergent inverse-square geometry
- Date: 2026-03-10T10:58:04
- Config: D=1.0, r_max=120.0, n_grid=1600, source_width=1.2, source_strength=1.0, dt=2.0, time_steps=220, fit_window=[4.0, 20.0]
- Results: `data/20260310_105804_emergent_inverse_square_geometry.json`
- Plots: `plots/20260310_105804_inverse_square_geometry_field.png`, `plots/20260310_105804_inverse_square_geometry_force.png`, `plots/20260310_105804_inverse_square_geometry_relaxation.png`, `plots/20260310_105804_inverse_square_geometry_source.png`
- Observation: the steady radial field generated by a localized disturbance follows an inverse-distance law outside the source core
- Conclusion: inverse-square geometry emerges from flux-conserving 3D substrate propagation once a localized disturbance relaxes into a recoverable stationary field

## Kernel graph Green response
- Date: 2026-03-10T11:05:40
- Config: n_sides=[9, 11, 13], epsilon_coeffs=[0.5, 1.0], cutoff_factor=2.5, fit_r_max=0.55
- Results: `data/20260310_110539_kernel_graph_green_response.json`
- Plots: `plots/20260310_110539_kernel_graph_green_profiles.png`, `plots/20260310_110539_kernel_graph_green_exponents.png`, `plots/20260310_110539_kernel_graph_green_fit_quality.png`
- Observation: kernel-graph Green response approaches A + B/r on the cubic substrate scan
- Conclusion: the weighted interaction kernel induces a graph Laplacian whose far field is consistent with inverse-distance geometry; direct shell-derivative force estimates remain noisier than the field fit

## Kernel Laplacian convergence
- Date: 2026-03-10T11:15:23
- Config: n_sides=[9, 11, 13, 17, 21], epsilon_coeffs=[0.5, 1.0, 2.0], cutoff_factor=2.5
- Results: `data/20260310_111522_kernel_operator_convergence.json`
- Plots: `plots/20260310_111522_operator_error_vs_h.png`, `plots/20260310_111522_operator_profiles.png`
- Observation: the kernel-induced operator reproduces the quadratic test exactly and converges toward the continuum Laplacian for smooth test functions
- Conclusion: after discrete second-moment normalization, the graph operator shows clear Laplacian convergence on the cubic scan; broader kernels remain less accurate at finite resolution

## L1 defect transverse test
- Date: 2026-03-10T13:03:26
- Config: epsilon=0.2, sizes=[6, 8, 10], variants=['baseline', 'puncture', 'line_defect'], phase_modes=48, restricted_modes=12
- Results: `data/20260310_130326_L1_defect_scan.json`
- Plots: `plots/20260310_130326_divergence_curl_phase.png`, `plots/20260310_130326_restricted_transverse_spectrum.png`, `plots/20260310_130326_restricted_eigenvalue_vs_n.png`
- Observation: the lowest restricted transverse eigenvalue decreases with lattice size in all three substrate branches, and both defect branches lie below the baseline torus at every tested n
- Conclusion: the restricted transverse sector develops a descending low band under the tested puncture and line-defect substrate modifications

## L1 transverse scaling test
- Date: 2026-03-10T13:23:23
- Config: epsilon=0.2, sizes=[6, 8, 10, 12, 14, 16], variants=['baseline', 'puncture', 'line_defect'], phase_modes=18, restricted_modes=6
- Results: `data/20260310_132323_L1_transverse_scaling.json`
- Plots: `plots/20260310_132323_transverse_floor_scaling_loglog.png`, `plots/20260310_132323_transverse_floor_ipr_vs_n.png`, `plots/20260310_132323_transverse_mode_profiles.png`, `plots/20260310_132323_divergence_curl_phase_scaling.png`
- Observation: the restricted floor keeps descending while the lowest-mode IPR drops with lattice size in every branch, and the defect branches do not retain a growing fraction of norm near the defect
- Conclusion: the descending restricted floor is behaving more like a continuum transverse-band candidate than a defect-pinned low mode in the tested range

## L1 transverse band test
- Date: 2026-03-10T13:33:22
- Config: epsilon=0.2, sizes=[6, 8, 10, 12, 14, 16], variants=['baseline', 'puncture', 'line_defect', 'flux_tube'], restricted_modes=20, flux_tube_phase=1.5707963267948966
- Results: `data/20260310_133322_L1_transverse_band_scan.json`
- Plots: `plots/20260310_133322_transverse_band_collapse.png`, `plots/20260310_133322_transverse_floor_vs_n.png`, `plots/20260310_133322_transverse_ipr_vs_n.png`, `plots/20260310_133322_divergence_curl_phase_band.png`, `plots/20260310_133322_transverse_mode_profiles.png`
- Observation: the first twenty restricted transverse eigenvalues show a partial n^2 spectral collapse while the lowest-mode IPR and defect concentration both fall with lattice size
- Conclusion: the restricted transverse sector is organizing as a stable low spectral band consistent with continuum scaling in the tested range

## Milestone: transverse-band formation test (2026-03-10)
- Kernel-induced substrate now supports:
  - Laplacian scalar sector
  - restricted transverse edge sector
- Numerical evidence:
  - descending restricted floor with lattice size
  - decreasing IPR and defect concentration
  - partial n^2 spectral collapse of first modes
- Interpretation boundary:
  - operator-level evidence for a low transverse spectral band consistent with continuum scaling in the tested range
  - no claim of gauge dynamics or emergent photons

## kernel width sensitivity test
- Date: 2026-03-10T14:45:13
- Config: n=16, c_epsilon_values=[0.25, 0.5, 1.0, 1.5], phase_modes=20, restricted_modes=20
- Results: `data/20260310_144513_kernel_width_sensitivity_test.json`
- Plots: `plots/20260310_144513_transverse_band_vs_kernel_width.png`, `plots/20260310_144513_divergence_curl_phase_kernel_width.png`
- Observation: the restricted transverse spectrum remains present across the tested kernel-width range, while the lowest eigenvalue, first-band spread, and lowest-mode IPR shift smoothly with c_epsilon
- Conclusion: within the tested local-kernel range, the restricted transverse band is robust to moderate kernel-width changes

## L1 geometric disorder test
- Date: 2026-03-10T14:46:21
- Config: epsilon=0.2, sizes=[12, 16, 20], sigma_factor=0.05, phase_modes=20, restricted_modes=20
- Results: `data/20260310_144621_L1_geometric_disorder_test.json`
- Plots: `plots/20260310_144621_transverse_band_disorder.png`, `plots/20260310_144621_divergence_curl_phase_disorder.png`
- Observation: the restricted transverse spectrum persists under weak geometric disorder, with the ordered and disordered scaled bands remaining close and the lowest restricted mode staying divergence-free
- Conclusion: in the tested weak-disorder regime, loss of exact lattice symmetry does not destroy the low restricted transverse band

## L1 large-n scaling test
- Date: 2026-03-10T14:50:17
- Config: epsilon=0.2, sizes=[12, 16, 20, 24, 28, 32], phase_modes=30, restricted_modes=30, variant=baseline
- Results: `data/20260310_145017_L1_large_n_scaling_test.json`
- Plots: `plots/20260310_145017_transverse_floor_vs_n_large.png`, `plots/20260310_145017_transverse_band_collapse_large.png`, `plots/20260310_145017_transverse_ipr_vs_n_large.png`, `plots/20260310_145017_divergence_curl_phase_large.png`
- Observation: the lowest restricted transverse eigenvalue keeps decreasing on the baseline torus, the lowest-mode IPR continues to fall, and the first thirty scaled eigenvalues show partial n^2 collapse
- Conclusion: the large-n baseline scan remains consistent with a continuum Laplacian-type transverse band in the tested range

## transverse mode structure
- Date: 2026-03-10T14:50:35
- Config: epsilon=0.2, n=16, variant=baseline, phase_modes=20, restricted_modes=10
- Results: `data/20260310_145035_transverse_mode_structure.json`
- Plots: `plots/20260310_145035_transverse_vector_field_profiles.png`, `plots/20260310_145035_divergence_curl_phase_mode_structure.png`
- Observation: the lowest restricted transverse modes remain numerically divergence-free while their spatial support extends over the periodic lattice rather than concentrating at a single location
- Conclusion: the inspected low restricted modes behave like extended circulation fields in the tested baseline branch

## transverse continuum comparison
- Date: 2026-03-10T15:00:06
- Config: epsilon=0.2, sizes=[12, 16, 20], variants=['baseline', 'puncture', 'line_defect', 'flux_tube'], restricted_modes=20
- Results: `data/20260310_150006_transverse_continuum_comparison.json`
- Plots: `plots/20260310_150006_transverse_continuum_comparison.png`, `plots/20260310_150006_transverse_mode_spacing.png`, `plots/20260310_150006_divergence_curl_phase_continuum.png`
- Observation: after n^2 rescaling, the low restricted transverse spectrum aligns with the continuum transverse mode ordering and the mode-spacing pattern remains stable across the tested branches
- Conclusion: for the largest baseline case, the first-ten relative spectrum error is 0.219 and the first-ten spacing error is 0.000, consistent with a low continuum transverse operator in the tested range

## transverse PDE reconstruction
- Date: 2026-03-10T15:48:31
- Config: epsilon=0.2, sizes=[12, 16], restricted_modes=3
- Results: `data/20260310_154831_transverse_pde_reconstruction.json`
- Plots: `plots/20260310_154831_pde_reconstruction_residuals.png`, `plots/20260310_154831_transverse_field_residual_maps.png`, `plots/20260310_154831_divergence_curl_phase_pde_reconstruction.png`
- Observation: after coarse-graining to a periodic vector field, the lowest restricted modes remain nearly divergence-free and the curl-curl residual decreases with lattice size
- Conclusion: for the lowest mode, the relative curl-curl residual decreases from 0.051 at n=12 to 0.029 at n=16, consistent with an emergent local coarse-grained PDE in the tested range

## transverse real-space correlations
- Date: 2026-03-10T15:51:53
- Config: epsilon=0.2, sizes=[12, 16], pair_samples=8000, bins=24, restricted_modes=6
- Results: `data/20260310_155153_transverse_realspace_correlations.json`
- Plots: `plots/20260310_155153_transverse_correlation_decay.png`, `plots/20260310_155153_transverse_anisotropy.png`, `plots/20260310_155153_divergence_curl_phase_correlations.png`
- Observation: the lowest restricted transverse mode shows extended two-point correlations while directional anisotropy remains modest across the tested sizes
- Conclusion: the correlation length remains of the same order, from 0.094 at n=12 to 0.088 at n=16, while the anisotropy ratio stays moderate at 1.197 for the largest size, consistent with an extended continuum-like branch in the tested range

## transverse effective operator fit
- Date: 2026-03-10T15:57:13
- Config: epsilon=0.2, sizes=[12, 16, 20, 24], restricted_modes=24, fit_modes=2
- Results: `data/20260310_155713_transverse_effective_operator_fit.json`
- Plots: `plots/20260310_155713_effective_operator_fit.png`, `plots/20260310_155713_effective_gap_vs_n.png`, `plots/20260310_155713_effective_stiffness_vs_n.png`, `plots/20260310_155713_divergence_curl_phase_effective_fit.png`
- Observation: the low restricted spectrum on the clean periodic branch is well fit by a linear function of effective momentum squared across the tested lattice sizes
- Conclusion: the fitted gap stays small, from -8.275e-17 at n=12 to 0.000e+00 at n=24, while the effective stiffness varies only mildly across the tested sizes, consistent with a nearly gapless continuum transverse branch in the tested range

## scalar-transverse coupling test
- Date: 2026-03-10T16:40:41
- Config: epsilon=0.2, sizes=[12], variants=['baseline', 'puncture', 'line_defect'], families=['radial', 'anisotropic', 'bump'], etas=[0.02, 0.05, 0.1], restricted_modes=20
- Results: `data/20260310_164041_scalar_transverse_coupling_test.json`
- Plots: `plots/20260310_164041_scalar_transverse_eigenvalue_shifts.png`, `plots/20260310_164041_scalar_transverse_mode_overlap.png`, `plots/20260310_164041_scalar_transverse_deformation_response.png`, `plots/20260310_164041_divergence_curl_phase_scalar_transverse.png`
- Observation: small scalar/background deformations shift the restricted transverse spectrum smoothly, while the matched low-mode overlaps remain moderate to high depending on branch and deformation family
- Conclusion: for the baseline anisotropic branch at eta=0.10, the lowest restricted eigenvalue shifts by -0.000479 with mean matched overlap 0.4121, indicating a measurable but bounded operator-level spectral response

## scalar-transverse coupling matrix
- Date: 2026-03-10T16:42:16
- Config: epsilon=0.2, n_side=12, variant=baseline, family=anisotropic, etas=[0.02, 0.05, 0.1], matrix_modes=10
- Results: `data/20260310_164216_scalar_transverse_coupling_matrix.json`
- Plots: `plots/20260310_164216_transverse_coupling_matrix_heatmap.png`, `plots/20260310_164216_transverse_mode_mixing_vs_eta.png`
- Observation: the deformation-induced coupling matrix remains structured, but its weight is broadly distributed across the low-mode block rather than sharply diagonal-dominant
- Conclusion: at eta=0.10, the diagonal fraction is 0.3365 and the off-diagonal fraction is 0.6635, indicating structured but non-diagonal mode mixing in the low transverse sector

## transverse backreaction maps
- Date: 2026-03-10T16:44:15
- Config: epsilon=0.2, n_side=12, variants=['baseline', 'puncture'], family=bump, etas=[0.05, 0.1]
- Results: `data/20260310_164414_transverse_backreaction_maps.json`
- Plots: `plots/20260310_164414_transverse_backreaction_maps.png`, `plots/20260310_164414_transverse_sensitivity_vs_scalar_profile.png`
- Observation: localized scalar/background bumps produce nonuniform changes in the lowest restricted transverse mode rather than a flat redistribution across the lattice
- Conclusion: for the baseline branch at eta=0.10, the mean sensitivity is 0.008907 with profile correlation -0.0934, so the response is spatially structured but not simply proportional to the scalar bump profile

## transverse gap response
- Date: 2026-03-10T16:46:00
- Config: epsilon=0.2, sizes=[12], variant=baseline, families=['radial', 'anisotropic', 'bump'], etas=[0.0, 0.02, 0.05, 0.1]
- Results: `data/20260310_164600_transverse_gap_response.json`
- Plots: `plots/20260310_164600_transverse_gap_vs_deformation.png`, `plots/20260310_164600_transverse_stiffness_vs_deformation.png`
- Observation: small scalar/background deformations shift the effective transverse fit mainly through smooth stiffness renormalization, with no large unstable gap opening in the tested range
- Conclusion: for the anisotropic family at n=12 and eta=0.10, the fitted stiffness is 23.654496 and the fitted gap is 9.884579e-02, consistent with a structured operator-level response rather than chaotic spectral rearrangement

## transverse gap response
- Date: 2026-03-10T16:57:17
- Config: epsilon=0.2, sizes=[12, 16], variant=baseline, families=['radial', 'anisotropic', 'bump'], etas=[0.0, 0.02, 0.05, 0.1]
- Results: `data/20260310_165717_transverse_gap_response.json`
- Plots: `plots/20260310_165717_transverse_gap_vs_deformation.png`, `plots/20260310_165717_transverse_stiffness_vs_deformation.png`
- Observation: small scalar/background deformations shift the effective transverse fit mainly through smooth stiffness renormalization, with no large unstable gap opening in the tested range
- Conclusion: for the anisotropic family at n=16 and eta=0.10, the fitted stiffness is 38.561330 and the fitted gap is 5.394913e-05, consistent with a structured operator-level response rather than chaotic spectral rearrangement


## DH Hermiticity Test
- Date: 2026-03-10T17:32:41
- Config: epsilon=0.2, two_d_n=6, three_d_n=8, prototype=full chiral 4-spinor lift
- Results: `data/20260310_172558_DH_hermiticity_test.json`
- Plots: 
- Observation: the chiral Dirac-type lift is Hermitian to numerical precision on both the 2D and 3D prototype grids
- Conclusion: the Stage 5 prototype passes the Hermiticity gate and is admissible for the square-root and spectral tests

## DH Square-Root Test
- Date: 2026-03-10T17:32:41
- Config: epsilon=0.2, two_d_n=6, three_d_n=8, modes=48, prototype=full chiral 4-spinor lift
- Results: `data/20260310_172136_DH_square_root_test.json`
- Plots: `plots/20260310_172136_DH2_vs_L0_spectrum.png`
- Observation: the chiral Dirac prototype tracks the low scalar spectrum only qualitatively: eigenvalue correlations remain high, but the relative spectral and operator errors stay large on both prototype grids
- Conclusion: the current Stage 5 prototype does not satisfy the square-root criterion tightly enough; under the prompt rule, this lift is rejected in its present form

## DH Flux Response Test
- Date: 2026-03-10T17:32:41
- Config: epsilon=0.2, two_d_n=6, flux_quanta=[0, 1, 2, 3], spectral_flow_modes=12, prototype=full chiral 4-spinor lift
- Results: `data/20260310_172754_DH_flux_response_test.json`
- Plots: `plots/20260310_172754_DH_flux_spectral_flow.png`
- Observation: the low positive DH branch shifts smoothly as torus flux is increased in the 2D prototype
- Conclusion: the prototype has smooth 2D background-flux spectral flow, but it remains rejected overall because the square-root criterion against the scalar operator failed

## DH Spectrum Structure
- Date: 2026-03-10T17:32:41
- Config: epsilon=0.2, two_d_n=6, three_d_n=8, modes=50, prototype=full chiral 4-spinor lift
- Results: `data/20260310_173100_DH_spectrum_structure.json`
- Plots: `plots/20260310_173100_DH_spectrum.png`, `plots/20260310_173100_DH_ipr_vs_mode.png`
- Observation: the low DH spectrum is sign-symmetric and extended in the 2D prototype, while the 3D prototype remains extended but shows only partial near-zero pairing
- Conclusion: the prototype shows qualitative Dirac-like spectral structure, but this does not overcome the failed square-root test against the scalar operator

## DK 2D Cochain Infrastructure
- Date: 2026-03-10T17:49:55
- Config: epsilon=0.2, cochain_sizes=[8, 12, 16], cycle_phase_x=0.0, cycle_phase_y=0.0
- Results: `data/20260310_174955_DK_2D_cochain_infrastructure.json`
- Plots: 
- Observation: the weighted periodic 2D cochain complex satisfies the coboundary and adjoint-coboundary identities to numerical precision across the tested lattice sizes
- Conclusion: the Stage 6 2D cochain infrastructure is internally consistent and admissible for the Dirac-Kaehler square test

## DK Square Test
- Date: 2026-03-10T17:50:01
- Config: epsilon=0.2, cochain_sizes=[8, 12, 16], square_modes=40
- Results: `data/20260310_175001_DK_square_test.json`
- Plots: `plots/20260310_175001_DK_square_vs_Hodge.png`
- Observation: the graded Dirac-Kaehler square matches the graded Hodge Laplacian to numerical precision across the tested 2D lattice sizes
- Conclusion: the Stage 6 Dirac-Kaehler lift satisfies the defining square identity on the validated 2D cochain architecture

## DK Spectrum Structure
- Date: 2026-03-10T17:51:59
- Config: epsilon=0.2, n_side=12, modes=60
- Results: `data/20260310_175158_DK_spectrum_structure.json`
- Plots: `plots/20260310_175158_DK_spectrum.png`, `plots/20260310_175158_DK_mode_grading.png`, `plots/20260310_175158_DK_ipr_vs_mode.png`
- Observation: the low graded Dirac-Kaehler spectrum shows clean sign symmetry, low near-zero modes, and stable weight distribution across 0-, 1-, and 2-form sectors
- Conclusion: the Stage 6 2D prototype exhibits a stable graded spectral structure consistent with a Dirac-Kaehler lift on the validated cochain complex

## DK Flux Response Test
- Date: 2026-03-10T17:53:52
- Config: epsilon=0.2, n_side=12, cycle_phases=[0.0, 0.2617993877991494, 0.5235987755982988, 0.7853981633974483], spectral_flow_modes=12
- Results: `data/20260310_175351_DK_flux_response_test.json`
- Plots: `plots/20260310_175351_DK_flux_spectral_flow.png`, `plots/20260310_175351_DK_flux_square_residual.png`
- Observation: the low Dirac-Kaehler branch shifts smoothly under flat torus-cycle holonomy while the graded square residual remains at numerical precision
- Conclusion: the Stage 6 2D Dirac-Kaehler lift preserves its square identity under mild background cycle holonomy and shows stable spectral flow

## DK Flux Response Test
- Date: 2026-03-10T17:55:52
- Config: epsilon=0.2, n_side=12, cycle_phases=[0.0, 0.2617993877991494, 0.5235987755982988, 0.7853981633974483], spectral_flow_modes=12
- Results: `data/20260310_175551_DK_flux_response_test.json`
- Plots: `plots/20260310_175551_DK_flux_spectral_flow.png`, `plots/20260310_175551_DK_flux_square_residual.png`
- Observation: the low Dirac-Kaehler branch shifts smoothly under flat torus-cycle holonomy while the graded square residual remains at numerical precision
- Conclusion: the Stage 6 2D Dirac-Kaehler lift preserves its square identity under mild background cycle holonomy and shows stable spectral flow

## DK Flux Response Test
- Date: 2026-03-10T17:55:53
- Config: epsilon=0.2, n_side=12, cycle_phases=[0.0, 0.2617993877991494, 0.5235987755982988, 0.7853981633974483], spectral_flow_modes=12
- Results: `data/20260310_175553_DK_flux_response_test.json`
- Plots: `plots/20260310_175553_DK_flux_spectral_flow.png`, `plots/20260310_175553_DK_flux_square_residual.png`
- Observation: the low Dirac-Kaehler branch shifts smoothly under flat torus-cycle holonomy while the graded square residual remains at numerical precision
- Conclusion: the Stage 6 2D Dirac-Kaehler lift preserves its square identity under mild background cycle holonomy and shows stable spectral flow

## DK 3D Cochain Infrastructure
- Date: 2026-03-10T18:04:12
- Config: epsilon=0.2, cochain_sizes=[4, 6, 8], cycle_phase_x=0.0, cycle_phase_y=0.0, cycle_phase_z=0.0
- Results: `data/20260310_180412_DK_3D_cochain_infrastructure.json`
- Plots: 
- Observation: the weighted periodic 3D cochain complex satisfies the nilpotency and adjoint-nilpotency identities to numerical precision across the tested lattice sizes
- Conclusion: the Stage 7 3D cochain infrastructure is internally consistent and admissible for the graded square test

## DK 3D Square Test
- Date: 2026-03-10T18:04:19
- Config: epsilon=0.2, cochain_sizes=[4, 6, 8], square_modes=40
- Results: `data/20260310_180419_DK_3D_square_test.json`
- Plots: `plots/20260310_180419_DK_3D_square_vs_Hodge.png`
- Observation: the graded 3D Dirac-Kaehler square matches the graded Hodge Laplacian to numerical precision across the tested lattice sizes and on every form degree
- Conclusion: the Stage 7 3D Dirac-Kaehler lift satisfies the grade-by-grade square identity on the validated 3D cochain architecture

## DK 3D Spectrum Structure
- Date: 2026-03-10T18:04:33
- Config: epsilon=0.2, n_side=6, modes=60
- Results: `data/20260310_180433_DK_3D_spectrum_structure.json`
- Plots: `plots/20260310_180433_DK_3D_spectrum.png`, `plots/20260310_180433_DK_3D_mode_grading.png`, `plots/20260310_180433_DK_3D_ipr_vs_mode.png`
- Observation: the low graded 3D Dirac-Kaehler spectrum shows clean sign symmetry, near-zero modes, and stable weight distribution across 0-, 1-, 2-, and 3-form sectors
- Conclusion: the Stage 7 3D prototype exhibits a stable graded spectral structure consistent with the validated Dirac-Kaehler lift

## DK 3D Flux Response Test
- Date: 2026-03-10T18:04:36
- Config: epsilon=0.2, n_side=4, cycle_phases=[0.0, 0.2617993877991494, 0.5235987755982988, 0.7853981633974483], spectral_flow_modes=12
- Results: `data/20260310_180436_DK_3D_flux_response_test.json`
- Plots: `plots/20260310_180436_DK_3D_flux_spectral_flow.png`, `plots/20260310_180436_DK_3D_flux_square_residual.png`
- Observation: the low 3D Dirac-Kaehler branch shifts smoothly under flat x-cycle holonomy while the graded square residual remains at numerical precision
- Conclusion: the Stage 7 3D Dirac-Kaehler lift preserves its square identity under mild torus-cycle holonomy and shows stable low spectral flow
