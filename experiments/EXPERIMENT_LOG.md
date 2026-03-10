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
