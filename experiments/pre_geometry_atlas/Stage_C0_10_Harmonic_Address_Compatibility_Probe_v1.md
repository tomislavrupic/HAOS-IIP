# Stage C0.10 Harmonic Address Compatibility Probe v1

Timestamped JSON: `data/20260315_213225_stage_c0_10_harmonic_address_compatibility_probe.json`
Timestamped CSV: `data/20260315_213225_stage_c0_10_harmonic_address_compatibility_probe.csv`

Architecture notice: This branch keeps the clustered Phase III DK seed geometry, the projected-transverse sector, the fixed graph resolution, and the Stage 23 topology measurement stack in place. Only a discrete harmonic-address compatibility factor is added as a static relational multiplier on the combinatorial operator.

Stage summary: classification counts {'address_protected_braid': 2, 'address_sorted_encounter': 6, 'address_induced_smear': 1}; topology counts {'braid_like_exchange': 2, 'transfer_smeared': 7}.

Per-run summary:

- `C010_A_exact_symmetry_broken`
  - group: `Fully compatible`
  - topology class: `transfer_smeared`
  - classification: `address_sorted_encounter`
  - persistence gain: `0.0000`
  - flow concentration: `0.9050`
  - transfer smear index: `0.3267`
  - encounter sorting score: `0.5211`
  - address selectivity index: `0.8116`
- `C010_A_exact_uniform_00`
  - group: `Fully compatible`
  - topology class: `braid_like_exchange`
  - classification: `address_protected_braid`
  - persistence gain: `0.0000`
  - flow concentration: `0.8989`
  - transfer smear index: `0.3195`
  - encounter sorting score: `1.0000`
  - address selectivity index: `1.0000`
- `C010_A_exact_uniform_11`
  - group: `Fully compatible`
  - topology class: `braid_like_exchange`
  - classification: `address_protected_braid`
  - persistence gain: `0.0000`
  - flow concentration: `0.8989`
  - transfer smear index: `0.3195`
  - encounter sorting score: `1.0000`
  - address selectivity index: `1.0000`
- `C010_B_near_forward_01`
  - group: `Near compatible`
  - topology class: `transfer_smeared`
  - classification: `address_sorted_encounter`
  - persistence gain: `0.0000`
  - flow concentration: `0.8931`
  - transfer smear index: `0.3306`
  - encounter sorting score: `0.8071`
  - address selectivity index: `0.9299`
- `C010_B_near_reverse_10`
  - group: `Near compatible`
  - topology class: `transfer_smeared`
  - classification: `address_sorted_encounter`
  - persistence gain: `0.0000`
  - flow concentration: `0.8931`
  - transfer smear index: `0.3306`
  - encounter sorting score: `0.8071`
  - address selectivity index: `0.9299`
- `C010_B_near_shifted_12`
  - group: `Near compatible`
  - topology class: `transfer_smeared`
  - classification: `address_sorted_encounter`
  - persistence gain: `0.0000`
  - flow concentration: `0.8931`
  - transfer smear index: `0.3306`
  - encounter sorting score: `0.8071`
  - address selectivity index: `0.9299`
- `C010_C_incompatible_02`
  - group: `Incompatible`
  - topology class: `transfer_smeared`
  - classification: `address_sorted_encounter`
  - persistence gain: `0.0000`
  - flow concentration: `0.8852`
  - transfer smear index: `0.3411`
  - encounter sorting score: `0.6763`
  - address selectivity index: `0.8382`
- `C010_C_incompatible_13`
  - group: `Incompatible`
  - topology class: `transfer_smeared`
  - classification: `address_sorted_encounter`
  - persistence gain: `0.0000`
  - flow concentration: `0.8852`
  - transfer smear index: `0.3411`
  - encounter sorting score: `0.6763`
  - address selectivity index: `0.8382`
- `C010_C_randomized_control`
  - group: `Incompatible`
  - topology class: `transfer_smeared`
  - classification: `address_induced_smear`
  - persistence gain: `0.0000`
  - flow concentration: `0.8988`
  - transfer smear index: `0.3376`
  - encounter sorting score: `0.0881`
  - address selectivity index: `0.5990`

Interpretive note:
- Harmonic addressing behaves like a weak proto-selection rule on this clustered seed: exact-match labeling preserves the braid class, while near-match and incompatible sectors relax to smeared transfer. The visible signal is topology and sorting bias, not measurable persistence gain.

## Plots
- `plots/20260315_213225_stage_c0_10_C010_A_exact_uniform_00_topology_trajectory.png`
- `plots/20260315_213225_stage_c0_10_C010_A_exact_uniform_00_flow_concentration_trace.png`
- `plots/20260315_213225_stage_c0_10_C010_A_exact_uniform_00_grade_exchange_trace.png`
- `plots/20260315_213225_stage_c0_10_C010_A_exact_uniform_00_grade_weights.png`
- `plots/20260315_213225_stage_c0_10_C010_A_exact_uniform_11_topology_trajectory.png`
- `plots/20260315_213225_stage_c0_10_C010_A_exact_uniform_11_flow_concentration_trace.png`
- `plots/20260315_213225_stage_c0_10_C010_A_exact_uniform_11_grade_exchange_trace.png`
- `plots/20260315_213225_stage_c0_10_C010_A_exact_uniform_11_grade_weights.png`
- `plots/20260315_213225_stage_c0_10_C010_A_exact_symmetry_broken_topology_trajectory.png`
- `plots/20260315_213225_stage_c0_10_C010_A_exact_symmetry_broken_flow_concentration_trace.png`
- `plots/20260315_213225_stage_c0_10_C010_A_exact_symmetry_broken_grade_exchange_trace.png`
- `plots/20260315_213225_stage_c0_10_C010_A_exact_symmetry_broken_grade_weights.png`
- `plots/20260315_213225_stage_c0_10_C010_B_near_forward_01_topology_trajectory.png`
- `plots/20260315_213225_stage_c0_10_C010_B_near_forward_01_flow_concentration_trace.png`
- `plots/20260315_213225_stage_c0_10_C010_B_near_forward_01_grade_exchange_trace.png`
- `plots/20260315_213225_stage_c0_10_C010_B_near_forward_01_grade_weights.png`
- `plots/20260315_213225_stage_c0_10_C010_B_near_reverse_10_topology_trajectory.png`
- `plots/20260315_213225_stage_c0_10_C010_B_near_reverse_10_flow_concentration_trace.png`
- `plots/20260315_213225_stage_c0_10_C010_B_near_reverse_10_grade_exchange_trace.png`
- `plots/20260315_213225_stage_c0_10_C010_B_near_reverse_10_grade_weights.png`
- `plots/20260315_213225_stage_c0_10_C010_B_near_shifted_12_topology_trajectory.png`
- `plots/20260315_213225_stage_c0_10_C010_B_near_shifted_12_flow_concentration_trace.png`
- `plots/20260315_213225_stage_c0_10_C010_B_near_shifted_12_grade_exchange_trace.png`
- `plots/20260315_213225_stage_c0_10_C010_B_near_shifted_12_grade_weights.png`
- `plots/20260315_213225_stage_c0_10_C010_C_incompatible_02_topology_trajectory.png`
- `plots/20260315_213225_stage_c0_10_C010_C_incompatible_02_flow_concentration_trace.png`
- `plots/20260315_213225_stage_c0_10_C010_C_incompatible_02_grade_exchange_trace.png`
- `plots/20260315_213225_stage_c0_10_C010_C_incompatible_02_grade_weights.png`
- `plots/20260315_213225_stage_c0_10_C010_C_incompatible_13_topology_trajectory.png`
- `plots/20260315_213225_stage_c0_10_C010_C_incompatible_13_flow_concentration_trace.png`
- `plots/20260315_213225_stage_c0_10_C010_C_incompatible_13_grade_exchange_trace.png`
- `plots/20260315_213225_stage_c0_10_C010_C_incompatible_13_grade_weights.png`
- `plots/20260315_213225_stage_c0_10_C010_C_randomized_control_topology_trajectory.png`
- `plots/20260315_213225_stage_c0_10_C010_C_randomized_control_flow_concentration_trace.png`
- `plots/20260315_213225_stage_c0_10_C010_C_randomized_control_grade_exchange_trace.png`
- `plots/20260315_213225_stage_c0_10_C010_C_randomized_control_grade_weights.png`
- `plots/20260315_213225_stage_c0_10_topology_vs_address_heatmap.png`
- `plots/20260315_213225_stage_c0_10_persistence_distribution_panel.png`
- `plots/20260315_213225_stage_c0_10_address_selectivity_panel.png`
