# Control Integrity Audit

Implemented fact: this audit compares each matched control to the source target observation used by the CST toy benchmark.
Design choice: spectral density is represented only by the available low-mode summary proxy; no full eigenspectrum is reconstructed.

| Control | Preserved Properties | Destroyed Properties | Label Leakage | Validity Note |
| --- | --- | --- | --- | --- |
| degraded_signal_control | degree_distribution, event_chain_length, causal_depth | spectral_density_proxy, shell_ordering, invariant_marginals | high | component-specific validity only; preserves some tested signal |
| generic_graph_operator_control | spectral_density_proxy, shell_ordering, invariant_marginals | degree_distribution, event_chain_length, causal_depth | low | component-specific validity only; preserves some tested signal |
| label_permutation_control | degree_distribution, spectral_density_proxy, event_chain_length, causal_depth, shell_ordering, invariant_marginals | none | low | component-specific validity only; preserves some tested signal |
| parameter_matched_null_control | shell_ordering | degree_distribution, spectral_density_proxy, event_chain_length, causal_depth, invariant_marginals | low | component-specific validity only; preserves some tested signal |
| periodic_diagonal_augmented_control | degree_distribution, event_chain_length, partial:causal_depth | spectral_density_proxy, shell_ordering, invariant_marginals | high | strongest frozen matched control; still shares benchmark probe, seed grid, and address fields |
| randomized_edge_control | degree_distribution, spectral_density_proxy, event_chain_length, causal_depth, shell_ordering, invariant_marginals | none | high | invalid as a negative control for invariant, spectral, causal-depth, shell-ordering, and label-set-sensitive tests |
| shuffled_structure_control | degree_distribution, spectral_density_proxy, event_chain_length, causal_depth, shell_ordering, invariant_marginals | none | high | invalid as a negative control for invariant, spectral, causal-depth, shell-ordering, and label-set-sensitive tests |

Key hostile finding: randomized and shuffled controls retain the exact event-label set and most non-chain descriptor marginals.
A control that preserves a component's tested signal is not a valid negative control for that component.
