# Integrated Proto-Particle Feasibility Run

## Objective

Test one self-contained integrated feasibility question on top of the frozen Phase V--IX contracts: whether the branch keeps coherent large-scale scaling, supports at least one localized excitation candidate, and preserves that candidate under deterministic perturbations better than a deterministic altered-connectivity control.

## Frozen Inputs

- Phase VI operator hierarchy: `phase6-operator/phase6_operator_manifest.json`
- Phase VIII short-time trace contract: `phase8-trace/phase8_trace_manifest.json`
- Phase IX stabilized descriptor contract: `phase9-invariants/phase9_manifest.json`

## Scaling Coherence

- Branch extension reuses `n_side = [12, 24, 36, 48, 60]` and the frozen short-time window `[0.02, 0.03, 0.05, 0.075, 0.1]`.
- Branch short-window trace prediction error at the extension level is `0.00031931043` under the minimal `h^2` predictor.
- Branch ratio span across the extended hierarchy is `0.027957610451`, with R5 ratios `{'b1_over_b0': -3.866498178779, 'b2_over_b0': 7.207694668199}`.
- Control scaling remains deterministic but descriptor-shifted, with R5 ratios `{'b1_over_b0': -5.150258900129, 'b2_over_b0': 11.442568506369}`.

## Localized Excitation Scan

- Deterministic families scanned: `['low_mode_localized_wavepacket', 'mid_spectrum_localized_superposition', 'random_localized_seed']` on levels `[48, 60]` with tau grid `[0.0, 0.1, 0.2, 0.4, 0.8]`.
- Branch candidates selected at `n_side = 60`: `['low_mode_localized_wavepacket', 'mid_spectrum_localized_superposition']`.
- The strongest branch candidate is `low_mode_localized_wavepacket`, with aggregate persistent classification on both repeat levels `[48, 60]`.

## Persistence and Control

- Branch persistent candidates at the scan level: `['low_mode_localized_wavepacket']`.
- Control persistent candidate count at the scan level: `0`.
- Branch stable candidate set across refinement: `['low_mode_localized_wavepacket']`.
- The spectral truncation shift uses the deterministic cutoff `0.5`, and the control branch shows broader width growth and weaker concentration retention under the same evolution window.

## Bounded Interpretation

This integrated run supports a narrow feasibility statement only: within the frozen operator and spectral-trace contracts, one localized branch excitation family survives deterministic evolution and perturbation better than the altered-connectivity control while the large-scale trace descriptors stay coherent. No particle claim, continuum claim, or physical interpretation is asserted.

Integrated run establishes proto-particle feasibility for the frozen operator hierarchy.
