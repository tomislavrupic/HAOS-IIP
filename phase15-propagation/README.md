# Phase XV — Effective Propagation Structure and Causal-Scale Feasibility

Phase XV tests whether the frozen dilute collective sector from Phase XIII and the frozen collective evolution policy from Phase XIV support bounded propagation-level structure.

This phase is constrained by the existing frozen contracts:

- Phase V recovery readout definitions
- Phase VI frozen operator hierarchy `delta_h`
- Phase VII spectral feasibility bundle
- Phase VIII short-time trace window
- Phase IX coefficient stabilization bundle
- Phase X cautious continuum-bridge bundle
- Phase XI survival-basin bundle
- Phase XII interaction bundle
- Phase XIII multi-mode sector bundle
- Phase XIV collective dynamics bundle

Phase XV does not introduce new candidate objects, change normalization, or modify any earlier phase.

## Frozen Inputs

- Candidate label: `low_mode_localized_wavepacket`
- Refinement levels: `48, 60, 72, 84`
- Ensemble sizes: `3, 5, 7`
- Layout seeds: `1301, 1302, 1303`
- Tau grid: `0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 3.6, 4.2`
- Minimum physical separation: `0.12`
- Control hierarchy: `periodic_diagonal_augmented_control`

## Propagation Probes

Phase XV applies three deterministic disturbances on top of the frozen Phase XIV collective state:

1. `density_pulse`
2. `candidate_removal`
3. `bias_onset`

For each disturbance and each hierarchy label (`frozen_branch`, `periodic_diagonal_augmented_control`), the builder records:

- disturbance radius over time
- threshold-crossing latencies by candidate distance
- effective propagation speed estimates
- long-wavelength response descriptors from the density-profile response

## Outputs

- `phase15_propagation_ledger.csv`
- `phase15_effective_speed_ledger.csv`
- `phase15_influence_range_ledger.csv`
- `phase15_transport_descriptor_ledger.csv`
- `phase15_runs.json`
- `phase15_manifest.json`
- `phase15_summary.md`

Plots are written to `phase15-propagation/plots/`.

## Claim Boundary

Phase XV is limited to effective propagation-structure feasibility diagnostics on the frozen collective sector. It does not claim metric structure, spacetime, light-cones, continuum field equations, or physical correspondence.
