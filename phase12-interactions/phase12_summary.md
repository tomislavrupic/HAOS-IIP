# Phase XII - Two-Mode Interaction and Identity Preservation

## Objective

Test whether two copies of the frozen persistent localized candidate preserve identity and produce reproducible interaction classes under deterministic placement sweeps on the frozen branch hierarchy, with direct comparison against the deterministic control hierarchy.

## Frozen Inputs

- Candidate: `low_mode_localized_wavepacket` from Phase XI and the proto-particle manifest
- Refinement levels: `[48, 60, 72, 84]` with `h = 1 / n_side`
- Separation grid (lattice units): `[4, 8, 12, 16, 20, 24]`
- Placement schedules: `['axis_x_symmetric', 'diagonal_symmetric']`
- Evaluation tau: `0.8`
- Persistence tau grid: `[0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 3.6, 4.2]`
- Phase VIII short-time window retained by contract: `[0.02, 0.03, 0.05, 0.075, 0.1]`

## Key Results

- Branch survival outcomes: `26` of `48` scanned branch placements.
- Control survival outcomes: `26` of `48` scanned control placements.
- Branch merger outcomes: `16`; control merger outcomes: `7`.
- Branch decoherence outcomes: `6`; control decoherence outcomes: `15`.
- Mean branch successive regime-map distance: `0.194444444444`.
- Mean branch schedule distance: `0.083333333333`.
- Mean branch/control regime-map distance: `0.1875`.
- Branch onset physical-threshold relative span: `0.400000000002`.
- Control onset physical-threshold relative span: `0.400000000002`.
- Max branch survival-band size: `4` separations.
- Mean control survival-band size: `3.25` separations.
- Mean branch/control norm-retention gap at evaluation tau: `0.024632172441`.
- Mean control minus branch low-mode ratio gap at evaluation tau: `0.040266723284`.

## Bounded Interpretation

The branch is treated as feasible only if its interaction outcome maps remain reproducible across refinement and schedule, preserve identity over a non-trivial separation band, and remain distinguishable from the deterministic control. These results do not assert particles, forces, fields, or any physical interaction law.

Phase XII establishes two-mode interaction and identity-preservation feasibility for the frozen operator hierarchy.
