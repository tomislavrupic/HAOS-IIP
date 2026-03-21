# Phase XIV - Collective Dynamics and Mesoscopic Transport

Phase XIV tests whether the frozen dilute sector from Phase XIII supports bounded collective dynamics descriptors without introducing new candidate objects or modifying earlier freezes.

## Frozen Inputs

- Candidate: `low_mode_localized_wavepacket`
- Refinement hierarchy: `n_side in [48, 60, 72, 84]`
- Ensemble sizes: `[3, 5, 7]`
- Layout seeds: `[1301, 1302, 1303]`
- Minimum physical separation: `0.12`
- Tau grid: `[0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 3.6, 4.2]`
- Phase VIII short-time window: `[0.02, 0.03, 0.05, 0.075, 0.1]`

## Probe Set

- Weak density-gradient transport probe
- Collective relaxation-timescale extraction
- Fluctuation-spectrum proxy on coarse x-bin density profiles
- Weak phase-coherent driving with bias removal
- Spacing-density equation-of-state proxy from the frozen seeded layouts

## Outputs

- `runs/phase14_transport_ledger.csv`
- `runs/phase14_relaxation_times.csv`
- `runs/phase14_fluctuation_spectra.csv`
- `runs/phase14_density_response.csv`
- `runs/phase14_equation_of_state_proxy.csv`
- `runs/phase14_runs.json`
- `phase14_manifest.json`
- `phase14_summary.md`

## Constraint

Phase XIV may only derive collective observables from previously frozen quantities: survival fraction, localization width, spectral ensemble proxy, pair-distance statistics, and identity metrics. It does not assert hydrodynamics, thermodynamics, or continuum medium laws.
