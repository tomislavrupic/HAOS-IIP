# Phase XIII - Multi-Mode Statistical Sector Formation

Phase XIII tests whether dilute populations of the frozen localized candidate develop reproducible survival and spacing structure on the frozen branch hierarchy, with direct comparison against the deterministic altered-connectivity control.

Run:

```bash
python3 phase13-sector-formation/build_multi_mode_sector.py
python3 phase13-sector-formation/diagnostics/check_phase13_bundle.py
```

Artifacts:

- `runs/phase13_population_survival_ledger.csv`
- `runs/phase13_spacing_statistics_ledger.csv`
- `runs/phase13_cluster_metrics.csv`
- `runs/phase13_spectral_ensemble_proxy.csv`
- `runs/phase13_runs.json`
- `phase13_summary.md`
- `phase13_manifest.json`
- `plots/phase13_survival_fraction_vs_time.svg`
- `plots/phase13_pair_distance_histogram_vs_refinement.svg`
- `plots/phase13_nearest_neighbor_scale_vs_refinement.svg`
- `plots/phase13_cluster_count_vs_time.svg`
- `plots/phase13_ensemble_spectral_proxy_vs_population.svg`

Constraint boundary:

- Uses only the frozen Phase XI persistent candidate.
- Uses only the frozen branch-local operator hierarchy and the deterministic control hierarchy already validated in earlier phases.
- Reuses the frozen Phase VIII and Phase XI operational windows.
- Does not assert thermodynamics, particles, gases, or kinetic laws.
