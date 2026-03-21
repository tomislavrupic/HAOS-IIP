# Phase XII - Two-Mode Interaction and Identity Preservation

Phase XII tests whether two copies of the frozen persistent branch candidate preserve identity and produce structured interaction outcomes under deterministic placement sweeps.

Run:

```bash
python3 phase12-interactions/build_two_mode_interaction.py
python3 phase12-interactions/diagnostics/check_phase12_bundle.py
```

Artifacts:

- `runs/phase12_interaction_outcome_ledger.csv`
- `runs/phase12_identity_metrics_ledger.csv`
- `runs/phase12_interaction_thresholds.csv`
- `runs/phase12_runs.json`
- `phase12_summary.md`
- `phase12_manifest.json`
- `plots/phase12_persistence_time_vs_separation.svg`
- `plots/phase12_identity_metric_vs_time.svg`
- `plots/phase12_interaction_regime_map.svg`
- `plots/phase12_threshold_scaling_vs_refinement.svg`

Constraint boundary:

- Reuses the frozen `low_mode_localized_wavepacket` candidate definition from Phase XI.
- Reuses the frozen `delta_h` hierarchy and Phase VIII / Phase XI probe windows.
- Does not claim particles, forces, fields, or scattering laws.
