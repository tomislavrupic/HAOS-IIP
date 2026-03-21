# Phase 5 Readout

Phase V is a deterministic recovery-statistics rig over frozen Phase IV sector bundles.

- Scope: recovery histogram and scaling readout only
- Runner: `runs/run_phase5.py`
- Diagnostics module: `diagnostics/check_phase5_bundle.py`
- Modes: `frozen_sector`, `dummy_sector`
- Control classes: `stable_frozen_sector`, `degraded_sector_control`, `shuffled_null_control`
- Imports: `haos_core` only, no phase imports
- Outputs: `phase5_readout_bundle_latest.json`, `phase5_recovery_histogram_latest.json`, `phase5_scaling_fit_latest.json`
- Audit outputs: `phase5_reproducibility_ledger_latest.json`, `phase5_robustness_audit_latest.json`, `phase5_control_discrimination_latest.json`
- Authority freeze: `phase5_authoritative_summary.md`, `phase5_authoritative_manifest.json`, `phase5_summary_figures/`, `phase5_runs_ledger_latest.json`

Canonical Phase V input schema:

```json
{
  "sector_identifier": "...",
  "graph_reference_or_seed": {},
  "kernel_configuration_snapshot": {},
  "selector_configuration": {},
  "invariant_baseline_vector": {},
  "perturbation_policy_descriptor": {}
}
```

Frozen-sector mode accepts only a compact Phase IV bundle and rejects oversized payloads. The runner compresses that bundle into the schema above, samples incidence perturbations, computes one structural recovery score per trial, builds the histogram, and labels the histogram form descriptively.

Phase V-B uses the same runner with `--audit phase5b`. It performs:

- a 5-schedule reproducibility sweep for `stable_frozen_sector`
- narrow one-at-a-time policy variation over `flip_probability`, `max_phase_drift_fraction_of_pi`, and `amplitude_jitter`
- matched frozen-sector control discrimination across `stable_frozen_sector`, `degraded_sector_control`, and `shuffled_null_control`

All audit runs stay inside `phase5-readout/`, keep the existing score and schema unchanged, and write compact JSON-only ledgers.

Phase V-close uses the same runner with `--audit phase5close`. It refreshes the stable frozen baseline, regenerates the Phase V-B audit outputs, writes the authority summary and manifest, records the edge-case classifier flip explicitly, and emits three minimal SVG figures.
