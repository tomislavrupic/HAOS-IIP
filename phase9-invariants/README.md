# Phase IX Coefficient Stabilization and Invariant Tracking

This phase tests whether the coefficient-like structures extracted from the frozen Phase VIII short-time trace window stabilize across refinement and whether simple normalized descriptors behave invariant-like on the frozen branch.

The phase consumes only frozen authority artifacts:

- `phase6-operator/phase6_operator_manifest.json`
- `phase7-spectral/phase7_spectral_manifest.json`
- `phase8-trace/phase8_trace_manifest.json`
- `PHASE_VIII_TRACE_CONTRACT_NOTE.md`

Run:

```bash
python3 phase9-invariants/build_coefficient_stability.py
```

Artifacts:

- `runs/phase9_coefficient_stability_ledger.csv`
- `runs/phase9_invariant_ratios_ledger.csv`
- `phase9_summary.md`
- `phase9_manifest.json`
- `plots/phase9_coefficient_vs_refinement.svg`
- `plots/phase9_ratio_vs_refinement.svg`
- `plots/phase9_rescaled_invariant_vs_refinement.svg`
- `plots/phase9_refinement_distance_convergence.svg`
- `runs/phase9_runs.json`
