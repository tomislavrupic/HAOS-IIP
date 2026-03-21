# Phase VIII Trace Asymptotic Feasibility

This phase tests whether the frozen branch-local cochain-Laplacian hierarchy supports a reproducible short-time trace scaling regime across refinement.

The phase consumes only the frozen Phase V, Phase VI, and Phase VII manifests. It does not modify the operator hierarchy, the refinement definition, or the Phase V recovery observable.

Run:

```bash
python3 phase8-trace/build_trace_asymptotic_feasibility.py
```

Artifacts:

- `runs/phase8_trace_scaling_ledger.csv`
- `phase8_trace_summary.md`
- `phase8_trace_manifest.json`
- `plots/phase8_trace_vs_t.svg`
- `plots/phase8_log_trace_slope_vs_t.svg`
- `plots/phase8_coefficient_like_terms_vs_refinement.svg`
- `runs/phase8_runs.json`
