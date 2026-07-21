# EL-R4-OC-01 Operational Closure

This versioned experiment tests whether an observed internal event-transition
context predicts the next event in frozen Phase XVI chains better than external
metadata and position-only baselines.

The contract was frozen before execution. Construction uses refinement levels
48 and 60, validation uses 72, and holdout uses 84. Transition-shuffled
training, context ablation, and target-label permutation are routed through the
same prediction path.

Run:

```bash
uv run python experiments/emergence_ladder/rung4_operational_closure/run_operational_closure.py
```

Current classification: `NEGATIVE_RESULT`.

The aggregate internal advantage passes its numerical margin, but the run-level
95% interval includes zero. This closes the compressed event-chain mechanism.
It does not establish operational, causal, thermodynamic, biological,
cross-scale, or physical closure.
