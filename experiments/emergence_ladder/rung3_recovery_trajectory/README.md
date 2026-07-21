# EL-R3-RT-01 Recovery Trajectory

This experiment distinguishes persistence from post-perturbation recovery. It
uses the frozen public `build_graph` and `run_transport` APIs, compares a
perturbed trajectory with its unperturbed reference, and executes passive,
operator-only, topology-altered, and trivial-attractor controls.

Run:

```bash
uv run python experiments/emergence_ladder/rung3_recovery_trajectory/run_recovery_trajectory.py
```

Current classification: `NEGATIVE_RESULT`.

The grade-locked mechanism does not return toward the reference trajectory and
does not preserve the declared grade signature. This is a bounded synthetic
operator result, not a claim about universal recoverability or physics.
