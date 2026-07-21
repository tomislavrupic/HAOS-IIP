# EL-R3-RT-02 Independent Restorative Dynamics

This experiment tests local cochain constraint-error feedback with redundant
one-bit edge-orientation memory. The mechanism sees current endpoint values,
stored relation signs, local incidence, and current inequality violations. It
does not see the unperturbed future trajectory, node-state checkpoint, recovery
score, or evaluation classification.

The primary target is `identity-preserving organizational-regime recovery`.
Exact microscopic reversal is not required.

Run tests before the frozen final execution:

```bash
uv run python -m unittest discover experiments/emergence_ladder/rung3_recovery_trajectory_v2/tests -v
```

Run the frozen calibration, validation, and final evaluation once:

```bash
uv run python experiments/emergence_ladder/rung3_recovery_trajectory_v2/run_experiment.py
```

Check the generated bundle:

```bash
uv run python experiments/emergence_ladder/rung3_recovery_trajectory_v2/check_bundle.py
```

This bundle does not alter EL-R3-RT-01, frozen HAOS metrics, or phase manifests.
It makes no operational-closure, agency, cross-scale, continuum, or physics
claim.

## Frozen Result

- Classification: `PARTIAL_RECOVERY_ONLY`
- Target recovery rate: `0.000000`
- Median state-region gain: `0.997665`
- Median edge-sign identity: `1.000000`
- Median functional recovery: `0.000000`
- Result hash: `el_r3_rt_02_cade4e057a54d119e6af143c`

The mechanism repairs relational constraints and preserves identity, but does
not satisfy the full multidimensional recovery contract. Rung 3 remains
unsupported.
