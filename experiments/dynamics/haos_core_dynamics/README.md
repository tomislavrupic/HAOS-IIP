# HAOS Core Dynamics

This is a bounded sidecar dynamics rig for HAOS-IIP. It does not modify frozen
core APIs, phase artifacts, or canonical notes.

The first probe asks:

> Can a native HAOS-derived interaction graph recover a perturbed scalar field
> under local operator dynamics better than deterministic controls?

This is not yet full branch-local cochain-Laplacian hierarchy evolution. It is
a minimal native dynamics falsifier: graph operator, local propagation,
recoverability, invariant drift, causal spread, and controls.

## Run

```bash
uv run --with numpy --with matplotlib python experiments/dynamics/haos_core_dynamics/run_haos_dynamics.py
```

Useful options:

- `--steps`: number of discrete dynamics steps.
- `--mode scalar|rd`: scalar recovery dynamics or graph reaction-diffusion.
- `--dt`: time-step size.
- `--diffusion`: Laplacian-flow strength.
- `--recovery-gain`: attraction back toward the frozen reference field.
- `--rd-du`, `--rd-dv`, `--rd-feed`, `--rd-kill`: Gray-Scott-style graph
  reaction-diffusion parameters.
- `--perturbation-scale`: controlled local perturbation size.
- `--damage-fraction`: optional edge damage applied after the baseline graph is
  built.
- `--seed`: deterministic seed.

## Outputs

Generated outputs are ignored:

- `bridge_status.json`
- `dynamics_report.md`
- `dynamics_timeseries.csv`
- `recoverability_timeseries.png`
- `invariant_drift.png`
- `control_comparison.png`
- `observed_final_pattern.png`

Reaction-diffusion run:

```bash
uv run --with numpy --with matplotlib python experiments/dynamics/haos_core_dynamics/run_haos_dynamics.py \
  --mode rd \
  --steps 300 \
  --dt 0.4
```

## Status Semantics

- `PASS`: observed recoverability is high, invariant drift is bounded, and all
  controls are worse by the configured margin.
- `MARGINAL`: observed dynamics recover or preserve structure, but control
  contrast is incomplete.
- `FAIL`: controls match/exceed observed recovery, recoverability collapses, or
  invariant drift is too large.
- `OPEN_NO_DATA_SYNTHETIC`: no HAOS geometry config is available; the run is a
  plumbing check only.

These are dynamics diagnostics, not physics claims.
