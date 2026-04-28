# Microtubule Lattice Demo

This is Biology Line B.

It is an experimental application layer on top of frozen HAOS-IIP. It models a microtubule-inspired cylindrical lattice as a deterministic interaction graph and tests whether lightweight HAOS-style diagnostics detect early degradation under lattice perturbation.

The model uses the minimal grounded background that microtubules are cylindrical polymers made of tubulin dimers, and that a common eukaryotic microtubule has 13 protofilaments arranged laterally into a hollow tube. Local longitudinal and lateral interactions stabilize the coarse lattice.

This demo does not modify HAOS core. It makes no claims about consciousness, Orch-OR, quantum computation, or HAOS explaining microtubules.

Run from the repository root:

```bash
python experiments/biology/microtubule_lattice_demo/run_microtubule_lattice_demo.py
```

Outputs are written to `experiments/biology/microtubule_lattice_demo/outputs/`.

## Spectral Telemetry Prototype

This folder also contains an experimental spectral-dynamics telemetry bridge:

```bash
uv run --with numpy --with matplotlib python experiments/biology/microtubule_lattice_demo/run_microtubule_spectral_telemetry.py
```

The spectral prototype uses a `13 x 30` cylindrical lattice, a toy GTP-cap-like
reference field, spectral address dynamics, and the stricter null ladder from
the HAOS dynamics sidecar. It reports recoverability, protofilament identity
retention, delta persistence, safety margin, and control comparison. Outputs are
ignored under `spectral_outputs/`.

This is still a toy telemetry probe, not biological validation.
