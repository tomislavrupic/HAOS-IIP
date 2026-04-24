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
