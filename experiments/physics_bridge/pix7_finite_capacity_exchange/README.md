# PIX-7 finite-capacity relational exchange

Internal speculative sidecar, created 2026-09-08 at the user's request to construct a candidate locally. No web search or originality certification. Existing HAOS-IIP results and public APIs are unchanged.

The model extends the relational Dirichlet energy shape with two unit-vector internal states per node and a load/mediator coupling. It uses the zero-phase incidence operator returned by `haos_core.build_graph`. The 2D periodic graph is an explicit supplied test substrate, not emergent spacetime or the canonical 3D physics bridge.

## Frozen first-run protocol

- J=K=1; coupling g=2; three sizes 4×4, 6×6, 8×8, fixed edge weight exp(-1/2).
- Seed a cosine contrast with amplitude 0.001. Evolve to dimensionless t=12.
- Compare to the derived tangent system, the g=0 control, and tighter integration tolerances.
- Stress with three random sphere-state seeds (7,19,41), plus near-pole initial data.
- Check weak-amplitude errors at 0.04, 0.02, 0.01; no parameter selection from outcomes.
- Check mass, energy, unit norms, relabeling symmetry and the tangent Poisson residual.
- No clipping, post-step renormalization, or damping. Initial vectors are normalized only to define initial data.
- Thresholds in `run_probe.py` check numerical consistency only. PASS_SCOPED is not a PIX-7 gravity-gate pass.

Run from repository root:

```bash
MPLCONFIGDIR=/tmp/pix7-mpl .venv/bin/python experiments/physics_bridge/pix7_finite_capacity_exchange/run_probe.py
```

Outputs: `outputs/results.json`, `outputs/trajectory.npz`, `outputs/comparison.png`.

The full derivation and gate failures are in `docs/notes/foundations/PIX_7_Constructed_Finite_Capacity_Exchange_v1.md`.
