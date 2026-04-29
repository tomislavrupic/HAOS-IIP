# Phase 58 Spherical Harmonic Control Probe

This physics-bridge sidecar is the first bounded follow-up to the Phase 57
celestial boundary audit.

It tests a known target: low-order spherical harmonic mode organization on a
discrete sphere graph. It does not test celestial holography, BMS symmetry,
Virasoro structure, soft theorems, S-matrix recovery, or gravitational memory.

## Run

```bash
python3 experiments/physics_bridge/spherical_harmonic_control_probe/run_spherical_harmonic_probe.py
```

If your local `python3` does not have NumPy/Matplotlib:

```bash
uv run --with numpy --with matplotlib python experiments/physics_bridge/spherical_harmonic_control_probe/run_spherical_harmonic_probe.py
```

## Outputs

- `mode_diagnostics.csv`
- `bridge_status.json`
- `spherical_harmonic_probe.md`
- `figures/mode_band_overlap.png`
- `figures/refinement_drift.png`
- `figures/control_comparison.png`
- `figures/leakage_by_band.png`
- `figures/degeneracy_splitting.png`

## Controls

- `coordinate_permutation_control`: preserves the sphere graph and spectrum but
  breaks node-to-coordinate identity.
- `weight_shuffle_control`: preserves topology while destroying the
  distance-weight signature.
- `degree_rewire_control`: preserves degree-like graph statistics while breaking
  local spherical adjacency.
- `ring_smooth_control`: supplies a smooth low-spectrum graph that is not an S2
  boundary.

## Status Semantics

- `PASS`: the known sphere target is separable from controls on low-l subspace
  recovery, leakage, and organization score.
- `MARGINAL`: the probe detects spherical organization, but controls compete.
- `OPEN`: the telemetry cannot distinguish S2 mode organization from generic
  smoothness.

Even a `PASS` is only a boundary-geometry sanity check. It is not a celestial
holography result.
