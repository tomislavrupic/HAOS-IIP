# Phase 62 Multi-Pole Supertranslation Toy Probe

This sidecar extends the Phase 61 toy memory benchmark from one permanent
deformation to an ordered sequence of synthetic `l=2,3,4` mode shifts on a
sampled `S2` boundary.

The narrow question:

> Can HAOS-style spectral telemetry distinguish a structured multi-pole,
> supertranslation-like toy memory function from smooth decoys, wrong harmonic
> addresses, wrong event ordering, transient bursts, and stochastic drift?

It does **not** claim BMS supertranslation recovery, BMS charge recovery, real
gravitational memory, celestial holography, soft hair detection, or S-matrix
reconstruction.

## Run

```bash
python3 experiments/physics_bridge/multipole_supertranslation_probe/run_multipole_supertranslation_probe.py
```

If the default Python runtime lacks plotting dependencies:

```bash
uv run --with numpy --with matplotlib python experiments/physics_bridge/multipole_supertranslation_probe/run_multipole_supertranslation_probe.py
```

## Outputs

- `multipole_diagnostics.csv`
- `bridge_status.json`
- `multipole_supertranslation_report.md`
- `figures/control_comparison.png`
- `figures/metric_breakdown.png`
- `figures/event_projection_traces.png`
- `figures/band_power_signature.png`
- `figures/field_identity.png`

## Claim Gate

`PASS` means only that the synthetic ordered multi-pole target is separable from
the included toy controls. The established physics claims around BMS
supertranslations, gravitational memory, celestial holography, and soft theorems
remain `OPEN` under the Phase 60 claim gate.
