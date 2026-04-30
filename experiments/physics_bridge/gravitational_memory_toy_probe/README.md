# Phase 61 Gravitational Memory Toy Probe

This sidecar is a bounded celestial-facing toy benchmark. It asks whether
HAOS-style spectral telemetry can detect a synthetic permanent displacement
pattern on a sampled `S2` boundary better than controls that preserve smoothness,
time-series scale, or low-mode structure while breaking memory identity.

It does **not** claim real gravitational-memory recovery, BMS charge recovery,
soft-hair detection, celestial holography, or S-matrix reconstruction.

## Run

```bash
python3 experiments/physics_bridge/gravitational_memory_toy_probe/run_gravitational_memory_toy_probe.py
```

If the default Python runtime lacks plotting dependencies:

```bash
uv run --with numpy --with matplotlib python experiments/physics_bridge/gravitational_memory_toy_probe/run_gravitational_memory_toy_probe.py
```

## Outputs

- `memory_diagnostics.csv`
- `bridge_status.json`
- `gravitational_memory_toy_report.md`
- `figures/control_comparison.png`
- `figures/metric_breakdown.png`
- `figures/memory_projection_trace.png`
- `figures/memory_field_identity.png`

## Claim Gate

`PASS` means only that the synthetic target is separable from the included toy
controls. The established physics claim `gravitational_memory_observable` remains
`OPEN` under the Phase 60 claim gate.
