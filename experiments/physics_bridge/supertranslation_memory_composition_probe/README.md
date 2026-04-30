# Phase 63 Supertranslation + Memory Composition Toy Probe

This sidecar composes the previous celestial toy probes into a stricter
two-stage benchmark:

1. a synthetic multi-pole supertranslation-like shift on a sampled `S2`;
2. a later synthetic memory response induced by that shift through a fixed toy
   harmonic response map.

The narrow question:

> Can HAOS-style spectral telemetry recover shift, memory, temporal order, and
> the induced shift-to-memory relation better than controls that preserve only
> pieces of the construction?

It does **not** claim BMS supertranslation recovery, BMS charge recovery, real
gravitational memory, soft-hair detection, celestial holography, or S-matrix
reconstruction.

## Run

```bash
python3 experiments/physics_bridge/supertranslation_memory_composition_probe/run_supertranslation_memory_composition_probe.py
```

If the default Python runtime lacks plotting dependencies:

```bash
uv run --with numpy --with matplotlib python experiments/physics_bridge/supertranslation_memory_composition_probe/run_supertranslation_memory_composition_probe.py
```

## Outputs

- `composition_diagnostics.csv`
- `bridge_status.json`
- `supertranslation_memory_composition_report.md`
- `figures/control_comparison.png`
- `figures/metric_breakdown.png`
- `figures/composition_projection_traces.png`
- `figures/composition_relation.png`
- `figures/field_identity.png`

## Claim Gate

`PASS` means only that the synthetic two-stage target is separable from the
included toy controls. The established physics claims around BMS
supertranslations, BMS charges, gravitational memory, celestial holography, and
soft theorems remain `OPEN` under the Phase 60 claim gate.
