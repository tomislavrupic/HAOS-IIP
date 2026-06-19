# Literature Component Bridge

This sidecar turns the bridge-literature extraction into an executable
component calibration.

It combines:

- spectral identity diagnostics,
- Hodge / cochain sector diagnostics,
- graph curvature diagnostics,
- transport metrics,
- kernel two-sample diagnostics.

The input is a deterministic synthetic relational graph plus matched controls.
The output is a component-level bridge instrument, not a physical bridge claim.

## Run

```bash
uv run python experiments/physics_bridge/literature_component_bridge/run_literature_component_bridge.py
```

## Outputs

- `precommitment_contract.json`
- `source_manifest.json`
- `component_scores.csv`
- `control_results.csv`
- `literature_component_bridge_report.md`
- `literature_component_bridge_result.json`

## Current Result

The current frozen run reports:

- `PASS` for component calibration,
- `LABEL_INVARIANCE_PASS`,
- `COMPONENT_CONTROLS_PASS`,
- `PHYSICAL_BRIDGE_NOT_ESTABLISHED`,
- `CLAIM_GATED_OPERATIONAL_MAPPING`.

## What PASS Means

`PASS` means the component instrument behaves as intended on the synthetic
calibration case:

- self-distance is zero,
- label-preserving permutations do not create false separation,
- weight-shuffle controls move weighted spectral/transport metrics,
- topology-destroyed controls move spectral, curvature, and transport metrics,
- triangle-removed controls move Hodge-sector metrics while graph-only metrics
  remain stable.

## What PASS Does Not Mean

This does not establish:

- physical correspondence,
- continuum spacetime,
- quantum correlations,
- gravitational dynamics,
- field theory,
- empirical validation.

It is only an operational bridge-method calibration built from the extracted
literature components.
