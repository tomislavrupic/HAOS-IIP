# Synthetic Relational Geometry Calibration

This sidecar calibrates the semantic/refinement diagnostics on a synthetic
domain where the relation graph is known by construction.

It is not a physics bridge, not a Bell experiment, and not evidence that
HAOS-IIP derives quantum correlations.

## Run

```bash
uv run python experiments/physics_bridge/synthetic_relational_geometry_calibration/run_synthetic_calibration.py
```

## Conditions

- `known_positive`: known semantic relation graph with coherent refinements.
- `semantic_permutation_control`: destroys label semantics while preserving
  refinement stability.
- `topology_destroyed_control`: destroys the relation topology.
- `refinement_broken_control`: destroys cross-refinement continuity.
- `ambiguous_partial_preservation`: partial preservation expected to remain
  `OPEN`.

## Outputs

- `precommitment_contract.json`
- `synthetic_source_manifest.json`
- `semantic_relation_matrix.csv`
- `calibration_geometry_matrix.csv`
- `calibration_control_results.csv`
- `calibration_refinement_persistence.csv`
- `synthetic_calibration_report.md`
- `synthetic_calibration_result.json`

## Boundary

Passing this calibration means the measurement instrument can distinguish known
synthetic semantic/refinement structure from declared controls. It does not
authorize any Bell or physics claim.
