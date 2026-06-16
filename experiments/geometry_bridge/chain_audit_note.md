# Geometry Chain Audit Note

Current read:

- Geometry layer exists as a blind synthetic recovery benchmark.
- Transformation semantics layer exists as a frozen transport/holonomy
  benchmark.
- Probability-rule layer exists as a frozen observable-to-probability bridge.
- Observable prediction layer exists as a synthetic forecast on a held-out
  target.
- Chain audit layer now exists and freezes the dependency chain explicitly.

What is still open:

- the chain is synthetic calibration, not a Bell derivation
- the observable layer is informative at the pairwise level but still weak at
  the coarse class level
- none of the layers authorize a physical claim
- hidden geometry now transfers distance, orientation, and held-out relations,
  but transformation semantics still needed its own recovery rung
- hidden transformation recovery now exists as GEO-T1-01 and closes the
  algebra-of-change gap for the synthetic chain

Reproduction commands:

```bash
uv run python -m unittest discover experiments/geometry_bridge/synthetic_intrinsic_geometry_recovery/tests
uv run python -m unittest discover experiments/geometry_bridge/synthetic_transformation_semantics_recovery/tests
uv run python -m unittest discover experiments/geometry_bridge/synthetic_probability_rule_recovery/tests
uv run python -m unittest discover experiments/geometry_bridge/synthetic_observable_prediction/tests
uv run python experiments/geometry_bridge/synthetic_intrinsic_geometry_recovery/run_synthetic_geometry_recovery.py
uv run python experiments/geometry_bridge/synthetic_transformation_semantics_recovery/run_synthetic_transformation_semantics.py
uv run python experiments/geometry_bridge/synthetic_probability_rule_recovery/run_synthetic_probability_rule_recovery.py
uv run python experiments/geometry_bridge/synthetic_observable_prediction/run_synthetic_observable_prediction.py
uv run python experiments/geometry_bridge/run_geometry_chain_audit.py
```

Evidence pointers:

- [Geometry README](./README.md)
- [Geometry chain requirements](./geometry_chain_requirements.md)
- [Objective completion audit](./goal_completion_audit.md)
- [Geometry transfer summary](./geometry_transfer_summary.md)
- [Geometry transfer consistency](./geometry_transfer_consistency.md)
- [Geometry transfer ledger](./geometry_transfer_ledger.md)
- [Geometry gap map](./geometry_gap_map.md)
- [Geometry evidence manifest](./geometry_evidence_manifest.md)
- [Geometry current state snapshot](./geometry_current_state_snapshot.md)
- [Hidden transformation recovery README](./synthetic_hidden_transformation_recovery/README.md)
- [Hidden transformation result](./synthetic_hidden_transformation_recovery/hidden_transformation_result.json)
- [Hidden probability-rule recovery README](./synthetic_hidden_probability_rule_recovery/README.md)
- [Hidden probability-rule result](./synthetic_hidden_probability_rule_recovery/probability_rule_result.json)
- [GEO-01 result](./synthetic_intrinsic_geometry_recovery/geometry_recovery_result.json)
- [GEO-02 result](./synthetic_transformation_semantics_recovery/semantics_result.json)
- [GEO-03 result](./synthetic_probability_rule_recovery/probability_rule_result.json)
- [GEO-04 result](./synthetic_observable_prediction/observable_prediction_result.json)
- [GEO-P1-01 result](./synthetic_hidden_probability_rule_recovery/probability_rule_result.json)
- [GEO-05 result](./geometry_chain_audit_result.json)

Bottom line:

The requested chain is implemented only as synthetic calibration. The
observable rung now clears its frozen synthetic holdout, but the chain-level
audit still keeps the missing dependency explicit rather than silently
upgrading it into a Bell mechanism.
