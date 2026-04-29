# HAOS-IIP Prediction Workbench

This directory is a bounded prediction layer for HAOS-IIP.

It exists to turn HAOS-style intuitions into falsifiable prediction records:

- What system is being tested?
- What perturbation is applied?
- What recoverability observable should move first?
- What standard or visible failure marker should occur later?
- What would falsify the prediction?

This is not a new-physics claim engine. It is a disciplined place to state candidate predictions that may be absent from current field practice, then force them through measurable pass/fail criteria.

## Boundary

The workbench does not modify:

- `haos_core/`
- `telemetry/`
- frozen phase artifacts
- canonical notes
- previous experiment outputs

The workbench can reference prior artifacts as evidence, but it does not rewrite them.

## Files

- `PREDICTION_PROTOCOL.md` defines admission rules and evidence levels.
- `prediction_ledger.json` stores machine-readable prediction records.
- `validate_predictions.py` checks the ledger and writes compact summaries.
- `templates/prediction_record_template.md` provides the manual record shape.
- `outputs/` contains generated validation summaries.

## Run

```bash
python3 predictions/validate_predictions.py
```

Expected result:

```text
prediction_records: 5
validation_passed: True
outputs_written: predictions/outputs/
```

## Safe Claim

HAOS-IIP may propose testable recoverability predictions that are not standard observables in current physics, biology, or network fields.

It does not claim those predictions are true until they survive external perturbation tests, independent controls, and failure attempts.

