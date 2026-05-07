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
- `ferron_extension_testable_signatures.md` states the ferron-derived
  field-facing predictions in short-note form.
- `ferron_extension_testable_signatures.pdf` is the compiled short note.
- `magnon_address_stability_testable_signatures.md` states the spin-sector
  prediction layer seeded by ferron calibration and a user-provided alpha-Fe2O3
  audit anchor.
- `magnon_address_stability_testable_signatures.pdf` is the compiled short note.
- `validate_predictions.py` checks the ledger and writes compact summaries.
- `templates/prediction_record_template.md` provides the manual record shape.
- `outputs/` contains generated validation summaries.

## Run

```bash
python3 predictions/validate_predictions.py
```

Expected result:

```text
prediction_records: 13
validation_passed: True
outputs_written: predictions/outputs/
```

## Ferron Extension

The ferron extension records four bounded, falsifiable signatures seeded by
Materials Bridge Line A:

- polar-mode address selection in layered ferroelectrics
- thickness / boundary dependence of target-band propagation
- target-window re-addressing after secondary perturbation
- cross-dataset narrowband control discrimination

These are field-facing diagnostic candidates, not proof claims. The raw
STFT/time-frequency grid boundary from the ferron sidecar remains active:
`NO_STFT_DATA_FOUND`.

## Magnon Extension

The magnon extension records four bounded, falsifiable spin-sector signatures:

- crystal-orientation-dependent narrowband magnon address
- thickness / boundary-dependent group velocity
- target-address recovery after secondary perturbation
- hybrid magnon-phonon / cross-sector locking

These are field-facing diagnostic candidates. The user-provided alpha-Fe2O3
audit anchor is recorded as `persistence_proxy=0.902`, but the corresponding
local audit artifact was not present in this repo at note creation time.

## Safe Claim

HAOS-IIP may propose testable recoverability predictions that are not standard observables in current physics, biology, or network fields.

It does not claim those predictions are true until they survive external perturbation tests, independent controls, and failure attempts.
