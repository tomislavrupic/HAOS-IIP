# Geometry Chain Audit

- bridge id: GEO-05
- status: OPEN

## Chain Readout
- intrinsic_geometry: synthetic=OPEN / external=OPEN / unresolved=holdout transfer does not outperform the best baseline on the spiral family (holdout_pass=False / best_baseline_spearman=0.8298584298584298)
- transformation_semantics: synthetic=VALID / external=OPEN (synthetic semantics layer present; frozen control summary exists but no Bell relevance)
- hidden_probability_rule: synthetic=VALID / external=OPEN (holdout_accuracy=0.7083333333333334 / null_accuracy=0.625 / brier=0.1963532621393027)
- observable_prediction: synthetic=VALID / external=OPEN (pairwise_accuracy=1.0 / pairwise_null_accuracy=0.5 / accuracy=0.8333333333333334 / null_accuracy=0.25)

## Status Axes
- synthetic_calibration: CLOSED
- external_semantics: OPEN

## Boundary
Synthetic calibration only; no Bell scoring; no physical mechanism claim.
