# Continuum Sketch Summary

## Inputs

- Refinement slice: `n_side = {48, 60, 72, 84}`
- Branch seeds: `1302, 1303`
- Matched control hierarchy: `periodic_diagonal_augmented_control`
- Disturbance family for propagation / ordering / shell metrics: `bias_onset`

## Observable Conventions

- Low-mode dispersion proxy: frozen `spectral_gap` from the Phase XI persistence ledger.
- Effective propagation speed band: mean `effective_speed` on the Phase XV `bias_onset` / ensemble-7 / seed-pair slice.
- Temporal-ordering consistency score: `1 - event_time_distance / 4.2` on the Phase XVI primary event ledger for the same seed pair.
- Distance-surrogate shell slope: the Phase XVIII least-squares arrival-vs-depth slope, extended to `n_side = 48` using the same artifact-only reconstruction rule.

## Feasibility Checks

- Smooth branch trends identified: `3` of `5` observables.
- Branch propagation-speed coefficient of variation: `0.027656` vs control `0.051604`.
- Branch shell-slope max adjacent drift: `0.085714` vs control `0.584615`.
- Branch persistence-time band: `1.200000` to `3.600000`.
- Branch temporal-ordering score band: `0.780045` to `0.845805`.

## Fit Notes

- `dispersion_proxy`: branch relative span `1.115610`, control relative span `1.115610`, branch best fit `R^2=1.000000`, control best fit `R^2=1.000000`, branch log-slope `1.996967`, control log-slope `1.996967`.
- `effective_speed_band`: branch relative span `0.076033`, control relative span `0.136723`, branch best fit `R^2=0.789029`, control best fit `R^2=0.527198`, branch log-slope `0.115737`, control log-slope `-0.174273`.
- `persistence_time`: branch relative span `1.043478`, control relative span `1.230769`, branch best fit `R^2=0.977398`, control best fit `R^2=0.996844`, branch log-slope `-1.868387`, control log-slope `-2.819219`.
- `ordering_consistency_score`: branch relative span `0.079835`, control relative span `0.200418`, branch best fit `R^2=0.924506`, control best fit `R^2=0.149767`, branch log-slope `-0.143710`, control log-slope `0.099396`.
- `distance_surrogate_shell_slope`: branch relative span `0.091898`, control relative span `0.581341`, branch best fit `R^2=0.780069`, control best fit `R^2=0.065301`, branch log-slope `-0.159002`, control log-slope `-0.276169`.

## Branch vs Control

- `dispersion_proxy` mean branch/control gap: `0.007217`.
- `effective_speed_band` mean branch/control gap: `0.011985`.
- `persistence_time` mean branch/control gap: `1.000000`.
- `ordering_consistency_score` mean branch/control gap: `0.066893`.
- `distance_surrogate_shell_slope` mean branch/control gap: `0.182463`.

The minimal scaling sketch indicates bounded refinement-ordered behavior of key frozen observables. This is consistent with a proto-continuum effective description within the tested branch regime. No continuum ontology or universality is asserted.
