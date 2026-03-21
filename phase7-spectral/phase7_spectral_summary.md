# Phase VII Spectral Summary

Timestamp: `2026-03-21T14:55:04Z`.

## Track Selection

- Primary track: `spectral_feasibility` (Track A).
- Secondary track executed: `False`.
- Selection rationale: Phase VI already showed stable operator diagnostics across refinement, while Phase V contains a bounded robustness edge-case classifier flip that makes observable coupling the less clean primary authority statement.

## Frozen References

- Operator manifest: `phase6-operator/phase6_operator_manifest.json`.
- Phase V authority manifest: `phase5-readout/phase5_authoritative_manifest.json`.
- Frozen operator class: `cochain_laplacian`.
- Frozen stable Phase V class: `linear_like`.

## Per-Level Probe Ledger

| level | n_side | h | band 1 | band 2 | band 3 | radius | nullspace | trace(t=0.25) | trace(t=0.5) | trace(t=1.0) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| R1 | 12 | 0.083333333333 | 0.263337445092 | 0.526674890185 | 0.982788724620 | 7.862309796963 | 4 | 36.135312835762 | 32.768913274466 | 27.241830650454 |
| R2 | 24 | 0.041666666667 | 0.067853205626 | 0.135706411252 | 0.266788738673 | 7.965353020924 | 4 | 38.939069201182 | 37.917120873497 | 35.983284258013 |
| R3 | 36 | 0.027777777778 | 0.030325938407 | 0.060651876814 | 0.120382315335 | 7.984582776023 | 4 | 39.468496959778 | 38.947172738264 | 37.934120590456 |
| R4 | 48 | 0.020833333333 | 0.017091721482 | 0.034183442965 | 0.068074441836 | 7.991324152244 | 4 | 39.728129343185 | 39.458847149761 | 38.927927288909 |

## Low-Mode Scaling

- Band 1: exponent `1.972182228379`, R^2 `0.999979005926`, monotone `True`.
- Band 2: exponent `1.972182228365`, R^2 `0.999979005926`, monotone `True`.
- Band 3: exponent `1.924243971927`, R^2 `0.999835444451`, monotone `True`.
- Exponent span across the three fitted bands: `0.047938256452`.
- Minimum low-mode fit R^2: `0.999835444451`.

## Spectral Envelope and Kernel Stability

- Spectral radius values: `[7.862309796963, 7.965353020924, 7.984582776023, 7.991324152244]`.
- Spectral radius monotone: `True`.
- Nullspace estimate values: `[4, 4, 4, 4]`.
- Nullspace consistent: `True`.

## Normalized Density Proxy

- The density surrogate is the empirical CDF of the first 24 sampled positive low-mode eigenvalues after normalization by the first positive band.
- R1 vs R2: max normalized CDF distance `0.0`.
- R1 vs R3: max normalized CDF distance `0.041667`.
- R1 vs R4: max normalized CDF distance `0.0`.
- R2 vs R3: max normalized CDF distance `0.041667`.
- R2 vs R4: max normalized CDF distance `0.0`.
- R3 vs R4: max normalized CDF distance `0.041667`.
- Mean pairwise normalized CDF distance: `0.0208335`.
- Max pairwise normalized CDF distance: `0.041667`.

## Short-Time Trace Proxy

- The trace proxy is the truncated sum of `exp(-t lambda_i)` over the 40 sampled lowest modes.
- t_0_25: values `[36.135312835762, 38.939069201182, 39.468496959778, 39.728129343185]`, monotone `True`.
- t_0_5: values `[32.768913274466, 37.917120873497, 38.947172738264, 39.458847149761]`, monotone `True`.
- t_1_0: values `[27.241830650454, 35.983284258013, 37.934120590456, 38.927927288909]`, monotone `True`.

## Bounded Interpretation

- The frozen operator hierarchy supports a stable descriptive low-mode scaling regime across the declared refinement levels.
- Kernel size remains fixed, the sampled spectral envelope remains monotone, and the normalized low-mode proxy stays tightly collapsed after scaling.
- This is a branch-local feasibility result only. It is not a continuum, geometric, or physical-law claim.

## Non-Claims

- No continuum limit is established.
- No geometric invariants or universal coefficients are extracted.
- No observable-spectral coupling conclusion is made from this primary track.

Phase VII establishes spectral feasibility for the frozen operator hierarchy.
