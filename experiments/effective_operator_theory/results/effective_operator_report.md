# Effective Operator Expansion Report

## Status

- Result: `PASS`
- Classification: `SYNTHETIC_EFFECTIVE_OPERATOR_SCAFFOLD`
- Result hash: `effective_operator_a740ce933dd1dbd479931847`
- leading Laplacian coefficient: `0.17987548725`
- higher-order correction coefficient: `0.0136274925612`
- correction ratio: `0.0757606985231`
- fit R^2: `0.999999981611`
- max long-wavelength relative residual: `0.000673617871882`

## Labels

- `EFFECTIVE_OPERATOR_EXPANSION_BUILT`
- `EFT_DISCIPLINE_METHOD_IMPORTED`
- `DIFFUSION_LIKE_LEADING_TERM_DETECTED`
- `SUPPRESSED_CORRECTION_HIERARCHY_REPORTED`
- `CUTOFF_DECLARED`
- `CONTROL_RESPONSE_REPORTED`
- `EFT_QG_NOT_DERIVED`
- `PHYSICAL_FIELD_THEORY_NOT_ESTABLISHED`
- `PHYSICAL_BRIDGE_NOT_ESTABLISHED`
- `EFFECTIVE_OPERATOR_SCAFFOLD_PASS`

## Controls

| Control | Leading coeff | Correction coeff | R^2 | Status |
| --- | ---: | ---: | ---: | --- |
| identity_operator_control | -1.43539468358e-32 | 2.54604976057e-32 | 1 | `PASS` |
| unstable_sign_control | -0.180102930412 | -0.010426826667 | 0.999999987317 | `PASS` |

## Boundary

This benchmark imports EFT discipline: scale hierarchy, allowed terms, cutoff, matching, and controls.
It does not derive EFT quantum gravity, physical field theory, spacetime, Lorentz invariance, matter sectors, constants, or empirical physics.
