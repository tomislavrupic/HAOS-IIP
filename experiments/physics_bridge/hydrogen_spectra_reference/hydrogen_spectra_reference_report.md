# Hydrogen Spectra Computational Reference Probe

Implemented fact: this sidecar runs a deterministic gross Rydberg-series reference and component controls.
Design choice: gross line positions use the declared relation `sigma = R_H * (1/n_lower^2 - 1/n_upper^2)`.
Heuristic: the reference score is the minimum of line match, fitted scale, inverse-square residual, line order, and selection-rule scores.
Unverified hypothesis: no CST or HAOS hydrogen-spectra recovery hypothesis is tested here.

## Verdict Labels
- HYDROGEN_REFERENCE_SANITY_PASS
- CST_HYDROGEN_STATUS_OPEN
- HAOS_DERIVATION_NOT_TESTED
- FINE_STRUCTURE_OPEN_NOT_MODELED
- NO_PHYSICAL_EXPERIMENT

## Reference
- declared R_H: `10967758.340280000120 m^-1`
- transition count: `15`
- fitted R_H: `10967758.340280001983 m^-1`
- max transition relative error: `0`
- reference score: `1.000000000000`

## Controls
- arithmetic_spacing_control: score `0.000000000000`, line `0.000000000000`, scale `0.000000000000`, law `0.000000000000`, order `1.000000000000`, selection `1.000000000000`, status `CONTROL_REJECTED`
- wrong_exponent_control: score `0.000000000000`, line `0.000000000000`, scale `0.000000000000`, law `0.000000000000`, order `1.000000000000`, selection `1.000000000000`, status `CONTROL_REJECTED`
- scaled_rydberg_control: score `0.000000000000`, line `0.000000000000`, scale `0.000000000000`, law `1.000000000000`, order `1.000000000000`, selection `1.000000000000`, status `CONTROL_REJECTED`
- reversed_transition_identity_control: score `0.000000000000`, line `0.000000000000`, scale `0.232764583742`, law `0.000000000000`, order `0.000000000000`, selection `1.000000000000`, status `CONTROL_REJECTED`
- forbidden_selection_rule_control: score `0.000000000000`, line `1.000000000000`, scale `1.000000000000`, law `1.000000000000`, order `1.000000000000`, selection `0.000000000000`, status `CONTROL_REJECTED`

## Boundary
- This is not a laboratory hydrogen spectrum.
- This does not model fine structure, Lamb shift, hyperfine splitting, line intensities, or apparatus effects.
- This does not derive hydrogen spectra from CST or HAOS-IIP.
- This does not change the frozen CST v0.2.2 checkpoint.
- `HYDROGEN_REFERENCE_SANITY_PASS` means only that the computational reference harness behaves as expected.
