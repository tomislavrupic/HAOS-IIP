# Semiconductor Band-Structure Computational Reference Probe

Implemented fact: this sidecar runs a deterministic toy direct-gap two-band reference and component controls.
Design choice: the reference uses declared cosine bands on a reduced unitless k-grid.
Heuristic: the reference score is the minimum of gap, directness, curvature, band-shape, band-order, and symmetry scores.
Unverified hypothesis: no CST or HAOS semiconductor band-recovery hypothesis is tested here.

## Verdict Labels
- SEMICONDUCTOR_REFERENCE_SANITY_PASS
- CST_SEMICONDUCTOR_STATUS_OPEN
- HAOS_DERIVATION_NOT_TESTED
- TOY_BAND_STRUCTURE_ONLY
- NO_PHYSICAL_EXPERIMENT

## Reference
- reference gap: `1.200000000000 eV`
- k-point count: `81`
- measured band gap: `1.200000000000 eV`
- conduction edge k: `0.000000000000`
- valence edge k: `0.000000000000`
- reference score: `1.000000000000`

## Controls
- scaled_gap_control: score `0.000000000000`, gap `0.000000000000`, direct `1.000000000000`, curvature `1.000000000000`, shape `0.000000000000`, order `1.000000000000`, symmetry `1.000000000000`, status `CONTROL_REJECTED`
- metallic_overlap_control: score `0.000000000000`, gap `0.000000000000`, direct `1.000000000000`, curvature `1.000000000000`, shape `0.000000000000`, order `0.814814814815`, symmetry `1.000000000000`, status `CONTROL_REJECTED`
- indirect_gap_control: score `0.000000000000`, gap `1.000000000000`, direct `0.000000000000`, curvature `1.000000000000`, shape `0.000000000000`, order `1.000000000000`, symmetry `0.000000000000`, status `CONTROL_REJECTED`
- wrong_curvature_control: score `0.000000000000`, gap `1.000000000000`, direct `1.000000000000`, curvature `0.000000000000`, shape `0.000000000000`, order `1.000000000000`, symmetry `1.000000000000`, status `CONTROL_REJECTED`
- flat_band_control: score `0.000000000000`, gap `1.000000000000`, direct `0.000000000000`, curvature `0.000000000000`, shape `0.000000000000`, order `1.000000000000`, symmetry `1.000000000000`, status `CONTROL_REJECTED`
- band_label_swap_control: score `0.000000000000`, gap `0.000000000000`, direct `0.000000000000`, curvature `0.000000000000`, shape `0.000000000000`, order `0.000000000000`, symmetry `1.000000000000`, status `CONTROL_REJECTED`
- k_order_shuffle_control: score `0.000000000000`, gap `1.000000000000`, direct `0.000000000000`, curvature `0.000000000000`, shape `0.000000000000`, order `1.000000000000`, symmetry `0.308641975309`, status `CONTROL_REJECTED`

## Boundary
- This is not a physical semiconductor calculation.
- This does not model real materials, doping, defects, many-body effects, transport, devices, or carrier dynamics.
- This does not derive band structure from CST or HAOS-IIP.
- This does not change the frozen CST v0.2.2 checkpoint.
- `SEMICONDUCTOR_REFERENCE_SANITY_PASS` means only that the toy reference harness behaves as expected.
