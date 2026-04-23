# HAOS-IIP Localized Bump Response 53.2

## Purpose

`53.2` resolves the broad localized-bump `OPEN` row from `52.5` into a sharper threshold statement.

The original localized bump test produced maximum refinement drift `0.5124`, so it could not support a proto-particle / localized-excitation proxy.

The follow-up keeps the same scalar carrier and same bump shape, but changes the readout window:

> exclude the earliest source-core transient layer and test whether weak localized bump excitations recover a stable shell-native transport response.

This is still an external scalar-carrier runner. HAOS core is untouched.

## Construction

- bump shape: same localized scalar bump, `sigma = 0.18`
- refinements: `n = 13, 15`
- weak amplitudes: `eta = 0.01, 0.02, 0.03, 0.04, 0.05`
- stress amplitudes: `eta = 0.075, 0.10, 0.15`
- delayed fit window: `[0.010, 0.026]`
- thresholds: same current-closure thresholds used for `52.5`

## Result

The weak localized bump window passes:

- all `10/10` weak cases pass;
- maximum refinement drift: `0.1095`;
- maximum median relative error: `0.0943`;
- maximum p90 relative error: `0.2984`;
- maximum shell-kappa CV: `0.1248`;
- minimum metric-tracking correlation remains above threshold.

The weakest weak case is:

- `weak_bump|n=13|eta=0.050`
- `kappa / D_eff = 0.9699`
- median relative error `0.0943`
- p90 relative error `0.2984`
- shell-kappa CV `0.1248`
- `|metric-tracking corr| = 0.9752`

The stronger stress window remains open:

- maximum refinement drift: `0.3578`;
- maximum median relative error: `0.2784`;
- maximum p90 relative error: `0.5287`;
- maximum shell-kappa CV: `0.3406`.

The hardest stress case is:

- `stress_bump|n=13|eta=0.150`
- median relative error `0.2784`
- p90 relative error `0.5287`
- shell-kappa CV `0.3406`

## Correct Interpretation

The old statement was:

> localized bump response is open.

The new statement is:

> weak localized bump excitations close as a bounded metric-tracking transport response after the source-core transient is excluded, while stronger localized bumps remain outside the same closure.

This is a thresholded proto-particle-style result, not a particle model.

## What This Supports

`53.2` supports:

- a weak localized-excitation proxy on the scalar carrier;
- bounded metric-transport tracking for small localized bump amplitudes;
- a narrowed amplitude boundary between `eta = 0.05` and `eta = 0.075`;
- an improved bridge row from undifferentiated `OPEN` to weak `PASS` plus stress `OPEN`.

## What This Does Not Support

This note does not claim:

- particle ontology;
- arbitrary localized excitation closure;
- stable strong-source behavior;
- quantization;
- stress-energy coupling;
- curvature extraction.

## Authority and Artifacts

- script: `numerics/simulations/scalar_kernel_graph_localized_bump_response.py`
- operator note: `experiments/operators/Scalar_Kernel_Graph_Localized_Bump_Response_v1.md`
- result: `data/20260423_182539_scalar_kernel_graph_localized_bump_response.json`
- latest: `data/scalar_kernel_graph_localized_bump_response_latest.json`
- plots:
  - `plots/20260423_182539_scalar_kernel_graph_localized_bump_response_summary.png`
  - `plots/20260423_182539_scalar_kernel_graph_localized_bump_response_profiles.png`

## Next Honest Move

The next localized-excitation step should map the boundary between `eta = 0.05` and `eta = 0.075`, then test whether that threshold survives mild carrier disorder.
