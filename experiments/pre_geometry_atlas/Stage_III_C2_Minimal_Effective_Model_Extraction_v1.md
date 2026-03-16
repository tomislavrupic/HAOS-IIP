# Stage III-C2 Minimal Effective Model Extraction v1

This reduction stays inside the frozen HAOS (Harmonic Address Operating System) / IIP (Interaction Invariance Physics) Dirac-Kahler clustered braid window.

Timestamped JSON: `data/20260316_104146_stage23_11_dk_minimal_effective_model_extraction.json`
Timestamped CSV: `data/20260316_104146_stage23_11_dk_minimal_effective_model_extraction.csv`

effective_model_valid = `TRUE`

## Calibration note
The source brief listed more phase and detuning values than a 9-run matrix can contain. This implementation therefore uses the resolved 3 x 3 anchor lattice `phase in {0.375, 0.500, 0.625}` crossed with bounded clustered-seed grade detuning `delta in {0.0, 0.5, 0.75}`.
The extracted `X` coordinate is the signed separation envelope of the ordered peak tracks, because the unsigned clustered peak distance remains nearly constant on this branch.

## Fitted map
Use the reduced state
- `X_n`: signed separation envelope
- `Y_n`: net grade-exchange amplitude

`X_{n+1} = X_n + c0 + c1 Y_n + c2 Y_n^2 + c3 phi + c4 phi Y_n`
`Y_{n+1} = Y_n + d0 + d1 X_n + d2 X_n^2 + d3 delta + d4 delta X_n`

### X-map coefficients
- `c0` = `-0.005783`
- `c1` = `0.051961`
- `c2` = `-0.818464`
- `c3` = `-0.023173`
- `c4` = `0.504837`

### Y-map coefficients
- `d0` = `0.012768`
- `d1` = `-0.017323`
- `d2` = `1.351037`
- `d3` = `0.003961`
- `d4` = `0.074234`

## Learned classifier thresholds
- midpoint phase halfwidth: `0.0625`
- unresolved detuning threshold: `0.2500`
- smear peak threshold: `0.1783`

## Validation metrics
- topology agreement fraction: `1.0000`
- detuning threshold error: `0.0000`
- phase-band width error: `0.0000`
- single smeared band reproduced: `True`
- mean survival-time envelope mismatch: `0.0139`

## Per-run summary
- `S23_11_p0375_d000`: obs=`braid_like_exchange`, pred=`braid_like_exchange`, `x_rmse=0.0944`, `y_rmse=0.0404`
- `S23_11_p0500_d000`: obs=`transfer_smeared`, pred=`transfer_smeared`, `x_rmse=0.1088`, `y_rmse=0.0224`
- `S23_11_p0625_d000`: obs=`braid_like_exchange`, pred=`braid_like_exchange`, `x_rmse=0.0610`, `y_rmse=0.0279`
- `S23_11_p0375_d050`: obs=`unresolved`, pred=`unresolved`, `x_rmse=0.1013`, `y_rmse=0.0201`
- `S23_11_p0500_d050`: obs=`unresolved`, pred=`unresolved`, `x_rmse=0.0415`, `y_rmse=0.0167`
- `S23_11_p0625_d050`: obs=`unresolved`, pred=`unresolved`, `x_rmse=0.0303`, `y_rmse=0.0128`
- `S23_11_p0375_d075`: obs=`transfer_smeared`, pred=`transfer_smeared`, `x_rmse=0.0403`, `y_rmse=0.0320`
- `S23_11_p0500_d075`: obs=`transfer_smeared`, pred=`transfer_smeared`, `x_rmse=0.0401`, `y_rmse=0.0357`
- `S23_11_p0625_d075`: obs=`transfer_smeared`, pred=`transfer_smeared`, `x_rmse=0.0319`, `y_rmse=0.0417`

## Interpretation
The clustered DK braid / smear sector admits a compact two-state surrogate on this anchor lattice. The reduction is descriptive rather than ontic: it compresses the clustered texture without overturning the Stage 23.10 result that the braid does not generalize into a family-wide intrinsic mechanism.

## Plots
- `plots/20260316_104146_S23_11_p0375_d000_phase_plane.png`
- `plots/20260316_104146_S23_11_p0500_d000_phase_plane.png`
- `plots/20260316_104146_S23_11_p0625_d000_phase_plane.png`
- `plots/20260316_104146_S23_11_p0375_d050_phase_plane.png`
- `plots/20260316_104146_S23_11_p0500_d050_phase_plane.png`
- `plots/20260316_104146_S23_11_p0625_d050_phase_plane.png`
- `plots/20260316_104146_S23_11_p0375_d075_phase_plane.png`
- `plots/20260316_104146_S23_11_p0500_d075_phase_plane.png`
- `plots/20260316_104146_S23_11_p0625_d075_phase_plane.png`
- `plots/20260316_104146_stage23_11_effective_phase_diagram.png`
- `plots/20260316_104146_stage23_11_failure_map.png`
