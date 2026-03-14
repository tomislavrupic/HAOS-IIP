# Atlas-0 Baseline Morphology v1

Atlas-0 is the baseline morphology pass for Stage 10.

Purpose:
- fix the experimental grammar for pre-geometry classification
- keep observables consistent across runs
- classify packet behavior before interpretation

Current scaffold:
- explicit run sheet in `numerics/simulations/stage10_atlas0_runs.json`
- shared metrics and export helpers in `numerics/simulations/stage10_common.py`
- regime labeling in `numerics/simulations/stage10_regime_labels.py`
- Atlas-0 driver in `numerics/simulations/stage10_atlas0.py`

Interpretation boundary:
Atlas-0 is a morphology and constraint-classification framework only. It does not assert continuum physics or particle content.
