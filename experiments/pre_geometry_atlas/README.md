# Stage 10 Pre-Geometry Atlas

Atlas-0 is the baseline morphology pass for Stage 10. It reuses the frozen scalar and projected transverse operators and classifies packet behavior before interpretation.

Run the full Atlas-0 sheet:

`python numerics/simulations/stage10_atlas0.py`

Run a small sanity subset:

`python numerics/simulations/stage10_atlas0.py --max-runs 4`

Outputs:

- stamped JSON summary in `data/`
- stamped CSV run table in `data/`
- stamped plots in `plots/`
- summary note in `experiments/pre_geometry_atlas/`

Atlas-0 keeps regime labels descriptive and operator-level only. No physical interpretation is attached at this stage.
