# Geometry Emergence

This module explores when interaction systems acquire geometry-like structure.

It is an additive sandbox layer for HAOS-IIP. It does not modify the existing operator architecture and does not import continuum assumptions, field equations, or spacetime interpretations.

Canonical figures in this module follow a frozen quiet-numerics standard: single-column layout, serif text, major-tick grids only, monochrome structure with one perceptual gradient, and deterministic SVG/PDF/PNG export.

## Scope

The module measures only structural diagnostics:

- path distortion between interaction paths and transport arrival structure
- neighborhood persistence across transport steps
- emergent dimensionality from transport signatures
- flow bending as divergence deformation under transport
- recoverability coupling between stable transport and low-distortion structure

## Layout

```text
geometry_emergence/
    operators/
    metrics/
    experiments/
    diagnostics/
    notebooks/
    README.md
```

## Run

From the repository root:

```bash
python3 geometry_emergence/experiments/run_geometry_probe.py
```

With a config override:

```bash
python3 geometry_emergence/experiments/run_geometry_probe.py --config geometry_emergence/experiments/default_geometry_probe_config.json
```

Canonical transition panel:

```bash
python3 geometry_emergence/experiments/run_canonical_transition.py
```

Clustered micro-scan:

```bash
python3 geometry_emergence/experiments/run_clustered_micro_scan.py
```

Cluster-scale sweep:

```bash
python3 geometry_emergence/experiments/run_cluster_scale_sweep.py --config geometry_emergence/configs/cluster_scale_sweep_config.json
```

Canonical figure render:

```bash
python3 geometry_emergence/experiments/render_cluster_scale_canonical_figure.py
```

## Output

The experiment stores deterministic JSON output at:

```text
geometry_emergence/diagnostics/geometry_probe_results.json
```

The canonical transition panel stores:

```text
geometry_emergence/diagnostics/canonical_transition_summary.json
geometry_emergence/diagnostics/canonical_transition_primary_metrics.svg
```

The clustered micro-scan stores:

```text
geometry_emergence/diagnostics/clustered_micro_scan_summary.json
```

The cluster-scale sweep stores:

```text
geometry_emergence/diagnostics/cluster_scale_sweep_full.json
geometry_emergence/diagnostics/cluster_scale_sweep_summary.json
geometry_emergence/diagnostics/cluster_scale_transition.svg
geometry_emergence/diagnostics/cluster_scale_transition.pdf
geometry_emergence/diagnostics/cluster_scale_transition.png
```

The publication-grade canonical figure and bounded interpretation layer store:

```text
geometry_emergence/diagnostics/cluster_scale_transition_canonical.svg
geometry_emergence/diagnostics/cluster_scale_transition_canonical.pdf
geometry_emergence/diagnostics/cluster_scale_transition_canonical.png
geometry_emergence/diagnostics/cluster_scale_minimal_interpretation.md
papers/figures/cluster_scale_transition_canonical.pdf
papers/figures/cluster_scale_transition_canonical.png
papers/cluster_scale_dependence_of_coupled_transport_geometry_signatures_in_interaction_graphs.pdf
```

The default probe sweeps kernel width on a fixed clustered embedding and looks for the onset of a regime with:

- high recoverability
- low path distortion
- stable neighborhood relations

No interpretation beyond that is claimed.
