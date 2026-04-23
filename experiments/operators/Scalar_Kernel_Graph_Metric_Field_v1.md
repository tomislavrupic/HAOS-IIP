# Scalar Kernel Graph Metric Field v1

## Purpose

Upgrade the scalar-carrier geometry read from one global coordinate-stiffness tensor to a **local metric-like tensor field** on the same operator/substrate family.

The bounded question is:

> does the same scalar kernel-graph carrier support a positive, near-isotropic, spatially stable local metric-like field in the bulk, under the same mild disorder and bounded kernel-family window already used in the `51.3` robustness pass?

## Construction

The carrier stays fixed:

- scalar kernel graph on `3D` cubic point clouds
- local regime `c_epsilon = 0.5`
- same mild point-set disorder window
- same bounded kernel-family window

The local metric-field estimator is taken directly from the operator's row-local quadratic response. For each node `i`,

$$
G_{ab}(i) = \frac{1}{2} \sum_{j \neq i} \widehat{L}_{ij} (x_j^a - x_i^a)(x_j^b - x_i^b).
$$

This is the right bounded local estimator on the present carrier because it measures the operator's second-order response directly at each node, rather than importing a new geometric fitting rule.

The raw node-level field is kept as provenance, but the bounded physical read is a local coarse-grained field averaged over a bulk neighborhood radius of `2.5 h`. This is the same kind of move the repository already used in other places when raw local readouts were too stencil-sensitive to carry the physical statement alone.

Bulk readout excludes a boundary collar of `2` lattice layers.

## Artifacts

- script: `numerics/simulations/scalar_kernel_graph_metric_field.py`
- config: `config.json -> scalar_kernel_graph_metric_field`
- results: `data/20260420_222424_scalar_kernel_graph_metric_field.json`
- latest: `data/scalar_kernel_graph_metric_field_latest.json`
- plots:
  - `plots/20260420_222424_scalar_kernel_graph_metric_field_summary.png`
  - `plots/20260420_222424_scalar_kernel_graph_metric_field_tensors.png`

## Disorder pass

| case | bulk nodes | min eigenvalue | mean anisotropy | mean offdiag ratio | spatial trace CV | normalized mean-tensor error | status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `disorder|n=11|jitter=0.000` | `343` | `1.1439` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | PASS |
| `disorder|n=11|jitter=0.020` | `214` | `1.1235` | `0.0067` | `0.0036` | `0.0027` | `0.0008` | PASS |
| `disorder|n=11|jitter=0.040` | `214` | `1.1051` | `0.0139` | `0.0073` | `0.0052` | `0.0018` | PASS |
| `disorder|n=13|jitter=0.000` | `729` | `1.1197` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | PASS |
| `disorder|n=13|jitter=0.020` | `523` | `1.0953` | `0.0063` | `0.0031` | `0.0030` | `0.0008` | PASS |
| `disorder|n=13|jitter=0.040` | `523` | `1.0768` | `0.0129` | `0.0064` | `0.0059` | `0.0016` | PASS |

## Kernel-family pass

| case | bulk nodes | min eigenvalue | mean anisotropy | mean offdiag ratio | spatial trace CV | normalized mean-tensor error | status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `kernel|n=11|family=gaussian_local` | `343` | `1.1439` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | PASS |
| `kernel|n=11|family=gaussian_half` | `343` | `1.1899` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | PASS |
| `kernel|n=11|family=inverse_quadratic` | `343` | `1.1975` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | PASS |
| `kernel|n=13|family=gaussian_local` | `729` | `1.1197` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | PASS |
| `kernel|n=13|family=gaussian_half` | `729` | `1.1576` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | PASS |
| `kernel|n=13|family=inverse_quadratic` | `729` | `1.1638` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | PASS |

## Refinement drift of the normalized mean tensor

| pair | drift |
| --- | --- |
| `disorder|jitter=0.000` | `0.0000` |
| `disorder|jitter=0.020` | `0.0006` |
| `disorder|jitter=0.040` | `0.0023` |
| `kernel|family=gaussian_local` | `0.0000` |
| `kernel|family=gaussian_half` | `0.0000` |
| `kernel|family=inverse_quadratic` | `0.0000` |

## Interpretation

- observation: the same scalar carrier now supports a stable local metric-like tensor field: the operator-level row-local quadratic response stays positive, near-isotropic, and nearly constant across the tested bulk window under refinement, mild disorder, and bounded kernel-family variation
- conclusion: all 6 disorder cases and all 6 kernel-family cases pass the local metric-field thresholds, and the maximum normalized refinement drift is 0.0023; the weakest passing case is `disorder|n=11|jitter=0.040` with normalized mean-tensor error 0.0018, mean anisotropy 0.0139, and spatial trace CV 0.0052; this means the scalar-carrier geometry bridge now extends from one global coordinate-stiffness read to a bounded stable local metric-field statement on the tested window

## Current boundary

If this note is read positively, the correct bounded statement is:

> the scalar-carrier geometry bridge now extends to a stable local bulk metric-like tensor field on the tested window.

This note still does **not** claim:

- curvature extraction
- current conservation
- broad universality beyond the tested mild carrier window
- ontology or spacetime
