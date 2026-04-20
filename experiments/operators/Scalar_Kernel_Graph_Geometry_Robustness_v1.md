# Scalar Kernel Graph Geometry Robustness v1

## Purpose

Run the next honest move after the bounded scalar-carrier `CP4` bridge:

1. mild point-set disorder on the same carrier
2. bounded kernel-family variation on the same carrier
3. response / metric extraction on that same coarse reconstruction

This note does **not** introduce a new geometry story. It tests whether the existing scalar-carrier geometry bridge survives these bounded stresses.

## Construction

The carrier remains the same:

- scalar kernel graph on `3D` cubic point clouds
- local regime `c_epsilon = 0.5`
- same Green / heat / shell-arrival / first-shell low-mode observables
- same Euclidean reconstruction centered at the source

Two bounded extensions are tested.

### 1. Disorder pass

- `n = [11, 13]`
- jitter fractions `= [0.0, 0.02, 0.04]`
- kernel family fixed to `gaussian_local`

For jittered cases, shell-arrival bins are assigned to the nearest clean-grid shell at the same `n`, so the coarse reconstruction itself stays fixed while the point cloud is perturbed.

### 2. Kernel-family pass

- `n = [11, 13]`
- families `= ['gaussian_local', 'gaussian_half', 'inverse_quadratic']`

Families tested:

- `gaussian_local`: `exp(-|x_i-x_j|^2 / epsilon_k)`
- `gaussian_half`: `exp(-|x_i-x_j|^2 / (2 epsilon_k))`
- `inverse_quadratic`: `(1 + |x_i-x_j|^2 / epsilon_k)^(-2)`

### 3. Response / metric extraction

On every case, extract the coordinate-response stiffness tensor on the orthonormalized coordinate basis `span(x,y,z)`:

$$
K_{coord} = Q^T (-\widehat{L}_h) Q.
$$

This is the bounded metric-response read for the same carrier. In the isotropic case it should stay close to a scalar multiple of the identity.

## Artifacts

- script: `numerics/simulations/scalar_kernel_graph_geometry_robustness.py`
- config: `config.json -> scalar_kernel_graph_geometry_robustness`
- results: `data/20260420_131037_scalar_kernel_graph_geometry_robustness.json`
- latest: `data/scalar_kernel_graph_geometry_robustness_latest.json`
- plots:
  - `plots/20260420_131037_scalar_kernel_graph_geometry_robustness_disorder.png`
  - `plots/20260420_131037_scalar_kernel_graph_geometry_robustness_kernels.png`
  - `plots/20260420_131037_scalar_kernel_graph_geometry_robustness_metric_tensor.png`

## Disorder pass

| case | Green fit `R^2` | Green dim hint | Heat fit `R^2` | shell-arrival fit `R^2` | min low-mode cosine | metric anisotropy | metric offdiag ratio | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `disorder|n=11|jitter=0.000` | `0.9834` | `2.9673` | `0.9993` | `0.9797` | `0.9934` | `0.0000` | `0.0000` | PASS |
| `disorder|n=11|jitter=0.020` | `0.9807` | `2.9682` | `0.9993` | `0.9840` | `0.9933` | `0.0003` | `0.0005` | PASS |
| `disorder|n=11|jitter=0.040` | `0.9736` | `2.9684` | `0.9993` | `0.9834` | `0.9933` | `0.0005` | `0.0014` | PASS |
| `disorder|n=13|jitter=0.000` | `0.9881` | `2.9663` | `0.9993` | `0.9811` | `0.9932` | `0.0000` | `0.0000` | PASS |
| `disorder|n=13|jitter=0.020` | `0.9876` | `2.9666` | `0.9993` | `0.9818` | `0.9932` | `0.0004` | `0.0003` | PASS |
| `disorder|n=13|jitter=0.040` | `0.9856` | `2.9664` | `0.9993` | `0.9834` | `0.9930` | `0.0007` | `0.0006` | PASS |

## Kernel-family pass

| case | Green fit `R^2` | Green dim hint | Heat fit `R^2` | shell-arrival fit `R^2` | min low-mode cosine | metric anisotropy | metric offdiag ratio | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `kernel|n=11|family=gaussian_local` | `0.9834` | `2.9673` | `0.9993` | `0.9797` | `0.9934` | `0.0000` | `0.0000` | PASS |
| `kernel|n=11|family=gaussian_half` | `0.9856` | `2.9703` | `0.9990` | `0.9734` | `0.9932` | `0.0000` | `0.0000` | PASS |
| `kernel|n=11|family=inverse_quadratic` | `0.9857` | `2.9708` | `0.9990` | `0.9734` | `0.9932` | `0.0000` | `0.0000` | PASS |
| `kernel|n=13|family=gaussian_local` | `0.9881` | `2.9663` | `0.9993` | `0.9811` | `0.9932` | `0.0000` | `0.0000` | PASS |
| `kernel|n=13|family=gaussian_half` | `0.9892` | `2.9687` | `0.9991` | `0.9732` | `0.9931` | `0.0000` | `0.0000` | PASS |
| `kernel|n=13|family=inverse_quadratic` | `0.9893` | `2.9691` | `0.9991` | `0.9732` | `0.9931` | `0.0000` | `0.0000` | PASS |

## Interpretation

- observation: the scalar kernel-graph geometry bridge survives the first bounded robustness pass: mild point-set disorder, bounded kernel-family variation, and coordinate-response extraction all stay compatible with the same 3D-like metric read on the tested carrier
- conclusion: all 6 disorder cases and all 6 kernel-family cases pass the shared CP4 thresholds; the weakest disorder case is `disorder|n=11|jitter=0.000` with shell-arrival fit R^2 = 0.9797, and the weakest kernel case is `kernel|n=13|family=gaussian_half` with shell-arrival fit R^2 = 0.9732; this means the scalar-carrier geometry closure is now robust across the tested mild disorder and bounded kernel-family window, while remaining a bounded carrier-level statement rather than full universality

## Current boundary

If this note is read positively, the correct bounded statement is:

> the scalar-carrier geometry closure now survives the first mild disorder and bounded kernel-family robustness pass, while the same coordinate-response tensor remains close to isotropic on the tested window.

The remaining open step is still stronger than this:

- broader disorder families
- broader kernel families
- stronger response / metric extraction beyond the coordinate stiffness tensor
