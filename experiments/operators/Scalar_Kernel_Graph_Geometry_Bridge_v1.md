# Scalar Kernel Graph Geometry Bridge v1

## Purpose

Close the preflight blocker for geometry closure by rebuilding all four geometry-facing observables on **one common scalar kernel-graph family** instead of mixing the old 3D Green-response line with the frozen 2D branch-local torus line.

The question here is deliberately bounded:

> on the same 3D cubic scalar kernel graph, under the same Euclidean shell reconstruction, do Green response, heat behavior, shell-arrival structure, and low-mode organization point to one coherent metric-like geometry read?

## Shared carrier

Use the strongest current scalar operator family:

- cubic 3D grids with `n = 9, 11, 13`
- local kernel regime `c_epsilon = 0.5`
- cutoff factor `2.5`
- kernel weights

$$
w_{ij} = \exp\left(-\frac{|x_i-x_j|^2}{\epsilon_k}\right),
\qquad
\epsilon_k = c_\epsilon h^2
$$

- induced scalar operator

$$
\widehat{L}_h = -\frac{2}{\mu_2 h^2} (D-W)
$$

This is the same local scalar carrier already used for the strongest CP1 operator-control result.

## One common metric-like coarse reconstruction

All four channels use the same Euclidean coarse geometry centered at the source node:

- shell radius `r`
- shell second moment `r^2`
- coordinate subspace `span(x, y, z)`

No branch-local torus artifacts enter this bridge.

## Measured channels

1. **Green response**
   Fit the mean-zero Poisson field from `-\widehat{L}_h \phi = \delta - N^{-1} 1` against `A + B/r`.

2. **Heat behavior**
   Evolve `u_t = \widehat{L}_h u` from a point source and test whether the shell second moment `\langle r^2 \rangle(t)` is linear in `t`.

3. **Shell-arrival structure**
   On the same heat trajectory, define shell arrival as the first half-peak crossing of the shell-mean signal and fit arrival time against `r^2`.

4. **Low-mode organization**
   Compare the first positive eigen-triplet subspace of `-\widehat{L}_h` against the Euclidean coordinate subspace `span(x, y, z)` using principal angles. This avoids fake instability from arbitrary rotations inside the degenerate first shell.

## Artifacts

- script: `numerics/simulations/scalar_kernel_graph_geometry_bridge.py`
- config: `config.json -> scalar_kernel_graph_geometry_bridge`
- results: `data/20260419_234444_scalar_kernel_graph_geometry_bridge.json`
- latest: `data/scalar_kernel_graph_geometry_bridge_latest.json`
- plots:
  - `plots/20260419_234444_scalar_kernel_graph_geometry_bridge_summary.png`
  - `plots/20260419_234444_scalar_kernel_graph_geometry_bridge_detail.png`

## Results

| `n` | Green dim hint | Green fit `R^2` | `D_eff` | Heat fit `R^2` | arrival slope | arrival fit `R^2` | min principal cosine | subspace affinity |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `9` | `2.9758` | `0.9769` | `0.9551` | `0.9996` | `0.0775` | `0.9886` | `0.9937` | `0.9875` |
| `11` | `2.9673` | `0.9834` | `0.9594` | `0.9996` | `0.0725` | `0.9779` | `0.9934` | `0.9869` |
| `13` | `2.9663` | `0.9881` | `0.9616` | `0.9996` | `0.0762` | `0.9787` | `0.9932` | `0.9865` |

## Interpretation

- observation: the shared-family scalar kernel-graph carrier now gives one coherent 3D-like geometry read across all four channels: inverse-distance Green response, linear diffusion on shell second moments, shell-arrival times linear in radius squared, and a first low-mode triplet that matches the Euclidean coordinate subspace
- conclusion: the previous geometry-closure preflight mismatch is removed on the tested cubic scalar carrier: one bounded metric-like reconstruction now organizes Green response, heat behavior, shell-arrival structure, and low-mode organization on the same operator/substrate family. This is a carrier-level bridge, not yet a universality or full geometry-closure claim beyond the tested cubic local-kernel regime.

## Current verdict

This bridge fixes the **category mistake** identified in the geometry-closure preflight note: all four channels now live on the same scalar operator/substrate family.

The positive read should stay bounded:

- this supports a common 3D-like metric organization on the tested local cubic scalar carrier;
- it does **not** yet claim universality across other kernel families, point-set disorder families, or full geometry emergence beyond the tested bridge regime.
