# Kernel Laplacian Convergence v1

## Purpose

Test whether the interaction kernel induces a graph operator that converges to the Laplace operator in the continuum limit.

The question is strictly operator-level:

> does the kernel-weighted graph Laplacian reproduce `Delta f` for smooth test functions as the cubic lattice is refined?

No new physics assumptions are added here.

## Construction

### Substrates

Regular cubic 3D grids with

- `n = 9, 11, 13, 17, 21`
- lattice spacing `h = 1/(n-1)`

### Kernel weights

For grid points `x_i`,

$$
w_{ij} = \exp\left(-\frac{|x_i-x_j|^2}{\epsilon_k}\right),
$$

with

$$
\epsilon_k = c_\epsilon h^2,
\qquad
c_\epsilon \in \{0.5, 1.0, 2.0\}.
$$

The graph is cut off at radius

$$
r_c = 2.5 \sqrt{\epsilon_k}.
$$

### Graph Laplacian

$$
L = D - W,
\qquad
D_{ii} = \sum_j w_{ij}.
$$

### Induced operator

Because the scan uses a local stencil with `epsilon_k ~ h^2`, the correct continuum normalization is set by the discrete second moment of the actual kernel stencil, not by the continuum Gaussian moment.

Define

$$
\widehat{L}_h = -\frac{2}{\mu_2 h^2} L,
$$

where `mu_2` is the averaged discrete second moment of the stencil.

This is the operator tested against the analytic Laplacian.

## Test Functions

Three analytic functions were used:

$$
f_1 = x^2 + y^2 + z^2,
$$

$$
f_2 = \sin(\pi x)\sin(\pi y)\sin(\pi z),
$$

$$
f_3 = \exp(-(x^2+y^2+z^2)).
$$

Analytic Laplacians:

$$
\Delta f_1 = 6,
$$

$$
\Delta f_2 = -3\pi^2 f_2,
$$

$$
\Delta f_3 = (4r^2 - 6) f_3.
$$

Errors were measured on interior nodes only, excluding a boundary layer equal to the maximum stencil offset so that every retained node has a full symmetric neighborhood.

The discrete error norm is

$$
\|\widehat{L}_h f - \Delta f\|_{L^2}
=
\left(h^3 \sum_{\text{interior}} |\widehat{L}_h f - \Delta f|^2\right)^{1/2}.
$$

## Artifacts

- script: `numerics/simulations/kernel_operator_convergence.py`
- results: `data/20260310_111522_kernel_operator_convergence.json`
- latest: `data/kernel_operator_convergence_latest.json`
- plots:
  - `plots/20260310_111522_operator_error_vs_h.png`
  - `plots/20260310_111522_operator_profiles.png`

## Results

### 1. Quadratic test

The quadratic function is reproduced to machine precision across the entire scan.

Representative errors:

- `n=9`, `c_epsilon=0.5`: `1.97e-15`
- `n=21`, `c_epsilon=1.0`: `6.54e-14`
- `n=13`, `c_epsilon=2.0`: `9.69e-15`

So the normalized kernel operator reproduces the constant Laplacian of `f_1` essentially exactly on interior nodes.

### 2. Trigonometric test

Fitted convergence orders:

- `c_epsilon=0.5`: `p = 1.9825`
- `c_epsilon=1.0`: `p = 1.8579`
- `c_epsilon=2.0`: `p = 1.3770`

This is clear positive-order convergence, close to second order for the more local kernels.

### 3. Gaussian test

Fitted convergence orders:

- `c_epsilon=0.5`: `p = 1.8035`
- `c_epsilon=1.0`: `p = 1.2495`
- `c_epsilon=2.0`: `p = 0.3316`

Again the local kernels converge clearly. The broader kernel still improves with refinement, but much more slowly at the present resolutions.

## Interpretation

The operator-level verdict is:

> after discrete second-moment normalization, the kernel-induced graph operator converges toward the continuum Laplacian on the cubic scan.

The evidence is strongest for `c_epsilon <= 1.0`:

- exact quadratic reproduction
- strong convergence on the oscillatory and Gaussian tests
- profile agreement on the finest grid

The `c_epsilon = 2.0` branch is still convergent for the trigonometric test, but noticeably less accurate on the Gaussian profile at finite resolution. So broad kernels remain less reliable if the grid is not yet fine.

## Current Verdict

- yes, the interaction kernel induces a Laplacian-type operator in the continuum limit
- the conclusion is strongest in the local-kernel regime `epsilon_k = c_epsilon h^2` with `c_epsilon <= 1`
- broader kernels are still finite-resolution sensitive and should not be used as the main convergence evidence
