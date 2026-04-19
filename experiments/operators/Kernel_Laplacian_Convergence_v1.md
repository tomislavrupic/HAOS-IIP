# Kernel Laplacian Convergence v1

## Purpose

Test whether the scalar kernel-induced operator keeps a clean continuum limit only on the ideal cubic interior scan, or whether that control survives two stricter CP1 extensions:

- weakly perturbed 3D point sets
- explicit homogeneous boundary-condition families

The question remains strictly operator-level:

> does the same local kernel support a stable scalar Laplacian-type continuum control across refinement, mild geometric disorder, and boundary handling?

No new physics assumptions are added here.

## Construction

### 1. Cubic interior baseline

Regular cubic 3D grids with:

- `n = 9, 11, 13, 17, 21`
- lattice spacing `h = 1 / (n - 1)`
- kernel-width scan `c_epsilon in {0.5, 1.0, 2.0}`
- cutoff radius `r_c = 2.5 sqrt(epsilon_k)`

Kernel weights:

$$
w_{ij} = \exp\left(-\frac{|x_i-x_j|^2}{\epsilon_k}\right),
\qquad
\epsilon_k = c_\epsilon h^2
$$

Graph operator:

$$
L = D - W
$$

Interior scalar normalization:

$$
\widehat{L}_h = -\frac{2}{\mu_2 h^2} L
$$

where `mu_2` is the discrete second moment of the actual stencil.

### 2. Weakly perturbed point-set control

The same scalar test family is evaluated on weakly jittered copies of the cubic grid:

- `n = 11, 13, 15, 17`
- jitter amplitude `0.04 h`
- fixed seed `42`
- local kernel regime `c_epsilon = 0.5`

Boundary nodes stay fixed on the box and interior points are jittered.

Because the regular-grid stencil no longer applies directly, the operator is recovered from the same kernel weights by a weighted local quadratic fit around each point. The Laplacian is read off from the fitted quadratic coefficients.

This is still a kernel-induced scalar operator test. It is not a new phase contract and it does not change the frozen public API.

### 3. Boundary-condition controls

The same cubic scan is extended from interior-only error to full-domain boundary handling by reflected ghost values:

- odd reflection for homogeneous Dirichlet data
- even reflection for homogeneous Neumann data

Boundary-compatible test families:

$$
f_D = \sin(\pi x)\sin(\pi y)\sin(\pi z),
\qquad
\Delta f_D = -3\pi^2 f_D
$$

$$
f_N = \cos(\pi x)\cos(\pi y)\cos(\pi z),
\qquad
\Delta f_N = -3\pi^2 f_N
$$

## Test Functions

Scalar interior functions:

$$
f_1 = x^2 + y^2 + z^2,
\qquad
\Delta f_1 = 6
$$

$$
f_2 = \sin(\pi x)\sin(\pi y)\sin(\pi z),
\qquad
\Delta f_2 = -3\pi^2 f_2
$$

$$
f_3 = \exp(-(x^2+y^2+z^2)),
\qquad
\Delta f_3 = (4r^2 - 6) f_3
$$

For the cubic baseline and the jittered point-set control, errors are measured on the interior region supported by the local kernel.

For the boundary controls, errors are measured on the full domain.

In every case the discrete norm is

$$
\| \widehat{L}_h f - \Delta f \|_{L^2}
=
\left(h^3 \sum |\widehat{L}_h f - \Delta f|^2\right)^{1/2}.
$$

## Artifacts

- script: `numerics/simulations/kernel_operator_convergence.py`
- config: `config.json -> kernel_operator_convergence`
- results: `data/20260412_101326_kernel_operator_convergence.json`
- latest: `data/kernel_operator_convergence_latest.json`
- plots:
  - `plots/20260412_101326_operator_error_vs_h.png`
  - `plots/20260412_101326_operator_profiles.png`
  - `plots/20260412_101326_point_set_operator_error_vs_h.png`
  - `plots/20260412_101326_boundary_condition_error_vs_h.png`
  - `plots/20260412_101326_boundary_condition_profiles.png`

## Results

### 1. Cubic interior baseline

The original scalar result remains intact.

- quadratic reproduction stays at machine precision across the scan
- trigonometric convergence orders:
  - `c_epsilon = 0.5`: `p = 1.9825`
  - `c_epsilon = 1.0`: `p = 1.8579`
  - `c_epsilon = 2.0`: `p = 1.3770`
- Gaussian convergence orders:
  - `c_epsilon = 0.5`: `p = 1.8035`
  - `c_epsilon = 1.0`: `p = 1.2495`
  - `c_epsilon = 2.0`: `p = 0.3316`

Interpretation:

> after discrete second-moment normalization, the cubic interior scan remains a clean scalar continuum control, strongest in the local-kernel regime `c_epsilon <= 1`.

### 2. Weakly perturbed point-set control

In the bounded CP1 point-set window:

- `n = 11, 13, 15, 17`
- jitter amplitude `0.04 h`
- `c_epsilon = 0.5`

the point-set control now remains convergent.

Measured orders:

- quadratic reproduction stays near machine precision:
  - max `L^2` error across the scan: `5.96e-14`
- trigonometric convergence:
  - `c_epsilon = 0.5`: `p = 1.8223`
- Gaussian convergence:
  - `c_epsilon = 0.5`: `p = 1.0611`

Operator-quality diagnostics:

- valid stencil fraction stays above `0.975`
- average neighbor count stays near `20-21`

Interpretation:

> scalar continuum control survives mild geometric disorder in the strictly local kernel regime once the operator is recovered from the same kernel weights by local quadratic reproduction.

This is a bounded result. It does not say every wider kernel or every coarser jittered scan behaves equally well.

### 3. Boundary-condition controls

The reflected boundary controls pass cleanly on both tested families.

Dirichlet-compatible convergence:

- `c_epsilon = 0.5`: `p = 1.9825`
- `c_epsilon = 1.0`: `p = 1.9680`
- `c_epsilon = 2.0`: `p = 1.9405`

Neumann-compatible convergence:

- `c_epsilon = 0.5`: `p = 2.1902`
- `c_epsilon = 1.0`: `p = 2.1756`
- `c_epsilon = 2.0`: `p = 2.1482`

Interpretation:

> the same scalar kernel regime remains stable when the scan is extended from interior-only error to explicit homogeneous Dirichlet and Neumann controls.

## Current Verdict

The operator-level verdict is now stronger than the original March cubic-only note.

Shortest honest summary:

> the scalar branch has bounded continuum control beyond the clean interior cubic scan.

More precisely:

- the cubic interior baseline remains a clean Laplacian-control scan
- the local-kernel point-set control survives mild geometric disorder in a bounded higher-resolution window
- reflected Dirichlet and Neumann controls also remain convergent
- broader kernels are still the finite-resolution weak point and should not be treated as the main evidence

## CP1 Reading

For the scalar branch only, the repository now supports this bounded CP1 statement:

> in the local-kernel regime, HAOS-IIP's scalar operator admits refinement-stable continuum control across the cubic interior scan, a weakly perturbed 3D point-set control, and explicit homogeneous boundary-condition families.

This is still not:

- a metric claim
- a curvature claim
- a universality claim across arbitrary kernels or substrates
- a field-theory or spacetime claim

It is a stronger scalar operator-control result, and it is the right base from which to start `CP2`.
