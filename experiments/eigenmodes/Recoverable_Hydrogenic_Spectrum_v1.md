# Recoverable Hydrogenic Spectrum v1

## Purpose

Test whether a discrete hydrogen-like bound-state ladder can emerge from a continuous interaction substrate without postulating quantization as a primitive axiom.

The working assumptions are:

- a complex disturbance field `psi(x,t)` propagates on a continuous substrate
- the substrate has quadratic dispersion `omega(k) = D k^2`
- the substrate also carries a radial interaction structure `V(r) = -alpha / r`
- stable states are bounded, recoverable coherent modes

## Governing Equation

To match the stationary operator requested for bound-mode analysis, use the sign convention

$$
\partial_t \psi = i \left( D \nabla^2 - V(r) \right) \psi,
\qquad
V(r) = -\frac{\alpha}{r}.
$$

With the modal ansatz

$$
\psi(x,t) = u(x) e^{-i E t},
$$

the time-independent problem becomes

$$
-D \nabla^2 u + V(r) u = E u.
$$

For spherical separation

$$
u(r,\Omega) = R_{n\ell}(r) Y_{\ell m}(\Omega),
$$

and the reduced radial field

$$
\chi(r) = r R_{n\ell}(r),
$$

the radial equation is

$$
-D \chi''(r) + \left[ \frac{D \ell(\ell+1)}{r^2} - \frac{\alpha}{r} \right] \chi(r) = E \chi(r).
$$

Boundary conditions:

- regularity near `r = 0`
- square-integrability of `chi`
- decay as `r -> infinity`

## Structural Origin of Discreteness

For bound states set

$$
E = -D \kappa^2,
\qquad
\kappa > 0.
$$

Then

$$
\chi'' + \left[
-\kappa^2
+ \frac{\alpha}{D r}
- \frac{\ell(\ell+1)}{r^2}
\right] \chi = 0.
$$

Define `rho = 2 kappa r` and write

$$
\chi(\rho) = \rho^{\ell+1} e^{-\rho/2} f(\rho).
$$

The remaining equation is

$$
\rho f'' + (2\ell + 2 - \rho) f' + (\lambda - \ell - 1) f = 0,
\qquad
\lambda = \frac{\alpha}{2 D \kappa}.
$$

Normalizability requires the series for `f` to terminate rather than grow exponentially. That forces

$$
\lambda - \ell - 1 = n_r,
\qquad
n_r = 0,1,2,\dots
$$

so that

$$
\kappa = \frac{\alpha}{2 D (\ell + 1 + n_r)}.
$$

Defining the principal index

$$
n = \ell + 1 + n_r,
$$

gives the discrete ladder

$$
E_n = -\frac{\alpha^2}{4 D n^2}.
$$

So the discreteness comes from three ingredients only:

- quadratic dispersion, which yields a second-order radial operator
- inverse-distance geometry, which creates the `1/r` binding scale
- regularity and decay, which force series termination

No separate quantization postulate is required.

## Numerical Setup

Implementation:

- file: `numerics/simulations/recoverable_hydrogenic_spectrum.py`
- radial finite-difference discretization of the reduced equation
- Dirichlet cutoff on `chi(r)` at small `r_min` and large `r_max`

Parameters used in the first run:

- `D = 1.0`
- `alpha = 1.0`
- `r_min = 0.001`
- `r_max = 80.0`
- `n_grid = 1400`
- `ell = 0, 1, 2`
- `4` bound states per angular channel when available

Artifacts:

- results: `data/20260310_104607_recoverable_hydrogenic_spectrum.json`
- latest: `data/recoverable_hydrogenic_spectrum_latest.json`
- plots:
  - `plots/20260310_104607_recoverable_hydrogenic_energy_scaling.png`
  - `plots/20260310_104607_recoverable_hydrogenic_n2_energy.png`
  - `plots/20260310_104607_recoverable_hydrogenic_modes.png`

## Numerical Results

For the s-wave channel `ell = 0`, the first four bound energies are:

| `n` | numerical `E_n` | predicted `-1/(4 n^2)` | relative error |
| --- | --- | --- | --- |
| 1 | `-0.249466` | `-0.250000` | `0.21%` |
| 2 | `-0.062436` | `-0.062500` | `0.10%` |
| 3 | `-0.027759` | `-0.027778` | `0.07%` |
| 4 | `-0.015267` | `-0.015625` | `2.29%` |

The scaling diagnostic is nearly constant:

$$
n^2 |E_n| = \{0.249466,\; 0.249745,\; 0.249827,\; 0.244269\},
$$

with mean `0.248327` and standard deviation `0.002347`.

The higher-`ell` channels show the same principal-index structure in the resolved range:

- `ell = 1`, `n = 2`: `E = -0.062501`
- `ell = 1`, `n = 3`: `E = -0.027778`
- `ell = 2`, `n = 3`: `E = -0.027778`

So the numerical spectrum reproduces both discreteness and the expected `n`-based scaling to good accuracy before finite-box effects become visible at larger radii.

## Interpretation

The mechanism is structural rather than axiomatic.

The quadratic dispersion creates a recoverable second-order propagation law. The inverse-distance geometry creates a radial binding term with no preferred finite length scale. The regularity and decay conditions then select only a discrete subset of globally recoverable modes. In that sense the ladder is a consequence of bounded coherent wave dynamics on the substrate, not a separately imposed quantization rule.

This is the cleanest statement supported by the current run:

> discrete hydrogen-like energies can emerge from a continuous interaction substrate once recoverable bounded modes are required.

## Minimal Additional Exploration

The present derivation also indicates how the ladder should change under modifications:

- alternative dispersion:
  replacing `D k^2` with a different power changes the radial operator order or scaling and generally breaks the exact `1/n^2` law
- nonlinear recoverability terms:
  adding amplitude-dependent self-restoration shifts the bound-state energies and can destroy the exact principal-index degeneracy
- substrate curvature:
  a curved radial measure changes the effective centrifugal and volume terms, which should split the clean flat-substrate ladder
- coherence-restoration operators:
  extra stabilizing radial terms can preserve discreteness while moving the spectrum away from the pure inverse-distance scaling

These are natural follow-on tests, but they are not needed for the main result.

## Current Verdict

- the minimal substrate equation supports discrete bound eigenmodes
- the spectrum follows `E_n propto -1/n^2` in the resolved low-state range
- the discreteness is generated by dispersion, radial geometry, and boundary conditions rather than a separate quantization postulate
