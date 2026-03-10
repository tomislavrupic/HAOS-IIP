# Interaction Attractor Spectrum v1

## Purpose

Test whether hydrogen-like discrete bound states can appear as recoverable attractors of a continuous interaction field, rather than by imposing quantization as a primitive axiom.

The target question is:

> can a continuous inverse-distance interaction substrate produce a discrete `-1/n^2` ladder through attractor dynamics alone?

## Minimal Model

Use the reduced radial interaction operator

$$
H_\ell = -D \frac{d^2}{dr^2} + \frac{D \ell(\ell+1)}{r^2} - \frac{\alpha}{r},
$$

acting on the reduced radial field `chi(r) = r R(r)`.

Instead of treating bound states as static eigenmodes only, define a normalized interaction flow

$$
\frac{d\chi}{d\tau} = -\big(H_\ell \chi - \lambda[\chi] \chi \big),
\qquad
\lambda[\chi] = \frac{\langle \chi, H_\ell \chi \rangle}{\langle \chi, \chi \rangle}.
$$

This has three useful properties:

- the norm is preserved by renormalization
- the Rayleigh energy `lambda[chi]` decreases along the flow
- fixed points satisfy `H_\ell chi = E chi`

So stable coherent states are attractors of the interaction flow.

## Stability Criterion

A state counts as a recoverable attractor if it satisfies:

- bounded profile on the radial grid
- convergence under the flow:
  `chi(tau) -> chi_*`
- return after perturbation within the same allowed sector

Numerically the residual criterion is

$$
\|H_\ell \chi - \lambda[\chi]\chi\| \to 0.
$$

## Why Projection Is Needed

The unconstrained flow has a sharp bias:

- generic initial data always descends toward the lowest available energy state

So for the s-wave channel the unconstrained attractor is only the `n = 1` state.

To test whether a discrete ladder exists at all, the experiment uses projected recoverable sectors. After the first attractor is identified, the next run removes the lower sector and repeats the flow. In a finite modal basis this is implemented by allowing only coefficients above the already recovered states.

This matters:

- the discrete ladder is present in the operator
- attractor dynamics selects one state per allowed sector
- without sector separation, higher states are not generic attractors

## Numerical Implementation

File:

- `numerics/simulations/interaction_attractor_spectrum.py`

Method:

- finite-difference discretization of the radial operator
- modal reduction onto the first `12` radial modes
- normalized semigroup update in coefficient space:

$$
c_j \mapsto e^{-dt E_j} c_j,
$$

followed by projection and normalization

Parameters used:

- `D = 1.0`
- `alpha = 1.0`
- `r_min = 0.001`
- `r_max = 80.0`
- `n_grid = 1200`
- `states = 4`
- `modal_basis_size = 12`
- `dt = 6.0`
- `max_steps = 400`
- `recovery_noise = 0.12`

Artifacts:

- results: `data/20260310_105044_interaction_attractor_spectrum.json`
- latest: `data/interaction_attractor_spectrum_latest.json`
- plots:
  - `plots/20260310_105044_interaction_attractor_energy_descent.png`
  - `plots/20260310_105044_interaction_attractor_scaling.png`
  - `plots/20260310_105044_interaction_attractor_recovery.png`
  - `plots/20260310_105044_interaction_attractor_modes.png`

## Attractor Results

Recovered s-wave attractor energies:

| `n` | attractor energy | target `-1/(4 n^2)` | relative error | recovered overlap |
| --- | --- | --- | --- | --- |
| 1 | `-0.249450` | `-0.250000` | `0.22%` | `1.000000` |
| 2 | `-0.062435` | `-0.062500` | `0.10%` | `1.000000` |
| 3 | `-0.027758` | `-0.027778` | `0.07%` | `1.000000` |
| 4 | `-0.015267` | `-0.015625` | `2.29%` | `1.000000` |

The scaling diagnostic is nearly constant:

$$
n^2 |E_n| = \{0.249450,\;0.249742,\;0.249826,\;0.244268\}.
$$

Mean `0.248321`, standard deviation `0.002344`.

So the attractor energies follow the same resolved `-1/n^2` ladder as the static bound-state calculation.

## Mechanism

The discreteness comes from three ingredients:

1. radial inverse-distance interaction geometry
2. quadratic second-order propagation operator
3. boundary-constrained recoverability

The attractor flow does not invent a new spectrum. It stabilizes the discrete bound modes already permitted by the interaction geometry and the boundary conditions.

This leads to a more precise statement:

> discrete hydrogen-like states do emerge from continuous interaction dynamics, but higher members of the ladder are attractors only within appropriately separated recoverable sectors.

## Interpretation

This is a positive but limited result.

Positive:

- a continuous interaction field can support a discrete attractor ladder
- the resulting stability energies follow `E_n propto -1/n^2` in the resolved range
- perturb-and-recover tests return to the same attractor inside each projected sector

Limitation:

- generic unconstrained flow selects only the ground state
- multiple stable attractors do not appear as simultaneous global sinks of one unconstrained phase flow

So the strongest defensible claim is:

> the inverse-distance interaction substrate contains a discrete ladder of recoverable coherent states, and normalized interaction flow can realize them as attractors once the corresponding sectors are separated.

## Current Verdict

- yes, discrete hydrogen-like stability levels emerge from a continuous interaction substrate
- yes, the resolved ladder follows `-1/n^2`
- no, the full unconstrained attractor flow does not give many globally competing bound attractors; it selects the lowest sector unless additional sector separation is imposed
