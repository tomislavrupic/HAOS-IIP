# Periodic / Twisted L1 Experiment

## Purpose

Test whether the HAOS-IIP edge/Hodge branch develops a low-lying coexact or divergence-free family on a topologically nontrivial substrate, or whether the low edge spectrum is still dominated by scalar leakage from the exact sector.

The direct target is:

- open box vs periodic torus
- trivial flux vs nontrivial flux
- full `L1` spectrum vs coexact/divergence-free projected spectrum

## Operator Construction

The open baseline uses the existing discrete Hodge operator

$$
L_1 = d_0 d_0^* + d_1^* d_1
$$

on the weighted cubic complex.

For the periodic branch, the substrate is a cubic 3-torus with oriented nearest-neighbor edges in the `x`, `y`, and `z` directions. The edge weight is kept kernel-derived:

$$
w_e = \exp\left(-\frac{h^2}{2\varepsilon}\right), \qquad h = \frac{1}{n}.
$$

A background `U(1)` connection is introduced in a discrete Landau-gauge form with flux quantum `m` through each `xy` plaquette:

$$
\phi = \frac{2\pi m}{n^2},
$$

$$
U_x(i,j,k) =
\begin{cases}
1, & i < n-1, \\
\exp(-i \phi n j), & i = n-1,
\end{cases}
\qquad
U_y(i,j,k) = \exp(i \phi i),
\qquad
U_z(i,j,k) = 1.
$$

The covariant node-edge differential is

$$
(d_{0,A} \psi)_e = \sqrt{w_e}\,(U_e \psi_{\mathrm{head}(e)} - \psi_{\mathrm{tail}(e)}).
$$

A compatible face-edge operator `d_{1,A}` is built by transporting boundary edge values back to a chosen face basepoint before summing the oriented boundary. The periodic/twisted edge operator is then

$$
L_{1,A} = d_{0,A} d_{0,A}^* + d_{1,A}^* d_{1,A}.
$$

The coexact/divergence-free projected operator is obtained by restricting `L_{1,A}` to

$$
\ker d_{0,A}^*.
$$

Numerically, if `Q` is an orthonormal basis for `\ker d_{0,A}^*`, the projected operator is

$$
L_{1,A}^{\mathrm{div}} = Q^* L_{1,A} Q.
$$

This projection removes exact-gradient contamination but still includes harmonic torus modes when they exist.

## Numerical Setup

- `\varepsilon = 0.2`
- periodic sizes: `n = 4, 5`
- flux sectors: `m = 0, 1`
- open comparison: `n = 5`
- diagnostics per low mode:
  - eigenvalue
  - exact fraction
  - coexact/divergence-free fraction
  - `||d_0^* a||`
  - `||d_1 a||`
  - inverse participation ratio
  - qualitative support pattern

## Low-Mode Results

### Open baseline (`n = 5`)

The first three open-box modes are exactly gradient-like:

| mode | eigenvalue | exact | coexact | `||d_0^* a||` | `||d_1 a||` | pattern |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | 0.326713 | 1.000 | 0.000 | 0.5716 | 0.0000 | delocalized exact `x`-gradient |
| 1 | 0.326713 | 1.000 | 0.000 | 0.5716 | 0.0000 | delocalized exact `y`-gradient |
| 2 | 0.326713 | 1.000 | 0.000 | 0.5716 | 0.0000 | delocalized exact `z`-gradient |

This is the scalar-leakage baseline.

### Periodic, trivial flux (`n = 5`, `m = 0`)

The lowest three modes become exactly divergence-free and curl-free torus cycles:

| mode | eigenvalue | exact | coexact | `||d_0^* a||` | `||d_1 a||` | pattern |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | 0.000000 | 0.000 | 1.000 | 0.0000 | 0.0000 | delocalized harmonic `y`-cycle |
| 1 | 0.000000 | 0.000 | 1.000 | 0.0000 | 0.0000 | delocalized harmonic `x`-cycle |
| 2 | 0.000000 | 0.000 | 1.000 | 0.0000 | 0.0000 | delocalized harmonic `z`-cycle |

So the periodic topology creates a genuine low non-exact family that does not exist on the open box.

### Periodic, nontrivial flux (`n = 5`, `m = 1`)

The harmonic zero modes are lifted, but the lowest modes remain mostly divergence-free:

| mode | eigenvalue | exact | coexact | `||d_0^* a||` | `||d_1 a||` | pattern |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| 0 | 0.013057 | 0.025 | 0.975 | 0.0850 | 0.0764 | mixed `x`-biased support |
| 1 | 0.111287 | 0.194 | 0.806 | 0.2344 | 0.2373 | divergence-free `x`-biased circulation |
| 2 | 0.220336 | 0.000 | 1.000 | 0.0000 | 0.4694 | divergence-free `z`-biased circulation |
| 3 | 0.221528 | 0.097 | 0.903 | 0.3099 | 0.3542 | divergence-free `x`-biased circulation |

The coexact/divergence-free projected spectrum begins at

$$
0.027359,\ 0.220336,\ 0.243475,\ 0.440293, \dots
$$

This shows that the low divergence-free branch survives projection and is not an artifact of exact scalar gradients.

### Persistence across sizes

For `n = 4`, `m = 1`, the same behavior persists:

| mode | eigenvalue | exact | coexact | `||d_0^* a||` | `||d_1 a||` |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.014110 | 0.020 | 0.980 | 0.0906 | 0.0769 |
| 1 | 0.174639 | 0.131 | 0.869 | 0.2541 | 0.3317 |
| 3 | 0.319646 | 0.000 | 1.000 | 0.0000 | 0.5654 |

So the low divergence-free family is not confined to one lattice size.

## Direct Answer

### Does a low family remain mostly coexact?

Yes.

On the periodic substrate, the lowest modes are overwhelmingly divergence-free/coexact-dominated. In the trivial-flux case they appear as harmonic torus cycles. In the nontrivial-flux case they remain mostly coexact after the zero-mode degeneracy is lifted.

### Does it respond cleanly to flux?

Yes.

The trivial-flux harmonic zero modes move to small positive eigenvalues under nontrivial plaquette holonomy, and the projected divergence-free spectrum shifts accordingly.

### Does it persist across lattice sizes?

Yes, at least for the tested sizes `n = 4` and `n = 5`.

### Is it plausibly vector/gauge-like, or still scalar contamination?

It is plausibly vector-like at the operator level.

The open-box low modes are exactly gradient-like, while the periodic/twisted branch produces a low family that remains mostly divergence-free and survives the `\ker d_0^*` projection. That is explicit non-scalar structure.

However, this is not yet a clean propagating Maxwell-like family. In trivial flux the lowest modes are harmonic topological cycles, and in nontrivial flux the low family is still a mixture of near-harmonic and circulation-type modes rather than a fully isolated transverse continuum branch.

## Verdict

The periodic/twisted `L1` experiment clears the main falsification test for the current edge operator.

The edge sector is not just scalar leakage from the node branch.

What is present is:

- a robust low divergence-free family on the periodic substrate
- clean flux sensitivity of that family
- persistence across the tested lattice sizes

What is not present yet is:

- a clearly separated propagating Maxwell-like band
- an unambiguous massless vector continuum limit
- any basis for particle-level claims

## Current verdict

- what worked: periodic topology plus plaquette holonomy exposed a robust low divergence-free `L1` family that is absent on the open box
- what did not appear: a clean propagating transverse band distinct from harmonic torus modes and mixed circulation states
- what must be tried next: flux sweeps, defect or puncture backgrounds, and then only after that a Dirac-type branch coupled to the same edge connection
