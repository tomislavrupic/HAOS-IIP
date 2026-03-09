# OPERATOR_STRUCTURE

## Operator tower

The current HAOS-IIP architecture is:

```text
same substrate
same kernel
different operator lifts
```

This replaces the earlier compressed reading that one node Laplacian might explain all sectors directly.

## Scalar operator layer: `L0`

`[E]` The primary discrete scalar operator acts on node `0-forms`:

$$
L_0 = D - A.
$$

`[E]` It can be written as

$$
L_0 = B^T W B.
$$

`[P]` In the dense limit, the current program expects

$$
L_0 \longrightarrow \mathcal{L}_0 = -\Delta_g + V(x).
$$

This is a scalar operator layer.

## Gauge / vector layer: `L1`

`[P]` Gauge structure should live on oriented edge `1-forms`, not on node scalars.

With node-edge incidence matrix `B`, and eventually edge-face incidence matrix `C`, the discrete Hodge candidate is

$$
L_1 = B W_0 B^T + C^T W_2 C.
$$

`[P]` Before faces are built explicitly, a rough starter is

$$
L_1 \approx B B^T.
$$

This is the first operator where one can seriously ask about:

- transverse modes
- circulation
- flux sectors
- Maxwell-like candidates

## Phase and transport on the scalar branch

`[E]` Legacy operator notes already support:

- complex phase as the minimal interference-capable representation
- global phase redundancy
- translation generator proportional to `-i nabla`

`[E]` The same notes also state that first-derivative structure requires a background field; the minimal homogeneous scalar generator is Laplacian.

## Covariant graph Laplacian

`[P]` A candidate covariant scalar operator is

$$
(L_\theta \psi)_i = \sum_j A_{ij}(\psi_i - U_{ij}\psi_j),
\qquad
U_{ij}=e^{i\theta_{ij}}.
$$

This operator is appropriate for testing:

- connection sensitivity
- holonomy effects
- edge-current structure

But it still acts on node scalars. It is not by itself the full gauge/vector sector.

## Fermion / spinor layer: `D_H`

`[P]` A fermion sector requires a first-order operator.

The minimal block-form toy architecture is

$$
D_H =
\begin{bmatrix}
0 & B^T \\
B & 0
\end{bmatrix},
$$

or in a weighted / phase-dressed version

$$
D_H =
\begin{bmatrix}
0 & B^T U^T \\
U B & 0
\end{bmatrix}.
$$

This is the right direction for:

- sign-symmetric spectrum
- first-order propagation
- candidate chiral structure
- defect-bound modes

## Current bottleneck

`[O]` The repository now has a clear scalar branch, but does not yet contain a built:

- weighted edge/Hodge branch `L1`
- Dirac-type branch `D_H`
- topological defect operator with controlled fermion-like localization
