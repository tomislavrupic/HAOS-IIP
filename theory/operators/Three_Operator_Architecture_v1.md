# Three Operator Architecture

Date: March 9, 2026

## Core statement

The simplest serious HAOS-IIP architecture is no longer:

```text
one Laplacian explains everything
```

It is:

```text
one substrate
one kernel
three operator sectors
```

## Shared substrate and kernel

Start from

$$
A_{ij} = \exp\!\left(-\frac{|x_i-x_j|^2}{2\epsilon_k}\right).
$$

This defines a weighted discrete geometry. Different geometric objects on that same weighted geometry carry different operators.

## Sector 1: scalar / geometry

Acts on:

- node `0-forms`

Operator:

$$
L_0 = D - A = B^T W B.
$$

Role:

- smooth low modes
- diffusion
- scalar continuum limit
- curvature/geometry program

## Sector 2: gauge / vector

Acts on:

- edge `1-forms`

Operator:

$$
L_1 = B W_0 B^T + C^T W_2 C.
$$

Starter approximation before faces:

$$
L_1 \approx B B^T.
$$

Role:

- circulation
- flux
- longitudinal vs transverse structure
- Maxwell-like candidate branch

## Sector 3: fermion / spinor

Acts on:

- doubled or framed node variables
- spinorial toy fields

Minimal block form:

$$
D_H =
\begin{bmatrix}
0 & B^T \\
B & 0
\end{bmatrix}.
$$

Weighted / phase-dressed form:

$$
D_H =
\begin{bmatrix}
0 & B^T U^T \\
U B & 0
\end{bmatrix}.
$$

Role:

- first-order propagation
- sign-symmetric spectrum
- candidate chirality
- defect-bound modes

## Why this is cleaner

The unification target becomes:

```text
interaction fabric
    ->
weighted discrete geometry
    ->
operator hierarchy
```

This is stronger and more realistic than forcing all sectors out of one node Laplacian.

## Immediate consequence

The current 3D eigenmode study validates only stage 1:

```text
kernel -> L0 -> scalar low modes
```

The next explicit target is therefore:

```text
build L1
```

and only after that:

```text
build the minimal D_H branch
```
