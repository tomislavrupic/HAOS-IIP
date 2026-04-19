# HAOS_Hilbert_Complex_Upgrade_F4_v1

Status:

- structural refinement note
- purpose: close `F4` of the continuum-execution program by upgrading the finite-dimensional `T6` sector theorem to a closed-range Hilbert-complex statement while preserving the bounded interpretation of the repo

Status labels:

- `[D]` derived in this note under the stated assumptions
- `[A]` assumption imported from earlier notes

## 1. Statement

`T6` was already proved in finite-dimensional inner-product spaces:

$$
V_k
=
\operatorname{im} d_{k-1}
\oplus
\mathcal H^k
\oplus
\operatorname{im} d_k^*.
$$

The open `F4` task was to state the same result in cleaner Hilbert-complex language.

This note closes that in the strongest bounded form currently justified:

> on Hilbert spaces with bounded incidence maps and closed ranges, the same exact / harmonic / coexact decomposition holds as an orthogonal direct sum.

Without closed-range hypotheses, the same theorem still holds with closures on the image sectors.

## 2. Setup

Let

$$
V_{k-1},\ V_k,\ V_{k+1}
$$

be Hilbert spaces, and let

$$
d_{k-1}:V_{k-1}\to V_k,
\qquad
d_k:V_k\to V_{k+1}
$$

be bounded operators satisfying the chain condition

$$
d_k d_{k-1}=0.
$$

Define the harmonic sector

$$
\mathcal H^k := \ker d_k \cap \ker d_{k-1}^*.
$$

As in `T5`, define the grade-`k` recoverability operator

$$
\Delta_k := d_{k-1}d_{k-1}^* + d_k^* d_k.
$$

## 3. Basic Hilbert-space identities

### Proposition 3.1. Orthogonal-complement identities

`[D]`

For bounded operators on Hilbert spaces,

$$
(\overline{\operatorname{im} d_{k-1}})^\perp = \ker d_{k-1}^*,
\qquad
(\overline{\operatorname{im} d_k^*})^\perp = \ker d_k.
$$

This is the standard Hilbert-space adjoint identity.

### Proposition 3.2. Exact and coexact orthogonality

`[D]`

The image sectors are orthogonal:

$$
\overline{\operatorname{im} d_{k-1}}
\perp
\overline{\operatorname{im} d_k^*}.
$$

#### Proof

For `u\in V_{k-1}` and `\beta\in V_{k+1}`,

$$
\langle d_{k-1}u,d_k^*\beta\rangle
=
\langle d_k d_{k-1}u,\beta\rangle
=
0
$$

by the chain condition. Taking closures preserves orthogonality. `\square`

## 4. Weak Hilbert-complex decomposition

### Theorem 4.1. Decomposition with closures

`[D]`

Without any closed-range assumption,

$$
V_k
=
\overline{\operatorname{im} d_{k-1}}
\oplus
\mathcal H^k
\oplus
\overline{\operatorname{im} d_k^*}.
$$

#### Proof

First,

$$
V_k = \ker d_k \oplus \overline{\operatorname{im} d_k^*}
$$

by Proposition 3.1.

Inside `\ker d_k`, the chain condition implies

$$
\overline{\operatorname{im} d_{k-1}} \subseteq \ker d_k.
$$

Also,

$$
(\overline{\operatorname{im} d_{k-1}})^\perp \cap \ker d_k
=
\ker d_{k-1}^* \cap \ker d_k
=
\mathcal H^k.
$$

Therefore

$$
\ker d_k
=
\overline{\operatorname{im} d_{k-1}}
\oplus
\mathcal H^k.
$$

Combining the two splittings gives the stated decomposition. `\square`

This is the natural Hilbert-space upgrade of the finite-dimensional theorem.

## 5. Closed-range upgrade

### Assumption 5.1. Closed ranges

`[A]`

Assume

$$
\operatorname{im} d_{k-1}
\quad\text{and}\quad
\operatorname{im} d_k
$$

are closed.

Because bounded operators have closed range if and only if their adjoints do, this also gives closed range for `d_{k-1}^*` and `d_k^*`.

### Theorem 5.2. Closed-range Hilbert-complex decomposition

`[D]`

Under Assumption 5.1,

$$
V_k
=
\operatorname{im} d_{k-1}
\oplus
\mathcal H^k
\oplus
\operatorname{im} d_k^*.
$$

#### Proof

Apply Theorem 4.1 and remove the closures using closed range. `\square`

So the finite-dimensional `T6` theorem was not a one-off linear-algebra accident. It is the closed-range Hilbert-complex form of the same statement.

## 6. Harmonic sector and the recoverability operator

### Proposition 6.1. Harmonic sector equals `\ker \Delta_k`

`[D]`

In this Hilbert setting,

$$
\mathcal H^k = \ker \Delta_k.
$$

#### Proof

Exactly as in `T6`,

$$
\langle \omega,\Delta_k\omega\rangle
=
\|d_{k-1}^*\omega\|^2 + \|d_k\omega\|^2.
$$

Therefore `\Delta_k\omega=0` implies both channel norms vanish, which is equivalent to

$$
\omega\in \ker d_{k-1}^*\cap \ker d_k = \mathcal H^k.
$$

The reverse inclusion is immediate. `\square`

## 7. Why this preserves the bounded repo interpretation

This note does not turn the repo into an unbounded PDE program. It keeps the same bounded reading:

- bounded incidence maps
- bounded recoverability operator
- only one extra functional-analytic input: closed range, when one wants the exact image sectors without closures

So the price of the Hilbert-complex upgrade is explicit and modest:

- without closed range, one must write closures
- with closed range, one gets the exact direct-sum theorem familiar from Hodge theory

## 8. Consequence for `T6`

The strongest bounded current reading of `T6` is now:

1. in Hilbert-complex form one always has

$$
V_k
=
\overline{\operatorname{im} d_{k-1}}
\oplus
\mathcal H^k
\oplus
\overline{\operatorname{im} d_k^*}
$$

2. if the relevant incidence maps have closed range, then

$$
V_k
=
\operatorname{im} d_{k-1}
\oplus
\mathcal H^k
\oplus
\operatorname{im} d_k^*
$$

3. the harmonic sector is still exactly `\ker \Delta_k`

So `F4` is closed in bounded form: the finite-dimensional theorem is now embedded in standard closed-range Hilbert-complex language.

## 9. Relation to the stack

- `T5` remains the Hodge-defect theorem in [HAOS_Hodge_Defect_Theorem_T5_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/notes/foundations/HAOS_Hodge_Defect_Theorem_T5_v1.md)
- `T6` remains the finite-dimensional sector theorem in [HAOS_Sector_Decomposition_Theorem_T6_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/notes/foundations/HAOS_Sector_Decomposition_Theorem_T6_v1.md)
- this note executes `F4` by giving the closed-range Hilbert-complex upgrade of the same decomposition

## 10. Shortest honest conclusion

`[D]`

The exact / harmonic / coexact split is not confined to the finite-dimensional toy setting. In Hilbert-complex language one always has the decomposition with closed images replaced by closures, and under closed-range hypotheses one recovers the exact orthogonal direct sum

$$
V_k
=
\operatorname{im} d_{k-1}
\oplus
\mathcal H^k
\oplus
\operatorname{im} d_k^*.
$$

That is the bounded `F4` upgrade.
