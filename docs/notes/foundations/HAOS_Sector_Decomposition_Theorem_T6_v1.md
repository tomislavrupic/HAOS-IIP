# HAOS_Sector_Decomposition_Theorem_T6_v1

Status:

- structural derivation note
- purpose: execute `T6` of the HAOS-to-harmonic derivation program by deriving the exact / harmonic / coexact decomposition from the incidence maps and the `T5` recoverability operator, without assuming Hodge decomposition up front

Status labels:

- `[D]` derived in this note under the stated assumptions
- `[A]` assumption introduced for the theorem

## 1. Statement

We execute `T6` in the strongest bounded form already justified by `T4` and `T5`:

> the grade-`k` space decomposes as
>
> $$
> V_k
> =
> \operatorname{im} d_{k-1}
> \oplus
> \mathcal H^k
> \oplus
> \operatorname{im} d_k^*,
> $$
>
> where
>
> $$
> \mathcal H^k
> :=
> \ker d_k \cap \ker d_{k-1}^*.
> $$

Interpretation:

- exact sector = longitudinal or redundancy content
- harmonic sector = neutral persistent content
- coexact sector = transverse or circulation-like content

This is the point where the sector language becomes a theorem of the derived operator stack rather than an imported vocabulary.

## 2. Inputs from earlier steps

From `T4`:

$$
d_{k-1}:V_{k-1}\to V_k,
\qquad
d_k:V_k\to V_{k+1},
\qquad
d_k d_{k-1}=0
$$

in normalized compatibility coordinates.

From `T5`:

$$
\Delta_k = d_{k-1}d_{k-1}^* + d_k^*d_k
$$

and

$$
\langle \omega,\Delta_k\omega\rangle
=
\|d_{k-1}^*\omega\|^2 + \|d_k\omega\|^2.
$$

### Assumption 2.1. Finite-dimensional working setting

`[A]`

For this theorem we work in finite-dimensional inner-product spaces `V_k`.

Reason:

- finite dimension is enough for the current toy derivation program
- it lets orthogonal complements and direct-sum decompositions be used without bringing in extra functional analysis assumptions

## 3. Exact and coexact sectors are orthogonal

### Proposition 3.1. Orthogonality of `im d_{k-1}` and `im d_k^*`

`[D]`

The exact and coexact sectors are orthogonal:

$$
\operatorname{im} d_{k-1}\perp \operatorname{im} d_k^*.
$$

#### Proof

Take arbitrary `u\in V_{k-1}` and `\beta\in V_{k+1}`. Then

$$
\langle d_{k-1}u,d_k^*\beta\rangle
=
\langle d_k d_{k-1}u,\beta\rangle
=
0
$$

by the chain identity `d_k d_{k-1}=0`. `\square`

So the sum

$$
W_k := \operatorname{im} d_{k-1} + \operatorname{im} d_k^*
$$

is automatically a direct orthogonal sum:

$$
W_k = \operatorname{im} d_{k-1}\oplus \operatorname{im} d_k^*.
$$

## 4. The harmonic sector is the orthogonal complement

### Proposition 4.1. Orthogonal-complement characterization

`[D]`

The harmonic space satisfies

$$
\mathcal H^k = W_k^\perp.
$$

#### Proof

Let `\omega\in V_k`. Then `\omega\in W_k^\perp` exactly when

$$
\langle \omega,d_{k-1}u\rangle = 0
\quad \forall u\in V_{k-1}
$$

and

$$
\langle \omega,d_k^*\beta\rangle = 0
\quad \forall \beta\in V_{k+1}.
$$

By adjointness these are equivalent to

$$
d_{k-1}^*\omega = 0
\qquad\text{and}\qquad
d_k\omega = 0.
$$

Therefore

$$
W_k^\perp
=
\ker d_{k-1}^*\cap \ker d_k
=
\mathcal H^k.
$$

`\square`

Interpretation:

- harmonic states are precisely those grade-`k` states invisible to both adjacent failure channels
- so they are not added as a third primitive structure; they appear as what remains after exact and coexact channels are removed

## 5. Sector decomposition theorem

### Theorem 5.1. Exact / harmonic / coexact decomposition

`[D]`

In the finite-dimensional working setting,

$$
V_k
=
\operatorname{im} d_{k-1}
\oplus
\mathcal H^k
\oplus
\operatorname{im} d_k^*.
$$

This is an orthogonal direct sum.

#### Proof

From Proposition 3.1, `W_k = \operatorname{im} d_{k-1}\oplus \operatorname{im} d_k^*` is an orthogonal direct sum. Finite-dimensional inner-product spaces decompose as

$$
V_k = W_k \oplus W_k^\perp.
$$

By Proposition 4.1,

$$
W_k^\perp = \mathcal H^k.
$$

Substituting gives

$$
V_k
=
\operatorname{im} d_{k-1}
\oplus
\mathcal H^k
\oplus
\operatorname{im} d_k^*.
$$

Because `\mathcal H^k = W_k^\perp`, each summand is orthogonal to the other two. `\square`

## 6. Harmonic sector equals the kernel of the recoverability operator

### Proposition 6.1. Harmonic states are exactly zero-defect directions

`[D]`

The harmonic sector coincides with the kernel of the grade-`k` recoverability operator:

$$
\mathcal H^k = \ker \Delta_k.
$$

#### Proof

If `\omega\in \mathcal H^k`, then

$$
d_{k-1}^*\omega = 0,
\qquad
d_k\omega = 0.
$$

So

$$
\Delta_k\omega
=
d_{k-1}d_{k-1}^*\omega + d_k^*d_k\omega
=
0.
$$

Hence `\mathcal H^k \subseteq \ker \Delta_k`.

Conversely, if `\omega\in \ker \Delta_k`, then

$$
0
=
\langle \omega,\Delta_k\omega\rangle
=
\|d_{k-1}^*\omega\|^2 + \|d_k\omega\|^2.
$$

Both norms must therefore vanish, so

$$
d_{k-1}^*\omega = 0,
\qquad
d_k\omega = 0.
$$

Thus `\omega\in \mathcal H^k`, proving `\ker \Delta_k \subseteq \mathcal H^k`. `\square`

This is the operational meaning of the harmonic sector:

- harmonic states are exactly the grade-`k` directions with zero derived recoverability defect
- they are neutral under the minimal graded recovery law

## 7. Why this does not assume decomposition up front

This note does not begin by postulating three sectors and then naming them exact, harmonic, and coexact.

Instead the logic is:

```text
incidence maps from T4
    -> minimal graded operator from T5
    -> orthogonality of image sectors
    -> harmonic space as orthogonal complement
    -> full decomposition
```

So the decomposition is not a borrowed slogan. It is the linear-algebraic consequence of the chain identity plus the derived recoverability operator in the finite-dimensional toy setting.

## 8. Relation to the theorem stack

- `T4` is executed in [HAOS_Incidence_Theorem_T4_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/notes/foundations/HAOS_Incidence_Theorem_T4_v1.md)
- `T5` is executed in [HAOS_Hodge_Defect_Theorem_T5_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/notes/foundations/HAOS_Hodge_Defect_Theorem_T5_v1.md)
- the closed-range Hilbert-complex upgrade is now executed in [HAOS_Hilbert_Complex_Upgrade_F4_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/notes/foundations/HAOS_Hilbert_Complex_Upgrade_F4_v1.md)
- the earlier one-form toy-model realization remains in [HAOS_Relational_Scalar_Recoverability_Defect_and_First_Dirichlet_Theorem_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/notes/foundations/HAOS_Relational_Scalar_Recoverability_Defect_and_First_Dirichlet_Theorem_v1.md)

## 9. Shortest honest conclusion

`[D]`

Once `T4` gives the incidence maps and `T5` gives the minimal graded recoverability operator, the exact / harmonic / coexact split is forced in the finite-dimensional toy setting:

$$
V_k
=
\operatorname{im} d_{k-1}
\oplus
\ker \Delta_k
\oplus
\operatorname{im} d_k^*.
$$

That is the bounded `T6` result.
