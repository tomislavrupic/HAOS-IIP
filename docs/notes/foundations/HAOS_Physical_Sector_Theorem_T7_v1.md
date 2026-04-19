# HAOS_Physical_Sector_Theorem_T7_v1

Status:

- structural derivation note
- purpose: execute `T7` of the HAOS-to-harmonic derivation program by deriving the physically active sector from the `T6` exact / harmonic / coexact decomposition, rather than treating the coexact restriction as a post hoc modeling choice

Status labels:

- `[D]` derived in this note under the stated assumptions
- `[A]` assumption introduced for the theorem

## 1. Statement

We execute `T7` in the bounded form already justified by `T6`:

> if HAOS requires redundancy removal before a mode is called dynamically meaningful, then the physical grade-`k` sector is not the raw space `V_k` but the quotient by exact redundancy and harmonic neutral content.

In the finite-dimensional toy setting, that quotient has the canonical orthogonal representative

$$
P_k
:=
\ker(d_{k-1}^*) \cap (\mathcal H^k)^\perp
=
\operatorname{im} d_k^*,
$$

and the derived active operator on that sector is

$$
T_k
:=
\Delta_k\big|_{P_k}
=
d_k^*d_k\big|_{P_k}.
$$

For `1`-forms this becomes

$$
T
=
d_1^* d_1\big|_{\ker(d_0^*) \cap (\mathcal H^1)^\perp},
$$

which is exactly the coexact transverse restriction used in HAOS-IIP.

## 2. Inputs from earlier steps

From `T6`:

$$
V_k
=
\operatorname{im} d_{k-1}
\oplus
\mathcal H^k
\oplus
\operatorname{im} d_k^*
$$

as an orthogonal direct sum, with

$$
\mathcal H^k = \ker \Delta_k = \ker d_k \cap \ker d_{k-1}^*.
$$

From `T5`:

$$
\Delta_k = d_{k-1}d_{k-1}^* + d_k^* d_k.
$$

## 3. Physical meaning assumptions

### Assumption 3.1. Redundancy quotient principle

`[A]`

If two grade-`k` states differ by an exact component in `\operatorname{im} d_{k-1}`, they represent the same physical content for the purposes of the active sector.

Reason:

- exact directions are the longitudinal / redundancy directions identified in `T6`

### Assumption 3.2. Neutral harmonic exclusion

`[A]`

Harmonic directions may encode persistent topological or background information, but they are not counted as active propagating content in the present dynamical sector.

Reason:

- by `T6`, harmonic states are exactly `\ker \Delta_k`, so they are zero-defect neutral directions under the minimal graded recovery law

These two assumptions are exactly the conditional bridge already stated in the derivation program: physical meaning requires redundancy removal before raw operator content is called propagation.

## 4. The physical quotient

### Definition 4.1. Physical grade-`k` state space

`[D]`

Under Assumptions 3.1-3.2, the active physical sector is the quotient

$$
\mathcal P_k
:=
V_k / \bigl(\operatorname{im} d_{k-1} \oplus \mathcal H^k\bigr).
$$

Interpretation:

- quotient by `\operatorname{im} d_{k-1}` removes redundancy
- quotient by `\mathcal H^k` removes neutral zero-defect directions from the active spectrum

This is the honest physical object. The coexact subspace will appear next as its canonical representative, not as a substitute for it.

## 5. Canonical representative of the quotient

### Proposition 5.1. Unique coexact representative

`[D]`

Every class in `\mathcal P_k` has a unique representative in

$$
P_k
:=
\ker(d_{k-1}^*) \cap (\mathcal H^k)^\perp.
$$

#### Proof

By `T6`, every `\omega\in V_k` has a unique orthogonal decomposition

$$
\omega = \omega_{\rm ex} + \omega_{\rm har} + \omega_{\rm coex}
$$

with

$$
\omega_{\rm ex}\in \operatorname{im} d_{k-1},
\qquad
\omega_{\rm har}\in \mathcal H^k,
\qquad
\omega_{\rm coex}\in \operatorname{im} d_k^*.
$$

Passing to the quotient by `\operatorname{im} d_{k-1}\oplus \mathcal H^k` kills the first two pieces and leaves only `\omega_{\rm coex}`. Since the orthogonal decomposition is unique, this representative is unique. `\square`

### Proposition 5.2. Canonical representative equals the coexact sector

`[D]`

The canonical representative space satisfies

$$
P_k
=
\operatorname{im} d_k^*.
$$

#### Proof

First, `\operatorname{im} d_k^* \subseteq \ker(d_{k-1}^*)` because

$$
d_{k-1}^* d_k^* = (d_k d_{k-1})^* = 0.
$$

Also `\operatorname{im} d_k^* \perp \mathcal H^k` by the orthogonal decomposition from `T6`. Hence

$$
\operatorname{im} d_k^* \subseteq P_k.
$$

Conversely, let `\omega\in P_k`. By `T6`,

$$
\omega = \omega_{\rm ex} + \omega_{\rm har} + \omega_{\rm coex}.
$$

Because `\omega\perp \mathcal H^k`, we must have `\omega_{\rm har}=0`. Because `d_{k-1}^*\omega=0`, and because

$$
d_{k-1}^*\omega_{\rm har}=0,
\qquad
d_{k-1}^*\omega_{\rm coex}=0,
$$

we get

$$
d_{k-1}^*\omega_{\rm ex}=0.
$$

Write `\omega_{\rm ex}=d_{k-1}u`. Then

$$
0
=
\langle u,d_{k-1}^*d_{k-1}u\rangle
=
\|d_{k-1}u\|^2
=
\|\omega_{\rm ex}\|^2.
$$

So `\omega_{\rm ex}=0`. Therefore `\omega=\omega_{\rm coex}\in \operatorname{im} d_k^*`. Hence `P_k\subseteq \operatorname{im} d_k^*`. `\square`

So the coexact sector is not chosen after the fact. It is the unique orthogonal representative of the physical quotient.

## 6. Reduction of the recoverability operator

### Proposition 6.1. `P_k` is invariant under the active operator

`[D]`

For `\omega\in P_k`, one has `\Delta_k\omega\in P_k`.

#### Proof

By Proposition 5.2, `P_k = \operatorname{im} d_k^*`. So there exists `\beta\in V_{k+1}` with `\omega=d_k^*\beta`. Then

$$
\Delta_k\omega
=
d_{k-1}d_{k-1}^*\omega + d_k^*d_k\omega
=
0 + d_k^*d_k\omega
\in
\operatorname{im} d_k^*
=
P_k.
$$

`\square`

### Proposition 6.2. Operator reduction on the physical sector

`[D]`

On `P_k`, the grade-`k` recoverability operator reduces to

$$
T_k
:=
\Delta_k\big|_{P_k}
=
d_k^*d_k\big|_{P_k}.
$$

#### Proof

If `\omega\in P_k`, then `d_{k-1}^*\omega=0` by definition. Therefore

$$
\Delta_k\omega
=
d_{k-1}d_{k-1}^*\omega + d_k^*d_k\omega
=
d_k^*d_k\omega.
$$

`\square`

Interpretation:

- the exact sector is removed by quotient
- the harmonic sector is removed because it is neutral
- the surviving active dynamics is the coexact / transverse branch

## 7. Physical-sector theorem

### Theorem 7.1. Derived physical-sector restriction

`[D]`

Under Assumptions 3.1-3.2, the active physical grade-`k` sector is canonically represented by

$$
P_k
=
\ker(d_{k-1}^*) \cap (\mathcal H^k)^\perp
=
\operatorname{im} d_k^*,
$$

and its derived active operator is

$$
T_k
=
d_k^*d_k\big|_{\ker(d_{k-1}^*) \cap (\mathcal H^k)^\perp}.
$$

This is the physical-sector theorem promised by the derivation program.

### Corollary 7.2. The `1`-form transverse restriction

`[D]`

For the minimal one-form sector

$$
V_0 \xrightarrow{d_0} V_1 \xrightarrow{d_1} V_2,
$$

the active physical operator is

$$
T
=
d_1^* d_1\big|_{\ker(d_0^*) \cap (\mathcal H^1)^\perp}.
$$

This is exactly the restricted coexact transverse operator already used in the HAOS-IIP numerical architecture.

## 8. Why this is not post hoc

The logical order is now:

```text
incidence maps
    -> Hodge-type recoverability operator
    -> exact / harmonic / coexact decomposition
    -> quotient by redundancy and neutral modes
    -> canonical coexact representative
    -> restricted physical operator
```

So the transverse restriction is not an aesthetic projection added after seeing the numerics. It is the canonical representative of the physical quotient once the HAOS interpretation of redundancy and neutral sectors is imposed.

## 9. Relation to the theorem stack

- `T6` is executed in [HAOS_Sector_Decomposition_Theorem_T6_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/notes/foundations/HAOS_Sector_Decomposition_Theorem_T6_v1.md)
- `T5` is executed in [HAOS_Hodge_Defect_Theorem_T5_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/notes/foundations/HAOS_Hodge_Defect_Theorem_T5_v1.md)
- the larger one-form toy-model realization remains in [HAOS_Relational_Scalar_Recoverability_Defect_and_First_Dirichlet_Theorem_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/notes/foundations/HAOS_Relational_Scalar_Recoverability_Defect_and_First_Dirichlet_Theorem_v1.md)

## 10. Shortest honest conclusion

`[D]`

Once `T6` has separated exact, harmonic, and coexact content, the physical active sector is the quotient by exact redundancy and harmonic neutrality. Its canonical orthogonal representative is the coexact sector, and the induced active operator is the restricted coexact operator

$$
T_k = d_k^*d_k\big|_{\ker(d_{k-1}^*) \cap (\mathcal H^k)^\perp}.
$$

That is the bounded `T7` result.
