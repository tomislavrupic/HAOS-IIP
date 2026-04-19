# HAOS_Incidence_Theorem_T4_v1

Status:

- structural derivation note
- purpose: execute `T4` of the HAOS-to-harmonic derivation program by deriving adjacent-grade incidence maps from graded recoverability and identifying the no-skip condition that becomes the chain identity in normalized compatibility coordinates

Status labels:

- `[D]` derived in this note under the stated assumptions
- `[A]` assumption introduced for the theorem
- `[PS]` proof sketch / normalization step rather than a fully closed theorem

## 1. Statement

We execute `T4` in the strongest bounded form currently justified by the program:

> multi-grade coherence requires linear adjacency maps
>
> $$
> d_k : V_k \to V_{k+1},
> $$
>
> such that all adjacent-grade consistency constraints factor through them.

Under strict adjacent-grade locality and a normalized identification of the active compatibility channel on the middle grade, the no-skip condition becomes the standard chain / cochain identity

$$
d_{k+1}\circ d_k = 0.
$$

This is the point where discrete complex structure enters the derivation. No simplicial complex, embedding, metric, or continuum geometry is assumed.

## 2. Setup

Let

$$
\mathcal V = \bigoplus_k V_k
$$

be a graded real inner-product space.

Input from `T1`:

$$
D(\mathbf u)=\frac{1}{2}\langle \mathbf u,Q\mathbf u\rangle
$$

near a coherent state, with `Q` symmetric and positive semidefinite.

### Assumption 2.1. Strict adjacent-grade locality

`[A]`

`H4` and `H7` are sharpened in the minimal way needed for `T4`:

- primitive cross-grade defect terms occur only between adjacent grades `V_k` and `V_{k+1}`
- no primitive term directly couples `V_k` to `V_{k+2}` or beyond

Under this assumption, the quadratic form has block-tridiagonal structure

$$
\langle \mathbf u,Q\mathbf u\rangle
=
\sum_k \langle u_k,Q_k u_k\rangle
+
2\sum_k b_k(u_k,u_{k+1}),
$$

where each `b_k` is a continuous bilinear form

$$
b_k : V_k\times V_{k+1}\to \mathbb R.
$$

Interpretation:

- `Q_k` measures same-grade recoverability defect
- `b_k` measures failure of grade-`k` data to cohere with grade-`k+1` data

`T2` and `T3` are not algebraically required for the existence of the graded maps below, but they fix the programmatic meaning of these maps as higher-grade analogues of the scalar relational operator and its recovery flow.

## 3. Step 1: adjacent-grade consistency factors through linear maps

### Proposition 3.1. Existence of incidence maps

`[D]`

For each adjacent pair `(V_k,V_{k+1})` there exists a unique bounded linear map

$$
d_k : V_k\to V_{k+1}
$$

such that

$$
b_k(u_k,u_{k+1})=\langle d_k u_k,u_{k+1}\rangle_{V_{k+1}}
$$

for all `u_k\in V_k` and `u_{k+1}\in V_{k+1}`.

#### Proof

Fix `u_k`. The map

$$
u_{k+1}\mapsto b_k(u_k,u_{k+1})
$$

is a continuous linear functional on `V_{k+1}`. By the Riesz representation theorem there is a unique vector `d_k u_k\in V_{k+1}` representing that functional. Linearity in `u_k` follows from bilinearity of `b_k`, and boundedness follows from continuity of `b_k`. `\square`

### Corollary 3.2. Factorization of adjacent-grade consistency

`[D]`

Every adjacent-grade consistency constraint in the quadratic recoverability law factors through `d_k` and its adjoint `d_k^*`.

Interpretation:

- graded compatibility is not an extra geometry laid on top of HAOS
- it is the linear form taken by adjacent-grade recoverability constraints once the defect is quadratic

## 4. Step 2: eliminating the middle grade exposes the no-skip condition

Now consider three consecutive grades `V_k`, `V_{k+1}`, and `V_{k+2}`. The part of the defect involving the middle grade is

$$
\Phi(u_k,u_{k+1},u_{k+2})
=
\frac{1}{2}\langle u_{k+1},Q_{k+1}u_{k+1}\rangle
+
\langle d_k u_k,u_{k+1}\rangle
+
\langle u_{k+1},d_{k+1}^*u_{k+2}\rangle.
$$

### Assumption 4.1. Active-channel nondegeneracy

`[A]`

On the compatibility channel of `V_{k+1}` that is seen by the adjacent-grade terms, `Q_{k+1}` is invertible or at least pseudoinvertible with Moore-Penrose inverse `Q_{k+1}^+`.

This is the minimum needed to eliminate the middle grade without ambiguity on the active subspace.

### Proposition 4.2. Induced non-adjacent coupling

`[D]`

After minimizing over `u_{k+1}`, the effective defect acquires the mixed term

$$
-\langle d_k u_k,Q_{k+1}^+d_{k+1}^*u_{k+2}\rangle.
$$

#### Proof

The Euler-Lagrange equation for `u_{k+1}` on the active channel is

$$
Q_{k+1}u_{k+1}+d_k u_k+d_{k+1}^*u_{k+2}=0.
$$

Solving with the pseudoinverse gives the minimizing value

$$
u_{k+1}^{\rm opt}
=
-Q_{k+1}^+(d_k u_k+d_{k+1}^*u_{k+2})
$$

on that channel. Substituting back yields the Schur-complement correction

$$
-\frac{1}{2}\langle d_k u_k+d_{k+1}^*u_{k+2},Q_{k+1}^+(d_k u_k+d_{k+1}^*u_{k+2})\rangle,
$$

whose mixed component is exactly

$$
-\langle d_k u_k,Q_{k+1}^+d_{k+1}^*u_{k+2}\rangle.
$$

`\square`

### Corollary 4.3. The honest no-skip condition

`[D]`

If recoverability is strictly adjacent-grade in the sense of Assumption 2.1 even after optimal elimination of intermediate grades, then the induced mixed term must vanish identically. Therefore

$$
d_{k+1}Q_{k+1}^+d_k=0
$$

on the active compatibility channel.

This is the exact operator statement forced by no-skip recoverability before any normalization.

## 5. Step 3: normalized compatibility coordinates give chain form

The previous section already gives the real structural content of `T4`: adjacent-grade coherence factors through linear maps, and no-skip recoverability kills the two-step composition after passage through the middle-grade compatibility metric.

### Proposition 5.1. Chain identity in normalized coordinates

`[PS]`

Choose grade-wise compatibility coordinates so that the positive operator on the active middle channel is absorbed into the identification of `V_{k+1}` with its dual. In those normalized coordinates, Corollary 4.3 becomes

$$
d_{k+1}\circ d_k = 0.
$$

#### Proof sketch

Corollary 4.3 says that the genuine no-skip condition is

$$
d_{k+1}Q_{k+1}^+d_k=0.
$$

If the inner product on the active compatibility channel is replaced by the equivalent normalized inner product induced by `Q_{k+1}^+`, then the factor `Q_{k+1}^+` is absorbed into the grade-`k+1` identification between vectors and covectors. In that normalized picture, the composite takes the standard chain form and reads

$$
d_{k+1}d_k=0.
$$

This is a change of compatibility coordinates, not an extra physical axiom. `\square`

## 6. What is proved and what is not

### Fully proved in this note

- `[D]` graded recoverability induces adjacent-grade bilinear couplings
- `[D]` each such coupling is represented by a unique linear incidence map `d_k`
- `[D]` strict no-skip recoverability implies the operator condition `d_{k+1}Q_{k+1}^+d_k=0`

### Still a normalization step rather than a closed theorem

- `[PS]` the textbook chain identity `d_{k+1}d_k=0` requires choosing compatibility coordinates adapted to the active middle-grade quadratic law

That is still enough to justify the programmatic statement that a cochain-complex-type structure emerges from multi-grade recoverability rather than being inserted by hand.

## 7. Consequence for the next theorem

`T5` can now be stated honestly in the following form:

- once the incidence maps are present, the minimal graded defect is built from upward failure `\|d_k\omega\|^2` and downward failure `\|d_{k-1}^*\omega\|^2`
- after normalized compatibility coordinates are fixed, this becomes the usual Hodge-defect law

So `T4` is the exact point where HAOS recoverability first compels a discrete complex rather than merely suggesting one.
