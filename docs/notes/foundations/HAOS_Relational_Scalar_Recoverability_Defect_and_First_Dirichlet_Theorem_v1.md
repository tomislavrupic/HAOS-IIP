# HAOS_Relational_Scalar_Recoverability_Defect_and_First_Dirichlet_Theorem_v1

Status:

- toy derivation note
- purpose: define a concrete HAOS-aligned recoverability defect on relational scalar states, prove the first quadratic / Dirichlet reduction, and only then lift the same logic to graded incidence and Hodge structure

Status labels:

- `[E]` established in the current repository or in standard finite-dimensional linear algebra
- `[D]` derived in this note under the stated assumptions
- `[A]` assumption introduced for the toy derivation

## Abstract

This note executes the first real step of the HAOS-to-harmonic derivation program on a minimal toy model. We begin with scalar states on a finite relational substrate and define a recoverability defect that depends only on local interaction differences. Under explicit HAOS-aligned assumptions of relational invariance, local failure assembly, weak compositionality, smoothness near coherence, and reciprocity, we prove that the defect admits a quadratic leading term around coherent states and that this leading term is a weighted Dirichlet form. The corresponding operator is a weighted graph Laplacian

$$
L_0 = B^\top \Lambda B.
$$

Only after this scalar step is fixed do we lift the same logic to a graded incidence complex. There the minimal one-form recoverability defect becomes

$$
\mathcal R_1(\omega) = \frac{1}{2}\|d_0^\ast \omega\|^2 + \frac{1}{2}\|d_1 \omega\|^2,
$$

whose Euler-Lagrange operator is the Hodge one-form operator

$$
L_1 = d_0 d_0^\ast + d_1^\ast d_1.
$$

The result is not a full foundational proof of HAOS. It is the first concrete derivation showing how harmonic operator structure can arise from recoverability requirements rather than being assumed from the start.

The abstract quadratic step `T1` is now isolated separately in
[HAOS_Quadratic_Defect_Theorem_T1_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/notes/foundations/HAOS_Quadratic_Defect_Theorem_T1_v1.md).
The present note should therefore be read as the first scalar specialization of that abstract theorem.

## 1. Scalar substrate

Let `G = (V, E, w)` be a finite connected undirected weighted graph.

Data:

- vertex set `V = {1, ..., n}`
- edge set `E = {e_1, ..., e_m}`
- positive edge weights `w_e > 0`

Choose an arbitrary orientation on edges and let

$$
B \in \mathbb R^{m \times n}
$$

be the oriented incidence matrix. For a scalar state

$$
f \in \mathbb R^n,
$$

the relational mismatch on edges is

$$
r = Bf \in \mathbb R^m,
\qquad
r_e = f_{h(e)} - f_{t(e)}.
$$

Interpretation:

- scalar states live on vertices
- only edge differences are interaction-visible
- the graph orientation is auxiliary; all final expressions are orientation-independent

## 2. Coherent states and gauge redundancy

`[A]` A scalar coherent state is one for which every visible relational mismatch vanishes:

$$
Bf = 0.
$$

Since the graph is connected, this implies

$$
f = c \mathbf 1
$$

for some constant `c`.

`[D]` Constant shifts are therefore gauge-like redundancies for scalar coherence:

$$
f \sim f + c \mathbf 1.
$$

This is the toy scalar form of relational invariance. Only differences matter.

## 3. Local recoverability potentials

For each edge `e`, choose a local potential

$$
\phi_e : \mathbb R \to \mathbb R_{\ge 0}
$$

subject to the following assumptions.

### S1. Nonnegativity and exact minimum

`[A]`

$$
\phi_e(x) \ge 0,
\qquad
\phi_e(x) = 0 \iff x = 0.
$$

### S2. Reciprocity

`[A]`

$$
\phi_e(-x) = \phi_e(x).
$$

This says recoverability defect does not depend on the chosen orientation sign.

### S3. Smoothness near coherence

`[A]`

Each `phi_e` is `C^3` in a neighborhood of `0`.

### S4. Strict local curvature

`[A]`

$$
\kappa_e := \phi_e''(0) > 0.
$$

This excludes flat directions at the local coherent point.

The simplest examples are:

$$
\phi_e(x) = \frac{1}{2}x^2,
\qquad
\phi_e(x) = 1 - \cos x,
\qquad
\phi_e(x) = \sqrt{1+x^2} - 1.
$$

Each of these is even, nonnegative, minimized at `0`, and has positive second derivative at `0`.

## 4. The scalar recoverability defect

We now define the toy HAOS recoverability defect on scalar relational states:

$$
\mathcal R_0(f)
:=
\sum_{e \in E} w_e \, \phi_e\!\big((Bf)_e\big).
$$

Interpretation:

- recoverability failure is assembled locally on edges
- only relational mismatches enter
- the defect vanishes exactly on coherent scalar states

### Proposition 4.1. Gauge invariance

`[D]`

For every constant `c`,

$$
\mathcal R_0(f + c \mathbf 1) = \mathcal R_0(f).
$$

#### Proof

Since `B \mathbf 1 = 0`,

$$
B(f + c \mathbf 1) = Bf + c B \mathbf 1 = Bf.
$$

Therefore every edge argument of `phi_e` is unchanged, so the sum is unchanged. `\square`

### Proposition 4.2. Coherent-state characterization

`[D]`

If the graph is connected, then

$$
\mathcal R_0(f) = 0
\iff
f = c \mathbf 1
$$

for some constant `c`.

#### Proof

By `S1`, the sum vanishes iff every term vanishes, so

$$
\mathcal R_0(f)=0
\iff
(Bf)_e = 0 \text{ for all } e.
$$

Thus `Bf = 0`. On a connected graph the kernel of the incidence matrix is exactly the constant vectors, so `f = c \mathbf 1`. The converse is immediate because `B(c \mathbf 1)=0`. `\square`

## 5. First quadratic step

We now prove the first real derivation theorem.

### Theorem 5.1. Quadratic leading defect

`[D]`

Let `u` be any fluctuation orthogonal to constants,

$$
u \perp \mathbf 1,
\qquad
f = c \mathbf 1 + u.
$$

Then, for `u` sufficiently small,

$$
\mathcal R_0(c \mathbf 1 + u)
=
\frac{1}{2}\sum_{e \in E} w_e \kappa_e (Bu)_e^2
+
\sum_{e \in E} w_e \rho_e\!\big((Bu)_e\big),
$$

where each remainder satisfies

$$
\rho_e(x) = O(|x|^3)
\quad \text{as } x \to 0.
$$

In particular,

$$
\mathcal R_0(c \mathbf 1 + u)
=
\frac{1}{2}\langle u, L_0 u \rangle
+
O(\|Bu\|_3^3),
$$

with

$$
L_0 = B^\top \Lambda B,
\qquad
\Lambda = \operatorname{diag}(w_e \kappa_e).
$$

#### Proof

Because `phi_e` is even and `C^3` near `0`,

$$
\phi_e(0)=0,
\qquad
\phi_e'(0)=0,
\qquad
\phi_e''(0)=\kappa_e.
$$

Taylor expansion at `0` gives

$$
\phi_e(x)
=
\frac{1}{2}\kappa_e x^2 + \rho_e(x),
\qquad
\rho_e(x)=O(|x|^3).
$$

Now use gauge invariance:

$$
\mathcal R_0(c \mathbf 1 + u)
=
\sum_{e \in E} w_e \phi_e\!\big((Bu)_e\big).
$$

Substituting the Taylor expansion,

$$
\mathcal R_0(c \mathbf 1 + u)
=
\frac{1}{2}\sum_{e \in E} w_e \kappa_e (Bu)_e^2
+
\sum_{e \in E} w_e \rho_e\!\big((Bu)_e\big).
$$

The quadratic term can be written in matrix form:

$$
\frac{1}{2}\sum_{e \in E} w_e \kappa_e (Bu)_e^2
=
\frac{1}{2}(Bu)^\top \Lambda (Bu)
=
\frac{1}{2}u^\top B^\top \Lambda B u.
$$

Thus the leading operator is

$$
L_0 = B^\top \Lambda B.
$$

The remainder bound is inherited from the sum of the local `O(|x|^3)` terms. `\square`

### Corollary 5.2. Dirichlet form

`[D]`

The quadratic leading defect is exactly a weighted Dirichlet form:

$$
\mathcal E_0(u)
:=
\frac{1}{2}\langle u, L_0 u \rangle
=
\frac{1}{2}\sum_{e \in E} w_e \kappa_e (Bu)_e^2.
$$

Equivalently, written over unordered vertex pairs,

$$
\mathcal E_0(u)
=
\frac{1}{2}\sum_{\{i,j\} \in E}
a_{ij}(u_i-u_j)^2
$$

for suitable symmetric coefficients `a_ij > 0`.

Interpretation:

- the first nontrivial HAOS recoverability law is quadratic
- relational invariance turns it into a difference energy
- the scalar harmonic operator is therefore not assumed, but forced by the local recoverability structure of the toy model

### Corollary 5.3. Kernel of the leading operator

`[D]`

If the graph is connected, then

$$
\ker L_0 = \operatorname{span}\{\mathbf 1\}.
$$

#### Proof

If `L_0 u = 0`, then

$$
0 = \langle u, L_0 u \rangle = (Bu)^\top \Lambda (Bu).
$$

Since `Lambda` is positive diagonal, this implies `Bu = 0`. Connectedness then gives `u = c \mathbf 1`. `\square`

## 6. Recovery flow

The Dirichlet form immediately gives the linearized scalar recovery flow.

### Proposition 6.1. Linearized steepest-descent recovery

`[D]`

If scalar recovery is modeled by steepest descent of the quadratic leading defect, then

$$
\partial_t u = -L_0 u.
$$

#### Proof

The first variation of

$$
\mathcal E_0(u)=\frac{1}{2}\langle u, L_0 u\rangle
$$

is

$$
\delta \mathcal E_0(u)[v] = \langle L_0 u, v \rangle.
$$

Therefore the gradient of the defect is `L0 u`, and steepest descent gives

$$
\partial_t u = - \nabla \mathcal E_0(u) = -L_0 u.
$$

`\square`

Interpretation:

- low modes are slow recoverability modes
- constants are neutral coherent modes
- scalar Laplacian-type recovery appears before any Hodge lift is introduced

## 7. Only now: lift to graded incidence

The scalar step above is the real guardrail. Only after it is complete do we lift.

### 7.1 Graded data

Introduce a finite graded complex

$$
V_0 \xrightarrow{d_0} V_1 \xrightarrow{d_1} V_2
$$

with inner products on each grade and the compatibility identity

$$
d_1 d_0 = 0.
$$

Interpretation:

- `V0`: scalar states
- `V1`: relational circulation-like states
- `V2`: face / plaquette consistency data

This is the minimal graded setting needed for one-form Hodge structure.

### 7.2 One-form recoverability defect

`[A]` The minimal one-form recoverability defect is defined by penalizing:

1. downward inconsistency through `d_0^\ast \omega`
2. upward inconsistency through `d_1 \omega`

Thus define

$$
\mathcal R_1(\omega)
:=
\frac{1}{2}\|d_0^\ast \omega\|^2
+
\frac{1}{2}\|d_1 \omega\|^2.
$$

Interpretation:

- the first term removes exact leakage / unresolved source-like content
- the second term penalizes failure of circulation closure

## 8. Hodge lift theorem

### Theorem 8.1. One-form Euler-Lagrange operator

`[D]`

The first variation of `R1` is

$$
\delta \mathcal R_1(\omega)[\eta]
=
\langle (d_0 d_0^\ast + d_1^\ast d_1)\omega, \eta \rangle.
$$

Hence the one-form recoverability operator is

$$
L_1 = d_0 d_0^\ast + d_1^\ast d_1.
$$

#### Proof

Using the adjoint identities,

$$
\delta \mathcal R_1(\omega)[\eta]
=
\langle d_0^\ast \omega, d_0^\ast \eta \rangle
+
\langle d_1 \omega, d_1 \eta \rangle
$$

$$
=
\langle d_0 d_0^\ast \omega, \eta \rangle
+
\langle d_1^\ast d_1 \omega, \eta \rangle
$$

$$
=
\langle (d_0 d_0^\ast + d_1^\ast d_1)\omega, \eta \rangle.
$$

Since this holds for every test direction `eta`, the Euler-Lagrange operator is exactly `L1`. `\square`

## 9. Exact / harmonic / coexact split

### Proposition 9.1. Orthogonality of exact and coexact sectors

`[D]`

The subspaces `im d_0` and `im d_1^\ast` are orthogonal.

#### Proof

For any `u in V0` and `beta in V2`,

$$
\langle d_0 u, d_1^\ast \beta \rangle
=
\langle d_1 d_0 u, \beta \rangle
=
0
$$

because `d_1 d_0 = 0`. `\square`

### Proposition 9.2. Harmonic complement

`[D]`

Define the harmonic subspace

$$
\mathcal H^1 := \ker d_0^\ast \cap \ker d_1.
$$

Then

$$
V_1 = \operatorname{im} d_0 \oplus \mathcal H^1 \oplus \operatorname{im} d_1^\ast
$$

as an orthogonal direct sum in finite dimension.

#### Proof sketch

By Proposition 9.1, `im d_0` and `im d_1^\ast` are orthogonal. Their orthogonal complement consists of vectors `omega` satisfying

$$
\langle \omega, d_0 u \rangle = 0
\quad \forall u,
\qquad
\langle \omega, d_1^\ast \beta \rangle = 0
\quad \forall \beta.
$$

Equivalently,

$$
d_0^\ast \omega = 0,
\qquad
d_1 \omega = 0.
$$

So the orthogonal complement is exactly `H^1`. Finite-dimensional orthogonal decomposition then gives the result. `\square`

Interpretation:

- exact = longitudinal / redundancy sector
- harmonic = topological memory / neutral sector
- coexact = transverse circulation sector

## 10. Relation to the current HAOS-IIP operators

`[E]` The repository already uses the operator forms

$$
L_0 = D - A
$$

for the scalar node sector and

$$
L_1 = d_0 d_0^\ast + d_1^\ast d_1
$$

for the one-form edge sector, together with the exact / harmonic / coexact split and the restricted transverse operator.

`[D]` What this note contributes is not a new numerical result. It gives a toy derivation path explaining how those operator forms can arise from recoverability logic rather than being introduced as independent mathematical machinery.

## 11. What is actually proved here

`[D]` Under the toy assumptions of Sections 1-4, this note proves:

1. scalar recoverability defect is gauge-invariant under constant shifts
2. coherent scalar states are exactly the constant states on a connected substrate
3. the leading defect near coherence is quadratic
4. the quadratic leading defect is a weighted Dirichlet form
5. the induced scalar recovery operator is Laplacian-type
6. the first graded lift of the same logic yields the one-form Hodge operator
7. exact / harmonic / coexact decomposition follows from the graded incidence structure

## 12. What remains open

`[A]`

This is still a toy derivation, not yet a full HAOS proof, because several ingredients are assumed rather than derived from first principles:

- local potentials `phi_e`
- smoothness near coherence
- the graded incidence complex itself
- the one-form defect ansatz

The next real foundational step would be:

> derive why HAOS forces these assumptions, rather than merely showing that once they are accepted, harmonic operators follow.

## 13. Shortest honest conclusion

`[D]`

The first scalar step now closes cleanly:

> a HAOS-aligned recoverability defect on relational scalar states reduces near coherence to a weighted Dirichlet form, and its recovery generator is a Laplacian-type operator.

Only after that scalar step is fixed does the graded lift produce the Hodge operator and exact / harmonic / coexact sector structure.

That is the correct order of derivation.
