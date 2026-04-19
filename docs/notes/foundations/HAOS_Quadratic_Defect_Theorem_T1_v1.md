# HAOS_Quadratic_Defect_Theorem_T1_v1

Status:

- abstract derivation note
- purpose: execute `T1` of the HAOS-to-harmonic derivation program in an abstract state-space setting before any graph, Dirichlet, or Hodge specialization is introduced

Status labels:

- `[D]` derived in this note under the stated assumptions
- `[A]` assumption introduced for the theorem

## 1. Statement

This note proves the first theorem in the HAOS-to-harmonic derivation ladder.

Target statement:

> near a coherent state `u = 0`, the recoverability defect admits the expansion
>
> $$
> D(u) = \frac{1}{2}\langle u, Q u \rangle + o(\|u\|^2)
> $$
>
> where `Q` is a positive-semidefinite self-adjoint bounded linear operator on the state space.

This is the first genuinely harmonic step. It shows that the first nontrivial recoverability law is quadratic before any scalar Dirichlet form, graph Laplacian, incidence complex, or Hodge decomposition is introduced.

## 2. Setup

Let `V` be a real Banach space of states attached to the relational degrees of freedom of an interaction substrate.

Data:

- state space `V`
- norm `||.||`
- recoverability defect

$$
D : V \to [0, +\infty)
$$

such that the coherent state is represented by the origin:

$$
D(0) = 0.
$$

Interpretation:

- `D(u)` measures failure of recoverable coherence
- `u = 0` is the fully coherent reference state
- nonnegativity means recoverability failure is never negative

## 3. Assumptions used

The theorem uses the following HAOS-aligned assumptions.

### H2. Recoverability principle

`[A]`

The coherent state `u = 0` is not just any point. It is a state of recoverable coherence, and the defect functional is nonnegative everywhere:

$$
D(u) \ge 0
\quad \forall u \in V.
$$

### H4. Weak compositionality

`[A]`

For weakly coupled subsystems, recoverability defect is additive to leading order.

In the present theorem, this assumption is not needed to force the existence of the quadratic term itself. Its role is interpretive: it guarantees that the leading quadratic form is the correct leading-order coherence law rather than an arbitrary global artifact.

### H5. Smoothness near coherent states

`[A]`

The defect is twice Fr\'echet differentiable at `0`.

Equivalently, there exist:

- a continuous linear functional `ell in V*`
- a continuous symmetric bilinear form `B : V x V -> R`
- a remainder `r(u)`

such that

$$
D(u) = D(0) + \ell(u) + \frac{1}{2} B(u,u) + r(u),
$$

with

$$
\lim_{\|u\| \to 0} \frac{r(u)}{\|u\|^2} = 0.
$$

Since `D(0)=0`, this becomes

$$
D(u) = \ell(u) + \frac{1}{2} B(u,u) + r(u).
$$

## 4. Quadratic defect theorem

### Theorem 4.1. Quadratic leading defect

`[D]`

Under `H2`, `H4`, and `H5`, the recoverability defect satisfies

$$
D(u) = \frac{1}{2}\langle u, Q u \rangle + o(\|u\|^2)
$$

near `u = 0`, where `Q` is a bounded self-adjoint positive-semidefinite operator.

#### Proof

By `H5`, we already have the second-order expansion

$$
D(u) = \ell(u) + \frac{1}{2} B(u,u) + r(u),
$$

with

$$
\frac{r(u)}{\|u\|^2} \to 0
\quad \text{as } \|u\| \to 0.
$$

So the only missing steps are:

1. show the linear term vanishes
2. identify the quadratic operator
3. prove positivity

### Step 1. Vanishing of the linear term

Assume, for contradiction, that `ell != 0`.

Then there exists `v in V` with `||v|| = 1` such that

$$
\ell(v) < 0.
$$

For `t > 0` small, evaluate the defect on `u = tv`:

$$
D(tv) = t \ell(v) + \frac{t^2}{2} B(v,v) + r(tv).
$$

Divide by `t > 0`:

$$
\frac{D(tv)}{t}
=
\ell(v) + \frac{t}{2} B(v,v) + \frac{r(tv)}{t}.
$$

Because

$$
\frac{r(tv)}{t}
=
t \cdot \frac{r(tv)}{t^2}
\to 0
\quad \text{as } t \to 0^+,
$$

the last two terms vanish in the limit. Therefore for sufficiently small `t > 0`,

$$
D(tv) < 0.
$$

This contradicts nonnegativity of the recoverability defect.

Hence

$$
\ell \equiv 0.
$$

So the expansion reduces to

$$
D(u) = \frac{1}{2} B(u,u) + r(u)
=
\frac{1}{2} B(u,u) + o(\|u\|^2).
$$

### Step 2. Identification of the quadratic operator

If `V` is a Hilbert space, the Riesz representation theorem gives a unique bounded self-adjoint operator `Q` such that

$$
B(u,v) = \langle u, Q v \rangle
\quad \forall u,v \in V.
$$

In particular,

$$
B(u,u) = \langle u, Q u \rangle.
$$

Therefore

$$
D(u) = \frac{1}{2}\langle u, Q u \rangle + o(\|u\|^2).
$$

If one works only on a Banach space, then `B` is still a continuous symmetric bilinear form, and the theorem should be read at the bilinear-form level unless a specific identification of `V` with its dual has been chosen.

### Step 3. Positive semidefiniteness

For arbitrary `u in V` and `t > 0`,

$$
D(tu) \ge 0.
$$

Divide by `t^2/2`:

$$
\frac{2D(tu)}{t^2} \ge 0.
$$

Now use the expansion:

$$
\frac{2D(tu)}{t^2}
=
\frac{B(tu,tu)}{t^2}
+
\frac{2r(tu)}{t^2}
=
B(u,u) + \frac{2r(tu)}{t^2}.
$$

Taking `t -> 0`,

$$
B(u,u) \ge 0.
$$

Hence, in Hilbert-space form,

$$
\langle u, Q u \rangle \ge 0
\quad \forall u.
$$

So `Q` is positive semidefinite.

This completes the proof. `\square`

## 5. Role of H4

`[D]`

The abstract Taylor argument uses only nonnegativity and second differentiability to force:

- vanishing of the linear term
- existence of a quadratic leading defect
- positivity of the quadratic form

`H4` enters at the next interpretive layer.

Without `H4`, the quadratic form could still exist but be an arbitrary global bilinear object with no reason to respect the decomposition of the interaction substrate.

With `H4`, the leading quadratic law is required to assemble additively across weakly coupled subsystems. This is what later makes the quadratic form localizable and differential-like, and what allows the scalar specialization to become a Dirichlet form rather than a generic positive quadratic energy.

So the role of `H4` is:

```text
T1:
quadratic law exists and is positive

T2:
quadratic law localizes into a relational difference energy
```

## 6. What is and is not proved here

### Proved

`[D]`

Under the stated assumptions:

1. the defect has no linear term near coherence
2. the first nontrivial coherence law is quadratic
3. the quadratic law is positive semidefinite
4. in Hilbert-space form, the law is generated by a bounded self-adjoint operator `Q`

### Not yet proved

`[A]`

This theorem does not yet force:

- locality
- graph structure
- Dirichlet-form structure
- incidence maps
- Hodge operators

Those belong to later steps.

## 7. Exact place in the derivation ladder

This note closes only `T1`.

The next steps are:

1. specialize the quadratic form to relational scalar states and prove the Dirichlet reduction
2. only after that, introduce graded incidence
3. then derive the one-form Hodge operator and the exact / harmonic / coexact split

That is the correct derivation order.

## 8. Relation to the existing notes

- the abstract derivation ladder is in [HAOS_to_Harmonic_Operator_Derivation_Program_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/notes/foundations/HAOS_to_Harmonic_Operator_Derivation_Program_v1.md)
- the scalar specialization and first graded lift are in [HAOS_Relational_Scalar_Recoverability_Defect_and_First_Dirichlet_Theorem_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/notes/foundations/HAOS_Relational_Scalar_Recoverability_Defect_and_First_Dirichlet_Theorem_v1.md)

## 9. Shortest honest conclusion

`[D]`

The first theorem now closes cleanly:

> near a coherent state, a HAOS-aligned recoverability defect has a positive quadratic leading term.

That is the first real bridge from HAOS coherence language to harmonic operator language.
