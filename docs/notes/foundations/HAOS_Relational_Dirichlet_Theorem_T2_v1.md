# HAOS_Relational_Dirichlet_Theorem_T2_v1

Status:

- abstract structural derivation note
- purpose: execute `T2` of the HAOS-to-harmonic derivation program by specializing the abstract quadratic defect from `T1` to node-like relational states and deriving the scalar Dirichlet form

Status labels:

- `[D]` derived in this note under the stated assumptions
- `[A]` assumption introduced for the theorem

## 1. Statement

We prove the `T2` statement in the form needed by the derivation program:

> on node-like states `f`, the leading quadratic defect reduces to a weighted difference energy
>
> $$
> D_0(f) \sim \frac{1}{2}\sum_{\{i,j\}\in E} w_{ij}(f_i-f_j)^2,
> $$
>
> where the sum runs only over directly interacting pairs.

By representation this equals

$$
D_0(f) = \frac{1}{2}\langle f, L_0 f \rangle
$$

with `L0` a symmetric positive-semidefinite Laplacian-type operator on the discrete interaction substrate.

No metric, embedding, continuum limit, or geometry is assumed. Only the relational interaction pattern and the HAOS axioms `H1`, `H3`, `H6`, together with the quadratic-defect result of `T1`, are used.

## 2. Setup

Let the node-like state space be

$$
V_0 = \mathbb R^n.
$$

The interaction substrate provides only:

- a finite set of relational elements `1, ..., n`
- a symmetric interaction relation `E`

where

$$
\{i,j\}\in E
$$

means that `i` and `j` interact directly.

Input from `T1`:

$$
D_0(f) = \frac{1}{2}\langle f, Q f \rangle + o(\|f\|^2),
$$

where `Q` is symmetric and positive semidefinite.

The task is to identify the structure of `Q`.

## 3. Step 1: relational invariance forces constants into the kernel

### Proposition 3.1. Translation invariance

`[A]`

`H1` says that only interaction differences matter. Therefore the node-like defect must be invariant under uniform shifts:

$$
D_0(f + c\mathbf 1) = D_0(f)
\quad \forall c \in \mathbb R.
$$

### Corollary 3.2. Constants lie in `ker Q`

`[D]`

The quadratic generator satisfies

$$
Q\mathbf 1 = 0.
$$

#### Proof

Apply the quadratic leading term from `T1`:

$$
\langle f + c\mathbf 1,\, Q(f + c\mathbf 1)\rangle
=
\langle f,Qf\rangle
+
2c\langle f,Q\mathbf 1\rangle
+
c^2\langle \mathbf 1,Q\mathbf 1\rangle.
$$

Since the defect is invariant for every `c`, the coefficients of `c` and `c^2` must vanish. Hence

$$
Q\mathbf 1 = 0,
\qquad
\langle \mathbf 1,Q\mathbf 1\rangle = 0.
$$

`\square`

Interpretation:

- constants are invisible to relational coherence
- the quadratic law can only depend on differences

## 4. Step 2: local failure assembly forces pairwise support

`H3` says global coherence failure is assembled from local interaction mismatches.

### Assumption 4.1. Strict pairwise locality

`[A]`

In the node-like scalar sector, the leading local contributions are supported on single interacting pairs only.

Equivalently, no primitive quadratic term is allowed to depend simultaneously on a non-interacting triple or on a nonlocal pair.

This is the scalar form of local failure assembly.

The bounded `F1` bridge note
[HAOS_Strict_Scalar_Locality_From_Local_Failure_Assembly_F1_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/notes/foundations/HAOS_Strict_Scalar_Locality_From_Local_Failure_Assembly_F1_v1.md)
shows that this step is derived once the scalar audit substrate is specified as pair-primitive.

### Proposition 4.2. Pairwise decomposition

`[D]`

Under strict pairwise locality, the quadratic form decomposes as

$$
\langle f,Qf\rangle
=
\sum_{\{i,j\}\in E} q_{ij}(f_i,f_j),
$$

where each `q_ij` is a quadratic form on two scalar variables.

#### Proof

Because the defect is quadratic by `T1` and primitive contributions are pairwise local by Assumption 4.1, every leading local term must be supported on one interacting pair and no other coordinates. Summing over all interacting pairs gives the decomposition. `\square`

## 5. Step 3: reciprocity forces difference squares

Now fix one interacting pair `{i,j}`.

The most general quadratic form on two scalars is

$$
q_{ij}(x,y)=a_{ij}x^2+2b_{ij}xy+c_{ij}y^2.
$$

`H6` imposes reciprocity:

$$
q_{ij}(x,y)=q_{ij}(y,x).
$$

So

$$
a_{ij}=c_{ij}.
$$

Translation invariance from Step 1 imposes

$$
q_{ij}(x+c,y+c)=q_{ij}(x,y)
\quad \forall c.
$$

Expanding in `c`, the coefficient of `c^2` must vanish:

$$
a_{ij}+2b_{ij}+c_{ij}=0.
$$

Since `a_{ij}=c_{ij}`, this becomes

$$
b_{ij}=-a_{ij}.
$$

Therefore

$$
q_{ij}(x,y)
=
a_{ij}(x^2-2xy+y^2)
=
a_{ij}(x-y)^2.
$$

Define

$$
w_{ij}:=a_{ij}=w_{ji}\ge 0.
$$

Nonnegativity follows from positive semidefiniteness of `Q`.

### Theorem 5.1. Relational Dirichlet form

`[D]`

The leading node-like defect has the form

$$
D_0(f)
\sim
\frac{1}{2}\sum_{\{i,j\}\in E} w_{ij}(f_i-f_j)^2.
$$

This is exactly a weighted Dirichlet energy on the interaction graph.

## 6. Step 4: operator representation

Define the weighted Laplacian-type operator

$$
(L_0 f)_i
:=
\sum_{j:\{i,j\}\in E} w_{ij}(f_i-f_j).
$$

Then a direct expansion gives

$$
\langle f,L_0 f\rangle
=
\sum_{\{i,j\}\in E} w_{ij}(f_i-f_j)^2.
$$

Hence

$$
D_0(f)
\sim
\frac{1}{2}\langle f,L_0 f\rangle.
$$

Properties:

- `L0` is symmetric
- `L0` is positive semidefinite
- `L0 1 = 0`

So `L0` is exactly the minimal semibounded scalar recoverability generator.

## 7. Why this is the minimal relational scalar operator

Each axiom removes a class of alternatives.

- `H1` removes absolute-value dependence and forces constant-shift invariance
- `H3` removes nonlocal primitive couplings
- `H6` removes directed or antisymmetric local responses

What remains is only:

$$
\text{weighted pairwise difference squares}
\quad \Longrightarrow \quad
\text{Dirichlet form}
\quad \Longrightarrow \quad
\text{graph Laplacian}.
$$

So the scalar Laplacian is not inserted as geometry. It is the unique quadratic operator shape compatible with:

1. relational invariance
2. pairwise local failure assembly
3. symmetric reciprocity

## 8. Relation to the existing notes

- the abstract quadratic step is in [HAOS_Quadratic_Defect_Theorem_T1_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/notes/foundations/HAOS_Quadratic_Defect_Theorem_T1_v1.md)
- the concrete potential-based realization and first graded lift are in [HAOS_Relational_Scalar_Recoverability_Defect_and_First_Dirichlet_Theorem_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/notes/foundations/HAOS_Relational_Scalar_Recoverability_Defect_and_First_Dirichlet_Theorem_v1.md)

## 9. Shortest honest conclusion

`[D]`

The `T2` step now closes in structural form:

> once the recoverability defect is quadratic, translation-invariant, pairwise local, and reciprocal, the node-like sector is forced into a weighted Dirichlet form and therefore into a Laplacian-type operator.

That is the first strictly relational scalar operator theorem in the HAOS derivation stack.
