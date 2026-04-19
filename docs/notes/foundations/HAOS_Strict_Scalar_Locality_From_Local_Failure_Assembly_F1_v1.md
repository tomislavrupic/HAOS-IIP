# HAOS_Strict_Scalar_Locality_From_Local_Failure_Assembly_F1_v1

Status:

- structural bridge note
- purpose: execute `F1` in bounded form by tightening the scalar locality gap left open in `T2`

Status labels:

- `[D]` derived in this note under the stated assumptions
- `[A]` assumption introduced for the theorem
- `[O]` open residue that remains after the present argument

## 1. Statement

The scalar `T2` note currently assumes strict pairwise locality:

```text
the leading scalar primitive defect terms are supported on single interacting pairs only
```

This note tightens that step in the strongest form currently justified by HAOS:

> if scalar local-failure assembly is evaluated on a pair-primitive scalar audit substrate, then the leading quadratic scalar defect necessarily decomposes over interacting pairs.

In formulas,

$$
D_0(f)
\sim
\frac{1}{2}\sum_{\{i,j\}\in E} q_{ij}(f_i,f_j),
$$

where each primitive term is supported on one interacting pair.

This does not yet prove that every conceivable HAOS scalar substrate must be pair-primitive.
It proves that once the scalar audit layer is pair-primitive, strict scalar pair locality is no longer an independent extra assumption.

## 2. What is being closed

`T2` already used:

- `T1`: the leading defect is quadratic
- `H1`: relational invariance
- `H6`: reciprocity

But its Step 2 introduced a standalone assumption:

```text
strict pairwise scalar locality
```

The present note reduces that gap by deriving the same pairwise support from:

1. `H3`: global coherence failure is assembled from local mismatch primitives
2. a scalar audit-substrate assumption saying that the primitive local mismatch cells in the scalar sector are the directly interacting pairs

That second input is the smallest structural choice needed to stop `H3` from being compatible with irreducible triangle-like or larger scalar primitives.

## 3. Setup

Let

$$
V_0 = \mathbb R^n
$$

be the scalar node-state space on a finite interaction substrate with interaction relation `E`.

From `T1`,

$$
D_0(f)
=
\frac{1}{2}\langle f,Qf\rangle + o(\|f\|^2),
$$

with `Q` symmetric and positive semidefinite.

### Assumption 3.1. Pair-primitive scalar audit substrate

`[A]`

In the scalar sector, the primitive local audit cells are exactly the directly interacting pairs:

$$
\mathcal A_0 = \bigl\{\{i,j\} : \{i,j\}\in E\bigr\}.
$$

Meaning:

- there are no independent scalar primitive audit cells supported on non-interacting pairs
- there are no irreducible scalar primitive audit cells supported on triples or larger sets

This is weaker than assuming the final Dirichlet form.
It only specifies what counts as a primitive scalar locality cell before the quadratic form is identified.

### Assumption 3.2. Local-failure assembly on primitive cells

`[A]`

`H3` is taken in the scalar sector to mean:

- the leading scalar defect is assembled from primitive local mismatch contributions
- each primitive contribution is supported on one scalar audit cell in `\mathcal A_0`

Since `T1` already fixed the leading defect to be quadratic, every primitive leading contribution is a quadratic form on its supporting audit cell.

## 4. Primitive-support theorem

### Theorem 4.1. Pair support from local-failure assembly

`[D]`

Under Assumptions 3.1-3.2, the leading quadratic scalar defect has the form

$$
\langle f,Qf\rangle
=
\sum_{\{i,j\}\in E} q_{ij}(f_i,f_j),
$$

where each `q_{ij}` is a quadratic form supported on one interacting pair.

#### Proof

By `T1`, the leading scalar defect is quadratic. By Assumption 3.2, that leading defect is assembled from primitive local mismatch terms, each supported on one scalar audit cell. By Assumption 3.1, the scalar audit cells are exactly the interacting pairs `{i,j}` in `E`. Therefore each primitive quadratic contribution depends only on the two coordinates `(f_i,f_j)`. Summing those primitive contributions over all scalar audit cells yields

$$
\langle f,Qf\rangle
=
\sum_{\{i,j\}\in E} q_{ij}(f_i,f_j).
$$

`\square`

## 5. Relation to the scalar Dirichlet theorem

Theorem 4.1 is exactly the Step 2 input needed by `T2`.

Once pair support is in place, the remaining `T2` steps still follow:

1. `H1` forces constants into `\ker Q`
2. `H6` forces reciprocal pair terms
3. shift invariance on each pair forces each local quadratic to reduce to

$$
q_{ij}(x,y)=w_{ij}(x-y)^2
$$

with `w_{ij}=w_{ji}\ge 0`

So the scalar Dirichlet form becomes

$$
D_0(f)
\sim
\frac{1}{2}\sum_{\{i,j\}\in E} w_{ij}(f_i-f_j)^2.
$$

## 6. What this note does and does not prove

### 6.1 What is now derived

`[D]`

The old `T2` strict pairwise-locality step is no longer free-standing once the scalar audit substrate is specified as pair-primitive.

That is a real tightening because the theorem now separates:

- primitive scalar locality structure
- the later difference-form identification

instead of assuming both at once.

### 6.2 What remains open

`[O]`

The note does not prove that `H3` alone forces pair-primitive scalar audit cells in every possible HAOS substrate.

If HAOS were allowed to contain irreducible scalar primitive cells on triples or larger finite clusters, then `H3` would still permit a different leading scalar locality law.

So the remaining open residue is:

```text
does HAOS local-failure assembly itself force pair-primitive scalar audit cells,
or is pair-primitivity an additional structural choice for the scalar sector?
```

That is the honest remainder of `F1`.

## 7. Bounded conclusion

`[D]`

The strongest justified `F1` statement is:

> on pair-primitive scalar audit substrates, HAOS local-failure assembly forces the leading scalar defect to assemble over interacting pairs, so the pair-support step in `T2` is derived rather than separately assumed.

This is enough to tighten the scalar derivation program now, while still keeping the final substrate-independence question open.
