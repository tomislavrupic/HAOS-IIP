# EL-R3-THEORY-REVISION-01

Status: `FROZEN_THEORY_REVISION`

Date frozen: `2026-08-18`

## Decision

The next Rung 3 candidate must not try harder with one-bit local relational
identity. It must test a genuinely different information architecture:

> Functional recovery may require bounded, spatially distributed mesoscopic
> memory that retains amplitude and phase information while remaining
> provably incomplete with respect to the microscopic state.

The representation is a uniformly quantized field on a fixed `4 x 4`
partition of the physical domain. Refinement increases microscopic resolution
inside each physical cell; it does not increase the number of stored values.

| Microscopic field | Microscopic values per coarse cell | Stored field |
| --- | ---: | --- |
| `8 x 8` | `2 x 2 = 4` | `4 x 4 = 16` values |
| `12 x 12` | `3 x 3 = 9` | `4 x 4 = 16` values |
| `16 x 16` | `4 x 4 = 16` | `4 x 4 = 16` values |

This is a mesoscopic checkpoint. It is not described as non-checkpoint
memory. The protection against microscopic leakage is its fixed budget,
target-blind construction, and provable nullspace.

## Representation

Let `x` be an `n x n` scalar field with `n` divisible by four. Let `C_n` have
sixteen rows, one for each nonempty, disjoint physical coarse cell. Each row
averages the `m = (n/4)^2` microscopic values in its cell.

The stored payload is

```text
q_b = Q_b(C_n x),
```

where `Q_b` is a completely frozen scalar quantizer applied componentwise. No
per-sample scale, offset, normalization, codebook, or adaptive range is
available to the decoder.

The direct four functional phasors are a prohibited oracle because storing
them stores the evaluation target. The coarse field is constructed without
using the function, outcome labels, final trajectories, or recovery scores.

## Rank And Information Loss

The sixteen averaging rows have disjoint nonempty support. They are therefore
linearly independent:

```text
rank(C_n) = 16,
nullity(C_n) = n^2 - 16.
```

For equal-size cells,

```text
C_n C_n^T = (1/m) I_16.
```

The unquantized coarse map already leaves `n^2 - 16` microscopic degrees of
freedom unresolved. Quantization makes the representation more many-to-one;
it cannot improve microscopic reconstruction.

## Deterministic Decoder

Given disrupted state `x'`, restore the stored quantized field with

```text
x_hat = x' + C_n^T (C_n C_n^T)^-1 (q_b - C_n x').
```

Because the cells are equal and disjoint,

```text
C_n^T (C_n C_n^T)^-1 = m C_n^T.
```

For microscopic index `i` in coarse cell `S_j`, the update is simply

```text
x_hat_i = x'_i + q_b,j - mean(x'_k for k in S_j).
```

It follows exactly, up to declared floating-point tolerance, that

```text
C_n x_hat = q_b.
```

There is no optimizer, training, gain, iteration count, fitted parameter, or
full-state reference hidden in the decoder.

## Predicted Failure Channels

Let

```text
P_n = C_n^T (C_n C_n^T)^-1 C_n.
```

Then

```text
x_hat - x
  = (I - P_n)(x' - x)
  + C_n^T (C_n C_n^T)^-1 (q_b - C_n x).
```

The candidate has two predeclared failure channels:

1. disruption inside the unresolved within-cell subspace;
2. quantization or saturation error in the stored field.

For the frozen linear function, functional error inherits the same
decomposition. This supplies a natural failure basin and makes a negative
result mechanistically interpretable.

## Information Budget

The primary candidate stores sixteen codes of exactly `b` bits each:

```text
payload = 16 b bits.
```

The canonical row-major `4 x 4` ordering requires no per-sample addresses.
Partition geometry, quantizer definition, decoder rule, and code version are
constant schema overhead and must be hashed and reported separately. Any
state-dependent scale, offset, lookup table, missingness side channel, or
adaptive metadata would be additional payload and invalidates the declared
budget.

## Governance Closure

This revision satisfies the prior theory gate because it:

- introduces new amplitude-bearing mesoscopic information rather than
  re-encoding RT-02 orientation or RP-01 parity;
- fixes a bounded representation that cannot reconstruct the continuous
  microscopic state;
- supplies a deterministic restorative rule;
- predicts an explicit failure basin;
- defines destructive equal-budget controls before execution;
- preserves RT-01, RT-02, and RP-01 unchanged;
- makes no recovery claim by itself.

The authorized successor candidate is
`EL-R3-MESOSCOPIC-FIELD-01`. Its separate frozen precommitment governs all
execution and claim language.
