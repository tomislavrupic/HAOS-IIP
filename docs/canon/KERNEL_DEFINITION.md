# KERNEL_DEFINITION

## Baseline kernel

`[E]` The primary kernel in the current HAOS-IIP program is

$$
A_{ij} = \exp\!\left(-\frac{|x_i-x_j|^2}{2\epsilon_k}\right), \qquad A_{ii}=0.
$$

Interpretation:

- `x_i`: node coordinate or embedded feature vector
- `epsilon_k`: kernel width
- `A_ij`: interaction weight between nodes `i` and `j`

## Derived matrices

`[E]`

$$
D_{ii} = \sum_j A_{ij}, \qquad L = D - A.
$$

Optional normalized forms may be used for diagnostics, but the ordinary Laplacian remains primary unless stated otherwise.

## Recoverability link

`[E]` Low-spectrum behavior is tied to perturbation recovery:

$$
x(t)=e^{-Lt}x(0), \qquad e^{-\lambda_1 T}\delta < \epsilon_c.
$$

The symbol `epsilon_c` is reserved for coherence threshold, not kernel width.

## Practical conventions

`[E]` Current numerics choose `epsilon_k` from local spacing statistics, typically a multiple of mean nearest-neighbor squared distance.

`[P]` Kernel truncation by finite radius is allowed for numerics if the truncation scale is reported explicitly.

## Phase-dressed extension

`[P]` The minimal gauge-oriented extension is

$$
A_{ij}^{(\theta)} = A_{ij} e^{i\theta_{ij}}.
$$

This is not part of the baseline kernel definition. It belongs to the gauge program.
