# Kernel Graph Green Response v1

## Purpose

Test the remaining open assumption in the current operator chain:

$$
\text{interaction kernel}
\rightarrow
\text{graph Laplacian}
\rightarrow
\text{Green field } \sim \frac{1}{r}.
$$

The question is not whether a continuum Laplacian produces `1/r`; that is already known. The question is whether the actual interaction kernel used in the repo induces a graph Laplacian whose discrete Green response approaches

$$
\phi(r) \approx A + \frac{B}{r}.
$$

## Construction

### Substrate

Use a regular cubic 3D point cloud on the unit cube.

- cubic sides tested: `n = 9, 11, 13`
- total nodes: `729, 1331, 2197`

The source node is the central lattice point.

### Kernel graph

For embedded points `x_i`, define the weighted adjacency matrix

$$
W_{ij}
=
\exp\left(-\frac{|x_i-x_j|^2}{2\epsilon_k}\right)
$$

on a cutoff graph, with cutoff radius

$$
r_c = 2.5 \sqrt{\epsilon_k}.
$$

The kernel width is scaled with the lattice spacing `h = 1/(n-1)`:

$$
\epsilon_k = c_\epsilon h^2,
\qquad
c_\epsilon \in \{0.5, 1.0\}.
$$

This keeps the graph local as the lattice is refined.

### Graph Laplacian

Construct

$$
L = D - W,
$$

with `D_ii = sum_j W_ij`.

### Mean-zero Poisson problem

Because `L 1 = 0`, solve in the mean-zero subspace using a neutralized source:

$$
s = \delta_{\text{source}} - \frac{1}{N}\mathbf{1}.
$$

Then solve

$$
L \phi = s
$$

via

$$
(L + J)\phi = s,
\qquad
J = \frac{1}{N}\mathbf{1}\mathbf{1}^T,
$$

which fixes the constant-mode ambiguity and returns a mean-zero field.

## Measurement

For each run:

1. compute Euclidean distance from the source node
2. radially average the field on distance shells
3. fit the node-level response against

$$
\phi(r) \approx A + \frac{B}{r}
$$

over a clean midrange window

$$
r \in [\max(2h, 1.5\sqrt{\epsilon_k}), 0.55].
$$

Also estimate the shell-averaged gradient slope to test the inverse-square branch, with the caveat that discrete shell derivatives are noisier than the field fit itself.

## Artifacts

- script: `numerics/simulations/kernel_graph_green_response.py`
- results: `data/20260310_110539_kernel_graph_green_response.json`
- latest: `data/kernel_graph_green_response_latest.json`
- plots:
  - `plots/20260310_110539_kernel_graph_green_profiles.png`
  - `plots/20260310_110539_kernel_graph_green_exponents.png`
  - `plots/20260310_110539_kernel_graph_green_fit_quality.png`

## Results

Field-fit diagnostics:

| `n` | `c_epsilon` | slope of `|phi-A|` vs `r` | fit `R^2` |
| --- | --- | --- | --- |
| 9 | 0.5 | `-0.9830` | `0.9834` |
| 11 | 0.5 | `-0.9703` | `0.9856` |
| 13 | 0.5 | `-0.9687` | `0.9892` |
| 9 | 1.0 | `-0.9780` | `0.9801` |
| 11 | 1.0 | `-0.9663` | `0.9829` |
| 13 | 1.0 | `-0.9645` | `0.9870` |

These are the key numbers. Across the entire scan,

- the field exponent stays close to `-1`
- the `A + B/r` fit quality stays above `0.98`

So the primary Green-response test is positive.

Direct shell-derivative estimates are noisier:

| `n` | `c_epsilon` | slope of `|dphi/dr|` vs `r` |
| --- | --- | --- |
| 9 | 0.5 | `-2.5533` |
| 11 | 0.5 | `-2.8068` |
| 13 | 0.5 | `-2.6982` |
| 9 | 1.0 | `-2.3387` |
| 11 | 1.0 | `-2.6744` |
| 13 | 1.0 | `-2.6216` |

These derivative slopes decay faster than `r^-2`, but they are visibly less stable than the field fit. On the finite graph, the derivative estimate is more sensitive to shell irregularity than the field itself.

## Interpretation

The correct reading is:

> the weighted interaction kernel induces a graph Laplacian whose Green response is consistent with an inverse-distance field on the cubic 3D substrate.

That closes the main structural gap more cleanly than the raw force estimate does.

The graph experiment therefore supports the chain

$$
\text{interaction kernel}
\rightarrow
\text{graph Laplacian}
\rightarrow
\phi(r) \approx A + \frac{B}{r}.
$$

The derivative branch is directionally consistent with inverse-square decay, but the present shell-based measurement is not yet clean enough to use as the primary diagnostic.

## Current Verdict

- yes, the kernel-induced graph Laplacian produces a near-`A + B/r` Green response on the cubic scan
- yes, this is the right closure test for the current operator program
- the direct inverse-square derivative estimate is still noisier than the field fit, so the strongest claim should remain at the level of the Green field rather than the raw shell derivative

## Next Step

The natural follow-up is now explicit:

1. couple the measured graph Green field directly into the bound-state solver
2. test whether the full chain closes numerically:

$$
\text{interaction kernel}
\rightarrow
\text{graph Laplacian}
\rightarrow
\text{graph Green field}
\rightarrow
\text{bound ladder}
$$
