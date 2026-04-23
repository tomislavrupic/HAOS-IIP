# Scalar Kernel Graph Current Closure Shell-Native v1

## Purpose

Retest the open transient current-closure question on the same scalar carrier without changing the carrier, the `51.4` local metric field, or the clean refinement line.

The only change is the reconstruction:

- keep the same local scalar kernel graph
- keep the same coarse local metric field
- replace the nearest-shell profile fit with exact bulk shells on the clean scaffold
- read the constitutive law against shell density rather than raw node mass

This is a shell-native reconstruction audit, not a new substrate claim.

## Construction

For each clean cubic refinement:

1. evolve the same point-source heat trajectory;
2. build exact bulk shells from the clean radius labels;
3. convert shell mass to shell density via `rho = m_shell / (count * h^3)`;
4. compare empirical cumulative shell current against

$$
I_{const}(r,t) = - \kappa \; 4\pi r^2 g_{rr}(r) \; \partial_r \rho(r,t).
$$

The closure diagnostic is no longer the raw current profile alone. It is whether one bounded `\kappa` family stays close to the heat diffusivity and remains stable across radius and refinement.

## Artifacts

- script: `numerics/simulations/scalar_kernel_graph_current_closure_shell_native.py`
- config: `config.json -> scalar_kernel_graph_current_closure_shell_native`
- results: `data/20260422_204609_scalar_kernel_graph_current_closure_shell_native.json`
- latest: `data/scalar_kernel_graph_current_closure_shell_native_latest.json`
- plots:
  - `plots/20260422_204609_scalar_kernel_graph_current_closure_shell_native_summary.png`
  - `plots/20260422_204609_scalar_kernel_graph_current_closure_shell_native_profiles.png`

## Clean refinement scan

| case | `D_eff` | `kappa_fit` | `kappa_fit / D_eff` | median rel. error | p90 rel. error | shell `kappa` CV | status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `clean|n=13` | `1.0593` | `0.9397` | `0.8871` | `0.0705` | `0.2328` | `0.0967` | PASS |
| `clean|n=15` | `1.0466` | `0.9530` | `0.9106` | `0.0709` | `0.1984` | `0.0742` | PASS |
| `clean|n=17` | `1.0368` | `0.9604` | `0.9263` | `0.0739` | `0.2718` | `0.1176` | PASS |

## Normalized shell-kappa drift

| pair | drift |
| --- | --- |
| `13->15` | `0.0493` |
| `15->17` | `0.0966` |

## Interpretation

- observation: on the clean scalar carrier, a shell-native transient current reconstruction now closes onto one bounded constitutive family: exact bulk shells plus density-normalized gradients recover a stable effective conductivity close to the heat diffusivity
- conclusion: all 3 clean refinement cases pass, the maximum normalized shell-kappa drift is 0.0966, and the weakest passing case is `clean|n=17` with kappa/D_eff 0.9263, median relative error 0.0739, p90 error 0.2718, and shell-kappa CV 0.1176; this supports a bounded shell-native transient current-closure read on the same scalar carrier

## Current boundary

If this note is read positively, the correct bounded statement is:

> on the clean scalar carrier, a shell-native transient current reconstruction closes onto one bounded constitutive family whose fitted conductivity stays close to the heat diffusivity and whose shellwise profile remains refinement-stable.

This note still does **not** claim:

- disorder or kernel-family universality
- conserved-current closure
- curvature extraction
- spacetime or ontology
