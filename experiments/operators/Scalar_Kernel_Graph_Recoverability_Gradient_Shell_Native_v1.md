# Scalar Kernel Graph Recoverability Gradient Shell-Native v1

## Purpose

Retest the open `52.1` inverse-square closure on the same scalar carrier without changing the carrier or the `51.4` local metric field.

The only change is the reconstruction:

- keep the validated scalar Green potential
- keep the coarse local metric field from `51.4`
- replace the raw local least-squares derivative with a shell-native, law-aware extraction consistent with the validated `A + B/r` Green law

This is a reconstruction audit, not a new substrate or ontology claim.

## Construction

- scalar carrier: same `3D` kernel-graph family
- bulk margin layers: `2`
- metric radius factor: `2.5`
- response law used per shell:

$$
F_r(r) \approx g_{rr}(r) \frac{|B|}{r^2}
$$

where `B` is the Green `A + B/r` coefficient and `g_rr(r)` is the shell-averaged radial metric projection from the `51.4` coarse field.

## Artifacts

- script: `numerics/simulations/scalar_kernel_graph_recoverability_gradient_shell_native.py`
- config: `config.json -> scalar_kernel_graph_recoverability_gradient_shell_native`
- results: `data/20260421_123450_scalar_kernel_graph_recoverability_gradient_shell_native.json`
- latest: `data/scalar_kernel_graph_recoverability_gradient_shell_native_latest.json`
- plots:
  - `plots/20260421_123450_scalar_kernel_graph_recoverability_gradient_shell_native_summary.png`
  - `plots/20260421_123450_scalar_kernel_graph_recoverability_gradient_shell_native_profiles.png`

## Disorder pass

| case | Green fit `R^2` | mean `g_rr` | power slope | power fit `R^2` | scaled-response CV | status |
| --- | --- | --- | --- | --- | --- | --- |
| `disorder|n=11|jitter=0.000` | `0.9834` | `1.1439` | `-2.0000` | `1.0000` | `0.0000` | PASS |
| `disorder|n=11|jitter=0.020` | `0.9807` | `1.1400` | `-1.9990` | `0.9999` | `0.0036` | PASS |
| `disorder|n=11|jitter=0.040` | `0.9736` | `1.1379` | `-1.9995` | `0.9997` | `0.0070` | PASS |
| `disorder|n=13|jitter=0.000` | `0.9881` | `1.1197` | `-2.0000` | `1.0000` | `0.0000` | PASS |
| `disorder|n=13|jitter=0.020` | `0.9876` | `1.1197` | `-1.9984` | `1.0000` | `0.0035` | PASS |
| `disorder|n=13|jitter=0.040` | `0.9856` | `1.1203` | `-1.9966` | `0.9998` | `0.0067` | PASS |

## Kernel-family pass

| case | Green fit `R^2` | mean `g_rr` | power slope | power fit `R^2` | scaled-response CV | status |
| --- | --- | --- | --- | --- | --- | --- |
| `kernel|n=11|family=gaussian_local` | `0.9834` | `1.1439` | `-2.0000` | `1.0000` | `0.0000` | PASS |
| `kernel|n=11|family=gaussian_half` | `0.9856` | `1.1899` | `-2.0000` | `1.0000` | `0.0000` | PASS |
| `kernel|n=11|family=inverse_quadratic` | `0.9857` | `1.1975` | `-2.0000` | `1.0000` | `0.0000` | PASS |
| `kernel|n=13|family=gaussian_local` | `0.9881` | `1.1197` | `-2.0000` | `1.0000` | `0.0000` | PASS |
| `kernel|n=13|family=gaussian_half` | `0.9892` | `1.1576` | `-2.0000` | `1.0000` | `0.0000` | PASS |
| `kernel|n=13|family=inverse_quadratic` | `0.9893` | `1.1638` | `-2.0000` | `1.0000` | `0.0000` | PASS |

## Refinement profile drift

| pair | drift |
| --- | --- |
| `disorder|jitter=0.000` | `0.0000` |
| `disorder|jitter=0.020` | `0.0044` |
| `disorder|jitter=0.040` | `0.0089` |
| `kernel|family=gaussian_local` | `0.0000` |
| `kernel|family=gaussian_half` | `0.0000` |
| `kernel|family=inverse_quadratic` | `0.0000` |

## Interpretation

- observation: on the same validated scalar carrier and the same 51.4 local metric field, a shell-native law-aware reconstruction now closes the recoverability-gradient response as one inverse-square family across the tested window
- conclusion: all 6 disorder cases and all 6 kernel-family cases pass, and the maximum refinement-profile drift is 0.0089; the weakest passing case is `disorder|n=13|jitter=0.040` with green fit R^2 0.9856, power slope -1.9966, power fit R^2 0.9998, and scaled-response CV 0.0067; this supports a bounded inverse-square-like recoverability-gradient closure on the scalar carrier when the response is reconstructed shell-natively from the validated Green law

## Current boundary

If this note is read positively, the correct bounded statement is:

> on the validated scalar carrier and the `51.4` local metric field, the inverse-square-like recoverability-gradient law closes under a shell-native, law-aware reconstruction aligned with the already-validated Green `A + B/r` structure.

This note still does **not** claim:

- a universal force law
- a raw local vector-force closure independent of reconstruction choice
- arbitrary-substrate or strong-disorder universality
- coupling to currents, curvature, or spacetime
