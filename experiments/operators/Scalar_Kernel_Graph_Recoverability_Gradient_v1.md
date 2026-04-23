# Scalar Kernel Graph Recoverability Gradient v1

## Purpose

Build the first bounded `52`-line force-like observable on the validated scalar geometry carrier without importing a new stochastic recovery model.

The deterministic proxy used here is simple:

- recoverability potential: the same mean-zero scalar Green field already validated on the carrier
- local metric: the coarse local metric field from `51.4`
- effective recoverability-gradient field:

$$
J_{rec} = - G \nabla \phi
$$

This is not yet called a fundamental force. It is the first bounded operator-native response field on the same carrier.

## Construction

- scalar carrier: same `3D` kernel-graph family
- metric radius factor: `2.5`
- gradient radius factor: `2.5`
- bulk margin layers: `2`

The bounded test asks whether the coarse response field is:

1. strongly radial
2. close to inverse-square in shell profile
3. close to constant in `r^2 F_r`
4. stable under the same mild disorder and bounded kernel-family window

## Artifacts

- script: `numerics/simulations/scalar_kernel_graph_recoverability_gradient.py`
- config: `config.json -> scalar_kernel_graph_recoverability_gradient`
- results: `data/20260420_223721_scalar_kernel_graph_recoverability_gradient.json`
- latest: `data/scalar_kernel_graph_recoverability_gradient_latest.json`
- plots:
  - `plots/20260420_223721_scalar_kernel_graph_recoverability_gradient_summary.png`
  - `plots/20260420_223721_scalar_kernel_graph_recoverability_gradient_profiles.png`

## Disorder pass

| case | radial alignment | power slope | power fit `R^2` | scaled-flux CV | status |
| --- | --- | --- | --- | --- | --- |
| `disorder|n=11|jitter=0.000` | `0.9978` | `-0.8953` | `0.5966` | `0.3139` | OPEN |
| `disorder|n=11|jitter=0.020` | `0.9941` | `-0.7349` | `0.3220` | `0.2880` | OPEN |
| `disorder|n=11|jitter=0.040` | `0.9909` | `-0.8087` | `0.3736` | `0.2877` | OPEN |
| `disorder|n=13|jitter=0.000` | `0.9983` | `-1.2780` | `0.8090` | `0.2312` | OPEN |
| `disorder|n=13|jitter=0.020` | `0.9968` | `-1.4371` | `0.8278` | `0.1659` | OPEN |
| `disorder|n=13|jitter=0.040` | `0.9963` | `-1.4518` | `0.8214` | `0.1688` | OPEN |

## Kernel-family pass

| case | radial alignment | power slope | power fit `R^2` | scaled-flux CV | status |
| --- | --- | --- | --- | --- | --- |
| `kernel|n=11|family=gaussian_local` | `0.9978` | `-0.8953` | `0.5966` | `0.3139` | OPEN |
| `kernel|n=11|family=gaussian_half` | `0.9983` | `-0.9147` | `0.5978` | `0.3131` | OPEN |
| `kernel|n=11|family=inverse_quadratic` | `0.9984` | `-0.9181` | `0.5982` | `0.3132` | OPEN |
| `kernel|n=13|family=gaussian_local` | `0.9983` | `-1.2780` | `0.8090` | `0.2312` | OPEN |
| `kernel|n=13|family=gaussian_half` | `0.9986` | `-1.2890` | `0.8112` | `0.2302` | OPEN |
| `kernel|n=13|family=inverse_quadratic` | `0.9986` | `-1.2908` | `0.8116` | `0.2302` | OPEN |

## Refinement profile drift

| pair | drift |
| --- | --- |
| `disorder|jitter=0.000` | `0.5009` |
| `disorder|jitter=0.020` | `0.5249` |
| `disorder|jitter=0.040` | `0.5212` |
| `kernel|family=gaussian_local` | `0.5009` |
| `kernel|family=gaussian_half` | `0.5022` |
| `kernel|family=inverse_quadratic` | `0.5024` |

## Interpretation

- observation: the scalar Green potential supports a geometry carrier and a local metric field, but the recoverability-gradient response does not yet close as one stable inverse-square family across the tested window
- conclusion: open recoverability-gradient cases remain: failing_cases=['disorder|n=11|jitter=0.000', 'disorder|n=11|jitter=0.020', 'disorder|n=11|jitter=0.040', 'disorder|n=13|jitter=0.000', 'disorder|n=13|jitter=0.020', 'disorder|n=13|jitter=0.040', 'kernel|n=11|family=gaussian_local', 'kernel|n=11|family=gaussian_half', 'kernel|n=11|family=inverse_quadratic', 'kernel|n=13|family=gaussian_local', 'kernel|n=13|family=gaussian_half', 'kernel|n=13|family=inverse_quadratic'], max_profile_drift=0.5249; the correct boundary therefore stays at the scalar geometry carrier and local metric-field statement without yet promoting a positive force-like response claim

## Current boundary

The correct bounded statement for this note is negative:

> on the validated scalar carrier, the first raw local least-squares recoverability-gradient reconstruction remains open; it gives a strongly radial response field, but it does not yet close the inverse-square law on its own.

The positive follow-on belongs to the separate `52.2` shell-native / law-aware reconstruction, not to this raw-local note.

This note still does **not** claim:

- a universal force law
- a fundamental ontic force
- closure on arbitrary substrates or strong disorder
- coupling to currents, curvature, or spacetime
