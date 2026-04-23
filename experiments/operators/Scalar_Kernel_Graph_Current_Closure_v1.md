# Scalar Kernel Graph Current Closure v1

## Purpose

Push the scalar carrier from the static `52.2` inverse-square response law to the next stronger question:

> does the same scalar carrier support a transient response/current closure on the clean refinement line?

The carrier is unchanged:

- same local scalar kernel graph
- same `51.4` local metric field
- same bounded Euclidean shell reconstruction

The new ingredient is transient heat-flow current.

## Construction

For each clean cubic refinement:

1. evolve the point-source heat trajectory on the same scalar operator;
2. build shell-mean profiles on the clean shell scaffold;
3. infer empirical shell current from cumulative mass flow,
4. compare it against a constitutive shell current built from the same metric field and shell gradient:

$$
I_{const}(r,t) = - \kappa \; 4\pi r^2 g_{rr}(r) \; \partial_r u(r,t)
$$

where `\kappa` is fitted per refinement on the bounded time window.

This tests whether one transient current law closes on the same scalar carrier, not just whether the static Green response is inverse-square.

## Artifacts

- script: `numerics/simulations/scalar_kernel_graph_current_closure.py`
- config: `config.json -> scalar_kernel_graph_current_closure`
- results: `data/20260421_131559_scalar_kernel_graph_current_closure.json`
- latest: `data/scalar_kernel_graph_current_closure_latest.json`
- plots:
  - `plots/20260421_131559_scalar_kernel_graph_current_closure_summary.png`
  - `plots/20260421_131559_scalar_kernel_graph_current_closure_profiles.png`

## Clean refinement scan

| case | `D_eff` | `kappa_fit` | `kappa_fit / D_eff` | median rel. error | p90 rel. error | status |
| --- | --- | --- | --- | --- | --- | --- |
| `clean|n=13` | `1.0593` | `655.0391` | `618.3748` | `0.5925` | `0.6671` | OPEN |
| `clean|n=15` | `1.0466` | `1334.1709` | `1274.7134` | `0.5090` | `0.5997` | OPEN |
| `clean|n=17` | `1.0368` | `2117.9635` | `2042.7973` | `0.4850` | `0.6198` | OPEN |

## Normalized profile drift

| pair | drift |
| --- | --- |
| `13->15` | `0.1874` |
| `15->17` | `0.2073` |

## Interpretation

- observation: the scalar carrier and local metric field now support a bounded inverse-square response law, but the stronger transient response/current closure remains open on the clean refinement scan
- conclusion: open current-closure cases remain: failing_cases=['clean|n=13', 'clean|n=15', 'clean|n=17'], max_profile_drift=0.2073; the weakest case is `clean|n=13` with median relative error 0.5925, p90 error 0.6671, and conductivity/diffusivity ratio 618.3748; the correct boundary therefore stays at the scalar carrier, local metric field, and shell-native inverse-square response law without yet promoting a positive transient current-closure claim

## Current boundary

If this note is read negatively, the correct bounded statement is:

> the same scalar carrier now supports the static `52.2` inverse-square response law, but transient constitutive current closure remains open even on the clean refinement scan.

So the next honest move after this note is not a stronger claim. It is:

1. improve the transient shell-current reconstruction,
2. understand the scaling of the fitted constitutive coefficient,
3. then only afterward widen back toward disorder or kernel-family robustness.
