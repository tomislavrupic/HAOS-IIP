# Scalar Kernel Graph Power-Law Scaling v1

## Purpose

Resolve the broad physics-bridge `power-law scaling` OPEN row into two auditable readouts:

- raw local-gradient power scaling, which remains `OPEN`;
- shell-native law-aware power scaling, which passes as a bounded continuum-like proxy.

This is a post-processing audit over frozen scalar artifacts. It does not change HAOS core, the scalar carrier, or either recoverability-gradient runner.

## Source Artifacts

- raw local gradient: `data/scalar_kernel_graph_recoverability_gradient_latest.json`
- shell-native reconstruction: `data/scalar_kernel_graph_recoverability_gradient_shell_native_latest.json`
- result: `data/20260423_183612_scalar_kernel_graph_power_law_scaling.json`
- latest: `data/scalar_kernel_graph_power_law_scaling_latest.json`
- plots:
  - `plots/20260423_183612_scalar_kernel_graph_power_law_scaling_summary.png`
  - `plots/20260423_183612_scalar_kernel_graph_power_law_scaling_profiles.png`

## Summary

| readout | cases | status | min power-fit `R^2` | max `|slope + 2|` | max scaled CV | max profile drift |
| --- | --- | --- | --- | --- | --- | --- |
| raw local gradient | `12` | OPEN | `0.3220` | `1.2651` | `0.3139` | `0.5249` |
| shell native | `12` | PASS | `0.9997` | `0.0034` | `0.0070` | `0.0089` |

## Weakest Cases

- raw: `disorder|n=11|jitter=0.020`, slope `-0.7349`, `R^2 = 0.3220`, CV `0.2880`;
- shell-native: `disorder|n=13|jitter=0.040`, Green `R^2 = 0.9856`, slope `-1.9966`, power `R^2 = 0.9998`, CV `0.0067`.

## Interpretation

- observation: the broad power-law scaling row splits cleanly: the raw local-gradient read remains open, while the shell-native law-aware read passes the tested continuum-like scaling thresholds
- conclusion: raw local-gradient power scaling remains open with minimum R^2 0.3220, maximum |slope+2| 1.2651, and profile drift 0.5249; shell-native scaling passes with minimum power-fit R^2 0.9997, maximum |slope+2| 0.0034, scaled-response CV 0.0070, and profile drift 0.0089; this supports a bounded continuum-like power-law proxy only under the shell-native reconstruction

## Correct Boundary

The old broad statement was:

> power-law scaling remains open because raw power-fit R^2 is too low.

The new statement is:

> raw local-gradient power scaling remains open, but shell-native law-aware power scaling closes as a bounded continuum-like proxy on the tested scalar-carrier window.

This is not a continuum-limit proof and not a gravity claim. It is a reconstruction boundary.
