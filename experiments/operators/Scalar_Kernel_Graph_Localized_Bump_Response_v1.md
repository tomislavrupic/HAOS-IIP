# Scalar Kernel Graph Localized Bump Response v1

## Purpose

Resolve the broad `52.5` localized-bump OPEN row into a sharper threshold statement.

The original bump test used the same early transient window as the smooth radial branch. That made the localized source-core layer dominate the shellwise constitutive fit and produced a maximum refinement drift of `0.5124`.

This follow-up keeps the same scalar carrier and the same bump shape, but excludes the earliest source-core transient layer with fit window `[0.01, 0.026]`.

## Construction

- bump sigma: `0.18`
- refinements: `[13, 15]`
- weak eta values: `[0.01, 0.02, 0.03, 0.04, 0.05]`
- stress eta values: `[0.075, 0.1, 0.15]`
- thresholds: `{'max_median_relative_error': 0.1, 'max_p90_relative_error': 0.37, 'max_diffusivity_match_gap': 0.15, 'max_shell_kappa_cv': 0.14, 'max_refinement_profile_drift': 0.15, 'min_metric_tracking_abs_corr': 0.95}`

## Regime Summary

| regime | eta values | all cases pass | max drift | max median error | max p90 error | max shell-kappa CV |
| --- | --- | --- | --- | --- | --- | --- |
| `weak` | `[0.01, 0.02, 0.03, 0.04, 0.05]` | `True` | `0.1095` | `0.0943` | `0.2984` | `0.1248` |
| `stress` | `[0.075, 0.1, 0.15]` | `False` | `0.3578` | `0.2784` | `0.5287` | `0.3406` |

## Case Table

| case | `kappa / D_eff` | median error | p90 error | shell-kappa CV | metric corr | status |
| --- | --- | --- | --- | --- | --- | --- |
| `stress_bump|n=13|eta=0.075` | `0.9374` | `0.1404` | `0.3275` | `0.1545` | `0.9782` | OPEN |
| `stress_bump|n=15|eta=0.075` | `0.9089` | `0.1112` | `0.2777` | `0.1519` | `0.9903` | OPEN |
| `stress_bump|n=13|eta=0.100` | `0.8893` | `0.1717` | `0.3847` | `0.2079` | `0.9743` | OPEN |
| `stress_bump|n=15|eta=0.100` | `0.8626` | `0.1412` | `0.3498` | `0.2020` | `0.9909` | OPEN |
| `stress_bump|n=13|eta=0.150` | `0.7787` | `0.2784` | `0.5287` | `0.3406` | `0.9628` | OPEN |
| `stress_bump|n=15|eta=0.150` | `0.7523` | `0.2374` | `0.4845` | `0.3165` | `0.9853` | OPEN |
| `weak_bump|n=13|eta=0.010` | `1.0104` | `0.0465` | `0.1323` | `0.0478` | `0.9891` | PASS |
| `weak_bump|n=15|eta=0.010` | `1.0029` | `0.0541` | `0.1392` | `0.0621` | `0.9892` | PASS |
| `weak_bump|n=13|eta=0.020` | `1.0000` | `0.0593` | `0.1887` | `0.0729` | `0.9889` | PASS |
| `weak_bump|n=15|eta=0.020` | `0.9926` | `0.0603` | `0.1573` | `0.0713` | `0.9890` | PASS |
| `weak_bump|n=13|eta=0.030` | `0.9931` | `0.0659` | `0.1796` | `0.0735` | `0.9825` | PASS |
| `weak_bump|n=15|eta=0.030` | `0.9739` | `0.0719` | `0.1802` | `0.0873` | `0.9877` | PASS |
| `weak_bump|n=13|eta=0.040` | `0.9833` | `0.0658` | `0.2414` | `0.0968` | `0.9702` | PASS |
| `weak_bump|n=15|eta=0.040` | `0.9593` | `0.0807` | `0.2001` | `0.1022` | `0.9931` | PASS |
| `weak_bump|n=13|eta=0.050` | `0.9699` | `0.0943` | `0.2984` | `0.1248` | `0.9752` | PASS |
| `weak_bump|n=15|eta=0.050` | `0.9492` | `0.0902` | `0.2201` | `0.1142` | `0.9924` | PASS |

## Refinement Drift

| case | drift |
| --- | --- |
| `stress_bump|eta=0.075` | `0.1474` |
| `stress_bump|eta=0.100` | `0.2019` |
| `stress_bump|eta=0.150` | `0.3578` |
| `weak_bump|eta=0.010` | `0.0568` |
| `weak_bump|eta=0.020` | `0.0728` |
| `weak_bump|eta=0.030` | `0.0834` |
| `weak_bump|eta=0.040` | `0.0928` |
| `weak_bump|eta=0.050` | `0.1095` |

## Interpretation

- observation: after excluding the earliest source-core transient layer, weak localized bump excitations close as a bounded metric-tracking transport response, while stronger localized bumps remain outside the same closure
- conclusion: all 10 weak localized bump cases pass for eta values [0.01, 0.02, 0.03, 0.04, 0.05], with maximum refinement drift 0.1095; the weakest weak case is `weak_bump|n=13|eta=0.050` with kappa/D_eff 0.9699, median relative error 0.0943, p90 error 0.2984, shell-kappa CV 0.1248, and |metric-tracking corr| 0.9752; the stress window remains open with maximum drift 0.3578, and the hardest stress case is `stress_bump|n=13|eta=0.150` with median relative error 0.2784, p90 error 0.5287, and shell-kappa CV 0.3406

## Authority

- script: `numerics/simulations/scalar_kernel_graph_localized_bump_response.py`
- result: `data/20260423_182539_scalar_kernel_graph_localized_bump_response.json`
- latest: `data/scalar_kernel_graph_localized_bump_response_latest.json`
- plots:
  - `plots/20260423_182539_scalar_kernel_graph_localized_bump_response_summary.png`
  - `plots/20260423_182539_scalar_kernel_graph_localized_bump_response_profiles.png`
