# Scalar Kernel Graph Current Closure Radial Disorder Native Flux v1

## Purpose

Test the missing transport bridge after the scaffold-shell disorder failure:

> if smooth radial mild-disorder cases are read with a disorder-native interior-ball flux instead of scaffold-shell bookkeeping, does the bounded constitutive law return?

The carrier stays fixed:

- same scalar kernel graph family
- same smooth radial deformation family
- same mild-disorder window
- same constitutive target

What changes is only the transport readout:

- native bulk shells built from the actual disordered radii
- interior-ball cumulative mass on the actual bulk nodes
- local-linear radial density gradient
- delayed asymptotic fit window `[0.006, 0.02]` to exclude the earliest source-core transient layer

## Construction

- clean refinements: `[13, 15]`
- radial amplitudes: `[0.05, 0.1, 0.15]`
- disorder fractions of lattice spacing: `[0.02, 0.04]`
- seeds: `[0, 1, 2]`
- gradient mode: `local_linear`

Artifacts:

- script: `numerics/simulations/scalar_kernel_graph_current_closure_radial_disorder_native_flux.py`
- config: `config.json -> scalar_kernel_graph_current_closure_radial_disorder_native_flux`
- results: `data/20260422_215916_scalar_kernel_graph_current_closure_radial_disorder_native_flux.json`
- latest: `data/scalar_kernel_graph_current_closure_radial_disorder_native_flux_latest.json`
- plots:
  - `plots/20260422_215916_scalar_kernel_graph_current_closure_radial_disorder_native_flux_summary.png`
  - `plots/20260422_215916_scalar_kernel_graph_current_closure_radial_disorder_native_flux_profiles.png`

## Aggregated cases

| case | trial passes | `kappa / D_eff` | median rel. error | p90 rel. error | shell `kappa` CV | `|metric corr|` | status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `radial_native_flux|n=13|eta=0.050|sigma=0.020` | `3/3` | `0.9141` | `0.0514` | `0.2355` | `0.1030` | `0.9985` | PASS |
| `radial_native_flux|n=15|eta=0.050|sigma=0.020` | `3/3` | `0.9513` | `0.0368` | `0.2285` | `0.0837` | `0.9947` | PASS |
| `radial_native_flux|n=13|eta=0.100|sigma=0.020` | `3/3` | `0.9360` | `0.0583` | `0.2654` | `0.1101` | `0.9977` | PASS |
| `radial_native_flux|n=15|eta=0.100|sigma=0.020` | `3/3` | `0.9734` | `0.0515` | `0.2559` | `0.0953` | `0.9939` | PASS |
| `radial_native_flux|n=13|eta=0.150|sigma=0.020` | `3/3` | `0.9553` | `0.0801` | `0.3004` | `0.1277` | `0.9944` | PASS |
| `radial_native_flux|n=15|eta=0.150|sigma=0.020` | `3/3` | `1.0017` | `0.0637` | `0.2792` | `0.1024` | `0.9902` | PASS |
| `radial_native_flux|n=13|eta=0.050|sigma=0.040` | `3/3` | `0.9118` | `0.0677` | `0.2494` | `0.1178` | `0.9926` | PASS |
| `radial_native_flux|n=15|eta=0.050|sigma=0.040` | `3/3` | `0.9425` | `0.0520` | `0.2354` | `0.0890` | `0.9932` | PASS |
| `radial_native_flux|n=13|eta=0.100|sigma=0.040` | `3/3` | `0.9335` | `0.0770` | `0.2777` | `0.1286` | `0.9966` | PASS |
| `radial_native_flux|n=15|eta=0.100|sigma=0.040` | `3/3` | `0.9673` | `0.0675` | `0.2492` | `0.0947` | `0.9931` | PASS |
| `radial_native_flux|n=13|eta=0.150|sigma=0.040` | `2/3` | `0.9575` | `0.0911` | `0.2950` | `0.1369` | `0.9952` | PASS |
| `radial_native_flux|n=15|eta=0.150|sigma=0.040` | `3/3` | `0.9959` | `0.0736` | `0.2818` | `0.1064` | `0.9916` | PASS |

## Refinement drift

| case | drift |
| --- | --- |
| `sigma=0.020|eta=0.050` | `0.0988` |
| `sigma=0.020|eta=0.100` | `0.0902` |
| `sigma=0.020|eta=0.150` | `0.1061` |
| `sigma=0.040|eta=0.050` | `0.1140` |
| `sigma=0.040|eta=0.100` | `0.1233` |
| `sigma=0.040|eta=0.150` | `0.1278` |

## Interpretation

- observation: once the smooth-radial mild-disorder family is read with native bulk shells, interior-ball cumulative flux, and a delayed asymptotic fit window, the bounded constitutive closure returns on aggregated observables
- conclusion: all 12 aggregated mild-disorder cases pass, with maximum refinement drift 0.1278; the weakest aggregated case is `radial_native_flux|n=13|eta=0.150|sigma=0.040` with kappa/D_eff 0.9575, median relative error 0.0911, p90 error 0.2950, shell-kappa CV 0.1369, and |metric-tracking corr| 0.9952; the hardest seed-level trial is `radial_native_flux|n=13|eta=0.150|sigma=0.040|seed=1` with median relative error 0.1100, p90 error 0.2817, shell-kappa CV 0.1482, and pass state FALSE; this supports a bounded disorder-native transient closure on the smooth-radial branch, with the minimum seed-level pass count still only 2/3 in the hardest tested case
