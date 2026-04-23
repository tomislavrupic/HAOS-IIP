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
- delayed asymptotic fit window `[0.006, 0.024]` to exclude the earliest source-core transient layer

## Construction

- clean refinements: `[13, 15]`
- radial amplitudes: `[0.05, 0.1, 0.15]`
- disorder fractions of lattice spacing: `[0.02, 0.04]`
- seeds: `[0, 1, 2]`
- gradient mode: `local_linear`

Artifacts:

- script: `numerics/simulations/scalar_kernel_graph_current_closure_radial_disorder_native_flux.py`
- config: `config.json -> scalar_kernel_graph_current_closure_radial_disorder_native_flux`
- results: `data/20260423_184421_scalar_kernel_graph_current_closure_radial_disorder_native_flux.json`
- latest: `data/scalar_kernel_graph_current_closure_radial_disorder_native_flux_latest.json`
- plots:
  - `plots/20260423_184421_scalar_kernel_graph_current_closure_radial_disorder_native_flux_summary.png`
  - `plots/20260423_184421_scalar_kernel_graph_current_closure_radial_disorder_native_flux_profiles.png`

## Aggregated cases

| case | trial passes | `kappa / D_eff` | median rel. error | p90 rel. error | shell `kappa` CV | `|metric corr|` | status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `radial_native_flux|n=13|eta=0.050|sigma=0.020` | `3/3` | `0.9387` | `0.0587` | `0.2323` | `0.0985` | `0.9985` | PASS |
| `radial_native_flux|n=15|eta=0.050|sigma=0.020` | `3/3` | `0.9769` | `0.0418` | `0.2277` | `0.0853` | `0.9947` | PASS |
| `radial_native_flux|n=13|eta=0.100|sigma=0.020` | `3/3` | `0.9564` | `0.0658` | `0.2579` | `0.1045` | `0.9977` | PASS |
| `radial_native_flux|n=15|eta=0.100|sigma=0.020` | `3/3` | `0.9969` | `0.0545` | `0.2453` | `0.0925` | `0.9939` | PASS |
| `radial_native_flux|n=13|eta=0.150|sigma=0.020` | `3/3` | `0.9737` | `0.0821` | `0.2876` | `0.1177` | `0.9944` | PASS |
| `radial_native_flux|n=15|eta=0.150|sigma=0.020` | `3/3` | `1.0245` | `0.0617` | `0.2665` | `0.0998` | `0.9902` | PASS |
| `radial_native_flux|n=13|eta=0.050|sigma=0.040` | `3/3` | `0.9363` | `0.0680` | `0.2406` | `0.1056` | `0.9926` | PASS |
| `radial_native_flux|n=15|eta=0.050|sigma=0.040` | `3/3` | `0.9677` | `0.0567` | `0.2250` | `0.0899` | `0.9932` | PASS |
| `radial_native_flux|n=13|eta=0.100|sigma=0.040` | `3/3` | `0.9551` | `0.0713` | `0.2673` | `0.1178` | `0.9966` | PASS |
| `radial_native_flux|n=15|eta=0.100|sigma=0.040` | `3/3` | `0.9902` | `0.0728` | `0.2397` | `0.0962` | `0.9931` | PASS |
| `radial_native_flux|n=13|eta=0.150|sigma=0.040` | `3/3` | `0.9765` | `0.0796` | `0.2827` | `0.1242` | `0.9952` | PASS |
| `radial_native_flux|n=15|eta=0.150|sigma=0.040` | `3/3` | `1.0179` | `0.0789` | `0.2677` | `0.1057` | `0.9916` | PASS |

## Refinement drift

| case | drift |
| --- | --- |
| `sigma=0.020|eta=0.050` | `0.1024` |
| `sigma=0.020|eta=0.100` | `0.0926` |
| `sigma=0.020|eta=0.150` | `0.1059` |
| `sigma=0.040|eta=0.050` | `0.1111` |
| `sigma=0.040|eta=0.100` | `0.1217` |
| `sigma=0.040|eta=0.150` | `0.1249` |

## Interpretation

- observation: once the smooth-radial mild-disorder family is read with native bulk shells, interior-ball cumulative flux, and a delayed asymptotic fit window, the bounded constitutive closure returns on aggregated observables
- conclusion: all 12 aggregated mild-disorder cases pass, with maximum refinement drift 0.1249; the weakest aggregated case is `radial_native_flux|n=13|eta=0.150|sigma=0.020` with kappa/D_eff 0.9737, median relative error 0.0821, p90 error 0.2876, shell-kappa CV 0.1177, and |metric-tracking corr| 0.9944; the hardest seed-level trial is `radial_native_flux|n=13|eta=0.150|sigma=0.020|seed=1` with median relative error 0.0895, p90 error 0.2660, shell-kappa CV 0.1092, and pass state TRUE; this supports a bounded disorder-native transient closure on the smooth-radial branch, with the minimum seed-level pass count still only 3/3 in the hardest tested case
