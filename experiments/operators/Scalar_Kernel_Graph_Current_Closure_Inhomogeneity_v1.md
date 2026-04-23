# Scalar Kernel Graph Current Closure Inhomogeneity v1

## Purpose

Test the next honest follow-on after the clean `52.4` shell-native current closure:

> do the extracted metric and the shell-native transport law track a designed scalar inhomogeneity together on the same carrier?

The carrier stays fixed:

- same scalar kernel graph family
- same local metric-field construction
- same shell-native current law

The only new ingredient is a small deterministic deformation profile on the point cloud.

## Construction

- families: `['radial', 'bump']`
- amplitudes: `[0.05, 0.1, 0.15]`
- clean refinements: `[13, 15]`

For each case:

1. deform the same clean carrier by one designed radial profile;
2. rebuild the scalar operator on the deformed point cloud;
3. extract the coarse local metric field;
4. test whether the shell-native transient current closure still fits one bounded constitutive family;
5. compare the shellwise metric shift against the imposed deformation profile.

Artifacts:

- script: `numerics/simulations/scalar_kernel_graph_current_closure_inhomogeneity.py`
- config: `config.json -> scalar_kernel_graph_current_closure_inhomogeneity`
- results: `data/20260422_213513_scalar_kernel_graph_current_closure_inhomogeneity.json`
- latest: `data/scalar_kernel_graph_current_closure_inhomogeneity_latest.json`
- plots:
  - `plots/20260422_213513_scalar_kernel_graph_current_closure_inhomogeneity_summary.png`
  - `plots/20260422_213513_scalar_kernel_graph_current_closure_inhomogeneity_profiles.png`

## Smooth Radial Deformation

| case | `kappa / D_eff` | median rel. error | p90 rel. error | shell `kappa` CV | metric corr | status |
| --- | --- | --- | --- | --- | --- | --- |
| `radial|n=13|eta=0.050` | `0.9019` | `0.0704` | `0.3240` | `0.1193` | `-0.9885` | PASS |
| `radial|n=13|eta=0.100` | `0.9243` | `0.0665` | `0.3430` | `0.1265` | `-0.9939` | PASS |
| `radial|n=13|eta=0.150` | `0.9466` | `0.0568` | `0.3568` | `0.1310` | `-0.9840` | PASS |
| `radial|n=15|eta=0.050` | `0.9272` | `0.0539` | `0.2429` | `0.1072` | `-0.9872` | PASS |
| `radial|n=15|eta=0.100` | `0.9588` | `0.0496` | `0.2686` | `0.1178` | `-0.9935` | PASS |
| `radial|n=15|eta=0.150` | `0.9848` | `0.0562` | `0.2773` | `0.1215` | `-0.9846` | PASS |
## Localized Scalar Bump

| case | `kappa / D_eff` | median rel. error | p90 rel. error | shell `kappa` CV | metric corr | status |
| --- | --- | --- | --- | --- | --- | --- |
| `bump|n=13|eta=0.050` | `0.7857` | `0.2053` | `0.4658` | `0.2164` | `-0.9752` | OPEN |
| `bump|n=13|eta=0.100` | `0.6876` | `0.3113` | `0.6112` | `0.3331` | `-0.9743` | OPEN |
| `bump|n=13|eta=0.150` | `0.5889` | `0.4299` | `0.7346` | `0.4963` | `-0.9628` | OPEN |
| `bump|n=15|eta=0.050` | `0.7971` | `0.1569` | `0.3896` | `0.1944` | `-0.9924` | OPEN |
| `bump|n=15|eta=0.100` | `0.6929` | `0.2487` | `0.5368` | `0.3088` | `-0.9909` | OPEN |
| `bump|n=15|eta=0.150` | `0.5929` | `0.3565` | `0.6609` | `0.4453` | `-0.9853` | OPEN |

## Refinement drift

| case | drift |
| --- | --- |
| `radial|eta=0.050` | `0.1202` |
| `radial|eta=0.100` | `0.1384` |
| `radial|eta=0.150` | `0.1465` |
| `bump|eta=0.050` | `0.1922` |
| `bump|eta=0.100` | `0.3150` |
| `bump|eta=0.150` | `0.5124` |

## Interpretation

- observation: the extracted metric and the shell-native constitutive law track one another under smooth radial inhomogeneity on the same scalar carrier, while a more localized bump deformation remains too sharp for the same bounded transport closure
- conclusion: all 6 smooth radial cases pass, with maximum radial refinement drift 0.1465; the weakest radial passing case is `radial|n=13|eta=0.050` with kappa/D_eff 0.9019, median relative error 0.0704, p90 error 0.3240, and |metric-tracking corr| 0.9885; by contrast, the sharpest failing bump case is `bump|n=13|eta=0.150` with kappa/D_eff 0.5889, median relative error 0.4299, p90 error 0.7346, and shell-kappa CV 0.4963; this supports a bounded smooth-inhomogeneity closure on the scalar carrier, while keeping localized inhomogeneity explicitly open
