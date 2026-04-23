# Scalar Kernel Graph Current Closure Radial Disorder v1

## Purpose

Test the next honest robustness question after the bounded smooth-inhomogeneity closure `52.5`:

> if the passing smooth radial branch is perturbed by mild carrier disorder, does the shell-native transient closure remain bounded?

The carrier stays fixed:

- same scalar kernel graph family
- same smooth radial deformation family
- same `52.5` transport thresholds

The only new ingredient is mild random carrier disorder added on top of the already-passing smooth radial deformation.

## Construction

- clean refinements: `[13, 15]`
- radial amplitudes: `[0.05, 0.1, 0.15]`
- disorder fractions of lattice spacing: `[0.02, 0.04]`
- seeds: `[0, 1, 2]`

This probe uses a scaffold-shell read:

1. build the clean shell layout once on the undeformed cubic scaffold;
2. apply smooth radial deformation;
3. add mild random interior disorder while keeping the boundary fixed;
4. read shell transport on the clean scaffold-shell assignment so shell bookkeeping itself does not collapse under jitter;
5. test whether the transient constitutive law remains within the same bounded window.

Artifacts:

- script: `numerics/simulations/scalar_kernel_graph_current_closure_radial_disorder.py`
- config: `config.json -> scalar_kernel_graph_current_closure_radial_disorder`
- results: `data/20260422_214620_scalar_kernel_graph_current_closure_radial_disorder.json`
- latest: `data/scalar_kernel_graph_current_closure_radial_disorder_latest.json`
- plots:
  - `plots/20260422_214620_scalar_kernel_graph_current_closure_radial_disorder_summary.png`
  - `plots/20260422_214620_scalar_kernel_graph_current_closure_radial_disorder_profiles.png`

## Aggregated cases

| case | trial passes | median rel. error | p90 rel. error | diffusivity gap | shell `kappa` CV | `|metric corr|` | status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `radial_disorder|n=13|eta=0.050|sigma=0.020` | `0/3` | `1.3652` | `2.6666` | `1.2110` | `1.0198` | `0.9890` | OPEN |
| `radial_disorder|n=15|eta=0.050|sigma=0.020` | `0/3` | `1.0221` | `1.0610` | `1.0146` | `0.7471` | `0.9854` | OPEN |
| `radial_disorder|n=13|eta=0.100|sigma=0.020` | `0/3` | `1.4294` | `2.7874` | `1.2553` | `1.0839` | `0.9900` | OPEN |
| `radial_disorder|n=15|eta=0.100|sigma=0.020` | `0/3` | `1.0628` | `1.1929` | `1.0437` | `0.8801` | `0.9859` | OPEN |
| `radial_disorder|n=13|eta=0.150|sigma=0.020` | `0/3` | `1.4827` | `3.3771` | `1.3027` | `1.0805` | `0.9850` | OPEN |
| `radial_disorder|n=15|eta=0.150|sigma=0.020` | `0/3` | `1.1029` | `1.3601` | `1.0744` | `0.9520` | `0.9807` | OPEN |
| `radial_disorder|n=13|eta=0.050|sigma=0.040` | `0/3` | `1.3713` | `2.5173` | `1.2131` | `1.0367` | `0.9826` | OPEN |
| `radial_disorder|n=15|eta=0.050|sigma=0.040` | `0/3` | `1.0337` | `1.0934` | `1.0224` | `0.7645` | `0.9813` | OPEN |
| `radial_disorder|n=13|eta=0.100|sigma=0.040` | `0/3` | `1.4345` | `2.9166` | `1.2564` | `1.0834` | `0.9869` | OPEN |
| `radial_disorder|n=15|eta=0.100|sigma=0.040` | `0/3` | `1.0755` | `1.2449` | `1.0522` | `0.9006` | `0.9840` | OPEN |
| `radial_disorder|n=13|eta=0.150|sigma=0.040` | `0/3` | `1.4844` | `3.2989` | `1.3041` | `1.0799` | `0.9858` | OPEN |
| `radial_disorder|n=15|eta=0.150|sigma=0.040` | `0/3` | `1.1151` | `1.4192` | `1.0826` | `0.9583` | `0.9819` | OPEN |

## Refinement drift

| case | drift |
| --- | --- |
| `sigma=0.020|eta=0.050` | `0.8398` |
| `sigma=0.020|eta=0.100` | `0.9359` |
| `sigma=0.020|eta=0.150` | `0.9611` |
| `sigma=0.040|eta=0.050` | `0.8578` |
| `sigma=0.040|eta=0.100` | `0.9446` |
| `sigma=0.040|eta=0.150` | `0.9596` |

## Interpretation

- observation: the scaffold-shell metric-tracking signal survives mild disorder on the smooth radial branch, but the shell-native transient constitutive closure does not remain bounded under the same disorder window
- conclusion: all 12 aggregated smooth-radial disorder cases remain open; the weakest failure is `radial_disorder|n=15|eta=0.050|sigma=0.020` with median relative error 1.0221, p90 error 1.0610, diffusivity gap 1.0146, shell-kappa CV 0.7471, and |metric-tracking corr| 0.9854; maximum refinement drift is 0.9611; this keeps disorder-robust transient closure explicitly open and points to a disorder-native flux reconstruction as the next honest step
