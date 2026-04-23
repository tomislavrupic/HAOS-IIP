# HAOS-IIP Scalar-Carrier Smooth-Inhomogeneity Current Closure 52.5

## Purpose

Freeze the bounded `52.5` follow-on after the positive clean-line shell-native current closure `52.4`.

`52.4` showed that, on the undeformed scalar carrier, the shell-native transient reconstruction closes onto one bounded constitutive family.

The next honest question is stricter:

> if the same carrier is deformed in a controlled way, do the extracted metric field and the shell-native transport law co-deform together?

So `52.5` keeps the carrier, the `51.4` local metric field construction, and the `52.4` shell-native current law fixed.
It changes only the geometry of the point cloud through small deterministic inhomogeneity profiles.

## Construction

The carrier is unchanged:

- same local scalar kernel graph
- same coarse local metric-field extraction
- same shell-native constitutive law

The controlled deformation families are:

- smooth radial deformation
- localized scalar bump

Tested amplitudes:

- `eta = 0.05, 0.10, 0.15`

Clean refinement line:

- `n = 13, 15`

For each case:

1. deform the same clean scalar carrier by the chosen deterministic profile;
2. rebuild the scalar operator on the deformed point cloud;
3. extract the coarse shellwise metric field;
4. fit the shell-native constitutive current law against the transient shell flow;
5. test whether the metric shift and constitutive profile track the imposed deformation in a bounded way.

## Bounded result

The correct `52.5` result is **mixed but positive on the smooth branch**.

All six smooth radial cases pass:

- `radial|n=13|eta=0.050`: `kappa / D_eff = 0.9019`, median relative error `0.0704`, p90 relative error `0.3240`, shell-`kappa` CV `0.1193`, `|corr| = 0.9885`
- `radial|n=13|eta=0.100`: `0.9243`, `0.0665`, `0.3430`, `0.1265`, `0.9939`
- `radial|n=13|eta=0.150`: `0.9466`, `0.0568`, `0.3568`, `0.1310`, `0.9840`
- `radial|n=15|eta=0.050`: `0.9272`, `0.0539`, `0.2429`, `0.1072`, `0.9872`
- `radial|n=15|eta=0.100`: `0.9588`, `0.0496`, `0.2686`, `0.1178`, `0.9935`
- `radial|n=15|eta=0.150`: `0.9848`, `0.0562`, `0.2773`, `0.1215`, `0.9846`

Radial refinement stability also stays bounded:

- normalized shell-`kappa` drift `eta = 0.05`: `0.1202`
- normalized shell-`kappa` drift `eta = 0.10`: `0.1384`
- normalized shell-`kappa` drift `eta = 0.15`: `0.1465`

The localized bump branch stays open:

- `bump|n=13|eta=0.050`: `kappa / D_eff = 0.7857`, median relative error `0.2053`, p90 relative error `0.4658`, shell-`kappa` CV `0.2164`
- `bump|n=13|eta=0.100`: `0.6876`, `0.3113`, `0.6112`, `0.3331`
- `bump|n=13|eta=0.150`: `0.5889`, `0.4299`, `0.7346`, `0.4963`
- `bump|n=15|eta=0.050`: `0.7971`, `0.1569`, `0.3896`, `0.1944`
- `bump|n=15|eta=0.100`: `0.6929`, `0.2487`, `0.5368`, `0.3088`
- `bump|n=15|eta=0.150`: `0.5929`, `0.3565`, `0.6609`, `0.4453`

So the strongest bounded statement now supported is:

> on the same scalar carrier, the extracted local metric field and the shell-native transient constitutive law co-deform together under smooth radial inhomogeneity, while more localized inhomogeneity remains outside the same bounded closure.

## What 52.5 now supports

The scalar line now supports, on one common carrier:

- one bounded local metric field;
- one shell-native inverse-square response closure;
- one bounded clean-line transient current closure;
- one bounded smooth-inhomogeneity co-deformation closure between the metric read and the shell-native transport law.

This is closer to physics than `52.4` alone because the closure is no longer restricted to the undeformed clean line.
It now survives a controlled smooth geometric deformation on the same carrier.

## What 52.5 still does not support

This note still does **not** claim:

- closure for arbitrary or sharply localized inhomogeneity;
- disorder-universal transient transport;
- conserved-current closure;
- curvature extraction;
- spacetime closure;
- ontology or force-law language beyond the bounded scalar carrier read already earned.

So `52.4` remains the clean-line closure, and `52.5` is the bounded smooth-inhomogeneity follow-on.

## Authority and artifacts

- operator note: `experiments/operators/Scalar_Kernel_Graph_Current_Closure_Inhomogeneity_v1.md`
- script: `numerics/simulations/scalar_kernel_graph_current_closure_inhomogeneity.py`
- config key: `config.json -> scalar_kernel_graph_current_closure_inhomogeneity`
- frozen result: `data/20260422_213513_scalar_kernel_graph_current_closure_inhomogeneity.json`
- latest result: `data/scalar_kernel_graph_current_closure_inhomogeneity_latest.json`

## Next honest move

The next scalar step should again stay narrow:

1. test whether the smooth radial branch survives mild carrier disorder on top of the same deformation family;
2. measure the localization threshold where the bump branch stops behaving like the smooth branch;
3. only afterward revisit conserved-current or curvature-style language.
