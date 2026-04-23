# HAOS-IIP Scalar-Carrier Current Closure 52.4

## Purpose

Freeze the bounded `52.4` follow-on after the open `52.3` transient current probe.

`52.3` kept the carrier fixed and asked the right stronger question:

> does the same scalar carrier support a constitutive response/current closure?

Its answer was still open. The main failure was not the carrier or the `51.4` local metric field. It was the transient reconstruction:

- empirical current was read from cumulative shell-mass flow,
- constitutive current was read from a shell profile that still mixed shell-center and shell-boundary viewpoints,
- the fitted `kappa` scaled incorrectly with refinement.

So `52.4` keeps the carrier fixed and changes only the reconstruction.

## What changed from 52.3

The carrier is unchanged:

- same local scalar kernel graph
- same `51.4` coarse local metric field
- same clean refinement line `n = 13, 15, 17`

The reconstruction is tighter:

1. use exact bulk shells on the clean cubic scaffold;
2. keep empirical current as cumulative shell-mass flow;
3. convert shell mass to shell density

`rho(r,t) = m_shell(r,t) / (count_shell(r) h^3)`;

4. compare it against the shell-native constitutive law

`I_const(r,t) = - kappa 4 pi r^2 g_rr(r) partial_r rho(r,t)`.

So `52.4` is not a new substrate result. It is a reconstruction audit on the same carrier.

## Bounded result

The `52.4` shell-native transient current reconstruction is **positive** on the clean refinement line.

Clean passing cases:

- `clean|n=13`: `kappa / D_eff = 0.8871`, median relative error `0.0705`, p90 relative error `0.2328`, shell-`kappa` CV `0.0967`
- `clean|n=15`: `kappa / D_eff = 0.9106`, median relative error `0.0709`, p90 relative error `0.1984`, shell-`kappa` CV `0.0742`
- `clean|n=17`: `kappa / D_eff = 0.9263`, median relative error `0.0739`, p90 relative error `0.2718`, shell-`kappa` CV `0.1176`

Refinement stability:

- normalized shell-`kappa` drift `13 -> 15`: `0.0493`
- normalized shell-`kappa` drift `15 -> 17`: `0.0966`

So the strongest bounded statement is now:

> on the clean scalar carrier, a shell-native transient current reconstruction closes onto one bounded constitutive family whose fitted conductivity stays close to the heat diffusivity and whose shellwise profile remains refinement-stable.

## What 52.4 now supports

This strengthens the scalar line in a specifically physics-facing way.

The same carrier now supports:

- one shared-family scalar geometry bridge;
- one bounded local metric field;
- one shell-native inverse-square recoverability-gradient law;
- one bounded shell-native transient current-closure read on the clean refinement line.

This means the scalar carrier no longer supports only static geometry and static response language.
It now supports one bounded transport-law closure on the same geometry carrier.

## What 52.4 does not support

This note still does **not** claim:

- mild-disorder or kernel-family current-closure universality;
- a conserved-current theorem;
- curvature extraction;
- spacetime closure;
- ontology or fundamental-force language.

`52.3` therefore remains the honest open baseline for the first transient probe, while `52.4` is the corrected shell-native closure on the same clean carrier.

## Next honest move

The next step should again stay narrow:

1. test whether the same shell-native current closure survives controlled inhomogeneity on the same carrier;
2. then test mild disorder and bounded kernel-family variation for the transient law itself;
3. only afterward ask for conserved-current or curvature-style closure.
