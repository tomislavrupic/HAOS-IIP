# HAOS-IIP Scalar-Carrier Recoverability Gradient 52.2

## Purpose

Freeze the bounded `52.2` follow-on after the open `52.1` probe.

The carrier is unchanged:

- same scalar kernel-graph family
- same `51.4` local metric field
- same bounded disorder and kernel-family window

The only change is the response reconstruction.

Instead of a raw local least-squares derivative, `52.2` uses a shell-native, law-aware reconstruction aligned with the already-validated scalar Green law

`phi(r) ≈ A + B/r`.

## What changed from 52.1

`52.1` already showed that the response field is directionally coherent:

- radial alignment stayed around `0.99+`
- the carrier and local metric were not the weak point

What stayed open in `52.1` was the law extraction itself:

- shell slope
- fit quality
- scaled-flux constancy
- refinement-profile stability

`52.2` therefore keeps the carrier fixed and changes only the reconstruction:

`F_r(r) ≈ g_rr(r) |B| / r^2`

with:

- `B` from the validated Green `A + B/r` fit
- `g_rr(r)` from the shell-averaged radial projection of the `51.4` coarse local metric field

## Bounded result

The `52.2` shell-native reconstruction is **positive** on the tested window.

Passed window:

- all `6/6` mild-disorder cases
- all `6/6` bounded kernel-family cases

Worst passing case:

- `disorder|n=13|jitter=0.040`
- Green fit `R^2 = 0.9856`
- power slope `-1.9966`
- power fit `R^2 = 0.9998`
- scaled-response CV `0.0067`

Refinement stability:

- maximum normalized profile drift `0.0089`

So the correct bounded statement is now:

> on the validated scalar carrier and the `51.4` local metric field, the recoverability-gradient response closes as one inverse-square-like family under a shell-native, law-aware reconstruction aligned with the already-validated Green structure.

## What this does and does not mean

What it supports:

- the `52.1` failure was a reconstruction failure, not a carrier failure
- the scalar carrier now supports a bounded force-like inverse-square closure in a law-aware shell-native reading
- the same operator family still organizes the potential, the metric field, and the response law on one common carrier

What it does **not** support:

- a universal force law
- a raw local vector-force closure independent of reconstruction choice
- arbitrary-substrate or strong-disorder universality
- current, curvature, or spacetime closure

## Clean next move

The next honest step is no longer “find an inverse-square law at all.”

It is to ask whether the same scalar carrier can now support:

1. response/current closure on the same reconstructed geometry,
2. controlled inhomogeneity tests where the extracted metric and the response law track designed deformation together,
3. then only afterward curvature-like and conserved-current questions.
