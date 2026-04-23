# HAOS-IIP Scalar-Carrier Current Closure 52.3

## Purpose

Freeze the first bounded `52.3` step after the positive `52.2` inverse-square response closure.

`52.2` answered the static question:

- on the validated scalar carrier,
- with the `51.4` local metric field,
- a shell-native, law-aware reconstruction closes one inverse-square-like response family.

`52.3` asks the stronger transient question:

> does the same scalar carrier also support a constitutive response/current closure?

## Construction

The carrier stays fixed:

- same local scalar kernel graph
- same coarse local metric field
- same clean Euclidean shell scaffold

For each clean refinement `n = 13, 15, 17`:

1. evolve the point-source heat trajectory,
2. infer empirical shell current from cumulative shell-mass flow,
3. compare it to a constitutive shell current built from the same metric field and shell gradient,
4. fit one constitutive coefficient `kappa` per refinement on the bounded time window.

So this note is not asking whether a static inverse-square law exists. It is asking whether the transient flow itself closes onto one current law on the same carrier.

## Bounded result

The first `52.3` current-closure probe is **open**.

Representative clean refinement read:

- `n = 13`: median relative error `0.5925`, `p90 = 0.6671`
- `n = 15`: median relative error `0.5090`, `p90 = 0.5997`
- `n = 17`: median relative error `0.4850`, `p90 = 0.6198`

Normalized current-profile drift also remains above the bounded pass target:

- `13 -> 15`: `0.1874`
- `15 -> 17`: `0.2073`

So the strongest honest statement is:

> the scalar carrier now supports the static `52.2` inverse-square response law, but transient constitutive current closure is not yet closed even on the clean refinement scan.

## Meaning

This does **not** undo the positive scalar line.

What remains positive:

- scalar operator control
- shared-family scalar geometry bridge
- scalar-carrier geometry robustness
- local metric-field extraction
- shell-native inverse-square recoverability-gradient closure

What remains open:

- a transient current law that uses the same geometry and closes with bounded residual under refinement

So the present bottleneck is no longer the scalar carrier itself.
It is the transient shell-current reconstruction and constitutive closure.

## Next honest move

The next step should be smaller, not broader:

1. improve transient shell-current reconstruction on the same clean carrier,
2. understand the scaling of the fitted constitutive coefficient,
3. only then widen toward disorder or kernel-family robustness.
