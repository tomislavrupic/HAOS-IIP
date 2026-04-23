# HAOS-IIP Scalar-Carrier Recoverability Gradient 52.1

## Purpose

Freeze the first bounded `52`-line probe after `51.4`.

The scalar geometry carrier is now strong enough to support a first force-like response test. The question is whether the same carrier and the same local metric field already support one stable inverse-square-like **recoverability-gradient** family.

## Deterministic proxy used here

The first probe stays fully inside the existing scalar carrier:

- recoverability potential: the same mean-zero scalar Green field
- local metric: the `51.4` coarse local metric field
- effective response field:

`J_rec = - G ∇φ`

This is intentionally weaker than a claim about fundamental force. It is only the first operator-native response proxy on the validated scalar carrier.

## Bounded result

The first probe is **open**.

What already works:

- the response field is strongly radial on all tested cases
- radial alignment stays around `0.99+`

What does **not** close yet:

- inverse-square shell slope
- shell-profile fit quality
- scaled-flux constancy
- refinement-profile stability

Representative clean values:

- `n=11`: slope `-0.8953`, fit `R^2 = 0.5966`, flux CV `0.3139`
- `n=13`: slope `-1.2780`, fit `R^2 = 0.8090`, flux CV `0.2312`

So the carrier and local metric are already good enough to give a directionally coherent response field, but not yet one stable inverse-square family.

## What 52.1 now supports

The repository can now say:

> the scalar Green potential supports a coherent radial response field on the validated scalar carrier, but the first recoverability-gradient closure does not yet support a positive inverse-square law claim.

## What the failure means

This is not a failure of the scalar carrier itself.

It means the first local differentiation and shell-flux reconstruction are still too crude to close the response law, even though:

- the Green potential itself is already geometry-compatible;
- the local metric field is already stable;
- the response direction is already strongly radial.

## Practical next move

The next honest step is `52.2`:

- improve the response reconstruction, not the carrier;
- test shell-native or law-aware gradient extraction;
- only then ask again whether one stable inverse-square family closes.
