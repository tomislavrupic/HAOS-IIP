# Periodic / Twisted L1 Hodge Projection

## Purpose

Separate the periodic/twisted edge space explicitly into

- exact sector: `im d0`
- harmonic sector: `H1 = ker d0* intersect ker d1`
- coexact sector: `im d1*`

and test whether the restricted operator on `im d1*` develops a stable low band after harmonic torus cycles are removed.

## Setup

- substrate: periodic cubic lattice (`3`-torus)
- kernel weight: `w_e = exp(-h^2 / (2 epsilon))`
- `epsilon = 0.2`
- sizes: `[4, 5]`
- flux scan: `[0, 1, 2, 3, 4]`
- low modes analyzed: `6`

## Hodge Projection Construction

The decomposition is built from the untwisted periodic reference complex with the same lattice size and edge weights, then applied to the twisted operator `L1`. This keeps the subspace split topological while the operator carries the flux dependence:

- `im d0` from the image of the untwisted node-edge differential
- `H1` from the simultaneous nullspace of untwisted `d0*` and `d1`
- `im d1*` from the image of the untwisted face-edge adjoint

The coexact operator of interest is

$$
L_{1,\mathrm{coexact}} = P_{\mathrm{coexact}} L_1 P_{\mathrm{coexact}}.
$$

## Decomposition Verification

Reference untwisted periodic complex: `periodic_n5_flux1`

- exact dimension: `124`
- harmonic dimension: `3`
- coexact dimension: `248`
- chain error `||d1 d0||`: `2.175e-16`
- orthogonality error `exact-harmonic`: `8.535e-16`
- orthogonality error `exact-coexact`: `2.733e-15`
- orthogonality error `harmonic-coexact`: `1.496e-15`
- projector reconstruction error: `5.319e-15`
- sample reconstruction error: `1.888e-15`

These checks stay small across the scanned cases, so the explicit Hodge splitting is numerically stable in the tested range.

## Low Coexact Branch

### Size `n = 4`

| `m` | `phi` | first coexact `lambda` | first-mode coexact fraction | harmonic dim |
| --- | ---: | ---: | ---: | ---: |
| 0 | 0.000000 | 1.710691 | 1.0000 | 3 |
| 1 | 0.392699 | 0.679582 | 1.0000 | 3 |
| 2 | 0.785398 | 0.606733 | 1.0000 | 3 |
| 3 | 1.178097 | 0.771051 | 1.0000 | 3 |
| 4 | 1.570796 | 0.651564 | 1.0000 | 3 |

### Size `n = 5`

| `m` | `phi` | first coexact `lambda` | first-mode coexact fraction | harmonic dim |
| --- | ---: | ---: | ---: | ---: |
| 0 | 0.000000 | 1.250455 | 1.0000 | 3 |
| 1 | 0.251327 | 0.458457 | 1.0000 | 3 |
| 2 | 0.502655 | 0.426759 | 1.0000 | 3 |
| 3 | 0.753982 | 0.619483 | 1.0000 | 3 |
| 4 | 1.005310 | 0.657667 | 1.0000 | 3 |

Representative comparison at nonzero flux:

- `n = 5`, `m = 1`: full lowest `lambda = 0.013057`, harmonic-projected `lambda = 0.509808`, coexact-projected `lambda = 0.458457`
- `n = 4`, `m = 1`: full lowest `lambda = 0.014110`, harmonic-projected `lambda = 0.611897`, coexact-projected `lambda = 0.679582`

These comparisons are the core Maxwell test in the current range: once the harmonic torus cycles are removed, the coexact branch remains at moderate eigenvalue rather than forming the actual bottom band of the spectrum.

## Direct Answers

### Does a low coexact family persist?

Only as a projected sector. It survives algebraically, but not as the actual low branch once harmonic modes are removed.

### Does it respond smoothly to flux?

Only partially. The larger lattice is smoother than the smaller one.

### Does it scale consistently with lattice size?

Yes, at least across the tested sizes.

### Does it still look harmonic/topological, or more like a lattice vector band?

no clear coexact low band yet; the very low modes remain dominated by harmonic or mixed topological remnants

## Current verdict

- what worked: explicit projection onto `im d0`, `H1`, and `im d1*` separated the harmonic torus cycles from the coexact branch without breaking the current periodic/twisted operator path
- what did not appear: no clean continuum-like Maxwell band or fully smooth low coexact branch across the whole scan
- what must be tried next: larger periodic sizes and defect/puncture backgrounds to test whether the coexact branch localizes into genuine circulation modes or settles into a lattice band

## Artifacts

- results: `data/20260310_000308_periodic_twisted_l1_hodge_projection.json`
- plots: `plots/20260310_000308_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_full_vs_coexact.png`, `plots/20260310_000308_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_harmonic_fraction.png`, `plots/20260310_000308_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_coexact_fraction.png`, `plots/20260310_000308_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_coexact_flow.png`, `plots/20260310_000308_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_persistence.png`, `plots/20260310_000308_periodic_twisted_l1_hodge_projection_periodic_twisted_l1_hodge_projection_modes.png`
- timestamp: `20260310_000308`
