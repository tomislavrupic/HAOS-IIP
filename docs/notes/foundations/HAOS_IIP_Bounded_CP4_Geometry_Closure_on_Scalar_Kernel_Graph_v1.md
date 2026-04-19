# HAOS-IIP Bounded CP4 Geometry Closure on the Scalar Kernel-Graph v1

## Status

- bounded geometry-closure note
- scope: scalar kernel-graph carrier only
- authority boundary: positive on the tested local cubic carrier, not yet a repo-wide universality or curvature claim

## Purpose

Freeze the first clean positive `CP4` result after the geometry-closure preflight blocker was identified.

The narrow question is:

> on one common scalar kernel-graph family, under one common metric-like coarse reconstruction, do Green response, heat behavior, shell-arrival structure, and low-mode organization support one coherent effective geometry read?

## Preflight blocker

The preflight note [HAOS_IIP_Geometry_Closure_Preflight_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/notes/foundations/HAOS_IIP_Geometry_Closure_Preflight_v1.md) showed that the previous strongest ingredients were not commensurate:

- Green response lived on a `3D cubic kernel-graph` line
- heat, low-mode, and shell-arrival surrogates lived on the old `2D periodic branch-local torus` line

So the original repo state could not honestly support a geometry-closure claim.

## Executed shared-family rebuild

That mismatch is now removed by the shared scalar-kernel-graph bridge:

- script: [scalar_kernel_graph_geometry_bridge.py](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/numerics/simulations/scalar_kernel_graph_geometry_bridge.py)
- experiment note: [Scalar_Kernel_Graph_Geometry_Bridge_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/experiments/operators/Scalar_Kernel_Graph_Geometry_Bridge_v1.md)
- stamped artifact: [20260419_234444_scalar_kernel_graph_geometry_bridge.json](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/data/20260419_234444_scalar_kernel_graph_geometry_bridge.json)

The rebuild keeps one carrier fixed:

- cubic `3D` scalar kernel graph
- local kernel regime `c_epsilon = 0.5`
- same induced scalar operator `\widehat{L}_h`
- same Euclidean shell reconstruction centered at the source

All four channels now live on that one carrier:

1. Green field fit `A + B/r`
2. shell second moment diffusion `\langle r^2 \rangle(t)`
3. shell-arrival times against `r^2`
4. first positive low-mode triplet against the coordinate subspace `span(x,y,z)`

## Quantitative read

On the tested `n = 9, 11, 13` window:

- Green effective-dimension hint: `2.9758, 2.9673, 2.9663`
- Green fit quality: `R^2 = 0.9769, 0.9834, 0.9881`
- effective diffusivity from shell second moment: `0.9551, 0.9594, 0.9616`
- shell-arrival fit quality against `r^2`: `R^2 = 0.9886, 0.9779, 0.9787`
- minimum first-shell / coordinate principal cosine: `0.9937, 0.9934, 0.9932`

That is the first bounded repo state in which all four geometry-facing channels point in the same direction on one common operator/substrate family.

## Bounded conclusion

The correct positive statement is:

> `CP4` now closes in bounded scalar-carrier form: on the tested local cubic scalar kernel graph, one metric-like reconstruction organizes Green response, heat behavior, shell-arrival structure, and low-mode organization.

The correct limits are just as important:

- this is **not** yet a full geometry closure for every active line in the repository
- this is **not** yet a universality statement across kernel families, point-set disorder, or alternative coarse reconstructions
- this is **not** yet a curvature, conserved-current, or spacetime claim

## What is now closed

The following narrow failure mode is closed:

- the repo can no longer be accused of stitching incompatible `3D` and `2D` substrates together to fake a geometry read

The following broader questions remain open:

- does the same geometry-style organization survive mild substrate disorder on the scalar carrier?
- does it survive bounded kernel-family variation?
- can metric-like structure, response, and downstream geometry-like quantities be extracted from the same carrier without changing the coarse reconstruction rule?

## Practical next move

The next honest step after this note is not another declaration. It is controlled robustness testing on the same scalar carrier:

1. mild point-set disorder on the same geometry bridge
2. bounded kernel-family variation
3. response / metric extraction from the same coarse reconstruction

Only if those survive should the repo say more than the bounded carrier-level `CP4` statement above.
