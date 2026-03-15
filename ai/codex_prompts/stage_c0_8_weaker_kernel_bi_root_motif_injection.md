# Stage C0.8 Weaker Kernel and Bi-Root Multi-Layer Motif Injection

HAOS-IIP Pre-Geometry Atlas - Combinatorial Kernel Line

Purpose

Stage C0.7 established the current negative boundary: the no-distance graph-shell branch supports protected ordering under positive mismatch forcing, but explicit kernel-update consistency failures remain zero under adjacency versus graph-shell alternation.

Stage C0.8 tests the minimal next control:

- weaken the nonlocal combinatorial kernel from rigid shell bins to a bounded-path-length weighted kernel
- extend local motif forcing from single-root motifs to bi-root and multi-layer motifs
- keep the protected topology-stability constraint active
- alternate adjacency and the weaker bounded-path kernel once mismatch is live

The question is whether explicit kernel inconsistency can now be forced without importing any coordinate metric or packet-law oracle.

Frozen architecture

Keep fixed from C0.1-C0.7:

- graph size
- representative encounter families
- projected-transverse state sector
- shell-weight measurement stack
- protected topology-stability constraint
- no Euclidean distance
- no Gaussian kernel
- no adaptive memory law

Only change:

- replace the rigid graph-shell survivor kernel with a weaker bounded-path-length kernel
- replace the C0.4 single-root motif ladder with a bi-root and multi-layer motif ladder

Kernel definition

The weaker kernel is purely combinatorial and bounded in hop count.

For hop distance `h`:

`w(h) = 1 / (h + 1)`

with finite support up to a fixed hop cap.

Allowed ingredients:

- adjacency
- bounded hop distance
- degree
- local motif incidence
- shell counts

Not allowed:

- coordinates
- Euclidean distance
- Gaussian falloff
- embedded geometry

Forcing protocols

Run a 9-run matrix:

1. incidence-noise forcing
2. bi-root / multi-layer motif forcing
3. combined noise plus motif forcing

across:

1. clustered composite anchor
2. counter-propagating corridor
3. phase-ordered symmetric triad

Motif ladder

Cycle the bi-root motif family through:

- bi-root triangle cluster
- bi-root diamond
- bi-root hub micro-star
- two-layer stacked motif

Observables

Retain the C0.6-C0.7 protected sequencing stack and add:

- kernel_update_consistency_failure_count
- protected_arrow_step_under_weaker_inconsistency
- bi_root_motif_correlation

Expected honest outcomes

Positive:

- consistency failures become positive
- protected arrow under weaker inconsistency appears in at least 6 of 9 runs
- positive mismatch remains reproducible

This would support a native inconsistency-mediated sequencing mechanism on the no-distance branch.

Negative:

- consistency failures stay zero
- or arrow forcing collapses under the weaker kernel

This would mean the current incidence structure is still too rigid even after softening the kernel family and enriching motif forcing.

Deliverables

- stamped JSON
- stamped CSV
- summary note
- Python runner
- trace and summary plots

Interpretation boundary

Do not claim emergent geometry or physical irreversibility.

At most claim:

the weaker bounded-path combinatorial kernel with bi-root and multi-layer motif forcing does or does not generate explicit kernel-update inconsistency under protected sequencing constraints on the no-distance edge-graph control branch.
