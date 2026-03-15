# Stage C0.7 Kernel Inconsistency Forcing v1

Timestamped summary: `data/20260315_205157_stage_c0_7_kernel_inconsistency_forcing.json`
Timestamped run table: `data/20260315_205157_stage_c0_7_kernel_inconsistency_forcing.csv`

Architecture notice: This branch keeps graph size, adjacency and graph-shell kernels, the clustered / corridor / triad representative families, the projected-transverse sector, and the shell-weight measurement vocabulary fixed. Only the post-mismatch update rule changes: survivor updates are deliberately alternated across adjacency and graph-shell kernels.
Baseline notice: C0.7 keeps the corrected C0.6 forcing schedules and the protected topology-stability rule active, but adds an explicit composition-level inconsistency probe. The key quantity is the disagreement between adjacency-to-shell and shell-to-adjacency survivor refinements on the same forced graph.

Stage summary: explicit consistency failures in `0/9` runs; protected arrow under inconsistency in `0/9` runs.

Per-run summary:

- `Noise-seeded inconsistency forcing` / `clustered_composite_anchor`
  - protected topology label: `trapped_local`
  - sequencing class: `ordered_without_inconsistency`
  - final topology label: `trapped_local`
  - consistency failure total: `0`
  - consistency failure max: `0`
  - protected arrow step under inconsistency: `-1`
  - positive mismatch density: `0.0802`
- `Motif-seeded inconsistency forcing` / `clustered_composite_anchor`
  - protected topology label: `trapped_local`
  - sequencing class: `ordered_without_inconsistency`
  - final topology label: `trapped_local`
  - consistency failure total: `0`
  - consistency failure max: `0`
  - protected arrow step under inconsistency: `-1`
  - positive mismatch density: `0.0064`
- `Hybrid seeded inconsistency forcing` / `clustered_composite_anchor`
  - protected topology label: `trapped_local`
  - sequencing class: `ordered_without_inconsistency`
  - final topology label: `trapped_local`
  - consistency failure total: `0`
  - consistency failure max: `0`
  - protected arrow step under inconsistency: `-1`
  - positive mismatch density: `0.0802`
- `Noise-seeded inconsistency forcing` / `counter_propagating_corridor`
  - protected topology label: `smeared_transfer`
  - sequencing class: `rigid_positive_mismatch`
  - final topology label: `smeared_transfer`
  - consistency failure total: `0`
  - consistency failure max: `0`
  - protected arrow step under inconsistency: `-1`
  - positive mismatch density: `0.0816`
- `Motif-seeded inconsistency forcing` / `counter_propagating_corridor`
  - protected topology label: `smeared_transfer`
  - sequencing class: `rigid_positive_mismatch`
  - final topology label: `smeared_transfer`
  - consistency failure total: `0`
  - consistency failure max: `0`
  - protected arrow step under inconsistency: `-1`
  - positive mismatch density: `0.0158`
- `Hybrid seeded inconsistency forcing` / `counter_propagating_corridor`
  - protected topology label: `smeared_transfer`
  - sequencing class: `rigid_positive_mismatch`
  - final topology label: `smeared_transfer`
  - consistency failure total: `0`
  - consistency failure max: `0`
  - protected arrow step under inconsistency: `-1`
  - positive mismatch density: `0.0816`
- `Noise-seeded inconsistency forcing` / `phase_ordered_symmetric_triad`
  - protected topology label: `split_channel_exchange`
  - sequencing class: `ordered_without_inconsistency`
  - final topology label: `split_channel_exchange`
  - consistency failure total: `0`
  - consistency failure max: `0`
  - protected arrow step under inconsistency: `-1`
  - positive mismatch density: `0.0802`
- `Motif-seeded inconsistency forcing` / `phase_ordered_symmetric_triad`
  - protected topology label: `split_channel_exchange`
  - sequencing class: `ordered_without_inconsistency`
  - final topology label: `split_channel_exchange`
  - consistency failure total: `0`
  - consistency failure max: `0`
  - protected arrow step under inconsistency: `-1`
  - positive mismatch density: `0.0041`
- `Hybrid seeded inconsistency forcing` / `phase_ordered_symmetric_triad`
  - protected topology label: `split_channel_exchange`
  - sequencing class: `ordered_without_inconsistency`
  - final topology label: `split_channel_exchange`
  - consistency failure total: `0`
  - consistency failure max: `0`
  - protected arrow step under inconsistency: `-1`
  - positive mismatch density: `0.0802`

Interpretation boundary:
- this scan asks whether explicit adjacency to graph-shell composition failure can coexist with protected arrow forcing on the no-distance branch
- topology labels remain persistence-class observables and the protected topology-stability rule stays active
- if consistency failures become positive and the protected arrow still fires broadly, the branch supports an irreversibility mechanism tied directly to kernel-update inconsistency

Series boundary:

“The C0 line (C0.1–C0.7) demonstrates reproducible incidence-driven mismatch accumulation and protected ordering on the pure graph-shell branch, but kernel-update consistency failures remain zero under forced alternation. The branch therefore supports incidence-forced protected sequencing while remaining too rigid for explicit kernel-mediated irreversibility. This sets the precise boundary before any Phase II back-reaction claim.”
