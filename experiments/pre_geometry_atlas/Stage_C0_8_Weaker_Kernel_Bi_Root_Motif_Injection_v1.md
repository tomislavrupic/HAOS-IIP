# Stage C0.8 Weaker Kernel Bi-Root Motif Injection v1

Timestamped summary: `data/20260315_210832_stage_c0_8_weaker_kernel_bi_root_motif_injection.json`
Timestamped run table: `data/20260315_210832_stage_c0_8_weaker_kernel_bi_root_motif_injection.csv`

Architecture notice: This branch keeps graph size, the clustered / corridor / triad representative families, the projected-transverse sector, and the shell-weight measurement vocabulary fixed. Only the survivor kernel and motif forcing family change: the graph-shell kernel is softened to a bounded-path combinatorial kernel and single-root motifs are replaced by bi-root and two-layer injections.
Baseline notice: C0.8 keeps the corrected C0.6 positive-mismatch forcing logic and the protected topology-stability rule active, but replaces the graph-shell survivor map with a bounded-path map using hop weights 1 / (h + 1) up to a fixed hop cap, then alternates that weaker kernel against adjacency once mismatch is live.

Stage summary: explicit consistency failures in `5/9` runs; protected arrow under weaker inconsistency in `2/9` runs.

Per-run summary:

- `Noise-seeded weaker-kernel inconsistency` / `clustered_composite_anchor`
  - protected topology label: `trapped_local`
  - sequencing class: `weaker_inconsistency_stable_without_arrow`
  - final topology label: `trapped_local`
  - consistency failure total: `2`
  - consistency failure max: `2`
  - protected arrow step under weaker inconsistency: `-1`
  - positive mismatch density: `0.0802`
  - bi-root motif correlation: `0.0000`
- `Bi-root motif-seeded weaker-kernel inconsistency` / `clustered_composite_anchor`
  - protected topology label: `trapped_local`
  - sequencing class: `ordered_without_weaker_inconsistency`
  - final topology label: `trapped_local`
  - consistency failure total: `0`
  - consistency failure max: `0`
  - protected arrow step under weaker inconsistency: `-1`
  - positive mismatch density: `0.0068`
  - bi-root motif correlation: `0.4582`
- `Hybrid noise plus motif weaker-kernel inconsistency` / `clustered_composite_anchor`
  - protected topology label: `trapped_local`
  - sequencing class: `weaker_inconsistency_ordered_persistent`
  - final topology label: `trapped_local`
  - consistency failure total: `2`
  - consistency failure max: `2`
  - protected arrow step under weaker inconsistency: `7`
  - positive mismatch density: `0.0864`
  - bi-root motif correlation: `0.2505`
- `Noise-seeded weaker-kernel inconsistency` / `counter_propagating_corridor`
  - protected topology label: `smeared_transfer`
  - sequencing class: `weaker_inconsistency_stable_without_arrow`
  - final topology label: `smeared_transfer`
  - consistency failure total: `8`
  - consistency failure max: `4`
  - protected arrow step under weaker inconsistency: `-1`
  - positive mismatch density: `0.0816`
  - bi-root motif correlation: `0.0000`
- `Bi-root motif-seeded weaker-kernel inconsistency` / `counter_propagating_corridor`
  - protected topology label: `smeared_transfer`
  - sequencing class: `ordered_without_weaker_inconsistency`
  - final topology label: `smeared_transfer`
  - consistency failure total: `0`
  - consistency failure max: `0`
  - protected arrow step under weaker inconsistency: `-1`
  - positive mismatch density: `0.0160`
  - bi-root motif correlation: `0.8163`
- `Hybrid noise plus motif weaker-kernel inconsistency` / `counter_propagating_corridor`
  - protected topology label: `smeared_transfer`
  - sequencing class: `weaker_inconsistency_stable_without_arrow`
  - final topology label: `smeared_transfer`
  - consistency failure total: `2`
  - consistency failure max: `2`
  - protected arrow step under weaker inconsistency: `-1`
  - positive mismatch density: `0.0828`
  - bi-root motif correlation: `0.2761`
- `Noise-seeded weaker-kernel inconsistency` / `phase_ordered_symmetric_triad`
  - protected topology label: `split_channel_exchange`
  - sequencing class: `ordered_without_weaker_inconsistency`
  - final topology label: `split_channel_exchange`
  - consistency failure total: `0`
  - consistency failure max: `0`
  - protected arrow step under weaker inconsistency: `-1`
  - positive mismatch density: `0.0802`
  - bi-root motif correlation: `0.0000`
- `Bi-root motif-seeded weaker-kernel inconsistency` / `phase_ordered_symmetric_triad`
  - protected topology label: `split_channel_exchange`
  - sequencing class: `ordered_without_weaker_inconsistency`
  - final topology label: `split_channel_exchange`
  - consistency failure total: `0`
  - consistency failure max: `0`
  - protected arrow step under weaker inconsistency: `-1`
  - positive mismatch density: `0.0041`
  - bi-root motif correlation: `0.5767`
- `Hybrid noise plus motif weaker-kernel inconsistency` / `phase_ordered_symmetric_triad`
  - protected topology label: `split_channel_exchange`
  - sequencing class: `weaker_inconsistency_ordered_persistent`
  - final topology label: `split_channel_exchange`
  - consistency failure total: `8`
  - consistency failure max: `4`
  - protected arrow step under weaker inconsistency: `2`
  - positive mismatch density: `0.0853`
  - bi-root motif correlation: `0.2485`

Interpretation boundary:
- this scan asks whether a softer bounded-path combinatorial kernel can open explicit inconsistency under the same protected sequencing rule that stayed rigid in C0.7
- topology labels remain persistence-class observables and the protected topology-stability rule stays active
- if consistency failures become positive and the protected arrow still fires broadly, the branch supports an irreversibility mechanism tied directly to adjacency versus weaker-kernel inconsistency