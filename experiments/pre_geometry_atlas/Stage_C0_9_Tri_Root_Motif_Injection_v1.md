# Stage C0.9 Tri-Root Motif Injection v1

Timestamped summary: `data/20260315_212227_stage_c0_9_tri_root_motif_injection.json`
Timestamped run table: `data/20260315_212227_stage_c0_9_tri_root_motif_injection.csv`

Architecture notice: This branch keeps graph size, the clustered / corridor / triad representative families, the projected-transverse sector, and the shell-weight measurement vocabulary fixed. Only the forcing family changes: the weaker bounded-path survivor kernel from C0.8 is retained, while the motif ladder is upgraded from bi-root and two-layer motifs to tri-root and three-layer constructions.
Baseline notice: C0.9 inherits the corrected C0.6 positive-mismatch forcing logic and the C0.8 bounded-path kernel, keeps the protected topology-stability rule active, and asks whether tri-root / 3-layer motifs can scale the partial C0.8 inconsistency signal into broad protected irreversibility under adjacency-versus-bounded-path alternation.

Stage summary: explicit consistency failures in `4/9` runs; tri-root-specific failure forcing in `2/9` runs; protected arrow under weaker inconsistency in `0/9` runs; protected arrow under tri-root inconsistency in `0/9` runs; irreversibility scaling factor `0.00` relative to the C0.8 `2/9` baseline.

Per-run summary:

- `Noise-seeded tri-root inconsistency` / `clustered_composite_anchor`
  - protected topology label: `trapped_local`
  - sequencing class: `weaker_inconsistency_stable_without_arrow`
  - final topology label: `trapped_local`
  - consistency failure total: `2`
  - tri-root failure total: `0`
  - protected arrow step under weaker inconsistency: `-1`
  - protected arrow step under tri-root inconsistency: `-1`
  - positive mismatch density: `0.0802`
  - tri-root motif correlation: `0.0000`
- `Tri-root motif-seeded inconsistency` / `clustered_composite_anchor`
  - protected topology label: `trapped_local`
  - sequencing class: `ordered_without_tri_root_inconsistency`
  - final topology label: `trapped_local`
  - consistency failure total: `0`
  - tri-root failure total: `0`
  - protected arrow step under weaker inconsistency: `-1`
  - protected arrow step under tri-root inconsistency: `-1`
  - positive mismatch density: `0.0027`
  - tri-root motif correlation: `0.6321`
- `Hybrid noise plus tri-root inconsistency` / `clustered_composite_anchor`
  - protected topology label: `trapped_local`
  - sequencing class: `rigid_positive_mismatch`
  - final topology label: `trapped_local`
  - consistency failure total: `0`
  - tri-root failure total: `0`
  - protected arrow step under weaker inconsistency: `-1`
  - protected arrow step under tri-root inconsistency: `-1`
  - positive mismatch density: `0.0820`
  - tri-root motif correlation: `0.2761`
- `Noise-seeded tri-root inconsistency` / `counter_propagating_corridor`
  - protected topology label: `smeared_transfer`
  - sequencing class: `weaker_inconsistency_stable_without_arrow`
  - final topology label: `smeared_transfer`
  - consistency failure total: `8`
  - tri-root failure total: `0`
  - protected arrow step under weaker inconsistency: `-1`
  - protected arrow step under tri-root inconsistency: `-1`
  - positive mismatch density: `0.0816`
  - tri-root motif correlation: `0.0000`
- `Tri-root motif-seeded inconsistency` / `counter_propagating_corridor`
  - protected topology label: `smeared_transfer`
  - sequencing class: `weaker_inconsistency_stable_without_arrow`
  - final topology label: `smeared_transfer`
  - consistency failure total: `2`
  - tri-root failure total: `2`
  - protected arrow step under weaker inconsistency: `-1`
  - protected arrow step under tri-root inconsistency: `-1`
  - positive mismatch density: `0.0027`
  - tri-root motif correlation: `0.7903`
- `Hybrid noise plus tri-root inconsistency` / `counter_propagating_corridor`
  - protected topology label: `smeared_transfer`
  - sequencing class: `weaker_inconsistency_stable_without_arrow`
  - final topology label: `smeared_transfer`
  - consistency failure total: `2`
  - tri-root failure total: `2`
  - protected arrow step under weaker inconsistency: `-1`
  - protected arrow step under tri-root inconsistency: `-1`
  - positive mismatch density: `0.0828`
  - tri-root motif correlation: `0.2761`
- `Noise-seeded tri-root inconsistency` / `phase_ordered_symmetric_triad`
  - protected topology label: `split_channel_exchange`
  - sequencing class: `ordered_without_tri_root_inconsistency`
  - final topology label: `split_channel_exchange`
  - consistency failure total: `0`
  - tri-root failure total: `0`
  - protected arrow step under weaker inconsistency: `-1`
  - protected arrow step under tri-root inconsistency: `-1`
  - positive mismatch density: `0.0802`
  - tri-root motif correlation: `0.0000`
- `Tri-root motif-seeded inconsistency` / `phase_ordered_symmetric_triad`
  - protected topology label: `split_channel_exchange`
  - sequencing class: `ordered_without_tri_root_inconsistency`
  - final topology label: `split_channel_exchange`
  - consistency failure total: `0`
  - tri-root failure total: `0`
  - protected arrow step under weaker inconsistency: `-1`
  - protected arrow step under tri-root inconsistency: `-1`
  - positive mismatch density: `0.0019`
  - tri-root motif correlation: `0.7903`
- `Hybrid noise plus tri-root inconsistency` / `phase_ordered_symmetric_triad`
  - protected topology label: `split_channel_exchange`
  - sequencing class: `ordered_without_tri_root_inconsistency`
  - final topology label: `split_channel_exchange`
  - consistency failure total: `0`
  - tri-root failure total: `0`
  - protected arrow step under weaker inconsistency: `-1`
  - protected arrow step under tri-root inconsistency: `-1`
  - positive mismatch density: `0.0816`
  - tri-root motif correlation: `0.2761`

Interpretation boundary:
- this scan asks whether tri-root and 3-layer motif forcing can scale the weaker-kernel inconsistency signal opened in C0.8 into a broadly reproducible protected irreversibility mechanism
- topology labels remain persistence-class observables and the protected topology-stability rule stays active
- if consistency failures and protected arrow events both cross the `6/9` bar, the branch supports a broad irreversibility claim; otherwise the correct outcome is a stronger boundary statement, not a forced Phase II interpretation