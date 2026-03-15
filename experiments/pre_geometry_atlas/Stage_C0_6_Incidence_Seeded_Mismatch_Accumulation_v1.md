# Stage C0.6 Incidence-Seeded Mismatch Accumulation v1

Timestamped summary: `data/20260315_204114_stage_c0_6_incidence_seeded_mismatch_accumulation.json`
Timestamped run table: `data/20260315_204114_stage_c0_6_incidence_seeded_mismatch_accumulation.csv`

Architecture notice: This branch keeps graph size, adjacency and graph-shell kernels, the clustered / corridor / triad representative families, the projected-transverse sector, and the shell-weight measurement vocabulary fixed. Only the protected sequencing step is externally forced by incidence-noise or motif injections.
Baseline notice: C0.6 keeps the C0.5 protected topology-stability rule active, but deliberately injects the C0.3 and C0.4 combinatorial perturbation menus to test whether positive mismatch density and protected arrow forcing can be made reproducible.

Stage summary: positive mismatch in `9/9` runs; protected arrow forcing in `6/9` runs.
Protected arrow step distribution: `{1: 2, 2: 4, 3: 6, 4: 5, 5: 5, 6: 5, 7: 5, 8: 5}`

Per-run summary:

- `Noise-seeded protected sequencing` / `clustered_composite_anchor`
  - protected topology label: `trapped_local`
  - sequencing class: `protected_ordered_persistent`
  - final topology label: `trapped_local`
  - positive mismatch density: `0.0802`
  - positive mismatch steps: `8`
  - protected arrow step distribution: `[3, 4, 5, 6, 7, 8]`
  - kernel consistency max: `0`
- `Motif-seeded protected sequencing` / `clustered_composite_anchor`
  - protected topology label: `trapped_local`
  - sequencing class: `protected_ordered_persistent`
  - final topology label: `trapped_local`
  - positive mismatch density: `0.0064`
  - positive mismatch steps: `8`
  - protected arrow step distribution: `[1, 2, 3, 4, 5, 6, 7, 8]`
  - kernel consistency max: `0`
- `Hybrid seeded protected sequencing` / `clustered_composite_anchor`
  - protected topology label: `trapped_local`
  - sequencing class: `protected_ordered_persistent`
  - final topology label: `trapped_local`
  - positive mismatch density: `0.0802`
  - positive mismatch steps: `8`
  - protected arrow step distribution: `[3]`
  - kernel consistency max: `0`
- `Noise-seeded protected sequencing` / `counter_propagating_corridor`
  - protected topology label: `smeared_transfer`
  - sequencing class: `stable_without_ordering`
  - final topology label: `smeared_transfer`
  - positive mismatch density: `0.0816`
  - positive mismatch steps: `8`
  - protected arrow step distribution: `[]`
  - kernel consistency max: `0`
- `Motif-seeded protected sequencing` / `counter_propagating_corridor`
  - protected topology label: `smeared_transfer`
  - sequencing class: `stable_without_ordering`
  - final topology label: `smeared_transfer`
  - positive mismatch density: `0.0158`
  - positive mismatch steps: `6`
  - protected arrow step distribution: `[]`
  - kernel consistency max: `0`
- `Hybrid seeded protected sequencing` / `counter_propagating_corridor`
  - protected topology label: `smeared_transfer`
  - sequencing class: `stable_without_ordering`
  - final topology label: `smeared_transfer`
  - positive mismatch density: `0.0816`
  - positive mismatch steps: `8`
  - protected arrow step distribution: `[]`
  - kernel consistency max: `0`
- `Noise-seeded protected sequencing` / `phase_ordered_symmetric_triad`
  - protected topology label: `split_channel_exchange`
  - sequencing class: `protected_ordered_persistent`
  - final topology label: `split_channel_exchange`
  - positive mismatch density: `0.0802`
  - positive mismatch steps: `8`
  - protected arrow step distribution: `[2, 3, 4, 5, 6, 7, 8]`
  - kernel consistency max: `0`
- `Motif-seeded protected sequencing` / `phase_ordered_symmetric_triad`
  - protected topology label: `split_channel_exchange`
  - sequencing class: `protected_ordered_persistent`
  - final topology label: `split_channel_exchange`
  - positive mismatch density: `0.0041`
  - positive mismatch steps: `8`
  - protected arrow step distribution: `[1, 2, 3, 4, 5, 6, 7, 8]`
  - kernel consistency max: `0`
- `Hybrid seeded protected sequencing` / `phase_ordered_symmetric_triad`
  - protected topology label: `split_channel_exchange`
  - sequencing class: `protected_ordered_persistent`
  - final topology label: `split_channel_exchange`
  - positive mismatch density: `0.0802`
  - positive mismatch steps: `8`
  - protected arrow step distribution: `[2, 3, 4, 5, 6, 7, 8]`
  - kernel consistency max: `0`

Interpretation boundary:
- this scan asks whether external combinatorial forcing can produce positive mismatch density while the protected topology-stability rule remains active
- topology labels are still treated as persistence-class observables on the no-distance branch
- if positive mismatch remains absent or isolated, the branch is still a rigid persistence-class control rather than an incidence-driven sequencing mechanism