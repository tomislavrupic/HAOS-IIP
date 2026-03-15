# Stage C0.5 Emergent Sequencing v1

Timestamped summary: `data/20260315_202621_stage_c0_5_emergent_sequencing.json`
Timestamped run table: `data/20260315_202621_stage_c0_5_emergent_sequencing.csv`

Architecture notice: This branch keeps graph size, adjacency and graph-shell kernels, the clustered / corridor / triad representative families, the projected-transverse sector, and the shell-weight measurement vocabulary fixed. Only the sequencing rule changes.
Baseline notice: C0.5 replaces the assumed packet-evolution clock with a combinatorial support-update rule. Sequencing is derived from incidence-mismatch accumulation, kernel-update consistency, and protected topology-stability blocks.

Per-run summary:

- `Adjacency sequencing` / `clustered_composite_anchor`
  - protected topology label: `trapped_local`
  - sequencing class: `static_regime`
  - sequencing depth: `0`
  - mismatch quiet step: `0`
  - protected arrow step: `-1`
  - persistence step: `1`
  - final topology label: `trapped_local`
  - max mismatch density: `0.0000`
- `Graph-shell sequencing` / `clustered_composite_anchor`
  - protected topology label: `trapped_local`
  - sequencing class: `static_regime`
  - sequencing depth: `0`
  - mismatch quiet step: `0`
  - protected arrow step: `-1`
  - persistence step: `1`
  - final topology label: `trapped_local`
  - max mismatch density: `0.0000`
- `Protected alternating sequencing` / `clustered_composite_anchor`
  - protected topology label: `trapped_local`
  - sequencing class: `static_regime`
  - sequencing depth: `0`
  - mismatch quiet step: `0`
  - protected arrow step: `-1`
  - persistence step: `1`
  - final topology label: `trapped_local`
  - max mismatch density: `0.0000`
- `Adjacency sequencing` / `counter_propagating_corridor`
  - protected topology label: `smeared_transfer`
  - sequencing class: `static_regime`
  - sequencing depth: `0`
  - mismatch quiet step: `0`
  - protected arrow step: `-1`
  - persistence step: `1`
  - final topology label: `smeared_transfer`
  - max mismatch density: `0.0000`
- `Graph-shell sequencing` / `counter_propagating_corridor`
  - protected topology label: `smeared_transfer`
  - sequencing class: `static_regime`
  - sequencing depth: `0`
  - mismatch quiet step: `0`
  - protected arrow step: `-1`
  - persistence step: `1`
  - final topology label: `smeared_transfer`
  - max mismatch density: `0.0000`
- `Protected alternating sequencing` / `counter_propagating_corridor`
  - protected topology label: `smeared_transfer`
  - sequencing class: `static_regime`
  - sequencing depth: `0`
  - mismatch quiet step: `0`
  - protected arrow step: `-1`
  - persistence step: `1`
  - final topology label: `smeared_transfer`
  - max mismatch density: `0.0000`
- `Adjacency sequencing` / `phase_ordered_symmetric_triad`
  - protected topology label: `split_channel_exchange`
  - sequencing class: `static_regime`
  - sequencing depth: `0`
  - mismatch quiet step: `0`
  - protected arrow step: `-1`
  - persistence step: `1`
  - final topology label: `split_channel_exchange`
  - max mismatch density: `0.0000`
- `Graph-shell sequencing` / `phase_ordered_symmetric_triad`
  - protected topology label: `split_channel_exchange`
  - sequencing class: `static_regime`
  - sequencing depth: `0`
  - mismatch quiet step: `0`
  - protected arrow step: `-1`
  - persistence step: `1`
  - final topology label: `split_channel_exchange`
  - max mismatch density: `0.0000`
- `Protected alternating sequencing` / `phase_ordered_symmetric_triad`
  - protected topology label: `split_channel_exchange`
  - sequencing class: `protected_ordered_persistent`
  - sequencing depth: `0`
  - mismatch quiet step: `0`
  - protected arrow step: `2`
  - persistence step: `1`
  - final topology label: `split_channel_exchange`
  - max mismatch density: `0.0000`

Interpretation boundary:
- this scan derives sequencing from incidence-mismatch accumulation on combinatorial support updates
- topology labels are treated as persistence-class observables on the no-distance branch
- it does not claim physical time, geometry emergence, particles, or universal irreversibility