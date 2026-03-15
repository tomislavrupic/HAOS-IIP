# Stage 19B Mesoscopic Substrate Reinforcement v1

Timestamped summary: `data/20260315_105216_stage19b_mesoscopic.json`
Timestamped run table: `data/20260315_105216_stage19b_mesoscopic.csv`

Output-label counts: {'no mesoscopic substrate effect': 9}

This Stage 19B pilot tests whether bounded basin-residence reinforcement on the kernel edge weights can create the first mesoscopic substrate-assisted persistence effect.

Base comparisons:

- `clustered_composite_anchor` / `null_control`
  - composite lifetime delta: `0.0000`
  - binding persistence delta: `0.0000`
  - coarse persistence delta: `0.0000`
  - reinforcement field norm: `0.0000`
  - kernel deviation norm: `0.0000`
  - basin dwell fraction: `0.0000`
  - gate met: `0`
  - label: `no mesoscopic substrate effect`
- `clustered_composite_anchor` / `very_weak_basin_residence`
  - composite lifetime delta: `0.0000`
  - binding persistence delta: `0.0000`
  - coarse persistence delta: `0.0000`
  - reinforcement field norm: `0.3902`
  - kernel deviation norm: `0.0039`
  - basin dwell fraction: `0.9866`
  - gate met: `0`
  - label: `no mesoscopic substrate effect`
- `clustered_composite_anchor` / `weak_basin_residence`
  - composite lifetime delta: `0.0000`
  - binding persistence delta: `0.0000`
  - coarse persistence delta: `0.0000`
  - reinforcement field norm: `0.3888`
  - kernel deviation norm: `0.0078`
  - basin dwell fraction: `0.9867`
  - gate met: `0`
  - label: `no mesoscopic substrate effect`
- `phase_ordered_symmetric_triad` / `null_control`
  - composite lifetime delta: `0.0000`
  - binding persistence delta: `0.0000`
  - coarse persistence delta: `0.0000`
  - reinforcement field norm: `0.0000`
  - kernel deviation norm: `0.0000`
  - basin dwell fraction: `0.0000`
  - gate met: `0`
  - label: `no mesoscopic substrate effect`
- `phase_ordered_symmetric_triad` / `very_weak_basin_residence`
  - composite lifetime delta: `0.0000`
  - binding persistence delta: `0.0000`
  - coarse persistence delta: `0.0000`
  - reinforcement field norm: `0.2179`
  - kernel deviation norm: `0.0022`
  - basin dwell fraction: `1.0000`
  - gate met: `0`
  - label: `no mesoscopic substrate effect`
- `phase_ordered_symmetric_triad` / `weak_basin_residence`
  - composite lifetime delta: `0.0000`
  - binding persistence delta: `0.0000`
  - coarse persistence delta: `0.0000`
  - reinforcement field norm: `0.2179`
  - kernel deviation norm: `0.0044`
  - basin dwell fraction: `1.0000`
  - gate met: `0`
  - label: `no mesoscopic substrate effect`
- `counter_propagating_corridor` / `null_control`
  - composite lifetime delta: `0.0000`
  - binding persistence delta: `0.0000`
  - coarse persistence delta: `0.0000`
  - reinforcement field norm: `0.0000`
  - kernel deviation norm: `0.0000`
  - basin dwell fraction: `0.0000`
  - gate met: `0`
  - label: `no mesoscopic substrate effect`
- `counter_propagating_corridor` / `very_weak_basin_residence`
  - composite lifetime delta: `0.0000`
  - binding persistence delta: `0.0000`
  - coarse persistence delta: `0.0000`
  - reinforcement field norm: `0.2451`
  - kernel deviation norm: `0.0025`
  - basin dwell fraction: `1.0000`
  - gate met: `0`
  - label: `no mesoscopic substrate effect`
- `counter_propagating_corridor` / `weak_basin_residence`
  - composite lifetime delta: `0.0000`
  - binding persistence delta: `0.0000`
  - coarse persistence delta: `0.0000`
  - reinforcement field norm: `0.2451`
  - kernel deviation norm: `0.0049`
  - basin dwell fraction: `1.0000`
  - gate met: `0`
  - label: `no mesoscopic substrate effect`

Promoted follow-up: `none`

Interpretation boundary:
- any positive effect is read only as bounded mesoscopic substrate-assisted persistence
- no geometric or force-law claim is made here
