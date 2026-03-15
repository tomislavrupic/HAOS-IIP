# Stage 19B / Branch 1 - Basin-Residence Reinforcement Negative Snapshot

Timestamped summary:
- `data/20260315_105216_stage19b_mesoscopic.json`
- `data/20260315_105216_stage19b_mesoscopic.csv`

Summary result:
- `9/9` runs labeled `no mesoscopic substrate effect`
- promoted follow-up: `none`

What was tested:
- Stage 19B Branch 1 replaced local substrate reinforcement with bounded mesoscopic basin-residence reinforcement.
- The substrate law remained weak, bounded, symmetric, and support-preserving.
- The representative set stayed fixed:
  - `clustered_composite_anchor`
  - `phase_ordered_symmetric_triad`
  - `counter_propagating_corridor`

What stayed unchanged:
- composite lifetime
- binding persistence
- coarse basin persistence
- ordering capture / stabilized class

Why this is a real negative:
- the branch was active, not dead
- reinforcement field norms were clearly nonzero
- kernel deformation was bounded but measurable
- no regime-level observable moved

Representative deformation scale:
- `clustered_composite_anchor`
  - reinforcement field norm: `0.3902` / `0.3888`
  - kernel deviation norm: `0.0039` / `0.0078`
- `phase_ordered_symmetric_triad`
  - reinforcement field norm: `0.2179` / `0.2179`
  - kernel deviation norm: `0.0022` / `0.0044`
- `counter_propagating_corridor`
  - reinforcement field norm: `0.2451` / `0.2451`
  - kernel deviation norm: `0.0025` / `0.0049`

Core conclusion:

Basin-residence reinforcement produces bounded, nonzero mesoscopic substrate deformation, but does not improve composite lifetime, binding persistence, coarse basin persistence, or ordering capture at the tested strengths.

Interpretation boundary:
- this is a negative result for occupation-history reinforcement at mesoscopic scale
- it does not rule out other substrate laws
- the next coherent branch is transport-corridor reinforcement, not stronger `delta`
