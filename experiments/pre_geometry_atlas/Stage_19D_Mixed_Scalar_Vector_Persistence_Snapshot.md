# Stage 19D Persistence Snapshot

Hypothesis tested:

> A shallow scalar substrate well plus a weak directional transport bias may jointly stabilize trajectories that neither scalar nor vector reinforcement can stabilize alone.

Parameter envelope:
- representatives: `clustered_composite_anchor`, `phase_ordered_symmetric_triad`, `counter_propagating_corridor`
- conditions:
  - `null`
  - `balanced_low` with `eps_s = 0.01`, `eps_v = 0.01`
  - `balanced_moderate` with `eps_s = 0.02`, `eps_v = 0.02`
- memory class: `single_horizon_exponential_decay`
- temporal horizon: `tau = 0.3`
- spatial smoothing:
  - scalar `sigma = 2.0`
  - vector `sigma = 2.0`
- kernel clamp: `0.05`

Persistence outcome:
- `9/9` runs returned `no mixed substrate effect`
- no composite lifetime gain
- no binding persistence gain
- no coarse basin persistence gain
- no corridor dwell or locking gain
- no promoted `12 -> 24` follow-up

Why this still counts as a real pilot:
- scalar and vector channels were both active in every non-null run
- mixed kernel deformation was bounded but clearly nonzero
- the strongest mixed response remained mesoscopic rather than numerically trivial
- the corridor case showed deflection texture, but no persistence consequence

Representative mixed deformation scale:
- clustered composite anchor:
  - mixed field norm up to `6.50e-3`
  - kernel deviation norm up to `6.50e-3`
  - max local deformation up to `3.47e-2`
- phase-ordered symmetric triad:
  - mixed field norm up to `9.42e-3`
  - kernel deviation norm up to `9.42e-3`
  - max local deformation up to `2.83e-2`
- counter-propagating corridor:
  - mixed field norm up to `6.70e-3`
  - kernel deviation norm up to `6.70e-3`
  - max local deformation up to `3.38e-2`

Structural interpretation:

> Mixed scalar-vector mesoscopic reinforcement is insufficient at the tested couplings.

The combined law enters the substrate more strongly than the pure directional branch and at a comparable mesoscopic scale to the stronger scalar branch, but still fails to create persistence capture, basin stabilization, or corridor locking.

Recommendation:
- do not widen the same reinforcement class by coefficient escalation
- either escalate symmetry search cleanly
- or terminate reinforcement-class exploration and move to `Stage 20` effective-law search
