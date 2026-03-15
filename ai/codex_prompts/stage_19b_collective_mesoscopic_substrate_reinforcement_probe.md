Stage 19B — Master Prompt
(Collective / mesoscopic substrate reinforcement probe)

Purpose:
Test whether geometry-like persistence can emerge only when kernel back-reaction is not purely local, but instead responds to mesoscopic collective structure (cluster-level or basin-level signals), while remaining strictly bounded and HAOS-compliant.

This stage does not strengthen ε.
It changes the class of substrate law.

⸻

Conceptual shift from 19A → 19B

Stage 19A tested:

local pair-density → tiny kernel deformation

Result:
bounded, real, but regime-inert

Stage 19B tests:

collective morphology → slow kernel reinforcement field

Key hypothesis:

A substrate may only start behaving “geometry-like” when reinforcement depends on coarse structure, not microscopic collisions.

⸻

Allowed kernel evolution class

Kernel becomes weakly time-dependent:

K(x,y,t) = K₀(x,y) · (1 + δ·S(x,y,t))

with constraints:
	•	δ ≪ 1
	•	S bounded and low-frequency
	•	no explicit metric construction
	•	no long-range instantaneous forcing
	•	update slower than packet time scale

This preserves:
	•	HAOS recoverability
	•	non-absorbing null
	•	bounded composition

⸻

Reinforcement signal candidates (test branches)

You run three distinct substrate laws, not parameter sweeps.

Branch A — cluster-density reinforcement
Kernel strengthens slightly where:
	•	persistent multi-packet clusters exist
	•	cluster membership stable over window τ

Example signal:

S_A(x,y,t) ∝ smoothed cluster co-occupancy

Interpretation:
substrate “remembers” regions of collective occupation.

⸻

Branch B — basin-residence reinforcement
Kernel responds to:
	•	long dwell inside coarse dynamical basins

Signal:

S_B(x,y,t) ∝ local time-averaged persistence score

Interpretation:
substrate aligns with stable attractor geography

⸻

Branch C — transport-corridor reinforcement
Kernel evolves along:
	•	repeatedly used propagation corridors

Signal:

S_C(x,y,t) ∝ directional flux memory

Interpretation:
substrate begins encoding preferred transport topology

⸻

Representatives (same as Stage 17–19A)

Use exactly:
	•	phase_ordered_symmetric_triad
	•	counter_propagating_corridor
	•	clustered_composite_anchor

No new morphology classes.

⸻

Strength ladder

Per branch:
	•	null (δ = 0)
	•	very weak
	•	weak

Do not exceed weak.

Goal is detection of regime change, not forcing one.

⸻

Observables (promotion gate)

A run promotes only if at least one appears:
	1.	composite lifetime increase
	2.	binding persistence increase
	3.	coarse basin persistence increase
	4.	new stabilized ordering class
	5.	transport channel locking

Additionally require:
	•	sustained effect
	•	not transient
	•	not proxy-sensitive
	•	survives morphology cross-check

⸻

Negative freeze rule

Freeze Stage 19B negative if:
	•	reinforcement field remains bounded
	•	kernel deformation non-zero but tiny
	•	no regime-level persistence gain
	•	no ordering capture
	•	no basin stabilization
	•	no refinement promotion

This is a valid physical result.

It would imply:

collective substrate memory at this scale is still insufficient to seed geometry-like structure.

⸻

Expected diagnostic outputs

For each run:
	•	mean / median persistence metrics
	•	cluster dwell statistics
	•	reinforcement field norm
	•	max local deformation
	•	morphology comparison vs null
	•	regime label

And global:
	•	response matrix
	•	deformation summary
	•	persistence comparison plots

⸻

Strategic interpretation axis

Stage 19B is the first true “proto-geometry” test in the atlas.

Possible reads:

If negative:
Architecture likely requires
→ thresholded reinforcement
→ multiscale coupling
→ or explicit constraint law

If weak positive:
You have first hint of
→ emergent effective metric memory

If strong positive (unlikely at weak δ):
You discovered a genuine regime boundary.
