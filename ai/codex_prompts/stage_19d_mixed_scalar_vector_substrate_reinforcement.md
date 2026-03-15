Stage 19D — Mixed Scalar + Vector Substrate Reinforcement
Master Prompt (Persistence-First Pilot)

Role

You are operating inside the frozen HAOS-IIP numerical architecture.
Do not modify operators, discretization, normalization conventions, or promotion rules.

Your task is to run a minimal persistence-first probe testing whether simultaneous scalar (density-like) and vector (flux-like) substrate feedback can produce any measurable gain in:
	•	composite lifetime
	•	binding persistence
	•	coarse basin residence
	•	corridor dwell / locking

This is a law-search precursor, not a tuning exercise.

No parameter escalation.
No architectural changes.
Record bounded nulls cleanly.

⸻

Hypothesis

Previous branches established:
	•	Frozen substrate → no binding
	•	Scalar back-reaction → no persistence gain
	•	Vector / directional back-reaction → bounded deformation, no persistence gain

Stage 19D tests whether combined feedback symmetry introduces a qualitatively new mesoscopic effect:

A shallow scalar well plus weak directional bias may jointly stabilize trajectories that neither mechanism can stabilize alone.

⸻

Substrate Update Structure

At each reinforcement update step:
	1.	Scalar channel (density memory)
	•	Stress functional: smoothed amplitude density
	•	Produces symmetric kernel modulation
	•	Acts as shallow mass-like basin shaping
	2.	Vector channel (flux memory)
	•	Stress functional: smoothed transverse flux magnitude or signed current proxy
	•	Produces weak directional kernel asymmetry
	•	Acts as synthetic transport bias
	3.	Mixed reinforcement rule

K_{ij}(t+1)
=
K_{ij}^{(0)}
\Big[
1
+
\epsilon_s\, S_{ij}(t)
\Big]
\,
\exp\!\big(
\epsilon_v\, A_{ij}(t)
\big)

where:
	•	S_{ij} = scalar reinforcement field
	•	A_{ij} = antisymmetric directional reinforcement field
	•	both fields are weak, smoothed, bounded memory signals

Clamp total kernel deviation per update.

No runaway growth allowed.

⸻

Numerical Discipline

Use:
	•	single-horizon exponential memory
	•	Gaussian spatial smoothing
	•	strict kernel update clamp
	•	reinforcement decay in absence of signal

Reinforcement must remain mesoscopic and reversible.

No permanent substrate plasticity.

⸻

Representatives

Run only the canonical three:
	1.	clustered_composite_anchor
	2.	phase_ordered_symmetric_triad
	3.	counter_propagating_corridor

No new geometries.

⸻

Run Structure

Minimal pilot:
	•	scalar strength grid: low → moderate
	•	vector strength grid: low → moderate
	•	include one balanced mixed case
	•	fixed temporal horizon
	•	fixed refinement rule: promote only on earned persistence signal

Total runs: small controlled matrix (e.g. 9).

⸻

Observables

Primary persistence metrics:
	•	composite lifetime
	•	binding persistence
	•	coarse basin persistence
	•	corridor dwell / locking

Secondary diagnostics:
	•	reinforcement field norms
	•	kernel deviation norms
	•	maximum local deformation
	•	directional deflection proxy
	•	energy dispersion shift

Do not over-interpret secondary signals.

Persistence metrics dominate decision.

⸻

Success Criteria

Promotion only if any representative shows:
	•	statistically clear positive delta in persistence metric
	•	stable behavior across full horizon
	•	no numerical instability signature

Otherwise classify as:

“Mixed scalar-vector mesoscopic reinforcement insufficient at tested coupling.”

⸻

Output Requirements

For each run:
	•	structured JSON metrics record
	•	persistence comparison vs frozen baseline
	•	reinforcement amplitude summary
	•	kernel deviation bounds
	•	short qualitative regime label

After pilot completion:

Produce a Stage 19D persistence snapshot note containing:
	•	hypothesis tested
	•	parameter envelope
	•	persistence outcome
	•	structural interpretation
	•	recommendation:
	•	escalate symmetry search
	•	or terminate reinforcement class and move to Stage 20 law-search

⸻

Stop Rule

If all runs return bounded null:
	•	freeze branch cleanly
	•	do not widen coupling
	•	do not introduce new memory channels

Return control for architectural decision.
