CODEX MASTER PROMPT — STAGE 12

Emergent Coarse-Grained Field Structure Atlas

Role

You are implementing Stage 12 of the pre-geometric atlas program.

All previous stages are frozen.
You must not modify earlier operators, thresholds, observables, or classification logic.

Your task is to build a coarse-grained morphology layer that studies whether multi-packet superposition fields exhibit:
	•	persistent envelope structures
	•	basin-like dominance regions
	•	scale-dependent regime stability

This stage is still operator-level and diagnostic only.

No spacetime interpretation.
No particle or force language.
No nonlinear claims.

⸻

Conceptual objective

Stage 11 established:
	•	transient binding corridors
	•	metastable composite envelopes
	•	geometry-driven interaction strength

Stage 12 asks:

When we observe the collective field at larger spatial and temporal scales,
do stable coarse structures emerge that organize the local morphology?

We therefore introduce:
	•	spatial coarse-graining
	•	envelope-dominance diagnostics
	•	basin persistence measures

⸻

Frozen baseline

All runs must:
	•	use the frozen kernel-induced cochain operator
	•	use the projected transverse sector
	•	use periodic boundary
	•	inherit packet family from strongest Stage 11 seed:
	•	bandwidth 0.1
	•	carrier cosine
	•	central k = 1.0
	•	amplitude scale 0.5

⸻

Run grammar

Create a small pilot atlas with 5 canonical coarse-field runs:
	1.	tight clustered pair
	2.	three-packet chain
	3.	symmetric wide triad
	4.	asymmetric density cluster
	5.	counter-propagating composite corridor

Axes varied per run:
	•	packet count
	•	geometric packing density
	•	relative phase pattern
	•	amplitude symmetry

No other parameter variation.

⸻

New Stage 12 observables

Add coarse-field diagnostics computed from the evolving field:

1. Envelope field

Gaussian-smoothed magnitude field:

E_sigma(x,t) = G_sigma * |ψ(x,t)|

Test at fixed smoothing scales:
	•	σ = 2 grid units
	•	σ = 4 grid units

2. Basin-dominance map

Define connected regions where:

E_sigma > α · max(E_sigma)

with α = 0.6.

Track:
	•	basin count over time
	•	dominant basin area fraction
	•	basin lifetime

3. Coarse persistence score

Measure stability of the dominant basin:

persistence = time_fraction( dominant basin identity stable )

4. Coarse-scale regime label

Classify each run into:
	•	basin-dominated regime
	•	multi-basin fluctuating regime
	•	diffuse coarse field regime

Label rule uses only:
	•	basin lifetime
	•	dominant area fraction
	•	envelope variance

⸻

Outputs

For each run produce:
	•	field snapshot
	•	envelope snapshot
	•	basin segmentation plot
	•	dominant-basin area trace
	•	persistence trace

Export:
	•	stamped JSON summary
	•	stamped CSV run table
	•	summary coarse-regime map

⸻

Constraints

You must preserve:
	•	linear dynamics interpretation
	•	superposition-based reading
	•	no emergent force claims

You must explicitly state:

Stage 12 identifies coarse-grained organizational morphology
and does not imply nonlinear effective interactions.

⸻

Pilot size

Exactly 5 runs.

Do not scale atlas yet.

Goal is:
	•	detect whether basin-dominated envelopes are reproducible
	•	test if coarse regimes align with Stage 11 composite strength

⸻

Completion criteria

Stage 12 pilot is considered successful if:
	•	at least one run shows high persistence basin dominance
	•	at least one run remains diffuse
	•	classifier separation is visually interpretable

If successful:
	•	freeze Stage 12 pilot
	•	then design Stage 13 resolution-interaction coupling.
