PROJECT: STAGE 10C — ATLAS-2 TRANSITION STRUCTURE
REPO CONTEXT: HAOS/IIP frozen architecture
MODE: paired-axis perturbation pilot
PRIORITY: transition topology before scaling

Core philosophy
Atlas-0 mapped baseline morphology.
Atlas-1 mapped single-axis resilience.
Atlas-2 maps how regime labels deform when controlled stresses interact.

Working motto
Map transitions, not just states.

NON-NEGOTIABLE CONSTRAINTS
	1.	Do not modify the frozen operator architecture.
	2.	Use the same observable set and regime classifier as Atlas-0 and Atlas-1.
	3.	Vary exactly two perturbation axes per run.
	4.	Keep perturbation strengths fixed across the pilot.
	5.	Do not introduce physical interpretation language.
	6.	Preserve full reproducibility (explicit seeds, deterministic naming).
	7.	Keep implementation minimal and inspectable.

HIGH-LEVEL GOAL

Construct a pilot transition map showing how baseline regime labels change under paired perturbations.

This stage studies interaction structure in regime space, not refinement, convergence, or continuum interpretation.

PILOT SEED SET

Reuse the curated Atlas-1 baseline seeds:
	•	3 Ballistic coherent
	•	3 Ballistic dispersive
	•	2 Diffusive
	•	1 Chaotic or irregular
	•	1 Fragmenting

Total baseline seeds: 10.

PERTURBATION PAIRS

Test only two paired-axis stresses:

Pair A
	•	edge-weight noise = 0.03
	•	anisotropic bias = 1.05

Pair B
	•	defect density = 0.03
	•	kernel-width jitter = 0.05

Total runs:
10 seeds × 2 pairs = 20 runs.

MINIMAL OBSERVABLE SET

Use the exact Atlas-0 observable set:
	•	packet center trajectory
	•	packet width trace
	•	spectral centroid and spread
	•	sector leakage norm
	•	norm or energy drift
	•	constraint norm
	•	anisotropy score
	•	coherence score
	•	regime label
	•	sector label

NEW FIELD

Add:

transition_type

Allowed values:
	•	stable
	•	regime_shift
	•	regime_softening
	•	fragmentation_induced
	•	chaos_induced

TRANSITION CLASSIFICATION RULE

Define transition_type only from:
	•	baseline regime label
	•	perturbed regime label
	•	width growth
	•	coherence change
	•	fragmentation indicators
	•	chaos threshold crossing

No new theoretical metrics.

ATLAS-2 QUESTIONS
	•	Do regime labels compose smoothly under paired perturbations?
	•	Are some baseline regimes structurally stable under combined stress?
	•	Do certain perturbation pairs produce systematic fragmentation corridors?
	•	Does the transverse sector remain more persistent under stress interaction?
	•	Are there clean transition clusters or scattered outcomes?

RUN RECORD FORMAT

Extend the Atlas CSV schema with:
	•	perturbation_pair
	•	baseline_regime
	•	perturbed_regime
	•	transition_type

VISUAL OUTPUT

For each run:
	•	field snapshot
	•	spectrum snapshot
	•	constraint trace
	•	width trace

Additionally:
	•	transition matrix plot
	•	regime transition map

ENGINEERING TARGET

Deliverables:
	•	stage10_atlas2.py
	•	atlas2_runs.json
	•	minimal transition labeling utility
	•	stamped CSV and JSON output
	•	short Atlas-2 pilot note

STOP CONDITIONS
	•	fixed final time
	•	optional early stop on constraint blow-up
	•	optional norm collapse threshold

INTERPRETATION BOUNDARY

Atlas-2 establishes:
	•	transition structure in regime space
	•	persistence patterns under paired stress
	•	morphology-level topology

It does not establish:
	•	geometric reconstruction
	•	particle interpretation
	•	continuum limits
	•	gauge or relativistic claims

SUCCESS CRITERION

Atlas-2 pilot is successful if:
	•	regime transitions are interpretable
	•	transition clusters are visible
	•	classifier behavior remains stable
	•	results justify or postpone Stage 11 refinement work
