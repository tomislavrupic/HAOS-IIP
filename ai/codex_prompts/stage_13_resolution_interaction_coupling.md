CODEX MASTER PROMPT — STAGE 13

Resolution-Interaction Coupling Atlas

Role

You are implementing Stage 13 of the pre-geometric atlas program.

All previous stages are frozen.
You must not modify earlier operators, thresholds, observables, or classification logic.

Your task is to test how collective interaction morphology and coarse-grained envelope structure behave when resolution scale becomes part of the experiment itself.

This stage is still operator-level and diagnostic only.

No continuum claims.
No spacetime interpretation.
No force or particle language.
No nonlinear claims.

⸻

Conceptual objective

Stage 11 established:
	•	transient binding corridors
	•	metastable composite envelopes
	•	overlap geometry as the dominant separator

Stage 12 established:
	•	basic coarse-grained basin structure
	•	both compact basin-dominated and diffuse coarse regimes
	•	one multi-basin fluctuating case

Stage 13 asks:

How do interaction motifs change when the resolution scale itself becomes part of the dynamics?

More precisely:
Do regime labels and basin structures persist when
	•	packet interaction scale
	•	coarse-graining scale
	•	lattice resolution
are varied together in a controlled way?

Stage 13 therefore studies scale coupling, not a new phenomenon class.

⸻

Why this stage exists

You now know that:
	•	coherent packets can form metastable composites
	•	composites coarse-grain into basin envelopes
	•	basin topology depends strongly on geometry

You do not yet know whether:
	•	those envelopes are scale-robust
	•	regime transitions are resolution-stable
	•	interaction strength is partly a discretization artifact

Stage 13 tests exactly that.

⸻

Frozen baseline

All runs must:
	•	use the frozen kernel-induced cochain operator
	•	use the projected transverse sector
	•	use periodic boundary
	•	inherit packet family from the strongest coherent Stage 11 / Stage 12 seed:
	•	bandwidth 0.1
	•	carrier cosine
	•	central k = 1.0
	•	amplitude scale 0.5

You must reuse the Stage 11 and Stage 12 label logic unchanged.

⸻

Minimal pilot design

Use only four representative frozen cases:
	1.	tight clustered pair
	2.	asymmetric density cluster
	3.	three-packet chain
	4.	counter-propagating corridor

These four cover:
	•	strong composite / basin cases
	•	one fluctuating chain case
	•	one diffuse corridor case

Do not expand beyond this pilot.

⸻

Coupled axis sets

Test exactly three coupled scale changes.

Axis set A
	•	resolution: N → 2N
	•	coarse sigma: σ → 2σ

Axis set B
	•	resolution: N → 2N
	•	packet width: w → w/2

Axis set C
	•	resolution: N → 2N
	•	packet separation: d → d/2

Do not vary these axes independently in Stage 13.
The point of this stage is coupled scale behavior.

⸻

Reused observables

Reuse all Stage 11 and Stage 12 observables, including:
	•	pairwise centroid distance
	•	overlap integral
	•	phase-lock indicator
	•	composite lifetime
	•	exchange / merger flag
	•	coarse basin count
	•	dominant basin area fraction
	•	coarse persistence
	•	envelope variance
	•	constraint norm
	•	sector leakage

⸻

New Stage 13 observables

Add only resolution-coupling diagnostics.

1. Basin persistence across resolution

Track:
	•	interaction label at each resolution
	•	coarse label at each resolution

2. Resolution drift type

Allowed labels:
	•	stable under refinement
	•	softening under refinement
	•	class flip under refinement
	•	coarse-only drift
	•	interaction-only drift

3. Metric drift summary

For each case track:
	•	composite lifetime scaling factor
	•	overlap integral scaling
	•	centroid drift scaling
	•	dominant basin area scaling
	•	coarse persistence scaling
	•	constraint norm scaling

4. Recoverability flag

Mark whether a case becomes:
	•	more stable with refinement
	•	less stable with refinement
	•	unchanged with refinement

⸻

Stage 13 questions
	•	Do Stage 11 metastable composites remain metastable under coupled scale changes?
	•	Do Stage 12 basin-dominated cases remain basin-dominated when resolution and smoothing co-vary?
	•	Is the multi-basin fluctuating chain structural or discretization-sensitive?
	•	Do diffuse corridor cases sharpen or remain diffuse under scale coupling?
	•	Are interaction-level and coarse-level labels coupled or decoupled under refinement?

⸻

Outputs

For each case and coupled-axis set produce:
	•	field snapshot
	•	overlap / pair-distance trace
	•	envelope snapshot
	•	basing segmentation plot

Additionally produce:
	•	resolution-transition matrix
	•	interaction / coarse persistence table
	•	metric-drift summary plot

Export:
	•	stamped JSON summary
	•	stamped CSV run table
	•	short Stage 13 pilot note

⸻

Output labels

Keep the vocabulary minimal:
	•	scale-robust composite
	•	scale-fragile composite
	•	scale-induced diffusion
	•	resolution-dependent morphology

These are summary descriptors for the pilot read, not replacements for the frozen Stage 11 or Stage 12 labels.

⸻

Interpretation boundary

Stage 13 establishes only:
	•	refinement stability or instability of interaction morphology
	•	refinement stability or instability of coarse envelope classes
	•	how the two are coupled under scale variation

It does not establish:
	•	continuum limits
	•	emergent geometry
	•	nonlinear forces
	•	particle structure

⸻

Success criterion

Stage 13 pilot is successful if:
	•	at least one composite class remains stable under coupled scaling
	•	at least one class visibly drifts under coupled scaling
	•	interaction-level and coarse-level stability can be compared cleanly
	•	the pilot clarifies whether Stage 14 should become continuum probing or further atlas refinement
