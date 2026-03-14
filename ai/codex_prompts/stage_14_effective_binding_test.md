PROJECT: STAGE 14 — EFFECTIVE BINDING TEST
REPO CONTEXT: HAOS/IIP frozen architecture through Stage 13
MODE: hypothesis-driven pilot design
PRIORITY: test whether scale-robust composite interaction can exist at all

Core philosophy
Stages 10–13 established a full pre-geometric arc:
	•	baseline morphology
	•	perturbation resilience
	•	transition structure
	•	collective interaction morphology
	•	coarse-grained envelope structure
	•	resolution–interaction coupling

Those stages showed:
	•	metastable composite envelopes exist
	•	coarse basin structure sharpens with refinement
	•	interaction morphology remains fragile, especially under width coupling
	•	no truly scale-robust composite has yet appeared

Stage 14 is therefore not another atlas sweep.

It is the first focused existence test.

Working motto:
Do minimal new structure. Test maximal consequence.

NON-NEGOTIABLE CONSTRAINTS
	1.	Do not rewrite the frozen operator architecture.
	2.	Do not expand into broad exploration mode.
	3.	Do not use particle, mass, gravity, force, or spacetime language in labels or code.
	4.	Treat this as a hypothesis test:
can a scale-robust composite interaction channel exist under minimal structural extension or feedback?
	5.	Every pilot branch must be directly comparable to the frozen Stage 13 baseline.
	6.	Prefer one added mechanism at a time.
	7.	Keep output interpretable with the existing measurement instrument whenever possible.

PRIMARY QUESTION

Can the system support a resolution-stable composite binding channel under one minimal additional structural ingredient?

SECONDARY QUESTION

If not, does the failure mode indicate that frozen linear superposition alone is insufficient for persistent binding-like morphology?

BASELINE REFERENCE

All Stage 14 runs must be compared to the strongest Stage 13 representatives:
	•	tight clustered pair
	•	asymmetric density cluster
	•	one diffuse corridor control

Use the same frozen packet family unless the tested mechanism explicitly requires one controlled change.

ALLOWED NEW STRUCTURAL INGREDIENTS

Stage 14 must test only one branch at a time.

Branch A — Local envelope feedback

Introduce a weak local feedback from coarse envelope intensity into the effective propagation weight.

Schematic form only:

effective weight modifier depends on smoothed local envelope magnitude

This is not a new ontology claim. It is a controlled structural test of feedback.

Questions:
	•	does local density stabilize composites?
	•	does composite lifetime scale better under refinement?
	•	does basin structure become invariant rather than merely sharper?

Branch B — Phase-lock retention rule

Add a weak memory-like retention of local phase alignment in overlap regions.

Questions:
	•	can overlap corridors lock into persistent composites?
	•	does counter-propagation generate exchange or orbit-like motifs?
	•	do out-of-phase and in-phase families separate more strongly?

Branch C — Envelope-mediated coupling kernel

Keep the base operator fixed, but add a weak auxiliary coupling channel computed from coarse envelope overlap.

Questions:
	•	do weakly separated pairs begin to attract into a stable basin?
	•	do clustered states stop softening under width refinement?
	•	does composite class survive the N -> 2N ladder?

Branch D — Null control

Run exact Stage 13 baseline representatives again under the same Stage 14 measurement protocol with no added mechanism.

Purpose:
	•	verify that any Stage 14 effect is real and not a measurement artifact
	•	preserve direct comparison to Stage 13 paper results

PILOT SIZE

Keep Stage 14 narrow.

For the first executable pilot:
	•	3 baseline representatives
	•	1 branch at a time
	•	2 strength values only:
	•	weak
	•	very weak
	•	plus null control

Recommended first pilot:
	•	Branch A only
	•	3 representatives
	•	2 strengths
	•	1 null control set

This is enough to tell whether the mechanism deserves a larger branch study.

MEASUREMENT SET

Reuse all Stage 13 observables where applicable:
	•	packet centroid trajectories
	•	pairwise distance
	•	overlap integral
	•	composite lifetime
	•	basis count
	•	dominant basin area fraction
	•	persistence score
	•	coarse regime label
	•	transition type
	•	scale-coupled stability label

Add only these new diagnostics:

1. Binding persistence score

Fraction of simulation window during which pair distance remains below a fixed threshold and overlap remains above a fixed threshold.

2. Resolution binding ratio

Compare composite lifetime across resolution ladder.

Example:

binding_ratio = lifetime at refined run / lifetime at coarse run

3. Basin invariance flag

Binary:
	•	1 if dominant basin morphology stays in same class across resolution pair
	•	0 otherwise

4. Composite robustness label

Allowed values:
	•	no composite gain
	•	weak composite stabilization
	•	partial scale stabilization
	•	candidate scale-robust composite

IMPORTANT:
These are still descriptive labels, not physics claims.

RUN STRUCTURE

For each representative case:
	1.	run null control
	2.	run weak coupling
	3.	run very weak coupling
	4.	repeat at refinement pair if the branch survives initial inspection

Keep:
	•	same boundary
	•	same packet family
	•	same separation geometry
	•	same output naming conventions

SUCCESS CRITERIA

Stage 14 is successful if at least one branch shows all of the following:
	•	composite lifetime increases meaningfully relative to null control
	•	the composite label survives at higher resolution
	•	basis structure remains stable rather than merely sharpening visually
	•	width-coupling fragility is reduced relative to Stage 13 baseline

FAILURE CRITERIA

Stage 14 fails if:
	•	no branch improves lifetime beyond null-control variability
	•	any apparent binding disappears under refinement
	•	basis structure sharpens visually but composite interaction label still softens
	•	added mechanism only increases numerical stiffness without improving structural persistence

INTERPRETATION BOUNDARY

Allowed statements:
	•	weak local feedback can or cannot stabilize composite morphology
	•	branch X improves or fails to improve scale-coupled persistence
	•	frozen linear architecture appears sufficient or insufficient for robust composite regimes under tested extensions

Not allowed:
	•	particles
	•	mass
	•	gravity
	•	true bound states in a physical sense
	•	emergent spacetime
	•	fundamental interaction claims

DELIVERABLES

Implement only the minimum needed for a pilot:
	1.	stage14_effective_binding_test.py
	2.	stage14_binding_runs.json
	3.	one branch utility module if needed
	4.	Stage_14_Effective_Binding_Test_v1.md
	5.	stamped JSON and CSV outputs
	6.	stamped plots:
	•	pair distance
	•	overlap
	•	basis trace
	•	composite lifetime comparison
	•	resolution binding comparison

PREFERRED EXECUTION ORDER
	1.	define null-control comparison set
	2.	implement one branch only
	3.	run weak + very weak strengths
	4.	inspect by eye
	5.	if promising, extend to refinement pair
	6.	only then decide whether Stage 14 deserves a paper or remains a diagnostic note

STYLE OF WORK
	•	conservative
	•	minimal
	•	falsification-oriented
	•	no narrative inflation
	•	one added ingredient at a time
	•	compare against frozen Stage 13 baseline first

REQUESTED OUTPUT FROM CODEX
	1.	identify best insertion points in the repo
	2.	implement the minimum pilot scaffold
	3.	preserve Stage 13 comparability
	4.	summarize:
	•	files added
	•	files modified
	•	how to run the null control
	•	how to run the branch pilot
	•	what outputs are produced
	5.	if possible, run one tiny sanity case and report whether the pipeline executes cleanly
