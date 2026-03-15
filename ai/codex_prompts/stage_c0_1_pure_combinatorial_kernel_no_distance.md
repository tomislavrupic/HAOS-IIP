# Stage C0.1 - Combinatorial Kernel Austerity Probe

## Purpose

Run a sharp control branch that answers one question cleanly:

How much of the observed structure survives when coordinate-distance weighting is removed from the interaction kernel?

This stage does not replace the current program. It is an austerity probe against hidden metric dependence.

The target is not novelty. The target is clarification.

If a claimed effect survives without geometric distance in the kernel, that is strong evidence of graph-structural content.

If it collapses, that is also useful:
- the Gaussian kernel was doing material ontological work
- the current model should be described more carefully as weighted-graph mesoscopic pre-geometry rather than strict combinatorial austerity

## Core Rule

In the new branch, interaction weights must be built with no coordinate distance at all.

Allowed ingredients:
- adjacency
- incidence
- node degree
- shortest-path distance in edge count only
- local motif counts
- triangle participation
- square or 2-cell participation
- discrete shell number from a seed or source node
- orientation or incidence sign already present in the complex

Forbidden ingredients:
- Euclidean coordinates inside the kernel builder
- `||x_i - x_j||`
- Gaussian radius in physical or embedded space
- anisotropic kernel stretch tied to coordinates
- any continuous embedding distance
- any hidden lookup table derived from geometric separation

Important boundary:
- existing coordinates may still be used only where already required to instantiate the frozen substrate, seed the same representative initial conditions, or produce the usual plots
- they must not be consumed by the new combinatorial kernel law

## Frozen Comparison Policy

Stage C0.1 is a control, so everything except kernel class stays as fixed as possible.

Do not change:
- representative initial conditions
- packet family
- operator sector under test
- boundary condition
- measurement definitions
- classifier thresholds
- resolution, at first pass
- memory or reinforcement laws

Do not add:
- delayed memory
- substrate back-reaction
- registration gates
- new packet species
- new DK couplings
- coefficient escalation

The only intentional intervention is:

Gaussian metric-coded kernel -> combinatorial kernel family

## Kernel Classes for C0.1

Test exactly three kernels in the first pass.

### 1. Gaussian Baseline

Use the current default kernel unchanged as the control condition.

Purpose:
- preserve direct comparability with prior stages
- measure what is lost, preserved, or deformed under austerity replacements

### 2. Pure Adjacency Kernel

Maximally austere no-distance baseline.

Canonical form:

`K_ij = 1` if `i` and `j` are adjacent, else `0`

Optional normalized form:

`K_ij = 1 / sqrt(d_i d_j)` if `i` and `j` are adjacent, else `0`

Constraints:
- no beyond-nearest-neighbor support
- no geometric bandwidth parameter
- if normalization is used, choose it once globally and keep it fixed for all runs

This is the primary no-smuggling baseline.

### 3. Graph-Shell Kernel

Use shortest-path distance in edge count only.

Canonical form:

`K_ij = f(l_ij)`

where `l_ij` is graph distance and `f` has finite support with no Gaussian profile.

Recommended first form:
- shell 1 weight `a`
- shell 2 weight `b`
- shell 3 weight `c`
- zero beyond shell 3

Example:
- `a = 1.0`
- `b = 0.5`
- `c = 0.25`

Constraints:
- choose one shell profile once before the scan
- keep it fixed across all representatives
- no tuning by run outcome
- no conversion from shell number to continuous radius

This tests whether discrete nonlocality alone preserves any of the reported structure.

## Not Yet Included

Do not include motif-weighted kernels in C0.1.

They are valid for a later branch because they still remain combinatorial, but they are not part of the first austerity control. The first pass should stay narrow:
- Gaussian baseline
- pure adjacency
- graph-shell

If signal survives C0.1, motif-weighted kernels become the natural C0.2 follow-up.

## Representative Initial Conditions

Use the same small representative set already established in the persistence-first program:
1. `clustered_composite_anchor`
2. `counter_propagating_corridor`
3. `phase_ordered_symmetric_triad`

These are enough for the first pass because they already probe distinct structural questions:
- anchor or composite retention
- transport corridor behavior
- symmetric multi-packet ordering

Do not launch a large atlas.

## Minimal Experiment Design

Stage C0.1 should be a minimal 9-run control matrix:

- 3 kernel classes
- 3 representative initial conditions
- 1 fixed base resolution

Recommended first pass:
- keep the same base resolution used in recent Stage 19-22 pilots
- no refinement ladder unless a non-Gaussian kernel shows a repeatable positive signal

Matrix:
- Gaussian baseline x clustered composite anchor
- Gaussian baseline x counter-propagating corridor
- Gaussian baseline x symmetric triad
- pure adjacency x clustered composite anchor
- pure adjacency x counter-propagating corridor
- pure adjacency x symmetric triad
- graph-shell x clustered composite anchor
- graph-shell x counter-propagating corridor
- graph-shell x symmetric triad

## Observables

Use the same observables already trusted in the current branch.

Primary persistence observables:
- composite lifetime
- binding persistence
- coarse basin persistence

Transport and ordering observables:
- corridor dwell time
- transport span
- exchange or collision class
- envelope ordering stability

Topology observables, only if the tested sector already supports them:
- topology label
- braid-like exchange label
- channel count
- loop count
- recirculation score

Implementation rule:
- do not invent new success metrics just to help the austerity branch look positive
- compare on the same measurement layer already used for the Gaussian branch

## Comparison Logic

The output should answer the critique directly, not narratively.

For each representative, report:
- Gaussian outcome
- pure adjacency outcome
- graph-shell outcome
- whether the qualitative class changed
- whether persistence metrics materially degraded, survived, or improved

The key table should be:

`representative x kernel class -> persistence and topology outcome`

## Interpretation Rules

### Case A - Structure Survives Under Adjacency and Graph-Shell

This is the strongest outcome.

Interpretation:
- the reported structure is not merely an artifact of embedded Gaussian distance
- at least part of the dynamics is genuinely graph-structural
- topology-like behavior survives a strict no-distance kernel test

### Case B - Structure Fails Under Pure Adjacency but Survives Under Graph-Shell

This is still informative.

Interpretation:
- nearest-neighbor austerity is too severe
- discrete nonlocality in graph shells carries real structural content
- the effect depends on combinatorial path structure, not necessarily on Euclidean embedding

### Case C - Structure Collapses Outside the Gaussian Baseline

This is a legitimate negative result.

Interpretation:
- the Gaussian kernel is doing material structural work
- the current program should not be described as metric-free in the strong sense
- the correct description is closer to weighted-graph mesoscopic pre-geometry

This is not failure. It is a clearer map of what is actually carrying the phenomenon.

## Promotion Rule

Only promote beyond the base 9-run matrix if at least one non-Gaussian kernel shows one of the following relative to Gaussian null expectations:
- preserved composite lifetime class
- preserved binding persistence class
- preserved coarse basin residence structure
- preserved exchange or braid-like classifier
- preserved corridor transport ordering

If no non-Gaussian kernel preserves any regime-level structure:
- freeze C0.1 negative
- do not broaden the search immediately
- write the conclusion plainly

## Scientific Boundary

Do not oversell either outcome.

Do not claim:
- full removal of geometry from the entire architecture, if only the kernel law changed
- emergence of topology from nothing
- proof that geometry is irrelevant

Do claim, if supported:
- whether the observed structure is robust or fragile under removal of coordinate-distance weighting in the kernel

That is the exact point of Stage C0.1.

## Deliverables

Produce:
- a 9-run runsheet or equivalent experiment specification
- one result table comparing all three kernels across the three representatives
- the usual per-run JSON outputs
- a short conclusion in one of two forms:

Either:

`Key phenomena survive removal of coordinate-distance weighting.`

Or:

`Key phenomena depend materially on metric-coded kernel structure.`

## Natural Next Step

If C0.1 shows survival under non-Gaussian kernels, the next branch should be:

Stage C0.2 - Motif-Weighted Combinatorial Kernel

If C0.1 freezes negative, the next step is not immediate escalation. The next step is honest reframing of the current architecture.
