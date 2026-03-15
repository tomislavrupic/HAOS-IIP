# Stage 23.1 - Dirac-Kaehler Collision Geometry Scan

Phase III focused opener - graded collision baseline

## Purpose

Run a minimal collision scan in the Dirac-Kaehler sector to test whether graded first-order propagation produces collision morphologies not seen in the scalar / transverse baseline.

This stage is not a full atlas.
It is a focused geometry-of-encounter probe.

The question is narrow:

When two or three Dirac-Kaehler packets meet, do they:
- pass through as before,
- reflect,
- deflect,
- exchange grade,
- form longer-lived composite motion?

No memory laws.
No adaptive kernel.
No reinforcement.
No gating.
No extra fields.

Only the operator sector changes.

## Context

Phases I-II established a strong negative boundary for weak bosonic-style persistence mechanisms.

Stage 23 opens Phase III by changing operator class.

The first lawful move is to test collision geometry, not broad emergence claims.

## Frozen Architecture Assumptions

Keep fixed:
- graph substrate family already validated
- kernel base structure
- normalization conventions
- persistence metrics
- base packet family where possible
- resolution baseline first at `n = 12`

Change only:
- propagator from scalar / transverse to Dirac-Kaehler
- packet state from scalar packet to minimal graded packet

## State Representation

Use minimal graded packet state:
- 0-form component
- 1-form component

Optional:
- 2-form component only if already implemented cleanly

Norm:

\[
\|\Psi\|^2 = \|\psi^{(0)}\|^2 + \|\psi^{(1)}\|^2
\]

Track grade weights over time.

## Propagation Rule

Use first-order Dirac-Kaehler evolution:

\[
\Psi_{t+\Delta t} = \Psi_t - i \Delta t\, D_{DK}\Psi_t
\]

with

\[
D_{DK} = d + \delta
\]

No nonlinear terms.
No feedback terms.
No sector mixing beyond Dirac-Kaehler itself.

## Collision Geometries To Scan

Run only a minimal set.

### 1. Head-On Symmetric Pair

Two identical Dirac-Kaehler packets launched toward each other.

Question:
- pass-through
- reflection
- graded exchange
- transient bound oscillation

### 2. Counter-Propagating Corridor Pair

Two packets in a corridor-like geometry.

Question:
- exclusion-like avoidance
- stabilized channel occupation
- grade-protected passage

### 3. Offset Glancing Collision

Two packets with small transverse offset.

Question:
- deflection
- rotational exchange
- asymmetric grade transfer

### 4. Tight Clustered Pair

Two initially close packets.

Question:
- metastable composite
- rapid dispersal
- grade beating

### 5. Symmetric Triad Collision

Three packets in a symmetric encounter pattern.

Question:
- central trapping
- exchange symmetry breaking
- collective grade oscillation

## Minimal Parameter Set

Keep the scan small.

For each geometry, vary only:
- relative phase: in-phase / out-of-phase
- one packet width class: baseline coherent width
- one amplitude class: baseline amplitude
- one graph type first: regular lattice
- one resolution first: `n = 12`

Optional second pass only if needed:
- `n = 24`

Total first-pass target:
- 5 geometries
- 2 phase relations

= 10 Dirac-Kaehler runs

## Primary Observables

Use persistence-first metrics plus Dirac-Kaehler-specific ones.

Persistence metrics:
- composite lifetime
- binding persistence
- corridor dwell
- coarse basin persistence

Collision metrics:
- minimum separation reached
- post-collision separation trend
- reflection / pass-through classification
- deflection angle proxy
- encounter dwell time

Dirac-Kaehler-specific metrics:
- grade weight evolution
- grade-transfer amplitude
- sign / symmetry exchange behavior
- localized circulation proxy if available

## Output Labels

Each run must receive one collision label:
- pass-through dispersive
- pass-through with grade exchange
- reflective / exclusion-like
- deflective / glancing
- metastable composite
- oscillatory grade-trapped
- unresolved / mixed

And one persistence label:
- no persistence gain
- weak persistence gain
- candidate bounded collision regime

## Positive Signal Criteria

A Stage 23.1 run is interesting only if at least one of the following appears relative to the scalar / transverse baseline:
- longer encounter dwell
- reflection or bounce not seen before
- stable or quasi-stable post-collision boundedness
- reproducible grade-exchange structure
- new collision class requiring a new label

Do not count visually interesting but persistence-null events as major positives.

## Negative Freeze Rule

Freeze Stage 23.1 negative if:
- all runs remain pass-through or rapidly dispersive
- no bounded post-collision regime appears
- grade exchange occurs but does not affect persistence
- collision classes do not differ materially from the scalar / transverse baseline

That would mean:

Dirac-Kaehler changes texture, not persistence.

## Promotion Rule

Promote only a run family that shows:
- nontrivial collision-class change
- persistence increase
- reproducibility under phase flip or small perturbation
- clean behavior at `n = 24`

Only then expand to broader Stage 23 atlas work.

## Deliverables

Create:
- `stage23_1_dk_collision_runs.json`
- `stage23_1_dk_collision_geometry_scan.py`
- `Stage_23_1_DK_Collision_Geometry_Scan_v1.md`

Outputs per run:
- stamped JSON / CSV
- collision trajectory plot
- separation trace
- grade-weight trace
- optional field snapshots before / during / after collision

Summary outputs:
- collision class matrix
- persistence comparison panel
- grade-transfer summary

## Interpretation Boundary

Do not claim:
- fermions
- Pauli exclusion
- matter formation
- emergent binding law

At most, claim:

graded first-order propagation changes or fails to change collision geometry in ways relevant to persistence.
