# Stage 23.8 - Dirac-Kahler Topology-Protection Mechanism Scan

HAOS (Harmonic Address Operating System) · IIP (Interaction Invariance Physics)
Pre-Geometry Atlas · Phase III continuation

## Purpose

Stages 23.5–23.7 established a reproducible braid-like exchange topology in clustered Dirac-Kahler encounters, together with a broad smeared midpoint band in phase space.

Stage 23.8 asks a new question:

Which structural mechanisms preserve the braid family, and which ones destabilize it into the smeared band?

This is a topology-protection scan, not a persistence or binding scan.

The goal is to identify whether the braid and smeared families are controlled by:
- internal grade-balance structure
- phase-order protection
- local symmetry of the clustered geometry
- or weak anisotropy / mismatch breaking

## Frozen architecture constraints

Keep fixed:
- DK operator sector
  `D_DK = d + delta`
- 2D periodic complex
- static Gaussian kernel
- minimal graded packet state with 0- and 1-form support
- clustered encounter family inherited from 23.5–23.7
- no substrate memory
- no mediator field
- no adaptive kernel law
- no persistence promotion logic

This stage is about topology class robustness, not composite stabilization.

## Structural hypothesis

The braid-like exchange family may be protected by one or more of:
1. grade-balance protection
   balanced 0-/1-form participation keeps the exchange structured
2. phase-order protection
   braid exchange exists only inside a finite phase corridor
3. geometric symmetry protection
   width/separation symmetry stabilizes the flow class
4. anisotropy breaking
   weak directional bias may selectively destroy braid coherence

## Branches to test

Use four mechanism branches, each isolated.

### Branch A - Grade-balance perturbation
Vary the initial ratio of 0-form to 1-form support while holding geometry and phase fixed.

Purpose:
- test whether braid topology depends on balanced grade participation

Examples:
- balanced baseline
- 0-form enhanced
- 1-form enhanced

### Branch B - Phase-band edge protection
Probe the braid corridor just inside and just outside the known smeared band.

Purpose:
- test whether topology changes sharply or smoothly at the band edge

Examples:
- lower braid-side anchor
- midpoint-band sample
- upper braid-recovery anchor

### Branch C - Width/symmetry breaking
Introduce small width asymmetry around the braid baseline.

Purpose:
- test whether braid topology is geometrically protected or only symmetry-tuned

Examples:
- symmetric baseline
- A wider
- B wider

### Branch D - Weak anisotropy breaking
Apply a mild kernel or launch-direction anisotropy without changing operator class.

Purpose:
- test whether braid topology survives directional bias

Examples:
- isotropic baseline
- weak x-stretch
- weak y-stretch

## Minimal run design

Keep it tight:
- 4 branches
- 3 runs per branch

Total: 12 runs

If you want a cheaper first pass, run Branches A and B only first.

## Primary observables

For each run compute:
- topology label
- braid_like_exchange
- transfer_smeared
- fragmented_exchange
- unresolved_mixed_topology
- flow concentration index
- grade exchange coherence
- channel count
- recirculation score
- grade asymmetry index
- phase alignment metric

## Derived protection labels

Assign one protection read per branch:
- protected
- edge-sensitive
- symmetry-sensitive
- anisotropy-sensitive
- unprotected / fragile

## Classification logic

A mechanism counts as topology-protecting if:
- the braid label survives across the branch perturbation set
- flow concentration remains high
- coherence remains within the braid envelope
- no collapse into smeared or fragmented class occurs

A mechanism counts as topology-breaking if:
- the braid class collapses systematically
- coherence drops below the braid corridor
- label transition is reproducible across the branch

## Output requirements

Per run:
- stamped JSON / CSV
- topology trajectory plot
- coherence trace
- flow-concentration trace

Summary outputs:
1. branch-by-branch topology protection matrix
2. braid-vs-smeared branch comparison panel
3. coherence robustness heatmap
4. short interpretation note

Note target:

`Stage_23_8_DK_Topology_Protection_Mechanisms_v1.md`

## Freeze logic

Stage 23.8 is positive if at least one protection mechanism is identified clearly.

That means:
- one branch preserves braid robustly
- or one branch reliably breaks braid into the smeared family

Stage 23.8 is negative / inconclusive if:
- all branches mix labels without stable structure
- no reproducible protection or breaking rule appears

## Interpretation boundary

Do not claim:
- topological matter
- bound states
- confinement
- emergent particles

At most claim:

the DK clustered collision family exhibits or fails to exhibit reproducible topology-protection structure under controlled local perturbations.

## Suggested filenames
- `stage_23_8_dk_topology_protection_mechanisms.md`
- `stage23_8_dk_topology_protection_runs.json`
- `stage23_8_dk_topology_protection_scan.py`
- `Stage_23_8_DK_Topology_Protection_Mechanisms_v1.md`

## Why this stage matters

23.5–23.7 showed that the braid family is real.
23.8 asks the next honest question:

Is the braid family an accident of tuning, or does it have a real protection structure?

That is the right next move.
