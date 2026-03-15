# Stage 23.7 - Dirac-Kahler Midpoint-Notch Topology Scan

## Role

You are extending the Phase III Dirac-Kahler collision atlas.
Your task is not to search for persistence or binding.
Your task is to resolve the internal topology of the midpoint detuning window identified in Stage 23.6.

Previous result:
- braid-like exchange is robust across sampled phase space
- a single smeared topology window exists near the midpoint detuning
- this window is not yet geometrically classified

Stage 23.7 is a local notch scan of that window.

## Objective

Determine whether the midpoint corridor contains:
- a narrow topology transition boundary
- a finite smeared band
- multiple micro-regimes
- or numerical sensitivity without geometric meaning

The output must be topology classification only.

No persistence claims.
No stabilization laws.
No kernel modification.

## Fixed Architecture
- periodic 2D complex
- same lattice size and metric conventions as Stage 23.6
- same Dirac-Kahler propagator
- graded state with minimal 0- and 1-form support
- clustered in-phase encounter baseline geometry
- kernel static and symmetric

## Baseline Packet Parameters

Unless varied by the run:
- clustered separation baseline
- equal geometric-mean width
- zero skew
- beta = 0.02
- identical amplitude normalization
- symmetric launch direction

## Scan Axis

Refine only:

`phase_offset_fraction_of_pi`

This parameter measures detuning along the in-phase -> anti-phase corridor.

The suspected smeared region lies near:

`approx 0.5`

Stage 23.7 resolves this region with a narrow symmetric notch.

## Required Run Matrix

Use a minimal 9-run notch scan:
- coarse anchors
- dense midpoint sampling
- symmetric control points

Do not introduce any additional variation.

## Diagnostics Per Run

Compute:
- topology_label
- braid_like_exchange
- transfer_smeared
- fragmented_exchange
- unresolved_mixed_topology
- flow_concentration_index
- grade_exchange_coherence
- effective channel_count
- recirculation_score

## Required Outputs

Per run:
- stamped JSON
- CSV row
- collision trajectory plot
- instantaneous grade-weight trace
- local flow-field visualization

Summary outputs:
- topology classification vs phase notch plot
- flow concentration panel
- coherence vs detuning curve
- phase-corridor topology table

## Interpretation Rules

A result is geometrically meaningful only if:
- topology classification changes reproducibly across adjacent notch samples
- or a narrow stable regime boundary appears
- or a finite smeared band is confirmed

A result is numerical sensitivity if:
- labels fluctuate without consistent structure
- or topology metrics vary without classification change

## Deliverable Note

Write:

`Stage_23_7_DK_Midpoint_Notch_Topology_Scan_v1.md`

It must state clearly:
- whether the midpoint corridor is
- sharp
- broad
- multi-layered
- or inconclusive

No speculation about binding or physical ontology.

## Continuation hint

Likely next step if the notch resolves cleanly:
- Stage 23.8 topology-protection mechanism scan
