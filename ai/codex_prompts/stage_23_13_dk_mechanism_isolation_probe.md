# Stage 23.13 -- DK Mechanism Isolation Probe

## Role

This stage is the final Phase III mechanism-level probe on the frozen HAOS (Harmonic Address Operating System) / IIP (Interaction Invariance Physics) Dirac-Kahler collision sector.

It is not a new atlas pass.
It is a finite ablation / isolation experiment.

## Purpose

Phase III already established:

- braid-like exchange topology
- smeared transfer topology
- a finite phase-band structure
- no persistence promotion
- no metastable transport family

The remaining question is:

which smallest dynamical ingredient actually controls the braid / smear split?

## Core question

Test whether braid-like exchange is primarily controlled by:

1. grade transfer itself
2. the phase-alignment corridor
3. clustered geometry / overlap
4. operator cross-coupling between the 0-form and 1-form sectors

## Frozen constraints

Keep fixed:

- DK operator sector
- graph family
- Gaussian kernel class
- no memory
- no adaptive kernel
- no selector logic
- no combinatorial austerity element
- same late-Phase III clustered seed convention

Do not add new stabilizers, fields, or persistence logic.

## Reference seed

Use the strongest late-Phase III clustered braid-supporting seed as the control family:

- symmetric clustered pair
- width ratio `1.0`
- separation `0.08`
- phase offset `0.625 pi`
- weak graded coupling `beta = 0.02`

## Branches

Run four branches with three runs each:

### A. Grade-transfer suppression

- control
- weak bounded suppression of net cross-grade redistribution
- moderate bounded suppression of net cross-grade redistribution

Implementation note:
because the frozen DK branch does not expose a literal transfer-off switch without changing operator class, use bounded blockwise grade-fraction anchoring as a transfer-suppression surrogate while leaving the spatial propagator and timestep fixed.

### B. Phase-corridor flattening

- control
- weak flattening of the braid-supporting phase offset toward the midpoint band
- moderate flattening to the midpoint band

### C. Geometry-overlap weakening

- control
- weak overlap weakening near the braid boundary
- moderate overlap weakening near the braid boundary

Implement this only through small clustered separation / width adjustments.

### D. Operator cross-coupling suppression

- control
- weak suppression of the `0 <-> 1` block coupling inside the DK operator
- moderate suppression of the same block coupling

Preserve the rest of the operator and the same weak graded modulation term.

## Observables

For each run record:

- `topology_class`
- `braid_like_exchange`
- `transfer_smeared`
- `localized_encounter`
- `unresolved_mixed`
- `braid_survival_time`
- `grade_exchange_coherence`
- `flow_concentration_index`
- `separation_oscillation_indicator`
- `topology_return_error`
- `mechanism_sensitivity_score`

Also save:

- topology trace
- braid-indicator trace
- coherence trace
- ablation-strength log

## Branch interpretation

Treat a branch as discriminative only if:

- braid collapses reproducibly under that ablation
- while the branch control remains braid-like
- and the collapse stays in a bounded dynamical regime rather than dissolving into trivial numerical failure

## Mechanism labels

Use only:

- `grade_transfer_primary`
- `phase_corridor_primary`
- `geometry_overlap_primary`
- `operator_coupling_primary`
- `coupled_mechanism_no_single_primary`
- `no_clear_mechanism_signal`

## Closure rule

Phase III mechanism isolation succeeds only if:

- one branch dominates the collapse structure
- or a reproducible coupled mechanism appears across more than one non-identical branch

At most conclude:

the DK braid / smear split is primarily controlled by `[X]`, with `[Y]` acting as secondary support.

Do not claim:

- persistence law
- bound-state mechanism
- emergent particles
- confinement

## Outputs

Per run:

- stamped JSON / CSV payload
- trace panel

Summary outputs:

- mechanism branch comparison matrix
- braid survival vs ablation panel
- coherence collapse panel
- short mechanism conclusion note

Write the note to:

`Stage_III_C4_Mechanism_Isolation_Probe_v1.md`
