# Stage 24.1 / IV-A1 -- Smeared-Sector Observable Ledger

Context:
HAOS = Harmonic Address Operating System.
IIP = Interaction Invariance Physics.

Purpose

Define and validate the minimal observable basis of the surviving smeared-transfer sector on the frozen clustered Dirac-Kahler branch.

This is the opening stage of Phase IV.
It is not a braid-rescue stage.
It is not a family-expansion stage.
It is not a new mechanism hunt.

Its purpose is to make the surviving smeared sector exact enough to support later reduction, quasi-invariant testing, and bounded extension analysis.

Frozen scientific basis

Treat the following Phase III outputs as authoritative:
- `23.11 / III-C2` -- valid clustered-sector two-state surrogate
- `23.12 / III-C3` -- one-sided smeared-dominant closed phase diagram
- `23.13 / III-C4` -- phase corridor primary, geometry overlap and `0 <-> 1` cross-coupling secondary
- `23.14 / III-C5` -- final consolidation and freeze

Preserve the following negatives:
- no family-wide intrinsic braid mechanism
- no stable braid phase
- no localized encounter phase
- no grade-transfer-primary interpretation

Primary scientific question

What is the smallest operational observable ledger that characterizes the closed smeared-transfer regime on the frozen clustered DK branch?

Restricted object of study

The object is the surviving smeared sector identified in Stage 23.12:
- clustered DK branch only
- one-sided smeared-transfer closure only
- stable connected region at `width_ratio = 1.35`
- sampled phase corridor `0.375 -> 0.575`
- internal weak-coupling stability probe `beta in {0.01, 0.02}` as already used in the closure pass

Do not expand to non-clustered families.
Do not search for new braid windows.
Do not change operator class.

Allowed work

You may:
- extract and reconcile observables from the frozen Stage 23.11 through 23.14 outputs,
- rerun only bookkeeping-level measurements on already-closed smeared cells if a needed observable is absent from the frozen outputs,
- derive compact secondary quantities from already-accepted primitive observables.

You may not:
- run exploratory family scans,
- introduce new memory laws or adaptive kernels,
- alter DK grading rules,
- alter kernel functional form,
- or broaden the phase-width lattice beyond what Phase III already closed.

Required observable candidates

Start from the following candidate ledger:
- `phase_corridor_position`
- `phase_corridor_width`
- `topology_survival_time`
- `transfer_asymmetry`
- `refinement_stability_flag`
- `reverse_stability_flag`
- `coherence_decay_rate`
- `flow_concentration_index`
- `topology_return_error`
- `smeared_regime_membership_flag`

If two quantities are redundant on the closed smeared region, collapse them.
If a quantity fails to survive the closure pass in a reproducible way, demote it.

The deliverable is not a long list.
The deliverable is the minimal sufficient observable basis.

Required outputs

Produce:
1. a note:
   - `Stage_IV_A1_Smeared_Sector_Observable_Ledger_v1.md`
2. a script:
   - `stage24_1_smeared_sector_observable_ledger.py`
3. a runsheet:
   - `stage24_1_smeared_sector_observable_ledger_runs.json`
4. a stamped JSON summary:
   - `*_stage24_1_smeared_sector_observable_ledger.json`
5. a stamped CSV ledger:
   - `*_stage24_1_smeared_sector_observable_ledger.csv`
6. summary visuals:
   - `*_stage24_1_smeared_sector_observable_table.png`
   - `*_stage24_1_smeared_sector_correlation_panel.png`
   - `*_stage24_1_smeared_sector_stability_panel.png`

Core task

Build a minimal observable architecture for the closed smeared regime by answering:
1. which observables are genuinely needed to characterize smeared-sector membership,
2. which observables are needed to characterize smeared-sector stability,
3. which observables are redundant or purely descriptive,
4. which observables survive refinement and reverse checks cleanly enough to carry into Phase IV reduction work.

Suggested analysis discipline

Use a three-layer classification:
- `primary_observable`
- `secondary_observable`
- `redundant_or_derived`

An observable counts as `primary_observable` only if it:
- remains well-defined across the closed smeared region,
- participates directly in stable regime characterization,
- and is not reducible to another kept quantity without losing Phase III closure information.

Validation criteria

The Stage 24.1 ledger is valid only if:
- every retained primary observable is supported across the closed smeared region,
- the retained basis is smaller than the candidate list,
- the basis can distinguish smeared-sector membership from transient mixed cells on the frozen closure lattice,
- and no retained observable contradicts the frozen Phase III ledger.

Required final read

The final note must end with a statement of the form:

On the frozen clustered DK branch, the closed smeared-transfer regime is minimally characterized by [X, Y, Z, ...], while [A, B, ...] remain secondary or redundant descriptors.

Interpretation boundary

Do not claim:
- a universal transport law
- a conserved quantity
- a continuum correspondence
- a particle interpretation
- or a family-wide mechanism

At most claim:

the surviving smeared-transfer sector admits a minimal operational observable basis on the frozen clustered DK branch.
