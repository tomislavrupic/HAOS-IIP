# Stage 23.6 - Dirac-Kahler Phase-Sensitivity and Braid-Exchange Stability Scan

## Purpose

Stage 23.5 established a repeatable braid-like local energy-flow topology in clustered Dirac-Kahler (DK) encounters.
Stage 23.6 isolates the phase axis as the primary control variable and tests whether:
- braid-like exchange is structurally robust across phase space,
- there exists a phase corridor that sharpens coherence into a new regime,
- or the braid family collapses into mixed / dispersive topology under detuning.

This is a topology-refinement probe, not a persistence or binding claim.

## Architecture (frozen)

Use the existing Phase III DK substrate:
- operator: Dirac-Kahler propagator
  `D_DK = d + delta`
- graded packet state with 0-form + 1-form support
- kernel:
  - static Gaussian interaction kernel
  - no temporal memory
  - no adaptive deformation
- geometry:
  - 2D periodic complex
  - clustered baseline configuration

No new stabilizers, mediator fields, or constraint sectors.

## Probe principle

Hold all clustered-family parameters fixed and vary only:
- relative phase offset
- minor phase-gradient perturbation
- weak grade-coupling modulation (single-run sensitivity check)

We are mapping:

topology response surface in phase space.

## Clustered baseline (reference)
- mean width sigma approx 0.08
- symmetric amplitudes
- separation: tight clustered regime
- in-phase condition previously produced the strongest braid-like exchange

## Phase-scan axes
1. global phase offset Delta phi
2. local phase skew (small gradient across packet support)
3. weak grade-coupling modulation beta (single exploratory axis)

## Expected classification outputs

Per run compute:
- topology label
- braid_like_exchange
- transfer_smeared
- mixed_unresolved
- dispersive_pass
- flow_concentration_index
- grade_exchange_coherence
- channel_count
- loop_count
- recirculation_score
- phase_alignment_metric

No persistence promotion logic.

## Minimal run design

Use a 9-run matrix:

Phase corridor sweep
- in-phase baseline
- quarter-cycle offset
- half-cycle offset
- three-quarter-cycle offset
- full anti-phase

Fine sensitivity
- slight positive phase skew
- slight negative phase skew

Weak coupling check
- beta = 0.01
- beta = 0.02

All other parameters frozen.

## Numerical safety rules
- bounded grade-coupling modulation
- phase perturbations <= 10 percent spatial gradient
- fixed time-step and grid resolution from Stage 23.5
- terminate run if:
  - packet norm loss > threshold
  - numerical instability indicator triggers

## Deliverables

Per run:
- stamped JSON / CSV
- topology trajectory plot
- flow concentration trace
- grade-exchange trace

Summary:
- topology-class matrix vs phase
- coherence heatmap
- braid-family robustness panel

Note target:

`Stage_23_6_DK_Phase_Sensitivity_Braid_Exchange_v1.md`

## Interpretation rule

Stage 23.6 is positive if:
- braid-like exchange persists across a finite phase corridor
- or topology sharpens into a new repeatable family

Stage 23.6 is negative if:
- braid structure collapses monotonically with detuning
- or topology fragments into mixed / dispersive classes

## Clean objective

Determine whether:

DK braid-exchange is a phase-protected topology,
or merely a fine-tuned clustered coincidence.
