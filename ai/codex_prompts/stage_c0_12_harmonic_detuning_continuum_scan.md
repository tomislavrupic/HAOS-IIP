# Stage C0.12 -- Harmonic Detuning Continuum Scan

HAOS-IIP Pre-Geometry Atlas / Combinatorial Branch

## Purpose

Stage C0.12 asks whether harmonic-address compatibility behaves like:

- a continuous topology selector
- a thresholded regime switch
- or a narrow resonance ladder

C0.10 established weak proto-selection on the clustered DK seed:

- exact harmonic matches favored braid-like exchange
- near matches relaxed toward smeared transfer
- incompatible assignments destroyed topology coherence

This stage now scans a dense detuning ladder and measures topology-class stability under harmonic detuning.

## Frozen baseline

Keep fixed:

- the same clustered DK seed used in C0.10
- no memory kernel
- no stabilizer or curvature feedback
- identical combinatorial kernel weights
- same graph size and motif placement

Only the harmonic-address assignment is modified.

## Control parameter

Define a scalar detuning parameter:

`delta in [0, 1]`

Interpretation:

- `delta = 0` exact address match
- `delta = 0.25` mild detuning
- `delta = 0.5` mid-detuning corridor
- `delta = 0.75` strong detuning
- `delta = 1` maximally incompatible assignment

Detuning must be implemented by harmonic-address relabeling or harmonic compatibility weights only.
No geometric deformation or amplitude rescaling is allowed.

## Required observables

Per run compute:

- `topology_class`
- `flow_concentration_index`
- `address_selectivity_index`
- `address_sorting_score`
- `topology_survival_time`
- `transfer_asymmetry`

Persistence metrics may be recorded, but they must not drive classification.

## Run logic

1. initialize the same clustered packet pair
2. apply detuning `delta`
3. evolve with the frozen DK propagator
4. detect topology class and topology-vs-time trace
5. record classification and detuning-sensitive indices

## Summary outputs

Produce:

1. topology phase diagram
2. selectivity decay curve
3. sorting degradation curve
4. braid survival panel
5. continuity assessment

Global classification must be one of:

- `smooth_continuum_decay`
- `sharp_transition_threshold`
- `multi_plateau_resonance_ladder`
- `incoherent_response`

## Note target

`Stage_C0_12_Harmonic_Detuning_Continuum_Scan_v1.md`
