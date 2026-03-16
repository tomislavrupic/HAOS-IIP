# Stage C0.15 - Harmonic-Address Hysteresis and Return-Error Probe

Master prompt (HAOS-IIP combinatorial no-distance branch)

## Purpose

Stages C0.10-C0.14 established weak harmonic-address topology selection, a sharp detuning threshold, and a locking null result. Stage C0.15 asks the sharper irreversibility question:

Can the harmonic-address branch fail to retrace when the detuning/forcing schedule is reversed?

This stage does not add latch rules or metric structure. It tests path dependence directly.

## Frozen architecture

Keep fixed:
- no coordinate embedding
- no Euclidean distance
- same clustered Dirac-Kahler projected-transverse seed from C0.10-C0.14
- same graph size, operator sector, packet widths, amplitudes, and update law
- same harmonic-address compatibility rule
- no memory kernel
- no new locking bonus

Only vary:
- a forward/reverse detuning schedule
- one bounded combinatorial perturbation protocol at a time

## Run design

Use a 9-run matrix:
- 3 seed corridors
- 3 perturbation protocols

Seed corridors:
1. exact-match corridor: starts and ends in the braid-selected sector
2. threshold corridor: starts and ends near the protected-to-smeared boundary
3. fragile corridor: starts and ends in the incompatible/smeared sector

Perturbation protocols:
1. plain detuning cycle
2. detuning cycle with mild degree-skew bias
3. detuning cycle with local motif bias

Each run uses a symmetric forward/reverse schedule. Example:
- exact-match corridor: `0.00 -> 0.25 -> 0.50 -> 0.75 -> 0.50 -> 0.25 -> 0.00`
- threshold corridor: `0.625 -> 0.6875 -> 0.75 -> 0.8125 -> 0.75 -> 0.6875 -> 0.625`
- fragile corridor: `1.00 -> 0.875 -> 0.75 -> 0.625 -> 0.75 -> 0.875 -> 1.00`

## Required observables

Per run compute:
- initial_topology_class
- final_topology_class
- topology_return_error
- flow_return_error
- selectivity_return_error
- coherence_return_error
- hysteresis_area_flow_concentration
- hysteresis_area_selectivity
- hysteresis_area_coherence
- irreversible_witness
- retrace_class

Derived classes:
- `reversible_retrace`
- `hysteretic_retrace`
- `irreversible_witness`

## Classification rule

A run is `reversible_retrace` if:
- initial and final topology classes match
- return errors stay below tolerance
- hysteresis areas stay below tolerance

A run is `hysteretic_retrace` if:
- topology returns to the initial class
- but one or more hysteresis areas remain visibly nonzero

A run is `irreversible_witness` if:
- the final topology class differs from the initial class
- or return error / hysteresis exceeds the irreversibility threshold

## Required outputs

Produce:
- stamped JSON
- stamped CSV
- per-run forward/reverse topology trace
- per-run forward/reverse selectivity trace
- hysteresis panel
- return-error summary panel
- short note

Note target:

`Stage_C0_15_Harmonic_Hysteresis_Return_Error_Probe_v1.md`

## Interpretation boundary

Do not claim:
- physical time emergence
- thermodynamic irreversibility
- bound states
- topology locking

At most claim:

the harmonic-address no-distance branch does or does not show path-dependent return failure under forward/reverse combinatorial detuning schedules.
