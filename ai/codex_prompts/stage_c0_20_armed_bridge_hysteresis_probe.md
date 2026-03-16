# Stage C0.20 - Armed-Bridge Hysteresis Probe

## Purpose

Push the late C0 bridge line toward the irreversibility side once locking has been tested fairly.

This stage starts from the real armed bridge braid window recovered in Stage C0.19 and asks:

`Does the bridge retrace when harmonic detuning is ramped up and then back down, or does it keep a return error even after the selector is restored?`

## Frozen Architecture

- same signed C0-to-DK bridge scaffold
- same no-distance combinatorial primitives
- same representative families
- same harmonic-address dressing
- no metric distance
- no Gaussian kernel

## Run Matrix

Exactly 9 runs:

- 3 representative families
- 3 hysteresis protocols

Protocols:

1. `plain_forward_reverse`
2. `frozen_combined_forward_reverse`
3. `frozen_hold_forward_reverse`

## Required Logic

For each family:

1. arm under exact harmonic match
2. recover the strongest braid window if one exists
3. start the hysteresis cycle from that armed state
4. ramp detuning `0 -> 1 -> 0`
5. optionally apply a frozen bridge-local latch and a maximum-detuning hold
6. measure return error at restored exact match

## Observables

- `armed_topology_class`
- `armed_window_length`
- `final_return_topology_class`
- `topology_return_error`
- `return_state_overlap`
- `hysteresis_area_flow`
- `hysteresis_area_loop`
- `braid_fraction_during_cycle`
- `irreversibility_label`

## Classification

- `irreversible_witness`
- `hysteretic_retrace`
- `reversible_retrace`
- `no_armed_braid_window`

## Honest Boundary

If the bridge still fails to retrace from a real armed braid state after the selector returns to exact match, then the late C0 line gains a stronger irreversibility witness even if locking remains absent.
