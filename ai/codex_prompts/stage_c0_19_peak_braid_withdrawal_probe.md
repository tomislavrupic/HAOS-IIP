# Stage C0.19 - Peak-Braid Withdrawal Probe

## Purpose

Push the late C0 bridge line toward a fairer locking test.

Stages C0.17 and C0.18 showed that bridge-local latches did not retain topology after selector withdrawal, but those runs withdrew from fixed arming cutoffs that did not always coincide with the real bridge braid window.

Stage C0.19 corrects that.

It first arms each bridge family under the exact-match harmonic selector, identifies the strongest or longest recovered `braid_like_exchange` window, and only then applies selector withdrawal.

## Frozen Architecture

- same signed C0-to-DK bridge scaffold
- same representative families
- same no-distance combinatorial primitives
- same harmonic-address dressing
- no metric distance
- no Gaussian kernel
- no external mediator field

## Run Matrix

Exactly 9 runs:

- 3 representative families
- 3 withdrawal protocols

Protocols:

1. `peak_abrupt_release`
2. `peak_frozen_combined_release`
3. `peak_ramped_combined_release`

## Required Logic

For each family:

1. arm the bridge under exact harmonic match
2. detect the longest recovered `braid_like_exchange` window
3. extract the terminal state of that window
4. withdraw the selector using one of the three protocols
5. classify retention

If no braid window appears, report that explicitly rather than forcing a retention claim.

## Observables

- `armed_topology_class`
- `armed_window_start_step`
- `armed_window_end_step`
- `armed_window_length`
- `armed_peak_loop_score`
- `post_withdrawal_topology_class`
- `withdrawal_braid_dwell_time`
- `loop_retention_score`
- `path_reuse_retention_score`
- `selector_off_recurrence_indicator`
- `locking_label`

## Classification

- `locked_after_withdrawal`
- `quasi_locked_retention`
- `topology_afterglow`
- `selector_dependent_collapse`
- `no_armed_braid_window`

## Honest Boundary

If the clustered bridge family still collapses immediately after withdrawal when started from a real braid window, the late C0 locking null becomes substantially stronger.
