# Stage C0.17 - Bridge Cycle-Latch Selector Withdrawal Probe

## Purpose

Push the late C0 line toward a stronger locking test on the derived C0-to-DK bridge scaffold.

This stage uses the signed combinatorial Dirac-Kahler bridge candidate from the C0 bridge validation note and asks:

`Can a bounded bridge-local latch built from loop occupancy or path reuse retain topology after the harmonic selector is withdrawn?`

The stage remains strictly no-distance:

- no Euclidean coordinates
- no Gaussian kernel
- no external stabilizer field
- no imported continuous geometry

## Frozen architecture

Keep fixed:

- representative families:
  - `clustered_composite_anchor`
  - `counter_propagating_corridor`
  - `phase_ordered_symmetric_triad`
- signed C0-to-DK bridge complex
- deterministic harmonic-address dressing
- derived `D_H` evolution on the bridge complex
- graph size and combinatorial motif support

Only add one bounded latch rule at a time.

## Latch rules

### Lock A - Cycle occupancy latch

During an arming window under exact-match harmonic compatibility, accumulate 2-cell loop occupancy and feed a bounded temporary weight back onto the active triangle boundary.

Interpretation:
- repeated closed loop support leaves a local cycle-weight trace

### Lock B - Edge-path reuse latch

During arming, accumulate edge support on the dominant exchange route and feed a bounded temporary weight back onto the re-used boundary edges.

Interpretation:
- repeated exchange prefers the route it actually used

### Lock C - Combined cycle-plus-path latch

Apply both bounded latches together with reduced strength so neither dominates by construction.

Interpretation:
- tests whether bridge-local loop support and path support need to act together before retention appears

## Protocol

For each run:

1. arm the family under exact harmonic match for a fixed arming window
2. accumulate the selected latch variable only during arming
3. freeze the latch at its bounded final amplitude
4. withdraw the selector by switching to the incompatible detuning sector
5. continue the same bridge evolution under the latched operator
6. classify whether topology survives or collapses

## Primary observables

Per run compute:

- `armed_topology_class`
- `post_withdrawal_topology_class`
- `withdrawal_braid_dwell_time`
- `loop_retention_score`
- `path_reuse_retention_score`
- `latch_duty_cycle`
- `flow_concentration_index_post`
- `grade_exchange_coherence_post`
- `selector_off_recurrence_indicator`

Derived labels:

- `locked_after_withdrawal`
- `quasi_locked_retention`
- `topology_afterglow`
- `selector_dependent_collapse`
- `no_locking_after_withdrawal`

## Run matrix

Use exactly 9 runs:

- 3 families
- 3 latch rules

## Positive criterion

The stage is positive if at least one latch family:

- keeps the clustered bridge family in `braid_like_exchange` after selector withdrawal
or
- produces a reproducible `quasi_locked_retention` band with nonzero dwell in more than one family

## Negative criterion

The stage is negative if:

- all armed braid cases collapse immediately after selector withdrawal
- all dwell times remain effectively zero
- latch scores activate but do not change class outcomes

## Required outputs

- stamped JSON
- stamped CSV
- topology trace plots
- latch trace plots
- clustered withdrawal summary panel
- note target:
  - `Stage_C0_17_Bridge_Cycle_Latch_Withdrawal_Probe_v1.md`
