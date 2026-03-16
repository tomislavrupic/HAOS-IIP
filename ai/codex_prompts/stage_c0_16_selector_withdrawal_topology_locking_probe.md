# Stage C0.16 - Selector-Withdrawal Topology Retention Probe

Master prompt (HAOS-IIP combinatorial no-distance branch)

## Purpose

Stage C0.14 showed that bounded latch rules do not create topology locking. Stage C0.16 asks the cleaner locking question:

If topology is first selected by harmonic-address compatibility and the selector is then removed, does the topology persist anyway?

This is a selector-withdrawal test, not a stronger forcing test.

## Frozen architecture

Keep fixed:
- no coordinate embedding
- no Euclidean distance
- same clustered Dirac-Kahler projected-transverse seed from C0.10-C0.14
- same packet width, amplitude, graph size, and update law
- same no-distance combinatorial operator family
- no memory kernel
- no new latch or reinforcement law

Only vary:
- the seed family used during the arming phase
- the bounded perturbation protocol applied after selector withdrawal

## Run design

Use a 9-run matrix:
- 3 armed seed families
- 3 withdrawal protocols

Armed seed families:
1. exact-match braid-selected seed
2. near-threshold protected seed
3. fragile incompatible control

Withdrawal protocols:
1. selector off, no extra perturbation
2. selector off plus mild degree-skew perturbation
3. selector off plus local motif perturbation

Each run has two phases:
1. arming phase under harmonic-address selection
2. withdrawal phase under the base combinatorial operator with the selector removed

## Required observables

Per run compute:
- armed_topology_class
- post_withdrawal_topology_class
- withdrawal_braid_dwell_time
- withdrawal_class_retention
- post_withdrawal_selectivity
- post_withdrawal_coherence
- post_withdrawal_flow_concentration
- selector_off_recurrence
- locking_after_withdrawal_label

Derived labels:
- `locked_after_withdrawal`
- `quasi_locked_retention`
- `topology_afterglow`
- `selector_dependent_collapse`
- `no_locking_after_withdrawal`

## Classification rule

A run counts as `locked_after_withdrawal` if:
- the armed phase reaches a braid-selected topology
- the withdrawal phase retains that topology through most of the selector-off window
- coherence and flow concentration remain inside the armed-phase class envelope

A run counts as `topology_afterglow` if:
- the withdrawal phase briefly retains the selected topology
- but the topology decays before the end of the selector-off window

A run counts as `selector_dependent_collapse` if:
- the armed phase reaches a selected topology
- and the topology collapses promptly once the selector is removed

## Required outputs

Produce:
- stamped JSON
- stamped CSV
- per-run topology trace
- per-run selector-off selectivity trace
- withdrawal dwell summary panel
- locking-after-withdrawal class matrix
- short note

Note target:

`Stage_C0_16_Selector_Withdrawal_Topology_Retention_Probe_v1.md`

## Interpretation boundary

Do not claim:
- true bound states
- confinement
- metric trapping
- broad irreversibility

At most claim:

the selected harmonic-address topology does or does not persist after the selector is removed on the no-distance combinatorial branch.
