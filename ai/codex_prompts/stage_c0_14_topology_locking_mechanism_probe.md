# Stage C0.14 -- Topology-Locking Mechanism Probe

HAOS / IIP combinatorial no-distance branch

## Purpose

Stage C0.14 asks whether a purely combinatorial locking rule can convert topology selection into a longer-lived exchange regime without introducing metric trapping or external stabilizers.

Previous C0 stages established:

- exchange-like topology can arise without coordinate distance
- harmonic-address compatibility can bias the topology class
- a detuning corridor separates braid-preserving and smeared-transfer sectors
- selection exists, but persistent locking has not yet appeared

This stage tests three bounded local locking rules:

1. exchange-path reinforcement
2. harmonic-address phase latch
3. motif occupancy latch

## Frozen constraints

Keep fixed:

- no coordinate embedding
- no Euclidean or Gaussian distance
- same clustered combinatorial packet seed family
- same graph size class and operator sector
- same harmonic-address weighting infrastructure from C0.10-C0.13
- no mediator field
- no delayed memory kernel
- no external persistence bonus

Only one local locking rule is added per run.

## Representative families

Use three compact families:

1. exact-match braid-protected clustered seed
2. near-boundary detuned smeared seed
3. fragile randomized harmonic-address control

## Required observables

Per run compute:

- `topology_class`
- `braid_like_exchange`
- `transfer_smeared`
- `locked_braid`
- `quasi_bounded_exchange`
- `unresolved`
- `topology_survival_time`
- `recurrence_indicator`
- `local_path_reuse_score`
- `harmonic_latch_duty_cycle`
- `motif_anchor_score`
- `flow_concentration_index`
- `exchange_coherence`
- `transport_span`
- `bounded_dwell_proxy`

Each run also receives a derived locking label:

- `no_locking_effect`
- `topology_stabilized_no_capture`
- `quasi_bounded_exchange`
- `motif_localized_lock`
- `address_latched_exchange`

## Minimal design

Run the 3 x 3 matrix:

- 3 locking rules
- 3 representative families

Each run is compared internally to a no-lock baseline for the same family.

## Outputs

Produce:

- stamped JSON and CSV
- per-run topology trace
- per-run latch or lock score trace
- per-run flow concentration trace
- locking-rule x family class matrix
- topology survival comparison panel
- latch score / anchor score summary
- short interpretation note

## Note target

`Stage_C0_14_Topology_Locking_Mechanism_Probe_v1.md`
