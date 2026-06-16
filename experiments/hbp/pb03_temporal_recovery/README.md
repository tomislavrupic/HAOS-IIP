# PB-03 Temporal Recovery Under Damage

PB-03 is the next HBP benchmark branch. It is a frozen precommitment draft and
does not claim any empirical result yet.

## Question

Can HAOS predict future recoverability of a system better than topology-only
methods?

## Candidate Domains

- temporal networks
- infrastructure evolution
- biological development sequences
- versioned software dependency graphs
- dynamic communication networks

## Task

Observe a network at time `t0`, introduce damage at `t1`, and predict recovery
potential at `t2`.

## Freeze Criteria

- dataset frozen
- split frozen
- metrics frozen
- baselines frozen

## Current Dataset Direction

See [`dataset_selection.md`](dataset_selection.md) for the current PowerGraph
cascade candidate family.

See [`source_manifest.json`](source_manifest.json) and
[`execution_readiness.md`](execution_readiness.md) for the current freeze
boundary.

See [`split_manifest.json`](split_manifest.json) and
[`execution_contract.json`](execution_contract.json) for the first frozen
implementation preconditions.

See [`metrics_manifest.json`](metrics_manifest.json) and
[`baselines_manifest.json`](baselines_manifest.json) for the frozen evaluation
choices.

See [`data_schema.json`](data_schema.json) and
[`control_manifest.json`](control_manifest.json) for the frozen data and
control rules.

## Non-Claims

- no empirical validation
- no physical mechanism claim
- no domain transfer claim
- no hidden tuning after holdout inspection
