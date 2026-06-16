# PB-04 Cross-Domain Structure Transfer

PB-04 is the harder HBP benchmark branch. It is a frozen precommitment draft
and does not claim any empirical result yet.

## Question

Can HAOS recognize the same recovery structure in completely different domains?

## Source / Target Domains

Source domain:

- biological tissue

Possible holdout domains:

- infrastructure networks
- social networks
- ecological systems

## Task

Learn nothing domain-specific. Use only HAOS structural descriptors, then
attempt to rank recoverability on a different domain.

## Freeze Criteria

- source domain frozen
- target domain frozen
- transfer protocol frozen
- evaluation frozen

## Current Dataset Direction

See [`dataset_selection.md`](dataset_selection.md) for the current biology-line
proxy source and PowerGraph target candidate pairing.

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

- no domain-bridging result yet
- no transfer success claim
- no physical mechanism claim
- no hidden tuning against the holdout
