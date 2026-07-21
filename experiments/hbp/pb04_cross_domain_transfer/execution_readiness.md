# PB-04 Execution Readiness

Audit note (2026-07-11): this file records the pre-run state. A runner and
stored result now exist, but the baseline and HAOS candidate currently execute
the same model path, while the seed-repeat and topology controls do not enact
their declared perturbations. The stored result is not authoritative for bridge
promotion until a new precommitted repair run closes those gaps.

PB-04 remains contract-first, but the source and target anchors are now pinned.

## Frozen Anchors

Source proxy:

- `experiments/biology/gene_network_demo`

Target candidate:

- `DATA/Powergraph/dataset_pf_opf`

## Why this is the right anchor

- The benchmark asks for cross-domain transfer.
- The source proxy is in-repo, deterministic, and already validated.
- The target is real infrastructure data, not a synthetic mirror of the source.
- The pairing keeps the benchmark hard enough that a win would matter.

## Explicit limitation

The biology-line artifact is a synthetic gene-network demo, not biological
tissue. That limitation is intentional and must stay visible.

## Still unfrozen

- transfer protocol
- evaluation freeze
- exact feature restrictions for the transfer model
- runner implementation

## Next implementation gate

Before a runner exists, PB-04 still needs:

1. a transfer protocol manifest
2. a frozen evaluation rule
3. a target-domain split manifest
4. a deterministic holdout runner
