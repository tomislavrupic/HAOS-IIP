# PB-03 Execution Readiness

PB-03 remains contract-first, but the source family is now pinned.

## Frozen Source Anchor

- `DATA/Powergraph/dataset_cascades`
- selected cases:
  - `ieee24`
  - `ieee39`
  - `ieee118`
  - `uk`

## Why this is the right anchor

- The benchmark question is temporal recovery under damage.
- The cascades dataset exposes explicit disturbance structure.
- This keeps the benchmark harder than a static graph task.
- It avoids reopening the already frozen PB-02 power-grid holdout path.

## Still unfrozen

- split assignment
- evaluation metrics
- baseline list
- runner implementation

## Next implementation gate

Before a runner exists, PB-03 still needs:

1. a split manifest
2. metric freeze
3. baseline freeze
4. deterministic holdout runner

