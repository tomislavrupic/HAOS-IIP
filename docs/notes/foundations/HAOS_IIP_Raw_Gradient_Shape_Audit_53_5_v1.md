# HAOS-IIP Raw Gradient Shape Audit 53.5

## Purpose

`53.5` is an external sidecar audit over the raw local-gradient boundary left open by `53.3`.

The question is narrower than a new bridge release:

> can the frozen raw recoverability-gradient artifact support a bounded shape-only closure under external source-core exclusion, light smoothing, and per-case amplitude normalization, while keeping the canonical raw bridge row unchanged?

This note does not modify HAOS core, frozen telemetry, bridge thresholds, or the existing `53.0` bridge rows. The canonical raw local-gradient row remains `OPEN`.

## Construction

The audit uses only:

- the frozen raw artifact `data/scalar_kernel_graph_recoverability_gradient_latest.json`;
- the same raw thresholds used in `53.3`;
- external post-processing only.

Three readouts are compared:

1. `canonical_raw_local_gradient`
   The stored `53.3` raw local-gradient row.
2. `windowed_raw_local_gradient`
   Source-core exclusion plus light local smoothing, with no amplitude normalization.
3. `shape_normalized_windowed_raw_local_gradient`
   The same source-core exclusion and smoothing, plus per-case normalization by the mean absolute scaled flux, used only to test refinement-stable profile shape.

The chosen external audit window is:

- minimum radius `r >= 0.28`;
- smoothing window `5`.

## Result

The audit splits the raw boundary into three readable layers:

| Readout | Scope | Status | Min `R^2` | Max `|slope+2|` | Max CV | Max drift |
| --- | --- | --- | --- | --- | --- | --- |
| canonical raw | full observable | `OPEN` | `0.3220` | `1.2651` | `0.3139` | `0.5249` |
| windowed raw | full observable | `OPEN` | `0.9751` | `0.0953` | `0.0436` | `0.5012` |
| shape-normalized windowed raw | shape-only observable | `PASS` | `0.9751` | `0.0953` | `0.0436` | `0.0582` |

The strongest bounded statement is therefore:

> the full raw local-gradient observable remains open, but a source-core-excluded, lightly smoothed, amplitude-normalized shape audit passes the existing raw thresholds as a shape-only proxy.

## Correct Interpretation

This does **not** close full raw inverse-square or power-law recovery.

What it does show is narrower and still useful:

- source-core exclusion and smoothing are enough to recover the raw slope and fit structure;
- the remaining full-observable failure is dominated by refinement-amplitude drift;
- once amplitude is normalized out, the raw profile shape becomes refinement-stable inside the existing raw thresholds.

The correct classification is therefore:

- canonical raw local-gradient: still `OPEN`;
- audited windowed raw: shape improved, amplitude drift still `OPEN`;
- audited shape-normalized windowed raw: bounded `PASS` as a shape-only sidecar proxy.

## What This Supports

`53.5` supports:

- a cleaner decomposition of the raw local-gradient failure mode;
- a bounded shape-only closure on frozen raw artifacts;
- the claim that the raw boundary is not purely a slope-fit problem, but mostly an amplitude-stability problem under refinement.

## What This Does Not Support

This note does not claim:

- full raw inverse-square closure;
- promotion of the shape-normalized read into the canonical bridge table;
- a replacement for the shell-native readout from `53.3`;
- continuum ontology or a gravity claim.

## Authority and Artifacts

- generator: `experiments/physics_bridge/raw_gradient_shape_audit.py`
- CSV: `experiments/physics_bridge/results/raw_gradient_shape_audit.csv`
- JSON: `experiments/physics_bridge/results/raw_gradient_shape_audit.json`
- Markdown summary: `experiments/physics_bridge/results/raw_gradient_shape_audit_summary.md`
- reproduction hook: `python3 examples/quick_reproduce.py --raw-gradient-audit-only`

## Next Honest Move

The next raw-gradient step should stay external:

1. test whether the same shape-only closure survives one larger refinement or one additional kernel-graph family;
2. if it does, probe whether a shell-anchored amplitude correction can close more than shape alone;
3. only then revisit whether the raw inverse-square boundary should split further.
