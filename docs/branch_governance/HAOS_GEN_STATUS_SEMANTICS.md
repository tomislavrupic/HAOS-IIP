# HAOS-GEN Evidence Status Semantics

HAOS-GEN separates three questions:

1. Did the experiment execute correctly?
2. Was the tested claim supported?
3. Is the research line open, terminal, or awaiting replication?

These questions must not share one red/green flag.

## Display contract

| Tone | Meaning |
|---|---|
| Green | Valid bounded positive or verified infrastructure |
| Blue | Verified terminal negative; a boundary was established |
| Amber | Verified open negative; tested advantage not supported |
| Red | Invalid execution, broken evidence contract, failed reproducibility, or corrupt artifact |
| Gray | Unclassified or not evaluated |

`DEGRADATION` may remain red at the individual intervention-row level because
that specific modification damaged admitted structure. It does not make the
experiment or framework red.

## Current line

- V0.2: green — frozen judge working as designed.
- V0.3: blue — terminal negative for mode sampling in the current regime.
- V0.4 Class 1: amber — valid open negative for sparse local edge intervention.

The machine-readable projection is `HAOS_GEN_STATUS.json`. Scientific status
strings remain unchanged; this layer changes interpretation and presentation,
not evidence.
