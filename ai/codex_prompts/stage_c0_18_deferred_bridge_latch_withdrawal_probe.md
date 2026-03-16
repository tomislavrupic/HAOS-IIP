# Stage C0.18 - Deferred Bridge Latch Withdrawal Probe

## Purpose

Correct the main fairness issue in C0.17.

In C0.17 the bridge-local latch acted during the arming window itself, which could smear the very braid topology the stage was supposed to test for retention.

C0.18 therefore uses a deferred-freeze protocol:

1. arm under pure exact-match harmonic compatibility with no latch feedback
2. accumulate the loop/path ledger passively during arming
3. freeze the latch only at the end of arming
4. withdraw the selector
5. continue under the frozen latched operator

The question is narrower and cleaner:

`If the latch is not allowed to disturb the arming braid, can it retain topology after selector withdrawal?`

## Frozen architecture

- same signed C0-to-DK bridge scaffold as the bridge validation note
- same three representative families
- same harmonic-address dressing
- same bounded latch families as C0.17

## Latch families

- `cycle_occupancy_latch`
- `edge_path_reuse_latch`
- `combined_cycle_path_latch`

## Positive criterion

Positive if at least one latch family keeps the clustered bridge family in `braid_like_exchange` after withdrawal or produces a reproducible nonzero braid dwell band.

## Negative criterion

Negative if the clustered family arms into braid correctly but still collapses immediately after withdrawal under all three latch rules.

## Required outputs

- stamped JSON
- stamped CSV
- note target:
  - `Stage_C0_18_Deferred_Bridge_Latch_Withdrawal_Probe_v1.md`
