# PIX-7 result infographics

These plates visualize the frozen `EL-R3-MESOSCOPIC-FIELD-01` contract and its
terminal `PARTIAL_RECOVERY_ONLY` result. They are explanatory graphics, not
additional evidence and not part of the frozen result hash.

## 01 — The 128-bit Contract

`01-the-128-bit-contract-pix7.png` shows the fixed `4 x 4` physical partition,
sixteen 8-bit quantized means, deterministic cell-local repair, rank `16`, and
nullity `n^2 - 16`. The prohibition seals mark the three central exclusions:
direct target phasors, full-state reconstruction, and adaptive range.

## 02 — The Blind Channel

`02-the-blind-channel-pix7.png` shows why the terminal recovery rate is `2/3`:
coarse-cell offsets and smooth physical twists are visible to `C_16` and were
recovered, while the within-cell balanced gradient lies in its nullspace and
was not recovered. The refinement strip preserves the same memory while the
microscopic nullity grows from `48` to `128` to `240`.

Both images were generated in built-in image-generation mode on 2026-08-18,
using the prior PIX-7 cathedral infographic as a style reference. Exact numeric
claims come from `../final/aggregate_result.json`.
