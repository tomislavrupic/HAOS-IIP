# Migration Map

## Active Runtime

| Old location | New location | Notes |
| --- | --- | --- |
| `numerics/simulations/stage23_14_phaseIII_final_consolidation_and_freeze.py` | `phase3-stability/runs/run_phase3.py` | New active Phase III closure runner. |
| `numerics/simulations/stage24_5_phaseIV_consolidation_and_freeze.py` | `phase4-sector-freeze/runs/run_phase4.py` | New active Phase IV closure runner. |
| none | `phase5-readout/runs/run_phase5.py` | Fresh Phase V scaffold with `dummy` and `frozen` modes. |
| none | `phase3-stability/diagnostics/check_phase3_bundle.py` | Independent bundle checker. |
| none | `phase4-sector-freeze/diagnostics/check_phase4_bundle.py` | Independent bundle checker. |
| none | `phase5-readout/diagnostics/check_phase5_bundle.py` | Independent bundle checker. |
| `archive-pre-refactor/numerics/simulations/L1_transverse_band_test.py` | `numerics/simulations/L1_transverse_band_test.py` | Active transverse-band utility reimplemented locally so stage-3/4 numerics no longer rely on archive-only runtime imports. |
| `archive-pre-refactor/numerics/simulations/L1_transverse_scaling_test.py` | `numerics/simulations/L1_transverse_scaling_test.py` | Active real-valued transverse scaling shim now wraps current stage-8 primitives so stage-2/3 numerics no longer depend on archive-only imports. |

## Core Extraction

| Primitive source quarry | New location | Notes |
| --- | --- | --- |
| `numerics/simulations/stage8_common.py` | `haos_core/core.py` | Periodic graph and operator primitives lifted without plotting/orchestration. |
| `numerics/simulations/stage9_common.py` | `haos_core/core.py` | Transport evolution primitives lifted into frozen APIs. |
| `numerics/simulations/stage9b_common.py` | `haos_core/core.py` | DK transport helpers lifted into reusable transport state helpers. |
| `numerics/simulations/stage10_common.py` | `haos_core/core.py` | Invariant and topology helpers lifted without reporting code. |
| `numerics/simulations/L1_stage2_common.py` | `haos_core/core.py` | Builder patterns kept only as reusable primitive precedents. |
| `numerics/simulations/DH_stage5_common.py` | `haos_core/core.py` | Shared invariant ideas retained only when structurally primitive. |
| `numerics/simulations/DK_stage6_common.py` | `haos_core/core.py` | DK primitive transport pieces retained only when phase-agnostic. |
| `numerics/simulations/DK_stage7_common.py` | `haos_core/core.py` | Selector and invariant helpers retained only when phase-agnostic. |
| `numerics/simulations/*_runs.json` for `23.10` to `24.5` | `phase3-stability/configs/`, `phase4-sector-freeze/configs/`, `phase5-readout/configs/` | Copied into phase-local configs so each phase can run independently. |

## Authority Preservation

- The authoritative note, JSON, CSV, and manifest paths remain unchanged.
- The old `stage23_10` to `stage24_5` simulation scripts remain readable as provenance anchors while the new phase folders become the active execution path.
- The new phase bundles carry compact authority references instead of embedding full upstream JSON payloads.
- Foundations cleanup release `50.2` closes the post-`50.1` note-level caveats `F2-F4` without changing numerical authority paths or promoting unfinished threshold work into frozen artifacts.
- Post-foundations validation release `50.3` freezes the transition from theorem-ladder cleanup to validation order, but still does not promote unfinished threshold or active-branch closure work into frozen numerical authority.
- Active-branch validation closure release `51.1` freezes the executed bounded `V1-V6b` tranche, including bounded universality and divergence-aware effective-law closure on the validated mild-disorder branch, while preserving the open clean-baseline identity boundary as an explicit next-program item rather than promoting it into a closed claim.
- Clean-baseline identity resolution release `51.2` freezes the bounded `CB1-lite`, `CB2`, `CB4`, and `CB3` closure path, promoting the clean vector identity issue from an open branch-survival question to a resolved symmetry-degenerate first-shell family statement without claiming a full exact-zero irrep package.

## Archive Policy

- `archive-pre-refactor/` is for exploratory and legacy runners only.
- `haos_core/` must not absorb plotting, report writing, CSV rendering, orchestration, or experiment narratives.
- No file under `phase3-stability/`, `phase4-sector-freeze/`, or `phase5-readout/` may import another phase script.
