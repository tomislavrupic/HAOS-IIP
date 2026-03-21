# Root Execution Guide

- Run a checker with `python3 run_phase.py <phase> --check` from repo root.
- A phase authority bundle means: builder output passes its checker and has `manifest + summary + runs.json`.
- Read manifests first, not phase source, unless you are fixing that phase locally.
- Regime-defining phases are `6`, `8`, `10`, `11`, `14`, `16`, `17`, and `18`.
- The current top-level anchor is `PROGRAM_STATE_MILESTONE_18.md`.
- The callable surface is frozen in `API_CONTRACT.md`.
- Never retroactively edit: `phase6_operator_manifest.json`, `phase8_trace_manifest.json`, `phase10_manifest.json`, `phase11_manifest.json`, `phase14_manifest.json`, `phase16_manifest.json`, `phase17_manifest.json`, `phase18_manifest.json`.
- Never change Phase V readout definitions, the Phase VI `delta_h` hierarchy, or the Phase VIII short-time window.
- New phases should consume frozen artifacts, not older phase scripts.
- Prefer `telemetry/frozen_metrics.py` over re-deriving metrics inside new builders.
- If a result conflicts with a frozen manifest, shrink scope and write a bounded failure note.
- Next open question: can the frozen mesoscopic-to-proto-geometric stack be compressed into one bounded emergent-kinematics statement without breaking contract?
