# Phase 4 Sector Freeze

This folder is the narrow authoritative Phase IV closure path only.

- Scope: frozen consolidation around `24.1` to `24.5`
- Runner: `runs/run_phase4.py`
- Diagnostics module: `diagnostics/check_phase4_bundle.py`
- Imports: `haos_core` only, no local phase imports
- Outputs: JSON bundles under `runs/`

This phase consumes the frozen authoritative Phase III and Phase IV artifacts and re-freezes them into a small JSON-first execution path.
