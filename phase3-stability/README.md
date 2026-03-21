# Phase 3 Stability

This folder is the narrow authoritative Phase III closure path only.

- Scope: frozen consolidation around `23.10` to `23.14`
- Runner: `runs/run_phase3.py`
- Diagnostics module: `diagnostics/check_phase3_bundle.py`
- Imports: `haos_core` only, no local phase imports
- Outputs: JSON bundles under `runs/`

This phase does not regenerate historical plots or notes. It freezes the authoritative Phase III input chain into a small JSON-first execution path.
