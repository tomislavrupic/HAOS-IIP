# Phase 5 Readout

This folder is a fresh JSON-first scaffold over frozen Phase III and Phase IV outputs.

- Scope: readout only, no new theory and no plots
- Runner: `runs/run_phase5.py`
- Diagnostics module: `diagnostics/check_phase5_bundle.py`
- Modes: `dummy` for standalone validation, `frozen` for authoritative Phase III/IV inputs
- Imports: `haos_core` only, no phase imports
- Outputs: JSON bundles under `runs/`

This scaffold prefers the new Phase III and Phase IV bundles when they exist, then falls back to the frozen authoritative `23.14` and `24.5` JSONs.
