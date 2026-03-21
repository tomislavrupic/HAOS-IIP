# AGENTS.md

- This is an experimental physics style codebase.
- Simplicity beats feature growth.
- `haos_core` public APIs are frozen: `build_graph`, `run_transport`, `apply_selector`, `compute_invariants`.
- Phase folders are measurement rigs and may import only from `haos_core`.
- Do not add cross-phase imports, shared globals, or hidden utilities.
- Prefer small deterministic functions and JSON-first outputs.
- Do not re-import or revive code from `archive-pre-refactor/`.
- Preserve authoritative data, manifests, and notes; use `MIGRATION_MAP.md` when logic moves.
