# External Validation

This folder holds the repo-local copy of the external HAOS-IIP validation layer.

Boundary rule:

- `HAOS-IIP/` core and frozen telemetry stay authoritative
- this folder is a runner layer above that authority surface
- nothing here modifies, redefines, or extends HAOS-IIP core

Current contents:

- `data/` graph artifacts consumed by the runner
- `runners/` the external filtration runner, HAOS-IIP sensor adapter, and prompt/contract notes
- `results/` deterministic CSV and summary outputs generated from the current artifact set

Toy backend run:

```bash
python3 external_validation/runners/run_validation.py
```

HAOS-IIP-backed run:

```bash
python3 external_validation/runners/run_validation.py --sensor haos_iip
```

Note:

- the `haos_iip` backend requires Python with `numpy` available
- in the Codex desktop environment, the bundled runtime was used for that path
