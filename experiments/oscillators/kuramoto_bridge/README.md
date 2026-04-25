# Kuramoto Oscillator Bridge

This is a sidecar oscillator bridge for HAOS-IIP. It does not modify frozen
core APIs, phase artifacts, theory notes, or Biology Line files.

The bridge asks a bounded question:

> Do HAOS-style recoverability, persistence, and identity-retention proxies
> detect oscillator coherence before or more robustly than the classic
> Kuramoto order parameter `R`?

## Run

The local Python environment may not include numerical packages. Use `uv` for a
reproducible run:

```bash
uv run --with numpy --with matplotlib python experiments/oscillators/kuramoto_bridge/run_kuramoto_bridge.py
```

To use a specific graph:

```bash
uv run --with numpy --with matplotlib python experiments/oscillators/kuramoto_bridge/run_kuramoto_bridge.py \
  --graph path/to/graph.json
```

Useful options:

- `--k-min`, `--k-max`, `--k-count`: critical-coupling scan range.
- `--steps`, `--dt`: integration horizon.
- `--frequency-mode gaussian|custom`: natural-frequency source.
- `--frequency-file`: newline or CSV numeric frequencies for custom mode.
- `--higher-order`: enable deterministic triangle coupling.
- `--max-nodes`: cap large HAOS artifacts for fast bridge runs.
- `--output-dir`: write outputs somewhere else.

## Outputs

Generated files go to ignored `outputs/`:

- `bridge_status.json`
- `probe_comparison.md`
- `probe_comparison.csv`
- `synchronization_curves.png`
- `proxy_evolution.png`
- `control_comparison.png`
- `failure_analysis.md`

## Data Loading

The loader first tries an explicit `--graph`, then common HAOS paths under
`data/`, `telemetry/`, `geometry_emergence/`, and `experiments/`.

It accepts JSON adjacency matrices, edge lists, Laplacians, and a few common
wrapped forms. If no serialized graph is found, it reconstructs the configured
`geometry_emergence` graph from its committed config and marks the source as a
HAOS-derived reconstruction. If even that is unavailable, it falls back to a
deterministic ring and reports `OPEN_NO_DATA_SYNTHETIC`.

## Status Semantics

- `PASS`: observed HAOS proxies detect coherence earlier than standard `R`, and
  deterministic controls do not reproduce the same early-detection signature.
- `MARGINAL`: observed signal is coherent but early detection or control
  contrast is incomplete.
- `FAIL`: proxies underperform or controls match the observed result.
- `OPEN_NO_DATA_SYNTHETIC`: no HAOS graph source was available; synthetic output
  is a plumbing check only.

These are bridge diagnostics, not physics claims.
