# Public Reproduction Spine

This reproduces a bounded mesoscopic-to-proto-geometric feasibility arc. It does not claim a continuum limit or physical ontology.

Run one command from repo root:

```bash
python3 examples/quick_reproduce.py
```

This path is intentionally narrow.

- It runs frozen Phase XV-XVIII bundle checks.
- It rebuilds the artifact-only continuum sketch.
- It writes a compact public table to `examples/output/quick_reproduce_output.csv`.
- It writes the corresponding figure to `examples/output/quick_reproduce_plot.svg`.
- It verifies both against the frozen public baselines in `examples/expected_output.csv` and `examples/expected_plot.svg`.

Expected terminal transcript:

```text
$ python3 examples/quick_reproduce.py
Reproduction spine passed.
Table: /Volumes/Samsung T5/2026/HAOS/HAOS DOCS/HAOS-IIP/examples/output/quick_reproduce_output.csv
Plot: /Volumes/Samsung T5/2026/HAOS/HAOS DOCS/HAOS-IIP/examples/output/quick_reproduce_plot.svg
Statement: Frozen propagation -> ordering -> causality -> distance-surrogate reproduces from stored artifacts, and the artifact-only continuum sketch stays bounded and branch-distinct. No new dynamics or continuum ontology are introduced.
```

Modest statement:

Frozen propagation -> ordering -> causality -> distance-surrogate reproduces from stored artifacts, and the artifact-only continuum sketch stays bounded and branch-distinct. No new dynamics or continuum ontology are introduced.

## Stability Demo

Run:

```bash
python3 -m haos_iip.demo stability
```

Machine-readable output:

```bash
python3 -m haos_iip.demo stability baseline --json
```

Legacy example wrapper:

```bash
python3 examples/stability_demo.py
```

This demo loads three synthetic scenarios from `examples/scenarios/`, evaluates them with the frozen telemetry layer, and writes a compact CSV plus an SVG metric panel to `examples/output/`.

Generator mode:

```bash
python3 -m haos_iip.demo stability baseline --noise 0.05 --connectivity-drop 0.2 --cluster-split
```

Batch scan mode:

```bash
python3 -m haos_iip.demo stability --scan noise=0.00:0.10:0.05 connectivity=0.00:0.20:0.10
```

Friendly metric aliases:

- `structural_retention` = persistence
- `temporal_consistency` = ordering
- `causal_deformation` = depth drift
- `geometric_integrity` = distance coherence

Scope:

It is a practical telemetry demonstration, not a new phase and not a physics claim. The distance score is a lightweight synthetic proxy so the frozen persistence, ordering, and causal diagnostics can be shown in a small public example.
