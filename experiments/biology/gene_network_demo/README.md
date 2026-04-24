# Gene Network Demo

This is the first HAOS-IIP biology-line experiment.

It is an experimental application layer only. It does not modify frozen HAOS-IIP core code, theory, telemetry, phase artifacts, authority bundles, canonical documents, paper files, or already frozen experiment outputs.

The demo tests whether lightweight HAOS-style stability diagnostics can detect early degradation in a synthetic gene regulatory network under controlled perturbation. It does not claim to explain biology, validate a biological theory, or prove HAOS-IIP core behavior.

Run from the repository root:

```bash
python experiments/biology/gene_network_demo/run_gene_network_demo.py
```

Primary outputs are written to `experiments/biology/gene_network_demo/outputs/`.
