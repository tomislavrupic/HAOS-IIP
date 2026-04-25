# HAOS-IIP Biology Line

This directory contains the HAOS-IIP Biology Line application artifacts.

The biology line is a reproducibility and recoverability testbed. It does not modify frozen HAOS-IIP core, theory, telemetry, phase artifacts, canonical notes, or authority bundles.

It does not claim that HAOS-IIP explains biology. It tests whether bounded HAOS-style recoverability diagnostics can be made deterministic, auditable, and falsifiable on toy biological interaction systems and public external datasets.

## Line Index

| Line | Artifact | Question | Status |
| --- | --- | --- | --- |
| A | `gene_network_demo/` | Does a synthetic gene regulatory network lose recoverable coherence before visible regulatory failure? | PASS |
| B | `microtubule_lattice_demo/` | Does a microtubule-inspired cylindrical lattice lose recoverable interaction structure before visible lattice failure? | PASS |
| C | `tissue_3d_demo/` | Does a 3D multicellular interaction structure lose recoverable coherence before visible tissue failure? | PASS |
| D | `morphogenesis_4d_demo/` | Does a 3D tissue-like structure evolving through developmental time lose recoverable coherence before visible developmental failure? | PASS |
| E | `external_dataset_bridge/` | Can public stress-response expression datasets be audited with bounded recoverability diagnostics, controls, and standard comparators? | MIXED |

## One-Command External Audit

The lab-facing external bridge can audit current Line E outputs without credentials, network access, numpy, or matplotlib:

```bash
python3 experiments/biology/external_dataset_bridge/quick_reproduce.py --check
```

Current expected summary:

```text
quick_reproduce_check: PASS
runs_checked: 3
checks_passed: 41/41
outputs_written: experiments/biology/external_dataset_bridge/outputs/
```

This audit passing does not mean every dataset is a positive probe. It means the bounded expected outcomes reproduce:

- `GDS112` heat shock remains an honest negative: `PROBE_FAIL_NO_EARLY_DETECTION`.
- `GDS20` hyper-osmotic shock remains a controlled external pass.
- `GDS20` plus direct `GO:0006970` remains a narrow sourced gene-set comparator pass.

## Toy Line Runners

Run toy lines from the repository root:

```bash
python experiments/biology/gene_network_demo/run_gene_network_demo.py
python experiments/biology/microtubule_lattice_demo/run_microtubule_lattice_demo.py
python experiments/biology/tissue_3d_demo/run_tissue_3d_demo.py
python experiments/biology/morphogenesis_4d_demo/run_morphogenesis_4d_demo.py
```

If local Python lacks `numpy` or `matplotlib`, use:

```bash
uv run --with numpy --with matplotlib python <runner>
```

## Numbered Paper

The current numbered Biology Line paper is:

`papers/pdf_releases/54.2 HAOS-IIP Biology Line A-E Toy Recoverability Ladder and External Dataset Bridge.pdf`

It supersedes the earlier A-D-only biology paper for the full biology-line overview. The A-D paper remains available as paper `54.1`.

## Limits

- Toy lines A-D are proof-of-instrument demos, not biological validation.
- Line E uses expression-pattern proxies, controls, and comparators; it is not a phenotype or wet-lab validation.
- Visible failure in Line E is a proxy, not cell death, organismal failure, or biological dysfunction.
- No consciousness, Orch-OR, quantum biology, or QATC claims are made here.
- If a bridge dataset fails, the failure should remain visible and documented.
