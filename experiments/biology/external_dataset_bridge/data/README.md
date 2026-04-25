# External Data Staging

Raw external data are not committed by default.

Use this directory for local staging:

```text
data/raw/
```

`quick_reproduce.py` can download the public GEO SOFT files listed in `../input_manifest.json` into `data/raw/`. These files are ignored by git.

Recommended first target:

```text
data/raw/GDS112.soft.gz
```

Optional known-gene-set target:

```text
data/raw/known_esr_genes.txt
```

Use one gene identifier per line. The bridge does not ship a curated ESR list by default, because that list should be sourced and cited externally rather than reconstructed from memory.

The bridge can also read CSV files from any local path.

Committed reference gene sets live under:

```text
data/reference/
```
