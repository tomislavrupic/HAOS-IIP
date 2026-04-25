# Biology Line E: External Dataset Bridge

This is Biology Line E for HAOS-IIP.

It is an external validation bridge, not another toy biological model.

The bridge asks whether a real public biological perturbation dataset can be tested with the same bounded HAOS-style diagnostic form used in Biology Lines A-D:

> Does recoverability loss appear before a stricter visible or conventional failure proxy?

## Boundary

This layer does not modify frozen HAOS-IIP core, theory, telemetry, phase artifacts, canonical notes, or Biology Lines A-D.

It does not claim:

- HAOS-IIP explains biology
- HAOS-IIP proves new physics
- gene expression stress response is a literal collapse process
- a proxy visible-failure threshold is the same as biological death or dysfunction

## First External Target

Primary target for v0.1:

- NCBI GEO DataSet `GDS112`
- Reference series: `GSE18`
- Organism: `Saccharomyces cerevisiae`
- Description: heat shock from 30 C to 37 C time course
- Samples: 0, 5, 15, 30, and 60 minutes
- Value type: log ratio
- GEO page: <https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS112>
- GEO SOFT download: <https://ftp.ncbi.nlm.nih.gov/geo/datasets/GDSnnn/GDS112/soft/GDS112.soft.gz>

The raw external data are not committed by default. Place downloaded data under `data/raw/` or pass a local path at runtime.

## Quick Lab Reproduction

Use the lab-facing reproducer when the goal is fast auditability rather than interactive exploration.

Audit the committed/generated outputs without downloading data:

```bash
python3 experiments/biology/external_dataset_bridge/quick_reproduce.py --check
```

Expected audit summary:

```text
quick_reproduce_check: PASS
runs_checked: 3
checks_passed: 41/41
outputs_written: experiments/biology/external_dataset_bridge/outputs/
```

Run the full public-data reproduction:

```bash
python3 experiments/biology/external_dataset_bridge/quick_reproduce.py
```

The full run downloads only public no-credential GEO SOFT files into ignored `data/raw/`, records input digests, runs the bounded GDS112/GDS20 probes, refreshes `outputs/probe_comparison.*`, and writes `outputs/quick_reproduce_check.*`.

Useful bounded options:

- `--check`: audit existing outputs only; no network and no numpy/matplotlib required.
- `--only gds112`: audit or rerun one bounded target.
- `--no-download`: require raw files to already exist under `data/raw/`.

Digitized input definitions live in `input_manifest.json`.

## Run Without Data

```bash
python3 experiments/biology/external_dataset_bridge/run_external_dataset_bridge.py
```

This writes an honest `OPEN_NO_DATA` validation report.

## Run With GEO SOFT Data

```bash
uv run --with numpy --with matplotlib python experiments/biology/external_dataset_bridge/run_external_dataset_bridge.py \
  --soft experiments/biology/external_dataset_bridge/data/raw/GDS112.soft.gz
```

With an externally curated known response gene set:

```bash
uv run --with numpy --with matplotlib python experiments/biology/external_dataset_bridge/run_external_dataset_bridge.py \
  --soft experiments/biology/external_dataset_bridge/data/raw/GDS112.soft.gz \
  --known-gene-set experiments/biology/external_dataset_bridge/data/raw/known_esr_genes.txt
```

## Current v0.1 External Probe Result

The first GDS112 bridge run is intentionally recorded as a negative/open result:

```text
bridge_status: PROBE_FAIL_NO_EARLY_DETECTION
baseline_stable: True
k_star_time: 5.0
first_visible_failure_time: 5.0
early_detection: False
contrast_pass: False
control_early_detection_count: 0
gene_set_status: OPEN_GENE_SET_NOT_SUPPLIED
```

Interpretation: the current proxy does not detect recoverability loss before the visible expression-pattern proxy on GDS112. This does not validate Line E; it exposes the next refinement target.

Post-mortem:

- recoverability collapse and visible expression-pattern proxy both occur at 5 minutes
- the bridge is sensitive to the strong heat-shock transition
- it does not beat direct expression-pattern movement for early detection on this dataset
- next iteration should use curated ESR gene-set ranking, slower stress dynamics, or matched impaired-response datasets

## Repair Attempt Result

Two fixes were added after the first GDS112 failure:

- GEO log-ratio datasets that start after time 0 now receive an inferred zero baseline when appropriate.
- `trajectory_identity_retention` is added as a standard time-series comparator, not a HAOS metric. This makes feature-shuffled gene identity visibly fail instead of allowing an aggregate expression-magnitude false positive.

Exploratory GDS20 hyper-osmotic shock result:

```text
bridge_status: PROBE_PASS_WITH_CONTROLS
baseline_stable: True
k_star_time: 5.0
first_visible_failure_time: 15.0
early_detection: True
contrast_pass: True
control_early_detection_count: 0
gene_set_status: OPEN_GENE_SET_NOT_SUPPLIED
failure_mode: OBSERVED_EARLY_DETECTION_WINDOW
```

Interpretation: the repaired bridge can produce an external recoverability-loss probe pass on a different public Gasch stress-response dataset, but the biology-specific gene-set claim remains open until a curated external response-gene list is supplied.

## GO:0006970 Gene-Set Comparator

A sourced direct SGD/GO gene set was added:

- file: `data/reference/GO_0006970_response_to_osmotic_stress_sgd_genes.txt`
- source: <https://current.geneontology.org/annotations/sgd.gaf.gz>
- term: `GO:0006970` response to osmotic stress
- filter: direct `GO:0006970` annotations, excluding `NOT`
- gene count: 24

GDS20 with the sourced GO comparator:

```text
gene_set_status: GENE_SET_EVALUATED
known_gene_set_size: 24
matched_gene_count: 10
haos_average_precision: 0.14409276508511915
fold_change_average_precision: 0.04289856779600831
average_precision_margin: 0.10119419728911083
haos_beats_fold_change: True
```

Interpretation: on the matched direct GO:0006970 subset, the HAOS-weighted feature ranking beats simple maximum fold-change average precision. This is a narrow sourced GO-term comparator result, not a full ESR or HOG-pathway validation.

## Run With A CSV

```bash
uv run --with numpy --with matplotlib python experiments/biology/external_dataset_bridge/run_external_dataset_bridge.py \
  --csv path/to/expression_time_series.csv
```

CSV expectation:

- one feature or gene per row
- one or more leading identifier columns
- numeric sample/time columns

The bridge detects numeric columns and treats the first numeric column as baseline.

## Controls

The supplied-data path runs three deterministic negative controls:

- `time_shuffle_control`: keeps baseline first and scrambles later sample order
- `feature_shuffle_control`: preserves each sample distribution while breaking gene identity across samples
- `row_permutation_control`: preserves each feature's value range while breaking its trajectory

A probe pass requires:

- observed data show early detection
- all deterministic controls do not show the same early-detection signature

If controls also pass, Line E is `PROBE_FAIL_CONTROL_MATCH`, not a success.

Current failure labels:

- `PROBE_FAIL_NO_EARLY_DETECTION`: observed data do not show k_star before visible proxy
- `PROBE_FAIL_CONTROL_MATCH`: observed data show early detection, but at least one deterministic control also matches
- `PROBE_PASS_WITH_CONTROLS`: observed data show early detection and deterministic controls do not match

## Standard Comparator

The bridge also writes ordinary expression movement summaries:

- mean absolute change from baseline
- median absolute change from baseline
- fraction of features with absolute change >= 0.5
- fraction of features with absolute change >= 1.0

These are not HAOS metrics. They are included so the HAOS-style signal can be compared against plain expression-change behavior.

## Bounded Bio-Metrics

The lab bridge keeps the biological claims bounded to auditable quantities:

- `baseline_stable`: observed baseline recoverability is high.
- `early_detection`: observed `k_star_index` precedes the visible-failure proxy index for pass runs.
- `control_contrast`: deterministic negative controls do not reproduce the observed early-detection signature.
- `honest_negative`: a primary dataset can remain useful as a reproduced negative result.
- `gene_set_comparator`: with a curated gene set, HAOS-weighted feature ranking must beat simple fold-change by the declared margin.

## Known Gene-Set Comparator

If `--known-gene-set` is supplied, the bridge writes `feature_rankings.csv` and compares:

- HAOS-weighted feature loss ranking
- simple maximum absolute fold-change ranking

The comparison uses average precision and top-k enrichment against the supplied known gene set.

This is intentionally strict:

- if no gene set is supplied, the biology-specific claim is `OPEN_GENE_SET_NOT_SUPPLIED`
- if supplied genes do not match feature identifiers, it is `OPEN_NO_MATCHED_GENES`
- if HAOS ranking does not beat simple fold-change by the declared margin, the biology-specific ranking claim fails

Do not treat a recoverability/control pass as proof of yeast biology specificity.

## Outputs

- `outputs/probe_comparison.csv`
- `outputs/probe_comparison.md`
- `outputs/gds112/` primary failed probe artifacts
- `outputs/gds20/` exploratory repaired-pass probe artifacts
- `outputs/gds20_go_0006970/` sourced GO gene-set comparator artifacts
- `external_dataset_bridge_validation.md`

Each per-dataset output folder contains `bridge_status.json`, `results.csv`, `standard_expression_metrics.csv`, `feature_rankings.csv`, control summaries, failure analysis, and plots.
