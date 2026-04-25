# Biology Line E Dataset Targets

## Primary v0.1 Target: GDS112

Source:

- NCBI GEO DataSet: `GDS112`
- GEO series: `GSE18`
- GEO page: <https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS112>
- GEO SOFT file: <https://ftp.ncbi.nlm.nih.gov/geo/datasets/GDSnnn/GDS112/soft/GDS112.soft.gz>
- GEO series page: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18>
- SGD dataset page: <https://www.yeastgenome.org/dataset/GSE18>
- PubMed: <https://pubmed.ncbi.nlm.nih.gov/11102521/>
- Study: Gasch AP, Spellman PT, Kao CM, Carmel-Harel O, Eisen MB, Storz G, Botstein D, Brown PO. Genomic expression programs in the response of yeast cells to environmental changes. Mol Biol Cell. 2000.

Why this target:

- real public gene expression perturbation data
- simple time course
- stress-response context
- small enough to serve as a first bridge target
- public GEO metadata expose sample count, time-course framing, and log-ratio value type

Known scope:

- This is not a tissue or microscopy dataset.
- It is not a visible-failure dataset.
- Any `visible_failure` value in this bridge is a strict expression-pattern proxy, not a biological failure label.

Required controls:

- time-shuffled ordered samples
- feature-shuffled sample distributions
- row-permuted feature trajectories
- ordinary absolute expression-change summaries
- optional curated known response gene set for head-to-head ranking against simple fold-change

Failure conditions:

- observed data do not show early detection
- controls show the same early-detection signature
- the result cannot be separated from plain expression-change movement
- HAOS feature ranking fails to beat simple fold-change against a supplied curated gene set
- the report implies biological failure where only an expression proxy exists

## Lab Reproduction Contract

Biology Line E now ships a bounded one-command reproducer:

```bash
python3 experiments/biology/external_dataset_bridge/quick_reproduce.py --check
```

`--check` audits existing outputs without network access, credentials, numpy, or matplotlib. A full run without `--check` downloads public GEO SOFT files into ignored `data/raw/`, records input digests, reruns the bounded probes, and audits the refreshed outputs.

Digitized inputs and expected bounded outcomes are declared in `input_manifest.json`.

The audit is allowed to pass while GDS112 remains a failed probe because the expected result is an honest negative:

- GDS112 must reproduce `PROBE_FAIL_NO_EARLY_DETECTION`.
- GDS20 must reproduce `PROBE_PASS_WITH_CONTROLS`.
- GDS20 plus direct `GO:0006970` must reproduce the narrow gene-set comparator pass.

## Current v0.1 Probe Status

The first GDS112 run is a negative/open external result:

- bridge_status: `PROBE_FAIL_NO_EARLY_DETECTION`
- baseline_stable: `True`
- k_star_time: `5.0`
- first_visible_failure_time: `5.0`
- early_detection: `False`
- contrast_pass: `False`
- control_early_detection_count: `0`
- gene_set_status: `OPEN_GENE_SET_NOT_SUPPLIED`

Interpretation:

The current Line E proxy does not show recoverability-loss detection before the strict expression-pattern failure proxy on GDS112. This is not a validation pass.

## Repair Probe: GDS20

GDS20 is a public hyper-osmotic shock time course from the same Gasch stress-response program.

After adding an inferred zero baseline for log-ratio time courses that begin after time 0, and adding `trajectory_identity_retention` as a standard control-facing comparator, the GDS20 exploratory run produced:

- bridge_status: `PROBE_PASS_WITH_CONTROLS`
- baseline_stable: `True`
- k_star_time: `5.0`
- first_visible_failure_time: `15.0`
- early_detection: `True`
- contrast_pass: `True`
- control_early_detection_count: `0`
- gene_set_status: `OPEN_GENE_SET_NOT_SUPPLIED`

Interpretation:

This is an external probe pass for the recoverability-loss ordering test, not a biology-specific validation claim. The gene-set ranking claim remains open.

## GO:0006970 Comparator Result

The direct SGD/GO `GO:0006970` response-to-osmotic-stress gene set was extracted from the current GO Consortium SGD GAF:

- source: <https://current.geneontology.org/annotations/sgd.gaf.gz>
- direct `GO:0006970` annotations only
- `NOT` qualifiers excluded
- gene count: 24

GDS20 head-to-head result:

- gene_set_status: `GENE_SET_EVALUATED`
- matched_gene_count: `10`
- haos_average_precision: `0.14409276508511915`
- fold_change_average_precision: `0.04289856779600831`
- average_precision_margin: `0.10119419728911083`
- haos_beats_fold_change: `True`

Interpretation:

This supports a narrow biology-specific ranking claim for the matched direct GO:0006970 set. It does not validate the full environmental stress response, HOG pathway dynamics, or all osmotic-response biology.

## Candidate Later Targets

- Public organoid or tissue imaging time series with morphology labels.
- Public microtubule dynamics microscopy datasets with catastrophe, rescue, or lattice-damage annotations.
- Public cell stress-response RNA-seq time series with explicit viability or phenotype readouts.
- Public perturb-seq or CRISPR perturbation datasets where perturbation strength and downstream expression recovery can be separated.

## Acceptance Rule

A dataset can become a Line E validation target only if it has:

- public source and citation
- clear perturbation or stress axis
- repeated or ordered samples
- measurable feature state
- explicit proxy or ground-truth failure marker
- honest limitations when the marker is only a proxy
