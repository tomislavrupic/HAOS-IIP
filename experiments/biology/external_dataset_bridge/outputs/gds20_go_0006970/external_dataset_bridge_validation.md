# External Dataset Bridge Validation

## System
Source: GDS20.soft.gz
Features analyzed: 2000
Samples analyzed: 7
Primary v0.1 target: public yeast heat-shock expression time course if GDS112 is supplied.

## Metrics
recoverability is a weighted proxy combining expression_match, step_coherence, and support_retention.
delta_persistence is the change in recoverability between ordered samples.
k_star is the first sustained recoverability crossing below 0.70.
safety_margin is distance in time or ordered-sample units to k_star.
visible_failure is a strict expression-pattern or trajectory-identity proxy, not a biological failure label.
standard expression metrics report ordinary absolute expression movement at 0.5 and 1.0 log-ratio thresholds.
trajectory_identity_retention is a standard time-series comparator used to catch identity-breaking controls.
feature_rankings.csv ranks genes by HAOS-weighted loss score and simple fold-change comparators.

## Negative Controls
The bridge runs deterministic time-shuffle, feature-shuffle, and row-permutation controls.
A probe pass requires observed early detection and zero matching early-detection hits in these controls.
- time_shuffle_control: early_detection=False, k_star_time=15.0, first_visible_failure_time=60.0
- feature_shuffle_control: early_detection=False, k_star_time=5.0, first_visible_failure_time=5.0
- row_permutation_control: early_detection=False, k_star_time=5.0, first_visible_failure_time=5.0

## Known Gene-Set Comparator
- gene_set_status: GENE_SET_EVALUATED
- matched_gene_count: 10
- haos_average_precision: 0.14409276508511915
- fold_change_average_precision: 0.04289856779600831
- haos_beats_fold_change: True

## Result
- bridge_status: PROBE_PASS_WITH_CONTROLS
- baseline_stable: True
- k_star_time: 5.0
- first_visible_failure_time: 15.0
- early_detection: True
- contrast_pass: True
- control_early_detection_count: 0

## Failure Analysis
- failure_mode: OBSERVED_EARLY_DETECTION_WINDOW
- interpretation: Observed data have k_star before the visible-failure proxy; controls decide whether this is a probe pass.
- recoverability_expression_match_correlation: 0.9966996164074223
- first_standard_large_change_time: None

## Explicit Failure Conditions
- The bridge fails if observed data do not show early detection.
- The bridge fails if deterministic controls show the same early-detection signature.
- A biology-specific claim fails if HAOS feature ranking does not beat simple fold-change ranking against a supplied curated gene set.
- The bridge remains open if visible_failure is only a proxy and no external phenotype label is supplied.

## Limitations
- External expression data do not by themselves provide visible biological failure labels.
- This bridge tests diagnostic ordering, not biological mechanism.
- This is not proof that HAOS-IIP explains biology.
- This does not modify or prove HAOS-IIP core.
