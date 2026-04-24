# External Dataset Bridge Validation

## System
Source: GDS112.soft.gz
Features analyzed: 2000
Samples analyzed: 5
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
- time_shuffle_control: early_detection=False, k_star_time=15.0, first_visible_failure_time=15.0
- feature_shuffle_control: early_detection=False, k_star_time=5.0, first_visible_failure_time=5.0
- row_permutation_control: early_detection=False, k_star_time=5.0, first_visible_failure_time=5.0

## Known Gene-Set Comparator
- gene_set_status: OPEN_GENE_SET_NOT_SUPPLIED
- matched_gene_count: None
- haos_average_precision: None
- fold_change_average_precision: None
- haos_beats_fold_change: None

## Result
- bridge_status: PROBE_FAIL_NO_EARLY_DETECTION
- baseline_stable: True
- k_star_time: 5.0
- first_visible_failure_time: 5.0
- early_detection: False
- contrast_pass: False
- control_early_detection_count: 0

## Failure Analysis
- failure_mode: NO_PRE_VISIBLE_POST_BASELINE_SAMPLE
- interpretation: The visible-failure proxy fires at the first post-baseline sample, so this dataset has no measured pre-visible window for early detection.
- recoverability_expression_match_correlation: 0.991822696285575
- first_standard_large_change_time: 15.0

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
