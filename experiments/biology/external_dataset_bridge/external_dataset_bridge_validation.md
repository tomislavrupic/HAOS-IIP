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
feature_rankings.csv ranks genes by HAOS-weighted loss score and simple fold-change comparators.
trajectory_identity_retention is a standard time-series comparator used to catch identity-breaking controls.

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
- recoverability_expression_match_correlation: 0.9930193053326482
- first_standard_large_change_time: 5.0

## Post-Mortem
The GDS112 probe failed because recoverability collapse and the visible expression-pattern proxy both occur at 5 minutes. The bridge therefore did not provide an earlier warning than direct expression-pattern movement.

This is a useful negative result. The current proxy is sensitive to the strong heat-shock transition, but GDS112 does not provide a measured post-baseline, pre-visible window under the default proxy.

## Repair Attempt
Two implementation corrections were added after the first failure:

- GEO log-ratio time courses that start after time 0 can receive an inferred zero baseline.
- `trajectory_identity_retention` was added as a standard time-series comparator so feature-shuffled gene identities count as visibly broken controls.

Exploratory GDS20 hyper-osmotic shock result:

- bridge_status: PROBE_PASS_WITH_CONTROLS
- baseline_stable: True
- k_star_time: 5.0
- first_visible_failure_time: 15.0
- early_detection: True
- contrast_pass: True
- control_early_detection_count: 0
- gene_set_status: OPEN_GENE_SET_NOT_SUPPLIED

This repaired pass does not erase the GDS112 failure. It shows that the bridge can pass on a different external stress-response time course when a measured pre-visible window exists and controls do not match.

## GO:0006970 Comparator
A sourced direct SGD/GO gene set was added for `GO:0006970` response to osmotic stress.

GDS20 with the GO comparator:

- gene_set_status: GENE_SET_EVALUATED
- known_gene_set_size: 24
- matched_gene_count: 10
- haos_average_precision: 0.14409276508511915
- fold_change_average_precision: 0.04289856779600831
- average_precision_margin: 0.10119419728911083
- haos_beats_fold_change: True

This is a narrow GO-term comparator pass. It is not a full ESR, HOG pathway, or systems-biology validation claim.

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
