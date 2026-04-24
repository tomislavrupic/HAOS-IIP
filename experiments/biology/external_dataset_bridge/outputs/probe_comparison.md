# External Dataset Bridge Probe Comparison

## Summary

Biology Line E now records both a failed primary probe and a repaired exploratory pass.

| dataset | source | bridge_status | k_star_time | first_visible_failure_time | early_detection | contrast_pass | gene_set_status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GDS112 | Heat shock from 30 C to 37 C | PROBE_FAIL_NO_EARLY_DETECTION | 5.0 | 5.0 | False | False | OPEN_GENE_SET_NOT_SUPPLIED |
| GDS20 | Hyper-osmotic shock time course | PROBE_PASS_WITH_CONTROLS | 5.0 | 15.0 | True | True | OPEN_GENE_SET_NOT_SUPPLIED |
| GDS20_GO_0006970 | Hyper-osmotic shock with sourced GO:0006970 comparator | PROBE_PASS_WITH_CONTROLS | 5.0 | 15.0 | True | True | GENE_SET_EVALUATED |

## Reading

GDS112 remains a real external failure for early detection because the visible expression-pattern proxy fires at the first post-baseline sample.

GDS20 becomes an exploratory external probe pass after adding an inferred zero baseline for log-ratio time courses and a trajectory-identity comparator that prevents feature-shuffled controls from passing as valid early-detection signals.

The sourced direct `GO:0006970` comparator evaluates 10 matched genes from a 24-gene SGD/GO response-to-osmotic-stress set:

- haos_average_precision: `0.14409276508511915`
- fold_change_average_precision: `0.04289856779600831`
- average_precision_margin: `0.10119419728911083`
- haos_beats_fold_change: `True`

This is a narrow GO-term comparator pass, not a full osmotic-stress regulon validation.
