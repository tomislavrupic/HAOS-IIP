# External Dataset Bridge Probe Comparison

| dataset | bridge_status | k_star_time | first_visible_failure_time | early_detection | controls early | gene_set_status | AP margin |
| --- | --- | --- | --- | --- | --- | --- | --- |
| gds112 | PROBE_FAIL_NO_EARLY_DETECTION | 5.0 | 5.0 | False | 0 | OPEN_GENE_SET_NOT_SUPPLIED | None |
| gds20 | PROBE_PASS_WITH_CONTROLS | 5.0 | 15.0 | True | 0 | OPEN_GENE_SET_NOT_SUPPLIED | None |
| gds20_go_0006970 | PROBE_PASS_WITH_CONTROLS | 5.0 | 15.0 | True | 0 | GENE_SET_EVALUATED | 0.10119419728911083 |

GDS112 is retained as an honest negative result. GDS20 is the current external pass with deterministic controls.
The GO:0006970 row is a narrow sourced gene-set comparator, not a full osmotic-stress biology validation.
