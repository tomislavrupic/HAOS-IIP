# Phase 56.2 FMO Pathway-Strengthening Variant Sweep

This sweep tests bounded address-rule variants after the Phase 56.1 FMO falsifier.
It should be read as a diagnostic sweep, not as biological validation.

| variant | runs | pass_rate | recoverability | site_identity | pathway_identity | active_null_z | max_controls |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| sink_spectral | 50 | 0.080 | 0.827528 +/- 0.037506 | 0.988291 +/- 0.004913 | 0.622057 +/- 0.089951 | 2.659739 +/- 0.507998 | 2 |
| spectral_baseline | 50 | 0.040 | 0.803724 +/- 0.039498 | 0.984934 +/- 0.005867 | 0.546597 +/- 0.094994 | 2.472024 +/- 0.574870 | 2 |
| environment_assisted | 50 | 0.060 | 0.491183 +/- 0.222524 | 0.884084 +/- 0.085967 | 0.128198 +/- 0.226975 | 0.901116 +/- 1.224950 | 2 |
| hybrid_sink | 50 | 0.040 | 0.378582 +/- 0.252680 | 0.830612 +/- 0.111975 | 0.091584 +/- 0.205845 | 0.262384 +/- 1.427271 | 2 |
| hybrid | 50 | 0.040 | 0.238918 +/- 0.263122 | 0.740080 +/- 0.144044 | 0.063999 +/- 0.182341 | -0.307723 +/- 1.459756 | 1 |
| local_only | 50 | 0.000 | 0.053467 +/- 0.126066 | 0.406772 +/- 0.208194 | 0.000000 +/- 0.000000 | -1.125053 +/- 0.859859 | 0 |
