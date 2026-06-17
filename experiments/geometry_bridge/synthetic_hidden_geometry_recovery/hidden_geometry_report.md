# Synthetic Hidden Geometry Recovery Report

- bridge id: GEO-HIDDEN-01
- version: synthetic-hidden-geometry-recovery-v0.1
- verdict: BENCHMARK_OPEN, MIXED_OPEN
- terminal labels: BENCHMARK_OPEN, TRANSFORMATION_RECOVERY_BOUNDARY_OPEN

This benchmark asks whether a frozen observer can recover distance, orientation, transformation class, and held-out relations from a hidden synthetic geometry.
- distance spearman: 0.972892
- orientation accuracy: 1.000000
- transform accuracy: 0.500000
- fiedler transform accuracy: 0.250000
- fiedler sign stability: 0.513889
- relation accuracy: 1.000000
- diagnosis: distance, orientation, and relations recover; normalized low-mode diagnostics remain below the frozen transformation threshold on holdout.
- fiedler note: low-mode diagnostics are recorded with the normalized Laplacian path and do not override the open verdict.
