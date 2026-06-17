# Post-67.1 Consolidated Status Snapshot

Status: continuity snapshot, not a new phase and not a new claim ceiling.

This page consolidates the current reading across the 66.5 scale-bridge
baseline, the HBP registry, the geometry bridge, and the equivalences layer.
It is a coordination note, not a mechanism statement.

## 1. Hardened 66.5 Baseline

- 66.5 remains the continuity baseline for the scale-bridge / foundations
  track.
- The current public smoke test remains:
  `uv run python examples/quick_reproduce.py`
- The scale-bridge hardening pass remains bounded: it clarifies where the
  ladder stops, but it does not reopen the claim boundary.

## 2. HBP Registry

The HBP registry remains live and bounded. The current terminal readings are:

| Bridge | Status | Terminal reading |
| --- | --- | --- |
| PB-01 | `PREDICTION_NOT_DISTINCT_FROM_BASELINES` | control-invalid synthetic calibration |
| PB-02 | `PRECOMMITMENT_DRAFT` | frozen draft for external PowerGraph holdout |
| PB-03 | `PREDICTION_NOT_DISTINCT_FROM_BASELINES` | temporal recovery boundary open |
| PB-04 | `PREDICTION_NOT_DISTINCT_FROM_BASELINES` | cross-domain transfer boundary open |

Reference labels:

- PB-01: `CONTROL_INVALID`, `HOLDOUT_TRANSFER_PASS`, `MIXED_OPEN`,
  `PHYSICAL_MECHANISM_NOT_ESTABLISHED`,
  `PREDICTION_NOT_DISTINCT_FROM_BASELINES`
- PB-02: `OPERATIONAL_MAPPING_PARTIAL`, `OPERATIONAL_MAPPING_VALID`,
  `PHYSICAL_MECHANISM_NOT_ESTABLISHED`
- PB-03: `DIMENSIONAL_ANALYSIS_FAIL`, `FORMAL_CORRESPONDENCE_ONLY`,
  `OPERATIONAL_MAPPING_PARTIAL`, `PHYSICAL_MECHANISM_NOT_ESTABLISHED`
- PB-04: `CROSS_DOMAIN_TRANSFER_BOUNDARY_OPEN`, `DIMENSIONAL_ANALYSIS_FAIL`,
  `FORMAL_CORRESPONDENCE_ONLY`, `OPERATIONAL_MAPPING_PARTIAL`,
  `PHYSICAL_MECHANISM_NOT_ESTABLISHED`

Implemented fact:
- PB-01 is still the synthetic calibration benchmark.
- PB-02 remains a frozen precommitment draft.
- PB-03 and PB-04 remain boundary cases, not upgraded bridges.

## 3. Geometry Bridge

The geometry bridge remains a synthetic calibration scaffold. Current chain
reading:

- `GEO-01`: open / non-dominant on holdout
- `GEO-02`: open
- `GEO-03`: open
- `GEO-04`: valid
- `GEO-HIDDEN-01`: open
- `GEO-T1-01`: valid
- `GEO-P1-01`: valid

The current sharpened GEO-HIDDEN-01 reading is:

- distance and orientation recover well;
- transformation recovery remains the weak point;
- normalized low-mode diagnostics improved;
- Cheeger conductance records a bottleneck signal;
- `TRANSFORMATION_RECOVERY_BOUNDARY_OPEN` remains the correct terminal
  reading.

Frozen spectral / curvature numbers:

- `transform_accuracy`: `0.500000`
- `fiedler_transform_accuracy`: `0.250000`
- `fiedler_sym_accuracy`: `0.250000`
- `fiedler_rw_accuracy`: `0.312500`
- `fiedler_sign_stability`: `0.513889`
- `cheeger_conductance`: `0.304962`

Related notes:

- [Geometry Bridge Chain](../../geometry_bridge/README.md)
- [Spectral Diagnostics Summary](../../geometry_bridge/spectral_diagnostics_summary.md)

## 4. Equivalences Layer

The equivalences layer is now mature enough to read as one compact structural
toolkit:

1. Spectral Address and Laplacian Eigenmodes
2. Current Closure and Discrete Conservation Balance
3. Hodge Decomposition and Discrete Sector Splitting
4. Hodge Laplacian Extensions and Higher-Order Sector Recovery
5. Discrete Curvature Models and Local-to-Global Diagnostics
6. Discrete Ricci Flow on Graphs as a Diagnostic Probe

The current reading is that these notes form a layered operator story:

- spectral address names recoverable identity in the harmonic / eigenmode
  neighborhood;
- current closure names flux-like balance in the exact / coexact neighborhood;
- Hodge decomposition names the decomposition rule for those neighborhoods on
  discrete cochain complexes;
- curvature and Ricci-flow-like updates now sit beside that picture as local
  diagnostics for why the global structure may be hard to recover.

Related note index:

- [Equivalences Layer Summary](../equivalences/equivalences_layer_summary.md)
- [Equivalences Index](../equivalences/equivalences_index.md)

## 5. Current Reading

Implemented:

- the 66.5 baseline remains intact;
- the HBP registry remains bounded and live;
- the geometry bridge now has a sharper synthetic bottleneck reading;
- the equivalences layer now includes spectral, Hodge, curvature, and
  Ricci-flow-like diagnostics.

What is still open:

- PB-02 is still a draft, not a scored bridge;
- GEO-HIDDEN-01 still does not recover transformation classes robustly;
- no note here upgrades any claim ceiling.

## Boundary

This snapshot is a coordination artifact. It preserves the separation between
frozen diagnostics and unsupported physics claims.
