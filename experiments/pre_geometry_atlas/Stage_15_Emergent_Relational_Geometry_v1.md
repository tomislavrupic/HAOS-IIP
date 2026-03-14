# Stage 15 Emergent Relational Geometry v1

This file freezes the minimal pilot grammar for Stage 15 before any implementation work begins.

## Scope

Stage 15 is a diagnostic probe, not a spacetime construction stage.

Question:
- can frozen collective interaction data induce a stable relational ordering that behaves metric-like under refinement?

Interpretation boundary:
- geometry is treated only as a compression of relational ordering
- no spacetime, gravitational, or physical metric claims are allowed here

## Pilot Size

The first pilot is intentionally small:
- 3 cases
- 2 resolutions (`12 -> 24`)
- 3 induced distance candidates

Total analysis cells:
- `3 x 2 x 3 = 18`

## Chosen Cases

1. `asymmetric_density_cluster`
- role: clustered composite anchor
- reason: compact triplet composite with non-degenerate triangle and ball-growth diagnostics

2. `counter_propagating_corridor`
- role: transport control
- reason: strongest transport-dominated control for directionality and shell-order tests

3. `symmetric_wide_triad`
- role: diffuse triad control
- reason: weakly structured field for testing whether metric-like ordering fails cleanly

Note:
- the two-packet `tight_clustered_pair` is intentionally left out of the first pilot because triangle and loop tests would be degenerate there
- it can be added later as a sanity cross-check if the first pilot is readable

## Distance Candidates

### 1. Persistence Distance

Definition:
- `d_p(i,j) = 1 / (epsilon + <overlap_ij>_T)`
- uses the full sampled time window
- symmetric by construction

Purpose:
- tests whether relational closeness can be induced from overlap persistence

### 2. Phase-Correlation Distance

Definition:
- `d_phi(i,j) = 1 - |<phase_lock_ij>_T|`

Purpose:
- tests whether synchronization induces a stable relational ordering

### 3. Transport-Delay Distance

Definition:
- `d_tau(i,j) = t_arrival(i->j)`
- arrival is detected when the smoothed envelope response at packet `j` crosses a fixed threshold after excitation at `i`

Purpose:
- tests whether influence ordering produces cone-like relational structure

## Metric-Like Tests

For every case and every distance candidate, compute:
- mean and max symmetry deviation
- triangle-inequality violation rate where packet count allows it
- pair-ordering persistence across `12 -> 24`

Threshold intent:
- high violation or high ordering drift -> no coherent ordering
- moderate but bounded violations -> weak ordering
- low violations with refinement drift -> scale-fragile metric-like ordering
- low violations with stable refinement behavior -> scale-stable metric-like ordering

## Higher-Level Probes

Only if the basic diagnostics are interpretable:
- persistence-distance ball growth for effective dimensionality
- transport-delay arrival shells for cone ordering
- loop inconsistency and geodesic-deviation analogs for triadic cases

These are diagnostics only.

## Plot Layout

Recommended output figures:
1. `distance_matrix_grid.png`
- base and refined matrix heatmaps for each distance candidate

2. `metric_violation_summary.png`
- symmetry, triangle, and ordering metrics side by side

3. `ordering_flip_heatmap.png`
- pair-rank changes between `12` and `24`

4. `effective_dimensionality_growth.png`
- relational ball-growth curves and fitted exponents

5. `transport_delay_shells.png`
- arrival shells and shell distortion summary

6. `loop_defect_summary.png`
- loop inconsistency and geodesic-deviation analog statistics

## Decision Gate

Stage 15 is positive only if at least one distance candidate shows:
- bounded symmetry deviation
- bounded triangle-inequality violations where applicable
- ordering stability across refinement
- stable or improving higher-level structure

Otherwise the correct conclusion is:
- geometry has not condensed from frozen interaction persistence in the tested pilot

## State

This is a frozen design layer only.
No implementation or results are included here.
