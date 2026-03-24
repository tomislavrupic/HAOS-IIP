# Continuum Sketch — Minimal Scaling Closure Protocol (Post-Phase XVIII)

## Program context

You are working inside the frozen HAOS-IIP program after Phase XVIII.

The following are already frozen and must not be modified:
- Phase VI cochain-Laplacian operator hierarchy
- Phase X proto-continuum extension law
- Phase XI localized-mode persistence diagnostics
- Phase XV propagation observables
- Phase XVI temporal-ordering observable
- Phase XVII causal-graph extraction
- Phase XVIII distance-surrogate construction
- frozen telemetry formulas (`telemetry/frozen_metrics.py`)

This protocol does not introduce new dynamics, new ontology, or new observables.

It tests only whether a small subset of existing observables exhibits simple refinement-ordered scaling.

## 0. Objective

Determine whether the frozen branch regime admits a proto-continuum scaling closure in the sense that:
- selected observables vary smoothly with refinement
- simple scaling exponents remain bounded
- an effective propagation bound remains stable
- causal-distance growth remains coherent

No continuum interpretation is claimed.

## 1. Refinement slice

Use a minimal deterministic hierarchy:

`n_side = { 48, 60, 72, 84 }`

Use the same frozen branch seeds already validated in Phases XVII–XVIII.
Use the matched altered-connectivity control hierarchy.

No new refinement levels may be added.

## 2. Observable set (fixed)

Evaluate only the following five quantities:
1. low-mode dispersion proxy  
   (e.g. sampled spectral radius or equivalent frozen Phase VII metric)
2. effective propagation speed band  
   (Phase XV disturbance-front observable)
3. localized-mode persistence time  
   (Phase XI metric)
4. temporal-ordering consistency score  
   (Phase XVI metric)
5. distance-surrogate shell slope  
   (Phase XVIII refinement descriptor)

No additional observables may be introduced.

## 3. Computation rules

For each refinement level:
- compute the five observables on the frozen branch
- compute the same observables on the matched control
- store results in a single structured table

Then:
- perform simple log-log or linear-in-h fits
- estimate drift bands or exponent spans
- evaluate propagation-speed variance across refinement
- check monotonicity of the distance-surrogate slope

No high-order fitting, symbolic modeling, or parameter scanning.

## 4. Required artifact set

Produce exactly:
1. `continuum_sketch_table.csv`  
   rows: refinement × observable × hierarchy (branch/control)
2. `continuum_sketch_plot_family.svg`  
   one compact multi-panel figure showing observable vs refinement
3. `continuum_sketch_summary.md`

No additional plots, logs, or diagnostics.

## 5. Success criteria

The sketch is considered feasible if:
- observable drift remains bounded across refinement
- at least two observables show smooth scaling trends
- effective propagation speed variance stays within a compact band
- branch/control separation persists at mesoscopic scale

The sketch is considered inconclusive if:
- scaling behavior is irregular
- propagation band broadens strongly
- ordering or distance metrics lose refinement coherence

## 6. Closure paragraph (must end exactly with)

If feasibility holds:

The minimal scaling sketch indicates bounded refinement-ordered behavior of key frozen observables. This is consistent with a proto-continuum effective description within the tested branch regime. No continuum ontology or universality is asserted.

If not:

The minimal scaling sketch does not yet show stable refinement-ordered behavior sufficient for a proto-continuum effective description within the tested branch regime.

## 7. Working style

Keep the execution extremely light:
- reuse frozen data paths whenever possible
- avoid rebuilding full histories
- avoid new stochastic probes
- prefer deterministic slices already validated

This protocol is a compression test, not a new phase.
