# Stage C0.13 -- Harmonic-Address Topology Phase-Boundary Refinement

## Context

Previous combinatorial-kernel austerity stages established:

- no coordinate distance
- motif, degree, and spectral constraints as the active structural variables
- harmonic-address compatibility as a weak topology-selection bias
- a harmonic detuning continuum with a finite braid-to-smear transition corridor

Stage C0.13 resolves the structure of that transition boundary.
This stage does not test persistence.
It tests phase-topology classification stability and boundary sharpness.

## Objective

Construct a high-resolution scan of harmonic detuning near the observed topology transition band and determine:

1. whether the braid-to-smear transition is sharp, diffuse, or multi-windowed
2. whether the boundary location depends on motif protection strength, degree skew, or weak grade coupling
3. whether a stable topology-phase diagram can be extracted

## Frozen architecture

- discrete combinatorial complex
- frozen interaction-invariant operator sector
- no metric embedding
- no distance kernel
- no memory or delayed reinforcement terms
- no transport stabilizers
- only harmonic-address tagging and internal graded propagation

## Control parameters

Primary axis:

`harmonic_detuning` in a narrow band around the observed transition.

For this stage the corridor is centered on the measured C0.12 threshold near `delta = 0.750`.

Secondary one-axis refinement toggles:

- `motif_protection_strength in {baseline, weak_plus, weak_minus}`
- `degree_skew in {symmetric, mild, asymmetric}`
- `grade_coupling_amplitude in {beta0, betaweak}`

Only one secondary family varies at a time.

## Observables

Per run compute:

- `topology_class`
- `topology_stability_index`
- `flow_concentration_index`
- `address_selectivity_index`
- `grade_exchange_asymmetry`
- `local_channel_count`
- `recirculation_score`

Derived diagnostics:

- `boundary_sharpness_metric`
- adjacent detuning classification variance
- topology-protection sensitivity summaries

## Required outputs

Per-run artifacts:

- stamped JSON
- stamped CSV
- topology vs detuning trace
- selectivity vs detuning trace

Summary artifacts:

- topology phase-diagram panel
- boundary-sharpness heatmap
- protection-sensitivity table

Target note:

`Stage_C0_13_Topology_Phase_Boundary_Refinement_v1.md`

## Boundary rule

A topology phase boundary is considered resolved if:

- at least 3 consecutive detuning samples share an identical topology class
- and adjacent samples outside that band show a different dominant class

Otherwise classify as a diffuse crossover.
