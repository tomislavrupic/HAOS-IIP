# Stage C0.3 - Topology-Exchange Robustness Under Incidence Noise

HAOS-IIP Pre-Geometry Atlas - Combinatorial Kernel Line

## Purpose

Stage C0.3 tests whether the no-distance combinatorial exchange surrogate remains stable when the incidence graph is perturbed while global degree statistics remain only weakly affected.

This stage is strictly combinatorial:

- no coordinates
- no Euclidean distance
- no Gaussian kernel
- no continuous embedding metric

The operational question is:

How much of the exchange-like topology survives controlled incidence noise alone?

## Frozen architecture

Freeze:

- field update rule
- graph-shell combinatorial kernel functional form
- damping and drive parameters
- time stepping
- graph size
- packet initialization logic
- grade coupling structure

Only the incidence structure of the edge-graph operator may vary.

## Honest baseline rule

The current C0.1-C0.2 combinatorial branch does not yet contain a clean `braid_like_exchange` baseline.

Therefore C0.3 is anchored to the nearest exchange-like protected surrogate from the no-distance branch:

- representative: `counter_propagating_corridor`
- branch status in C0.2: stable `transfer_smeared` class across uniform, mild, and strong degree families

This keeps the test falsifiable and avoids inventing a braid seed that the branch has not actually produced.

## Perturbation axis

Exactly one axis varies per run:

- incidence noise level and incidence noise pattern

Allowed noise operators:

- random endpoint rewiring
- deletion plus insertion
- localized defect injections
- clustered patch rewiring
- degree-biased incidence perturbation
- connectivity-stress perturbation

All perturbations must preserve:

- graph size
- exact edge count
- connectedness

## Minimal run matrix

Exactly nine runs:

1. baseline no-noise control
2. ultra-weak rewiring
3. weak rewiring
4. moderate rewiring
5. moderate deletion-insertion
6. clustered rewiring patch
7. random localized defect pair
8. degree-biased perturbation
9. connectivity-edge stress test

## Observables

Per run compute:

- `topology_class`
- `braid_like_exchange`
- `smeared_transfer`
- `unresolved`
- `flow_concentration_index`
- `exchange_coherence`
- `channel_count`
- `loop_count`
- `degree_spectrum_shift`
- `clustering_coefficient_shift`
- `recurrence_indicator`

Secondary diagnostics may be retained if already part of the branch.

## Classification logic

A run counts as topology-protected braid only if:

- `topology_class = braid_like_exchange`
- flow concentration remains within a tolerance band of the no-noise baseline
- exchange coherence remains in the protected band
- `loop_count = 0`
- recurrence remains near zero

Otherwise classify as:

- `smeared_transfer`
- or `unresolved`

## Interpretation boundary

This stage may claim only:

combinatorial exchange-like transport is or is not robust to controlled incidence perturbations.

This stage must not claim:

- geometry emergence
- gravity
- particles
- confinement
- universality

## Expected clean reads

A. Strong combinatorial protection

- exchange class survives most perturbations

B. Threshold protection

- exchange survives weak perturbations and fails beyond a noise band

C. Fragile

- exchange surrogate collapses under minimal incidence perturbation

## Required outputs

- stamped JSON result table
- stamped CSV summary
- topology classification summary
- flow concentration panel
- degree spectrum drift panel
- incidence perturbation log
- short atlas note
