# Stage C0.2 - Degree-Spectrum Protection Test

## Purpose

Test whether transport, exchange, or topology-like flow classes in the combinatorial no-distance branch are protected, suppressed, or reshaped by the degree spectrum of the underlying graph.

This stage remains fully combinatorial:
- no Euclidean distance
- no Gaussian kernel
- no embedded metric oracle
- no dynamic rewiring during the run
- no reinforcement memory

The core question is:

Are the observed flow classes protected by degree regularity, or do they depend on degree defects, hubs, or shell inhomogeneity?

## Frozen Ingredients

Keep fixed from C0.1:
- packet update law
- projected-transverse packet family
- packet amplitude and phase initialization conventions
- topology and coherence observables
- no adaptive kernel evolution

Only change:
- graph degree-spectrum structure

## Graph Families

Use a minimal degree-spectrum set:
1. degree-uniform baseline
2. mild degree spread
3. strong degree heterogeneity

Optional later extension:
4. bimodal degree graph

## Kernel Constraint

Keep the combinatorial kernel family fixed while scanning the degree spectrum.

Allowed ingredients:
- adjacency
- degree normalization
- bounded hop-radius shell weights
- motif or shared-neighbor structure

Forbidden:
- coordinates
- Euclidean distance
- Gaussian embedding distance

## Representatives

Use the smallest useful inherited set:
1. clustered encounter
2. corridor-like counterflow
3. symmetric exchange-pair surrogate

Do not open a new operator sector in C0.2.

## Primary Observables

For each run compute:
- topology label
- flow concentration index
- transport span
- coherence metric
- overlap persistence
- degree-localization correlation
- defect pinning score
- hub attraction or avoidance score

Where relevant, distinguish among:
- transfer_smeared
- trapped
- dispersive
- defect_pinned
- channel_split

## Protection Logic

A topology counts as degree-protected if:
- it survives across the degree-uniform and mild-spread graphs
- classification remains stable under small degree perturbations
- flow concentration and coherence stay inside the same class envelope

A topology counts as degree-sensitive if:
- it collapses under mild spread
- it changes class near degree defects
- transport is redirected, smeared, or pinned by hubs or low-degree pockets

A topology counts as degree-induced if:
- it does not exist in the uniform graph
- but appears reproducibly only in heterogeneous degree structures

## Minimal Run Design

Use the compact first-pass matrix:
- 3 graph families
- 3 representative encounters

Total:
- 9 runs

No parameter sweep beyond this in the first pass.

## Required Outputs

Per run:
- stamped JSON
- CSV row
- graph-flow visualization
- node-degree overlay
- topology classification
- coherence and span traces

Summary outputs:
1. topology vs degree-spectrum matrix
2. flow concentration by graph family
3. degree-localization correlation panel
4. hub and defect pinning summary
5. short interpretation note

Target note:

`Stage_C0_2_Degree_Spectrum_Protection_Test_v1.md`

## Interpretation Boundary

Do not claim:
- geometry emergence
- gravity
- physical particle trapping
- universality

At most claim:

combinatorial transport topology is or is not protected by graph degree-spectrum structure.
