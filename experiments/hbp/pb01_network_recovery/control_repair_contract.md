# PB-01 Control Repair Contract

Status: frozen repair contract

Objective

Repair PB-01 controls so the benchmark can reach `CONTROL_VALID` without
changing the frozen holdout result or chasing a higher HAOS score.

Frozen failure to repair

- leakage positive control is detected
- destructive controls are too weak
- target remains predictive even when topology, degree, weights, or null
  structure are destroyed under the current contract

Component-specific control contracts

## 1. Topology-destroyed graph

Preserve

- node count
- weight scale range
- perturbation metadata

Destroy

- adjacency topology
- shortest-path structure
- connectivity structure
- local neighborhood order

Expected metric response

- recovery-quality shift should be large enough to exceed the frozen absolute
  movement threshold
- degree-based and spectral baselines should no longer track the target closely

Invalid if

- the control leaves recoverability nearly unchanged
- the control still preserves the tested neighborhood signal

## 2. Degree-preserving rewiring

Preserve

- degree histogram
- node count
- weight scale

Destroy

- edge placement
- adjacency motifs
- long-range recovery pathways

Expected metric response

- the target should separate from degree-only predictors if topology is
  actually contributing
- recovery-quality shift should exceed the frozen threshold

Invalid if

- rewiring behaves like the original graph on recovery ranking
- degree-only structure still explains the signal

## 3. Weight-shuffled graph

Preserve

- adjacency support
- node count
- perturbation metadata

Destroy

- edge-weight ordering
- weighted recovery pathways
- spectral weight structure

Expected metric response

- performance should drop if the HAOS proxy is really responding to weighted
  structure rather than adjacency only

Invalid if

- recovery quality barely changes
- weighted-signal duplicates remain visible

## 4. Parameter-matched null

Preserve

- graph size
- perturbation family
- rough severity range

Destroy

- graph-specific recovery structure
- topology-specific ranking signal

Expected metric response

- should be the weakest of the destructive controls
- still must move enough to qualify as destructive under the frozen contract

Invalid if

- the null remains too close to the target
- the null is indistinguishable from the original benchmark

## 5. Perturbation-free baseline

Preserve

- original graph and dynamics

Destroy

- none

Purpose

- sanity check that the system can recover in the absence of perturbation

Invalid if

- the baseline is worse than a perturbed case

## 6. Seed repeat

Preserve

- all construction choices

Destroy

- none

Purpose

- verify determinism

Invalid if

- reruns change the case count or the control ordering

## 7. Intentional leakage positive control

Preserve

- direct access to the target signal

Destroy

- blindness

Purpose

- verify the leak detector

Invalid if

- the control is not detected as leakage

Repair checks required before any new PB-01 score claim

- compare HAOS ranks against degree, spectral gap, closeness, and shortest-path
  baselines on the same frozen split
- verify whether perturbation severity can be increased without breaking the
  leak detector
- verify that graph-family shortcuts are not carrying the target
- rerun controls only after the destruction criteria are frozen

Boundary

This contract changes control validity criteria only. It does not change the
benchmark target, thresholds, or claim ceiling.
