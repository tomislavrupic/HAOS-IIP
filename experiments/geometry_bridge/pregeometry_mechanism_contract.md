# GEO-05 Pregeometry Mechanism Contract

## Status

OPEN

## Purpose

Freeze the missing dependency between intrinsic geometry, transformation semantics,
probability rule, and observable prediction before any Bell interpretation is
attempted.

## Dependency Chain

`intrinsic_relations -> relational_geometry -> transformation_semantics -> probability_rule -> observable_prediction`

## What Must Hold

- Intrinsic relations are recovered from operators, transport, or closure behavior.
- Relational geometry exposes distance, orientation, adjacency, and covariance.
- Transformation semantics gives settings a native meaning and composition rule.
- Probability rule maps structure to bounded predictions without target import.
- Observable prediction improves on null baselines on a frozen holdout.

## What Is Frozen

- Source artifacts from GEO-01 through GEO-04
- The repository geometry bridge notes and README
- No Bell scoring
- No physical mechanism claim
- No target-curve import
- No post hoc geometry selection

## Gate Order

1. GEO-01 intrinsic geometry recovery
2. GEO-02 transformation semantics
3. GEO-03 probability-rule recovery
4. GEO-04 observable prediction

## Promotion Rule

- If holdout fails, status remains OPEN.
- If only pairwise observable structure succeeds, status remains OPEN.
- If coarse observable target beats null, status remains OPEN.
- Only if all synthetic gates pass is the chain marked complete.

## Boundary

This mechanism is a frozen synthetic dependency contract. It is the missing
intermediate layer between geometry and downstream observables, not a Bell
mechanism and not a physical derivation.
