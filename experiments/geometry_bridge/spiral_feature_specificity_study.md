# Spiral Feature Specificity Study

Status: frozen precommitment

Purpose

Determine which geometric invariant is necessary for spiral-family transfer
and whether it is already latent in HAOS-IIP or would need to be added as
genuinely new structure.

Frozen residual from GEO-01

`holdout transfer does not outperform the best baseline on the spiral family`

Candidate invariants to precommit

1. Curvature / turning rate
2. Geodesic ordering
3. Orientation transport
4. Winding / path history
5. Multiscale neighborhood structure

Decisive comparison

- best baseline
- vs best baseline + precommitted HAOS invariant

Matched controls

- degree-preserving rewiring
- parameter-matched nulls
- topology destruction
- label permutation
- baseline + candidate invariant
- holdout spiral transfer

What this study will not do

- it will not change the GEO-01 thresholds
- it will not tune the current HAOS geometry score against spiral holdout
- it will not claim external physics meaning
- it will not add features before the invariant is frozen
- it will not optimize the spiral score directly

Success criteria

- the specific missing invariant is isolated under a frozen precommitment
- the added invariant improves holdout spiral transfer against the best
  baseline and the matched control destroys that increment
- if the invariant is already latent, the study records the operational mapping
  without upgrading it into a physics claim
- if no candidate produces stable holdout improvement, freeze the result as
  `SPIRAL_SPECIFIC_STRUCTURE_NOT_RECOVERED`

Failure criteria

- the study quietly re-optimizes GEO-01
- the study adds features without precommitting them
- the study upgrades synthetic specificity into physical semantics
- the study chases score changes without a matched-control check
