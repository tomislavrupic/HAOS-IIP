# CST Phase 2.2 Discriminative Instrument Repair

Implemented fact: this repair layer consumes the precommitted contract before generating repaired outputs.
Design choice: hashes name representations only; strict component-distance equivalence determines identity.
Heuristic: component-specific controls are deterministic transforms over existing CST signatures, not new frozen physics runs.
Unverified hypothesis: none added; target victory is not a success criterion.

## Verdict Labels
- INSTRUMENT PARTIALLY DISCRIMINATIVE
- UNDERPOWERED
- TARGET NOT DISTINCT
- MIXED / OPEN

## Power Gate
- status: UNDERPOWERED
- actual independent seeds: 2
- required independent seeds: 4
- actual target-target pairs: 3
- required target-target pairs: 6

## Control Validation
- all predeclared component-control expectations passed

## Signature Collision V2
- coarse collision count: 17 (grouping-only)
- unexplained fine collision count: 0

## Threshold Boundary Sensitivity
- official strict boundary-sensitive pair count: 13
- epsilon: 0.02

## Limitations
- This does not alter frozen HAOS-IIP phases or haos_core APIs.
- This does not add fields because they favor the target; fields come from the static contract.
- The official verdict does not use scalar compression or exploratory weight families.
