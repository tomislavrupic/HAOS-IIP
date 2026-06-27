# GEO-MM-01 Failure Conditions

## Status

This ledger defines how the Moonshine / Betti geometry-bridge diagnostic fails
or shrinks. It is part of the recoverability contract, not a pessimistic
appendix.

## Threshold Failure

The Betti diagnostic fails if `Betti_0 = 4` is a knife-edge threshold artifact.

Current executable check:

```bash
uv run python experiments/geometry_bridge/monstrous_moonshine_diagnostic/run_betti_threshold_sweep.py
```

Current status:

- `Betti_0` stability band: `[8.0, 11.5]`
- exact edge-signature band: `[8.0, 8.5]`
- verdict: `LOCAL SIGNAL STABLE WITH ROBUSTNESS BAND`, but still local to the declared graph rule

Failure condition:

- `Betti_0 = 4` appears only at a single threshold or disappears under a small
  predeclared threshold perturbation.

## Null-Impostor Failure

The diagnostic fails or weakens if deterministic null ensembles frequently
reproduce the composite signature:

```text
S = {
  Betti_0,
  edge_count,
  component_size_distribution,
  D4/relabel preservation,
  destructive-control response profile
}
```

Current executable check:

```bash
uv run python experiments/geometry_bridge/monstrous_moonshine_diagnostic/run_betti_null_ensemble.py
```

Current status:

- null false-positive rate for the coarse definition: `0.110000`
- verdict: `BETTI_ROBUSTNESS_OPEN`

Failure condition:

- many impostors reproduce the same composite signature, or the false-positive
  rate increases under better-matched null families.

## D4 Symmetry Failure

The Betti diagnostic fails if declared D4 controls alter the canonical edge
signature or Betti vector.

Current executable check:

```bash
uv run python experiments/geometry_bridge/monstrous_moonshine_diagnostic/run_betti_vector.py
```

Failure condition:

- `d4_rotate_90`, `d4_reflect_real`, or graph-isomorphism relabeling changes
  `Betti_0`, `Betti_1`, or canonical edge signature.

## Destructive-Control Failure

The Betti diagnostic is too coarse if destructive controls do not move the
Betti vector or edge signature.

Failure condition:

- `threshold_mutation_control`, `topology_destroyed_control`, or
  `support_replacement_control` returns the reference Betti vector and near-zero
  edge distance.

## Moonshine-Betti Decoupling Failure

The bridge claim fails if Moonshine-channel support perturbations and
Betti-channel support perturbations decouple.

Current executable check:

```bash
uv run python experiments/geometry_bridge/monstrous_moonshine_diagnostic/run_moonshine_betti_bridge.py
```

Allowed phrase:

```text
shared-support perturbation covariance
```

Forbidden upgrade:

```text
Moonshine-Betti mechanism
```

Failure condition:

- support replacement changes only one channel;
- self-recovery moves either channel;
- channel-specific controls are compressed into a derivation claim.

## Formalization Failure

The formalization claim fails if Lean definitions are not local and checkable.

Current executable check:

```bash
uv run python experiments/geometry_bridge/monstrous_moonshine_diagnostic/run_formal_lean_targets.py
```

Current status:

- formal target scaffold: `OPEN`
- local Lean project: not present
- Lean check: not run

Failure condition:

- the repo claims a Lean theorem without local Lean source and a successful
  local check;
- target-only theorem comments are presented as proof;
- placeholder definitions are used to prove a ghost theorem.

## Source-Provenance Failure

The arithmetic diagnostic fails if embedded constants drift from pinned source
claims or fail local arithmetic identities.

Current executable check:

```bash
uv run python experiments/geometry_bridge/monstrous_moonshine_diagnostic/arithmetic_source_validation.py
```

Failure condition:

- supersingular prime list changes without a source update;
- Monster-order exponents drift;
- j-coefficient witnesses no longer sum exactly;
- embedded irrep dimensions lose declared factor-support behavior;
- Gaussian-prime representatives fail the stated residue-class rule.

## Physics Claim Failure

The physics claim is absent until an interaction model exists.

Immediate failure condition:

- any report claims continuum, quantum, gravity, field-theory, or physical
  ontology from the Moonshine / Betti diagnostic suite.

Current allowed classification:

```text
LOCAL SIGNAL STABLE
MATHEMATICALLY INTERESTING
ENGINEERING FEASIBLE
FORMALLY PARTIAL
GLOBAL CLAIM OPEN
```

Current disallowed classification:

```text
Moonshine proof
Monster action recovered
physical bridge established
quantum/gravity/continuum derivation
```
