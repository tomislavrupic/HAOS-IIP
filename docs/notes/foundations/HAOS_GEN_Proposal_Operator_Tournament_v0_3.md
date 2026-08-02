# HAOS-GEN V0.3 — Proposal Operator Tournament

## Frozen judge

V0.3 does not change candidate acceptance. The V0.2 judge is frozen at:

```text
haos_gen/hostile_validation.py  7409b6b96cc213fcaf30237bb5686bcd146be0c0effd4ee89084890a2fa56d8c
haos_gen/spectral_tuning.py     e69eb0123159fb2a94913e47cad6472a04f340157d384d7841f7a89ab33df063
telemetry/frozen_metrics.py     daa7a759ad8f0ed67249f194a7e7e0753c78a4eb4c0a8b44d47e2e7a0d121f87
```

The test suite fails if that byte hash changes. V0.3 imports the frozen
recovery, amplification, random-baseline, and hostile-null helpers and adapts
arbitrary candidate addresses to the same V0/V0.2 gates.

## Shared tournament contract

Each operator receives the same integer proposal budget for every seed and
regime. Insufficient distinct proposals are a hard error; the tournament never
silently lowers a budget.

The canonical demo declares:

- five admitted seed addresses;
- four proposals per operator, seed, and regime;
- standard and hard regimes;
- identical tuning stress for all;
- sealed long-range held-out stress in the standard regime;
- stronger mixed-range disorder in the hard regime;
- V0 local resonance and non-overlapping random eigenmode baselines;
- the frozen coordinate, degree-class, and spectrum-preserving null suite.

## Operator A: non-local spectral mixing

Version: `non_local_spectral_mixing_v0_3`.

1. Compute the seed Rayleigh coordinate in the declared symmetric operator.
2. Admit eigenmodes whose normalized spectral separation exceeds
   `min_spectral_separation`.
3. Rank modes by tuning-family mean recovery, then spectral separation.
4. Prefer modes whose mean recovery is at least the seed mean; fall back to the
   full distant set only when none qualify.
5. Generate normalized mixtures using the fixed alpha set `(0.3, 0.5, 0.7)`.

No locality term enters selection.

## Operator B: multi-seed interference

Version: `multi_seed_interference_v0_3`.

1. Select fixed-size seed groups with deterministic identifier order.
2. Apply the versioned `equal_weight_signed` combination rule.
3. Normalize the result.
4. Reject combinations whose maximum overlap with any source seed exceeds the
   declared input-overlap limit.
5. Deduplicate up to sign before taking the fixed budget.

The canonical experiment uses `k=2`.

## Operator C: recovery-guided graph walk

Version: `recovery_guided_graph_walk_v0_3`.

1. Infer weighted adjacency from the off-diagonal Laplacian entries.
2. Evaluate every node delta-address on the tuning family.
3. Define the local recovery field as:

```text
0.75 * normalized mean recovery
+ 0.25 * normalized persistence
```

4. Start from seed-support nodes in descending amplitude order.
5. At each step choose the non-revisited neighbor maximizing:

```text
walk_bias_recovery_weight * local recovery
+ (1 - walk_bias_recovery_weight) * normalized edge weight
```

6. At fixed stopping steps, expand support by the declared radius and lift the
   first non-constant restricted-operator eigenvector into the full space.

The walk is deterministic. There is no undeclared sampling.

## Uniqueness gate

A directed held-out survivor is unique only when, against every surviving V0
local and random baseline address for the same seed and regime:

1. normalized Rayleigh-coordinate distance is greater than
   `address_distance_threshold`; and
2. absolute address overlap is below `subspace_overlap_threshold`.

If either condition fails, the ledger records `BASELINE_RECOVERABLE`.

## Primary survivor-discovery ledger

Every row records:

- operator, seed, regime, and candidate;
- tuning and held-out decisions;
- recovery, persistence, and `k_star` changes;
- coordinate, degree, and spectrum-null margins;
- baseline match;
- uniqueness;
- explicit failure notes.

Aggregate status is positive only when at least one directed operator produces
unique survivors across the declared minimum number of seeds and in the hard
regime. Otherwise:

```text
OPEN_NO_UNIQUE_DIRECTED_ADVANTAGE
```

## Canonical V0.3 result

The first five-seed, two-regime run closes:

```text
OPEN_NO_UNIQUE_DIRECTED_ADVANTAGE
```

Observed:

- V0 local: 22/40 survivors;
- random eigenmodes: 35/40 survivors;
- non-local mixing: 3/40 survivors, one unique in standard only;
- multi-seed interference: 1/40 survivor, unique in standard only;
- recovery-guided walks: 0/40 survivors;
- no directed operator produced a unique hard-regime survivor across multiple
  seeds.

These are bounded deterministic demo results, not population estimates.

## Remaining gaps

- The uniqueness distance uses a normalized Rayleigh coordinate plus
  single-vector overlap. Degenerate subspaces require principal-angle tracking.
- “Admitted seeds” are declared synthetic addresses in the demo; a later
  authority bundle must prove admission from independent ledgers.
- Random comparison is one deterministic equal-budget draw per seed/regime.
  Repeated draws and confidence intervals remain OPEN.
- Non-local mixing still uses base-operator eigenmodes rather than learned
  operator modifications.
- Multi-seed interference uses pairwise equal weighting only.
- The graph walk uses node-delta recovery as a local proxy. That proxy may be
  too weak or misaligned with extended spectral recovery.
- The hard regime is stronger and structurally mixed, but it is still synthetic.
- Status: **OPEN proposal tournament with no unique directed advantage**.
