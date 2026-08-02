# HAOS-GEN V0.4 — Structural Intervention

## Class 1: local edge-weight intervention

V0.3 is terminal-negative for mode sampling in the current spectral regime:

```text
TERMINAL_NEGATIVE_MODE_SAMPLING_CURRENT_REGIME
```

V0.4 changes the object being generated. It applies a sparse, reversible update
to the substrate and asks whether the modified operator supports genuinely new
held-out recoverable addresses without degrading its admitted set.

The V0.2 address judge remains frozen.

## Frozen intervention contract

The canonical run fixes:

- edge budget `B = 2`;
- maximum per-edge change `delta_w_max = 0.08`;
- total L1 budget `<= 0.16`;
- retention tolerance `epsilon_ret = 0.02`;
- reversal tolerance `epsilon_rev = 1e-10`;
- subspace-equivalence overlap `theta_eq = 0.92`;
- four matched random interventions per substrate;
- three distinct weighted-path kernel realizations;
- standard and hard held-out regimes.

For an edge `(u, v)` with signed weight change `delta`, the Laplacian update is:

```text
delta * (e_u - e_v)(e_u - e_v)^T
```

The exact update records are sufficient to reverse the intervention.

## Directed selector

The selector sees tuning data only.

For every existing graph edge it evaluates both:

- reinforcement by `+delta_w_max`;
- suppression by `-min(delta_w_max, current_weight)`.

For each trial it evaluates the declared admitted seeds with frozen
`recovery_score` and `persistence_time`. The tuning objective is:

```text
mean(seed mean_recovery + 0.25 * normalized seed persistence)
```

The better sign is retained for that edge. Edges are ranked by objective gain,
with deterministic endpoint tie-breaking, and the top `B` are selected.

No held-out profile enters this computation.

## Address discovery and retention

Pre- and post-intervention address sets are discovered by running the frozen
V0.2 hostile validator from the same admitted seeds.

Post-intervention perturbation operators receive the same sparse Laplacian
update as the substrate. Stress coordinates and perturbation patterns remain
unchanged.

Existing-address retention requires:

- every pre-intervention V0.2 survivor to pass the V0.2-compatible post gate;
- no mean-recovery or persistence regression beyond `epsilon_ret`;
- no declared admitted seed to regress beyond `epsilon_ret`.

New addresses receive credit only when they:

- survive post-intervention V0.2 tuning, held-out, null, and amplification gates;
- have overlap below `theta_eq` with every pre-intervention survivor.

Candidate identifiers alone never establish novelty.

## Matched random baseline

Each random trial:

- touches exactly `B` graph edges;
- uses the identical signed update-magnitude multiset as the directed run;
- samples edges with a recorded deterministic random seed;
- faces the same discovery, retention, subspace, and held-out gates.

A random run that creates an address while degrading the old set is not counted
as a valid random structural expansion.

## Reversibility

The exact update matrix is subtracted from:

- the modified substrate;
- every modified tuning operator;
- every modified held-out operator.

V0.4 then re-runs address discovery. Matrix equality, address count, mean
recovery, mean persistence, and mean `k_star` coordinate must return within
`epsilon_rev`.

## Outcome labels

- `STRUCTURAL_EXPANSION`: retained original set, at least one new
  subspace-distinct V0.2 survivor, strict advantage over valid matched random
  interventions, and successful reversal.
- `STABILIZATION`: no new address, but strict recovery, persistence, or `k_star`
  improvement beyond valid random interventions.
- `RELABELING / NULL`: no net structural gain, failed reversal, equivalent
  subspace, or random-matched result.
- `DEGRADATION`: admitted seeds or existing survivors regress, or post address
  density falls.

Aggregate confirmation requires structural expansion on multiple substrates
and on those substrates in the hard regime. Otherwise:

```text
OPEN_NO_STRUCTURAL_INTERVENTION_ADVANTAGE
```

## Canonical result

The first run closes:

```text
OPEN_NO_STRUCTURAL_INTERVENTION_ADVANTAGE
```

All six directed substrate-regime cells are classified `DEGRADATION`.

- Two cells produced one subspace-distinct post survivor each, but existing
  addresses regressed, so no expansion credit was allowed.
- The remaining four cells lost address density or produced no unique survivor.
- Exact reversal passed in all six cells.
- No matched random trial simultaneously retained the old set and produced a
  unique new address.

This is a bounded synthetic negative, not a theorem that edge intervention
cannot work.

## Remaining gaps

- The selector optimizes admitted-seed conditioning, not expected new-address
  density directly.
- Only one edge budget and update magnitude are tested.
- Weighted paths are distinct kernel realizations but remain one graph family.
- V0.2 mode discovery limits the observable address pool.
- Degenerate subspaces still require principal-angle rather than vector-overlap
  equivalence.
- Degree-preserving rewiring, sparse low-rank updates, and cochain edits remain
  outside V0.4 Class 1.
- Status: **OPEN with no structural-intervention advantage**.
