# HAOS-GEN V0.2 — Hostile Validation

## Purpose

V0 demonstrated that a directed spectral generator could propose candidates and
reject regressions. V0.2 tests whether that selection survives outside the
perturbation family used during tuning.

The governing rule is:

> Search may occur on the tuning family. Retention requires amplification on a
> sealed, structurally distinct held-out family and separation from every
> declared hostile null.

This is still an experimental operator audit. It is not evidence for a physical
generative law.

## Two-family protocol

The caller declares two perturbation suites before the run:

- **tuning family**: used by the V0 generator and initial selection gate;
- **held-out family**: never used to rank or construct candidates.

The reference demo uses nearest-edge weight/diagonal perturbations for tuning
and next-nearest-neighbor/oscillatory-diagonal perturbations for held-out
validation. They are deterministic but structurally different.

A candidate must:

1. pass the V0 tuning gate;
2. pass recovery, persistence, and `k_star` gates on held-out stress;
3. beat the worst held-out mean recovery among all declared null families;
4. show at least one positive held-out amplification:
   - mean recovery gain,
   - persistence gain, or
   - later/no sustained collapse relative to the seed.

Non-regression alone is no longer sufficient.

## Hostile nulls

Each candidate is compared under the same held-out stress coordinates against:

1. **coordinate-permutation null** — preserves coefficient multiset and norm;
2. **degree-class permutation null** — permutes coefficients among nodes with
   matching operator diagonal/degree class;
3. **spectrum-preserving basis-scramble null** — preserves every operator
   eigenvalue while destroying the original eigenbasis/address relation.

The candidate must exceed the strongest null by the predeclared hostile-control
margin. These nulls are intentionally different attacks; none is treated as a
complete specificity proof.

## Equal-budget random baseline

V0.2 draws a deterministic seeded sample of pure random eigenmodes with the same
count as the directed proposal set. Modes already present in the directed set
are excluded from this pool. Each random mode receives the same:

- tuning non-regression gate;
- held-out non-regression gate;
- hostile null suite;
- amplification requirement.

The report records directed and random survival yields separately. A negative
directed-yield advantage is a valid result and must not be hidden. Probe status
is fail-closed:

- `FAIL_NO_HOSTILE_SURVIVORS`
- `OPEN_NO_DIRECTED_YIELD_ADVANTAGE`
- `PASS_BOUNDED_DIRECTED_YIELD_ADVANTAGE`

## Yield and cost accounting

Every run reports:

- proposed count;
- tuning survivors;
- held-out survivors;
- final hostile-validation survivors;
- cumulative yield after each proposal;
- proposals per final survivor;
- equal-budget random survivors;
- directed minus random yield.

This is the first explicit generation-cost proxy. It measures search yield, not
compute cost, energy, or universal creative efficiency.

## Integration

- implementation: `haos_gen/hostile_validation.py`;
- V0 metrics and proposal operator: `haos_gen/spectral_tuning.py`;
- frozen telemetry: `telemetry/frozen_metrics.py`;
- executable example: `examples/haos_gen_hostile_demo.py`;
- tests: `tests/test_haos_gen.py`.

V0.2 remains downstream of Phase XIX and outside the frozen `haos_core` API.
No existing phase manifest is edited. Promotion to a numbered measurement phase
would require a new authority bundle, checker, runs ledger, and claim gate.

## Failure modes

- `TUNING_GATE_FAILED`
- `HELDOUT_RECOVERABILITY_GATE_FAILED`
- `HOSTILE_NULL_SEPARATION_FAILED`
- `NO_HELDOUT_AMPLIFICATION`

All failures remain in the output ledger.

## Remaining gaps

- The held-out family is sealed by API contract and reporting discipline, not
  cryptographic commitment or a separate evaluator process.
- The degree-class null uses the operator diagonal as a graph-degree proxy. It
  does not preserve all graph motifs, community structure, or weighted-degree
  correlations.
- The spectrum-preserving null preserves eigenvalues but produces dense
  operators. Sparse spectrum-preserving graph nulls remain OPEN.
- Random eigenmodes test proposal directionality only against one seeded sample.
  Multi-seed confidence intervals and repeated equal-budget trials remain OPEN.
- Yield counts proposals, not wall time, memory, or perturbation-evaluation cost.
- Candidate generation is still local to the seed-supported eigensystem.
- Degenerate eigenspaces and eigenvalue crossings still require subspace
  tracking.
- A passing held-out run is bounded to the declared operator and two stress
  families. It does not establish generalization beyond them.
- Status: **OPEN hostile-validation prototype**.
