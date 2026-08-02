# Expansion Manifesto

HAOS-IIP learned how to distrust structure well. It detects degradation,
identifies sustained collapse, preserves frozen regimes, and refuses to turn a
surviving diagnostic into an ontology. That discipline remains the authority.
What was missing was an equally explicit rule for producing new testable
structure.

HAOS-GEN adds that rule: expand from what is already recoverable, then make the
expansion earn retention under the same hostile telemetry as its parent.
Generation is neither acceptance nor discovery. It is the controlled production
of candidates. Selection belongs to HAOS.

QATC contributes the intuition of tuning and selective amplification. Here those
words are reduced to finite operators: resonance is normalized overlap; tuning
is directed traversal of a declared spectral neighborhood; amplification is
retention of candidates that maintain or improve recovery under a perturbation
ladder. Nothing in this layer licenses a quantum, physical, creative-field, or
ontological claim.

# Technical Specification

## Generative principle

Given a frozen or declared symmetric operator \(L\), an already admitted seed
address \(a_0\), and a predeclared perturbation suite
\(\{L(\epsilon_i)\}\), HAOS-GEN:

1. diagonalizes \(L\) without changing it;
2. proposes nearby eigenmode addresses \(v_j\), ordered by
   \(R_j \times Q_j\), where \(R_j=|\langle a_0,v_j\rangle|\) is resonance and
   \(Q_j=\sum_n |v_j(n)|^4\) is a locality/concentration preference;
3. reidentifies every proposed address in each perturbed operator by maximum
   overlap;
4. evaluates the reconstructed trajectory with the frozen
   `recovery_score` and `persistence_time`, plus the declared first sustained
   collapse index `k_star`;
5. runs the same audit on a matched coordinate-permutation null preserving the
   candidate coefficient multiset and norm;
6. retains only candidates with no material recovery regression, no persistence
   regression, no earlier `k_star`, and declared separation from the matched
   null.

The output is an expanded set of *candidate spectral coordinates*. It is not a
new phase authority, a modified frozen operator, or evidence for new physics.

## Contracts

Input:

- one finite symmetric operator;
- one nonzero seed address of matching dimension;
- an ordered, equal-shape perturbation ladder and stress coordinates;
- frozen `SurvivalThresholds`;
- predeclared generation and selection thresholds.

Output:

- baseline recovery profile;
- candidate and matched-null recovery profiles;
- resonance, locality, novelty, `recovery_score` series, `persistence_time`,
  and `k_star`;
- accepted/rejected status with explicit failure labels;
- retained address vectors;
- no mutation of inputs, telemetry, manifests, or frozen phase artifacts.

## Failure modes

- `RECOVERABILITY_REGRESSION`: mean recovery falls beyond tolerance.
- `PERSISTENCE_REGRESSION`: the candidate fails the frozen survival gate sooner.
- `EARLIER_SUSTAINED_COLLAPSE`: candidate `k_star` precedes baseline `k_star`.
- `NO_MATCHED_CONTROL_SEPARATION`: the tuned candidate does not beat its null.
- `NO_RECOVERABLE_EXPANSION`: neither improvement nor stable novelty is present.
- hard input failure: non-square/asymmetric operators, dimension mismatch,
  empty perturbation suite, or unordered stress coordinates.

## Integration

- `haos_gen/spectral_tuning.py` is a downstream experimental layer.
- It imports frozen functions from `telemetry/frozen_metrics.py`; it does not
  duplicate or edit their formulas.
- It does not enter `haos_core`, whose public API remains frozen.
- Phase XIX supplies the bounded term *spectral address*. HAOS-GEN supplies the
  first experimental dynamics attached to that term.
- A later numbered measurement phase may consume HAOS-GEN output only through a
  new manifest, checker, runs ledger, and claim gate. No existing phase manifest
  should be retroactively edited.
- `examples/haos_gen_demo.py` is the minimal reproducible path:
  generation -> perturbation ladder -> frozen telemetry -> null comparison ->
  selection -> expanded address set.

# Remaining Gaps

- The first operator proposes 0-cochain/eigenvector addresses only. General
  cochain and Hodge-sector proposals are OPEN.
- The matched null preserves coefficient values and norm, but not every graph
  statistic. Degree-, locality-, band-, and spectrum-matched null families are
  still required before stronger specificity language.
- Maximum-overlap eigenvector tracking can fail at degenerate or crossing
  eigenspaces. Subspace tracking and degeneracy-aware addresses are OPEN.
- Candidate ranking uses overlap and inverse-participation-style locality. It is
  a deterministic heuristic, not an optimality theorem.
- Stress suites are supplied by the caller. Results are conditional on their
  coverage and can be overfit if generation and final evaluation reuse the same
  ladder. A held-out hostile perturbation suite is required for serious claims.
- `k_star` is implemented locally because the current frozen telemetry module
  does not export one canonical cross-phase function. This implementation is
  explicit and tested, but it is not promoted into frozen telemetry.
- Temporal consistency, causal deformation, and geometric integrity are not
  silently synthesized for a purely spectral address. They remain
  not-applicable until a candidate has declared temporal, causal, or geometric
  observables.
- No empirical advantage, continuum result, physical emergence, quantum
  mechanism, or general theory of creativity is claimed.
- Status: **OPEN experimental layer**. A passing demo establishes executable
  bookkeeping, not scientific validation.
