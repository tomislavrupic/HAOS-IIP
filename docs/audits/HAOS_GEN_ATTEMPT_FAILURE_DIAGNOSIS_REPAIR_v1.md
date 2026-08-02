# HAOS-GEN Attempt → Failure → Diagnosis → Repair

## Status

```text
VERIFIED_NEGATIVE_SYNTHETIC_GENERATIVE_LINE_PAUSED
```

This audit closes the current synthetic HAOS-GEN development loop. It does not
close HAOS-IIP as a telemetry or falsification framework.

## 1. Attempt

Attempted claim:

> A sparse, recovery-informed, reversible edge-weight intervention can create
> additional held-out recoverable spectral addresses without degrading the
> admitted set.

Protocol:

- three weighted-path kernel realizations;
- standard and hard held-out regimes;
- two edge changes of magnitude at most `0.08`;
- tuning-only directed selection;
- four equal-budget random interventions per substrate;
- frozen V0.2 address judge;
- old-address retention, new-address uniqueness, null separation, and exact
  reversal.

Reproduction command:

```bash
python3 examples/haos_gen_v04_edge_intervention.py
```

## 2. Failure

Observed:

- six of six directed substrate-regime cells: `DEGRADATION`;
- zero valid directed structural expansions;
- zero valid matched-random structural expansions;
- exact reversal: six of six;
- aggregate:
  `OPEN_NO_STRUCTURAL_INTERVENTION_ADVANTAGE`.

Two cells produced one subspace-distinct post-intervention survivor, but both
failed existing-address retention. No expansion credit was allowed.

The earliest visible warning occurred on `weighted_path_00`: both selected edge
updates had tuning score `0.0`. The fixed intervention budget nevertheless
forced two `+0.08` updates.

## 3. Diagnosis

### D1 — Forced action without positive evidence

The experiment contract required exactly `B` touched edges. It did not include
an abstention gate when no edge showed positive tuning evidence. The system
therefore converted “no informative preference” into an intervention.

### D2 — Selector/acceptance target mismatch

The selector optimized:

```text
mean(admitted-seed recovery + 0.25 * normalized persistence)
```

The acceptance gate protected:

- every pre-intervention V0.2 survivor;
- every admitted seed;
- post-intervention address density;
- subspace-distinct novelty;
- held-out and null separation.

Small gains on the selector proxy did not imply retention of the larger
discovered address set. This mismatch explains why positive tuning scores on
two substrates still produced retention failures.

### D3 — Synthetic search had no external task

V0.3 and V0.4 searched for abstract generative advantage inside synthetic
spectral systems. The framework had a rigorous judge but no external observable
whose recovery mattered independently. This made continued operator invention
technically measurable but weakly motivated.

### D4 — Negative evidence was visually conflated with broken execution

The runs were reproducible, reversible, and contract-valid, but their negative
scientific outcomes were presented like system failures. Status semantics have
now separated verified negative evidence from invalid execution.

## 4. New approach

Pause the current synthetic generative line.

Do not implement another proposal or intervention class until an application
passes the following admission gate:

1. **External task** — a concrete observable failure or recovery problem.
2. **Independent target** — success is defined outside HAOS terminology.
3. **Standard baseline** — comparison against established domain metrics.
4. **Matched controls** — nuisance variables and trivial explanations declared.
5. **False-positive budget** — harm from incorrect flags is measured.
6. **Held-out data** — selection and evaluation are separated.
7. **HAOS contribution** — a predeclared reason recoverability telemetry could
   add information beyond the standard baseline.

One possible application is musical structural-degradation measurement under
masking, compression, resynthesis, stem separation, or edit perturbation.
“AI-music detection” is not admitted at this stage because attribution requires
independent labeled corpora and much stronger false-positive controls.

## 5. Proof draft

### Draft claim — rejected

> HAOS-GEN demonstrates a generative principle that creates new recoverable
> structures.

The evidence does not support this claim. V0.3 found no unique directed
mode-sampling advantage. V0.4 Class 1 produced only degradation.

### Bounded claim — candidate

> In the declared deterministic synthetic suites, the frozen HAOS-GEN judge
> reproducibly distinguishes candidate survival from valid structural
> expansion, rejects interventions that trade existing recovery for apparent
> novelty, and verifies exact reversal of sparse edge updates.

Evidence:

- V0.2 hostile judge is hash-frozen;
- V0.3 rejects baseline-recoverable or seed-specific novelty;
- V0.4 refuses expansion credit when old addresses regress;
- all canonical V0.4 updates reverse to the original measured state;
- repository tests pass reproducibly.

## 6. Adversarial audit

Attack A — **The metrics may only restate their own thresholds.**

Correct. The evidence establishes internal protocol behavior, not external
scientific utility.

Attack B — **Synthetic weighted paths are too narrow.**

Correct. No graph-family generality is claimed.

Attack C — **The admitted seeds are declared, not independently validated.**

Correct. Their role is test-fixture input. “Admitted” does not mean externally
real or physically privileged.

Attack D — **Exact reversal is trivial because the stored matrix update is
subtracted.**

Partly correct. Reversal proves bookkeeping and recovery-profile reproducibility,
not a deep reversibility law.

Attack E — **HAOS may add nothing beyond standard robustness analysis.**

Open. No superiority or uniqueness claim is licensed without an external task
and established baselines.

Attack F — **An abstention gate could avoid degradation but would not create
structure.**

Correct. Abstention is a safety repair, not a generative success.

## 7. Repaired claim

Final claim:

> HAOS-IIP currently provides a reproducible, fail-closed audit protocol for
> recoverability experiments on declared discrete systems. In the tested
> synthetic spectral regime, directed mode sampling showed no unique advantage,
> and sparse local edge interventions failed existing-structure retention. The
> framework has therefore not demonstrated generative advantage. Further
> generative work is paused until a concrete external task passes an
> application-admission gate.

Allowed interpretation:

- the audit machinery is functioning;
- the tested generative hypotheses were unsupported;
- the negative results establish bounded research boundaries;
- an application-grounded experiment may be proposed later.

Forbidden interpretation:

- HAOS has proved a new physical law;
- HAOS generates structure in general;
- synthetic negative results prove all future applications impossible;
- HAOS can identify AI-generated music from degradation telemetry alone.

## 8. Repair applied

- Preserve V0.2 as the frozen judge.
- Preserve V0.3 as
  `TERMINAL_NEGATIVE_MODE_SAMPLING_CURRENT_REGIME`.
- Preserve V0.4 Class 1 as a verified open negative.
- Pause additional synthetic HAOS-GEN versions.
- Route any next experiment through the application-admission gate.
- Use blue/amber for verified negative evidence; reserve red for invalid runs.
