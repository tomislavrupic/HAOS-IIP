# HAOS-IIP Falsification Engine v1

Status:

- Phase 66 artifact
- failure-gate checklist
- no new result
- no stronger physical claim

Purpose:

Define public failure conditions for HAOS-IIP claims, dashboards, bridge
observables, and phase releases.

Core rule:

```text
A HAOS-IIP claim is not mature until it states how it can fail.
```

This artifact turns critique into executable failure gates. It does not weaken
the program. It prevents a claim from staying vague enough to avoid contact with
tests.

## 1. Claim Classes

Every claim should declare exactly one primary class. Secondary classes can be
listed only after the primary failure gates are satisfied.

| Class | Claim type | Short definition | Minimum evidence class |
| --- | --- | --- | --- |
| A | Internal numerical claim | A result passes inside a declared frozen computational regime | Reproducible run, fixed metrics, fixed thresholds |
| B | Structural stability claim | A diagnostic structure survives bounded perturbation | Perturbation sweep plus matched controls |
| C | Scaling / refinement claim | A structure remains stable or converges across size/refinement | Multi-scale run plus drift / convergence test |
| D | Geometry / metric surrogate claim | Operator response behaves like a distance or metric surrogate | Positivity, locality, isotropy, convergence, controls |
| E | Current / closure claim | Flux-like or conservation-like residual remains bounded | Closure residuals, shell / region checks, null controls |
| F | Physics-bridge observable | Internal metric maps to a physics-adjacent toy/proxy/external observable | Claim gate, null comparison, allowed-language boundary |
| G | Biology telemetry claim | Recoverability telemetry applies to biological or bio-adjacent data/model | Biological null model, data provenance, no mechanism overclaim |
| H | Materials bridge claim | Recoverability telemetry applies to public materials data | Provenance, parsing, hashes, external-data metrics, no embodiment claim |

## 2. Required Failure Gate Per Claim

Every claim must specify:

- what is measured
- what baseline is used
- what perturbation is applied
- what null model is compared
- what threshold counts as failure
- what result would downgrade the claim
- what result would kill the claim
- what reproduction command or artifact verifies the result

If one of these is missing, the claim remains immature.

## 3. General Failure Gates

These gates can apply to any claim class:

- fails under seed variation
- fails under bounded perturbation
- fails under refinement
- fails against matched controls
- cannot beat a declared null model
- does not reproduce from committed commands and artifacts
- depends on hand-picked parameters
- loses positivity under the claimed regime
- loses isotropy under the claimed regime
- loses closure under the claimed regime
- bridge observable is indistinguishable from generic graph behavior
- external-data result disappears when standard analysis or null models are used
- claim language exceeds the evidence class

## 4. Claim-Specific Failure Gates

| Claim class | Downgrade if | Kill if |
| --- | --- | --- |
| A. Internal numerical claim | thresholds are met only narrowly or controls are incomplete | reproduction fails or thresholds are altered after seeing results |
| B. Structural stability claim | stability is seed-specific or parameter-sensitive | structure collapses under the declared bounded perturbation |
| C. Scaling / refinement claim | trend is noisy or only two scales are tested | structure does not converge or relabels under refinement |
| D. Geometry / metric surrogate claim | metric-like behavior holds only locally or anisotropically | positivity, locality, isotropy, or convergence fails in the claimed regime |
| E. Current / closure claim | residuals are larger than expected or only one branch closes | closure disappears under matched controls, refinement, or alternate discretization |
| F. Physics-bridge observable | toy/proxy controls remain competitive | observable cannot beat generic graph/null behavior or claim gate is exceeded |
| G. Biology telemetry claim | biological metadata are incomplete or effect is weak | no advantage over standard biological null models or mechanism is unsupported |
| H. Materials bridge claim | data parse is partial or metadata are incomplete | no real data, no provenance, parse failure, or no advantage over standard materials analysis |

## 5. Claim Downgrade Ladder

Use this ladder when a result weakens under audit.

```text
Confirmed inside frozen regime
  ↓
Numerically suggestive
  ↓
Analogy only
  ↓
Unresolved
  ↓
Failed under current tests
```

Downgrading a claim is not a failure of discipline. It is the discipline working.

## 6. Mandatory Claim Box

Full claim-gating artifact:

- [HAOS_IIP_Claim_Gating_Template_v1.md](/Volumes/Samsung%20T5/2026/HAOS/HAOS%20DOCS/HAOS-IIP/docs/notes/foundations/HAOS_IIP_Claim_Gating_Template_v1.md)

Every new release should include:

```text
Claim class:
- A / B / C / D / E / F / G / H

Numerically shown:
- ...

What is measured:
- ...

Baseline:
- ...

Perturbation / stressor:
- ...

Null model / matched control:
- ...

Failure condition:
- ...

Downgrade condition:
- ...

Kill condition:
- ...

Bridge status:
- internal / toy / proxy / real-data-loaded / externally validated

Not claimed:
- ...

Reproduction path:
- ...
```

## 7. Bridge Status Vocabulary

Use these bridge labels exactly:

| Status | Meaning |
| --- | --- |
| internal | Result exists only inside a HAOS-IIP frozen computational regime |
| toy | Synthetic benchmark with known constructed target |
| proxy | Physics-adjacent or biology/materials-adjacent diagnostic, not direct validation |
| real-data-loaded | Public external data were loaded, parsed, and audited |
| externally validated | Independent data or experiment supports a prediction against declared null models |
| failed | Claim did not pass its declared gates |
| unresolved | Data, controls, or reproduction are insufficient |

## 8. Null-Model Requirements

A null model should be close enough to be dangerous.

Weak nulls do not mature a claim.

Preferred nulls:

- matched random graph preserving degree, spectrum, or locality where relevant
- seed variation with identical thresholds
- shuffled labels preserving marginal distributions
- wrong-band or off-target spectral windows
- smoothness-preserving controls
- partial-statistic controls
- baseline physical / biological / materials model when external data are used
- same parser and metrics on a negative-control dataset

Null-model rule:

```text
If the null is easy to beat, it is not enough.
```

## 9. Release Review Checklist

Before a new release is merged:

- [ ] claim class declared
- [ ] standard-language translation present
- [ ] claim box present
- [ ] null model or matched control stated
- [ ] failure condition stated
- [ ] downgrade condition stated
- [ ] kill condition stated
- [ ] reproduction path stated
- [ ] bridge status label used correctly
- [ ] non-claims are explicit
- [ ] public summary does not exceed evidence class

## 10. Authority Boundary

The falsification engine defines how claims fail. It does not decide that a
claim is true.

Passing a failure gate means:

- the claim survived that test

It does not mean:

- the claim survived all possible tests
- the claim became established physics
- the claim no longer needs null models

## 11. Next Work

- convert this Markdown checklist into a machine-readable claim-box template
- add a small validator for new releases that checks for required claim-box fields
- apply the checklist to the next canonical entry paper draft
- backfill claim boxes into the most visible public-facing releases
- build a dashboard table summarizing claim class, bridge status, and current downgrade level
- use the framework comparison matrix to add framework-specific null models
