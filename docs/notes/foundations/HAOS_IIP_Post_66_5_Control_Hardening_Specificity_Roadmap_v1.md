# HAOS-IIP Post-66.5 Control Hardening and Specificity Roadmap v1

This document is a Continuum Physics Scaling Roadmap artifact, not a new phase.

Status: next-cycle roadmap artifact.

Layer: CP ladder / scale-bridge legitimacy / Phase 66 audit discipline.

Phase 67 remains parked.

## 1. Current Verdict

The post-66.5 roadmap audit executed cleanly, but the major gates remain `OPEN`.

This is the correct outcome.

The repo is stronger because the audit located the exact failure surfaces instead of upgrading language prematurely.

## 2. Next Cycle Target

Next cycle target:

```text
control hardening + specificity
```

The next cycle should not add new phase concepts.

It should harden the controls and improve diagnostic specificity inside the existing CP ladder.

Core rule:

```text
next cycle is not expansion
next cycle is specificity
```

## 3. Priority 1 - CP2b: Control-Hardened Same-Surrogate Recovery

Problem:

CP2 showed strong admissible recovery, but matched controls also recovered too strongly.

Question:

```text
Is the recovery branch-specific, or is the current round-trip test too easy?
```

Required controls:

- eigenmode-order shuffled control
- matched-spectrum / broken-eigenvector control
- random orthogonal basis control
- kernel-weight permutation control
- degree-preserving graph rewiring
- phase-shuffled surrogate control

Success condition:

Admissible recovery stays high while hardened controls fall clearly below the declared threshold or show much worse round-trip error.

Failure condition:

Hardened controls recover at admissible levels.

Claim boundary:

Until CP2b separates admissible recovery from hardened controls, same-surrogate recovery remains useful but not scale-bridge closure.

## 4. Priority 2 - CP3b: Effective-Law Specificity Audit

Problem:

Some CP3 rows pass, but rescaled-invariant flow does not separate.

Question:

```text
Which effective-law candidates are genuinely branch-specific?
```

Required work:

- isolate coefficient-flow rows that pass
- isolate shell-slope rows that pass
- isolate propagation-speed rows that pass
- treat rescaled-invariant flow as `OPEN` / failed separator
- compare each candidate law against matched controls separately

Success condition:

At least one compact law family fits admissible branches and clearly outperforms controls under refinement.

Failure condition:

Every candidate law requires retuning or appears in controls.

Claim boundary:

Passing CP3b may support bounded effective-law behavior inside the tested regime. It does not derive physical field equations.

## 5. Priority 3 - Comparative Diagnostic Specificity

Problem:

The Kuramoto sidecar detected early-warning behavior but failed edge-signature specificity.

Question:

```text
Does HAOS recoverability add anything beyond generic early-warning detection?
```

Required work:

- separate early-warning detection from edge-signature specificity
- add null oscillator controls
- add shuffled coupling controls
- add graph Laplacian baseline
- report where HAOS is useful and where it is not

Success condition:

HAOS-IIP diagnostics discriminate structural degradation better than at least one standard baseline under defined conditions.

Failure condition:

HAOS-IIP detects only generic early-warning behavior.

Claim boundary:

Comparative diagnostic only. No superiority claim unless measured against declared controls and baselines.

## 6. Priority 4 - CP5b: Minimal Universality Grid

Problem:

CP5 has seed/refinement evidence, but kernel-width and substrate-family coverage remain insufficient.

Question:

```text
Does the observed behavior survive controlled variation?
```

Minimum grid:

- 2-3 seeds
- 2 kernel widths
- 2 substrate families
- one mild disorder condition
- matched controls for each

Success condition:

The same diagnostic pattern survives across the grid without retuning.

Failure condition:

Each substrate/kernel change rewrites the result.

Claim boundary:

Only after CP5b begins to pass should stronger continuum-feasibility language be considered, and even then only inside the tested operator regime.

## 7. Cycle Rules

- no score upgrade until at least one `OPEN` gate becomes `PASS` under hardened controls
- no Phase 67 activation
- no stronger continuum language
- no placeholder numbers
- no new phase concepts
- all scale-bridge work must include matched controls
- all continuum-facing language must use the Scale-Bridge Status Box

## 8. Immediate Start

Start with:

```text
CP2b control-hardened same-surrogate recovery
```

Recommended artifact:

```text
HAOS_IIP_CP2b_Control_Hardened_Same_Surrogate_Recovery_v1.md
```

Best placement:

```text
docs/notes/foundations/
```

If numerical code is added, keep it in the existing continuum / scale-bridge / simulation structure rather than creating Phase 67.

## 9. Final Boundary

Scale-bridge evidence is useful.

Scale-bridge evidence is not continuum proof.
