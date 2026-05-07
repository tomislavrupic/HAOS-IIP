# HAOS-IIP Claim-Gating Template v1

Status:

- Phase 66 artifact
- mandatory claim-boundary template
- no new result
- no stronger physical claim

Purpose:

Provide a mandatory claim-boundary template for HAOS-IIP papers, dashboards,
Zenodo releases, README updates, and public summaries.

Core rule:

```text
Internal PASS is not external proof.
```

## Mandatory Claim Box

Copy this box into every new release, paper, dashboard, or public summary.

```text
Claim title:
- ...

Claim class:
- A internal numerical / B structural stability / C scaling-refinement /
  D geometry-metric surrogate / E current-closure / F physics bridge /
  G biology telemetry / H materials bridge

Bridge status:
- internal diagnostic only / numerically suggestive / analogy only /
  bridge hypothesis / externally testable candidate / externally supported /
  failed under current tests

Numerically shown inside HAOS-IIP frozen regime:
- ...

Operational meaning:
- ...

Standard-language translation:
- ...

Analogy to known physics / math:
- ...

Bridge hypothesis:
- ...

Required null model:
- ...

Failure condition:
- ...

Downgrade condition:
- ...

Not claimed:
- ...

Reproduction path:
- ...
```

## Bridge Status Labels

Use one primary label.

| Label | Meaning | Allowed public language |
| --- | --- | --- |
| Internal diagnostic only | Result exists inside a declared HAOS-IIP frozen regime | "Inside the frozen regime, this diagnostic passes/fails." |
| Numerically suggestive | Result is positive but controls, scales, or evidence are incomplete | "This suggests a direction for further tests." |
| Analogy only | Result resembles known mathematics or physics but is not a bridge claim | "This is structurally analogous to..." |
| Bridge hypothesis | Result proposes an external observable or dataset target | "This motivates a testable bridge hypothesis." |
| Externally testable candidate | Prediction or metric is specified enough to be tested outside HAOS-IIP | "This can be tested by..." |
| Externally supported | External data or experiment supports the claim against declared nulls | "This external dataset supports the bounded claim..." |
| Failed under current tests | The claim did not pass its declared gates | "Under current tests, this claim fails or is downgraded." |

Do not use "externally supported" unless provenance, null comparison, and
failure conditions are present.

## Public Summary Rule

Every public statement should be reducible to:

```text
We tested X under Y conditions.
It passed/failed Z diagnostic.
This means A internally.
It does not yet mean B externally.
```

Example:

```text
We tested scalar field-like graph response under bounded smooth inhomogeneity.
It passed the declared current-closure diagnostic on the smooth radial branch.
This means the frozen scalar-carrier regime produced a bounded internal closure
signal. It does not yet mean a physical current law has been derived.
```

## Minimal Claim Examples

### Internal Numerical Claim

```text
Numerically shown inside HAOS-IIP frozen regime:
- The declared metric crossed the PASS threshold under the fixed runner.

Operational meaning:
- A diagnostic remained stable under the declared perturbation.

Not claimed:
- External physical validity.
```

### Physics-Bridge Claim

```text
Bridge hypothesis:
- A HAOS-style recoverability metric may map to a measurable external proxy.

Required null model:
- Matched toy/proxy controls and, where available, public external data.

Not claimed:
- Established physics recovery.
```

### Materials-Bridge Claim

```text
Numerically shown:
- Public materials data were loaded, parsed, and audited.

Bridge status:
- real-data-loaded, unless independent prediction against nulls has also passed.

Not claimed:
- The material proves or embodies HAOS-IIP.
```

## Forbidden Upgrades

Do not upgrade:

- internal PASS -> external proof
- toy PASS -> established physics
- proxy PASS -> real-world validation
- analogy -> derivation
- data loaded -> prediction confirmed
- one dataset -> universal claim
- stable metric -> physical ontology

## Release Review Checklist

Before publication or archive:

- [ ] claim box included
- [ ] bridge status label selected
- [ ] standard-language translation appears before HAOS-IIP shorthand
- [ ] operational meaning is stated
- [ ] null model or matched control is stated
- [ ] failure condition is stated
- [ ] downgrade condition is stated
- [ ] not-claimed section is explicit
- [ ] reproduction path is listed
- [ ] public summary passes the four-sentence rule

## Authority Boundary

This template controls claim language. It does not validate claims.

The strongest permitted interpretation is always bounded by the selected bridge
status label and the declared failure gates.

## Next Work

- backfill this claim box into the canonical entry paper draft
- add the template to future numbered paper boilerplate
- add dashboard support for bridge status labels
- create a machine-readable claim-box schema for release validation
