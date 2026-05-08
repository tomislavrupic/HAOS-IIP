# HAOS-IIP Next Roadmap - Post-66.5

Status: future roadmap artifact.

Layer: foundations / scale-bridge audit continuation.

This is not Phase 67, not a new mechanism, and not a physical correspondence claim.

## 1. Current State

HAOS-IIP is in strong methodological condition after Release 66.5.

Phase 66 has established:

- canonical entry paper
- hostile-audit layer
- translation table
- falsification engine
- framework comparison matrix
- claim-gating template
- scale-bridge legitimacy audit
- paper spine index
- `uv`-locked reproducibility

Current strongest unresolved bottleneck:

```text
scale-bridge legitimacy
```

Internal audit score:

```text
8.9 / 10
```

This does not mean the repo is weak. It means the last hard bridge gates are not fully closed yet.

## 2. Next Target

Same-surrogate coarse-graining recovery.

Core question:

```text
Can the same spectral address / metric surrogate survive a refinement round trip?

fine -> coarse -> fine

And does it survive better than matched controls?
```

## 3. Core Rule

A surrogate that cannot survive coarse-graining is not a scale bridge.

Rules:

- do not change the surrogate definition mid-test
- do not use placeholder numbers
- do not upgrade language beyond generated evidence

## 4. Immediate Artifact

Create one of these:

```text
HAOS_IIP_Same_Surrogate_Coarse_Graining_Recovery_v1.md
```

or:

```text
HAOS_IIP_Coarse_Graining_Recovery_Audit_66_6_v1.md
```

Best placement:

```text
docs/notes/foundations/
```

If numerical code is added, keep it under the existing continuum / scale-bridge / simulation structure rather than creating Phase 67.

## 5. Required Measurements

The audit should produce:

1. admissible branch recovery %
2. matched-control recovery %
3. round-trip error
4. spectral-address persistence
5. metric-surrogate ordering preservation
6. triangle-inequality or bounded-violation check, where applicable
7. `PASS` / `OPEN` / `FAIL` status
8. explicit claim boundary

## 6. Success Threshold

Declare the threshold before interpreting results.

Suggested starting threshold:

- admissible recovery >= 95%
- matched control < 15%
- same surrogate used before and after projection
- metric-like properties preserved within declared tolerance
- no physical-continuum claim

## 7. Coarse-Graining Recovery Table Schema

| Surrogate type | Refinement round trip | Recovery % admissible | Recovery % matched control | Round-trip error | Metric-like check | Status | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Spectral address separation | fine -> coarse -> fine | generated only | generated only | generated only | generated only | `PASS` / `OPEN` / `FAIL` | no placeholders |
| Local metric surrogate | fine -> coarse -> fine | generated only | generated only | generated only | generated only | `PASS` / `OPEN` / `FAIL` | no placeholders |
| Scalar-carrier norm | fine -> coarse -> fine | generated only | generated only | generated only | generated only | `PASS` / `OPEN` / `FAIL` | no placeholders |
| Phase XVIII distance surrogate | fine -> coarse -> fine | generated only | generated only | generated only | generated only | `PASS` / `OPEN` / `FAIL` | no placeholders |

## 8. Claim Boundary

This audit may show:

```text
A spectral / metric surrogate preserves identity under refinement round trip better than matched controls.
```

It does not show:

- proven continuum limit
- derived spacetime
- physical metric
- real curvature
- physical correspondence
- replacement of standard physics

## 9. After This Gate

If same-surrogate coarse-graining recovery passes cleanly:

1. update the Scale-Bridge Status and Limits paragraph
2. recompile the next PDF release, likely `66.6`
3. update the hostile-audit score from `8.9` toward `9.x` only as an internal audit score
4. then consider the bounded predictive bridge experiment

## 10. Next After Coarse-Graining

Bounded predictive bridge experiment.

Candidate targets:

- ferron coherence sidecar
- alpha-Fe2O3 magnon sidecar
- scalar-carrier physics-bridge observable

Core requirement:

```text
Make a prediction before the full audit.
Then compare against data or held-out structure.
```

## 11. Still Parked

Phase 67 remains parked.

No delayed-memory / interaction-memory expansion until the scale-bridge gate is cleaner.

Current priority:

```text
finish the scale bridge before expanding the system
```
