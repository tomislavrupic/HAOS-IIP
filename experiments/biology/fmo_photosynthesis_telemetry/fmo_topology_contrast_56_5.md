# Phase 56.5 FMO Topology Contrast and Lessons

This closes the initial FMO telemetry arc as a diagnostic module rather than forcing a PASS.
The comparison is bounded: microtubule and FMO are toy telemetry substrates, not biological validation.

## Summary Table

| system | runs | pass_rate | recoverability | identity | pathway_identity | active_null_z | max_controls | status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Microtubule 55.2 | 50 | 1.000000 | 0.916789 | 0.871810 | n/a | 28.705018 | 0 | ROBUST PASS |
| FMO 56.1 baseline | 50 | 0.040000 | 0.803724 | 0.984934 | 0.546597 | 2.472024 | 2 | DIAGNOSTIC FAIL |
| FMO 56.2 sink | 50 | 0.080000 | 0.827528 | 0.988291 | 0.622057 | 2.659739 | 2 | DIAGNOSTIC FAIL |
| FMO 56.3 flux | 50 | 0.000000 | 0.930908 | 0.998190 | 0.906155 | 3.227554 | 2 | CONTROL MATCH |
| FMO 56.4 intrinsic | 50 | 0.160000 | 0.848259 | 0.991098 | 0.671540 | 2.612894 | 2 | PRODUCTIVE FALSIFIER |

## Topology Sensitivity

The microtubule lattice and the FMO-like network expose different strengths and limits of the same telemetry instrument.

- Microtubule 55.2 is a robust external toy PASS: structured cylindrical topology, repeated protofilament identity, and multi-scale local coupling remain specific under strict controls.
- FMO preserves site identity strongly, but compact delocalized transfer topology makes pathway identity and control discrimination much harder.
- Explicit flux restoration in 56.3 proves pathway retention can be forced, but it also shows that engineered pathway terms can be copied by controls.
- Intrinsic pathway dynamics in 56.4 reduces the obviousness of the term and slightly improves pass rate, but sacrifices pathway fidelity.

## FMO Lesson

FMO should be kept as a productive falsifier. It does not invalidate the spectral dynamics instrument; it shows that the instrument is topology-sensitive.

The current bottleneck is not simple recoverability or site identity. The bottleneck is branch-specific pathway discrimination on compact transfer graphs.

## Closure Decision

Do not keep tuning FMO until it turns green. Freeze the initial FMO arc as an honest diagnostic result:

- 56.1: baseline falsifier;
- 56.2: sink improves pathway retention but not specificity;
- 56.3: explicit flux solves retention but creates control leakage;
- 56.4: intrinsic directed dynamics improves pass rate modestly but loses retention;
- 56.5: topology sensitivity closure.

Future FMO work should require a new transfer-state model or a branch-specific control-discrimination metric, not just stronger scalar address terms.

## Non-Claims

- No claim is made about molecular FMO dynamics.
- No claim is made about quantum biology or photosynthesis.
- No claim is made that HAOS spectral dynamics universally solve biological telemetry.
- The positive result remains the microtubule toy lattice; the FMO result remains diagnostic.
