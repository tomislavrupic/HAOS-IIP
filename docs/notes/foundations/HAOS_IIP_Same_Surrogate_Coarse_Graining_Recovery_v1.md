# HAOS-IIP Same-Surrogate Coarse-Graining Recovery v1

This document is a Continuum Physics Scaling Roadmap artifact, not a new phase.

Status: CP2 audit run, `OPEN`.

Phase 67 remains parked.

## 1. Purpose

This audit executes the first same-surrogate coarse-graining recovery check after Release 66.5.

Core question:

```text
Can the same spectral address / metric surrogate survive a refinement round trip?

fine -> coarse -> fine

And does it survive better than matched controls?
```

Core rule:

```text
A surrogate that cannot survive coarse-graining is not a scale bridge.
```

This run uses frozen ledgers only. It introduces no new dynamics, no new phase number, and no physical-continuum claim.

## 2. Source Artifacts

Executable extractor:

- `continuum-sketch/build_same_surrogate_coarse_graining_recovery.py`

Generated outputs:

- `continuum-sketch/same_surrogate_coarse_graining_recovery.csv`
- `continuum-sketch/same_surrogate_coarse_graining_recovery_summary.md`

Frozen source ledgers:

- `phase18-distance-surrogate/runs/phase18_refinement_scaling.csv`
- `phase18-distance-surrogate/runs/phase18_triangle_violation_rate.csv`
- `phase10-bridge/runs/phase10_coarse_grain_summary.csv`

## 3. Predeclared Thresholds

Thresholds were declared before interpreting the generated results:

- admissible branch recovery >= `95.0%`
- matched-control recovery < `15.0%`
- slope round-trip relative error <= `0.05`
- causal-depth drift <= `0.10`
- triangle violation <= `0.05`
- no physical-continuum claim

## 4. CP2 Recovery Table

| Surrogate type | Refinement round trip | Recovery % admissible | Recovery % matched control | Round-trip error admissible | Round-trip error control | Spectral-address persistence | Metric-like check | Status | Claim boundary |
| --- | --- | ---: | ---: | ---: | ---: | --- | --- | --- | --- |
| Phase XVIII distance surrogate | `84 -> 72 -> 84` | 100.000000 | 50.000000 | 0.004717 | 0.231481 | not applicable to metric-surrogate row | PASS | `OPEN` | metric-like surrogate round-trip diagnostic only; no continuum proof |
| Phase X low-mode spectral projection | `R5 -> projected low-mode subspace -> R5 readout` | 100.000000 | 100.000000 | 0.037338 | 0.034210 | 100.000000 | PASS | `OPEN` | spectral projection recovery bookkeeping only; no branch-specific scale-bridge proof |

## 5. Interpretation

The Phase XVIII distance surrogate recovers cleanly on the admissible branch under the declared checks:

- branch recovery: `100.000000%`
- matched-control recovery: `50.000000%`
- branch round-trip error: `0.004717`
- matched-control round-trip error: `0.231481`
- triangle violation: `0.000000`

This is encouraging, but it does not close CP2 because the matched control still recovers `50.000000%`, above the strict `< 15%` threshold.

The Phase X low-mode spectral projection also recovers cleanly:

- branch recovery: `100.000000%`
- matched-control recovery: `100.000000%`
- branch trace-window error: `0.037338`
- matched-control trace-window error: `0.034210`
- spectral-address persistence proxy: `100.000000%`

This supports projection-recovery bookkeeping, but it does not establish branch-specific scale-bridge separation because the matched control also recovers.

## 6. CP2 Status

```text
CP2 same-surrogate coarse-graining recovery: OPEN
```

Reason:

- admissible branch recovery passes
- same surrogate is used before and after projection
- metric-like checks pass where applicable
- matched controls do not fail strongly enough under the predeclared recovery threshold

## 7. Claim Boundary

This audit may show:

```text
A spectral / metric surrogate can preserve identity under refinement round trip inside frozen ledgers.
```

It does not show:

- proven continuum limit
- derived spacetime
- physical metric
- real curvature
- physical correspondence
- replacement of standard physics

## 8. Next Step

The next CP2 pass should tighten matched-control separation rather than changing the surrogate definition.

Recommended next move:

```text
Run the same CP2 recovery extractor on additional substrate / kernel / seed families, or add a control-specific recovery channel that is predeclared before interpretation.
```

No Phase 67 activation.
