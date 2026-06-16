# PB-01 Repair Hypotheses

Status: diagnostic draft

Objective

Move PB-01 from `CONTROL_INVALID` toward `CONTROL_VALID` without chasing a
higher HAOS score.

Observed control pattern

- leakage positive control is detected
- shuffled-label control degrades as expected
- topology-destroyed, degree-preserving, weight-shuffled, and parameter-matched
  null controls all fail to move recovery quality enough to satisfy the frozen
  threshold

Repair hypotheses

1. Feature duplication

- The current HAOS proxy may be too close to degree, shortest-path, or
  spectral summaries.
- Next check: compare HAOS feature ranks against degree, spectral gap,
  closeness, and shortest-path baselines on the same holdout split.
- Expected outcome if true: HAOS-only and baseline+HAOS collapse toward the
  same control-insensitive signal.

2. Dynamics too easy

- The synthetic recovery process may leave too much signal after destructive
  perturbations.
- Next check: increase perturbation severity, widen topology destruction, and
  verify that baseline signal drops before any HAOS comparison is rerun.
- Expected outcome if true: control shifts become larger without changing the
  leakage detector.

3. Target leakage through graph family or perturbation magnitude

- The target may still be partially recoverable from family- or severity-level
  metadata.
- Next check: hold graph family fixed, permute perturbation severity labels,
  and compare against parameter-matched nulls that preserve only obvious
  confounds.
- Expected outcome if true: the target score degrades when family/severity
  shortcuts are broken.

4. Control contract is too weak

- The destroyed controls may need a stronger destruction rule than the current
  absolute-shift threshold.
- Next check: define component-specific invalidation criteria for topology,
  degree, spectral, and weight controls before rerunning any score comparison.
- Expected outcome if true: control-invalid classification becomes component
  specific rather than a single coarse shift threshold.

What this note is not

- not a claim that PB-01 is wrong
- not a claim that the HAOS score is meaningless
- not a control repair
- not an external-validation claim

Boundary

This note exists to make the PB-01 repair path testable before PB-02 is run.
