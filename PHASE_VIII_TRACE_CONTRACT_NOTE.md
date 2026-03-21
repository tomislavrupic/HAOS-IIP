# Phase VIII Spectral-Trace Contract Note

Phase VIII freezes the spectral-trace baseline for all later short-time trace work on the current branch.

The following objects are now part of the frozen spectral-trace contract:

- the operational short-time window `[0.02, 0.03, 0.05, 0.075, 0.1]`
- the master deterministic probe grid `[0.01, 0.015, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0]`
- the exact trace proxy `T_h(t) = Tr(exp(-t * delta_h))`
- the recorded truncation checkpoints `full_exact`, `lambda_le_7_5`, and `lambda_le_7_8`
- the window-scoped slope, kernel-sensitivity, and coefficient-like diagnostics frozen in `phase8-trace/phase8_trace_manifest.json`

Future phases may use these objects, extend them explicitly, or supersede them only through a new recorded authority step.

Future phases may not silently:

- redefine the operational short-time regime
- change the probe grid while retaining the same contract language
- reinterpret later diagnostics as if they were produced on a different trace baseline

Any later phase that requires a different probe grid or a different operational window must declare that change explicitly and record the superseding authority artifact.
