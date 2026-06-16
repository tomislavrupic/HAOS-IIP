# Same-Surrogate Control Integrity Report

Status: control-hardening diagnostic, not a continuum claim.

Frozen thresholds: admissible recovery >= 95.0%, matched-control recovery < 15.0%, slope error <= 0.05, causal-depth drift <= 0.10, triangle violation <= 0.05.

The available CP2 evidence contains two frozen surrogate rows:

| Surrogate | Admissible recovery % | Matched-control recovery % | Admissible error | Control error | Status |
| --- | ---: | ---: | ---: | ---: | --- |
| Phase XVIII distance surrogate | 100.000000 | 50.000000 | 0.004717 | 0.231481 | FAIL |
| Phase X low-mode spectral projection | 100.000000 | 100.000000 | 0.037338 | 0.034210 | FAIL |

## Diagnostic Readout

- Branch recovery distribution: [100.0, 100.0]
- Control recovery distribution: [50.0, 100.0]
- Branch mean recovery: 100.000000%
- Control mean recovery: 75.000000%
- Effect size: 25.000000 percentage points
- Overlap coefficient: 0.500000
- Branch range: [100.000000, 100.000000]
- Control range: [50.000000, 100.000000]

## Control Integrity

- The matched control does not fall below the frozen strict threshold.
- One control row equals admissible recovery, so branch/control separation is not established.
- No additional hardening control family is available in the frozen same-surrogate artifact slice.

## Terminal Classification

- `CP2_CONTROL_INVALID`
- `CP2_SAME_SURROGATE_CLOSED`
