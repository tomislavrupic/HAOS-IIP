# Same-Surrogate Coarse-Graining Recovery Summary

Status: Continuum Physics Scaling Roadmap artifact, not a new phase.

Thresholds declared before interpretation:

- admissible recovery >= 95.0%
- matched-control recovery < 15.0%
- slope round-trip relative error <= 0.05
- causal-depth drift <= 0.10
- triangle violation <= 0.05

| Surrogate type | Round trip | Admissible recovery % | Matched-control recovery % | Admissible error | Control error | Status |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Phase XVIII distance surrogate | 84 -> 72 -> 84 | 100.000000 | 50.000000 | 0.004717 | 0.231481 | FAIL |
| Phase X low-mode spectral projection | R5 -> projected low-mode subspace -> R5 readout | 100.000000 | 100.000000 | 0.037338 | 0.034210 | FAIL |

Interpretation:

- Phase XVIII distance surrogate: admissible branch recovery is complete under declared checks, but matched-control recovery is not below the strict threshold, so the CP2 gate is `CP2_CONTROL_INVALID`.
- Phase X low-mode spectral projection: branch and control both recover, so the row supports projection bookkeeping but not branch-specific control separation.

Claim boundary: this is same-surrogate recovery bookkeeping inside frozen ledgers. It is not a continuum limit, physical metric, or physical correspondence claim.
