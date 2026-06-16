# Synthetic Hidden Geometry Recovery

This benchmark tests whether a frozen observer can recover:

- intrinsic distance
- orientation / handedness
- transformation class
- held-out relations

from a hidden synthetic geometry.

Current state: `BENCHMARK_OPEN`.

Synthetic distance, orientation, and held-out relation recovery pass. The
residual open criterion is transformation recovery: the frozen observer does
not yet recover the hidden transformation class robustly on holdout.

It is synthetic calibration only and does not authorize Bell or physical
mechanism claims.
