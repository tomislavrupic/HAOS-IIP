# HAOS-AF-AR-01: Bounded Adaptive Recoverability

This sidecar tests whether HAOS-IIP's perturbation and recoverability discipline
can support a falsifiable operational criterion for antifragility.

The tested system is a finite weighted ring with local distance-two backup
interactions. Exposure damages a contiguous block of primary interactions. The
target sees only the maximum fractional interaction loss and persistently
reallocates a frozen fraction of the same total interaction-weight budget from
primary edges to uniformly distributed local backup edges. It receives no
held-out location, family, score, future state, or full-state checkpoint.

The later test is held out by location and includes three families:

- displaced primary-block damage;
- distributed primary damage;
- backup-channel damage.

The invariant audit contract compares post-exposure recovery capacity with the
unadapted system under the identical later perturbation. Recovery capacity is
the harmonic mean of normalized algebraic connectivity and normalized global
effective conductance. Identity, primary-capacity retention, total budget,
tail behavior, collapse ordering, dose response, and matched controls are
separate gates.

Controls are passive no-update, equal-budget random reallocation, reversed
update, and sham no-stress. A statically hardened graph is reported only as a
noncausal ceiling reference: it shows that the adapted endpoint can be designed
in advance and therefore cannot establish a unique or optimal adaptation law.

## Reproduce

Run tests and the implementation-only smoke partition:

```bash
uv run python -m unittest experiments.antifragility.adaptive_recoverability_v1.tests.test_adaptive_recoverability
uv run python experiments/antifragility/adaptive_recoverability_v1/run_experiment.py \
  --partition smoke --output-root /tmp/haos-af-ar-01-smoke
```

The committed final bundle is validated with:

```bash
uv run python experiments/antifragility/adaptive_recoverability_v1/check_bundle.py
```

## Claim Boundary

The maximum positive classification is `BOUNDED_TOY_ANTIFRAGILITY_SUPPORTED`.
It means only that the frozen stress-triggered update improved later held-out
recoverability on this finite graph schedule while passing the declared
identity, budget, tail, collapse, and control gates.

It does not establish general antifragility, a physical law, biological
adaptation, economic advantage, autonomy, agency, universality, a continuum
limit, new physics, or a higher emergence rung.
