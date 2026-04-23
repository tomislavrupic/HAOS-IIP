# HAOS-IIP Physics Bridge 53.0

This directory is an external post-processing layer over frozen scalar-carrier artifacts.

It does not modify HAOS core, change scalar operators, re-run the scalar simulations, or claim physical ontology. Its job is narrower: translate already-frozen HAOS-IIP scalar measurements into physics-facing proxy observables with explicit `PASS`, `OPEN`, and `WATCH` boundaries.

## Run

```bash
python3 experiments/physics_bridge/physics_observables.py
```

or through the public reproduction entry point:

```bash
python3 examples/quick_reproduce.py --physics-observables-only
```

## Outputs

- `experiments/physics_bridge/results/physics_observables.csv`
- `experiments/physics_bridge/results/physics_observables.json`
- `experiments/physics_bridge/results/physics_observables_summary.md`

## Discipline

- `PASS` means the proxy passes its existing frozen scalar threshold.
- `OPEN` means the proxy fails a threshold and must not be promoted.
- `WATCH` means the aggregate read passes, but a narrower margin remains unresolved.

The bridge dictionary is a test protocol, not a theory feature.

The inverse-square row is intentionally split:

- raw recoverability-gradient inverse-square closure remains `OPEN`;
- shell-native law-aware inverse-square closure is `PASS` under its stated reconstruction thresholds.
