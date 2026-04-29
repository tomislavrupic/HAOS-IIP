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

Any root-note or harmonic-partial reading must remain external to this directory and cannot change the bridge rows.

The inverse-square row is intentionally split:

- raw recoverability-gradient inverse-square closure remains `OPEN`;
- shell-native law-aware inverse-square closure is `PASS` under its stated reconstruction thresholds.

The power-law scaling row is intentionally split:

- raw local-gradient power scaling remains `OPEN`;
- shell-native law-aware power scaling is `PASS` as a bounded continuum-like proxy.

The localized-bump row is also intentionally split:

- weak localized bump response is `PASS` after excluding the earliest source-core transient layer;
- stronger localized bump response remains `OPEN` under the same transport thresholds.

The disorder-native row is also intentionally split:

- aggregated mild-disorder closure is `PASS`;
- seed-universal disorder closure is separately tracked and is currently `PASS` on the tested delayed fit window `[0.006, 0.024]`.

## Optional Sidecar Audit

An additional external sidecar audit is available for the raw local-gradient boundary:

```bash
python3 experiments/physics_bridge/raw_gradient_shape_audit.py
```

or through the public reproduction entry point:

```bash
python3 examples/quick_reproduce.py --raw-gradient-audit-only
```

This sidecar audit does not change the canonical bridge rows. Its strongest bounded result is narrower: the full raw local-gradient observable remains `OPEN`, while a source-core-excluded, amplitude-normalized shape-only read can pass the same raw thresholds as a sidecar proxy.

## Future Celestial Boundary Audit

A queued post-Phase-56 roadmap is saved at:

`docs/notes/foundations/HAOS_IIP_Celestial_Holography_Audit_Roadmap_v1.md`

That roadmap treats celestial holography as a high-standard boundary audit for
HAOS-IIP physics-facing language. It is not a core change and not a claim of
BMS, Virasoro, S-matrix, soft-theorem, or gravitational-memory recovery.
