# A Minimal Structural-Stability Oracle Based on Frozen HAOS-IIP Telemetry

## 1. Problem

The immediate practical problem is not new phase discovery. It is deciding whether an evolving graph-like trajectory remains coherent, becomes weakly deformed, or collapses into fragmentation.

This note defines a small public oracle for that task.

## 2. Frozen Telemetry Basis

The oracle uses only the frozen telemetry layer in:

- `telemetry/frozen_metrics.py`

Core diagnostics are:

- persistence / structural retention
- ordering / temporal consistency
- depth drift / causal deformation
- distance coherence / geometric integrity

The first three are direct frozen telemetry computations. The fourth is a lightweight synthetic proxy used only inside the public demo layer so the same diagnostic ladder can be illustrated on small scenarios.

## 3. Synthetic Scenario Generator

The public demo consumes three baseline synthetic scenarios:

- `baseline`
- `perturbed`
- `fragmented`

It also supports deterministic generator modifiers:

- `--noise`
- `--connectivity-drop`
- `--cluster-split`

These modifiers do not touch any frozen phase bundle. They only perturb the synthetic example trajectories in a reproducible way.

## 4. CLI Interface

Primary entrypoint:

```bash
python3 -m haos_iip.demo stability
```

Machine-readable mode:

```bash
python3 -m haos_iip.demo stability baseline --json
```

Batch scan mode:

```bash
python3 -m haos_iip.demo stability --scan noise=0.00:0.10:0.05 connectivity=0.00:0.20:0.10
```

Generated variant example:

```bash
python3 -m haos_iip.demo stability baseline --noise 0.05 --connectivity-drop 0.2 --cluster-split
```

## 5. Stability Ladder Results

On the default deterministic bundle the classifier yields:

- `baseline -> stable`
- `perturbed -> marginal`
- `fragmented -> unstable`

This demonstrates that the frozen telemetry layer can separate coherent, weakly perturbed, and strongly fragmented regimes in a compact public-facing workflow.

## 6. Limits And Non-Claims

This oracle is a structural-stability tool. It is not a claim about continuum limits, physical spacetime emergence, or ontology.

The demo should be interpreted as:

- a reusable diagnostic layer
- a callable tool surface
- a small parameter-space explorer

It should not be interpreted as a new phase or as a replacement for the frozen authority chain.
