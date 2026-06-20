# Gaussian Prime Norm-Lift Geometry

This sidecar turns the Gaussian integer lattice into a bounded HAOS-IIP-style
arithmetic geometry calibration.

It tests whether:

- Gaussian-prime ramified / inert / split classes are recoverable as shell
  telemetry,
- the norm lift `N(a + bi) = a^2 + b^2` induces stable weighted graph
  Laplacian telemetry,
- cochain-like edge-flow summaries respond when lattice topology is destroyed,
- matched controls move the intended components.

The bundle is synthetic / arithmetic calibration only. It does not claim a
physical bridge, Monster moonshine result, continuum limit, quantum result,
gravity result, or field-theory derivation.

## Run

```bash
uv run python experiments/geometry_bridge/gaussian_prime_norm_lift/run_gaussian_prime_norm_lift.py
```

## Outputs

- `precommitment_contract.json`
- `source_manifest.json`
- `gaussian_prime_nodes.csv`
- `component_scores.csv`
- `control_results.csv`
- `gaussian_prime_norm_lift_report.md`
- `gaussian_prime_norm_lift_result.json`

## Controls

- `known_positive`: self-recovery.
- `rotation_invariance_control`: lattice rotation must preserve telemetry.
- `class_shuffle_control`: prime-shell semantic alignment should degrade.
- `norm_shuffle_control`: norm-lift and weighted spectral telemetry should move.
- `topology_destroyed_control`: spectral and cochain telemetry should move.

## Boundary

The Monster / supersingular-prime material is not included in v0.1. It remains
a future visualization or source-table hook until a trusted, inspectable data
table is added with its own precommitment.

