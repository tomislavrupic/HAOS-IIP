# Monstrous Moonshine Arithmetic Diagnostic

This sidecar records a small, bounded Monstrous Moonshine support diagnostic for
the HAOS-IIP geometry bridge.

It includes:

- the 15 supersingular primes,
- the Monster group order prime support and exponents,
- a small set of low j-coefficient decomposition witnesses,
- factor support for a small embedded subset of Monster irreducible dimensions,
- Gaussian-prime ramified / inert / split classes over the same prime support.

It does **not** prove Monstrous Moonshine, construct the Moonshine module,
include the complete 194-dimensional irrep catalog, or establish a physical
bridge.

## Run

```bash
uv run python experiments/geometry_bridge/monstrous_moonshine_diagnostic/run_monstrous_moonshine_diagnostic.py
```

## Outputs

- `precommitment_contract.json`
- `source_manifest.json`
- `supersingular_prime_table.csv`
- `moonshine_witnesses.csv`
- `component_scores.csv`
- `control_results.csv`
- `monstrous_moonshine_diagnostic_report.md`
- `monstrous_moonshine_diagnostic_result.json`

## Visual Sidecar

An external animated SVG visualization is stored under:

- `visuals/anim_supersingular_betti.svg`
- `visuals/visual_source_manifest.json`

The HAOS-IIP interpretation is documented in:

- `betti_diagram_haos_connection.md`

This visual sidecar is review material only. Its source text references a
Lean-certified Betti diagram, but this repository does not independently check
that Lean proof unless the proof source is added and run locally. The executable
authority for this folder remains `run_monstrous_moonshine_diagnostic.py`.

## Controls

- `known_positive`: self-recovery.
- `nonsupersingular_prime_control`: replaces `71` with `73`.
- `exponent_shuffle_control`: keeps support but moves Monster-order exponents.
- `decomposition_break_control`: breaks one j-coefficient witness.
- `gaussian_class_swap_control`: swaps one inert/split class relation.

## Boundary

This is an arithmetic support diagnostic only. The Gaussian-prime norm-lift
bundle remains a narrative and structural neighbor, not a derivation of
moonshine.
