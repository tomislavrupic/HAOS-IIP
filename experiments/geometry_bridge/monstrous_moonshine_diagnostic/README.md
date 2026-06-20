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

## Betti_0 Sidecar

The local Betti_0 / component-count diagnostic is available as:

```bash
uv run python experiments/geometry_bridge/monstrous_moonshine_diagnostic/run_betti_component_count.py
```

Generated outputs:

- `betti_relation_graph_spec.json`
- `betti_relation_graph_nodes.csv`
- `betti_relation_graph_edges.csv`
- `betti_control_results.csv`
- `betti_component_count_report.md`
- `betti_component_count_result.json`

Current reference result:

- `Betti_0`: `4`
- relation edges: `33`
- result hash: `betti_component_1454dd195fe727984643d0dc`
- labels include `BETTI_CONTROLS_PASS`, `LEAN_THEOREM_NOT_INCLUDED`,
  and `PHYSICAL_BRIDGE_NOT_ESTABLISHED`

Threshold sensitivity is checked by:

```bash
uv run python experiments/geometry_bridge/monstrous_moonshine_diagnostic/run_betti_threshold_sweep.py
uv run python experiments/geometry_bridge/monstrous_moonshine_diagnostic/run_betti_null_ensemble.py
```

Generated outputs:

- `betti_threshold_sweep.csv`
- `betti_threshold_sweep_report.md`
- `betti_null_ensemble_results.csv`
- `betti_null_ensemble_report.md`
- `betti_threshold_sweep_result.json`

Current sweep result:

- `Betti_0` stability band: `[8.0, 11.5]`
- exact edge-signature band: `[8.0, 8.5]`
- edge-neighborhood band: `[8.0, 10.0]`
- null false-positive rate: `0.110000`
- result hash: `betti_sweep_6d8f185a0d3a8ad5130832c4`

## Visual Sidecar

An external animated SVG visualization is stored under:

- `visuals/anim_supersingular_betti.svg`
- `visuals/visual_source_manifest.json`

The HAOS-IIP interpretation is documented in:

- `betti_diagram_haos_connection.md`

This visual sidecar is review material only. Its source text references a
Lean-certified Betti diagram, but this repository does not independently check
that Lean proof unless the proof source is added and run locally. The executable
authority for this folder remains `run_monstrous_moonshine_diagnostic.py` and
`run_betti_component_count.py`.

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
