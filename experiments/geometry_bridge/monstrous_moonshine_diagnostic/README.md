# Monstrous Moonshine Arithmetic Diagnostic

Lifecycle note: GEO-MM-01 is retained as a supporting arithmetic calibration.
The bounded successor `GEO-MM-02-FORMAL-INVARIANCE` is deferred while the
emergence ladder works at its first unsupported rung. Any later implementation
still requires a new precommitment and remains outside Moonshine or physics
mechanism claims. See the
[branch lifecycle registry](../../../docs/branch_governance/branch_lifecycle_summary.md).

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

## Source Validation

The pinned source and local arithmetic validation pass is available as:

```bash
uv run python experiments/geometry_bridge/monstrous_moonshine_diagnostic/arithmetic_source_validation.py
```

Generated outputs:

- `source_manifest_pinned.json`
- `source_validation_report.md`
- `arithmetic_source_validation_result.json`

Current validation result:

- status: `PASS`
- classification: `PINNED_SOURCE_VALIDATION_ONLY`
- result hash: `source_validation_90b451e478499b9587bb4a5d`
- pinned manifest hash: `source_manifest_pinned_4184b4b53867192d35ed6f01`
- validates supersingular-prime support, Monster-order exponents,
  j-coefficient witnesses, embedded irrep dimensions, and Gaussian-prime
  representatives
- labels include `SOURCE_MANIFEST_PINNED_BUILT`,
  `MOONSHINE_PROOF_NOT_ESTABLISHED`, and
  `PHYSICAL_BRIDGE_NOT_ESTABLISHED`

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

The local Betti vector diagnostic extends the same declared graph to the
graph-native cycle count:

```bash
uv run python experiments/geometry_bridge/monstrous_moonshine_diagnostic/run_betti_vector.py
```

Generated outputs:

- `betti_vector_control_results.csv`
- `betti_vector_report.md`
- `betti_vector_result.json`

Current reference signature:

- `Betti_0`: `4`
- `Betti_1`: `22`
- nodes: `15`
- edges: `33`
- result hash: `betti_vector_1d5373fa671fd513fe2d5ac0`
- labels include `BETTI_VECTOR_CONTROLS_PASS`, `LEAN_THEOREM_NOT_INCLUDED`,
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

## Moonshine-Betti Shared-Support Bridge

The bounded bridge manifest is available as:

```bash
uv run python experiments/geometry_bridge/monstrous_moonshine_diagnostic/run_moonshine_betti_bridge.py
```

Generated outputs:

- `moonshine_betti_bridge_manifest.json`
- `moonshine_betti_bridge_covariance.csv`
- `moonshine_betti_bridge_report.md`
- `moonshine_betti_bridge_result.json`

Current bridge result:

- bridge type: `shared-support diagnostic coupling`
- useful phrase: `shared-support perturbation covariance`
- dangerous phrase: `Moonshine-Betti mechanism`
- shared-support replacement moves both channels:
  - Moonshine total distance: `0.522389629805`
  - Betti vector distance: `7`
  - Betti edge distance: `0.175`
- result hash: `moonshine_betti_bridge_a7d2b62c0c8be6d8c6b757e9`
- labels include `SHARED_SUPPORT_COVARIANCE_PASS`,
  `MOONSHINE_BETTI_MECHANISM_NOT_ESTABLISHED`, and
  `PHYSICAL_BRIDGE_NOT_ESTABLISHED`

## Formal Lean Targets

The local formal target scaffold is available as:

```bash
uv run python experiments/geometry_bridge/monstrous_moonshine_diagnostic/run_formal_lean_targets.py
```

Generated outputs:

- `formal_lean_targets_manifest.json`
- `formal_lean_targets_report.md`
- `lean_graph_invariance_targets.lean`

Current formal target result:

- status: `OPEN`
- classification: `FORMAL_TARGET_SCAFFOLD_ONLY`
- local Lean project: `not present`
- Lean check: `not run`
- first target after definitions: `graph_iso_preserves_component_count`
- D4 target: `d4_preserves_component_count`
- Betti-vector target: `d4_preserves_betti_vector`
- result hash: `formal_lean_targets_ebc40a27ab6bd9db12cc27d5`
- labels include `LEAN_TARGET_SCAFFOLD_BUILT`, `LEAN_THEOREM_NOT_INCLUDED`,
  and `PHYSICAL_BRIDGE_NOT_ESTABLISHED`

## Failure Ledger and Composite Dashboard

Failure conditions are recorded in:

- `failure_conditions.md`

The composite recoverability dashboard is available as:

```bash
uv run python experiments/geometry_bridge/monstrous_moonshine_diagnostic/run_geometry_bridge_recoverability_dashboard.py
```

Generated outputs:

- `geometry_bridge_recoverability_report.md`
- `geometry_bridge_recoverability_result.json`

Current dashboard result:

- status: `PASS_WITH_FRAGILITY`
- classification: `COMPOSITE_DIAGNOSTIC_DASHBOARD_ONLY`
- result hash: `geometry_bridge_dashboard_dc14cbaf7f66529e66a57f49`
- open / partial channels remain visible:
  - Gaussian-prime norm-lift support: `PARTIAL`
  - null ensemble rarity: `OPEN`
  - formal Lean targets: `OPEN`
- best-case classification remains:
  `LOCAL_SIGNAL_STABLE`, `MATHEMATICALLY_INTERESTING`,
  `ENGINEERING_FEASIBLE`, `FORMALLY_PARTIAL`, `GLOBAL_CLAIM_OPEN`

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
the local Betti runners.

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
