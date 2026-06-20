# Betti Vector Diagnostic Report

## Status

- Result: `PASS`
- Classification: `TOPOLOGICAL_DIAGNOSTIC_SIDECAR`
- Result hash: `betti_vector_1d5373fa671fd513fe2d5ac0`
- Reference `Betti_0`: `4`
- Reference `Betti_1`: `22`
- Reference nodes: `15`
- Reference edges: `33`
- Reference component sizes: `[12, 1, 1, 1]`

## Labels

- `BETTI_VECTOR_BUILT`
- `BETTI_0_COMPONENT_COUNT_AVAILABLE`
- `BETTI_1_CYCLE_COUNT_AVAILABLE`
- `D4_SYMMETRY_CONTROLS_AVAILABLE`
- `DESTRUCTIVE_CONTROLS_AVAILABLE`
- `LEAN_THEOREM_NOT_INCLUDED`
- `PHYSICAL_BRIDGE_NOT_ESTABLISHED`
- `CLAIM_GATED_GRAPH_TOPOLOGY`
- `BETTI_VECTOR_CONTROLS_PASS`

## Control Results

| Condition | Nodes | Edges | Betti_0 | Betti_1 | Betti vector distance | Edge distance | Status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| known_positive | 15 | 33 | 4 | 22 | 0 | 0 | `PASS` |
| d4_rotate_90 | 15 | 33 | 4 | 22 | 0 | 0 | `PASS` |
| d4_reflect_real | 15 | 33 | 4 | 22 | 0 | 0 | `PASS` |
| isomorphism_relabel | 15 | 33 | 4 | 22 | 0 | 0 | `PASS` |
| threshold_mutation_control | 15 | 18 | 6 | 9 | 15 | 0.454545454545 | `PASS` |
| topology_destroyed_control | 15 | 33 | 2 | 20 | 4 | 0.8 | `PASS` |
| support_replacement_control | 15 | 40 | 3 | 28 | 7 | 0.175 | `PASS` |

## Boundary

`Betti_1` is computed here only as the independent cycle count of the declared finite undirected graph:

```text
Betti_1 = E - V + C
```

This is graph topology, not Moonshine topology. It does not prove Monstrous Moonshine, verify the external Lean SVG claim, or establish any physical bridge.
