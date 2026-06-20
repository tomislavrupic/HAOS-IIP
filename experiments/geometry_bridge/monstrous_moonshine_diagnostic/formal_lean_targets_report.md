# Formal Lean Target Scaffold Report

## Status

- Result: `OPEN`
- Classification: `FORMAL_TARGET_SCAFFOLD_ONLY`
- Result hash: `formal_lean_targets_ebc40a27ab6bd9db12cc27d5`
- Betti evidence hash: `betti_vector_1d5373fa671fd513fe2d5ac0`

## Labels

- `LEAN_TARGET_SCAFFOLD_BUILT`
- `GRAPH_ISO_COMPONENT_COUNT_TARGET_DECLARED`
- `D4_COMPONENT_COUNT_TARGET_DECLARED`
- `BETTI_VECTOR_INVARIANCE_TARGET_DECLARED`
- `LOCAL_LEAN_PROJECT_NOT_PRESENT`
- `LEAN_CHECK_NOT_RUN`
- `LEAN_THEOREM_NOT_INCLUDED`
- `MOONSHINE_PROOF_NOT_ESTABLISHED`
- `PHYSICAL_BRIDGE_NOT_ESTABLISHED`
- `CLAIM_GATED_FORMAL_TARGET_ONLY`

## Required Local Lean Definitions

| Name | Status | Executable source |
| --- | --- | --- |
| `GaussianPrimeNode` | `MISSING` | `run_betti_component_count.py:base_nodes` |
| `RelationGraph` | `MISSING` | `run_betti_component_count.py:relation_edges` |
| `D4Action` | `MISSING` | `run_betti_component_count.py:transform_nodes` |
| `edgeRelation` | `MISSING` | `run_betti_component_count.py:relation_edges` |
| `componentCount` | `MISSING` | `run_betti_component_count.py:component_count` |
| `bettiOne` | `MISSING` | `run_betti_vector.py:betti_signature` |

## Target Ladder

| ID | Target | Status |
| --- | --- | --- |
| `L1` | `finite_graph_component_count_exists` | `TARGET_ONLY_NOT_LEAN_CHECKED` |
| `L2` | `graph_iso_preserves_component_count` | `TARGET_ONLY_NOT_LEAN_CHECKED` |
| `L3` | `d4_action_induces_graph_iso` | `TARGET_ONLY_NOT_LEAN_CHECKED` |
| `L4` | `d4_preserves_component_count` | `TARGET_ONLY_NOT_LEAN_CHECKED` |
| `L5` | `graph_iso_preserves_betti_one` | `TARGET_ONLY_NOT_LEAN_CHECKED` |
| `L6` | `d4_preserves_betti_vector` | `TARGET_ONLY_NOT_LEAN_CHECKED` |

## Boundary

This is a formal target scaffold, not a Lean-certified theorem.
The checked evidence remains executable Python diagnostics over the declared graph.
The future Lean work should start with graph isomorphism preserving component count before any D4 specialization.
