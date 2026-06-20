#!/usr/bin/env python3
"""Betti vector diagnostic for the moonshine sidecar.

This script extends the local Betti_0 component-count diagnostic with the
graph-native cycle count Betti_1 = E - V + C. It reuses the declared relation
graph from run_betti_component_count.py and makes no Moonshine, Lean, or
physical-bridge claim.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from run_betti_component_count import (  # noqa: E402
    BettiComponentConfig,
    CONDITIONS,
    canonical_edge_signature,
    component_count,
    condition_edges,
    condition_nodes,
    expected_response,
    jaccard_distance,
    load_spec,
    relation_edges,
    stable_hash,
    write_csv,
    write_json,
)


CONTROL_RESULTS_PATH = ROOT / "betti_vector_control_results.csv"
REPORT_PATH = ROOT / "betti_vector_report.md"
RESULT_PATH = ROOT / "betti_vector_result.json"

CONTROL_FIELDS = [
    "condition",
    "expected_response",
    "nodes",
    "edges",
    "betti0",
    "betti1",
    "betti_vector_distance",
    "edge_jaccard_distance",
    "observed_status",
]


@dataclass(frozen=True)
class BettiVectorConfig:
    version: str = "betti-vector-v0.1"
    threshold: float = 8.0
    mutated_threshold: float = 4.0
    invariance_max: float = 1.0e-12
    destructive_betti_vector_min: float = 1.0
    destructive_edge_min: float = 0.10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Betti_0/Betti_1 vector diagnostic.")
    parser.add_argument("--output-dir", type=Path, default=ROOT)
    return parser.parse_args()


def component_sizes(node_count: int, edges: list[tuple[int, int, float]]) -> list[int]:
    neighbors: list[list[int]] = [[] for _ in range(node_count)]
    for left, right, _ in edges:
        neighbors[left].append(right)
        neighbors[right].append(left)
    seen = [False] * node_count
    sizes: list[int] = []
    for start in range(node_count):
        if seen[start]:
            continue
        stack = [start]
        seen[start] = True
        size = 0
        while stack:
            node = stack.pop()
            size += 1
            for neighbor in neighbors[node]:
                if not seen[neighbor]:
                    seen[neighbor] = True
                    stack.append(neighbor)
        sizes.append(size)
    return sorted(sizes, reverse=True)


def betti_signature(nodes: list[dict[str, Any]], edges: list[tuple[int, int, float]]) -> dict[str, Any]:
    node_count = len(nodes)
    edge_count = len(edges)
    betti0 = component_count(node_count, edges)
    betti1 = edge_count - node_count + betti0
    return {
        "Betti_0": betti0,
        "Betti_1": betti1,
        "nodes": node_count,
        "edges": edge_count,
        "component_size_distribution": component_sizes(node_count, edges),
    }


def betti_vector_distance(reference: dict[str, Any], observed: dict[str, Any]) -> int:
    return abs(int(reference["Betti_0"]) - int(observed["Betti_0"])) + abs(
        int(reference["Betti_1"]) - int(observed["Betti_1"])
    )


def observed_status(
    condition: str,
    vector_distance: float,
    edge_distance: float,
    config: BettiVectorConfig,
) -> str:
    if condition in {"known_positive", "d4_rotate_90", "d4_reflect_real", "isomorphism_relabel"}:
        if vector_distance <= config.invariance_max and edge_distance <= config.invariance_max:
            return "PASS"
        return "FAIL_SYMMETRY_STABILITY"
    if vector_distance >= config.destructive_betti_vector_min or edge_distance >= config.destructive_edge_min:
        return "PASS"
    return "FAIL_DESTRUCTIVE_CONTROL"


def as_component_config(config: BettiVectorConfig) -> BettiComponentConfig:
    return BettiComponentConfig(
        threshold=config.threshold,
        mutated_threshold=config.mutated_threshold,
        invariance_max=config.invariance_max,
        destructive_edge_min=config.destructive_edge_min,
    )


def run_betti_vector(config: BettiVectorConfig, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    component_config = as_component_config(config)
    reference_nodes = condition_nodes("known_positive")
    reference_edges = relation_edges(reference_nodes, config.threshold)
    reference_signature = canonical_edge_signature(reference_nodes, reference_edges)
    reference_betti = betti_signature(reference_nodes, reference_edges)

    control_rows: list[dict[str, Any]] = []
    condition_results: dict[str, Any] = {}
    for condition in CONDITIONS:
        nodes = condition_nodes(condition)
        edges = condition_edges(condition, nodes, component_config, len(reference_edges))
        observed_betti = betti_signature(nodes, edges)
        edge_signature = canonical_edge_signature(nodes, edges)
        vector_distance = betti_vector_distance(reference_betti, observed_betti)
        edge_distance = jaccard_distance(reference_signature, edge_signature)
        status = observed_status(condition, vector_distance, edge_distance, config)
        row = {
            "condition": condition,
            "expected_response": expected_response(condition),
            "nodes": observed_betti["nodes"],
            "edges": observed_betti["edges"],
            "betti0": observed_betti["Betti_0"],
            "betti1": observed_betti["Betti_1"],
            "betti_vector_distance": f"{vector_distance:.12g}",
            "edge_jaccard_distance": f"{edge_distance:.12g}",
            "observed_status": status,
        }
        control_rows.append(row)
        condition_results[condition] = {
            "signature": observed_betti,
            "expected_response": expected_response(condition),
            "betti_vector_distance": vector_distance,
            "edge_jaccard_distance": edge_distance,
            "observed_status": status,
        }

    controls_pass = all(result["observed_status"] == "PASS" for result in condition_results.values())
    labels = [
        "BETTI_VECTOR_BUILT",
        "BETTI_0_COMPONENT_COUNT_AVAILABLE",
        "BETTI_1_CYCLE_COUNT_AVAILABLE",
        "D4_SYMMETRY_CONTROLS_AVAILABLE",
        "DESTRUCTIVE_CONTROLS_AVAILABLE",
        "LEAN_THEOREM_NOT_INCLUDED",
        "PHYSICAL_BRIDGE_NOT_ESTABLISHED",
        "CLAIM_GATED_GRAPH_TOPOLOGY",
    ]
    labels.append("BETTI_VECTOR_CONTROLS_PASS" if controls_pass else "BETTI_VECTOR_CONTROLS_PARTIAL")

    result = {
        "bridge_id": "betti_vector_v0_1",
        "status": "PASS" if controls_pass else "MIXED_OPEN",
        "classification": "TOPOLOGICAL_DIAGNOSTIC_SIDECAR",
        "labels": labels,
        "spec": load_spec(),
        "config": asdict(config),
        "reference_signature": {
            "Betti_0": reference_betti["Betti_0"],
            "Betti_1": reference_betti["Betti_1"],
            "nodes": reference_betti["nodes"],
            "edges": reference_betti["edges"],
        },
        "reference_component_size_distribution": reference_betti["component_size_distribution"],
        "condition_results": condition_results,
        "claim_ceiling": "Graph-topology diagnostic only; no Moonshine proof, local Lean theorem, or physical bridge claim.",
    }
    result["result_hash"] = stable_hash("betti_vector", result)

    write_csv(output_dir / CONTROL_RESULTS_PATH.name, control_rows, CONTROL_FIELDS)
    write_json(output_dir / RESULT_PATH.name, result)
    write_report(output_dir / REPORT_PATH.name, result, control_rows)
    return result


def write_report(path: Path, result: dict[str, Any], control_rows: list[dict[str, Any]]) -> None:
    reference = result["reference_signature"]
    lines = [
        "# Betti Vector Diagnostic Report",
        "",
        "## Status",
        "",
        f"- Result: `{result['status']}`",
        f"- Classification: `{result['classification']}`",
        f"- Result hash: `{result['result_hash']}`",
        f"- Reference `Betti_0`: `{reference['Betti_0']}`",
        f"- Reference `Betti_1`: `{reference['Betti_1']}`",
        f"- Reference nodes: `{reference['nodes']}`",
        f"- Reference edges: `{reference['edges']}`",
        f"- Reference component sizes: `{result['reference_component_size_distribution']}`",
        "",
        "## Labels",
        "",
        *[f"- `{label}`" for label in result["labels"]],
        "",
        "## Control Results",
        "",
        "| Condition | Nodes | Edges | Betti_0 | Betti_1 | Betti vector distance | Edge distance | Status |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in control_rows:
        lines.append(
            "| {condition} | {nodes} | {edges} | {betti0} | {betti1} | {betti_vector_distance} | {edge_jaccard_distance} | `{observed_status}` |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "`Betti_1` is computed here only as the independent cycle count of the declared finite undirected graph:",
            "",
            "```text",
            "Betti_1 = E - V + C",
            "```",
            "",
            "This is graph topology, not Moonshine topology. It does not prove Monstrous Moonshine, verify the external Lean SVG claim, or establish any physical bridge.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    result = run_betti_vector(BettiVectorConfig(), args.output_dir)
    print(json.dumps({"status": result["status"], "result_hash": result["result_hash"]}, sort_keys=True))


if __name__ == "__main__":
    main()
