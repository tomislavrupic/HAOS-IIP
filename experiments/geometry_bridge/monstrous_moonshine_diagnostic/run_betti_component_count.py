#!/usr/bin/env python3
"""Betti_0 component-count diagnostic for the moonshine sidecar.

This script implements the local graph rule promised by
betti_diagram_haos_connection.md. It treats Betti_0 as connected-component
count for a declared relation graph over Gaussian-prime representatives of the
15 supersingular primes.

It is a topological diagnostic sidecar only. It does not prove Monstrous
Moonshine and does not verify the external Lean/SVG claim.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import deque
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parent
REPO_ROOT = Path(__file__).resolve().parents[3]

SPEC_PATH = ROOT / "betti_relation_graph_spec.json"
NODE_TABLE_PATH = ROOT / "betti_relation_graph_nodes.csv"
EDGE_TABLE_PATH = ROOT / "betti_relation_graph_edges.csv"
CONTROL_RESULTS_PATH = ROOT / "betti_control_results.csv"
REPORT_PATH = ROOT / "betti_component_count_report.md"
RESULT_PATH = ROOT / "betti_component_count_result.json"

SUPERSINGULAR_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71)
CONDITIONS = (
    "known_positive",
    "d4_rotate_90",
    "d4_reflect_real",
    "isomorphism_relabel",
    "threshold_mutation_control",
    "topology_destroyed_control",
    "support_replacement_control",
)

NODE_FIELDS = ["node_id", "prime", "x", "y", "gaussian_class"]
EDGE_FIELDS = ["source", "target", "distance"]
CONTROL_FIELDS = [
    "condition",
    "expected_response",
    "component_count",
    "betti0_distance",
    "edge_jaccard_distance",
    "observed_status",
]


@dataclass(frozen=True)
class BettiComponentConfig:
    version: str = "betti-component-count-v0.1"
    threshold: float = 8.0
    mutated_threshold: float = 4.0
    invariance_max: float = 1.0e-12
    destructive_betti_min: float = 1.0
    destructive_edge_min: float = 0.10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Betti_0 component-count diagnostic.")
    parser.add_argument("--output-dir", type=Path, default=ROOT)
    return parser.parse_args()


def stable_hash(prefix: str, payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(encoded).hexdigest()[:24]}"


def stable_unit(*parts: Any) -> float:
    digest = hashlib.sha256("|".join(str(part) for part in parts).encode("utf-8")).hexdigest()
    return int(digest[:16], 16) / float(16**16 - 1)


def repo_rel(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def gaussian_class(prime: int) -> str:
    if prime == 2:
        return "ramified"
    if prime % 4 == 1:
        return "split"
    if prime % 4 == 3:
        return "inert"
    return "unknown"


def split_prime_pair(prime: int) -> tuple[int, int]:
    for a in range(1, int(math.sqrt(prime)) + 1):
        for b in range(a, int(math.sqrt(prime)) + 1):
            if a * a + b * b == prime:
                return (a, b)
    raise ValueError(f"no split representative found for {prime}")


def representative_for_prime(prime: int) -> tuple[int, int]:
    if prime == 2:
        return (1, 1)
    if prime % 4 == 3:
        return (prime, 0)
    if prime % 4 == 1:
        return split_prime_pair(prime)
    raise ValueError(f"unsupported prime {prime}")


def base_nodes(primes: Iterable[int] = SUPERSINGULAR_PRIMES) -> list[dict[str, Any]]:
    nodes = []
    for prime in primes:
        x, y = representative_for_prime(prime)
        nodes.append({"prime": prime, "x": x, "y": y, "gaussian_class": gaussian_class(prime)})
    return nodes


def transform_nodes(nodes: list[dict[str, Any]], condition: str) -> list[dict[str, Any]]:
    if condition == "d4_rotate_90":
        return [{**node, "x": -node["y"], "y": node["x"]} for node in nodes]
    if condition == "d4_reflect_real":
        return [{**node, "x": node["x"], "y": -node["y"]} for node in nodes]
    if condition == "isomorphism_relabel":
        order = sorted(range(len(nodes)), key=lambda idx: stable_unit("relabel", idx))
        return [nodes[idx] for idx in order]
    return nodes


def condition_nodes(condition: str) -> list[dict[str, Any]]:
    if condition == "support_replacement_control":
        primes = [prime for prime in SUPERSINGULAR_PRIMES if prime != 71] + [73]
        return base_nodes(primes)
    return transform_nodes(base_nodes(), condition)


def distance(left: dict[str, Any], right: dict[str, Any]) -> float:
    return math.hypot(float(left["x"]) - float(right["x"]), float(left["y"]) - float(right["y"]))


def relation_edges(nodes: list[dict[str, Any]], threshold: float) -> list[tuple[int, int, float]]:
    edges = []
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            dist = distance(nodes[i], nodes[j])
            if dist <= threshold:
                edges.append((i, j, dist))
    return edges


def topology_destroyed_edges(node_count: int, edge_count: int) -> list[tuple[int, int, float]]:
    candidates = [(i, j) for i in range(node_count) for j in range(i + 1, node_count)]
    ranked = sorted(candidates, key=lambda edge: stable_unit("betti_topology_destroy", edge[0], edge[1]))
    return [(i, j, 1.0) for i, j in ranked[:edge_count]]


def condition_edges(condition: str, nodes: list[dict[str, Any]], config: BettiComponentConfig, reference_edge_count: int) -> list[tuple[int, int, float]]:
    if condition == "threshold_mutation_control":
        return relation_edges(nodes, config.mutated_threshold)
    if condition == "topology_destroyed_control":
        return topology_destroyed_edges(len(nodes), reference_edge_count)
    return relation_edges(nodes, config.threshold)


def canonical_edge_signature(nodes: list[dict[str, Any]], edges: list[tuple[int, int, float]]) -> set[tuple[int, int]]:
    prime_by_index = [int(node["prime"]) for node in nodes]
    signatures = set()
    for i, j, _ in edges:
        left = prime_by_index[i]
        right = prime_by_index[j]
        signatures.add((min(left, right), max(left, right)))
    return signatures


def component_count(node_count: int, edges: list[tuple[int, int, float]]) -> int:
    neighbors: list[list[int]] = [[] for _ in range(node_count)]
    for i, j, _ in edges:
        neighbors[i].append(j)
        neighbors[j].append(i)
    seen = [False] * node_count
    count = 0
    for start in range(node_count):
        if seen[start]:
            continue
        count += 1
        queue: deque[int] = deque([start])
        seen[start] = True
        while queue:
            node = queue.popleft()
            for neighbor in neighbors[node]:
                if not seen[neighbor]:
                    seen[neighbor] = True
                    queue.append(neighbor)
    return count


def jaccard_distance(left: set[tuple[int, int]], right: set[tuple[int, int]]) -> float:
    union = left | right
    if not union:
        return 0.0
    return float(1.0 - len(left & right) / len(union))


def expected_response(condition: str) -> str:
    return {
        "known_positive": "self Betti_0 and edge signature recover exactly",
        "d4_rotate_90": "D4 rotation preserves Betti_0 and edge signature",
        "d4_reflect_real": "D4 reflection preserves Betti_0 and edge signature",
        "isomorphism_relabel": "graph-isomorphism relabel preserves Betti_0 and edge signature",
        "threshold_mutation_control": "threshold mutation should move Betti_0 or relation edges",
        "topology_destroyed_control": "rewiring should move Betti_0 or relation edges",
        "support_replacement_control": "support replacement should move Betti_0 or relation edges",
    }[condition]


def observed_status(condition: str, betti0_distance: float, edge_distance: float, config: BettiComponentConfig) -> str:
    if condition in {"known_positive", "d4_rotate_90", "d4_reflect_real", "isomorphism_relabel"}:
        if betti0_distance <= config.invariance_max and edge_distance <= config.invariance_max:
            return "PASS"
        return "FAIL_SYMMETRY_STABILITY"
    if betti0_distance >= config.destructive_betti_min or edge_distance >= config.destructive_edge_min:
        return "PASS"
    return "FAIL_DESTRUCTIVE_CONTROL"


def node_rows(nodes: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "node_id": idx,
            "prime": node["prime"],
            "x": node["x"],
            "y": node["y"],
            "gaussian_class": node["gaussian_class"],
        }
        for idx, node in enumerate(nodes)
    ]


def edge_rows(edges: list[tuple[int, int, float]]) -> list[dict[str, Any]]:
    return [{"source": i, "target": j, "distance": f"{dist:.12g}"} for i, j, dist in edges]


def load_spec() -> dict[str, Any]:
    return json.loads(SPEC_PATH.read_text(encoding="utf-8"))


def run_betti(config: BettiComponentConfig, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    reference_nodes = condition_nodes("known_positive")
    reference_edges = relation_edges(reference_nodes, config.threshold)
    reference_betti0 = component_count(len(reference_nodes), reference_edges)
    reference_signature = canonical_edge_signature(reference_nodes, reference_edges)

    control_rows = []
    condition_results: dict[str, Any] = {}
    for condition in CONDITIONS:
        nodes = condition_nodes(condition)
        edges = condition_edges(condition, nodes, config, len(reference_edges))
        betti0 = component_count(len(nodes), edges)
        edge_signature = canonical_edge_signature(nodes, edges)
        betti0_distance = abs(betti0 - reference_betti0)
        edge_distance = jaccard_distance(reference_signature, edge_signature)
        status = observed_status(condition, betti0_distance, edge_distance, config)
        condition_results[condition] = {
            "component_count": betti0,
            "edge_count": len(edges),
            "betti0_distance": betti0_distance,
            "edge_jaccard_distance": edge_distance,
            "expected_response": expected_response(condition),
            "observed_status": status,
        }
        control_rows.append(
            {
                "condition": condition,
                "expected_response": expected_response(condition),
                "component_count": betti0,
                "betti0_distance": f"{betti0_distance:.12g}",
                "edge_jaccard_distance": f"{edge_distance:.12g}",
                "observed_status": status,
            }
        )

    controls_pass = all(result["observed_status"] == "PASS" for result in condition_results.values())
    labels = [
        "BETTI_RELATION_GRAPH_SPEC_BUILT",
        "BETTI_0_COMPONENT_COUNT_AVAILABLE",
        "D4_SYMMETRY_CONTROLS_AVAILABLE",
        "DESTRUCTIVE_CONTROLS_AVAILABLE",
        "LEAN_THEOREM_NOT_INCLUDED",
        "PHYSICAL_BRIDGE_NOT_ESTABLISHED",
        "CLAIM_GATED_TOPOLOGICAL_DIAGNOSTIC",
    ]
    labels.append("BETTI_CONTROLS_PASS" if controls_pass else "BETTI_CONTROLS_PARTIAL")

    result_payload = {
        "bridge_id": "betti_component_count_v0_1",
        "status": "PASS" if controls_pass else "MIXED_OPEN",
        "classification": "TOPOLOGICAL_DIAGNOSTIC_SIDECAR",
        "labels": labels,
        "spec": load_spec(),
        "config": asdict(config),
        "reference": {
            "component_count": reference_betti0,
            "edge_count": len(reference_edges),
            "node_count": len(reference_nodes),
        },
        "condition_results": condition_results,
        "claim_ceiling": "No Moonshine proof, local Lean theorem, or physical bridge claim.",
    }
    result_payload["result_hash"] = stable_hash("betti_component", result_payload)

    write_csv(output_dir / NODE_TABLE_PATH.name, node_rows(reference_nodes), NODE_FIELDS)
    write_csv(output_dir / EDGE_TABLE_PATH.name, edge_rows(reference_edges), EDGE_FIELDS)
    write_csv(output_dir / CONTROL_RESULTS_PATH.name, control_rows, CONTROL_FIELDS)
    write_json(output_dir / RESULT_PATH.name, result_payload)
    write_report(output_dir / REPORT_PATH.name, result_payload, control_rows)
    return result_payload


def write_report(path: Path, result: dict[str, Any], control_rows: list[dict[str, Any]]) -> None:
    lines = [
        "# Betti Component-Count Diagnostic Report",
        "",
        "## Status",
        "",
        f"- Result: `{result['status']}`",
        f"- Classification: `{result['classification']}`",
        f"- Result hash: `{result['result_hash']}`",
        f"- Reference `Betti_0`: `{result['reference']['component_count']}`",
        f"- Reference edges: `{result['reference']['edge_count']}`",
        "",
        "## Labels",
        "",
        *[f"- `{label}`" for label in result["labels"]],
        "",
        "## Control Results",
        "",
        "| Condition | Component count | Betti_0 distance | Edge distance | Status |",
        "| --- | ---: | ---: | ---: | --- |",
    ]
    for row in control_rows:
        lines.append(
            "| {condition} | {component_count} | {betti0_distance} | {edge_jaccard_distance} | `{observed_status}` |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "This script implements a local Betti_0 / component-count diagnostic for a declared arithmetic relation graph.",
            "It does not prove Monstrous Moonshine, verify the external Lean SVG claim, or establish any physical bridge.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    result = run_betti(BettiComponentConfig(), args.output_dir)
    print(json.dumps({"status": result["status"], "result_hash": result["result_hash"]}, sort_keys=True))


if __name__ == "__main__":
    main()

