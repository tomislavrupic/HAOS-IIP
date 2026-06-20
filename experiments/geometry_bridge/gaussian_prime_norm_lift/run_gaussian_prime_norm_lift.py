#!/usr/bin/env python3
"""Gaussian-prime norm-lift geometry diagnostic.

This sidecar turns the Gaussian integer lattice into a bounded HAOS-IIP style
geometry calibration. It measures whether norm shells, Gaussian-prime classes,
and graph Laplacian telemetry remain recoverable under matched controls.

It is arithmetic/synthetic geometry instrumentation only. It does not establish
physics, Monster moonshine, or a continuum bridge.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, deque
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np


ROOT = Path(__file__).resolve().parent
REPO_ROOT = Path(__file__).resolve().parents[3]

PRECOMMITMENT_PATH = ROOT / "precommitment_contract.json"
SOURCE_MANIFEST_PATH = ROOT / "source_manifest.json"
NODE_TABLE_PATH = ROOT / "gaussian_prime_nodes.csv"
COMPONENT_SCORES_PATH = ROOT / "component_scores.csv"
CONTROL_RESULTS_PATH = ROOT / "control_results.csv"
REPORT_PATH = ROOT / "gaussian_prime_norm_lift_report.md"
RESULT_PATH = ROOT / "gaussian_prime_norm_lift_result.json"

CONDITIONS = (
    "known_positive",
    "rotation_invariance_control",
    "class_shuffle_control",
    "norm_shuffle_control",
    "topology_destroyed_control",
)

CLASS_ORDER = (
    "origin",
    "unit",
    "ramified_norm2",
    "inert_axis_prime",
    "split_norm_prime",
    "composite_or_nonprime",
)

PRIME_CLASSES = {"ramified_norm2", "inert_axis_prime", "split_norm_prime"}

COMPONENT_SCORE_FIELDS = ["condition", "component", "metric", "value"]
CONTROL_RESULT_FIELDS = [
    "condition",
    "expected_response",
    "shell_distance",
    "spectral_distance",
    "norm_lift_distance",
    "cochain_distance",
    "observed_status",
]
NODE_FIELDS = ["node_id", "a", "b", "norm", "gaussian_class"]


@dataclass(frozen=True)
class GaussianPrimeNormLiftConfig:
    version: str = "gaussian-prime-norm-lift-v0.1"
    radius: int = 8
    eigenvalue_count: int = 12
    invariance_max: float = 1.0e-9
    class_shuffle_shell_min: float = 0.20
    norm_shuffle_lift_min: float = 0.10
    norm_shuffle_spectral_min: float = 0.02
    topology_spectral_min: float = 0.05
    topology_cochain_min: float = 0.05


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Gaussian-prime norm-lift geometry diagnostic.")
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


def is_rational_prime(value: int) -> bool:
    if value < 2:
        return False
    if value == 2:
        return True
    if value % 2 == 0:
        return False
    limit = int(math.sqrt(value))
    for divisor in range(3, limit + 1, 2):
        if value % divisor == 0:
            return False
    return True


def gaussian_prime_class(a: int, b: int) -> str:
    norm = a * a + b * b
    if norm == 0:
        return "origin"
    if norm == 1:
        return "unit"
    if norm == 2 and abs(a) == 1 and abs(b) == 1:
        return "ramified_norm2"
    if a == 0 or b == 0:
        axis_value = abs(a if b == 0 else b)
        if is_rational_prime(axis_value) and axis_value % 4 == 3:
            return "inert_axis_prime"
        return "composite_or_nonprime"
    if is_rational_prime(norm):
        return "split_norm_prime"
    return "composite_or_nonprime"


def build_target_nodes(radius: int) -> list[dict[str, Any]]:
    nodes: list[dict[str, Any]] = []
    for a in range(-radius, radius + 1):
        for b in range(-radius, radius + 1):
            norm = a * a + b * b
            nodes.append(
                {
                    "a": a,
                    "b": b,
                    "norm": norm,
                    "gaussian_class": gaussian_prime_class(a, b),
                }
            )
    return sorted(nodes, key=lambda row: (row["a"], row["b"]))


def condition_nodes(condition: str, config: GaussianPrimeNormLiftConfig) -> list[dict[str, Any]]:
    nodes = build_target_nodes(config.radius)
    if condition in {"known_positive", "topology_destroyed_control"}:
        return nodes
    if condition == "rotation_invariance_control":
        rotated = [
            {
                "a": -row["b"],
                "b": row["a"],
                "norm": row["norm"],
                "gaussian_class": gaussian_prime_class(-row["b"], row["a"]),
            }
            for row in nodes
        ]
        return sorted(rotated, key=lambda row: (row["a"], row["b"]))
    if condition == "class_shuffle_control":
        shuffled = [dict(row) for row in nodes]
        classes = [row["gaussian_class"] for row in nodes]
        order = sorted(range(len(classes)), key=lambda idx: stable_unit("class_shuffle", idx))
        shuffled_classes = [classes[idx] for idx in order]
        for row, label in zip(shuffled, shuffled_classes):
            row["gaussian_class"] = label
        return shuffled
    if condition == "norm_shuffle_control":
        shuffled = [dict(row) for row in nodes]
        norms = [row["norm"] for row in nodes]
        order = sorted(range(len(norms)), key=lambda idx: stable_unit("norm_shuffle", idx))
        shuffled_norms = [norms[idx] for idx in order]
        for row, norm in zip(shuffled, shuffled_norms):
            row["norm"] = norm
        return shuffled
    raise ValueError(f"unknown condition {condition}")


def grid_edges(nodes: list[dict[str, Any]]) -> list[tuple[int, int]]:
    index = {(row["a"], row["b"]): idx for idx, row in enumerate(nodes)}
    edges: list[tuple[int, int]] = []
    for idx, row in enumerate(nodes):
        for neighbor in ((row["a"] + 1, row["b"]), (row["a"], row["b"] + 1)):
            if neighbor in index:
                edges.append((idx, index[neighbor]))
    return edges


def topology_destroyed_edges(node_count: int, edge_count: int) -> list[tuple[int, int]]:
    candidates = [(i, j) for i in range(node_count) for j in range(i + 1, node_count)]
    ranked = sorted(candidates, key=lambda edge: stable_unit("topology_destroy", edge[0], edge[1]))
    return ranked[:edge_count]


def condition_edges(condition: str, nodes: list[dict[str, Any]]) -> list[tuple[int, int]]:
    target_edges = grid_edges(build_target_nodes(int((math.sqrt(len(nodes)) - 1) / 2)))
    if condition == "topology_destroyed_control":
        return topology_destroyed_edges(len(nodes), len(target_edges))
    return grid_edges(nodes)


def weighted_adjacency(nodes: list[dict[str, Any]], edges: list[tuple[int, int]]) -> np.ndarray:
    adjacency = np.zeros((len(nodes), len(nodes)), dtype=float)
    norms = np.array([row["norm"] for row in nodes], dtype=float)
    for i, j in edges:
        weight = 1.0 / (1.0 + abs(norms[i] - norms[j]))
        adjacency[i, j] = adjacency[j, i] = weight
    return adjacency


def normalized_laplacian(adjacency: np.ndarray) -> np.ndarray:
    degree = np.sum(adjacency, axis=1)
    inv_sqrt = np.zeros_like(degree)
    mask = degree > 1.0e-12
    inv_sqrt[mask] = 1.0 / np.sqrt(degree[mask])
    return np.eye(adjacency.shape[0]) - (inv_sqrt[:, None] * adjacency * inv_sqrt[None, :])


def shell_component(nodes: list[dict[str, Any]]) -> dict[str, Any]:
    counts = Counter(row["gaussian_class"] for row in nodes)
    class_vector = [counts.get(label, 0) for label in CLASS_ORDER]
    prime_rows = [row for row in nodes if row["gaussian_class"] in PRIME_CLASSES]
    alignment_count = sum(1 for row in prime_rows if prime_class_matches_shell(row))
    shell_counts = Counter(row["norm"] for row in prime_rows)
    return {
        "class_counts": class_vector,
        "prime_count": len(prime_rows),
        "distinct_prime_shells": len(shell_counts),
        "prime_shell_alignment": alignment_count / len(prime_rows) if prime_rows else 0.0,
        "prime_shell_entropy": entropy(list(shell_counts.values())),
    }


def prime_class_matches_shell(row: dict[str, Any]) -> bool:
    label = row["gaussian_class"]
    norm = int(row["norm"])
    a = int(row["a"])
    b = int(row["b"])
    if label == "ramified_norm2":
        return norm == 2
    if label == "inert_axis_prime":
        axis_value = abs(a if b == 0 else b)
        return (a == 0 or b == 0) and is_rational_prime(axis_value) and axis_value % 4 == 3
    if label == "split_norm_prime":
        return a != 0 and b != 0 and is_rational_prime(norm) and norm % 4 == 1
    return False


def entropy(counts: list[int]) -> float:
    total = sum(counts)
    if total == 0:
        return 0.0
    result = 0.0
    for count in counts:
        if count:
            probability = count / total
            result -= probability * math.log(probability, 2)
    return float(result)


def spectral_component(adjacency: np.ndarray, config: GaussianPrimeNormLiftConfig) -> dict[str, Any]:
    eigenvalues = np.linalg.eigvalsh(normalized_laplacian(adjacency))
    eigenvalues = np.sort(eigenvalues)
    selected = eigenvalues[: config.eigenvalue_count]
    return {
        "eigenvalues": [float(value) for value in selected],
        "spectral_gap": float(eigenvalues[1] - eigenvalues[0]) if len(eigenvalues) > 1 else 0.0,
        "zero_mode_count": int(np.sum(np.abs(eigenvalues) < 1.0e-8)),
        "low_mode_sum": float(np.sum(selected)),
    }


def norm_lift_component(nodes: list[dict[str, Any]], edges: list[tuple[int, int]]) -> dict[str, Any]:
    norms = np.array([row["norm"] for row in nodes], dtype=float)
    norm_gradients = [abs(norms[i] - norms[j]) for i, j in edges]
    origin_index = next(idx for idx, row in enumerate(nodes) if row["a"] == 0 and row["b"] == 0)
    graph_distance = bfs_distances(len(nodes), edges, origin_index)
    finite = np.isfinite(graph_distance)
    radial = np.sqrt(norms)
    correlation = pearson(graph_distance[finite], radial[finite]) if np.any(finite) else 0.0
    return {
        "edge_norm_gradient_mean": mean(norm_gradients),
        "edge_norm_gradient_std": std(norm_gradients),
        "origin_graph_radial_correlation": correlation,
        "finite_origin_reachability": float(np.mean(finite)),
    }


def cochain_component(nodes: list[dict[str, Any]], edges: list[tuple[int, int]]) -> dict[str, Any]:
    norms = np.array([row["norm"] for row in nodes], dtype=float)
    divergence = np.zeros(len(nodes), dtype=float)
    flow_values = []
    for i, j in edges:
        flow = norms[j] - norms[i]
        divergence[i] -= flow
        divergence[j] += flow
        flow_values.append(abs(flow))
    return {
        "edge_flow_abs_mean": mean(flow_values),
        "edge_flow_abs_std": std(flow_values),
        "divergence_l2": float(np.linalg.norm(divergence) / max(len(nodes), 1)),
        "divergence_abs_mean": float(np.mean(np.abs(divergence))) if len(nodes) else 0.0,
    }


def bfs_distances(node_count: int, edges: list[tuple[int, int]], source: int) -> np.ndarray:
    neighbors: list[list[int]] = [[] for _ in range(node_count)]
    for i, j in edges:
        neighbors[i].append(j)
        neighbors[j].append(i)
    distances = np.full(node_count, np.inf, dtype=float)
    distances[source] = 0.0
    queue: deque[int] = deque([source])
    while queue:
        node = queue.popleft()
        for neighbor in neighbors[node]:
            if not np.isfinite(distances[neighbor]):
                distances[neighbor] = distances[node] + 1.0
                queue.append(neighbor)
    return distances


def pearson(left: np.ndarray, right: np.ndarray) -> float:
    if len(left) < 2:
        return 0.0
    a = left - float(np.mean(left))
    b = right - float(np.mean(right))
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a, b) / denom) if denom > 1.0e-12 else 0.0


def mean(values: Iterable[float]) -> float:
    items = list(values)
    return float(sum(items) / len(items)) if items else 0.0


def std(values: Iterable[float]) -> float:
    items = list(values)
    if not items:
        return 0.0
    avg = mean(items)
    return float(math.sqrt(sum((value - avg) ** 2 for value in items) / len(items)))


def vector_distance(left: Iterable[float], right: Iterable[float]) -> float:
    a = np.array(list(left), dtype=float)
    b = np.array(list(right), dtype=float)
    denom = max(float(np.linalg.norm(a)), 1.0e-12)
    return float(np.linalg.norm(a - b) / denom)


def shell_distance(reference: dict[str, Any], observed: dict[str, Any]) -> float:
    class_part = vector_distance(reference["class_counts"], observed["class_counts"])
    fields = ["prime_count", "distinct_prime_shells", "prime_shell_alignment", "prime_shell_entropy"]
    scalar_part = vector_distance([reference[field] for field in fields], [observed[field] for field in fields])
    return float(class_part + scalar_part)


def spectral_distance(reference: dict[str, Any], observed: dict[str, Any]) -> float:
    scalar = vector_distance(
        [reference["spectral_gap"], reference["zero_mode_count"], reference["low_mode_sum"]],
        [observed["spectral_gap"], observed["zero_mode_count"], observed["low_mode_sum"]],
    )
    return float(vector_distance(reference["eigenvalues"], observed["eigenvalues"]) + scalar)


def norm_lift_distance(reference: dict[str, Any], observed: dict[str, Any]) -> float:
    fields = [
        "edge_norm_gradient_mean",
        "edge_norm_gradient_std",
        "origin_graph_radial_correlation",
        "finite_origin_reachability",
    ]
    return vector_distance([reference[field] for field in fields], [observed[field] for field in fields])


def cochain_distance(reference: dict[str, Any], observed: dict[str, Any]) -> float:
    fields = ["edge_flow_abs_mean", "edge_flow_abs_std", "divergence_l2", "divergence_abs_mean"]
    return vector_distance([reference[field] for field in fields], [observed[field] for field in fields])


def expected_response(condition: str) -> str:
    return {
        "known_positive": "all component distances are zero",
        "rotation_invariance_control": "Gaussian-prime and norm-lift telemetry remains invariant under lattice rotation",
        "class_shuffle_control": "prime-shell semantic alignment degrades while graph geometry is mostly preserved",
        "norm_shuffle_control": "norm-lift and weighted spectral telemetry move under shell mismatch",
        "topology_destroyed_control": "spectral and cochain telemetry move under lattice-neighborhood destruction",
    }[condition]


def observed_status(condition: str, distances: dict[str, float], config: GaussianPrimeNormLiftConfig) -> str:
    if condition in {"known_positive", "rotation_invariance_control"}:
        return "PASS" if max(distances.values()) <= config.invariance_max else "FAIL_INVARIANCE"
    if condition == "class_shuffle_control":
        return "PASS" if distances["shell"] >= config.class_shuffle_shell_min else "FAIL_CLASS_SHUFFLE"
    if condition == "norm_shuffle_control":
        if distances["norm_lift"] >= config.norm_shuffle_lift_min and distances["spectral"] >= config.norm_shuffle_spectral_min:
            return "PASS"
        return "FAIL_NORM_SHUFFLE"
    if condition == "topology_destroyed_control":
        if distances["spectral"] >= config.topology_spectral_min and distances["cochain"] >= config.topology_cochain_min:
            return "PASS"
        return "FAIL_TOPOLOGY_DESTROYED"
    raise ValueError(f"unknown condition {condition}")


def component_rows(condition: str, components: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for component, metrics in components.items():
        for metric, value in metrics.items():
            if isinstance(value, list):
                for index, item in enumerate(value):
                    rows.append({"condition": condition, "component": component, "metric": f"{metric}_{index}", "value": f"{float(item):.12g}"})
            else:
                rows.append({"condition": condition, "component": component, "metric": metric, "value": f"{float(value):.12g}"})
    return rows


def node_rows(nodes: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "node_id": idx,
            "a": row["a"],
            "b": row["b"],
            "norm": row["norm"],
            "gaussian_class": row["gaussian_class"],
        }
        for idx, row in enumerate(nodes)
    ]


def source_manifest() -> dict[str, Any]:
    return {
        "source": "Gaussian integer lattice synthetic arithmetic geometry sidecar",
        "code": repo_rel(ROOT / "run_gaussian_prime_norm_lift.py"),
        "mathematical_objects": [
            "Gaussian integers Z[i]",
            "norm N(a+bi)=a^2+b^2",
            "Gaussian-prime ramified/inert/split classification",
            "weighted grid graph Laplacian",
        ],
        "external_data": "none",
        "monster_data": "not included in v0.1; future visualization hook only",
    }


def precommitment_contract(config: GaussianPrimeNormLiftConfig) -> dict[str, Any]:
    return {
        "name": "Gaussian Prime Norm-Lift Geometry v0.1",
        "status": "PRECOMMITTED_ARITHMETIC_GEOMETRY_CALIBRATION",
        "purpose": "Test whether discrete arithmetic invariants form stable norm-lift and spectral telemetry under matched controls.",
        "source_artifacts": source_manifest(),
        "conditions": {condition: expected_response(condition) for condition in CONDITIONS},
        "thresholds": asdict(config),
        "components": [
            "gaussian_prime_shells",
            "weighted_laplacian_spectrum",
            "norm_lift_graph_telemetry",
            "cochain_edge_flow_telemetry",
        ],
        "claim_ceiling": "SYNTHETIC_ARITHMETIC_GEOMETRY_CALIBRATION",
        "non_claims": [
            "not a physical bridge",
            "not Monster moonshine",
            "not a continuum limit",
            "not a quantum, gravity, or field-theory derivation",
            "not empirical validation",
            "not a change to frozen HAOS-IIP phases",
        ],
    }


def run_condition(condition: str, config: GaussianPrimeNormLiftConfig) -> dict[str, Any]:
    nodes = condition_nodes(condition, config)
    edges = condition_edges(condition, nodes)
    adjacency = weighted_adjacency(nodes, edges)
    return {
        "nodes": nodes,
        "edges": edges,
        "components": {
            "shell": shell_component(nodes),
            "spectral": spectral_component(adjacency, config),
            "norm_lift": norm_lift_component(nodes, edges),
            "cochain": cochain_component(nodes, edges),
        },
    }


def run_bridge(config: GaussianPrimeNormLiftConfig, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    reference = run_condition("known_positive", config)
    score_rows: list[dict[str, Any]] = []
    control_rows: list[dict[str, Any]] = []
    condition_results: dict[str, Any] = {}

    for condition in CONDITIONS:
        observed = run_condition(condition, config)
        components = observed["components"]
        distances = {
            "shell": shell_distance(reference["components"]["shell"], components["shell"]),
            "spectral": spectral_distance(reference["components"]["spectral"], components["spectral"]),
            "norm_lift": norm_lift_distance(reference["components"]["norm_lift"], components["norm_lift"]),
            "cochain": cochain_distance(reference["components"]["cochain"], components["cochain"]),
        }
        status = observed_status(condition, distances, config)
        condition_results[condition] = {
            "components": components,
            "distances": distances,
            "expected_response": expected_response(condition),
            "observed_status": status,
        }
        score_rows.extend(component_rows(condition, components))
        control_rows.append(
            {
                "condition": condition,
                "expected_response": expected_response(condition),
                "shell_distance": f"{distances['shell']:.12g}",
                "spectral_distance": f"{distances['spectral']:.12g}",
                "norm_lift_distance": f"{distances['norm_lift']:.12g}",
                "cochain_distance": f"{distances['cochain']:.12g}",
                "observed_status": status,
            }
        )

    controls_pass = all(result["observed_status"] == "PASS" for result in condition_results.values())
    labels = [
        "GAUSSIAN_PRIME_NORM_LIFT_BUILT",
        "ARITHMETIC_GEOMETRY_TELEMETRY_AVAILABLE",
        "SPECTRAL_TELEMETRY_AVAILABLE",
        "COCHAIN_FLOW_TELEMETRY_AVAILABLE",
        "PHYSICAL_BRIDGE_NOT_ESTABLISHED",
        "MONSTER_MOONSHINE_NOT_TESTED",
        "CLAIM_GATED_SYNTHETIC_CALIBRATION",
    ]
    labels.append("COMPONENT_CONTROLS_PASS" if controls_pass else "COMPONENT_CONTROLS_PARTIAL")

    result_payload = {
        "bridge_id": "gaussian_prime_norm_lift_v0_1",
        "status": "PASS" if controls_pass else "MIXED_OPEN",
        "classification": "SYNTHETIC_ARITHMETIC_GEOMETRY_CALIBRATION",
        "labels": labels,
        "condition_results": condition_results,
        "claim_ceiling": "No physics or Monster claim; arithmetic geometry calibration only.",
    }
    result_payload["result_hash"] = stable_hash("gaussian_norm_lift", result_payload)

    write_json(output_dir / PRECOMMITMENT_PATH.name, precommitment_contract(config))
    write_json(output_dir / SOURCE_MANIFEST_PATH.name, source_manifest())
    write_csv(output_dir / NODE_TABLE_PATH.name, node_rows(reference["nodes"]), NODE_FIELDS)
    write_csv(output_dir / COMPONENT_SCORES_PATH.name, score_rows, COMPONENT_SCORE_FIELDS)
    write_csv(output_dir / CONTROL_RESULTS_PATH.name, control_rows, CONTROL_RESULT_FIELDS)
    write_json(output_dir / RESULT_PATH.name, result_payload)
    write_report(output_dir / REPORT_PATH.name, result_payload, control_rows)
    return result_payload


def write_report(path: Path, result: dict[str, Any], control_rows: list[dict[str, Any]]) -> None:
    lines = [
        "# Gaussian Prime Norm-Lift Geometry Report",
        "",
        "## Status",
        "",
        f"- Result: `{result['status']}`",
        f"- Classification: `{result['classification']}`",
        f"- Result hash: `{result['result_hash']}`",
        "",
        "## Labels",
        "",
        *[f"- `{label}`" for label in result["labels"]],
        "",
        "## Control Results",
        "",
        "| Condition | Shell | Spectral | Norm lift | Cochain | Status |",
        "| --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in control_rows:
        lines.append(
            "| {condition} | {shell_distance} | {spectral_distance} | {norm_lift_distance} | {cochain_distance} | `{observed_status}` |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "This sidecar uses Gaussian integers as an arithmetic geometry calibration case.",
            "It does not establish a physical bridge, continuum limit, Monster moonshine connection, quantum result, gravity result, or field-theory derivation.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    result = run_bridge(GaussianPrimeNormLiftConfig(), args.output_dir)
    print(json.dumps({"status": result["status"], "result_hash": result["result_hash"]}, sort_keys=True))


if __name__ == "__main__":
    main()

