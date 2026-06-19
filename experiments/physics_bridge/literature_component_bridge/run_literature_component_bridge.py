#!/usr/bin/env python3
"""Literature-derived HAOS-IIP bridge component calibration.

This sidecar turns the bridge-literature extraction into an executable
instrument test. It combines spectral, Hodge, curvature, transport, and kernel
comparison components on one deterministic synthetic relational graph plus
matched controls.

It does not establish a physical bridge.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy.optimize import linprog
from scipy.sparse.csgraph import shortest_path


ROOT = Path(__file__).resolve().parent
REPO_ROOT = Path(__file__).resolve().parents[3]

PRECOMMITMENT_PATH = ROOT / "precommitment_contract.json"
SOURCE_MANIFEST_PATH = ROOT / "source_manifest.json"
COMPONENT_SCORES_PATH = ROOT / "component_scores.csv"
CONTROL_RESULTS_PATH = ROOT / "control_results.csv"
REPORT_PATH = ROOT / "literature_component_bridge_report.md"
RESULT_PATH = ROOT / "literature_component_bridge_result.json"

CONDITIONS = (
    "known_positive",
    "label_permutation_control",
    "weight_shuffle_control",
    "topology_destroyed_control",
    "hodge_triangle_removed_control",
)

COMPONENT_SCORE_FIELDS = ["condition", "component", "metric", "value"]
CONTROL_RESULT_FIELDS = [
    "condition",
    "expected_response",
    "spectral_distance",
    "hodge_distance",
    "curvature_distance",
    "transport_kernel_distance",
    "observed_status",
]


@dataclass(frozen=True)
class BridgeComponentConfig:
    version: str = "literature-component-bridge-v0.1"
    spectral_distance_min: float = 0.12
    hodge_distance_min: float = 0.20
    curvature_distance_min: float = 0.08
    transport_distance_min: float = 0.02
    label_invariance_max: float = 1.0e-9
    sinkhorn_epsilon: float = 0.2
    mmd_bandwidth: float = 0.75


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the literature-derived bridge component sidecar.")
    parser.add_argument("--output-dir", type=Path, default=ROOT)
    return parser.parse_args()


def repo_rel(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def stable_hash(prefix: str, payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(encoded).hexdigest()[:24]}"


def stable_unit(*parts: Any) -> float:
    digest = hashlib.sha256("|".join(str(part) for part in parts).encode("utf-8")).hexdigest()
    return int(digest[:16], 16) / float(16**16 - 1)


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


def target_graph() -> tuple[np.ndarray, np.ndarray, bool]:
    """Return target adjacency, cluster labels, and whether 2-cells are present."""
    labels = np.repeat(np.arange(3), 4)
    n = len(labels)
    adjacency = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            if labels[i] == labels[j]:
                weight = 1.0 - 0.03 * ((i + j) % 3)
                adjacency[i, j] = adjacency[j, i] = weight
    bridges = ((3, 4, 0.22), (7, 8, 0.21), (11, 0, 0.20), (2, 6, 0.16), (5, 9, 0.15))
    for i, j, weight in bridges:
        adjacency[i, j] = adjacency[j, i] = weight
    return adjacency, labels, True


def condition_graph(condition: str) -> tuple[np.ndarray, np.ndarray, bool]:
    adjacency, labels, include_triangles = target_graph()
    if condition == "known_positive":
        return adjacency, labels, include_triangles
    if condition == "label_permutation_control":
        order = np.array([4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3])
        return adjacency[np.ix_(order, order)], labels[order], include_triangles
    if condition == "hodge_triangle_removed_control":
        return adjacency, labels, False
    if condition == "weight_shuffle_control":
        shuffled = np.zeros_like(adjacency)
        edges = edge_list(adjacency)
        weights = [adjacency[i, j] for i, j in edges]
        order = sorted(range(len(weights)), key=lambda idx: stable_unit("weight_shuffle", idx))
        for edge, weight_idx in zip(edges, order):
            i, j = edge
            shuffled[i, j] = shuffled[j, i] = weights[weight_idx]
        return shuffled, labels, include_triangles
    if condition == "topology_destroyed_control":
        n = adjacency.shape[0]
        destroyed = np.zeros_like(adjacency)
        edge_count = len(edge_list(adjacency))
        candidates = [(i, j) for i in range(n) for j in range(i + 1, n)]
        ranked = sorted(candidates, key=lambda edge: stable_unit("topology", *edge))
        source_weights = sorted((adjacency[i, j] for i, j in edge_list(adjacency)), reverse=True)
        for edge, weight in zip(ranked[:edge_count], source_weights):
            i, j = edge
            destroyed[i, j] = destroyed[j, i] = weight
        return destroyed, labels, include_triangles
    raise ValueError(f"unknown condition {condition}")


def edge_list(adjacency: np.ndarray) -> list[tuple[int, int]]:
    return [(i, j) for i in range(adjacency.shape[0]) for j in range(i + 1, adjacency.shape[1]) if adjacency[i, j] > 0]


def triangle_list(adjacency: np.ndarray, include_triangles: bool) -> list[tuple[int, int, int]]:
    if not include_triangles:
        return []
    n = adjacency.shape[0]
    triangles: list[tuple[int, int, int]] = []
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if adjacency[i, j] > 0 and adjacency[i, k] > 0 and adjacency[j, k] > 0:
                    triangles.append((i, j, k))
    return triangles


def normalized_laplacian(adjacency: np.ndarray) -> np.ndarray:
    degree = np.sum(adjacency, axis=1)
    inv_sqrt = np.zeros_like(degree)
    inv_sqrt[degree > 1.0e-12] = 1.0 / np.sqrt(degree[degree > 1.0e-12])
    return np.eye(adjacency.shape[0]) - (inv_sqrt[:, None] * adjacency * inv_sqrt[None, :])


def spectral_component(adjacency: np.ndarray, labels: np.ndarray) -> dict[str, Any]:
    laplacian = normalized_laplacian(adjacency)
    eigenvalues, eigenvectors = np.linalg.eigh(laplacian)
    order = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    embedding = eigenvectors[:, 1:4]
    same_distances: list[float] = []
    different_distances: list[float] = []
    for i in range(adjacency.shape[0]):
        for j in range(i + 1, adjacency.shape[0]):
            distance = float(np.linalg.norm(embedding[i] - embedding[j]))
            if labels[i] == labels[j]:
                same_distances.append(distance)
            else:
                different_distances.append(distance)
    ordering_score = pairwise_order_score(same_distances, different_distances, lower_is_better=True)
    return {
        "eigenvalues": [float(value) for value in eigenvalues[:8]],
        "spectral_gap": float(eigenvalues[1] - eigenvalues[0]),
        "low_mode_cluster_ordering": ordering_score,
        "fiedler_sign_balance": float(np.mean(eigenvectors[:, 1] > 0)),
    }


def pairwise_order_score(left: list[float], right: list[float], lower_is_better: bool) -> float:
    total = 0
    good = 0.0
    for a in left:
        for b in right:
            total += 1
            if a == b:
                good += 0.5
            elif (a < b) == lower_is_better:
                good += 1.0
    return good / float(total) if total else 0.0


def incidence_matrices(adjacency: np.ndarray, include_triangles: bool) -> tuple[np.ndarray, np.ndarray, list[tuple[int, int]]]:
    edges = edge_list(adjacency)
    edge_index = {edge: idx for idx, edge in enumerate(edges)}
    b1 = np.zeros((adjacency.shape[0], len(edges)), dtype=float)
    for col, (i, j) in enumerate(edges):
        b1[i, col] = -1.0
        b1[j, col] = 1.0
    triangles = triangle_list(adjacency, include_triangles)
    b2 = np.zeros((len(edges), len(triangles)), dtype=float)
    for col, (i, j, k) in enumerate(triangles):
        b2[edge_index[(j, k)], col] = 1.0
        b2[edge_index[(i, j)], col] = 1.0
        b2[edge_index[(i, k)], col] = -1.0
    return b1, b2, edges


def project_energy(matrix: np.ndarray, vector: np.ndarray) -> tuple[np.ndarray, float]:
    if matrix.size == 0 or matrix.shape[1] == 0:
        return np.zeros_like(vector), 0.0
    coeffs, *_ = np.linalg.lstsq(matrix, vector, rcond=None)
    projection = matrix @ coeffs
    return projection, float(np.dot(projection, projection))


def hodge_component(adjacency: np.ndarray, labels: np.ndarray, include_triangles: bool) -> dict[str, Any]:
    b1, b2, edges = incidence_matrices(adjacency, include_triangles)
    edge_flow = np.array([(labels[j] - labels[i]) * adjacency[i, j] for i, j in edges], dtype=float)
    exact_projection, exact_energy = project_energy(b1.T, edge_flow)
    coexact_projection, coexact_energy = project_energy(b2, edge_flow)
    residual = edge_flow - exact_projection - coexact_projection
    residual_energy = float(np.dot(residual, residual))
    total_energy = exact_energy + coexact_energy + residual_energy
    l1 = b1.T @ b1 + b2 @ b2.T
    hodge_eigenvalues = np.linalg.eigvalsh(l1) if l1.size else np.array([], dtype=float)
    harmonic_dim = int(np.sum(np.abs(hodge_eigenvalues) < 1.0e-8))
    return {
        "edge_count": len(edges),
        "triangle_count": len(triangle_list(adjacency, include_triangles)),
        "exact_fraction": exact_energy / total_energy if total_energy else 0.0,
        "coexact_fraction": coexact_energy / total_energy if total_energy else 0.0,
        "harmonic_fraction": residual_energy / total_energy if total_energy else 0.0,
        "harmonic_dim": harmonic_dim,
        "lowest_l1_eigenvalues": [float(value) for value in hodge_eigenvalues[:8]],
    }


def all_pairs_cost(adjacency: np.ndarray) -> np.ndarray:
    cost = np.zeros_like(adjacency)
    positive = adjacency > 0
    cost[positive] = 1.0 / np.maximum(adjacency[positive], 1.0e-12)
    cost[~positive] = np.inf
    np.fill_diagonal(cost, 0.0)
    return shortest_path(cost, directed=False, unweighted=False)


def lazy_measure(adjacency: np.ndarray, node: int, alpha: float = 0.5) -> tuple[list[int], np.ndarray]:
    neighbors = [idx for idx, weight in enumerate(adjacency[node]) if weight > 0]
    support = [node] + neighbors
    measure = np.zeros(len(support), dtype=float)
    measure[0] = alpha
    weight_sum = float(np.sum(adjacency[node, neighbors])) if neighbors else 0.0
    if weight_sum > 0:
        for idx, neighbor in enumerate(neighbors, start=1):
            measure[idx] = (1.0 - alpha) * adjacency[node, neighbor] / weight_sum
    else:
        measure[0] = 1.0
    return support, measure


def wasserstein_one(left_support: list[int], left: np.ndarray, right_support: list[int], right: np.ndarray, distances: np.ndarray) -> float:
    cost = np.array([[distances[i, j] for j in right_support] for i in left_support], dtype=float)
    objective = cost.reshape(-1)
    row_count, col_count = cost.shape
    equalities = []
    values = []
    for row in range(row_count):
        constraint = np.zeros(row_count * col_count, dtype=float)
        constraint[row * col_count : (row + 1) * col_count] = 1.0
        equalities.append(constraint)
        values.append(left[row])
    for col in range(col_count):
        constraint = np.zeros(row_count * col_count, dtype=float)
        constraint[col::col_count] = 1.0
        equalities.append(constraint)
        values.append(right[col])
    result = linprog(objective, A_eq=np.vstack(equalities), b_eq=np.array(values), bounds=(0.0, None), method="highs")
    return float(result.fun) if result.success else float("nan")


def forman_curvatures(adjacency: np.ndarray) -> list[float]:
    degree = np.sum(adjacency > 0, axis=1)
    values: list[float] = []
    for i, j in edge_list(adjacency):
        common = sum(1 for k in range(adjacency.shape[0]) if adjacency[i, k] > 0 and adjacency[j, k] > 0)
        values.append(float(4.0 - degree[i] - degree[j] + common))
    return values


def ollivier_proxy_curvatures(adjacency: np.ndarray) -> list[float]:
    distances = all_pairs_cost(adjacency)
    values: list[float] = []
    for i, j in edge_list(adjacency):
        left_support, left_measure = lazy_measure(adjacency, i)
        right_support, right_measure = lazy_measure(adjacency, j)
        w1 = wasserstein_one(left_support, left_measure, right_support, right_measure, distances)
        edge_distance = distances[i, j]
        values.append(float(1.0 - w1 / edge_distance) if edge_distance > 1.0e-12 else 0.0)
    return values


def curvature_component(adjacency: np.ndarray) -> dict[str, Any]:
    forman = forman_curvatures(adjacency)
    ollivier = ollivier_proxy_curvatures(adjacency)
    return {
        "forman_mean": mean(forman),
        "forman_std": std(forman),
        "forman_negative_fraction": fraction_negative(forman),
        "ollivier_mean": mean(ollivier),
        "ollivier_std": std(ollivier),
        "ollivier_negative_fraction": fraction_negative(ollivier),
    }


def edge_feature_matrix(adjacency: np.ndarray) -> np.ndarray:
    forman_by_edge = dict(zip(edge_list(adjacency), forman_curvatures(adjacency)))
    ollivier_by_edge = dict(zip(edge_list(adjacency), ollivier_proxy_curvatures(adjacency)))
    rows = []
    for i, j in edge_list(adjacency):
        rows.append([adjacency[i, j], forman_by_edge[(i, j)], ollivier_by_edge[(i, j)]])
    return np.array(rows, dtype=float)


def sinkhorn_distance(left: np.ndarray, right: np.ndarray, epsilon: float) -> float:
    cost = squared_distance_matrix(left, right)
    kernel = np.exp(-cost / max(epsilon, 1.0e-12))
    a = np.full(left.shape[0], 1.0 / left.shape[0])
    b = np.full(right.shape[0], 1.0 / right.shape[0])
    u = np.ones_like(a)
    v = np.ones_like(b)
    for _ in range(200):
        u = a / np.maximum(kernel @ v, 1.0e-300)
        v = b / np.maximum(kernel.T @ u, 1.0e-300)
    plan = (u[:, None] * kernel) * v[None, :]
    return float(np.sum(plan * cost))


def sliced_wasserstein(left: np.ndarray, right: np.ndarray) -> float:
    directions = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
        ],
        dtype=float,
    )
    values = []
    for direction in directions:
        unit = direction / np.linalg.norm(direction)
        left_projection = np.sort(left @ unit)
        right_projection = np.sort(right @ unit)
        count = min(len(left_projection), len(right_projection))
        values.append(float(np.mean(np.abs(left_projection[:count] - right_projection[:count]))))
    return mean(values)


def mmd_rbf(left: np.ndarray, right: np.ndarray, bandwidth: float) -> float:
    gamma = 1.0 / (2.0 * bandwidth * bandwidth)
    k_xx = np.exp(-gamma * squared_distance_matrix(left, left))
    k_yy = np.exp(-gamma * squared_distance_matrix(right, right))
    k_xy = np.exp(-gamma * squared_distance_matrix(left, right))
    return float(np.mean(k_xx) + np.mean(k_yy) - 2.0 * np.mean(k_xy))


def transport_kernel_component(adjacency: np.ndarray, reference_features: np.ndarray, config: BridgeComponentConfig) -> dict[str, Any]:
    features = edge_feature_matrix(adjacency)
    sinkhorn_xy = sinkhorn_distance(reference_features, features, config.sinkhorn_epsilon)
    sinkhorn_xx = sinkhorn_distance(reference_features, reference_features, config.sinkhorn_epsilon)
    sinkhorn_yy = sinkhorn_distance(features, features, config.sinkhorn_epsilon)
    return {
        "sinkhorn_distance": sinkhorn_xy,
        "sinkhorn_divergence": max(0.0, sinkhorn_xy - 0.5 * sinkhorn_xx - 0.5 * sinkhorn_yy),
        "sliced_wasserstein": sliced_wasserstein(reference_features, features),
        "mmd_rbf": mmd_rbf(reference_features, features, config.mmd_bandwidth),
    }


def squared_distance_matrix(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return np.sum((left[:, None, :] - right[None, :, :]) ** 2, axis=2)


def mean(values: list[float]) -> float:
    return float(sum(values) / len(values)) if values else 0.0


def std(values: list[float]) -> float:
    if not values:
        return 0.0
    avg = mean(values)
    return float(math.sqrt(sum((value - avg) ** 2 for value in values) / len(values)))


def fraction_negative(values: list[float]) -> float:
    return float(sum(1 for value in values if value < 0) / len(values)) if values else 0.0


def vector_distance(left: list[float], right: list[float]) -> float:
    a = np.array(left, dtype=float)
    b = np.array(right, dtype=float)
    denominator = max(float(np.linalg.norm(a)), 1.0e-12)
    return float(np.linalg.norm(a - b) / denominator)


def spectral_distance(reference: dict[str, Any], observed: dict[str, Any]) -> float:
    return vector_distance(reference["eigenvalues"], observed["eigenvalues"])


def hodge_distance(reference: dict[str, Any], observed: dict[str, Any]) -> float:
    fields = ["exact_fraction", "coexact_fraction", "harmonic_fraction", "harmonic_dim", "triangle_count"]
    return vector_distance([reference[field] for field in fields], [observed[field] for field in fields])


def curvature_distance(reference: dict[str, Any], observed: dict[str, Any]) -> float:
    fields = [
        "forman_mean",
        "forman_std",
        "forman_negative_fraction",
        "ollivier_mean",
        "ollivier_std",
        "ollivier_negative_fraction",
    ]
    return vector_distance([reference[field] for field in fields], [observed[field] for field in fields])


def transport_kernel_distance(observed: dict[str, Any]) -> float:
    return float(observed["sinkhorn_divergence"] + observed["sliced_wasserstein"] + observed["mmd_rbf"])


def expected_response(condition: str) -> str:
    return {
        "known_positive": "all component distances near zero",
        "label_permutation_control": "all component distances near zero if labels do not leak",
        "weight_shuffle_control": "transport/kernel and weighted spectral diagnostics should move",
        "topology_destroyed_control": "spectral, curvature, and transport diagnostics should move",
        "hodge_triangle_removed_control": "Hodge sector diagnostics should move while graph-only components remain stable",
    }[condition]


def observed_status(condition: str, distances: dict[str, float], config: BridgeComponentConfig) -> str:
    if condition in {"known_positive", "label_permutation_control"}:
        if max(distances.values()) <= config.label_invariance_max:
            return "PASS"
        return "FAIL_LABEL_OR_SELF_INVARIANCE"
    if condition == "weight_shuffle_control":
        if (
            distances["transport_kernel"] >= config.transport_distance_min
            and distances["spectral"] >= config.spectral_distance_min * 0.25
        ):
            return "PASS"
        return "FAIL_WEIGHT_CONTROL"
    if condition == "topology_destroyed_control":
        if (
            distances["spectral"] >= config.spectral_distance_min
            and distances["curvature"] >= config.curvature_distance_min
            and distances["transport_kernel"] >= config.transport_distance_min
        ):
            return "PASS"
        return "FAIL_TOPOLOGY_CONTROL"
    if condition == "hodge_triangle_removed_control":
        if (
            distances["hodge"] >= config.hodge_distance_min
            and distances["spectral"] <= config.label_invariance_max
            and distances["curvature"] <= config.label_invariance_max
        ):
            return "PASS"
        return "FAIL_HODGE_CONTROL"
    raise ValueError(f"unknown condition {condition}")


def component_rows(condition: str, components: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for component, metrics in components.items():
        for metric, value in metrics.items():
            if isinstance(value, list):
                for index, item in enumerate(value):
                    rows.append({"condition": condition, "component": component, "metric": f"{metric}_{index}", "value": f"{item:.12g}"})
            else:
                rows.append({"condition": condition, "component": component, "metric": metric, "value": f"{float(value):.12g}"})
    return rows


def source_manifest() -> dict[str, Any]:
    return {
        "source": "literature_component_bridge synthetic sidecar",
        "literature_extraction": repo_rel(REPO_ROOT / "papers/bridge_literature/extraction/haos_iip_bridge_extraction.md"),
        "methods_spine": repo_rel(REPO_ROOT / "papers/bridge_literature/extraction/bridge_methods_spine.md"),
        "pdf_manifest": repo_rel(REPO_ROOT / "papers/bridge_literature/extraction/pdf_extraction_manifest.csv"),
        "code": repo_rel(ROOT / "run_literature_component_bridge.py"),
    }


def precommitment_contract(config: BridgeComponentConfig) -> dict[str, Any]:
    return {
        "name": "Literature Component Bridge v0.1",
        "status": "PRECOMMITTED_COMPONENT_CALIBRATION",
        "purpose": "Combine literature-extracted spectral, Hodge, curvature, transport, and kernel components in a bounded HAOS-IIP bridge instrument.",
        "source_artifacts": source_manifest(),
        "components": [
            "spectral_identity",
            "hodge_sector_split",
            "graph_curvature",
            "transport_metrics",
            "kernel_two_sample",
        ],
        "conditions": {condition: expected_response(condition) for condition in CONDITIONS},
        "thresholds": asdict(config),
        "mapping_status": {
            "spectral_identity": "HEURISTIC",
            "hodge_sector_split": "HEURISTIC",
            "graph_curvature": "HEURISTIC",
            "transport_metrics": "HEURISTIC",
            "kernel_two_sample": "HEURISTIC",
        },
        "claim_ceiling": "OPERATIONAL_MAPPING",
        "non_claims": [
            "not a physical bridge",
            "not a continuum limit",
            "not a quantum, gravity, or field-theory derivation",
            "not empirical validation",
            "not a change to frozen HAOS-IIP phases",
        ],
    }


def run_bridge(config: BridgeComponentConfig, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    reference_adjacency, reference_labels, reference_triangles = condition_graph("known_positive")
    reference_features = edge_feature_matrix(reference_adjacency)
    reference_components = {
        "spectral": spectral_component(reference_adjacency, reference_labels),
        "hodge": hodge_component(reference_adjacency, reference_labels, reference_triangles),
        "curvature": curvature_component(reference_adjacency),
        "transport_kernel": transport_kernel_component(reference_adjacency, reference_features, config),
    }
    score_rows: list[dict[str, Any]] = []
    control_rows: list[dict[str, Any]] = []
    condition_results: dict[str, Any] = {}

    for condition in CONDITIONS:
        adjacency, labels, include_triangles = condition_graph(condition)
        components = {
            "spectral": spectral_component(adjacency, labels),
            "hodge": hodge_component(adjacency, labels, include_triangles),
            "curvature": curvature_component(adjacency),
            "transport_kernel": transport_kernel_component(adjacency, reference_features, config),
        }
        distances = {
            "spectral": spectral_distance(reference_components["spectral"], components["spectral"]),
            "hodge": hodge_distance(reference_components["hodge"], components["hodge"]),
            "curvature": curvature_distance(reference_components["curvature"], components["curvature"]),
            "transport_kernel": transport_kernel_distance(components["transport_kernel"]),
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
                "spectral_distance": f"{distances['spectral']:.12g}",
                "hodge_distance": f"{distances['hodge']:.12g}",
                "curvature_distance": f"{distances['curvature']:.12g}",
                "transport_kernel_distance": f"{distances['transport_kernel']:.12g}",
                "observed_status": status,
            }
        )

    control_pass = all(result["observed_status"] == "PASS" for result in condition_results.values())
    labels = [
        "BRIDGE_COMPONENTS_BUILT",
        "SPECTRAL_COMPONENT_AVAILABLE",
        "HODGE_COMPONENT_AVAILABLE",
        "CURVATURE_COMPONENT_AVAILABLE",
        "TRANSPORT_KERNEL_COMPONENT_AVAILABLE",
        "PHYSICAL_BRIDGE_NOT_ESTABLISHED",
        "CLAIM_GATED_OPERATIONAL_MAPPING",
    ]
    labels.append("COMPONENT_CONTROLS_PASS" if control_pass else "COMPONENT_CONTROLS_PARTIAL")
    if condition_results["label_permutation_control"]["observed_status"] == "PASS":
        labels.append("LABEL_INVARIANCE_PASS")

    result_payload = {
        "bridge_id": "literature_component_bridge_v0_1",
        "status": "PASS" if control_pass else "MIXED_OPEN",
        "classification": "OPERATIONAL_MAPPING",
        "labels": labels,
        "condition_results": condition_results,
        "claim_ceiling": "No physical bridge claim; component calibration only.",
    }
    result_payload["result_hash"] = stable_hash("literature_bridge", result_payload)

    contract = precommitment_contract(config)
    manifest = source_manifest()
    write_json(output_dir / PRECOMMITMENT_PATH.name, contract)
    write_json(output_dir / SOURCE_MANIFEST_PATH.name, manifest)
    write_csv(output_dir / COMPONENT_SCORES_PATH.name, score_rows, COMPONENT_SCORE_FIELDS)
    write_csv(output_dir / CONTROL_RESULTS_PATH.name, control_rows, CONTROL_RESULT_FIELDS)
    write_json(output_dir / RESULT_PATH.name, result_payload)
    write_report(output_dir / REPORT_PATH.name, result_payload, control_rows)
    return result_payload


def write_report(path: Path, result: dict[str, Any], control_rows: list[dict[str, Any]]) -> None:
    lines = [
        "# Literature Component Bridge Report",
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
        "| Condition | Spectral | Hodge | Curvature | Transport/kernel | Status |",
        "| --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in control_rows:
        lines.append(
            "| {condition} | {spectral_distance} | {hodge_distance} | {curvature_distance} | {transport_kernel_distance} | `{observed_status}` |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "This sidecar builds bridge instrumentation from the literature-extracted components.",
            "It does not establish a physical bridge, continuum limit, quantum result, gravity result, or field-theory derivation.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    result = run_bridge(BridgeComponentConfig(), args.output_dir)
    print(json.dumps({"status": result["status"], "result_hash": result["result_hash"]}, sort_keys=True))


if __name__ == "__main__":
    main()
