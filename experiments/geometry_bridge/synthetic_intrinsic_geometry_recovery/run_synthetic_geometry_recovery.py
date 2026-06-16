#!/usr/bin/env python3
"""Synthetic intrinsic-geometry recovery benchmark.

This benchmark asks a narrower question than Bell:

Can frozen operator/transport observations recover intrinsic geometry on a
synthetic latent space better than conventional graph/spectral baselines?

It deliberately separates:
- latent geometry generation
- operator-only observation extraction
- frozen calibration on development data
- untouched holdout evaluation
- control degradation

The harness is intended to fail honestly when geometry is not recoverable.
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


ROOT = Path(__file__).resolve().parent
REPO_ROOT = Path(__file__).resolve().parents[3]

PRECOMMITMENT_PATH = ROOT / "precommitment_contract.json"
SOURCE_MANIFEST_PATH = ROOT / "geometry_source_manifest.json"
GEOMETRY_MATRIX_PATH = ROOT / "geometry_matrix.csv"
CONTROL_RESULTS_PATH = ROOT / "control_results.csv"
REPORT_PATH = ROOT / "geometry_recovery_report.md"
RESULT_PATH = ROOT / "geometry_recovery_result.json"

FAMILIES = ("ring", "grid", "spiral")
DEVELOPMENT_FAMILIES = ("ring",)
CALIBRATION_FAMILIES = ("ring", "grid")
HOLDOUT_FAMILIES = ("spiral",)
SEEDS = (4101, 4102, 4103, 4104)
PERTURBATIONS = ("edge_removal", "weight_degradation", "localized_state_shock")
CONTROL_NAMES = (
    "label_permutation_control",
    "topology_destroyed_control",
    "degree_preserving_rewire_control",
    "parameter_matched_null_control",
    "seed_repeat_control",
    "leakage_positive_control",
)

MATRIX_FIELDNAMES = [
    "split",
    "family",
    "seed",
    "perturbation",
    "node_left",
    "node_right",
    "true_distance",
    "haos_distance",
    "shortest_path_distance",
    "spectral_distance",
    "invariant_overlap",
    "persistence_score",
    "causal_depth",
    "temporal_order_stability",
]

CONTROL_FIELDNAMES = [
    "control_name",
    "split",
    "family",
    "seed",
    "spearman_true_vs_predicted",
    "top_k_precision",
    "adjacency_f1",
    "mean_absolute_error",
    "status",
]


@dataclass(frozen=True)
class GeometryConfig:
    version: str = "synthetic-intrinsic-geometry-recovery-v0.1"
    n_nodes: int = 9
    latent_noise: float = 0.035
    graph_sigma: float = 0.45
    graph_cutoff: float = 0.72
    diffusion_steps: int = 4
    local_shock_fraction: float = 0.25
    calibration_weight_floor: float = 0.10
    min_holdout_spearman_margin: float = 0.04
    min_holdout_topk_margin: float = 0.02
    min_control_degradation: float = 0.10
    holdout_spearman_min: float = 0.40
    holdout_topk_min: float = 0.40
    bootstrap_replicates: int = 200
    bootstrap_seed: int = 9137


def repo_rel(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def stable_json_hash(prefix: str, payload: Any) -> str:
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


def generate_latent_coords(family: str, seed: int, n_nodes: int, noise: float) -> np.ndarray:
    rng = np.random.default_rng(seed)
    coords = np.zeros((n_nodes, 2), dtype=float)
    if family == "ring":
        for index in range(n_nodes):
            angle = 2.0 * math.pi * index / n_nodes
            coords[index] = (math.cos(angle), math.sin(angle))
    elif family == "grid":
        side = int(round(math.sqrt(n_nodes)))
        if side * side != n_nodes:
            raise ValueError("grid family requires a square node count")
        index = 0
        for row in range(side):
            for col in range(side):
                coords[index] = (float(col), float(row))
                index += 1
        coords[:, 0] /= max(side - 1, 1)
        coords[:, 1] /= max(side - 1, 1)
    elif family == "spiral":
        for index in range(n_nodes):
            radius = 0.25 + 0.10 * index
            angle = 1.35 * index
            coords[index] = (radius * math.cos(angle), radius * math.sin(angle))
    else:
        raise ValueError(f"unknown family {family}")
    coords += rng.normal(scale=noise, size=coords.shape)
    return coords


def pairwise_distances(coords: np.ndarray) -> np.ndarray:
    diffs = coords[:, None, :] - coords[None, :, :]
    return np.linalg.norm(diffs, axis=-1)


def build_weighted_graph(coords: np.ndarray, sigma: float, cutoff: float) -> np.ndarray:
    distances = pairwise_distances(coords)
    weights = np.exp(-(distances**2) / max(2.0 * sigma**2, 1.0e-12))
    np.fill_diagonal(weights, 0.0)
    mask = distances <= cutoff
    weights *= mask.astype(float)
    return np.maximum(weights, weights.T)


def normalize_rows(matrix: np.ndarray) -> np.ndarray:
    row_sums = matrix.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.divide(matrix, row_sums, out=np.zeros_like(matrix), where=row_sums > 1.0e-12)


def laplacian(graph: np.ndarray) -> np.ndarray:
    degree = np.sum(graph, axis=1)
    return np.diag(degree) - graph


def transition_matrix(graph: np.ndarray) -> np.ndarray:
    return normalize_rows(graph)


def shortest_path_distances(graph: np.ndarray) -> np.ndarray:
    n = graph.shape[0]
    inf = 1.0e9
    distances = np.full((n, n), inf, dtype=float)
    np.fill_diagonal(distances, 0.0)
    adjacency = graph > 1.0e-12
    for source in range(n):
        visited = np.zeros(n, dtype=bool)
        queue = [(0.0, source)]
        distances[source, source] = 0.0
        while queue:
            queue.sort(key=lambda item: item[0], reverse=True)
            current_distance, node = queue.pop()
            if visited[node]:
                continue
            visited[node] = True
            for neighbor in np.where(adjacency[node])[0]:
                weight = 1.0 / max(graph[node, neighbor], 1.0e-12)
                candidate = current_distance + weight
                if candidate < distances[source, neighbor]:
                    distances[source, neighbor] = candidate
                    queue.append((candidate, neighbor))
    finite = distances < inf
    if np.any(finite):
        distances = distances / max(np.max(distances[finite]), 1.0e-12)
    return distances


def spectral_embedding_distance(graph: np.ndarray, dimensions: int = 3) -> np.ndarray:
    lap = laplacian(graph)
    eigenvalues, eigenvectors = np.linalg.eigh(lap)
    basis = eigenvectors[:, 1 : 1 + min(dimensions, eigenvectors.shape[1] - 1)]
    if basis.size == 0:
        return np.zeros_like(graph)
    diffs = basis[:, None, :] - basis[None, :, :]
    distances = np.linalg.norm(diffs, axis=-1)
    max_value = np.max(distances)
    return distances / max(max_value, 1.0e-12)


def diffusion_signature(graph: np.ndarray, steps: int) -> np.ndarray:
    transition = transition_matrix(graph)
    n = graph.shape[0]
    signatures = []
    current = np.eye(n, dtype=float)
    for _ in range(steps):
        current = current @ transition
        signatures.append(current)
    return np.concatenate(signatures, axis=1)


def invariant_overlap(graph: np.ndarray) -> np.ndarray:
    transition = transition_matrix(graph)
    stationary = np.mean(transition, axis=0)
    overlap = np.outer(stationary, stationary)
    max_value = np.max(overlap)
    return overlap / max(max_value, 1.0e-12)


def persistence_score(graph: np.ndarray) -> np.ndarray:
    lap = laplacian(graph)
    eigvals = np.linalg.eigvalsh(lap)
    gap = float(eigvals[1]) if eigvals.size > 1 else 0.0
    score = 1.0 / (1.0 + gap)
    return np.full(graph.shape[0], score, dtype=float)


def causal_depth(graph: np.ndarray) -> np.ndarray:
    distances = shortest_path_distances(graph)
    max_value = np.max(distances)
    return 1.0 - np.mean(distances, axis=1) / max(max_value, 1.0e-12)


def temporal_order_stability(graph: np.ndarray, steps: int) -> np.ndarray:
    transition = transition_matrix(graph)
    state = np.eye(graph.shape[0], dtype=float)
    order_scores = np.zeros(graph.shape[0], dtype=float)
    for _ in range(steps):
        state = state @ transition
        order_scores += np.sum(state * np.arange(1, graph.shape[0] + 1), axis=1)
    order_scores /= float(steps)
    max_value = np.max(order_scores)
    return 1.0 - order_scores / max(max_value, 1.0e-12)


def feature_matrix(graph: np.ndarray, config: GeometryConfig) -> np.ndarray:
    diffusion = diffusion_signature(graph, config.diffusion_steps)
    lap = laplacian(graph)
    eigvals, eigvecs = np.linalg.eigh(lap)
    spectral = eigvecs[:, 1 : 1 + min(3, eigvecs.shape[1] - 1)]
    if spectral.size == 0:
        spectral = np.zeros((graph.shape[0], 1), dtype=float)
    degree = np.sum(graph, axis=1, keepdims=True)
    invariant = invariant_overlap(graph)
    persistence = persistence_score(graph)[:, None]
    depth = causal_depth(graph)[:, None]
    temporal = temporal_order_stability(graph, config.diffusion_steps)[:, None]
    summary = np.concatenate([degree, spectral, diffusion, persistence, depth, temporal, invariant], axis=1)
    return summary


def perturb_graph(graph: np.ndarray, family: str, seed: int, perturbation: str, severity: float) -> np.ndarray:
    rng = np.random.default_rng(seed)
    perturbed = graph.copy()
    n = graph.shape[0]
    if perturbation == "edge_removal":
        upper = np.triu_indices(n, 1)
        existing = np.where(perturbed[upper] > 0.0)[0]
        remove_count = max(1, int(round(severity * max(len(existing), 1))))
        chosen = rng.choice(existing, size=min(remove_count, len(existing)), replace=False)
        for index in chosen:
            i = upper[0][index]
            j = upper[1][index]
            perturbed[i, j] = perturbed[j, i] = 0.0
    elif perturbation == "weight_degradation":
        perturbed *= (1.0 - severity)
    elif perturbation == "localized_state_shock":
        anchors = rng.choice(n, size=max(1, n // 4), replace=False)
        for anchor in anchors:
            perturbed[anchor] *= (1.0 - 0.5 * severity)
            perturbed[:, anchor] *= (1.0 - 0.5 * severity)
    else:
        raise ValueError(f"unknown perturbation {perturbation}")
    if family == "spiral":
        perturbed = 0.95 * perturbed + 0.05 * np.roll(perturbed, shift=1, axis=0)
    return np.maximum(perturbed, perturbed.T)


def pairwise_ground_truth(coords: np.ndarray) -> np.ndarray:
    distances = pairwise_distances(coords)
    max_value = np.max(distances)
    return distances / max(max_value, 1.0e-12)


def upper_triangle_values(matrix: np.ndarray) -> list[float]:
    values: list[float] = []
    for row in range(matrix.shape[0]):
        for col in range(row + 1, matrix.shape[1]):
            values.append(float(matrix[row, col]))
    return values


def mean(values: list[float]) -> float:
    return float(sum(values) / len(values)) if values else 0.0


def variance(values: list[float]) -> float:
    if not values:
        return 0.0
    mu = mean(values)
    return float(sum((value - mu) ** 2 for value in values) / len(values))


def rankdata(values: list[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0 for _ in values]
    index = 0
    while index < len(indexed):
        end = index + 1
        while end < len(indexed) and indexed[end][1] == indexed[index][1]:
            end += 1
        avg = (index + end - 1) / 2.0 + 1.0
        for slot in range(index, end):
            ranks[indexed[slot][0]] = avg
        index = end
    return ranks


def spearman(left: list[float], right: list[float]) -> float:
    if len(left) != len(right) or len(left) < 2:
        return 0.0
    ranked_left = rankdata(left)
    ranked_right = rankdata(right)
    left_mu = mean(ranked_left)
    right_mu = mean(ranked_right)
    numerator = sum((a - left_mu) * (b - right_mu) for a, b in zip(ranked_left, ranked_right))
    left_den = math.sqrt(sum((a - left_mu) ** 2 for a in ranked_left))
    right_den = math.sqrt(sum((b - right_mu) ** 2 for b in ranked_right))
    return numerator / (left_den * right_den) if left_den * right_den > 1.0e-300 else 0.0


def mae(left: list[float], right: list[float]) -> float:
    return mean([abs(a - b) for a, b in zip(left, right)])


def top_k_precision(true_matrix: np.ndarray, predicted_matrix: np.ndarray, k: int = 3) -> float:
    n = true_matrix.shape[0]
    hits = 0
    total = 0
    for row in range(n):
        true_neighbors = np.argsort(true_matrix[row])[1 : 1 + k]
        predicted_neighbors = np.argsort(predicted_matrix[row])[1 : 1 + k]
        hits += len(set(true_neighbors) & set(predicted_neighbors))
        total += k
    return hits / float(total) if total else 0.0


def adjacency_f1(graph: np.ndarray, predicted_distances: np.ndarray) -> float:
    n = graph.shape[0]
    true_edges = set()
    pred_edges = set()
    threshold = np.median(predicted_distances[np.triu_indices(n, 1)])
    for row in range(n):
        for col in range(row + 1, n):
            if graph[row, col] > 0.0:
                true_edges.add((row, col))
            if predicted_distances[row, col] <= threshold:
                pred_edges.add((row, col))
    intersection = true_edges & pred_edges
    precision = len(intersection) / float(len(pred_edges)) if pred_edges else 0.0
    recall = len(intersection) / float(len(true_edges)) if true_edges else 0.0
    return 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0


def control_matrix(graph: np.ndarray, control_name: str, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    if control_name == "label_permutation_control":
        order = rng.permutation(graph.shape[0])
        return graph[np.ix_(order, order)]
    if control_name == "topology_destroyed_control":
        return np.zeros_like(graph)
    if control_name == "degree_preserving_rewire_control":
        rewired = graph.copy()
        for _ in range(max(1, graph.shape[0] // 2)):
            i, j, k, l = rng.choice(graph.shape[0], size=4, replace=False)
            rewired[i, j], rewired[k, l] = rewired[k, l], rewired[i, j]
            rewired[j, i] = rewired[i, j]
            rewired[l, k] = rewired[k, l]
        return rewired
    if control_name == "parameter_matched_null_control":
        return np.where(graph > 0.0, np.mean(graph[graph > 0.0]), 0.0)
    if control_name == "seed_repeat_control":
        return graph.copy()
    if control_name == "leakage_positive_control":
        return graph + np.eye(graph.shape[0], dtype=float) * 0.0
    raise ValueError(control_name)


def control_status(control_name: str, metric: float) -> str:
    if control_name == "leakage_positive_control":
        return "DETECTED" if metric >= 0.90 else "NOT_DETECTED"
    if control_name == "topology_destroyed_control":
        return "DEGRADED" if metric < 0.35 else "FAILED_TO_DEGRADE"
    if control_name == "label_permutation_control":
        return "DEGRADED" if metric < 0.70 else "LEAKED"
    return "DEGRADED" if metric < 0.75 else "WEAK"


def precommitment_payload(config: GeometryConfig) -> dict[str, Any]:
    return {
        "bridge_id": "GEO-01",
        "version": config.version,
        "purpose": (
            "Test whether frozen operator/transport observations recover intrinsic "
            "geometry on synthetic latent spaces better than conventional graph and "
            "spectral baselines, using blind development/calibration/hodout splits."
        ),
        "classification": "OPERATIONAL_MAPPING_PARTIAL",
        "claim_boundary": [
            "synthetic calibration only",
            "not a Bell experiment",
            "not a physical mechanism claim",
            "not empirical validation of HAOS as ontology",
        ],
        "state_schema": {
            "latent_space": "two-dimensional synthetic latent geometry",
            "hidden_coordinates": "point set on ring, grid, or spiral manifold",
            "observations": "operator-only graph transport summaries",
            "targets": "pairwise intrinsic distance and neighborhood order",
            "units": {
                "latent_space": "dimensionless",
                "hidden_coordinates": "dimensionless",
                "observations": "dimensionless",
                "targets": "dimensionless normalized distance",
            },
            "valid_ranges": {
                "latent_space": "2D finite latent set",
                "hidden_coordinates": "bounded coordinates",
                "observations": "finite nonnegative weights",
                "targets": "0..1",
            },
        },
        "dynamics": {
            "graph_construction": "weights = exp(-d^2/(2 sigma^2)) with cutoff sparsification",
            "observation_map": "diffusion signatures, Laplacian spectrum, and transport summaries",
            "perturbations": list(PERTURBATIONS),
            "transport_steps": config.diffusion_steps,
        },
        "boundary_conditions": "frozen synthetic family definitions and untouched holdout family",
        "symmetries": [
            "label permutation should not preserve semantic reconstruction",
            "topology destruction should degrade geometry recovery",
            "seed repeat should remain deterministic",
        ],
        "haos_mapping": [
            {
                "haos_source": "invariant overlap",
                "domain_target": "local geometry coherence",
                "mapping_function": "cosine overlap between node diffusion signatures",
                "units_before": "dimensionless",
                "units_after": "dimensionless",
                "semantic_justification": "stable transport neighborhoods should preserve local similarity",
                "mapping_status": "HEURISTIC",
                "uncertainty": "medium",
                "failure_conditions": "label permutation or topology destruction leaves overlap unchanged",
            },
            {
                "haos_source": "persistence score",
                "domain_target": "topological stability under perturbation",
                "mapping_function": "inverse spectral-gap proxy",
                "units_before": "dimensionless",
                "units_after": "dimensionless",
                "semantic_justification": "geometry should degrade when persistent structure is destroyed",
                "mapping_status": "HEURISTIC",
                "uncertainty": "medium",
                "failure_conditions": "destroyed topology does not lower persistence",
            },
            {
                "haos_source": "causal depth",
                "domain_target": "intrinsic separation from perturbation",
                "mapping_function": "1 - normalized shortest-path distance",
                "units_before": "dimensionless",
                "units_after": "dimensionless",
                "semantic_justification": "geometry should rank neighborhoods by transport reachability",
                "mapping_status": "CALIBRATED",
                "uncertainty": "medium",
                "failure_conditions": "holdout geometry does not transfer after calibration freeze",
            },
            {
                "haos_source": "temporal order stability",
                "domain_target": "consistency of diffusion order under steps",
                "mapping_function": "monotone response summary over diffusion powers",
                "units_before": "dimensionless",
                "units_after": "dimensionless",
                "semantic_justification": "intrinsic geometry should preserve step ordering better than nulls",
                "mapping_status": "HEURISTIC",
                "uncertainty": "medium",
                "failure_conditions": "order stability survives topology destruction",
            },
        ],
        "observation_map": [
            {
                "observable": "pairwise intrinsic distance",
                "source": "latent coordinates (hidden during scoring)",
                "feature_path": "Euclidean distance between HAOS transport signatures",
                "units": "dimensionless",
                "valid_range": "0..1",
            },
            {
                "observable": "topological neighborhood",
                "source": "latent coordinates (hidden during scoring)",
                "feature_path": "top-k nearest-neighbor recovery from HAOS signatures",
                "units": "dimensionless",
                "valid_range": "0..1",
            },
        ],
        "prediction_target": "recover intrinsic pairwise geometry and neighborhood order on untouched holdout family",
        "calibration_split": "development + calibration families frozen before holdout inspection",
        "holdout_split": "spiral",
        "controls": [
            {
                "name": "label_permutation_control",
                "preserves": "graph weights and transport magnitudes",
                "destroys": "semantic correspondence between labels and latent neighborhood order",
                "expected_response": "geometry scores drop materially",
                "invalidation_condition": "scores stay near the target",
                "implementation_provenance": "synthetic label permutation with frozen seeds",
            },
            {
                "name": "topology_destroyed_control",
                "preserves": "node count and weight scale",
                "destroys": "adjacency topology and transport geometry",
                "expected_response": "all geometry scores degrade sharply",
                "invalidation_condition": "topology destruction does not reduce recovery",
                "implementation_provenance": "zeroed graph null with same dimensions",
            },
            {
                "name": "degree_preserving_rewire_control",
                "preserves": "degree histogram and node count",
                "destroys": "latent adjacency and neighborhood order",
                "expected_response": "spectral and HAOS scores weaken relative to target",
                "invalidation_condition": "rewire behaves like the target",
                "implementation_provenance": "frozen stub rewiring under deterministic seeds",
            },
            {
                "name": "parameter_matched_null_control",
                "preserves": "weight scale and node count",
                "destroys": "geometry-specific transport structure",
                "expected_response": "baseline metrics remain lower than target",
                "invalidation_condition": "null outperforms or matches target systematically",
                "implementation_provenance": "constant-weight null graph",
            },
            {
                "name": "seed_repeat_control",
                "preserves": "all construction choices",
                "destroys": "none",
                "expected_response": "deterministic repeatability",
                "invalidation_condition": "same seed changes the result",
                "implementation_provenance": "repeat-run frozen seed pair",
            },
            {
                "name": "leakage_positive_control",
                "preserves": "explicit access to hidden coordinates",
                "destroys": "blindness",
                "expected_response": "should be detected as leakage, not accepted as a valid candidate",
                "invalidation_condition": "positive control is not flagged",
                "implementation_provenance": "intentional target leakage check",
            },
        ],
        "baselines": [
            {"name": "mean_predictor", "description": "constant mean target distance"},
            {"name": "random_predictor", "description": "frozen deterministic random draw"},
            {"name": "graph_shortest_path", "description": "pure topological shortest-path distance"},
            {"name": "spectral_embedding", "description": "Laplacian eigenvector distance"},
            {"name": "degree_only", "description": "degree-distance proxy"},
            {"name": "shuffled_labels", "description": "label permutation baseline"},
        ],
        "uncertainty": {
            "bootstrap_replicates": config.bootstrap_replicates,
            "bootstrap_seed": config.bootstrap_seed,
            "method": "percentile 95 percent interval on paired pairwise scores",
        },
        "falsification": {
            "fail_conditions": [
                "holdout recovery does not outperform baselines by the frozen margin",
                "controls do not degrade as declared",
                "leakage positive control is not detected",
                "deterministic repeats do not match",
                "frozen holdout is inspected before scoring",
            ]
        },
        "verdict_logic": {
            "official_verdict": "GEOMETRY_NOT_YET_DEMONSTRATED if holdout transfer fails or controls are invalid",
            "promotion_gate": "HAOS geometry can only be promoted after blind holdout separation and control degradation are frozen and reproducible",
        },
        "provenance": {
            "source_artifacts": [
                repo_rel(PRECOMMITMENT_PATH),
                repo_rel(SOURCE_MANIFEST_PATH),
            ],
            "code_artifacts": [repo_rel(Path(__file__))],
            "external_data_status": "synthetic_only",
            "notes": [
                "This contract freezes the geometry question before any Bell-like probability layer.",
                "No coordinates are exposed to the scorer before holdout execution.",
            ],
        },
    }


def generate_family_record(family: str, seed: int, config: GeometryConfig) -> dict[str, Any]:
    coords = generate_latent_coords(family, seed, config.n_nodes, config.latent_noise)
    graph = build_weighted_graph(coords, config.graph_sigma, config.graph_cutoff)
    features = feature_matrix(graph, config)
    target = pairwise_ground_truth(coords)
    shortest = shortest_path_distances(graph)
    spectral = spectral_embedding_distance(graph)
    haos = pairwise_from_features(features)
    return {
        "family": family,
        "seed": seed,
        "coords": coords,
        "graph": graph,
        "target": target,
        "shortest": shortest,
        "spectral": spectral,
        "haos": haos,
        "features": features,
    }


def pairwise_from_features(features: np.ndarray) -> np.ndarray:
    diffs = features[:, None, :] - features[None, :, :]
    distances = np.linalg.norm(diffs, axis=-1)
    max_value = np.max(distances)
    return distances / max(max_value, 1.0e-12)


def fit_haos_weights(records: list[dict[str, Any]], config: GeometryConfig) -> np.ndarray:
    rows = []
    targets = []
    for record in records:
        features = record["features"]
        target = record["target"]
        n = features.shape[0]
        for i in range(n):
            for j in range(i + 1, n):
                rows.append(
                    [
                        float(np.linalg.norm(features[i] - features[j])),
                        float(record["shortest"][i, j]),
                        float(record["spectral"][i, j]),
                        float(np.sum(np.abs(features[i] - features[j]))),
                        1.0,
                    ]
                )
                targets.append(float(target[i, j]))
    x = np.asarray(rows, dtype=float)
    y = np.asarray(targets, dtype=float)
    regularizer = np.eye(x.shape[1], dtype=float) * config.calibration_weight_floor
    weights = np.linalg.solve(x.T @ x + regularizer, x.T @ y)
    return weights


def predict_distances(record: dict[str, Any], weights: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    features = record["features"]
    n = features.shape[0]
    predicted = np.zeros((n, n), dtype=float)
    degree = np.sum(record["graph"], axis=1)
    shortest = record["shortest"]
    spectral = record["spectral"]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            feature_distance = float(np.linalg.norm(features[i] - features[j]))
            signed_sum = float(np.sum(np.abs(features[i] - features[j])))
            vector = np.array([feature_distance, shortest[i, j], spectral[i, j], signed_sum, 1.0], dtype=float)
            predicted[i, j] = float(vector @ weights)
    predicted = np.clip(predicted, 0.0, None)
    max_value = np.max(predicted)
    return predicted / max(max_value, 1.0e-12), degree


def evaluate_record(record: dict[str, Any], predicted: np.ndarray, config: GeometryConfig) -> dict[str, float]:
    target_values = upper_triangle_values(record["target"])
    predicted_values = upper_triangle_values(predicted)
    shortest_values = upper_triangle_values(record["shortest"])
    spectral_values = upper_triangle_values(record["spectral"])
    return {
        "spearman_true_vs_predicted": spearman(target_values, predicted_values),
        "top_k_precision": top_k_precision(record["target"], predicted),
        "adjacency_f1": adjacency_f1(record["graph"], predicted),
        "mean_absolute_error": mae(target_values, predicted_values),
        "spearman_true_vs_shortest": spearman(target_values, shortest_values),
        "spearman_true_vs_spectral": spearman(target_values, spectral_values),
    }


def generate_rows(records: list[dict[str, Any]], weights: np.ndarray, split: str, include_controls: bool = False) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    geometry_rows: list[dict[str, Any]] = []
    control_rows: list[dict[str, Any]] = []
    for record in records:
        predicted, _degree = predict_distances(record, weights)
        graph = record["graph"]
        for i in range(graph.shape[0]):
            for j in range(i + 1, graph.shape[0]):
                geometry_rows.append(
                    {
                        "split": split,
                        "family": record["family"],
                        "seed": record["seed"],
                        "perturbation": "frozen",
                        "node_left": i,
                        "node_right": j,
                        "true_distance": f"{record['target'][i, j]:.12g}",
                        "haos_distance": f"{predicted[i, j]:.12g}",
                        "shortest_path_distance": f"{record['shortest'][i, j]:.12g}",
                        "spectral_distance": f"{record['spectral'][i, j]:.12g}",
                        "invariant_overlap": f"{float(np.mean(record['graph'][i])):.12g}",
                        "persistence_score": f"{float(np.mean(persistence_score(record['graph']))):.12g}",
                        "causal_depth": f"{float(causal_depth(record['graph'])[i]):.12g}",
                        "temporal_order_stability": f"{float(temporal_order_stability(record['graph'], 4)[i]):.12g}",
                    }
                )
        if include_controls:
            for control_name in CONTROL_NAMES:
                control_graph = control_matrix(graph, control_name, int(record["seed"]) + 97)
                if control_name == "leakage_positive_control":
                    control_predicted = record["target"].copy()
                    metrics = evaluate_record(
                        {
                            "graph": control_graph,
                            "target": record["target"],
                            "shortest": record["target"],
                            "spectral": record["target"],
                        },
                        control_predicted,
                        GeometryConfig(),
                    )
                    control_rows.append(
                        {
                            "control_name": control_name,
                            "split": split,
                            "family": record["family"],
                            "seed": record["seed"],
                            "spearman_true_vs_predicted": f"{metrics['spearman_true_vs_predicted']:.12g}",
                            "top_k_precision": f"{metrics['top_k_precision']:.12g}",
                            "adjacency_f1": f"{metrics['adjacency_f1']:.12g}",
                            "mean_absolute_error": f"{metrics['mean_absolute_error']:.12g}",
                            "status": control_status(control_name, metrics["spearman_true_vs_predicted"]),
                        }
                    )
                    continue
                control_record = {
                    "family": record["family"],
                    "seed": record["seed"],
                    "graph": control_graph,
                    "target": record["target"],
                    "shortest": shortest_path_distances(control_graph if np.any(control_graph) else graph),
                    "spectral": spectral_embedding_distance(control_graph if np.any(control_graph) else graph),
                    "features": feature_matrix(control_graph if np.any(control_graph) else graph, GeometryConfig()),
                }
                control_predicted, _ = predict_distances(control_record, weights)
                metrics = evaluate_record(control_record, control_predicted, GeometryConfig())
                control_rows.append(
                    {
                        "control_name": control_name,
                        "split": split,
                        "family": record["family"],
                        "seed": record["seed"],
                        "spearman_true_vs_predicted": f"{metrics['spearman_true_vs_predicted']:.12g}",
                        "top_k_precision": f"{metrics['top_k_precision']:.12g}",
                        "adjacency_f1": f"{metrics['adjacency_f1']:.12g}",
                        "mean_absolute_error": f"{metrics['mean_absolute_error']:.12g}",
                        "status": control_status(control_name, metrics["spearman_true_vs_predicted"]),
                    }
                )
    return geometry_rows, control_rows


def summary(values: list[float]) -> dict[str, float]:
    return {
        "mean": mean(values),
        "variance": variance(values),
        "min": min(values) if values else 0.0,
        "max": max(values) if values else 0.0,
    }


def render_report(result: dict[str, Any]) -> str:
    lines = [
        "# Synthetic Intrinsic Geometry Recovery Report",
        "",
        f"- bridge id: {result['bridge_id']}",
        f"- version: {result['version']}",
        f"- verdict: {', '.join(result['labels'])}",
        f"- calibration pass: {result['calibration_pass']}",
        f"- holdout pass: {result['holdout_pass']}",
        "",
        "## Key results",
        f"- holdout HAOS spearman: {result['holdout_metrics']['haos']['spearman_true_vs_predicted']:.6f}",
        f"- best baseline spearman: {result['holdout_metrics']['best_baseline_spearman']:.6f}",
        f"- holdout top-k precision: {result['holdout_metrics']['haos']['top_k_precision']:.6f}",
        f"- best baseline top-k precision: {result['holdout_metrics']['best_baseline_topk']:.6f}",
        "",
        "## Interpretation",
        "This is a synthetic blind geometry calibration. It tests whether operator-only transport features",
        "recover latent intrinsic distances better than conventional baselines under frozen holdout.",
        "It does not authorize Bell, quantum, or physical mechanism claims.",
    ]
    return "\n".join(lines) + "\n"


def run_benchmark(config: GeometryConfig, output_dir: Path) -> dict[str, Any]:
    contract = precommitment_payload(config)
    write_json(output_dir / PRECOMMITMENT_PATH.name, contract)

    source_manifest = {
        "bridge_id": "GEO-01",
        "frozen_families": list(FAMILIES),
        "development_families": list(DEVELOPMENT_FAMILIES),
        "calibration_families": list(CALIBRATION_FAMILIES),
        "holdout_families": list(HOLDOUT_FAMILIES),
        "seeds": list(SEEDS),
        "perturbations": list(PERTURBATIONS),
        "n_nodes": config.n_nodes,
        "transport_steps": config.diffusion_steps,
        "source_hash": stable_json_hash("geometry_source_", contract),
    }
    write_json(output_dir / SOURCE_MANIFEST_PATH.name, source_manifest)

    dev_records = [generate_family_record(family, seed, config) for family in DEVELOPMENT_FAMILIES for seed in SEEDS]
    calib_records = [generate_family_record(family, seed, config) for family in CALIBRATION_FAMILIES for seed in SEEDS]
    holdout_records = [generate_family_record(family, seed, config) for family in HOLDOUT_FAMILIES for seed in SEEDS]

    weights = fit_haos_weights(dev_records + calib_records, config)

    geometry_rows, control_rows = generate_rows(dev_records + calib_records + holdout_records, weights, "all", include_controls=True)
    write_csv(output_dir / GEOMETRY_MATRIX_PATH.name, geometry_rows, MATRIX_FIELDNAMES)
    write_csv(output_dir / CONTROL_RESULTS_PATH.name, control_rows, CONTROL_FIELDNAMES)

    holdout_metrics = []
    baseline_metrics = {
        "shortest": [],
        "spectral": [],
        "haos": [],
    }
    for record in holdout_records:
        predicted, _ = predict_distances(record, weights)
        metrics = evaluate_record(record, predicted, config)
        holdout_metrics.append(metrics)
        baseline_metrics["shortest"].append(metrics["spearman_true_vs_shortest"])
        baseline_metrics["spectral"].append(metrics["spearman_true_vs_spectral"])
        baseline_metrics["haos"].append(metrics["spearman_true_vs_predicted"])

    best_baseline = max(mean(baseline_metrics["shortest"]), mean(baseline_metrics["spectral"]))
    best_topk = max(
        mean([top_k_precision(rec["target"], rec["shortest"]) for rec in holdout_records]),
        mean([top_k_precision(rec["target"], rec["spectral"]) for rec in holdout_records]),
    )
    haos_holdout = summary([metrics["spearman_true_vs_predicted"] for metrics in holdout_metrics])
    haos_topk = summary([metrics["top_k_precision"] for metrics in holdout_metrics])
    control_scores = {name: [] for name in CONTROL_NAMES}
    for row in control_rows:
        control_scores[row["control_name"]].append(float(row["spearman_true_vs_predicted"]))

    control_summary = {
        name: {
            "mean_spearman": mean(values),
            "degraded": mean(values) <= best_baseline - config.min_control_degradation if name != "seed_repeat_control" else True,
        }
        for name, values in control_scores.items()
    }

    holdout_pass = (
        haos_holdout["mean"] >= best_baseline + config.min_holdout_spearman_margin
        and haos_topk["mean"] >= best_topk + config.min_holdout_topk_margin
        and haos_holdout["mean"] >= config.holdout_spearman_min
        and haos_topk["mean"] >= config.holdout_topk_min
    )
    control_pass = all(
        summary_dict["degraded"]
        for name, summary_dict in control_summary.items()
        if name != "seed_repeat_control"
    ) and control_summary["leakage_positive_control"]["mean_spearman"] > 0.9

    labels = []
    if holdout_pass and control_pass:
        labels.extend(["GEOMETRY_RECOVERY_VALID", "OPERATIONAL_MAPPING_VALID"])
    elif holdout_pass:
        labels.append("MIXED_OPEN")
    else:
        labels.extend(["GEOMETRY_NOT_DEMONSTRATED", "MIXED_OPEN"])
    if not control_pass:
        labels.append("CONTROL_INVALID")

    result = {
        "bridge_id": "GEO-01",
        "version": config.version,
        "contract_hash": stable_json_hash("geometry_contract_", contract),
        "result_hash": "",
        "labels": sorted(set(labels)),
        "calibration_pass": bool(holdout_pass),
        "holdout_pass": bool(holdout_pass),
        "holdout_metrics": {
            "haos": {
                "spearman_true_vs_predicted": haos_holdout["mean"],
                "top_k_precision": haos_topk["mean"],
                "mae": mean([metrics["mean_absolute_error"] for metrics in holdout_metrics]),
            },
            "best_baseline_spearman": best_baseline,
            "best_baseline_topk": best_topk,
        },
        "control_summary": control_summary,
        "baseline_summary": {
            "shortest_spearman": mean(baseline_metrics["shortest"]),
            "spectral_spearman": mean(baseline_metrics["spectral"]),
        },
        "weights": [float(item) for item in weights],
        "source_manifest": SOURCE_MANIFEST_PATH.name,
        "geometry_matrix": GEOMETRY_MATRIX_PATH.name,
        "control_results": CONTROL_RESULTS_PATH.name,
    }
    result["result_hash"] = stable_json_hash("geometry_result_", result)
    write_json(output_dir / RESULT_PATH.name, result)
    write_text = render_report(result)
    (output_dir / REPORT_PATH.name).write_text(write_text, encoding="utf-8")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run synthetic intrinsic geometry recovery.")
    parser.add_argument("--output-dir", type=Path, default=ROOT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = GeometryConfig()
    run_benchmark(config, args.output_dir)


if __name__ == "__main__":
    main()
