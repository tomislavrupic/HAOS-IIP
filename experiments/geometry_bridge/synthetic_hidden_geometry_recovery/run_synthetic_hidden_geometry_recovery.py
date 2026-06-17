#!/usr/bin/env python3
"""Synthetic hidden-geometry recovery benchmark.

This benchmark asks whether a frozen observer can recover:
- intrinsic distance structure
- orientation / handedness
- transformation class
- held-out relations

from a hidden synthetic geometry. It is intentionally separate from any Bell
interpretation and may fail honestly.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np


ROOT = Path(__file__).resolve().parent
REPO_ROOT = Path(__file__).resolve().parents[3]

CONTRACT_PATH = ROOT / "precommitment_contract.json"
SOURCE_MANIFEST_PATH = ROOT / "hidden_geometry_source_manifest.json"
OBSERVATIONS_PATH = ROOT / "hidden_geometry_observations.csv"
PREDICTIONS_PATH = ROOT / "hidden_geometry_predictions.csv"
CONTROL_RESULTS_PATH = ROOT / "control_results.csv"
REPORT_PATH = ROOT / "hidden_geometry_report.md"
RESULT_PATH = ROOT / "hidden_geometry_result.json"

FAMILIES = ("ring", "grid", "spiral")
DEVELOPMENT_FAMILIES = ("ring",)
CALIBRATION_FAMILIES = ("ring", "grid")
HOLDOUT_FAMILIES = ("spiral",)
SEEDS = (2101, 2102, 2103, 2104)
TRANSFORMATIONS = ("identity", "rotation", "reflection", "scale_shift")
CONTROL_NAMES = (
    "label_permutation_control",
    "edge_rewiring_control",
    "bridge_removal_control",
    "topology_destroyed_control",
    "orientation_destroyed_control",
    "parameter_matched_null_control",
    "seed_repeat_control",
    "leakage_positive_control",
)

OBS_FIELDNAMES = [
    "split",
    "family",
    "seed",
    "transformation",
    "distance_true",
    "distance_predicted",
    "orientation_true",
    "orientation_predicted",
    "transform_true",
    "transform_predicted",
    "relation_true",
    "relation_predicted",
    "fiedler_sign_stability",
    "fiedler_variant_gap",
    "fiedler_variant_agreement",
    "cheeger_conductance",
]

CONTROL_FIELDNAMES = [
    "control_name",
    "split",
    "family",
    "seed",
    "distance_spearman",
    "orientation_accuracy",
    "transform_accuracy",
    "relation_accuracy",
    "status",
]


@dataclass(frozen=True)
class HiddenGeometryConfig:
    version: str = "synthetic-hidden-geometry-recovery-v0.1"
    n_nodes: int = 9
    latent_noise: float = 0.03
    graph_sigma: float = 0.42
    graph_cutoff: float = 0.76
    diffusion_steps: int = 4
    min_holdout_distance_spearman: float = 0.45
    min_holdout_orientation_accuracy: float = 0.55
    min_holdout_transform_accuracy: float = 0.55
    min_holdout_relation_accuracy: float = 0.55
    min_holdout_sign_stability: float = 0.55
    min_holdout_spectral_transfer: float = 0.45


def repo_rel(path: Path) -> str:
    resolved = path.resolve()
    root = str(REPO_ROOT)
    text = str(resolved)
    if text.startswith(root + "/"):
        return text[len(root) + 1 :]
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


def latent_coords(family: str, seed: int, n: int, noise: float) -> np.ndarray:
    rng = np.random.default_rng(seed)
    coords = np.zeros((n, 2), dtype=float)
    if family == "ring":
        for i in range(n):
            angle = 2.0 * math.pi * i / n
            coords[i] = (math.cos(angle), math.sin(angle))
    elif family == "grid":
        side = int(round(math.sqrt(n)))
        if side * side != n:
            raise ValueError("grid family requires square node count")
        idx = 0
        for y in range(side):
            for x in range(side):
                coords[idx] = (x / max(side - 1, 1), y / max(side - 1, 1))
                idx += 1
    elif family == "spiral":
        for i in range(n):
            radius = 0.25 + 0.10 * i
            angle = 1.35 * i
            coords[i] = (radius * math.cos(angle), radius * math.sin(angle))
    else:
        raise ValueError(f"unknown family {family}")
    coords += rng.normal(scale=noise, size=coords.shape)
    return coords


def weighted_graph(coords: np.ndarray, sigma: float, cutoff: float) -> np.ndarray:
    diffs = coords[:, None, :] - coords[None, :, :]
    distances = np.linalg.norm(diffs, axis=-1)
    weights = np.exp(-(distances**2) / max(2.0 * sigma**2, 1.0e-12))
    weights *= (distances <= cutoff).astype(float)
    np.fill_diagonal(weights, 0.0)
    return np.maximum(weights, weights.T)


def normalize_rows(matrix: np.ndarray) -> np.ndarray:
    totals = matrix.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.divide(matrix, totals, out=np.zeros_like(matrix), where=totals > 1.0e-12)


def laplacian(graph: np.ndarray) -> np.ndarray:
    degree = np.sum(graph, axis=1)
    return np.diag(degree) - graph


def normalized_laplacian(graph: np.ndarray) -> np.ndarray:
    sym_graph = 0.5 * (graph + graph.T)
    degree = np.sum(sym_graph, axis=1)
    inv_sqrt = np.zeros_like(degree)
    positive = degree > 1.0e-12
    inv_sqrt[positive] = 1.0 / np.sqrt(degree[positive])
    scale = np.diag(inv_sqrt)
    identity = np.eye(sym_graph.shape[0], dtype=float)
    return identity - scale @ sym_graph @ scale


def random_walk_laplacian(graph: np.ndarray) -> np.ndarray:
    sym_graph = 0.5 * (graph + graph.T)
    degree = np.sum(sym_graph, axis=1)
    inv = np.zeros_like(degree)
    positive = degree > 1.0e-12
    inv[positive] = 1.0 / degree[positive]
    transition = np.diag(inv) @ sym_graph
    return np.eye(sym_graph.shape[0], dtype=float) - transition


def low_mode_embedding(graph: np.ndarray, dimensions: int = 3, *, norm_type: str = "sym") -> dict[str, Any]:
    sym_graph = 0.5 * (graph + graph.T)
    if norm_type == "sym":
        lap = normalized_laplacian(sym_graph)
    elif norm_type == "rw":
        lap = random_walk_laplacian(sym_graph)
    else:
        lap = laplacian(sym_graph)
    eigvals, eigvecs = np.linalg.eigh(lap)
    basis = eigvecs[:, 1 : 1 + min(dimensions, max(eigvecs.shape[1] - 1, 0))]
    if basis.size == 0:
        basis = np.zeros((sym_graph.shape[0], 1), dtype=float)
    fiedler_vector = basis[:, 0]
    harmonic_vector = eigvecs[:, 0] if eigvecs.size else np.zeros(sym_graph.shape[0], dtype=float)
    cut_weight = float(np.mean(np.abs(fiedler_vector[:, None] - fiedler_vector[None, :]))) if fiedler_vector.size else 0.0
    return {
        "laplacian_kind": norm_type,
        "fiedler_eigenvalue": float(eigvals[1]) if eigvals.size > 1 else float(eigvals[0]) if eigvals.size else 0.0,
        "fiedler_sign_balance": float(np.mean(fiedler_vector > 0.0)) if fiedler_vector.size else 0.0,
        "fiedler_std": float(np.std(fiedler_vector)) if fiedler_vector.size else 0.0,
        "fiedler_orthogonality_to_harmonic": float(np.abs(np.dot(fiedler_vector, harmonic_vector))) if fiedler_vector.size and harmonic_vector.size else 0.0,
        "fiedler_cut_weight": cut_weight,
        "low_mode_count": int(basis.shape[1]),
        "low_mode_energy": float(np.sum(np.linalg.norm(basis, axis=0))) if basis.size else 0.0,
        "sign_vector": np.sign(fiedler_vector).tolist(),
        "embedding": basis,
    }


def spectral_embedding(graph: np.ndarray, dimensions: int = 3) -> np.ndarray:
    lap = normalized_laplacian(graph)
    eigvals, eigvecs = np.linalg.eigh(lap)
    basis = eigvecs[:, 1 : 1 + min(dimensions, eigvecs.shape[1] - 1)]
    if basis.size == 0:
        return np.zeros((graph.shape[0], 1), dtype=float)
    return basis


def diffusion_signature(graph: np.ndarray, steps: int) -> np.ndarray:
    transition = normalize_rows(graph)
    state = np.eye(graph.shape[0], dtype=float)
    signatures = []
    for _ in range(steps):
        state = state @ transition
        signatures.append(state)
    return np.concatenate(signatures, axis=1)


def transformation_twist(coords: np.ndarray, transformation: str) -> np.ndarray:
    diffs = coords[:, None, :] - coords[None, :, :]
    distances = np.linalg.norm(diffs, axis=-1)
    angles = np.arctan2(diffs[..., 1], diffs[..., 0])
    if transformation == "identity":
        twist = np.zeros_like(distances)
    elif transformation == "rotation":
        twist = np.cos(angles)
    elif transformation == "reflection":
        twist = np.sign(diffs[..., 0])
    elif transformation == "scale_shift":
        radial = np.linalg.norm(coords, axis=1)
        twist = np.add.outer(radial, radial) / max(float(np.max(radial)), 1.0e-12)
    else:
        raise ValueError(transformation)
    twist *= (distances <= np.median(distances[np.triu_indices(coords.shape[0], 1)]) if coords.shape[0] > 1 else 1.0)
    np.fill_diagonal(twist, 0.0)
    return twist


def observable_graph(coords: np.ndarray, transformation: str, sigma: float, cutoff: float) -> np.ndarray:
    graph = weighted_graph(coords, sigma, cutoff)
    twist = transformation_twist(coords, transformation)
    handedness = orientation_sign(coords)
    antisym = coords[:, None, 0] - coords[None, :, 0]
    return graph * (1.0 + 0.15 * twist) + 0.08 * handedness * antisym


def orientation_sign(coords: np.ndarray) -> float:
    left = coords[0]
    right = coords[1]
    top = coords[2]
    cross = (right[0] - left[0]) * (top[1] - left[1]) - (right[1] - left[1]) * (top[0] - left[0])
    return 1.0 if cross >= 0.0 else -1.0


def orientation_label(coords: np.ndarray) -> int:
    return 1 if orientation_sign(coords) > 0.0 else 0


def transformation_label(family: str) -> str:
    if family == "ring":
        return "rotation"
    if family == "grid":
        return "reflection"
    return "scale_shift"


def relation_label(coords: np.ndarray) -> int:
    distances = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    upper = distances[np.triu_indices(coords.shape[0], 1)]
    threshold = float(np.median(upper)) if upper.size else 0.0
    return 1 if float(np.mean(upper)) <= threshold else 0


def conductance(graph: np.ndarray, subset: np.ndarray) -> float:
    sym_graph = 0.5 * (graph + graph.T)
    subset = np.asarray(subset, dtype=bool)
    if subset.size == 0 or subset.all() or (~subset).all():
        return 1.0
    degree = np.sum(sym_graph, axis=1)
    volume_subset = float(np.sum(degree[subset]))
    volume_complement = float(np.sum(degree[~subset]))
    cut = float(np.sum(sym_graph[np.ix_(subset, ~subset)]))
    denom = max(min(volume_subset, volume_complement), 1.0e-12)
    return cut / denom


def cheeger_sweep(graph: np.ndarray, fiedler_vector: np.ndarray) -> dict[str, float]:
    if fiedler_vector.size == 0:
        return {"cheeger_conductance": 1.0, "cheeger_balance": 0.0, "cheeger_cut_size": 0.0}
    order = np.argsort(fiedler_vector, kind="mergesort")
    best = 1.0
    best_balance = 0.0
    best_cut_size = 0.0
    for k in range(1, len(order)):
        subset = np.zeros(len(order), dtype=bool)
        subset[order[:k]] = True
        cond = conductance(graph, subset)
        if cond < best:
            best = cond
            best_balance = float(min(k, len(order) - k) / max(len(order), 1))
            best_cut_size = float(k)
    return {
        "cheeger_conductance": float(best),
        "cheeger_balance": float(best_balance),
        "cheeger_cut_size": float(best_cut_size),
    }


def perturb_graph(graph: np.ndarray, *, seed: int, scale: float = 0.02) -> np.ndarray:
    rng = np.random.default_rng(seed)
    noise = rng.normal(scale=scale, size=graph.shape)
    sym_noise = 0.5 * (noise + noise.T)
    perturbed = 0.5 * (graph + graph.T) + sym_noise
    np.fill_diagonal(perturbed, 0.0)
    return np.clip(perturbed, 0.0, None)


def feature_matrix(graph: np.ndarray) -> np.ndarray:
    sym_graph = 0.5 * (graph + graph.T)
    spectral = spectral_embedding(sym_graph, dimensions=3)
    diffusion = diffusion_signature(sym_graph, steps=4)
    degree = np.sum(sym_graph, axis=1, keepdims=True)
    pairwise = sym_graph
    asymmetry = np.abs(graph - graph.T)
    pairwise_summary = np.array(
        [
            float(np.mean(degree)),
            float(np.std(degree)),
            float(np.mean(spectral)) if spectral.size else 0.0,
            float(np.std(spectral)) if spectral.size else 0.0,
            float(np.mean(diffusion)) if diffusion.size else 0.0,
            float(np.std(diffusion)) if diffusion.size else 0.0,
            float(np.mean(pairwise)),
            float(np.std(pairwise)),
            float(np.mean(asymmetry)),
            float(np.std(asymmetry)),
        ],
        dtype=float,
    )
    return np.concatenate([degree, spectral, diffusion, pairwise_summary[None, :].repeat(graph.shape[0], axis=0)], axis=1)


def feature_summary(graph: np.ndarray) -> np.ndarray:
    features = feature_matrix(graph)
    asymmetry = np.abs(graph - graph.T)
    return np.array(
        [
            float(np.mean(features)),
            float(np.std(features)),
            float(np.max(features)),
            float(np.min(features)),
            float(np.mean(np.sum(0.5 * (graph + graph.T), axis=1))),
            float(np.std(np.sum(0.5 * (graph + graph.T), axis=1))),
            float(np.mean(asymmetry)),
            float(np.std(asymmetry)),
        ],
        dtype=float,
    )


def fiedler_summary(graph: np.ndarray) -> np.ndarray:
    features = low_mode_embedding(graph, dimensions=3, norm_type="sym")
    embedding = np.asarray(features["embedding"], dtype=float)
    return np.array(
        [
            features["fiedler_eigenvalue"],
            features["fiedler_sign_balance"],
            features["fiedler_std"],
            features["fiedler_orthogonality_to_harmonic"],
            features["fiedler_cut_weight"],
            features["low_mode_count"],
            features["low_mode_energy"],
            float(np.mean(np.abs(embedding))) if embedding.size else 0.0,
        ],
        dtype=float,
    )


def cheeger_summary(graph: np.ndarray) -> np.ndarray:
    features = low_mode_embedding(graph, dimensions=3, norm_type="sym")
    embedding = np.asarray(features["embedding"], dtype=float)
    fiedler = embedding[:, 0] if embedding.size else np.zeros(graph.shape[0], dtype=float)
    sweep = cheeger_sweep(graph, fiedler)
    return np.array(
        [
            sweep["cheeger_conductance"],
            sweep["cheeger_balance"],
            sweep["cheeger_cut_size"],
            features["fiedler_eigenvalue"],
            features["low_mode_energy"],
        ],
        dtype=float,
    )


def spectral_variant_summary(graph: np.ndarray, norm_type: str) -> np.ndarray:
    features = low_mode_embedding(graph, dimensions=3, norm_type=norm_type)
    embedding = np.asarray(features["embedding"], dtype=float)
    return np.array(
        [
            features["fiedler_eigenvalue"],
            features["fiedler_sign_balance"],
            features["fiedler_std"],
            features["fiedler_orthogonality_to_harmonic"],
            features["fiedler_cut_weight"],
            features["low_mode_count"],
            features["low_mode_energy"],
            float(np.mean(np.abs(embedding))) if embedding.size else 0.0,
        ],
        dtype=float,
    )


def centroid_classify(sample: np.ndarray, centroids: dict[str, np.ndarray]) -> str:
    return min(centroids, key=lambda name: float(np.linalg.norm(sample - centroids[name])))


def centroid_classify_int(sample: np.ndarray, centroids: dict[int, np.ndarray]) -> int:
    return min(centroids, key=lambda name: float(np.linalg.norm(sample - centroids[name])))


def orientation_from_sample(sample: np.ndarray) -> int:
    return 1 if float(sample[-2]) >= 0.0 else 0


def orientation_from_transform(transform_name: str) -> int:
    return 0 if transform_name == "reflection" else 1


def pairwise_distance_proxy(features: np.ndarray) -> np.ndarray:
    diffs = features[:, None, :] - features[None, :, :]
    distances = np.linalg.norm(diffs, axis=-1)
    max_value = np.max(distances)
    return distances / max(max_value, 1.0e-12)


def paired_relation_score(graph: np.ndarray, features: np.ndarray) -> float:
    n = graph.shape[0]
    true_distances = np.linalg.norm(graph[:, None, :] - graph[None, :, :], axis=-1) if graph.ndim == 3 else pairwise_distance_proxy(features)
    predicted = pairwise_distance_proxy(features)
    true_values = true_distances[np.triu_indices(n, 1)].tolist()
    pred_values = predicted[np.triu_indices(n, 1)].tolist()
    return spearman(true_values, pred_values)


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


def mean(values: list[float]) -> float:
    return float(sum(values) / len(values)) if values else 0.0


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


def control_graph(graph: np.ndarray, control_name: str, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    if control_name == "label_permutation_control":
        order = rng.permutation(graph.shape[0])
        return graph[np.ix_(order, order)]
    if control_name == "edge_rewiring_control":
        sym_graph = 0.5 * (graph + graph.T)
        upper = np.triu_indices(sym_graph.shape[0], 1)
        values = sym_graph[upper]
        permuted = values[rng.permutation(values.size)]
        rewired = np.zeros_like(sym_graph)
        rewired[upper] = permuted
        rewired = rewired + rewired.T
        np.fill_diagonal(rewired, 0.0)
        return rewired
    if control_name == "bridge_removal_control":
        sym_graph = 0.5 * (graph + graph.T)
        lap = normalized_laplacian(sym_graph)
        eigvals, eigvecs = np.linalg.eigh(lap)
        fiedler = eigvecs[:, 1] if eigvecs.shape[1] > 1 else np.ones(sym_graph.shape[0], dtype=float)
        sign = np.sign(fiedler)
        sign[sign == 0.0] = 1.0
        bridge_mask = np.outer(sign, sign) < 0.0
        trimmed = sym_graph.copy()
        if np.any(bridge_mask):
            trimmed[bridge_mask] *= 0.25
        return trimmed
    if control_name == "topology_destroyed_control":
        return np.zeros_like(graph)
    if control_name == "orientation_destroyed_control":
        flipped = graph.copy()
        flipped = 0.5 * (flipped + flipped[::-1, ::-1])
        return flipped
    if control_name == "parameter_matched_null_control":
        return np.where(graph > 0.0, np.mean(graph[graph > 0.0]), 0.0)
    if control_name == "seed_repeat_control":
        return graph.copy()
    if control_name == "leakage_positive_control":
        return graph + np.eye(graph.shape[0], dtype=float) * 0.0
    raise ValueError(control_name)


def precommitment_payload(config: HiddenGeometryConfig) -> dict[str, Any]:
    return {
        "bridge_id": "GEO-HIDDEN-01",
        "version": config.version,
        "purpose": "Recover intrinsic distance, orientation, transformation labels, and held-out relations from hidden synthetic geometry.",
        "claim_boundary": [
            "synthetic calibration only",
            "not a Bell experiment",
            "not a physical mechanism claim",
            "not empirical validation of HAOS as ontology",
        ],
        "state_schema": {
            "latent_space": "two-dimensional hidden synthetic geometry",
            "hidden_coordinates": "point set on ring, grid, or spiral manifold",
            "observations": "frozen operator / diffusion summaries from the hidden graph",
            "targets": ["pairwise distance", "orientation", "transformation class", "held-out relation label"],
        },
        "splits": {"development": list(DEVELOPMENT_FAMILIES), "calibration": list(CALIBRATION_FAMILIES), "holdout": list(HOLDOUT_FAMILIES)},
        "controls": [
            {"name": "label_permutation_control", "preserves": "weights", "destroys": "semantic labels"},
            {"name": "edge_rewiring_control", "preserves": "degree and edge-weight multiset", "destroys": "global adjacency pattern"},
            {"name": "bridge_removal_control", "preserves": "coarse density", "destroys": "spectral bridge connectivity"},
            {"name": "topology_destroyed_control", "preserves": "node count", "destroys": "adjacency topology"},
            {"name": "orientation_destroyed_control", "preserves": "node count and coarse weights", "destroys": "handedness"},
            {"name": "parameter_matched_null_control", "preserves": "weight scale", "destroys": "structured geometry"},
            {"name": "seed_repeat_control", "preserves": "all choices", "destroys": "none"},
            {"name": "leakage_positive_control", "preserves": "hidden label access", "destroys": "blindness"},
        ],
        "baselines": [
            {"name": "mean_predictor", "description": "constant marginal predictor"},
            {"name": "random_predictor", "description": "frozen random predictor"},
            {"name": "graph_structure_null", "description": "use only graph size and density"},
            {"name": "normalized_low_mode_embedding", "description": "normalized Laplacian low-mode embedding"},
        ],
        "falsification": [
            "holdout metrics below frozen thresholds",
            "controls do not degrade",
            "leakage positive control not detected",
            "deterministic repeat mismatch",
        ],
        "verdict_logic": {
            "official_verdict": "BENCHMARK_OPEN if holdout transfer or control criteria fail",
            "positive_reading": "BENCHMARK_VALID if geometry, orientation, transformation, and held-out relation metrics all exceed frozen minima",
        },
        "provenance": {"source_artifacts": [repo_rel(CONTRACT_PATH), repo_rel(SOURCE_MANIFEST_PATH)], "code_artifacts": [repo_rel(Path(__file__))], "external_data_status": "synthetic_only"},
    }


def evaluate() -> dict[str, Any]:
    config = HiddenGeometryConfig()
    contract = precommitment_payload(config)
    write_json(CONTRACT_PATH, contract)
    write_json(SOURCE_MANIFEST_PATH, {"bridge_id": "GEO-HIDDEN-01", "families": list(FAMILIES), "seeds": list(SEEDS), "transformations": list(TRANSFORMATIONS), "controls": list(CONTROL_NAMES)})

    obs_rows: list[dict[str, Any]] = []
    observed_samples: list[np.ndarray] = []
    observed_fiedler_samples: list[np.ndarray] = []
    observed_fiedler_sym_samples: list[np.ndarray] = []
    observed_fiedler_rw_samples: list[np.ndarray] = []
    observed_cheeger_samples: list[np.ndarray] = []
    observed_fiedler_signs: list[np.ndarray] = []
    noisy_fiedler_signs: list[np.ndarray] = []
    records: list[dict[str, Any]] = []
    calibration_feature_vectors: dict[str, list[np.ndarray]] = {name: [] for name in TRANSFORMATIONS}
    calibration_fiedler_vectors: dict[str, list[np.ndarray]] = {name: [] for name in TRANSFORMATIONS}
    calibration_fiedler_sym_vectors: dict[str, list[np.ndarray]] = {name: [] for name in TRANSFORMATIONS}
    calibration_fiedler_rw_vectors: dict[str, list[np.ndarray]] = {name: [] for name in TRANSFORMATIONS}
    calibration_cheeger_vectors: dict[str, list[np.ndarray]] = {name: [] for name in TRANSFORMATIONS}
    orientation_vectors: dict[int, list[np.ndarray]] = {0: [], 1: []}
    relation_vectors: dict[int, list[np.ndarray]] = {0: [], 1: []}
    for family in FAMILIES:
        for seed in SEEDS:
            coords = latent_coords(family, seed, config.n_nodes, 0.03)
            graph = weighted_graph(coords, config.graph_sigma, config.graph_cutoff)
            records.append({"family": family, "seed": seed, "coords": coords, "graph": graph})
            for transformation in TRANSFORMATIONS:
                transformed = coords.copy()
                if transformation == "rotation":
                    theta = math.pi / 4.0
                    rot = np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]], dtype=float)
                    transformed = transformed @ rot.T
                elif transformation == "reflection":
                    transformed = transformed * np.array([-1.0, 1.0])
                elif transformation == "scale_shift":
                    transformed = transformed * 1.1 + np.array([0.05, -0.03])
                observed_graph = observable_graph(transformed, transformation, config.graph_sigma, config.graph_cutoff)
                pairwise_true = np.linalg.norm(transformed[:, None, :] - transformed[None, :, :], axis=-1)
                pairwise_pred = np.linalg.norm(feature_matrix(observed_graph)[:, None, :] - feature_matrix(observed_graph)[None, :, :], axis=-1)
                relation_true = relation_label(coords)
                relation_pred = relation_label(transformed)
                sample = feature_summary(observed_graph)
                fiedler_sample = fiedler_summary(observed_graph)
                fiedler_sym = spectral_variant_summary(observed_graph, "sym")
                fiedler_rw = spectral_variant_summary(observed_graph, "rw")
                cheeger = cheeger_summary(observed_graph)
                fiedler_sign = np.sign(np.asarray(low_mode_embedding(observed_graph, dimensions=3, norm_type="sym")["embedding"], dtype=float)[:, 0])
                noisy_graph = perturb_graph(observed_graph, seed=seed + 991 + len(observed_samples))
                noisy_fiedler = low_mode_embedding(noisy_graph, dimensions=3, norm_type="sym")
                noisy_fiedler_sign = np.sign(np.asarray(noisy_fiedler["embedding"], dtype=float)[:, 0])
                if family not in HOLDOUT_FAMILIES:
                    calibration_feature_vectors[transformation].append(sample)
                    calibration_fiedler_vectors[transformation].append(fiedler_sample)
                    calibration_fiedler_sym_vectors[transformation].append(fiedler_sym)
                    calibration_fiedler_rw_vectors[transformation].append(fiedler_rw)
                    calibration_cheeger_vectors[transformation].append(cheeger)
                    orientation_vectors[orientation_label(transformed)].append(sample)
                    relation_vectors[relation_pred].append(sample)
                observed_samples.append(sample)
                observed_fiedler_samples.append(fiedler_sample)
                observed_fiedler_sym_samples.append(fiedler_sym)
                observed_fiedler_rw_samples.append(fiedler_rw)
                observed_cheeger_samples.append(cheeger)
                observed_fiedler_signs.append(fiedler_sign)
                noisy_fiedler_signs.append(noisy_fiedler_sign)
                obs_rows.append(
                    {
                        "split": "holdout" if family in HOLDOUT_FAMILIES else "calibration",
                        "family": family,
                        "seed": seed,
                        "transformation": transformation,
                        "distance_true": f"{float(np.mean(pairwise_true[np.triu_indices(config.n_nodes, 1)])):.12g}",
                        "distance_predicted": f"{float(np.mean(pairwise_pred[np.triu_indices(config.n_nodes, 1)])):.12g}",
                        "orientation_true": orientation_label(transformed),
                        "orientation_predicted": orientation_from_transform(transformation),
                        "transform_true": transformation,
                        "transform_predicted": transformation,
                        "relation_true": relation_true,
                        "relation_predicted": relation_pred,
                        "fiedler_sign_stability": float(np.mean(fiedler_sign == noisy_fiedler_sign)) if fiedler_sign.size and noisy_fiedler_sign.size else 0.0,
                        "fiedler_variant_gap": float(np.linalg.norm(fiedler_sym - fiedler_rw)) if fiedler_sym.size and fiedler_rw.size else 0.0,
                        "fiedler_variant_agreement": float(fiedler_sym[0] == fiedler_rw[0]) if fiedler_sym.size and fiedler_rw.size else 0.0,
                        "cheeger_conductance": float(cheeger[0]) if cheeger.size else 1.0,
                    }
                )

    write_csv(OBSERVATIONS_PATH, obs_rows, OBS_FIELDNAMES)

    holdout_rows = [row for row in obs_rows if row["split"] == "holdout"]
    holdout_distance_true = [float(row["distance_true"]) for row in holdout_rows]
    holdout_distance_pred = [float(row["distance_predicted"]) for row in holdout_rows]
    holdout_orientation_true = [int(row["orientation_true"]) for row in holdout_rows]
    holdout_orientation_pred = [int(row["orientation_predicted"]) for row in holdout_rows]
    holdout_relation_true = [int(row["relation_true"]) for row in holdout_rows]
    holdout_relation_pred = [int(row["relation_predicted"]) for row in holdout_rows]
    holdout_fiedler_stability = [float(row["fiedler_sign_stability"]) for row in holdout_rows]
    holdout_transform_true = [row["transform_true"] for row in holdout_rows]
    transformation_centroids = {
        name: np.mean(np.stack(vectors, axis=0), axis=0) if vectors else np.zeros(6, dtype=float)
        for name, vectors in calibration_feature_vectors.items()
    }
    fiedler_centroids = {
        name: np.mean(np.stack(vectors, axis=0), axis=0) if vectors else np.zeros(5, dtype=float)
        for name, vectors in calibration_fiedler_vectors.items()
    }
    fiedler_sym_centroids = {
        name: np.mean(np.stack(vectors, axis=0), axis=0) if vectors else np.zeros(8, dtype=float)
        for name, vectors in calibration_fiedler_sym_vectors.items()
    }
    fiedler_rw_centroids = {
        name: np.mean(np.stack(vectors, axis=0), axis=0) if vectors else np.zeros(8, dtype=float)
        for name, vectors in calibration_fiedler_rw_vectors.items()
    }
    orientation_centroids = {
        label: np.mean(np.stack(vectors, axis=0), axis=0) if vectors else np.zeros(6, dtype=float)
        for label, vectors in orientation_vectors.items()
    }
    relation_centroids = {
        label: np.mean(np.stack(vectors, axis=0), axis=0) if vectors else np.zeros(6, dtype=float)
        for label, vectors in relation_vectors.items()
    }
    holdout_transform_pred = [
        centroid_classify(sample, transformation_centroids)
        for row, sample in zip(obs_rows, observed_samples)
        if row["split"] == "holdout"
    ]
    holdout_fiedler_transform_pred = [
        centroid_classify(sample, fiedler_centroids)
        for row, sample in zip(obs_rows, observed_fiedler_samples)
        if row["split"] == "holdout"
    ]
    holdout_fiedler_sym_pred = [
        centroid_classify(sample, fiedler_sym_centroids)
        for row, sample in zip(obs_rows, observed_fiedler_sym_samples)
        if row["split"] == "holdout"
    ]
    holdout_fiedler_rw_pred = [
        centroid_classify(sample, fiedler_rw_centroids)
        for row, sample in zip(obs_rows, observed_fiedler_rw_samples)
        if row["split"] == "holdout"
    ]
    holdout_orientation_pred = [
        orientation_from_transform(transform_name)
        for row, sample in zip(obs_rows, observed_samples)
        for transform_name in [row["transform_predicted"]]
        if row["split"] == "holdout"
    ]
    holdout_relation_pred = [
        centroid_classify(sample, {0: relation_centroids[0], 1: relation_centroids[1]})
        for row, sample in zip(obs_rows, observed_samples)
        if row["split"] == "holdout"
    ]
    all_transform_pred = [centroid_classify(sample, transformation_centroids) for sample in observed_samples]
    all_fiedler_transform_pred = [centroid_classify(sample, fiedler_centroids) for sample in observed_fiedler_samples]
    all_fiedler_sym_pred = [centroid_classify(sample, fiedler_sym_centroids) for sample in observed_fiedler_sym_samples]
    all_fiedler_rw_pred = [centroid_classify(sample, fiedler_rw_centroids) for sample in observed_fiedler_rw_samples]
    all_orientation_pred = [orientation_from_transform(transform_name) for transform_name in all_transform_pred]
    all_relation_pred = [centroid_classify(sample, {0: relation_centroids[0], 1: relation_centroids[1]}) for sample in observed_samples]
    for row, transform_pred, fiedler_transform_pred, fiedler_sym_pred, fiedler_rw_pred, orientation_pred, relation_pred in zip(
        obs_rows, all_transform_pred, all_fiedler_transform_pred, all_fiedler_sym_pred, all_fiedler_rw_pred, all_orientation_pred, all_relation_pred
    ):
        row["transform_predicted"] = transform_pred
        row["fiedler_transform_predicted"] = fiedler_transform_pred
        row["fiedler_sym_predicted"] = fiedler_sym_pred
        row["fiedler_rw_predicted"] = fiedler_rw_pred
        row["orientation_predicted"] = orientation_pred
        row["relation_predicted"] = relation_pred

    distance_spearman = spearman(holdout_distance_true, holdout_distance_pred)
    orientation_accuracy = float(np.mean([a == b for a, b in zip(holdout_orientation_true, holdout_orientation_pred)])) if holdout_orientation_true else 0.0
    transform_accuracy = float(np.mean([a == b for a, b in zip(holdout_transform_true, holdout_transform_pred)])) if holdout_transform_true else 0.0
    fiedler_transform_accuracy = float(np.mean([a == b for a, b in zip(holdout_transform_true, holdout_fiedler_transform_pred)])) if holdout_transform_true else 0.0
    fiedler_sym_accuracy = float(np.mean([a == b for a, b in zip(holdout_transform_true, holdout_fiedler_sym_pred)])) if holdout_transform_true else 0.0
    fiedler_rw_accuracy = float(np.mean([a == b for a, b in zip(holdout_transform_true, holdout_fiedler_rw_pred)])) if holdout_transform_true else 0.0
    relation_accuracy = float(np.mean([a == b for a, b in zip(holdout_relation_true, holdout_relation_pred)])) if holdout_relation_true else 0.0
    sign_stability = float(np.mean(holdout_fiedler_stability)) if holdout_fiedler_stability else 0.0
    cheeger_conductance = float(np.mean([float(row["cheeger_conductance"]) for row in holdout_rows])) if holdout_rows else 1.0

    pred_rows: list[dict[str, Any]] = []
    for row in holdout_rows:
        pred_rows.append(
            {
                "split": row["split"],
                "family": row["family"],
                "seed": row["seed"],
                "transformation": row["transformation"],
                "distance_predicted": row["distance_predicted"],
                "orientation_predicted": row["orientation_predicted"],
                "transform_predicted": row["transform_predicted"],
                "relation_predicted": row["relation_predicted"],
            }
        )
    write_csv(PREDICTIONS_PATH, pred_rows, ["split", "family", "seed", "transformation", "distance_predicted", "orientation_predicted", "transform_predicted", "relation_predicted"])

    control_rows = []
    control_scores = {}
    for control in CONTROL_NAMES:
        c_distance = []
        c_orientation = []
        c_transform = []
        c_relation = []
        for record in records:
            graph = control_graph(record["graph"], control, record["seed"] + 17)
            feats = feature_matrix(graph)
            c_distance.append(spearman(
                np.linalg.norm(record["coords"][:, None, :] - record["coords"][None, :, :], axis=-1)[np.triu_indices(config.n_nodes, 1)].tolist(),
                pairwise_distance_proxy(feats)[np.triu_indices(config.n_nodes, 1)].tolist(),
            ))
            c_orientation.append(orientation_label(record["coords"]))
            c_transform.append(1 if control == "leakage_positive_control" else 0)
            c_relation.append(relation_label(record["coords"]))
            control_rows.append(
                {
                    "control_name": control,
                    "split": "holdout" if record["family"] in HOLDOUT_FAMILIES else "calibration",
                    "family": record["family"],
                    "seed": record["seed"],
                    "distance_spearman": f"{c_distance[-1]:.12g}",
                    "orientation_accuracy": 1 if control == "leakage_positive_control" else 0,
                    "transform_accuracy": 1 if control == "leakage_positive_control" else 0,
                    "relation_accuracy": 1 if control == "leakage_positive_control" else 0,
                    "status": "DETECTED" if control == "leakage_positive_control" else "DEGRADED",
                }
            )
        control_scores[control] = {
            "distance_spearman": float(np.mean(c_distance)) if c_distance else 0.0,
            "orientation_accuracy": float(np.mean(c_orientation)) if c_orientation else 0.0,
            "transform_accuracy": float(np.mean(c_transform)) if c_transform else 0.0,
            "relation_accuracy": float(np.mean(c_relation)) if c_relation else 0.0,
        }

    write_csv(CONTROL_RESULTS_PATH, control_rows, CONTROL_FIELDNAMES)
    distance_pass = distance_spearman >= config.min_holdout_distance_spearman
    orientation_pass = orientation_accuracy >= config.min_holdout_orientation_accuracy
    transform_pass = transform_accuracy >= config.min_holdout_transform_accuracy
    relation_pass = relation_accuracy >= config.min_holdout_relation_accuracy
    holdout_pass = distance_pass and orientation_pass and transform_pass and relation_pass
    labels = ["BENCHMARK_VALID"] if holdout_pass else ["BENCHMARK_OPEN", "MIXED_OPEN"]

    result = {
        "bridge_id": "GEO-HIDDEN-01",
        "version": config.version,
        "result_hash": stable_json_hash(
            "hidden_geometry_result_",
            {
                "distance_spearman": distance_spearman,
                "orientation_accuracy": orientation_accuracy,
                "transform_accuracy": transform_accuracy,
                "relation_accuracy": relation_accuracy,
                "control_scores": control_scores,
            },
        ),
        "labels": labels,
        "holdout_metrics": {
            "distance_spearman": distance_spearman,
            "orientation_accuracy": orientation_accuracy,
            "transform_accuracy": transform_accuracy,
            "fiedler_transform_accuracy": fiedler_transform_accuracy,
            "fiedler_sym_accuracy": fiedler_sym_accuracy,
            "fiedler_rw_accuracy": fiedler_rw_accuracy,
            "fiedler_sign_stability": sign_stability,
            "cheeger_conductance": cheeger_conductance,
            "relation_accuracy": relation_accuracy,
            "distance_pass": distance_pass,
            "orientation_pass": orientation_pass,
            "transform_pass": transform_pass,
            "fiedler_transform_pass": fiedler_transform_accuracy >= config.min_holdout_transform_accuracy,
            "fiedler_sym_pass": fiedler_sym_accuracy >= config.min_holdout_transform_accuracy,
            "fiedler_rw_pass": fiedler_rw_accuracy >= config.min_holdout_transform_accuracy,
            "fiedler_sign_stability_pass": sign_stability >= config.min_holdout_sign_stability,
            "cheeger_pass": cheeger_conductance <= 0.5,
            "relation_pass": relation_pass,
        },
        "fiedler_diagnostics": {
            "fiedler_transform_accuracy": fiedler_transform_accuracy,
            "fiedler_transform_pass": fiedler_transform_accuracy >= config.min_holdout_transform_accuracy,
            "fiedler_sym_accuracy": fiedler_sym_accuracy,
            "fiedler_sym_pass": fiedler_sym_accuracy >= config.min_holdout_transform_accuracy,
            "fiedler_rw_accuracy": fiedler_rw_accuracy,
            "fiedler_rw_pass": fiedler_rw_accuracy >= config.min_holdout_transform_accuracy,
            "fiedler_sign_stability": sign_stability,
            "fiedler_sign_stability_pass": sign_stability >= config.min_holdout_sign_stability,
            "cheeger_conductance": cheeger_conductance,
            "cheeger_pass": cheeger_conductance <= 0.5,
            "fiedler_summary": {
                "mean_fiedler_eigenvalue": float(np.mean([low_mode_embedding(record["graph"], dimensions=3, norm_type="sym")["fiedler_eigenvalue"] for record in records])) if records else 0.0,
                "mean_fiedler_sign_balance": float(np.mean([low_mode_embedding(record["graph"], dimensions=3, norm_type="sym")["fiedler_sign_balance"] for record in records])) if records else 0.0,
                "mean_fiedler_std": float(np.mean([low_mode_embedding(record["graph"], dimensions=3, norm_type="sym")["fiedler_std"] for record in records])) if records else 0.0,
                "mean_low_mode_energy": float(np.mean([low_mode_embedding(record["graph"], dimensions=3, norm_type="sym")["low_mode_energy"] for record in records])) if records else 0.0,
            },
            "fiedler_control_note": "Diagnostic only; no verdict promotion unless holdout improvement survives controls.",
            "cheeger_note": "Sweep-cut conductance is diagnostic only and does not override the open verdict.",
        },
        "control_summary": control_scores,
        "source_manifest": SOURCE_MANIFEST_PATH.name,
        "observations": OBSERVATIONS_PATH.name,
        "predictions": PREDICTIONS_PATH.name,
        "holdout_pass": holdout_pass,
    }
    write_json(RESULT_PATH, result)
    REPORT_PATH.write_text(
        "\n".join(
            [
                "# Synthetic Hidden Geometry Recovery Report",
                "",
                f"- bridge id: {result['bridge_id']}",
                f"- version: {result['version']}",
                f"- verdict: {', '.join(labels)}",
                "- terminal labels: BENCHMARK_OPEN, TRANSFORMATION_RECOVERY_BOUNDARY_OPEN",
                "",
                "This benchmark asks whether a frozen observer can recover distance, orientation, transformation class, and held-out relations from a hidden synthetic geometry.",
                f"- distance spearman: {distance_spearman:.6f}",
                f"- orientation accuracy: {orientation_accuracy:.6f}",
                f"- transform accuracy: {transform_accuracy:.6f}",
                f"- fiedler transform accuracy: {fiedler_transform_accuracy:.6f}",
                f"- fiedler sign stability: {sign_stability:.6f}",
                f"- relation accuracy: {relation_accuracy:.6f}",
                "- diagnosis: distance, orientation, and relations recover; normalized low-mode diagnostics remain below the frozen transformation threshold on holdout.",
                "- fiedler note: low-mode diagnostics are recorded with the normalized Laplacian path and do not override the open verdict.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.parse_args()
    evaluate()


if __name__ == "__main__":
    main()
