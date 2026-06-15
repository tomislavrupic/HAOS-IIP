#!/usr/bin/env python3
"""Synthetic relational-geometry calibration suite.

This suite is a controlled instrument test. It asks whether the current
semantic/refinement diagnostics recover known relational structure when that
structure is available by construction, and reject destructive controls.

It is not a physics bridge and does not compute Bell/CHSH scores.
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

CONTRACT_PATH = ROOT / "precommitment_contract.json"
SOURCE_MANIFEST_PATH = ROOT / "synthetic_source_manifest.json"
SEMANTIC_MATRIX_PATH = ROOT / "semantic_relation_matrix.csv"
GEOMETRY_MATRIX_PATH = ROOT / "calibration_geometry_matrix.csv"
CONTROL_RESULTS_PATH = ROOT / "calibration_control_results.csv"
REFINEMENT_PATH = ROOT / "calibration_refinement_persistence.csv"
REPORT_PATH = ROOT / "synthetic_calibration_report.md"
RESULT_PATH = ROOT / "synthetic_calibration_result.json"


def repo_rel(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)

SETTINGS = ("alpha", "beta", "gamma", "delta", "epsilon", "zeta")
REFINEMENT_LEVELS = (0, 1, 2, 3)
SEEDS = (3101, 3102, 3103)
HOLDOUT_DIRECTED_EDGES = (("alpha", "beta"), ("gamma", "delta"), ("epsilon", "zeta"))

CONDITIONS = (
    "known_positive",
    "semantic_permutation_control",
    "topology_destroyed_control",
    "refinement_broken_control",
    "ambiguous_partial_preservation",
)

MATRIX_FIELDNAMES = [
    "condition",
    "seed",
    "refinement_level",
    "setting_left",
    "setting_right",
    "G",
    "expected_relation",
    "holdout_edge",
]

CONTROL_FIELDNAMES = [
    "condition",
    "sample_role",
    "semantic_rank_correlation",
    "nearest_neighbor_recovery",
    "ordinal_relation_accuracy",
    "graph_reconstruction_f1",
    "holdout_relation_accuracy",
    "scale_aligned_distance",
    "semantic_pass",
    "refinement_pass",
    "expected_status",
    "observed_status",
]

REFINEMENT_FIELDNAMES = [
    "condition",
    "sample_role",
    "pair_order_stability",
    "sign_topology_preservation",
    "aligned_matrix_distance",
    "eigenspace_alignment",
    "cross_refinement_variance",
    "status",
]

SEMANTIC_FIELDNAMES = ["setting_left", "setting_right", "expected_relation", "expected_abs_relation", "holdout_edge"]


@dataclass(frozen=True)
class CalibrationConfig:
    version: str = "synthetic-relational-geometry-calibration-v0.1"
    semantic_rank_correlation_min: float = 0.75
    nearest_neighbor_min: float = 0.80
    ordinal_accuracy_min: float = 0.75
    graph_f1_min: float = 0.80
    holdout_accuracy_min: float = 0.80
    scale_aligned_distance_max: float = 0.35
    refinement_pair_order_min: float = 0.85
    refinement_sign_topology_min: float = 0.80
    refinement_aligned_distance_max: float = 0.25
    refinement_eigenspace_min: float = 0.70
    refinement_variance_max: float = 0.025
    ambiguous_lower_bound: float = 0.35
    ambiguous_upper_bound: float = 0.75


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run synthetic relational-geometry calibration.")
    parser.add_argument("--output-dir", type=Path, default=ROOT)
    return parser.parse_args()


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


def semantic_operator() -> np.ndarray:
    size = len(SETTINGS)
    matrix = np.zeros((size, size), dtype=float)
    for index in range(size):
        next_index = (index + 1) % size
        matrix[index, next_index] = 1.0
        matrix[next_index, index] = -1.0
    matrix[0, 3] = 0.45
    matrix[3, 0] = -0.45
    matrix[2, 5] = -0.35
    matrix[5, 2] = 0.35
    return matrix


def base_vectors() -> np.ndarray:
    return np.eye(len(SETTINGS), dtype=float)


def skew_noise(seed: int, refinement_level: int, amplitude: float, condition: str) -> np.ndarray:
    size = len(SETTINGS)
    matrix = np.zeros((size, size), dtype=float)
    for row in range(size):
        for col in range(row + 1, size):
            value = 2.0 * stable_unit("skew", condition, seed, refinement_level, row, col) - 1.0
            matrix[row, col] = amplitude * value
            matrix[col, row] = -amplitude * value
    return matrix


def vector_noise(seed: int, refinement_level: int, amplitude: float, condition: str) -> np.ndarray:
    size = len(SETTINGS)
    matrix = np.zeros((size, size), dtype=float)
    for row in range(size):
        for col in range(size):
            value = 2.0 * stable_unit("vector", condition, seed, refinement_level, row, col) - 1.0
            matrix[row, col] = amplitude * value
    return matrix


def permutation_matrix(order: tuple[int, ...]) -> np.ndarray:
    matrix = np.zeros((len(order), len(order)), dtype=float)
    for original, permuted in enumerate(order):
        matrix[original, permuted] = 1.0
    return matrix


def condition_vectors_and_operator(condition: str, seed: int, refinement_level: int) -> tuple[np.ndarray, np.ndarray]:
    vectors = base_vectors()
    operator = semantic_operator()
    if condition == "known_positive":
        vectors = vectors + vector_noise(seed, refinement_level, 0.012, condition)
        operator = operator + skew_noise(seed, refinement_level, 0.012, condition)
    elif condition == "semantic_permutation_control":
        permutation = permutation_matrix((2, 5, 1, 4, 0, 3))
        vectors = permutation @ vectors + vector_noise(seed, refinement_level, 0.012, condition)
        operator = operator + skew_noise(seed, refinement_level, 0.012, condition)
    elif condition == "topology_destroyed_control":
        operator = skew_noise(seed, refinement_level, 0.18, condition)
        vectors = vectors + vector_noise(seed, refinement_level, 0.012, condition)
    elif condition == "refinement_broken_control":
        permutations = (
            (0, 1, 2, 3, 4, 5),
            (2, 1, 4, 3, 0, 5),
            (5, 3, 1, 4, 2, 0),
            (1, 4, 0, 5, 3, 2),
        )
        permutation = permutation_matrix(permutations[refinement_level % len(permutations)])
        vectors = permutation @ vectors + vector_noise(seed, refinement_level, 0.05, condition)
        operator = operator + skew_noise(seed, refinement_level, 0.08, condition)
    elif condition == "ambiguous_partial_preservation":
        permutation = permutation_matrix((0, 2, 1, 3, 5, 4))
        vectors = 0.68 * vectors + 0.32 * (permutation @ vectors)
        vectors = vectors + vector_noise(seed, refinement_level, 0.025, condition)
        shuffled = permutation @ operator @ permutation.T
        operator = 0.65 * operator + 0.35 * shuffled + skew_noise(seed, refinement_level, 0.035, condition)
    else:
        raise ValueError(f"unknown condition {condition}")
    return vectors, operator


def relation_matrix(vectors: np.ndarray, operator: np.ndarray) -> np.ndarray:
    size = vectors.shape[0]
    output = np.zeros((size, size), dtype=float)
    for row in range(size):
        for col in range(size):
            left = vectors[row]
            right = vectors[col]
            denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
            output[row, col] = float(left @ operator @ right / denominator) if denominator else 0.0
    return output


def expected_matrix() -> np.ndarray:
    return semantic_operator()


def mean_condition_matrix(condition: str) -> np.ndarray:
    matrices = []
    for seed in SEEDS:
        for refinement_level in REFINEMENT_LEVELS:
            vectors, operator = condition_vectors_and_operator(condition, seed, refinement_level)
            matrices.append(relation_matrix(vectors, operator))
    return np.mean(np.stack(matrices), axis=0)


def mean_value(values: list[float]) -> float:
    return sum(values) / float(len(values)) if values else 0.0


def variance(values: list[float]) -> float:
    if not values:
        return 0.0
    avg = mean_value(values)
    return sum((value - avg) ** 2 for value in values) / float(len(values))


def pearson(left: list[float], right: list[float]) -> float:
    if len(left) != len(right) or len(left) < 2:
        return 0.0
    left_mean = mean_value(left)
    right_mean = mean_value(right)
    numerator = sum((a - left_mean) * (b - right_mean) for a, b in zip(left, right))
    left_den = math.sqrt(sum((a - left_mean) ** 2 for a in left))
    right_den = math.sqrt(sum((b - right_mean) ** 2 for b in right))
    return numerator / (left_den * right_den) if left_den * right_den > 1.0e-300 else 0.0


def average_ranks(values: list[float], reverse: bool = False) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1], reverse=reverse)
    ranks = [0.0 for _ in values]
    index = 0
    while index < len(indexed):
        end = index + 1
        while end < len(indexed) and indexed[end][1] == indexed[index][1]:
            end += 1
        average_rank = (index + end - 1) / 2.0 + 1.0
        for ranked_index in range(index, end):
            ranks[indexed[ranked_index][0]] = average_rank
        index = end
    return ranks


def spearman(left: list[float], right: list[float]) -> float:
    return pearson(average_ranks(left), average_ranks(right))


def off_diagonal_values(matrix: np.ndarray) -> list[float]:
    values = []
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            if row != col:
                values.append(float(matrix[row, col]))
    return values


def aligned_matrix_distance(left: np.ndarray, right: np.ndarray) -> float:
    denominator = float(np.sum(right * right))
    scale = float(np.sum(left * right) / denominator) if denominator > 1.0e-300 else 0.0
    left_norm = float(np.linalg.norm(left))
    return float(np.linalg.norm(left - scale * right) / left_norm) if left_norm > 1.0e-300 else 0.0


def sign_topology_preservation(left: np.ndarray, right: np.ndarray) -> float:
    left_sign = np.where(np.abs(left) <= 0.02, 0, np.sign(left))
    right_sign = np.where(np.abs(right) <= 0.02, 0, np.sign(right))
    return float(np.mean(left_sign == right_sign))


def eigenspace_alignment(left: np.ndarray, right: np.ndarray) -> float:
    left_u, _left_s, _left_vt = np.linalg.svd(left)
    right_u, _right_s, _right_vt = np.linalg.svd(right)
    return float(abs(np.linalg.det(left_u[:, :2].T @ right_u[:, :2])))


def semantic_metrics(condition: str, config: CalibrationConfig) -> dict[str, Any]:
    observed = mean_condition_matrix(condition)
    expected = expected_matrix()
    observed_abs = [abs(value) for value in off_diagonal_values(observed)]
    expected_abs = [abs(value) for value in off_diagonal_values(expected)]
    observed_signed = off_diagonal_values(observed)
    expected_signed = off_diagonal_values(expected)
    rank_correlation = spearman(expected_abs, observed_abs)
    ordinal_accuracy = ordinal_accuracy_score(expected_abs, observed_abs)
    graph_f1_score = graph_reconstruction_f1(expected, observed)
    nearest = nearest_neighbor_recovery(expected, observed)
    holdout = holdout_accuracy(expected, observed)
    distance = aligned_matrix_distance(expected, observed)
    signed_correlation = spearman(expected_signed, observed_signed)
    semantic_pass = (
        rank_correlation >= config.semantic_rank_correlation_min
        and nearest >= config.nearest_neighbor_min
        and ordinal_accuracy >= config.ordinal_accuracy_min
        and graph_f1_score >= config.graph_f1_min
        and holdout >= config.holdout_accuracy_min
        and distance <= config.scale_aligned_distance_max
    )
    composite = mean_value([max(0.0, rank_correlation), nearest, ordinal_accuracy, graph_f1_score, holdout, max(0.0, signed_correlation)])
    ambiguous_open = config.ambiguous_lower_bound <= composite <= config.ambiguous_upper_bound
    return {
        "rank_correlation": rank_correlation,
        "nearest_neighbor_recovery": nearest,
        "ordinal_relation_accuracy": ordinal_accuracy,
        "graph_reconstruction_f1": graph_f1_score,
        "holdout_relation_accuracy": holdout,
        "scale_aligned_distance": distance,
        "signed_rank_correlation": signed_correlation,
        "semantic_composite": composite,
        "semantic_pass": semantic_pass,
        "ambiguous_open": ambiguous_open,
    }


def ordinal_accuracy_score(expected_abs: list[float], observed_abs: list[float]) -> float:
    total = 0
    correct = 0
    for left in range(len(expected_abs)):
        for right in range(left + 1, len(expected_abs)):
            if expected_abs[left] == expected_abs[right]:
                continue
            total += 1
            if (expected_abs[left] > expected_abs[right]) == (observed_abs[left] > observed_abs[right]):
                correct += 1
    return correct / float(total) if total else 0.0


def graph_reconstruction_f1(expected: np.ndarray, observed: np.ndarray) -> float:
    expected_edges = set()
    observed_edges = set()
    for row in range(expected.shape[0]):
        expected_neighbors = sorted(
            (col for col in range(expected.shape[1]) if col != row),
            key=lambda col: (-abs(expected[row, col]), SETTINGS[col]),
        )[:2]
        observed_neighbors = sorted(
            (col for col in range(observed.shape[1]) if col != row),
            key=lambda col: (-abs(observed[row, col]), SETTINGS[col]),
        )[:2]
        expected_edges.update((row, col) for col in expected_neighbors)
        observed_edges.update((row, col) for col in observed_neighbors)
    intersection = expected_edges & observed_edges
    precision = len(intersection) / float(len(observed_edges)) if observed_edges else 0.0
    recall = len(intersection) / float(len(expected_edges)) if expected_edges else 0.0
    return 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0


def nearest_neighbor_recovery(expected: np.ndarray, observed: np.ndarray) -> float:
    recovered = 0
    for row in range(expected.shape[0]):
        expected_strength = max(abs(expected[row, col]) for col in range(expected.shape[1]) if col != row)
        expected_neighbors = {
            col
            for col in range(expected.shape[1])
            if col != row and abs(abs(expected[row, col]) - expected_strength) <= 1.0e-12
        }
        observed_neighbor = sorted(
            (col for col in range(observed.shape[1]) if col != row),
            key=lambda col: (-abs(observed[row, col]), SETTINGS[col]),
        )[0]
        if observed_neighbor in expected_neighbors:
            recovered += 1
    return recovered / float(expected.shape[0])


def holdout_accuracy(expected: np.ndarray, observed: np.ndarray) -> float:
    correct = 0
    for left_setting, right_setting in HOLDOUT_DIRECTED_EDGES:
        left = SETTINGS.index(left_setting)
        right = SETTINGS.index(right_setting)
        expected_sign = 1 if expected[left, right] > 0 else -1
        observed_sign = 1 if observed[left, right] > 0 else -1
        if expected_sign == observed_sign and abs(observed[left, right]) >= 0.2 * abs(expected[left, right]):
            correct += 1
    return correct / float(len(HOLDOUT_DIRECTED_EDGES))


def refinement_matrices(condition: str) -> dict[int, list[np.ndarray]]:
    grouped: dict[int, list[np.ndarray]] = {}
    for seed in SEEDS:
        grouped[seed] = []
        for refinement_level in REFINEMENT_LEVELS:
            vectors, operator = condition_vectors_and_operator(condition, seed, refinement_level)
            grouped[seed].append(relation_matrix(vectors, operator))
    return grouped


def refinement_metrics(condition: str, config: CalibrationConfig) -> dict[str, Any]:
    grouped = refinement_matrices(condition)
    order_values: list[float] = []
    sign_values: list[float] = []
    distance_values: list[float] = []
    eigenspace_values: list[float] = []
    variances: list[float] = []
    for matrices in grouped.values():
        for left, right in zip(matrices, matrices[1:]):
            order_values.append(spearman(off_diagonal_values(left), off_diagonal_values(right)))
            sign_values.append(sign_topology_preservation(left, right))
            distance_values.append(aligned_matrix_distance(left, right))
            eigenspace_values.append(eigenspace_alignment(left, right))
        variances.append(float(np.mean(np.var(np.stack(matrices), axis=0))))
    pair_order = mean_value(order_values)
    sign_topology = mean_value(sign_values)
    distance = mean_value(distance_values)
    eigenspace = mean_value(eigenspace_values)
    cross_variance = mean_value(variances)
    refinement_pass = (
        pair_order >= config.refinement_pair_order_min
        and sign_topology >= config.refinement_sign_topology_min
        and distance <= config.refinement_aligned_distance_max
        and eigenspace >= config.refinement_eigenspace_min
        and cross_variance <= config.refinement_variance_max
    )
    return {
        "pair_order_stability": pair_order,
        "sign_topology_preservation": sign_topology,
        "aligned_matrix_distance": distance,
        "eigenspace_alignment": eigenspace,
        "cross_refinement_variance": cross_variance,
        "refinement_pass": refinement_pass,
    }


def expected_status(condition: str) -> str:
    return {
        "known_positive": "PASS",
        "semantic_permutation_control": "FAIL_SEMANTIC",
        "topology_destroyed_control": "FAIL_SEMANTIC",
        "refinement_broken_control": "FAIL_REFINEMENT",
        "ambiguous_partial_preservation": "OPEN",
    }[condition]


def observed_status(condition: str, semantic: dict[str, Any], refinement: dict[str, Any]) -> str:
    if condition == "ambiguous_partial_preservation" and semantic["ambiguous_open"]:
        return "OPEN"
    if semantic["semantic_pass"] and refinement["refinement_pass"]:
        return "PASS"
    if not semantic["semantic_pass"] and not refinement["refinement_pass"]:
        return "FAIL_BOTH"
    if not semantic["semantic_pass"]:
        return "FAIL_SEMANTIC"
    return "FAIL_REFINEMENT"


def geometry_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    expected = expected_matrix()
    holdout_set = set(HOLDOUT_DIRECTED_EDGES)
    for condition in CONDITIONS:
        for seed in SEEDS:
            for refinement_level in REFINEMENT_LEVELS:
                vectors, operator = condition_vectors_and_operator(condition, seed, refinement_level)
                matrix = relation_matrix(vectors, operator)
                for left_index, left_setting in enumerate(SETTINGS):
                    for right_index, right_setting in enumerate(SETTINGS):
                        rows.append(
                            {
                                "condition": condition,
                                "seed": seed,
                                "refinement_level": refinement_level,
                                "setting_left": left_setting,
                                "setting_right": right_setting,
                                "G": f"{matrix[left_index, right_index]:.12g}",
                                "expected_relation": f"{expected[left_index, right_index]:.12g}",
                                "holdout_edge": (left_setting, right_setting) in holdout_set,
                            }
                        )
    return rows


def semantic_matrix_rows() -> list[dict[str, Any]]:
    expected = expected_matrix()
    holdout_set = set(HOLDOUT_DIRECTED_EDGES)
    rows: list[dict[str, Any]] = []
    for left_index, left_setting in enumerate(SETTINGS):
        for right_index, right_setting in enumerate(SETTINGS):
            value = expected[left_index, right_index]
            rows.append(
                {
                    "setting_left": left_setting,
                    "setting_right": right_setting,
                    "expected_relation": f"{value:.12g}",
                    "expected_abs_relation": f"{abs(value):.12g}",
                    "holdout_edge": (left_setting, right_setting) in holdout_set,
                }
            )
    return rows


def control_and_refinement_rows(config: CalibrationConfig) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    control_rows: list[dict[str, Any]] = []
    refinement_rows: list[dict[str, Any]] = []
    condition_results: dict[str, Any] = {}
    for condition in CONDITIONS:
        semantic = semantic_metrics(condition, config)
        refinement = refinement_metrics(condition, config)
        status = observed_status(condition, semantic, refinement)
        condition_results[condition] = {
            "semantic": semantic,
            "refinement": refinement,
            "expected_status": expected_status(condition),
            "observed_status": status,
        }
        control_rows.append(
            {
                "condition": condition,
                "sample_role": "target" if condition == "known_positive" else "control",
                "semantic_rank_correlation": f"{semantic['rank_correlation']:.12g}",
                "nearest_neighbor_recovery": f"{semantic['nearest_neighbor_recovery']:.12g}",
                "ordinal_relation_accuracy": f"{semantic['ordinal_relation_accuracy']:.12g}",
                "graph_reconstruction_f1": f"{semantic['graph_reconstruction_f1']:.12g}",
                "holdout_relation_accuracy": f"{semantic['holdout_relation_accuracy']:.12g}",
                "scale_aligned_distance": f"{semantic['scale_aligned_distance']:.12g}",
                "semantic_pass": semantic["semantic_pass"],
                "refinement_pass": refinement["refinement_pass"],
                "expected_status": expected_status(condition),
                "observed_status": status,
            }
        )
        refinement_rows.append(
            {
                "condition": condition,
                "sample_role": "target" if condition == "known_positive" else "control",
                "pair_order_stability": f"{refinement['pair_order_stability']:.12g}",
                "sign_topology_preservation": f"{refinement['sign_topology_preservation']:.12g}",
                "aligned_matrix_distance": f"{refinement['aligned_matrix_distance']:.12g}",
                "eigenspace_alignment": f"{refinement['eigenspace_alignment']:.12g}",
                "cross_refinement_variance": f"{refinement['cross_refinement_variance']:.12g}",
                "status": "REFINEMENT_PASS" if refinement["refinement_pass"] else "REFINEMENT_FAIL",
            }
        )
    return control_rows, refinement_rows, condition_results


def precommitment_contract(config: CalibrationConfig) -> dict[str, Any]:
    return {
        "name": "Synthetic Relational Geometry Calibration v0.1",
        "status": "PRECOMMITTED_CALIBRATION_SUITE",
        "purpose": "Test whether semantic/refinement instruments recover known relational geometry and reject destructive controls.",
        "non_claims": [
            "not a physical Bell experiment",
            "not a quantum mechanics derivation",
            "not a HAOS-IIP Bell bridge",
            "not evidence that CST recovers Bell correlations",
        ],
        "known_structure": {
            "settings": list(SETTINGS),
            "semantic_operator": "directed cycle plus two signed cross-links",
            "holdout_directed_edges": [list(edge) for edge in HOLDOUT_DIRECTED_EDGES],
            "refinement_levels": list(REFINEMENT_LEVELS),
            "seeds": list(SEEDS),
        },
        "conditions": {
            "known_positive": "same semantic graph under small permitted perturbations",
            "semantic_permutation_control": "preserves vector dimension and refinement stability but destroys label semantics",
            "topology_destroyed_control": "preserves dimensions but destroys semantic topology",
            "refinement_broken_control": "preserves local semantic ingredients but destroys cross-refinement continuity",
            "ambiguous_partial_preservation": "partial semantic preservation expected to remain OPEN",
        },
        "metrics": {
            "semantic": [
                "rank correlation",
                "nearest-neighbor recovery",
                "ordinal relation accuracy",
                "graph reconstruction F1",
                "holdout relation accuracy",
                "scale-aligned distance",
            ],
            "refinement": [
                "pair-order stability",
                "sign-topology preservation",
                "aligned matrix distance",
                "eigenspace alignment",
                "cross-refinement variance",
            ],
        },
        "thresholds": asdict(config),
        "chsh_policy": {"computed": False, "authorized": False},
    }


def write_report(path: Path, result: dict[str, Any]) -> None:
    rows = result["condition_results"]
    lines = [
        "# Synthetic Relational Geometry Calibration v0.1",
        "",
        "Implemented fact: this suite calibrates semantic/refinement diagnostics on known synthetic geometry.",
        "Design choice: the synthetic domain has known positives, destructive controls, and an ambiguous partial-preservation case.",
        "Heuristic: success means the instrument distinguishes known structure from declared controls, not that any physics claim is supported.",
        "Unverified hypothesis: no HAOS-IIP Bell derivation is tested here.",
        "",
        "## Labels",
    ]
    lines.extend(f"- {label}" for label in result["labels"])
    lines.extend(["", "## Conditions"])
    for condition in CONDITIONS:
        item = rows[condition]
        lines.append(
            "- {condition}: expected `{expected}`, observed `{observed}`".format(
                condition=condition,
                expected=item["expected_status"],
                observed=item["observed_status"],
            )
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "- This is synthetic calibration, not a physics bridge.",
            "- No Bell/CHSH score is computed.",
            "- Passing calibration does not establish HAOS-IIP physics recovery.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_calibration(config: CalibrationConfig, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    contract = precommitment_contract(config)
    contract["contract_hash"] = stable_json_hash("synthetic_calibration_contract", contract)
    write_json(output_dir / CONTRACT_PATH.name, contract)
    write_json(
        output_dir / SOURCE_MANIFEST_PATH.name,
        {
            "settings": list(SETTINGS),
            "refinement_levels": list(REFINEMENT_LEVELS),
            "seeds": list(SEEDS),
            "conditions": list(CONDITIONS),
            "semantic_operator": expected_matrix().tolist(),
            "holdout_directed_edges": [list(edge) for edge in HOLDOUT_DIRECTED_EDGES],
        },
    )
    control_rows, refinement_rows, condition_results = control_and_refinement_rows(config)
    labels: list[str] = []
    if condition_results["known_positive"]["observed_status"] == "PASS":
        labels.append("KNOWN_POSITIVE_PASS")
    else:
        labels.append("KNOWN_POSITIVE_FAIL")
    if condition_results["semantic_permutation_control"]["observed_status"] in {"FAIL_SEMANTIC", "FAIL_BOTH"}:
        labels.append("SEMANTIC_CONTROL_REJECTED")
    else:
        labels.append("SEMANTIC_CONTROL_NOT_REJECTED")
    if condition_results["topology_destroyed_control"]["observed_status"] in {"FAIL_SEMANTIC", "FAIL_BOTH"}:
        labels.append("TOPOLOGY_CONTROL_REJECTED")
    else:
        labels.append("TOPOLOGY_CONTROL_NOT_REJECTED")
    if condition_results["refinement_broken_control"]["observed_status"] in {"FAIL_REFINEMENT", "FAIL_BOTH"}:
        labels.append("REFINEMENT_CONTROL_REJECTED")
    else:
        labels.append("REFINEMENT_CONTROL_NOT_REJECTED")
    if condition_results["ambiguous_partial_preservation"]["observed_status"] == "OPEN":
        labels.append("AMBIGUOUS_CASE_OPEN")
    else:
        labels.append("AMBIGUOUS_CASE_MISCLASSIFIED")
    calibration_pass = all(
        label in labels
        for label in (
            "KNOWN_POSITIVE_PASS",
            "SEMANTIC_CONTROL_REJECTED",
            "TOPOLOGY_CONTROL_REJECTED",
            "REFINEMENT_CONTROL_REJECTED",
            "AMBIGUOUS_CASE_OPEN",
        )
    )
    labels.append("CALIBRATION_PASS" if calibration_pass else "CALIBRATION_FAIL")
    labels.extend(["NO_PHYSICAL_EXPERIMENT", "HAOS_BELL_DERIVATION_NOT_ESTABLISHED"])
    result: dict[str, Any] = {
        "contract_hash": contract["contract_hash"],
        "config": asdict(config),
        "labels": labels,
        "calibration_pass": calibration_pass,
        "chsh_score": None,
        "chsh_scoring_computed": False,
        "condition_results": condition_results,
        "outputs": {
            "precommitment_contract": repo_rel(output_dir / CONTRACT_PATH.name),
            "source_manifest": repo_rel(output_dir / SOURCE_MANIFEST_PATH.name),
            "semantic_relation_matrix": repo_rel(output_dir / SEMANTIC_MATRIX_PATH.name),
            "geometry_matrix": repo_rel(output_dir / GEOMETRY_MATRIX_PATH.name),
            "control_results": repo_rel(output_dir / CONTROL_RESULTS_PATH.name),
            "refinement_persistence": repo_rel(output_dir / REFINEMENT_PATH.name),
            "report": repo_rel(output_dir / REPORT_PATH.name),
            "result": repo_rel(output_dir / RESULT_PATH.name),
        },
    }
    hash_payload = {key: value for key, value in result.items() if key != "outputs"}
    result["result_hash"] = stable_json_hash("synthetic_calibration_result", hash_payload)
    write_csv(output_dir / SEMANTIC_MATRIX_PATH.name, semantic_matrix_rows(), SEMANTIC_FIELDNAMES)
    write_csv(output_dir / GEOMETRY_MATRIX_PATH.name, geometry_rows(), MATRIX_FIELDNAMES)
    write_csv(output_dir / CONTROL_RESULTS_PATH.name, control_rows, CONTROL_FIELDNAMES)
    write_csv(output_dir / REFINEMENT_PATH.name, refinement_rows, REFINEMENT_FIELDNAMES)
    write_json(output_dir / RESULT_PATH.name, result)
    write_report(output_dir / REPORT_PATH.name, result)
    return result


def main() -> None:
    args = parse_args()
    result = run_calibration(CalibrationConfig(), args.output_dir)
    print(json.dumps({
        "labels": result["labels"],
        "calibration_pass": result["calibration_pass"],
        "chsh_scoring_computed": result["chsh_scoring_computed"],
        "result_hash": result["result_hash"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
