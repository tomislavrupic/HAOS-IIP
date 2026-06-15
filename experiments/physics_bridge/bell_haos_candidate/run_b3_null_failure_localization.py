#!/usr/bin/env python3
"""B3.2.1 null failure localization for the HAOS-IIP Bell geometry bridge.

This script freezes the B3.2 invariant and localizes why nulls preserve too
much of its sign-changing geometry. It does not change J_lambda, setting
vectors, settings, tolerances, or any Bell scoreboard path.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np


ROOT = Path(__file__).resolve().parent
REPO_ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
geometry = importlib.import_module("run_b3_geometry_audit")

CONTRACT_PATH = ROOT / "b3_2_1_precommitment_contract.json"
DISTRIBUTION_PATH = ROOT / "b3_2_1_distribution_comparison.csv"
MATRIX_PATH = ROOT / "b3_2_1_matrix_comparison.csv"
TOPOLOGY_PATH = ROOT / "b3_2_1_topology_dependence.csv"
OPERATOR_PATH = ROOT / "b3_2_1_operator_dependence.csv"
ORDERING_PATH = ROOT / "b3_2_1_pair_ordering.csv"
INVARIANT_LEDGER_PATH = ROOT / "b3_2_1_null_invariant_ledger.json"
REPORT_PATH = ROOT / "b3_2_1_report.md"
RESULT_PATH = ROOT / "b3_2_1_result.json"


def repo_rel(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)

LEFT_SETTINGS = ("a0", "a1", "h0", "h3")
RIGHT_SETTINGS = ("b0", "b1", "h1", "h2")
ANALYSIS_PAIRS = tuple((role, left, right) for role, pairs in geometry.ALL_PAIR_ROLES.items() for left, right in pairs)

CONTROL_IDS = (
    "target_geometry",
    "identity_pairing",
    "shuffled_pairing",
    "random_orthogonal_pairing",
    "topology_destroyed_pairing",
    "label_permuted_settings",
    "closure_vector_permutation",
    "refinement_broken_control",
)

OPERATOR_VARIANTS = (
    "target_J",
    "degree_matched_J",
    "spectrum_matched_J",
    "locality_preserving_randomized_J",
    "orientation_reversed_J",
)

DISTRIBUTION_FIELDNAMES = [
    "control_id",
    "sample_role",
    "count",
    "mean_G",
    "variance_G",
    "min_G",
    "max_G",
    "range_G",
    "positive_count",
    "negative_count",
    "sign_balance",
    "cohens_d_vs_target",
    "distribution_overlap_vs_target",
    "pair_order_spearman_vs_target",
    "pair_order_signature",
]

MATRIX_FIELDNAMES = [
    "control_id",
    "sample_role",
    "frobenius_norm",
    "scale_aligned_distance_to_target",
    "raw_frobenius_distance_to_target",
    "matrix_rank",
    "singular_values",
    "symmetric_norm",
    "antisymmetric_norm",
    "antisymmetric_fraction",
    "matrix_values",
]

TOPOLOGY_FIELDNAMES = [
    "control_id",
    "sample_role",
    "count",
    "pearson_absG_inverse_chain_distance",
    "pearson_absG_adjacency_fraction",
    "adjacent_mean_absG",
    "nonadjacent_mean_absG",
    "pair_order_spearman_vs_target",
    "interpretation",
]

OPERATOR_FIELDNAMES = [
    "operator_variant",
    "preserves",
    "destroys",
    "count",
    "mean_G",
    "variance_G",
    "positive_count",
    "negative_count",
    "pair_order_spearman_vs_target",
    "scale_aligned_distance_to_target",
    "singular_values",
    "matrix_rank",
    "status",
]

ORDERING_FIELDNAMES = [
    "control_id",
    "pair_role",
    "setting_left",
    "setting_right",
    "target_mean_G",
    "control_mean_G",
    "target_rank",
    "control_rank",
    "rank_delta",
]


@dataclass(frozen=True)
class NullLocalizationConfig:
    version: str = "b3-2-1-null-failure-localization-v0.1"
    source_limit: int = 6
    sign_tolerance: float = 0.02
    matrix_similarity_distance_threshold: float = 0.75
    rank_similarity_threshold: float = 0.60
    topology_correlation_attention_threshold: float = 0.20


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run B3.2.1 null failure localization.")
    parser.add_argument("--output-dir", type=Path, default=ROOT)
    return parser.parse_args()


def stable_json_hash(prefix: str, payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(encoded).hexdigest()[:24]}"


def hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def mean(values: list[float]) -> float:
    return sum(values) / float(len(values)) if values else 0.0


def variance(values: list[float]) -> float:
    if not values:
        return 0.0
    avg = mean(values)
    return sum((value - avg) ** 2 for value in values) / float(len(values))


def pearson(left: list[float], right: list[float]) -> float:
    if len(left) != len(right) or len(left) < 2:
        return 0.0
    left_mean = mean(left)
    right_mean = mean(right)
    numerator = sum((a - left_mean) * (b - right_mean) for a, b in zip(left, right))
    left_den = math.sqrt(sum((a - left_mean) ** 2 for a in left))
    right_den = math.sqrt(sum((b - right_mean) ** 2 for b in right))
    denominator = left_den * right_den
    return numerator / denominator if denominator > 1.0e-300 else 0.0


def average_ranks(values: list[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
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
    if len(left) != len(right) or len(left) < 2:
        return 0.0
    return pearson(average_ranks(left), average_ranks(right))


def distribution_overlap(left: list[float], right: list[float], bins: int = 12) -> float:
    if not left or not right:
        return 0.0
    lower = min(min(left), min(right))
    upper = max(max(left), max(right))
    if abs(upper - lower) <= 1.0e-300:
        return 1.0
    width = (upper - lower) / float(bins)
    left_counts = [0 for _ in range(bins)]
    right_counts = [0 for _ in range(bins)]
    for values, counts in ((left, left_counts), (right, right_counts)):
        for value in values:
            index = min(bins - 1, int((value - lower) / width))
            counts[index] += 1
    left_total = float(len(left))
    right_total = float(len(right))
    return sum(min(a / left_total, b / right_total) for a, b in zip(left_counts, right_counts))


def cohen_d(left: list[float], right: list[float]) -> float:
    if not left or not right:
        return 0.0
    var_left = variance(left)
    var_right = variance(right)
    pooled = math.sqrt((var_left + var_right) / 2.0)
    return (mean(right) - mean(left)) / pooled if pooled > 1.0e-300 else 0.0


def token_positions(source: Any) -> dict[str, int]:
    return {token: index for index, token in enumerate(str(source.fields["chain_signature"]).split(">"))}


def setting_chain_distance(source: Any, left_setting: str, right_setting: str) -> tuple[float, float]:
    positions = token_positions(source)
    left_tokens = geometry.SETTING_DEFINITIONS[left_setting]
    right_tokens = geometry.SETTING_DEFINITIONS[right_setting]
    distances: list[float] = []
    adjacent = 0
    for left_token in left_tokens:
        for right_token in right_tokens:
            distance = abs(positions.get(left_token, 0) - positions.get(right_token, 0))
            distances.append(float(distance))
            if distance == 1:
                adjacent += 1
    return mean(distances), adjacent / float(len(distances))


def matrix_operator_variant(source: Any, variant: str) -> list[list[float]]:
    target = geometry.chain_orientation_operator(source)
    if variant == "target_J":
        return target
    if variant == "spectrum_matched_J":
        return geometry.signed_transform_matrix(target, geometry.SIGNED_PERMUTATION)
    if variant == "orientation_reversed_J":
        return [[-value for value in row] for row in target]
    if variant == "locality_preserving_randomized_J":
        chain = str(source.fields["chain_signature"]).split(">")
        matrix = geometry.zero_matrix(len(geometry.TOKENS))
        index_by_token = {token: index for index, token in enumerate(geometry.TOKENS)}
        weights: list[float] = []
        for edge_index, (left_token, right_token) in enumerate(zip(chain, chain[1:])):
            if left_token in index_by_token and right_token in index_by_token:
                weights.append(abs(target[index_by_token[left_token]][index_by_token[right_token]]))
        rotated_weights = list(reversed(weights))
        for edge_index, (left_token, right_token) in enumerate(zip(chain, chain[1:])):
            if left_token not in index_by_token or right_token not in index_by_token:
                continue
            weight = rotated_weights[edge_index % len(rotated_weights)] if rotated_weights else 0.0
            sign = 1.0 if edge_index % 2 == 0 else -1.0
            left_index = index_by_token[left_token]
            right_index = index_by_token[right_token]
            matrix[left_index][right_index] = sign * weight
            matrix[right_index][left_index] = -sign * weight
        return matrix
    if variant == "degree_matched_J":
        degrees = [sum(abs(value) for value in row) for row in target]
        order = sorted(range(len(degrees)), key=lambda index: (-degrees[index], index))
        matrix = geometry.zero_matrix(len(geometry.TOKENS))
        for rank, left_index in enumerate(order):
            right_index = order[(rank + 2) % len(order)]
            if left_index == right_index:
                continue
            weight = min(degrees[left_index], degrees[right_index]) / 2.0
            sign = 1.0 if rank % 2 == 0 else -1.0
            matrix[left_index][right_index] += sign * weight
            matrix[right_index][left_index] -= sign * weight
        return matrix
    raise ValueError(f"unknown operator variant {variant}")


def g_with_operator(source: Any, left_setting: str, right_setting: str, operator: list[list[float]]) -> float:
    left = geometry.setting_vector(source, left_setting)
    right = geometry.setting_vector(source, right_setting)
    denominator = geometry.norm(left) * geometry.norm(right)
    if denominator <= 1.0e-300:
        return 0.0
    return geometry.dot(left, geometry.matvec(operator, right)) / denominator


def control_sources(control_id: str, target_sources: list[Any], refinement_sources: list[Any]) -> tuple[str, list[Any]]:
    if control_id == "refinement_broken_control":
        return "control", refinement_sources
    if control_id == "target_geometry":
        return "target", target_sources
    return "null", target_sources


def pair_values_for_control(control_id: str, target_sources: list[Any], refinement_sources: list[Any]) -> list[dict[str, Any]]:
    sample_role, sources = control_sources(control_id, target_sources, refinement_sources)
    rows: list[dict[str, Any]] = []
    for source in sources:
        for pair_role, left_setting, right_setting in ANALYSIS_PAIRS:
            value = geometry.relational_invariant_g(source, left_setting, right_setting, control_id)
            distance, adjacency = setting_chain_distance(source, left_setting, right_setting)
            rows.append(
                {
                    "control_id": control_id,
                    "sample_role": sample_role,
                    "source_id": source.source_id,
                    "pair_role": pair_role,
                    "setting_left": left_setting,
                    "setting_right": right_setting,
                    "pair_key": f"{pair_role}:{left_setting}_{right_setting}",
                    "G": value,
                    "abs_G": abs(value),
                    "chain_distance": distance,
                    "inverse_chain_distance": 1.0 / (1.0 + distance),
                    "adjacency_fraction": adjacency,
                }
            )
    return rows


def pair_values_for_operator(variant: str, sources: list[Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for source in sources:
        operator = matrix_operator_variant(source, variant)
        for pair_role, left_setting, right_setting in ANALYSIS_PAIRS:
            value = g_with_operator(source, left_setting, right_setting, operator)
            rows.append(
                {
                    "operator_variant": variant,
                    "source_id": source.source_id,
                    "pair_role": pair_role,
                    "setting_left": left_setting,
                    "setting_right": right_setting,
                    "pair_key": f"{pair_role}:{left_setting}_{right_setting}",
                    "G": value,
                }
            )
    return rows


def relation_matrix_for_control(control_id: str, target_sources: list[Any], refinement_sources: list[Any]) -> np.ndarray:
    _sample_role, sources = control_sources(control_id, target_sources, refinement_sources)
    matrices = []
    for source in sources:
        rows = []
        for left_setting in LEFT_SETTINGS:
            rows.append([
                geometry.relational_invariant_g(source, left_setting, right_setting, control_id)
                for right_setting in RIGHT_SETTINGS
            ])
        matrices.append(np.array(rows, dtype=float))
    return np.mean(np.stack(matrices), axis=0)


def relation_matrix_for_operator(variant: str, sources: list[Any]) -> np.ndarray:
    matrices = []
    for source in sources:
        operator = matrix_operator_variant(source, variant)
        rows = []
        for left_setting in LEFT_SETTINGS:
            rows.append([g_with_operator(source, left_setting, right_setting, operator) for right_setting in RIGHT_SETTINGS])
        matrices.append(np.array(rows, dtype=float))
    return np.mean(np.stack(matrices), axis=0)


def matrix_stats(matrix: np.ndarray, target_matrix: np.ndarray | None) -> dict[str, Any]:
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    symmetric = 0.5 * (matrix + matrix.T)
    antisymmetric = 0.5 * (matrix - matrix.T)
    frobenius = float(np.linalg.norm(matrix))
    symmetric_norm = float(np.linalg.norm(symmetric))
    antisymmetric_norm = float(np.linalg.norm(antisymmetric))
    raw_distance = 0.0
    aligned_distance = 0.0
    if target_matrix is not None:
        raw_distance = float(np.linalg.norm(matrix - target_matrix))
        denominator = float(np.sum(matrix * matrix))
        scale = float(np.sum(target_matrix * matrix) / denominator) if denominator > 1.0e-300 else 0.0
        target_norm = float(np.linalg.norm(target_matrix))
        aligned_distance = float(np.linalg.norm(target_matrix - scale * matrix) / target_norm) if target_norm > 1.0e-300 else 0.0
    return {
        "frobenius_norm": frobenius,
        "scale_aligned_distance_to_target": aligned_distance,
        "raw_frobenius_distance_to_target": raw_distance,
        "matrix_rank": int(np.linalg.matrix_rank(matrix, tol=1.0e-9)),
        "singular_values": [float(value) for value in singular_values],
        "symmetric_norm": symmetric_norm,
        "antisymmetric_norm": antisymmetric_norm,
        "antisymmetric_fraction": antisymmetric_norm / frobenius if frobenius > 1.0e-300 else 0.0,
        "matrix_values": matrix.tolist(),
    }


def grouped_values(rows: list[dict[str, Any]]) -> dict[str, float]:
    buckets: dict[str, list[float]] = {}
    for row in rows:
        buckets.setdefault(row["pair_key"], []).append(float(row["G"]))
    return {key: mean(values) for key, values in sorted(buckets.items())}


def rank_signature(pair_means: dict[str, float]) -> str:
    ordered = sorted(pair_means.items(), key=lambda item: (-item[1], item[0]))
    return ">".join(pair for pair, _value in ordered)


def distribution_rows(
    all_pair_rows: dict[str, list[dict[str, Any]]],
    config: NullLocalizationConfig,
) -> list[dict[str, Any]]:
    target_values = [row["G"] for row in all_pair_rows["target_geometry"]]
    target_pair_means = grouped_values(all_pair_rows["target_geometry"])
    target_order_values = [target_pair_means[key] for key in sorted(target_pair_means)]
    rows: list[dict[str, Any]] = []
    for control_id in CONTROL_IDS:
        values = [row["G"] for row in all_pair_rows[control_id]]
        pair_means = grouped_values(all_pair_rows[control_id])
        order_values = [pair_means[key] for key in sorted(target_pair_means)]
        positive_count = sum(1 for value in values if value > config.sign_tolerance)
        negative_count = sum(1 for value in values if value < -config.sign_tolerance)
        rows.append(
            {
                "control_id": control_id,
                "sample_role": all_pair_rows[control_id][0]["sample_role"],
                "count": len(values),
                "mean_G": f"{mean(values):.12g}",
                "variance_G": f"{variance(values):.12g}",
                "min_G": f"{min(values):.12g}",
                "max_G": f"{max(values):.12g}",
                "range_G": f"{(max(values)-min(values)):.12g}",
                "positive_count": positive_count,
                "negative_count": negative_count,
                "sign_balance": f"{((positive_count-negative_count)/float(len(values))):.12g}",
                "cohens_d_vs_target": f"{cohen_d(target_values, values):.12g}",
                "distribution_overlap_vs_target": f"{distribution_overlap(target_values, values):.12g}",
                "pair_order_spearman_vs_target": f"{spearman(target_order_values, order_values):.12g}",
                "pair_order_signature": rank_signature(pair_means),
            }
        )
    return rows


def matrix_rows(target_sources: list[Any], refinement_sources: list[Any]) -> list[dict[str, Any]]:
    target_matrix = relation_matrix_for_control("target_geometry", target_sources, refinement_sources)
    rows: list[dict[str, Any]] = []
    for control_id in CONTROL_IDS:
        matrix = relation_matrix_for_control(control_id, target_sources, refinement_sources)
        stats = matrix_stats(matrix, None if control_id == "target_geometry" else target_matrix)
        sample_role, _sources = control_sources(control_id, target_sources, refinement_sources)
        rows.append(
            {
                "control_id": control_id,
                "sample_role": sample_role,
                "frobenius_norm": f"{stats['frobenius_norm']:.12g}",
                "scale_aligned_distance_to_target": f"{stats['scale_aligned_distance_to_target']:.12g}",
                "raw_frobenius_distance_to_target": f"{stats['raw_frobenius_distance_to_target']:.12g}",
                "matrix_rank": stats["matrix_rank"],
                "singular_values": json.dumps(stats["singular_values"], separators=(",", ":")),
                "symmetric_norm": f"{stats['symmetric_norm']:.12g}",
                "antisymmetric_norm": f"{stats['antisymmetric_norm']:.12g}",
                "antisymmetric_fraction": f"{stats['antisymmetric_fraction']:.12g}",
                "matrix_values": json.dumps(stats["matrix_values"], separators=(",", ":")),
            }
        )
    return rows


def topology_rows(all_pair_rows: dict[str, list[dict[str, Any]]]) -> list[dict[str, Any]]:
    target_pair_means = grouped_values(all_pair_rows["target_geometry"])
    target_order_values = [target_pair_means[key] for key in sorted(target_pair_means)]
    rows: list[dict[str, Any]] = []
    for control_id in CONTROL_IDS:
        pair_rows = all_pair_rows[control_id]
        abs_values = [row["abs_G"] for row in pair_rows]
        inverse_distances = [row["inverse_chain_distance"] for row in pair_rows]
        adjacency_fractions = [row["adjacency_fraction"] for row in pair_rows]
        adjacent = [row["abs_G"] for row in pair_rows if row["adjacency_fraction"] > 0.0]
        nonadjacent = [row["abs_G"] for row in pair_rows if row["adjacency_fraction"] <= 0.0]
        pair_means = grouped_values(pair_rows)
        order_values = [pair_means[key] for key in sorted(target_pair_means)]
        distance_corr = pearson(abs_values, inverse_distances)
        adjacency_corr = pearson(abs_values, adjacency_fractions)
        interpretation = "TOPOLOGY_DEPENDENCE_WEAK_OR_GENERIC"
        if control_id == "target_geometry" and abs(distance_corr) >= 0.20:
            interpretation = "TOPOLOGY_DEPENDENCE_VISIBLE_IN_TARGET"
        elif control_id != "target_geometry" and abs(distance_corr) >= 0.20:
            interpretation = "TOPOLOGY_DEPENDENCE_SURVIVES_CONTROL"
        rows.append(
            {
                "control_id": control_id,
                "sample_role": pair_rows[0]["sample_role"],
                "count": len(pair_rows),
                "pearson_absG_inverse_chain_distance": f"{distance_corr:.12g}",
                "pearson_absG_adjacency_fraction": f"{adjacency_corr:.12g}",
                "adjacent_mean_absG": f"{mean(adjacent):.12g}",
                "nonadjacent_mean_absG": f"{mean(nonadjacent):.12g}",
                "pair_order_spearman_vs_target": f"{spearman(target_order_values, order_values):.12g}",
                "interpretation": interpretation,
            }
        )
    return rows


OPERATOR_LEDGER = {
    "target_J": {
        "preserves": "frozen chain-order adjacency, orientation, antisymmetry, closure-strength weighting",
        "destroys": "nothing; reference operator",
    },
    "degree_matched_J": {
        "preserves": "approximate weighted node degree and antisymmetry",
        "destroys": "actual chain adjacency and orientation order",
    },
    "spectrum_matched_J": {
        "preserves": "skew spectrum under signed permutation, norm, antisymmetry",
        "destroys": "semantic token orientation alignment",
    },
    "locality_preserving_randomized_J": {
        "preserves": "chain-local adjacency and antisymmetry",
        "destroys": "edge-weight order and alternating orientations",
    },
    "orientation_reversed_J": {
        "preserves": "chain adjacency, spectrum magnitude, locality, antisymmetry",
        "destroys": "orientation sign convention",
    },
}


def operator_rows(target_sources: list[Any]) -> list[dict[str, Any]]:
    target_pair_rows = pair_values_for_operator("target_J", target_sources)
    target_pair_means = grouped_values(target_pair_rows)
    target_order_values = [target_pair_means[key] for key in sorted(target_pair_means)]
    target_matrix = relation_matrix_for_operator("target_J", target_sources)
    rows: list[dict[str, Any]] = []
    for variant in OPERATOR_VARIANTS:
        pair_rows = pair_values_for_operator(variant, target_sources)
        values = [row["G"] for row in pair_rows]
        pair_means = grouped_values(pair_rows)
        order_values = [pair_means[key] for key in sorted(target_pair_means)]
        matrix = relation_matrix_for_operator(variant, target_sources)
        stats = matrix_stats(matrix, None if variant == "target_J" else target_matrix)
        positive_count = sum(1 for value in values if value > 0.02)
        negative_count = sum(1 for value in values if value < -0.02)
        status = "REFERENCE_OPERATOR" if variant == "target_J" else "OPERATOR_GEOMETRY_PERSISTS"
        if variant != "target_J" and stats["scale_aligned_distance_to_target"] >= 0.75 and abs(spearman(target_order_values, order_values)) <= 0.60:
            status = "OPERATOR_GEOMETRY_DEGRADED"
        rows.append(
            {
                "operator_variant": variant,
                "preserves": OPERATOR_LEDGER[variant]["preserves"],
                "destroys": OPERATOR_LEDGER[variant]["destroys"],
                "count": len(values),
                "mean_G": f"{mean(values):.12g}",
                "variance_G": f"{variance(values):.12g}",
                "positive_count": positive_count,
                "negative_count": negative_count,
                "pair_order_spearman_vs_target": f"{spearman(target_order_values, order_values):.12g}",
                "scale_aligned_distance_to_target": f"{stats['scale_aligned_distance_to_target']:.12g}",
                "singular_values": json.dumps(stats["singular_values"], separators=(",", ":")),
                "matrix_rank": stats["matrix_rank"],
                "status": status,
            }
        )
    return rows


def pair_ordering_rows(all_pair_rows: dict[str, list[dict[str, Any]]]) -> list[dict[str, Any]]:
    target_means = grouped_values(all_pair_rows["target_geometry"])
    target_ranks = {key: rank for key, rank in zip(sorted(target_means), average_ranks([target_means[key] for key in sorted(target_means)]))}
    rows: list[dict[str, Any]] = []
    for control_id in CONTROL_IDS:
        control_means = grouped_values(all_pair_rows[control_id])
        ordered_keys = sorted(target_means)
        control_ranks = {key: rank for key, rank in zip(ordered_keys, average_ranks([control_means[key] for key in ordered_keys]))}
        for pair_key in ordered_keys:
            pair_role, settings = pair_key.split(":")
            setting_left, setting_right = settings.split("_")
            rows.append(
                {
                    "control_id": control_id,
                    "pair_role": pair_role,
                    "setting_left": setting_left,
                    "setting_right": setting_right,
                    "target_mean_G": f"{target_means[pair_key]:.12g}",
                    "control_mean_G": f"{control_means[pair_key]:.12g}",
                    "target_rank": f"{target_ranks[pair_key]:.12g}",
                    "control_rank": f"{control_ranks[pair_key]:.12g}",
                    "rank_delta": f"{(control_ranks[pair_key]-target_ranks[pair_key]):.12g}",
                }
            )
    return rows


def null_invariant_ledger() -> dict[str, Any]:
    return {
        "frozen_status": {
            "B3_2": "FROZEN",
            "relational_geometry": "DETECTED",
            "haos_specificity": "NOT_ESTABLISHED",
            "null_rejection": "FAIL",
            "chsh_scoring": "NOT_AUTHORIZED",
        },
        "controls": {
            "identity_pairing": {
                "dimension": "preserved",
                "norm": "partially preserved",
                "orientation": "destroyed",
                "antisymmetry": "destroyed",
                "locality": "destroyed",
                "chain_adjacency": "destroyed",
                "refinement_continuity": "preserved",
            },
            "shuffled_pairing": {
                "dimension": "preserved",
                "norm": "preserved",
                "orientation": "partially preserved",
                "antisymmetry": "preserved",
                "locality": "destroyed",
                "chain_adjacency": "destroyed",
                "refinement_continuity": "preserved",
            },
            "random_orthogonal_pairing": {
                "dimension": "preserved",
                "norm": "preserved",
                "orientation": "randomized",
                "antisymmetry": "preserved",
                "locality": "destroyed",
                "chain_adjacency": "destroyed",
                "refinement_continuity": "preserved",
            },
            "topology_destroyed_pairing": {
                "dimension": "preserved",
                "norm": "destroyed",
                "orientation": "destroyed",
                "antisymmetry": "trivial_zero",
                "locality": "destroyed",
                "chain_adjacency": "destroyed",
                "refinement_continuity": "preserved",
            },
            "label_permuted_settings": {
                "dimension": "preserved",
                "norm": "preserved",
                "orientation": "preserved",
                "antisymmetry": "preserved",
                "locality": "preserved",
                "chain_adjacency": "preserved",
                "refinement_continuity": "preserved",
                "setting_semantics": "destroyed",
            },
            "closure_vector_permutation": {
                "dimension": "preserved",
                "norm": "preserved",
                "orientation": "preserved",
                "antisymmetry": "preserved",
                "locality": "preserved",
                "chain_adjacency": "preserved",
                "refinement_continuity": "preserved",
                "vector_token_alignment": "destroyed",
            },
            "refinement_broken_control": {
                "dimension": "preserved",
                "norm": "partially preserved",
                "orientation": "preserved_if_chain_signature_survives",
                "antisymmetry": "preserved",
                "locality": "preserved",
                "chain_adjacency": "preserved",
                "refinement_continuity": "destroyed",
            },
        },
        "interpretation": "Nulls preserving orientation algebra, antisymmetry, or setting-vector diversity can reproduce much of the observed geometry. Current evidence localizes the failure to target-specific provenance, not to sign-change detection.",
    }


def precommitment_contract(config: NullLocalizationConfig) -> dict[str, Any]:
    return {
        "name": "B3.2.1_NULL_FAILURE_LOCALIZATION_PRECOMMITMENT",
        "purpose": "Localize why B3.2 nulls preserve sign-changing relational geometry.",
        "frozen_inputs": {
            "invariant": "G(a,b,lambda)=<u_a,J_lambda u_b>/(||u_a||||u_b||)",
            "source_artifact": geometry.GeometryConfig().source_artifact,
            "source_hash": hash_file(REPO_ROOT / geometry.GeometryConfig().source_artifact),
            "settings": {key: list(value) for key, value in geometry.SETTING_DEFINITIONS.items()},
            "left_matrix_settings": list(LEFT_SETTINGS),
            "right_matrix_settings": list(RIGHT_SETTINGS),
            "analysis_pairs": [list(pair) for pair in ANALYSIS_PAIRS],
            "prohibited_changes": [
                "do not alter J_lambda",
                "do not alter closure vectors",
                "do not alter setting mappings",
                "do not alter sign convention",
                "do not alter B3.2 thresholds",
                "do not compute CHSH",
            ],
        },
        "required_analyses": [
            "target/null G distributions",
            "null invariant ledger",
            "topology dependence",
            "operator dependence",
            "relation-matrix geometry",
            "pair-order stability",
            "holdout behavior",
        ],
        "operator_variants": {key: OPERATOR_LEDGER[key] for key in OPERATOR_VARIANTS},
        "matrix_metrics": [
            "singular values",
            "rank",
            "symmetric norm",
            "antisymmetric norm",
            "scale-aligned distance",
            "pair-order Spearman correlation",
        ],
        "thresholds_for_localization_only": asdict(config),
        "chsh_policy": {
            "computed": False,
            "authorized": False,
            "reason": "B3.2 null rejection failed; this script performs failure localization only.",
        },
    }


def write_report(path: Path, result: dict[str, Any], distribution: list[dict[str, Any]], operator: list[dict[str, Any]]) -> None:
    labels = result["labels"]
    failing_nulls = result["failing_nulls"]
    target_distribution = next(row for row in distribution if row["control_id"] == "target_geometry")
    lines = [
        "# B3.2.1 Null Failure Localization",
        "",
        "Implemented fact: B3.2 remains frozen and CHSH is not computed.",
        "Design choice: this harness compares whole G distributions, relation matrices, topology dependence, and operator variants.",
        "Heuristic: matrix similarity and rank-order similarity localize whether nulls preserve geometry provenance.",
        "Unverified hypothesis: HAOS-specific Bell geometry remains unestablished.",
        "",
        "## Labels",
    ]
    lines.extend(f"- {label}" for label in labels)
    lines.extend(
        [
            "",
            "## Target Distribution",
            f"- mean G: `{target_distribution['mean_G']}`",
            f"- variance G: `{target_distribution['variance_G']}`",
            f"- positive/negative count: `{target_distribution['positive_count']}` / `{target_distribution['negative_count']}`",
            f"- pair-order signature: `{target_distribution['pair_order_signature']}`",
            "",
            "## Failing Nulls",
        ]
    )
    lines.extend(f"- {null_id}" for null_id in failing_nulls)
    lines.extend(["", "## Operator Dependence"])
    for row in operator:
        lines.append(
            "- {variant}: status `{status}`, pair-order rho `{rho}`, aligned distance `{distance}`".format(
                variant=row["operator_variant"],
                status=row["status"],
                rho=row["pair_order_spearman_vs_target"],
                distance=row["scale_aligned_distance_to_target"],
            )
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "- This is a failure-localization harness.",
            "- No Bell/CHSH score is computed.",
            "- No B3.2 field, threshold, setting, vector, or operator is changed.",
            "- HAOS Bell derivation is not established.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_localization(config: NullLocalizationConfig, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    contract = precommitment_contract(config)
    contract["contract_hash"] = stable_json_hash("b3_2_1_contract", contract)
    write_json(output_dir / CONTRACT_PATH.name, contract)

    geometry_config = geometry.GeometryConfig(source_limit=config.source_limit)
    target_sources = geometry.load_sources(geometry_config, geometry_config.target_hierarchy_label)
    refinement_sources = geometry.load_sources(geometry_config, geometry_config.refinement_control_label)
    all_pair_rows = {
        control_id: pair_values_for_control(control_id, target_sources, refinement_sources)
        for control_id in CONTROL_IDS
    }

    distribution = distribution_rows(all_pair_rows, config)
    matrices = matrix_rows(target_sources, refinement_sources)
    topology = topology_rows(all_pair_rows)
    operators = operator_rows(target_sources)
    ordering = pair_ordering_rows(all_pair_rows)
    ledger = null_invariant_ledger()

    failing_nulls: list[str] = []
    for row in matrices:
        if row["control_id"] == "target_geometry":
            continue
        aligned_distance = float(row["scale_aligned_distance_to_target"])
        distribution_row = next(item for item in distribution if item["control_id"] == row["control_id"])
        rho = abs(float(distribution_row["pair_order_spearman_vs_target"]))
        if aligned_distance < config.matrix_similarity_distance_threshold or rho > config.rank_similarity_threshold:
            failing_nulls.append(row["control_id"])

    labels = [
        "NULL_FAILURE_LOCALIZED",
        "GENERIC_BILINEAR_GEOMETRY_SUSPECTED",
        "TARGET_SPECIFIC_GEOMETRY_NOT_ESTABLISHED",
        "CHSH_SCORING_NOT_AUTHORIZED",
        "HAOS_BELL_DERIVATION_NOT_ESTABLISHED",
        "MIXED / OPEN",
    ]
    if "topology_destroyed_pairing" not in failing_nulls:
        labels.append("TOPOLOGY_ZERO_CONTROL_REJECTED")
    if "refinement_broken_control" in failing_nulls:
        labels.append("REFINEMENT_SPECIFICITY_NOT_ESTABLISHED")
    if "label_permuted_settings" in failing_nulls:
        labels.append("SETTING_SEMANTICS_NOT_ESTABLISHED")

    result = {
        "contract_hash": contract["contract_hash"],
        "config": asdict(config),
        "labels": labels,
        "failing_nulls": failing_nulls,
        "chsh_score": None,
        "chsh_scoring_computed": False,
        "chsh_scoring_authorized": False,
        "b3_2_frozen_status": ledger["frozen_status"],
        "outputs": {
            "precommitment_contract": repo_rel(output_dir / CONTRACT_PATH.name),
            "distribution_comparison": repo_rel(output_dir / DISTRIBUTION_PATH.name),
            "matrix_comparison": repo_rel(output_dir / MATRIX_PATH.name),
            "topology_dependence": repo_rel(output_dir / TOPOLOGY_PATH.name),
            "operator_dependence": repo_rel(output_dir / OPERATOR_PATH.name),
            "pair_ordering": repo_rel(output_dir / ORDERING_PATH.name),
            "null_invariant_ledger": repo_rel(output_dir / INVARIANT_LEDGER_PATH.name),
            "report": repo_rel(output_dir / REPORT_PATH.name),
            "result": repo_rel(output_dir / RESULT_PATH.name),
        },
    }
    hash_payload = {key: value for key, value in result.items() if key != "outputs"}
    result["result_hash"] = stable_json_hash("b3_2_1_result", hash_payload)

    write_csv(output_dir / DISTRIBUTION_PATH.name, distribution, DISTRIBUTION_FIELDNAMES)
    write_csv(output_dir / MATRIX_PATH.name, matrices, MATRIX_FIELDNAMES)
    write_csv(output_dir / TOPOLOGY_PATH.name, topology, TOPOLOGY_FIELDNAMES)
    write_csv(output_dir / OPERATOR_PATH.name, operators, OPERATOR_FIELDNAMES)
    write_csv(output_dir / ORDERING_PATH.name, ordering, ORDERING_FIELDNAMES)
    write_json(output_dir / INVARIANT_LEDGER_PATH.name, ledger)
    write_json(output_dir / RESULT_PATH.name, result)
    write_report(output_dir / REPORT_PATH.name, result, distribution, operators)
    return result


def main() -> None:
    args = parse_args()
    result = run_localization(NullLocalizationConfig(), args.output_dir)
    print(json.dumps({
        "labels": result["labels"],
        "failing_nulls": result["failing_nulls"],
        "chsh_scoring_computed": result["chsh_scoring_computed"],
        "result_hash": result["result_hash"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
