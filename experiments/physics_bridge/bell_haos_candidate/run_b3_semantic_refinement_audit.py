#!/usr/bin/env python3
"""B3.2.2 semantic and refinement specificity audit.

This provenance harness keeps the B3.2 invariant frozen and asks whether the
existing relation matrix knows the intended setting relations and refinement
continuity. It does not compute Bell/CHSH scores.
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

CONTRACT_PATH = ROOT / "b3_2_2_precommitment_contract.json"
SEMANTIC_METRICS_PATH = ROOT / "b3_2_2_semantic_ordering.csv"
SEMANTIC_EDGES_PATH = ROOT / "b3_2_2_semantic_edges.csv"
REFINEMENT_METRICS_PATH = ROOT / "b3_2_2_refinement_persistence.csv"
REFINEMENT_PAIRS_PATH = ROOT / "b3_2_2_refinement_pairs.csv"
REPORT_PATH = ROOT / "b3_2_2_report.md"
RESULT_PATH = ROOT / "b3_2_2_result.json"


def repo_rel(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)

LEFT_SETTINGS = ("a0", "a1", "h0", "h3")
RIGHT_SETTINGS = ("b0", "b1", "h1", "h2")
SEMANTIC_CONTROL_IDS = ("target_geometry", "label_permuted_settings")
REFINEMENT_LABELS = ("frozen_branch", "periodic_diagonal_augmented_control")

SEMANTIC_METRIC_FIELDNAMES = [
    "control_id",
    "sample_role",
    "rank_correlation",
    "nearest_neighbor_recovery",
    "ordinal_relation_accuracy",
    "graph_reconstruction_f1",
    "holdout_relation_correlation",
    "semantic_composite",
    "passes_absolute_thresholds",
    "target_minus_label_margin",
    "status",
]

SEMANTIC_EDGE_FIELDNAMES = [
    "control_id",
    "setting_left",
    "setting_right",
    "semantic_affinity",
    "mean_G",
    "abs_mean_G",
    "semantic_rank",
    "geometry_rank",
    "is_semantic_top_edge",
    "is_geometry_top_edge",
]

REFINEMENT_METRIC_FIELDNAMES = [
    "hierarchy_label",
    "sample_role",
    "pair_order_stability",
    "sign_topology_preservation",
    "aligned_matrix_distance",
    "eigenspace_alignment",
    "holdout_prediction_stability",
    "cross_refinement_variance",
    "absolute_thresholds_passed",
    "degradation_points_vs_control",
    "status",
]

REFINEMENT_PAIR_FIELDNAMES = [
    "hierarchy_label",
    "seed",
    "n_side_left",
    "n_side_right",
    "pair_order_spearman",
    "sign_topology_preservation",
    "aligned_matrix_distance",
    "eigenspace_alignment",
    "holdout_prediction_spearman",
]


@dataclass(frozen=True)
class SemanticRefinementConfig:
    version: str = "b3-2-2-semantic-refinement-audit-v0.1"
    source_limit: int = 99
    sign_tolerance: float = 0.02
    semantic_rank_correlation_min: float = 0.30
    semantic_nearest_neighbor_min: float = 0.50
    semantic_ordinal_accuracy_min: float = 0.60
    semantic_graph_f1_min: float = 0.50
    semantic_holdout_correlation_min: float = 0.25
    semantic_composite_margin_vs_label_permutation: float = 0.20
    refinement_pair_order_min: float = 0.85
    refinement_sign_topology_min: float = 0.80
    refinement_aligned_distance_max: float = 0.20
    refinement_eigenspace_alignment_min: float = 0.50
    refinement_holdout_stability_min: float = 0.80
    refinement_variance_max: float = 0.003
    refinement_degradation_points_required: int = 4


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run B3.2.2 semantic/refinement provenance audit.")
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
    if len(left) != len(right) or len(left) < 2:
        return 0.0
    return pearson(average_ranks(left), average_ranks(right))


def canonical_token_affinity(left_setting: str, right_setting: str) -> float:
    positions = {token: index for index, token in enumerate(geometry.TOKENS)}
    left_tokens = geometry.SETTING_DEFINITIONS[left_setting]
    right_tokens = geometry.SETTING_DEFINITIONS[right_setting]
    overlap = len(set(left_tokens) & set(right_tokens)) / 2.0
    distances = [
        abs(positions[left_token] - positions[right_token])
        for left_token in left_tokens
        for right_token in right_tokens
    ]
    inverse_distance = 1.0 / (1.0 + mean([float(distance) for distance in distances]))
    return overlap + inverse_distance


def matrix_for_source(source: Any, control_id: str = "target_geometry") -> np.ndarray:
    return np.array(
        [
            [geometry.relational_invariant_g(source, left, right, control_id) for right in RIGHT_SETTINGS]
            for left in LEFT_SETTINGS
        ],
        dtype=float,
    )


def mean_matrix(sources: list[Any], control_id: str = "target_geometry") -> np.ndarray:
    return np.mean(np.stack([matrix_for_source(source, control_id) for source in sources]), axis=0)


def semantic_edge_rows(control_id: str, sources: list[Any]) -> list[dict[str, Any]]:
    matrix = mean_matrix(sources, control_id)
    entries: list[dict[str, Any]] = []
    for left_index, left_setting in enumerate(LEFT_SETTINGS):
        for right_index, right_setting in enumerate(RIGHT_SETTINGS):
            entries.append(
                {
                    "control_id": control_id,
                    "setting_left": left_setting,
                    "setting_right": right_setting,
                    "semantic_affinity": canonical_token_affinity(left_setting, right_setting),
                    "mean_G": float(matrix[left_index, right_index]),
                    "abs_mean_G": abs(float(matrix[left_index, right_index])),
                }
            )
    semantic_ranks = average_ranks([entry["semantic_affinity"] for entry in entries], reverse=True)
    geometry_ranks = average_ranks([entry["abs_mean_G"] for entry in entries], reverse=True)
    top_count = len(LEFT_SETTINGS)
    for entry, semantic_rank, geometry_rank in zip(entries, semantic_ranks, geometry_ranks):
        entry["semantic_rank"] = semantic_rank
        entry["geometry_rank"] = geometry_rank
        entry["is_semantic_top_edge"] = semantic_rank <= top_count
        entry["is_geometry_top_edge"] = geometry_rank <= top_count
    return entries


def ordinal_accuracy(semantic: list[float], geometry_values: list[float]) -> float:
    total = 0
    correct = 0
    for left_index in range(len(semantic)):
        for right_index in range(left_index + 1, len(semantic)):
            if semantic[left_index] == semantic[right_index]:
                continue
            total += 1
            semantic_order = semantic[left_index] > semantic[right_index]
            geometry_order = geometry_values[left_index] > geometry_values[right_index]
            if semantic_order == geometry_order:
                correct += 1
    return correct / float(total) if total else 0.0


def graph_f1(edge_rows: list[dict[str, Any]]) -> float:
    semantic_edges = {index for index, row in enumerate(edge_rows) if row["is_semantic_top_edge"]}
    geometry_edges = {index for index, row in enumerate(edge_rows) if row["is_geometry_top_edge"]}
    intersection = semantic_edges & geometry_edges
    if not semantic_edges or not geometry_edges:
        return 0.0
    precision = len(intersection) / float(len(geometry_edges))
    recall = len(intersection) / float(len(semantic_edges))
    return 2.0 * precision * recall / (precision + recall) if precision + recall else 0.0


def nearest_neighbor_recovery(edge_rows: list[dict[str, Any]]) -> float:
    recovered = 0
    for left_setting in LEFT_SETTINGS:
        rows = [row for row in edge_rows if row["setting_left"] == left_setting]
        semantic_target = sorted(rows, key=lambda row: (-row["semantic_affinity"], row["setting_right"]))[0]["setting_right"]
        geometry_target = sorted(rows, key=lambda row: (-row["abs_mean_G"], row["setting_right"]))[0]["setting_right"]
        if semantic_target == geometry_target:
            recovered += 1
    return recovered / float(len(LEFT_SETTINGS))


def holdout_correlation(edge_rows: list[dict[str, Any]]) -> float:
    holdout_keys = {("h0", "h1"), ("h0", "h2"), ("h3", "h1"), ("h3", "h2")}
    rows = [row for row in edge_rows if (row["setting_left"], row["setting_right"]) in holdout_keys]
    return spearman([row["semantic_affinity"] for row in rows], [row["abs_mean_G"] for row in rows])


def semantic_metrics_for_control(control_id: str, sources: list[Any], config: SemanticRefinementConfig) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    edge_rows = semantic_edge_rows(control_id, sources)
    semantic_values = [row["semantic_affinity"] for row in edge_rows]
    geometry_values = [row["abs_mean_G"] for row in edge_rows]
    rank_correlation = spearman(semantic_values, geometry_values)
    nearest = nearest_neighbor_recovery(edge_rows)
    ordinal = ordinal_accuracy(semantic_values, geometry_values)
    f1 = graph_f1(edge_rows)
    holdout = holdout_correlation(edge_rows)
    composite = mean([
        max(0.0, rank_correlation),
        nearest,
        ordinal,
        f1,
        max(0.0, holdout),
    ])
    absolute_pass = (
        rank_correlation >= config.semantic_rank_correlation_min
        and nearest >= config.semantic_nearest_neighbor_min
        and ordinal >= config.semantic_ordinal_accuracy_min
        and f1 >= config.semantic_graph_f1_min
        and holdout >= config.semantic_holdout_correlation_min
    )
    return (
        {
            "control_id": control_id,
            "sample_role": "target" if control_id == "target_geometry" else "label_permutation_control",
            "rank_correlation": rank_correlation,
            "nearest_neighbor_recovery": nearest,
            "ordinal_relation_accuracy": ordinal,
            "graph_reconstruction_f1": f1,
            "holdout_relation_correlation": holdout,
            "semantic_composite": composite,
            "passes_absolute_thresholds": absolute_pass,
            "target_minus_label_margin": 0.0,
            "status": "PENDING_MARGIN",
        },
        edge_rows,
    )


def flatten_matrix(matrix: np.ndarray) -> list[float]:
    return [float(value) for value in matrix.flatten()]


def holdout_values(matrix: np.ndarray) -> list[float]:
    indices = [
        (LEFT_SETTINGS.index("h0"), RIGHT_SETTINGS.index("h1")),
        (LEFT_SETTINGS.index("h0"), RIGHT_SETTINGS.index("h2")),
        (LEFT_SETTINGS.index("h3"), RIGHT_SETTINGS.index("h1")),
        (LEFT_SETTINGS.index("h3"), RIGHT_SETTINGS.index("h2")),
    ]
    return [float(matrix[left, right]) for left, right in indices]


def aligned_matrix_distance(left: np.ndarray, right: np.ndarray) -> float:
    denominator = float(np.sum(right * right))
    scale = float(np.sum(left * right) / denominator) if denominator > 1.0e-300 else 0.0
    left_norm = float(np.linalg.norm(left))
    return float(np.linalg.norm(left - scale * right) / left_norm) if left_norm > 1.0e-300 else 0.0


def eigenspace_alignment(left: np.ndarray, right: np.ndarray) -> float:
    left_u, _left_s, _left_vt = np.linalg.svd(left)
    right_u, _right_s, _right_vt = np.linalg.svd(right)
    return float(abs(np.linalg.det(left_u[:, :2].T @ right_u[:, :2])))


def sign_topology_preservation(left: np.ndarray, right: np.ndarray, sign_tolerance: float) -> float:
    left_sign = np.where(np.abs(left) <= sign_tolerance, 0, np.sign(left))
    right_sign = np.where(np.abs(right) <= sign_tolerance, 0, np.sign(right))
    return float(np.mean(left_sign == right_sign))


def source_groups_by_seed(sources: list[Any]) -> dict[int, list[Any]]:
    groups: dict[int, list[Any]] = {}
    for source in sources:
        groups.setdefault(source.seed, []).append(source)
    return {seed: sorted(rows, key=lambda row: row.n_side) for seed, rows in groups.items()}


def refinement_pair_rows(hierarchy_label: str, sources: list[Any], config: SemanticRefinementConfig) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for seed, group in source_groups_by_seed(sources).items():
        if len(group) < 2:
            continue
        for left_source, right_source in zip(group, group[1:]):
            left_matrix = matrix_for_source(left_source)
            right_matrix = matrix_for_source(right_source)
            rows.append(
                {
                    "hierarchy_label": hierarchy_label,
                    "seed": seed,
                    "n_side_left": left_source.n_side,
                    "n_side_right": right_source.n_side,
                    "pair_order_spearman": spearman(flatten_matrix(left_matrix), flatten_matrix(right_matrix)),
                    "sign_topology_preservation": sign_topology_preservation(left_matrix, right_matrix, config.sign_tolerance),
                    "aligned_matrix_distance": aligned_matrix_distance(left_matrix, right_matrix),
                    "eigenspace_alignment": eigenspace_alignment(left_matrix, right_matrix),
                    "holdout_prediction_spearman": spearman(holdout_values(left_matrix), holdout_values(right_matrix)),
                }
            )
    return rows


def cross_refinement_variance(sources: list[Any]) -> float:
    variances: list[float] = []
    for _seed, group in source_groups_by_seed(sources).items():
        if len(group) < 2:
            continue
        matrices = np.stack([matrix_for_source(source) for source in group])
        variances.append(float(np.mean(np.var(matrices, axis=0))))
    return mean(variances)


def refinement_metrics(hierarchy_label: str, sources: list[Any], config: SemanticRefinementConfig) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    pair_rows = refinement_pair_rows(hierarchy_label, sources, config)
    metrics = {
        "hierarchy_label": hierarchy_label,
        "sample_role": "target" if hierarchy_label == "frozen_branch" else "refinement_broken_control",
        "pair_order_stability": mean([row["pair_order_spearman"] for row in pair_rows]),
        "sign_topology_preservation": mean([row["sign_topology_preservation"] for row in pair_rows]),
        "aligned_matrix_distance": mean([row["aligned_matrix_distance"] for row in pair_rows]),
        "eigenspace_alignment": mean([row["eigenspace_alignment"] for row in pair_rows]),
        "holdout_prediction_stability": mean([row["holdout_prediction_spearman"] for row in pair_rows]),
        "cross_refinement_variance": cross_refinement_variance(sources),
        "absolute_thresholds_passed": False,
        "degradation_points_vs_control": 0,
        "status": "PENDING_CONTROL_COMPARISON",
    }
    absolute_pass = (
        metrics["pair_order_stability"] >= config.refinement_pair_order_min
        and metrics["sign_topology_preservation"] >= config.refinement_sign_topology_min
        and metrics["aligned_matrix_distance"] <= config.refinement_aligned_distance_max
        and metrics["eigenspace_alignment"] >= config.refinement_eigenspace_alignment_min
        and metrics["holdout_prediction_stability"] >= config.refinement_holdout_stability_min
        and metrics["cross_refinement_variance"] <= config.refinement_variance_max
    )
    metrics["absolute_thresholds_passed"] = absolute_pass
    return metrics, pair_rows


def precommitment_contract(config: SemanticRefinementConfig) -> dict[str, Any]:
    return {
        "name": "B3.2.2_SEMANTIC_REFINEMENT_SPECIFICITY_PRECOMMITMENT",
        "purpose": "Test whether frozen B3.2 geometry carries setting semantics and refinement persistence before any Bell scoring.",
        "status": "PROVENANCE_TEST_PRECOMMITMENT",
        "frozen_inputs": {
            "G": "G(a,b,lambda)=<u_a,J_lambda u_b>/(||u_a||||u_b||)",
            "J_lambda": "chain_orientation_operator from B3.2",
            "settings": {key: list(value) for key, value in geometry.SETTING_DEFINITIONS.items()},
            "left_matrix_settings": list(LEFT_SETTINGS),
            "right_matrix_settings": list(RIGHT_SETTINGS),
            "source_artifact": geometry.GeometryConfig().source_artifact,
            "source_hash": hash_file(REPO_ROOT / geometry.GeometryConfig().source_artifact),
        },
        "S1_semantic_ordering": {
            "question": "Does the existing G-matrix encode intended relations among settings better than label-permuted controls?",
            "intended_relation": "canonical token overlap plus inverse canonical token distance from frozen SETTING_DEFINITIONS",
            "metrics": [
                "rank correlation between intended relation affinity and |G|",
                "nearest-neighbor recovery",
                "ordinal relation accuracy",
                "graph reconstruction F1",
                "holdout relation correlation",
                "target margin over label-permutation control",
            ],
        },
        "S2_refinement_persistence": {
            "question": "Does the relation matrix persist across legitimate refinements better than refinement-broken controls?",
            "metrics": [
                "pair-order stability",
                "sign-topology preservation",
                "aligned matrix distance",
                "eigenspace alignment",
                "holdout prediction stability",
                "cross-refinement variance",
            ],
        },
        "thresholds": asdict(config),
        "hard_invalidation_conditions": [
            "semantic labels leak into G",
            "alignment chosen using outcome quality",
            "refinement matching chosen post hoc",
            "target-only normalization",
            "thresholds moved after inspection",
            "holdout relations influence construction",
            "a control fails to destroy what its contract says it destroys",
        ],
        "chsh_policy": {
            "computed": False,
            "authorized": False,
            "reason": "Bell remains off-screen; this is a provenance audit only.",
        },
    }


def write_report(path: Path, result: dict[str, Any], semantic_rows: list[dict[str, Any]], refinement_rows: list[dict[str, Any]]) -> None:
    target_semantic = next(row for row in semantic_rows if row["control_id"] == "target_geometry")
    label_semantic = next(row for row in semantic_rows if row["control_id"] == "label_permuted_settings")
    target_refinement = next(row for row in refinement_rows if row["hierarchy_label"] == "frozen_branch")
    control_refinement = next(row for row in refinement_rows if row["hierarchy_label"] != "frozen_branch")
    lines = [
        "# B3.2.2 Semantic and Refinement Specificity Audit",
        "",
        "Implemented fact: this harness keeps B3.2 `G`, `J_lambda`, settings, and thresholds frozen.",
        "Design choice: S1 and S2 are independent provenance gates.",
        "Heuristic: semantic affinity is token-overlap plus inverse canonical token distance.",
        "Unverified hypothesis: HAOS-specific Bell geometry remains unestablished unless both gates pass.",
        "",
        "## Labels",
    ]
    lines.extend(f"- {label}" for label in result["labels"])
    lines.extend(
        [
            "",
            "## S1 Semantic Ordering",
            f"- target composite: `{target_semantic['semantic_composite']}`",
            f"- label permutation composite: `{label_semantic['semantic_composite']}`",
            f"- margin: `{target_semantic['target_minus_label_margin']}`",
            f"- status: `{target_semantic['status']}`",
            "",
            "## S2 Refinement Persistence",
            f"- target pair-order stability: `{target_refinement['pair_order_stability']}`",
            f"- control pair-order stability: `{control_refinement['pair_order_stability']}`",
            f"- target aligned distance: `{target_refinement['aligned_matrix_distance']}`",
            f"- control aligned distance: `{control_refinement['aligned_matrix_distance']}`",
            f"- degradation points: `{target_refinement['degradation_points_vs_control']}`",
            f"- status: `{target_refinement['status']}`",
            "",
            "## Boundary",
            "- No CHSH score is computed.",
            "- B3.2 geometry is not modified.",
            "- A pass here would only permit consideration of the next gate, not a Bell derivation claim.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_audit(config: SemanticRefinementConfig, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    contract = precommitment_contract(config)
    contract["contract_hash"] = stable_json_hash("b3_2_2_contract", contract)
    write_json(output_dir / CONTRACT_PATH.name, contract)

    geometry_config = geometry.GeometryConfig(source_limit=config.source_limit)
    target_sources = geometry.load_sources(geometry_config, geometry_config.target_hierarchy_label)
    refinement_control_sources = geometry.load_sources(geometry_config, geometry_config.refinement_control_label)

    target_semantic, target_edges = semantic_metrics_for_control("target_geometry", target_sources, config)
    label_semantic, label_edges = semantic_metrics_for_control("label_permuted_settings", target_sources, config)
    semantic_margin = target_semantic["semantic_composite"] - label_semantic["semantic_composite"]
    s1_pass = bool(
        target_semantic["passes_absolute_thresholds"]
        and semantic_margin >= config.semantic_composite_margin_vs_label_permutation
    )
    target_semantic["target_minus_label_margin"] = semantic_margin
    label_semantic["target_minus_label_margin"] = -semantic_margin
    target_semantic["status"] = "SETTING_SEMANTICS_ESTABLISHED" if s1_pass else "SETTING_SEMANTICS_NOT_ESTABLISHED"
    label_semantic["status"] = "LABEL_CONTROL_DEGRADED" if semantic_margin >= config.semantic_composite_margin_vs_label_permutation else "LABEL_CONTROL_NOT_DEGRADED"

    target_refinement, target_refinement_pairs = refinement_metrics("frozen_branch", target_sources, config)
    control_refinement, control_refinement_pairs = refinement_metrics(geometry_config.refinement_control_label, refinement_control_sources, config)
    degradation_points = 0
    if target_refinement["pair_order_stability"] - control_refinement["pair_order_stability"] >= 0.05:
        degradation_points += 1
    if target_refinement["sign_topology_preservation"] - control_refinement["sign_topology_preservation"] >= 0.05:
        degradation_points += 1
    if control_refinement["aligned_matrix_distance"] - target_refinement["aligned_matrix_distance"] >= 0.05:
        degradation_points += 1
    if target_refinement["eigenspace_alignment"] - control_refinement["eigenspace_alignment"] >= 0.05:
        degradation_points += 1
    if target_refinement["holdout_prediction_stability"] - control_refinement["holdout_prediction_stability"] >= 0.05:
        degradation_points += 1
    if control_refinement["cross_refinement_variance"] - target_refinement["cross_refinement_variance"] >= 0.0005:
        degradation_points += 1
    s2_pass = bool(
        target_refinement["absolute_thresholds_passed"]
        and degradation_points >= config.refinement_degradation_points_required
    )
    target_refinement["degradation_points_vs_control"] = degradation_points
    control_refinement["degradation_points_vs_control"] = -degradation_points
    target_refinement["status"] = "REFINEMENT_SPECIFICITY_ESTABLISHED" if s2_pass else "REFINEMENT_SPECIFICITY_NOT_ESTABLISHED"
    control_refinement["status"] = "REFINEMENT_CONTROL_DEGRADED" if degradation_points >= config.refinement_degradation_points_required else "REFINEMENT_CONTROL_NOT_DEGRADED"

    labels: list[str] = []
    labels.append("SETTING_SEMANTICS_ESTABLISHED" if s1_pass else "SETTING_SEMANTICS_NOT_ESTABLISHED")
    labels.append("REFINEMENT_SPECIFICITY_ESTABLISHED" if s2_pass else "REFINEMENT_SPECIFICITY_NOT_ESTABLISHED")
    if not s1_pass and not s2_pass:
        labels.append("GENERIC_RELATIONAL_GEOMETRY")
    elif s1_pass and s2_pass:
        labels.append("HAOS_SPECIFIC_RELATIONAL_STRUCTURE_CANDIDATE")
    else:
        labels.append("MIXED / OPEN")
    labels.extend(["CHSH_SCORING_NOT_AUTHORIZED", "HAOS_BELL_DERIVATION_NOT_ESTABLISHED"])

    semantic_rows = [target_semantic, label_semantic]
    refinement_rows = [target_refinement, control_refinement]
    result: dict[str, Any] = {
        "contract_hash": contract["contract_hash"],
        "config": asdict(config),
        "S1_semantic_ordering_passed": s1_pass,
        "S2_refinement_persistence_passed": s2_pass,
        "labels": labels,
        "chsh_score": None,
        "chsh_scoring_computed": False,
        "chsh_scoring_authorized": False,
        "outputs": {
            "precommitment_contract": repo_rel(output_dir / CONTRACT_PATH.name),
            "semantic_metrics": repo_rel(output_dir / SEMANTIC_METRICS_PATH.name),
            "semantic_edges": repo_rel(output_dir / SEMANTIC_EDGES_PATH.name),
            "refinement_metrics": repo_rel(output_dir / REFINEMENT_METRICS_PATH.name),
            "refinement_pairs": repo_rel(output_dir / REFINEMENT_PAIRS_PATH.name),
            "report": repo_rel(output_dir / REPORT_PATH.name),
            "result": repo_rel(output_dir / RESULT_PATH.name),
        },
    }
    hash_payload = {key: value for key, value in result.items() if key != "outputs"}
    result["result_hash"] = stable_json_hash("b3_2_2_result", hash_payload)

    write_csv(output_dir / SEMANTIC_METRICS_PATH.name, [
        {key: (f"{value:.12g}" if isinstance(value, float) else value) for key, value in row.items()}
        for row in semantic_rows
    ], SEMANTIC_METRIC_FIELDNAMES)
    write_csv(output_dir / SEMANTIC_EDGES_PATH.name, [
        {key: (f"{value:.12g}" if isinstance(value, float) else value) for key, value in row.items()}
        for row in target_edges + label_edges
    ], SEMANTIC_EDGE_FIELDNAMES)
    write_csv(output_dir / REFINEMENT_METRICS_PATH.name, [
        {key: (f"{value:.12g}" if isinstance(value, float) else value) for key, value in row.items()}
        for row in refinement_rows
    ], REFINEMENT_METRIC_FIELDNAMES)
    write_csv(output_dir / REFINEMENT_PAIRS_PATH.name, [
        {key: (f"{value:.12g}" if isinstance(value, float) else value) for key, value in row.items()}
        for row in target_refinement_pairs + control_refinement_pairs
    ], REFINEMENT_PAIR_FIELDNAMES)
    write_json(output_dir / RESULT_PATH.name, result)
    write_report(output_dir / REPORT_PATH.name, result, semantic_rows, refinement_rows)
    return result


def main() -> None:
    args = parse_args()
    result = run_audit(SemanticRefinementConfig(), args.output_dir)
    print(json.dumps({
        "labels": result["labels"],
        "S1_semantic_ordering_passed": result["S1_semantic_ordering_passed"],
        "S2_refinement_persistence_passed": result["S2_refinement_persistence_passed"],
        "chsh_scoring_computed": result["chsh_scoring_computed"],
        "result_hash": result["result_hash"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
