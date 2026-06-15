#!/usr/bin/env python3
"""B3.2 HAOS-IIP relational geometry audit.

This harness tests whether frozen HAOS-IIP structure yields a relational,
sign-changing invariant G(a,b,lambda) before any Bell score is computed.

It does not run CHSH scoring. Passing G1-G5 only authorizes a later, separate
scoreboard run under a new precommitment.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import inspect
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable


REPO_ROOT = Path(__file__).resolve().parents[3]
ROOT = Path(__file__).resolve().parent
PHASE18_SOURCE = REPO_ROOT / "phase18-distance-surrogate" / "runs" / "phase18_shell_ordering_metrics.csv"

CONTRACT_PATH = ROOT / "b3_2_precommitment_contract.json"
GEOMETRY_MATRIX_PATH = ROOT / "b3_2_geometry_matrix.csv"
COVARIANCE_PATH = ROOT / "b3_2_covariance_diagnostics.csv"
CONTROL_RESULTS_PATH = ROOT / "b3_2_control_results.csv"
REPORT_PATH = ROOT / "b3_2_geometry_report.md"
RESULT_PATH = ROOT / "b3_2_result.json"


def repo_rel(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)

TOKENS = (
    "low_k_half",
    "dispersion_half",
    "radius_half",
    "spectral_half",
    "width_half",
    "front_near",
    "front_far",
)

SETTING_DEFINITIONS = {
    "a0": ("radius_half", "spectral_half"),
    "a1": ("width_half", "front_near"),
    "b0": ("dispersion_half", "spectral_half"),
    "b1": ("low_k_half", "front_far"),
    "h0": ("radius_half", "front_near"),
    "h1": ("dispersion_half", "width_half"),
    "h2": ("spectral_half", "front_far"),
    "h3": ("low_k_half", "radius_half"),
}

CHSH_HELD_BACK_PAIRS = (("a0", "b0"), ("a0", "b1"), ("a1", "b0"), ("a1", "b1"))
HOLDOUT_PAIRS = (("h0", "h1"), ("h0", "h2"), ("h3", "h1"), ("h3", "h2"))
ALL_PAIR_ROLES = {"chsh_held_back": CHSH_HELD_BACK_PAIRS, "holdout": HOLDOUT_PAIRS}

LABEL_PERMUTATION = {
    "a0": "a1",
    "a1": "h0",
    "b0": "b1",
    "b1": "h1",
    "h0": "h2",
    "h1": "h3",
    "h2": "a0",
    "h3": "b0",
}

INDEX_PERMUTATION = (2, 4, 6, 1, 3, 5, 0)
SIGNED_PERMUTATION = ((3, 1), (5, -1), (0, 1), (6, -1), (2, 1), (4, -1), (1, 1))

FORBIDDEN_GEOMETRY_LEAKAGE_PATTERNS = (
    "sqrt(2)",
    "2*sqrt(2)",
    "2 * sqrt(2)",
    "-cos(",
    "cos(theta",
    "cos(2*delta",
    "cos(2 * delta",
    "quantum_expectations",
    "target_chsh_value",
    "target_s_value",
    "reference_outcomes",
    "adaptive_optimization",
)

GEOMETRY_FIELDNAMES = [
    "candidate_id",
    "control_id",
    "pair_role",
    "source_id",
    "hierarchy_label",
    "seed",
    "n_side",
    "setting_left",
    "setting_right",
    "G",
    "operator_kind",
    "vector_kind",
    "closure_strength",
    "chain_signature",
]

COVARIANCE_FIELDNAMES = [
    "source_id",
    "pair_role",
    "setting_left",
    "setting_right",
    "G_original",
    "G_transformed",
    "absolute_error",
    "tolerance",
    "status",
]

CONTROL_FIELDNAMES = [
    "control_id",
    "purpose",
    "preserves",
    "destroys",
    "expected_metric_response",
    "invalidation_condition",
    "count",
    "mean_G",
    "min_G",
    "max_G",
    "range_G",
    "variance_G",
    "positive_count",
    "negative_count",
    "geometry_score",
    "degradation_ratio",
    "status",
]


@dataclass(frozen=True)
class GeometryConfig:
    version: str = "b3-2-geometry-audit-v0.1"
    source_artifact: str = "phase18-distance-surrogate/runs/phase18_shell_ordering_metrics.csv"
    target_hierarchy_label: str = "frozen_branch"
    refinement_control_label: str = "periodic_diagonal_augmented_control"
    source_limit: int = 6
    relational_range_tolerance: float = 0.05
    sign_tolerance: float = 0.02
    covariance_tolerance: float = 1.0e-12
    holdout_range_tolerance: float = 0.05
    null_degradation_ratio: float = 0.70
    seeds: tuple[int, ...] = (7201, 7202, 7203)


@dataclass(frozen=True)
class SourceRow:
    source_id: str
    source_hash: str
    hierarchy_label: str
    seed: int
    n_side: int
    fields: dict[str, Any]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run B3.2 relational geometry audit.")
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


def detect_geometry_leakage(source_texts: Iterable[str]) -> tuple[bool, list[str]]:
    hits: list[str] = []
    normalized = [text.replace(" ", "").lower() for text in source_texts]
    for pattern in FORBIDDEN_GEOMETRY_LEAKAGE_PATTERNS:
        needle = pattern.replace(" ", "").lower()
        if any(needle in text for text in normalized):
            hits.append(pattern)
    return bool(hits), sorted(set(hits))


def generator_sources() -> list[str]:
    functions = [
        setting_vector,
        chain_orientation_operator,
        relational_invariant_g,
        evaluate_geometry_rows,
        covariance_rows,
    ]
    return [inspect.getsource(function) for function in functions]


def load_sources(config: GeometryConfig, hierarchy_label: str) -> list[SourceRow]:
    source_path = REPO_ROOT / config.source_artifact
    source_hash = hash_file(source_path)
    rows: list[SourceRow] = []
    with source_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row.get("hierarchy_label") != hierarchy_label:
                continue
            fields = {
                "hierarchy_label": row["hierarchy_label"],
                "n_side": int(row["n_side"]),
                "seed": int(row["seed"]),
                "probe_name": row["probe_name"],
                "shell_order_score": float(row["shell_order_score"]),
                "max_shell_overlap_fraction": float(row["max_shell_overlap_fraction"]),
                "near_mean_arrival": float(row["near_mean_arrival"]),
                "far_mean_arrival": float(row["far_mean_arrival"]),
                "phase17_mismatch_rate": float(row["phase17_mismatch_rate"]),
                "chain_signature": row["chain_signature"],
            }
            source_id = stable_json_hash("phase18_geometry_source", fields)
            rows.append(
                SourceRow(
                    source_id=source_id,
                    source_hash=source_hash,
                    hierarchy_label=hierarchy_label,
                    seed=int(row["seed"]),
                    n_side=int(row["n_side"]),
                    fields=fields,
                )
            )
            if len(rows) >= config.source_limit:
                break
    if not rows:
        raise ValueError(f"no rows found for hierarchy label {hierarchy_label}")
    return rows


def clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


def source_closure_strength(source: SourceRow) -> float:
    shell = float(source.fields["shell_order_score"])
    overlap = float(source.fields["max_shell_overlap_fraction"])
    mismatch = float(source.fields["phase17_mismatch_rate"])
    near = float(source.fields["near_mean_arrival"])
    far = float(source.fields["far_mean_arrival"])
    arrival_gap = clamp((far - near) / max(far, 1.0e-12), 0.0, 1.0)
    return clamp(shell * (1.0 - overlap) * (1.0 - mismatch) + 0.25 * arrival_gap, 0.05, 0.95)


def rank_weights(source: SourceRow) -> dict[str, float]:
    chain = str(source.fields["chain_signature"]).split(">")
    if len(chain) <= 1:
        return {token: 1.0 for token in TOKENS}
    ranks = {token: index / float(len(chain) - 1) for index, token in enumerate(chain)}
    return {token: 1.0 + (0.5 - ranks.get(token, 0.5)) for token in TOKENS}


def setting_vector(source: SourceRow, setting: str, vector_kind: str = "target") -> list[float]:
    mapped_setting = LABEL_PERMUTATION[setting] if vector_kind == "label_permuted" else setting
    plus_token, minus_token = SETTING_DEFINITIONS[mapped_setting]
    weights = rank_weights(source)
    vector = [0.0 for _ in TOKENS]
    vector[TOKENS.index(plus_token)] = weights[plus_token]
    vector[TOKENS.index(minus_token)] = -weights[minus_token]
    if vector_kind == "closure_vector_permuted":
        vector = permute_vector(vector, INDEX_PERMUTATION)
    return vector


def zero_matrix(size: int) -> list[list[float]]:
    return [[0.0 for _ in range(size)] for _ in range(size)]


def identity_matrix(size: int) -> list[list[float]]:
    matrix = zero_matrix(size)
    for index in range(size):
        matrix[index][index] = 1.0
    return matrix


def chain_orientation_operator(source: SourceRow) -> list[list[float]]:
    index_by_token = {token: index for index, token in enumerate(TOKENS)}
    matrix = zero_matrix(len(TOKENS))
    chain = str(source.fields["chain_signature"]).split(">")
    strength = source_closure_strength(source)
    for edge_index, (left, right) in enumerate(zip(chain, chain[1:])):
        if left not in index_by_token or right not in index_by_token:
            continue
        left_index = index_by_token[left]
        right_index = index_by_token[right]
        weight = strength / float(edge_index + 1)
        matrix[left_index][right_index] += weight
        matrix[right_index][left_index] -= weight
    return matrix


def permute_vector(vector: list[float], permutation: tuple[int, ...]) -> list[float]:
    output = [0.0 for _ in vector]
    for old_index, new_index in enumerate(permutation):
        output[new_index] = vector[old_index]
    return output


def signed_permute_vector(vector: list[float], permutation: tuple[tuple[int, int], ...]) -> list[float]:
    output = [0.0 for _ in vector]
    for old_index, (new_index, sign) in enumerate(permutation):
        output[new_index] = sign * vector[old_index]
    return output


def transform_matrix(matrix: list[list[float]], permutation: tuple[int, ...]) -> list[list[float]]:
    size = len(matrix)
    output = zero_matrix(size)
    for old_i, new_i in enumerate(permutation):
        for old_j, new_j in enumerate(permutation):
            output[new_i][new_j] = matrix[old_i][old_j]
    return output


def signed_transform_matrix(matrix: list[list[float]], permutation: tuple[tuple[int, int], ...]) -> list[list[float]]:
    size = len(matrix)
    output = zero_matrix(size)
    for old_i, (new_i, sign_i) in enumerate(permutation):
        for old_j, (new_j, sign_j) in enumerate(permutation):
            output[new_i][new_j] = sign_i * sign_j * matrix[old_i][old_j]
    return output


def operator_for_control(source: SourceRow, control_id: str) -> list[list[float]]:
    target = chain_orientation_operator(source)
    if control_id in {"target_geometry", "label_permuted_settings", "closure_vector_permutation", "refinement_broken_control"}:
        return target
    if control_id == "identity_pairing":
        return identity_matrix(len(TOKENS))
    if control_id == "shuffled_pairing":
        return transform_matrix(target, INDEX_PERMUTATION)
    if control_id == "random_orthogonal_pairing":
        return signed_transform_matrix(target, SIGNED_PERMUTATION)
    if control_id == "topology_destroyed_pairing":
        return zero_matrix(len(TOKENS))
    raise ValueError(f"unknown control id {control_id}")


def dot(left: list[float], right: list[float]) -> float:
    return sum(a * b for a, b in zip(left, right))


def matvec(matrix: list[list[float]], vector: list[float]) -> list[float]:
    return [sum(row[index] * vector[index] for index in range(len(vector))) for row in matrix]


def norm(vector: list[float]) -> float:
    return math.sqrt(sum(value * value for value in vector))


def relational_invariant_g(
    source: SourceRow,
    setting_left: str,
    setting_right: str,
    control_id: str = "target_geometry",
) -> float:
    vector_kind = "target"
    if control_id == "label_permuted_settings":
        vector_kind = "label_permuted"
    elif control_id == "closure_vector_permutation":
        vector_kind = "closure_vector_permuted"
    left = setting_vector(source, setting_left, vector_kind)
    right = setting_vector(source, setting_right, vector_kind)
    denominator = norm(left) * norm(right)
    if denominator <= 1.0e-300:
        return 0.0
    operator = operator_for_control(source, control_id)
    return dot(left, matvec(operator, right)) / denominator


def all_pairs() -> list[tuple[str, str, str]]:
    pairs: list[tuple[str, str, str]] = []
    for pair_role, role_pairs in ALL_PAIR_ROLES.items():
        for left, right in role_pairs:
            pairs.append((pair_role, left, right))
    return pairs


def evaluate_geometry_rows(
    control_id: str,
    sources: list[SourceRow],
    config: GeometryConfig,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for source in sources:
        for pair_role, setting_left, setting_right in all_pairs():
            value = relational_invariant_g(source, setting_left, setting_right, control_id)
            vector_kind = "target"
            if control_id == "label_permuted_settings":
                vector_kind = "label_permuted"
            elif control_id == "closure_vector_permutation":
                vector_kind = "closure_vector_permuted"
            rows.append(
                {
                    "candidate_id": "b3_2_relational_geometry",
                    "control_id": control_id,
                    "pair_role": pair_role,
                    "source_id": source.source_id,
                    "hierarchy_label": source.hierarchy_label,
                    "seed": source.seed,
                    "n_side": source.n_side,
                    "setting_left": setting_left,
                    "setting_right": setting_right,
                    "G": f"{value:.12g}",
                    "operator_kind": control_id,
                    "vector_kind": vector_kind,
                    "closure_strength": f"{source_closure_strength(source):.12g}",
                    "chain_signature": source.fields["chain_signature"],
                }
            )
    return rows


def values_from_rows(rows: list[dict[str, Any]], control_id: str, pair_role: str | None = None) -> list[float]:
    output: list[float] = []
    for row in rows:
        if row["control_id"] != control_id:
            continue
        if pair_role is not None and row["pair_role"] != pair_role:
            continue
        output.append(float(row["G"]))
    return output


def geometry_stats(values: list[float], sign_tolerance: float) -> dict[str, Any]:
    if not values:
        return {
            "count": 0,
            "mean_G": 0.0,
            "min_G": 0.0,
            "max_G": 0.0,
            "range_G": 0.0,
            "variance_G": 0.0,
            "positive_count": 0,
            "negative_count": 0,
            "geometry_score": 0.0,
        }
    mean = sum(values) / float(len(values))
    variance = sum((value - mean) ** 2 for value in values) / float(len(values))
    positive = sum(1 for value in values if value > sign_tolerance)
    negative = sum(1 for value in values if value < -sign_tolerance)
    range_g = max(values) - min(values)
    geometry_score = range_g * min(positive, negative) / float(len(values))
    return {
        "count": len(values),
        "mean_G": mean,
        "min_G": min(values),
        "max_G": max(values),
        "range_G": range_g,
        "variance_G": variance,
        "positive_count": positive,
        "negative_count": negative,
        "geometry_score": geometry_score,
    }


def covariance_rows(sources: list[SourceRow], config: GeometryConfig) -> tuple[list[dict[str, Any]], float]:
    rows: list[dict[str, Any]] = []
    max_error = 0.0
    for source in sources:
        operator = chain_orientation_operator(source)
        transformed_operator = transform_matrix(operator, INDEX_PERMUTATION)
        for pair_role, setting_left, setting_right in all_pairs():
            left = setting_vector(source, setting_left)
            right = setting_vector(source, setting_right)
            denominator = norm(left) * norm(right)
            original = dot(left, matvec(operator, right)) / denominator if denominator else 0.0
            transformed_left = permute_vector(left, INDEX_PERMUTATION)
            transformed_right = permute_vector(right, INDEX_PERMUTATION)
            transformed_denominator = norm(transformed_left) * norm(transformed_right)
            transformed = (
                dot(transformed_left, matvec(transformed_operator, transformed_right)) / transformed_denominator
                if transformed_denominator
                else 0.0
            )
            error = abs(original - transformed)
            max_error = max(max_error, error)
            rows.append(
                {
                    "source_id": source.source_id,
                    "pair_role": pair_role,
                    "setting_left": setting_left,
                    "setting_right": setting_right,
                    "G_original": f"{original:.12g}",
                    "G_transformed": f"{transformed:.12g}",
                    "absolute_error": f"{error:.12g}",
                    "tolerance": config.covariance_tolerance,
                    "status": "COVARIANCE_PASS" if error <= config.covariance_tolerance else "COVARIANCE_FAIL",
                }
            )
    return rows, max_error


CONTROL_CONTRACTS = {
    "identity_pairing": {
        "purpose": "remove oriented chain transport while preserving vector dimensionality",
        "preserves": "setting vectors, source rows, token basis",
        "destroys": "skew orientation transport J_lambda",
        "expected_metric_response": "geometry_score lower than target by predeclared ratio",
        "invalidation_condition": "identity pairing preserves target geometry score",
    },
    "shuffled_pairing": {
        "purpose": "preserve skew matrix norm but break chain-aligned topology",
        "preserves": "operator sparsity and source closure strength",
        "destroys": "token adjacency alignment",
        "expected_metric_response": "geometry_score lower than target by predeclared ratio",
        "invalidation_condition": "shuffled pairing preserves target geometry score",
    },
    "random_orthogonal_pairing": {
        "purpose": "preserve orthogonal norm while randomizing orientation basis",
        "preserves": "operator norm and antisymmetry",
        "destroys": "semantic token orientation",
        "expected_metric_response": "geometry_score lower than target by predeclared ratio",
        "invalidation_condition": "random orthogonal pairing preserves target geometry score",
    },
    "topology_destroyed_pairing": {
        "purpose": "destroy chain topology entirely",
        "preserves": "source row availability and setting labels",
        "destroys": "all pairwise orientation edges",
        "expected_metric_response": "G collapses near zero",
        "invalidation_condition": "nonzero geometry survives zero-topology operator",
    },
    "label_permuted_settings": {
        "purpose": "test whether geometry depends on declared setting semantics",
        "preserves": "source rows and J_lambda",
        "destroys": "setting-label correspondence",
        "expected_metric_response": "geometry_score lower than target by predeclared ratio",
        "invalidation_condition": "label permutation preserves target geometry score",
    },
    "closure_vector_permutation": {
        "purpose": "break closure-vector token alignment",
        "preserves": "source rows, J_lambda, vector norms",
        "destroys": "u_a/u_b token alignment with J_lambda",
        "expected_metric_response": "geometry_score lower than target by predeclared ratio",
        "invalidation_condition": "permuted vectors preserve target geometry score",
    },
    "refinement_broken_control": {
        "purpose": "replace frozen target rows with a frozen matched control hierarchy",
        "preserves": "artifact schema, n_side/seed style, chain-signature format",
        "destroys": "target hierarchy identity",
        "expected_metric_response": "geometry_score lower than target by predeclared ratio",
        "invalidation_condition": "matched control hierarchy preserves target geometry score",
    },
}


def control_result_rows(geometry_rows: list[dict[str, Any]], config: GeometryConfig) -> tuple[list[dict[str, Any]], bool, dict[str, Any]]:
    target_values = values_from_rows(geometry_rows, "target_geometry")
    target_stats = geometry_stats(target_values, config.sign_tolerance)
    target_score = float(target_stats["geometry_score"])
    rows: list[dict[str, Any]] = []
    null_pass = True
    for control_id, contract in CONTROL_CONTRACTS.items():
        values = values_from_rows(geometry_rows, control_id)
        stats = geometry_stats(values, config.sign_tolerance)
        score = float(stats["geometry_score"])
        ratio = score / target_score if target_score > 1.0e-300 else math.inf
        if control_id == "topology_destroyed_pairing":
            status = "NULL_REJECTION_PASS" if abs(float(stats["max_G"])) <= config.sign_tolerance and abs(float(stats["min_G"])) <= config.sign_tolerance else "NULL_REJECTION_FAIL"
        else:
            status = "NULL_REJECTION_PASS" if ratio <= config.null_degradation_ratio else "NULL_REJECTION_FAIL"
        if status != "NULL_REJECTION_PASS":
            null_pass = False
        rows.append(
            {
                "control_id": control_id,
                **contract,
                "count": stats["count"],
                "mean_G": f"{stats['mean_G']:.12g}",
                "min_G": f"{stats['min_G']:.12g}",
                "max_G": f"{stats['max_G']:.12g}",
                "range_G": f"{stats['range_G']:.12g}",
                "variance_G": f"{stats['variance_G']:.12g}",
                "positive_count": stats["positive_count"],
                "negative_count": stats["negative_count"],
                "geometry_score": f"{score:.12g}",
                "degradation_ratio": f"{ratio:.12g}",
                "status": status,
            }
        )
    return rows, null_pass, target_stats


def gate_results(
    geometry_rows: list[dict[str, Any]],
    covariance_max_error: float,
    null_pass: bool,
    config: GeometryConfig,
) -> dict[str, Any]:
    target_all = geometry_stats(values_from_rows(geometry_rows, "target_geometry"), config.sign_tolerance)
    target_holdout = geometry_stats(values_from_rows(geometry_rows, "target_geometry", "holdout"), config.sign_tolerance)
    g1 = float(target_all["range_G"]) >= config.relational_range_tolerance
    g2 = int(target_all["positive_count"]) > 0 and int(target_all["negative_count"]) > 0
    g3 = covariance_max_error <= config.covariance_tolerance
    g4 = (
        float(target_holdout["range_G"]) >= config.holdout_range_tolerance
        and int(target_holdout["positive_count"]) > 0
        and int(target_holdout["negative_count"]) > 0
    )
    g5 = null_pass
    gates = {
        "G1_RELATIONAL_SENSITIVITY": {
            "passed": g1,
            "label": "RELATIONAL_STRUCTURE_DETECTED" if g1 else "RELATIONAL_STRUCTURE_NOT_DETECTED",
            "stats": target_all,
            "tolerance": config.relational_range_tolerance,
        },
        "G2_SIGN_STRUCTURE": {
            "passed": g2,
            "label": "SIGN_CHANGING_GEOMETRY_DETECTED" if g2 else "SIGN_CHANGING_GEOMETRY_NOT_DETECTED",
            "stats": target_all,
            "sign_tolerance": config.sign_tolerance,
        },
        "G3_COVARIANCE": {
            "passed": g3,
            "label": "COVARIANCE_PASS" if g3 else "COVARIANCE_FAIL",
            "max_error": covariance_max_error,
            "tolerance": config.covariance_tolerance,
        },
        "G4_HOLDOUT_TRANSFER": {
            "passed": g4,
            "label": "HOLDOUT_TRANSFER_PASS" if g4 else "HOLDOUT_TRANSFER_FAIL",
            "stats": target_holdout,
            "tolerance": config.holdout_range_tolerance,
        },
        "G5_NULL_REJECTION": {
            "passed": g5,
            "label": "NULL_REJECTION_PASS" if g5 else "NULL_REJECTION_FAIL",
            "degradation_ratio": config.null_degradation_ratio,
        },
    }
    authorized = all(gate["passed"] for gate in gates.values())
    gates["G6_CHSH_SCOREBOARD"] = {
        "passed": authorized,
        "label": "CHSH_SCORING_AUTHORIZED" if authorized else "CHSH_SCORING_NOT_AUTHORIZED",
        "computed": False,
        "reason": "This harness authorizes a later scoreboard run only; it does not compute CHSH.",
    }
    labels = [gate["label"] for gate in gates.values()]
    labels.extend(["HAOS_BELL_DERIVATION_NOT_ESTABLISHED", "MIXED / OPEN"])
    return {
        "gates": gates,
        "labels": labels,
        "chsh_scoring_authorized": authorized,
        "chsh_scoring_computed": False,
    }


def precommitment_contract(config: GeometryConfig) -> dict[str, Any]:
    return {
        "name": "B3.2_PRECOMMITMENT_CONTRACT",
        "purpose": "Test whether frozen HAOS-IIP structure yields a relational, sign-changing invariant G(a,b,lambda) before CHSH is computed.",
        "status": "GEOMETRY_AUDIT_PRECOMMITMENT",
        "source_artifacts": {
            "phase18_shell_ordering_metrics": config.source_artifact,
            "source_hash": hash_file(REPO_ROOT / config.source_artifact),
            "target_hierarchy_label": config.target_hierarchy_label,
            "refinement_control_label": config.refinement_control_label,
            "source_limit": config.source_limit,
        },
        "setting_representation": {
            "tokens": list(TOKENS),
            "setting_definitions": {key: list(value) for key, value in SETTING_DEFINITIONS.items()},
            "held_back_chsh_pairs": [list(pair) for pair in CHSH_HELD_BACK_PAIRS],
            "holdout_pairs": [list(pair) for pair in HOLDOUT_PAIRS],
        },
        "closure_vectors": {
            "definition": "u_s(lambda) is a signed two-token contrast vector weighted by token rank in the frozen chain_signature.",
            "normalization": "G divides by ||u_a|| ||u_b||.",
            "sign_convention": "plus token contributes positively; minus token contributes negatively.",
        },
        "pairing_operator": {
            "name": "chain_orientation_operator",
            "definition": "J_lambda is a skew-symmetric nearest-neighbor transport matrix along the frozen chain_signature.",
            "edge_weight": "source_closure_strength(lambda)/(edge_index+1)",
            "semantic_justification": "orientation is taken from observed chain order, not from Bell target signs.",
        },
        "relational_invariant": "G(a,b,lambda)=<u_a,J_lambda u_b>/(||u_a||||u_b||)",
        "covariance_transformations": {
            "joint_token_permutation": list(INDEX_PERMUTATION),
            "expected": "G is invariant when u_a, u_b, and J_lambda are transformed together.",
        },
        "null_operator_contracts": CONTROL_CONTRACTS,
        "seeds": list(config.seeds),
        "tolerances": {
            "relational_range_tolerance": config.relational_range_tolerance,
            "sign_tolerance": config.sign_tolerance,
            "covariance_tolerance": config.covariance_tolerance,
            "holdout_range_tolerance": config.holdout_range_tolerance,
            "null_degradation_ratio": config.null_degradation_ratio,
        },
        "gate_logic": [
            "G1_RELATIONAL_SENSITIVITY must pass before promotion.",
            "G2_SIGN_STRUCTURE must pass before promotion.",
            "G3_COVARIANCE must pass before promotion.",
            "G4_HOLDOUT_TRANSFER must pass before promotion.",
            "G5_NULL_REJECTION must pass before promotion.",
            "G6_CHSH_SCOREBOARD is authorization only; CHSH is not computed here.",
        ],
        "hard_prohibitions": [
            "do not tune beta against CHSH output",
            "do not tune J_lambda against CHSH output",
            "do not tune setting mappings against CHSH output",
            "do not tune u_a/u_b extraction against CHSH output",
            "do not tune sign convention against CHSH output",
            "do not tune normalization against CHSH output",
            "do not use proxy objectives correlated with a frozen quantum table",
        ],
        "non_claims": [
            "not a physical Bell experiment",
            "not a loophole-free Bell test",
            "not a quantum mechanics derivation",
            "not evidence that HAOS-IIP derives Bell correlations",
        ],
    }


def write_report(path: Path, result: dict[str, Any], control_rows: list[dict[str, Any]]) -> None:
    gates = result["gate_results"]["gates"]
    lines = [
        "# B3.2 Relational Geometry Audit",
        "",
        "Implemented fact: this harness computes a frozen HAOS-IIP relational invariant before any CHSH score.",
        "Design choice: CHSH scoring is authorization-gated and not computed here.",
        "Heuristic: `J_lambda` is an oriented chain-transport pairing from Phase 18 chain signatures.",
        "Unverified hypothesis: HAOS-IIP Bell derivation remains unestablished.",
        "",
        "## Gate Results",
    ]
    for gate_name, gate in gates.items():
        lines.append(f"- {gate_name}: `{gate['label']}`")
    lines.extend(
        [
            "",
            "## Target Geometry Summary",
            f"- all-pair range: `{gates['G1_RELATIONAL_SENSITIVITY']['stats']['range_G']:.12g}`",
            f"- all-pair positives/negatives: `{gates['G2_SIGN_STRUCTURE']['stats']['positive_count']}` / `{gates['G2_SIGN_STRUCTURE']['stats']['negative_count']}`",
            f"- holdout range: `{gates['G4_HOLDOUT_TRANSFER']['stats']['range_G']:.12g}`",
            f"- covariance max error: `{gates['G3_COVARIANCE']['max_error']:.12g}`",
            "",
            "## Null Controls",
        ]
    )
    for row in control_rows:
        lines.append(
            "- {control}: score `{score}`, ratio `{ratio}`, status `{status}`".format(
                control=row["control_id"],
                score=row["geometry_score"],
                ratio=row["degradation_ratio"],
                status=row["status"],
            )
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "- This is a geometry audit harness, not a Bell-score optimizer.",
            "- No CHSH S value is computed here.",
            "- A failed earlier gate blocks later promotion.",
            "- HAOS Bell derivation is not established.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_audit(config: GeometryConfig, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    contract = precommitment_contract(config)
    leakage_detected, leakage_hits = detect_geometry_leakage(generator_sources())
    contract["target_leakage_guard"] = {"detected": leakage_detected, "hits": leakage_hits}
    contract["contract_hash"] = stable_json_hash("b3_2_contract", contract)
    write_json(output_dir / CONTRACT_PATH.name, contract)

    target_sources = load_sources(config, config.target_hierarchy_label)
    refinement_sources = load_sources(config, config.refinement_control_label)
    geometry_rows: list[dict[str, Any]] = []
    for control_id in (
        "target_geometry",
        "identity_pairing",
        "shuffled_pairing",
        "random_orthogonal_pairing",
        "topology_destroyed_pairing",
        "label_permuted_settings",
        "closure_vector_permutation",
    ):
        geometry_rows.extend(evaluate_geometry_rows(control_id, target_sources, config))
    geometry_rows.extend(evaluate_geometry_rows("refinement_broken_control", refinement_sources, config))

    covariance_diagnostics, covariance_max_error = covariance_rows(target_sources, config)
    control_rows, null_pass, target_stats = control_result_rows(geometry_rows, config)
    gates = gate_results(geometry_rows, covariance_max_error, null_pass and not leakage_detected, config)
    if leakage_detected:
        gates["labels"].append("TARGET_LEAKAGE_DETECTED")
        gates["gates"]["G5_NULL_REJECTION"]["passed"] = False
        gates["gates"]["G5_NULL_REJECTION"]["label"] = "NULL_REJECTION_FAIL"
        gates["gates"]["G6_CHSH_SCOREBOARD"]["passed"] = False
        gates["gates"]["G6_CHSH_SCOREBOARD"]["label"] = "CHSH_SCORING_NOT_AUTHORIZED"
        gates["chsh_scoring_authorized"] = False

    result: dict[str, Any] = {
        "contract_hash": contract["contract_hash"],
        "config": asdict(config),
        "source_hash": hash_file(REPO_ROOT / config.source_artifact),
        "gate_results": gates,
        "target_geometry_stats": target_stats,
        "target_leakage_guard": {"detected": leakage_detected, "hits": leakage_hits},
        "chsh_score": None,
        "chsh_scoring_computed": False,
        "outputs": {
            "precommitment_contract": repo_rel(output_dir / CONTRACT_PATH.name),
            "geometry_matrix": repo_rel(output_dir / GEOMETRY_MATRIX_PATH.name),
            "covariance_diagnostics": repo_rel(output_dir / COVARIANCE_PATH.name),
            "control_results": repo_rel(output_dir / CONTROL_RESULTS_PATH.name),
            "geometry_report": repo_rel(output_dir / REPORT_PATH.name),
            "result": repo_rel(output_dir / RESULT_PATH.name),
        },
    }
    hash_payload = {key: value for key, value in result.items() if key != "outputs"}
    result["result_hash"] = stable_json_hash("b3_2_result", hash_payload)

    write_csv(output_dir / GEOMETRY_MATRIX_PATH.name, geometry_rows, GEOMETRY_FIELDNAMES)
    write_csv(output_dir / COVARIANCE_PATH.name, covariance_diagnostics, COVARIANCE_FIELDNAMES)
    write_csv(output_dir / CONTROL_RESULTS_PATH.name, control_rows, CONTROL_FIELDNAMES)
    write_json(output_dir / RESULT_PATH.name, result)
    write_report(output_dir / REPORT_PATH.name, result, control_rows)
    return result


def main() -> None:
    args = parse_args()
    result = run_audit(GeometryConfig(), args.output_dir)
    print(
        json.dumps(
            {
                "labels": result["gate_results"]["labels"],
                "chsh_scoring_authorized": result["gate_results"]["chsh_scoring_authorized"],
                "chsh_scoring_computed": result["chsh_scoring_computed"],
                "result_hash": result["result_hash"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
