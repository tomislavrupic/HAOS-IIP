"""Biology Line E external dataset bridge.

This module adapts real expression time-series data into lightweight
HAOS-style recoverability diagnostics. It is intentionally conservative:
without real input data, the validation status is OPEN rather than PASS.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import gzip
import json
import math
import re
from pathlib import Path
from typing import Iterable

import numpy as np


COLLAPSE_THRESHOLD = 0.70
SUSTAIN_STEPS = 1
VISIBLE_FAILURE_MATCH_THRESHOLD = 0.50
TRAJECTORY_IDENTITY_FAILURE_THRESHOLD = 0.15
CONTROL_SEED = 20260424
COLLAPSE_THRESHOLD_GRID = (0.55, 0.60, 0.65, 0.70, 0.75, 0.80)
VISIBLE_MATCH_THRESHOLD_GRID = (0.30, 0.40, 0.50, 0.60)
IDENTIFIER_HEADERS = {
    "id",
    "id_ref",
    "identifier",
    "gene",
    "gene_id",
    "gene_symbol",
    "probe",
    "probe_id",
    "feature",
    "feature_id",
}
CONTROL_FEATURE_NAMES = {
    "EMPTY",
    "3XSSC",
    "BUFFER",
    "BLANK",
    "CONTROL",
}


@dataclass(frozen=True)
class ExpressionDataset:
    source_name: str
    feature_ids: list[str]
    sample_labels: list[str]
    time_values: list[float]
    matrix: np.ndarray


@dataclass(frozen=True)
class BridgeRun:
    label: str
    dataset: ExpressionDataset
    rows: list[dict[str, float | int | bool | str | None]]
    summary: dict[str, object]
    standard_metrics: list[dict[str, float | int | str]]


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _parse_float(value: str) -> float | None:
    text = value.strip()
    if text in {"", "NA", "NaN", "null", "NULL"}:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    if not math.isfinite(number):
        return None
    return number


def _infer_time(label: str, fallback: int) -> float:
    minute_match = re.search(r"(\d+(?:\.\d+)?)\s*(?:min|minute)", label, re.IGNORECASE)
    if minute_match:
        return float(minute_match.group(1))
    number_match = re.search(r"(\d+(?:\.\d+)?)", label)
    if number_match:
        return float(number_match.group(1))
    return float(fallback)


def _numeric_column_indices(rows: list[list[str]], header: list[str], min_fraction: float = 0.60) -> list[int]:
    numeric_indices: list[int] = []
    if not rows:
        return numeric_indices

    for idx, _name in enumerate(header):
        numeric_count = 0
        checked = 0
        for row in rows[: min(len(rows), 500)]:
            if idx >= len(row):
                continue
            checked += 1
            if _parse_float(row[idx]) is not None:
                numeric_count += 1
        if checked and numeric_count / checked >= min_fraction:
            numeric_indices.append(idx)
    return numeric_indices


def _normalized_header(value: str) -> str:
    return re.sub(r"[^a-z0-9_]+", "_", value.strip().lower()).strip("_")


def _sample_numeric_column_indices(rows: list[list[str]], header: list[str]) -> list[int]:
    identifier_indices = {
        idx
        for idx, name in enumerate(header)
        if _normalized_header(name) in IDENTIFIER_HEADERS
    }
    # The bridge expects one or more leading identifier columns. Excluding the
    # first column prevents numeric probe IDs from being treated as samples.
    identifier_indices.add(0)
    return [
        idx
        for idx in _numeric_column_indices(rows, header)
        if idx not in identifier_indices
    ]


def _feature_id_index(header: list[str]) -> int:
    preferred = ("identifier", "gene", "gene_symbol", "gene_id", "id_ref")
    normalized = [_normalized_header(name) for name in header]
    for name in preferred:
        if name in normalized:
            return normalized.index(name)
    return 0


def _is_valid_feature_id(feature_id: str) -> bool:
    text = feature_id.strip().upper()
    if not text:
        return False
    if text in CONTROL_FEATURE_NAMES:
        return False
    if text.startswith("GENOMIC"):
        return False
    if " " in text:
        return False
    return True


def load_csv_expression(path: Path) -> ExpressionDataset:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = [row for row in reader if row]

    numeric_indices = _sample_numeric_column_indices(rows, header)
    if len(numeric_indices) < 2:
        raise ValueError("CSV input must contain at least two numeric sample columns")

    id_index = _feature_id_index(header)
    feature_ids: list[str] = []
    values: list[list[float]] = []
    for row in rows:
        parsed_row: list[float] = []
        for idx in numeric_indices:
            value = _parse_float(row[idx] if idx < len(row) else "")
            parsed_row.append(np.nan if value is None else value)
        feature_id = row[id_index].strip() if id_index < len(row) and row[id_index].strip() else f"feature_{len(feature_ids)}"
        if np.isfinite(parsed_row).sum() >= 2 and _is_valid_feature_id(feature_id):
            feature_ids.append(feature_id)
            values.append(parsed_row)

    matrix, feature_ids = _clean_matrix_and_ids(np.array(values, dtype=float), feature_ids)
    sample_labels = [header[idx] for idx in numeric_indices]
    time_values = [_infer_time(label, i) for i, label in enumerate(sample_labels)]
    return ExpressionDataset(path.name, feature_ids, sample_labels, time_values, matrix)


def load_geo_soft_expression(path: Path) -> ExpressionDataset:
    sample_titles: list[str] = []
    sample_label_by_id: dict[str, str] = {}
    current_subset_description: str | None = None
    dataset_value_type = ""
    table_lines: list[str] = []
    in_table = False

    with _open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("!dataset_sample_title"):
                _key, value = line.split("=", 1)
                sample_titles.append(value.strip())
            elif line.startswith("!dataset_value_type"):
                _key, value = line.split("=", 1)
                dataset_value_type = value.strip().lower()
            elif line.startswith("!subset_description"):
                _key, value = line.split("=", 1)
                current_subset_description = value.strip()
            elif line.startswith("!subset_sample_id"):
                _key, value = line.split("=", 1)
                sample_id = value.strip()
                if current_subset_description:
                    sample_label_by_id[sample_id] = current_subset_description
            elif line.startswith("#") and "=" in line:
                key, value = line[1:].split("=", 1)
                sample_id = key.strip()
                if sample_id.startswith("GSM") and sample_id not in sample_label_by_id:
                    sample_label_by_id[sample_id] = value.strip()
            elif line == "!dataset_table_begin":
                in_table = True
            elif line == "!dataset_table_end":
                in_table = False
            elif in_table:
                table_lines.append(line)

    if not table_lines:
        raise ValueError("No GEO dataset table found in SOFT file")

    reader = csv.reader(table_lines, delimiter="\t")
    header = next(reader)
    rows = [row for row in reader if row]
    numeric_indices = _sample_numeric_column_indices(rows, header)
    if len(numeric_indices) < 2:
        raise ValueError("SOFT file must contain at least two numeric sample columns")

    id_index = _feature_id_index(header)
    feature_ids: list[str] = []
    values: list[list[float]] = []
    for row in rows:
        parsed_row: list[float] = []
        for idx in numeric_indices:
            value = _parse_float(row[idx] if idx < len(row) else "")
            parsed_row.append(np.nan if value is None else value)
        feature_id = row[id_index].strip() if id_index < len(row) and row[id_index].strip() else f"feature_{len(feature_ids)}"
        if np.isfinite(parsed_row).sum() >= 2 and _is_valid_feature_id(feature_id):
            feature_ids.append(feature_id)
            values.append(parsed_row)

    matrix, feature_ids = _clean_matrix_and_ids(np.array(values, dtype=float), feature_ids)
    header_labels = [header[idx] for idx in numeric_indices]
    if len(sample_titles) >= len(header_labels):
        sample_labels = sample_titles[: len(header_labels)]
    else:
        sample_labels = [sample_label_by_id.get(label, label) for label in header_labels]
    time_values = [_infer_time(label, i) for i, label in enumerate(sample_labels)]
    if "log ratio" in dataset_value_type and time_values and min(time_values) > 0.0:
        matrix = np.column_stack([np.zeros(matrix.shape[0]), matrix])
        sample_labels = ["0 minute inferred log-ratio baseline", *sample_labels]
        time_values = [0.0, *time_values]
    return ExpressionDataset(path.name, feature_ids, sample_labels, time_values, matrix)


def _clean_matrix_and_ids(matrix: np.ndarray, feature_ids: list[str]) -> tuple[np.ndarray, list[str]]:
    if matrix.ndim != 2 or matrix.shape[1] < 2:
        raise ValueError("expression matrix must be feature x sample with at least two samples")

    cleaned = matrix.astype(float, copy=True)
    col_means = np.nanmean(cleaned, axis=0)
    col_means = np.where(np.isfinite(col_means), col_means, 0.0)
    missing_rows, missing_cols = np.where(~np.isfinite(cleaned))
    cleaned[missing_rows, missing_cols] = col_means[missing_cols]
    keep = np.std(cleaned, axis=1) > 1e-12
    cleaned = cleaned[keep]
    cleaned_ids = [feature_id for feature_id, keep_row in zip(feature_ids, keep) if keep_row]
    if cleaned.shape[0] == 0:
        raise ValueError("expression matrix contains no varying features")
    return cleaned, cleaned_ids


def select_top_variable_features(dataset: ExpressionDataset, max_features: int = 2000) -> ExpressionDataset:
    if dataset.matrix.shape[0] <= max_features:
        return dataset
    variances = np.var(dataset.matrix, axis=1)
    keep = np.argsort(variances)[-max_features:]
    keep.sort()
    return ExpressionDataset(
        source_name=dataset.source_name,
        feature_ids=[dataset.feature_ids[i] for i in keep if i < len(dataset.feature_ids)],
        sample_labels=dataset.sample_labels,
        time_values=dataset.time_values,
        matrix=dataset.matrix[keep],
    )


def make_time_shuffle_control(dataset: ExpressionDataset) -> ExpressionDataset:
    """Deterministically scramble non-baseline sample order."""
    if dataset.matrix.shape[1] <= 2:
        return dataset
    rng = np.random.default_rng(CONTROL_SEED)
    tail = np.arange(1, dataset.matrix.shape[1])
    rng.shuffle(tail)
    order = np.concatenate(([0], tail))
    return ExpressionDataset(
        source_name=f"{dataset.source_name}__time_shuffle_control",
        feature_ids=dataset.feature_ids,
        sample_labels=[dataset.sample_labels[i] for i in order],
        time_values=[dataset.time_values[i] for i in order],
        matrix=dataset.matrix[:, order],
    )


def make_feature_shuffle_control(dataset: ExpressionDataset) -> ExpressionDataset:
    """Permute gene identities while preserving each sample distribution."""
    rng = np.random.default_rng(CONTROL_SEED + 1)
    shuffled = dataset.matrix.copy()
    for col in range(shuffled.shape[1]):
        shuffled[:, col] = shuffled[rng.permutation(shuffled.shape[0]), col]
    return ExpressionDataset(
        source_name=f"{dataset.source_name}__feature_shuffle_control",
        feature_ids=dataset.feature_ids,
        sample_labels=dataset.sample_labels,
        time_values=dataset.time_values,
        matrix=shuffled,
    )


def make_row_permutation_control(dataset: ExpressionDataset) -> ExpressionDataset:
    """Destroy each gene trajectory while preserving row-level value ranges."""
    rng = np.random.default_rng(CONTROL_SEED + 2)
    shuffled = dataset.matrix.copy()
    for row in range(shuffled.shape[0]):
        shuffled[row] = shuffled[row, rng.permutation(shuffled.shape[1])]
    return ExpressionDataset(
        source_name=f"{dataset.source_name}__row_permutation_control",
        feature_ids=dataset.feature_ids,
        sample_labels=dataset.sample_labels,
        time_values=dataset.time_values,
        matrix=shuffled,
    )


def compute_standard_expression_metrics(dataset: ExpressionDataset) -> list[dict[str, float | int | str]]:
    """Compute plain biology-style expression movement summaries.

    These are not HAOS metrics. They provide a simple external comparator:
    how many features cross absolute log-ratio thresholds at each sample?
    """
    matrix = dataset.matrix
    baseline = matrix[:, 0]
    delta = matrix - baseline[:, None]
    rows: list[dict[str, float | int | str]] = []
    for idx in range(matrix.shape[1]):
        abs_delta = np.abs(delta[:, idx])
        rows.append(
            {
                "sample_index": idx,
                "sample_label": dataset.sample_labels[idx],
                "time_value": float(dataset.time_values[idx]),
                "mean_abs_change": float(np.mean(abs_delta)),
                "median_abs_change": float(np.median(abs_delta)),
                "fraction_abs_change_ge_0_5": float(np.mean(abs_delta >= 0.5)),
                "fraction_abs_change_ge_1_0": float(np.mean(abs_delta >= 1.0)),
            }
        )
    return rows


def compute_trajectory_identity_retention(matrix: np.ndarray) -> np.ndarray:
    """Retain gene identity through the response trajectory.

    This is a standard time-series comparator, not a HAOS metric. It checks
    whether the same feature identities carry related response vectors between
    adjacent samples. Feature-shuffled controls should fail this comparator.
    """
    baseline = matrix[:, [0]]
    delta = matrix - baseline
    retention = np.ones(matrix.shape[1], dtype=float)
    for idx in range(1, matrix.shape[1]):
        neighbor = idx + 1 if idx + 1 < matrix.shape[1] else idx - 1
        left = delta[:, idx]
        right = delta[:, neighbor]
        if float(np.std(left)) <= 1e-12 or float(np.std(right)) <= 1e-12:
            retention[idx] = 1.0
        else:
            retention[idx] = max(0.0, float(np.corrcoef(left, right)[0, 1]))
    return retention


def load_gene_set(path: Path) -> set[str]:
    genes: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            genes.add(_normalize_gene_id(text.split()[0]))
    return genes


def _normalize_gene_id(value: str) -> str:
    return re.sub(r"[^A-Z0-9_-]", "", value.upper())


def _average_precision(labels: np.ndarray, scores: np.ndarray) -> float:
    positives = int(np.sum(labels))
    if positives == 0:
        return 0.0
    order = np.argsort(-scores)
    hits = 0
    precision_sum = 0.0
    for rank, idx in enumerate(order, start=1):
        if labels[idx]:
            hits += 1
            precision_sum += hits / rank
    return precision_sum / positives


def _top_k_enrichment(labels: np.ndarray, scores: np.ndarray, k: int) -> float:
    if labels.size == 0:
        return 0.0
    k = min(k, labels.size)
    if k <= 0:
        return 0.0
    baseline = max(float(np.mean(labels)), 1e-12)
    order = np.argsort(-scores)[:k]
    return float(np.mean(labels[order]) / baseline)


def compute_feature_rankings(
    dataset: ExpressionDataset,
    rows: list[dict[str, object]],
) -> list[dict[str, float | str]]:
    """Rank features by HAOS-weighted trajectory contribution and fold change."""
    matrix = dataset.matrix
    baseline = matrix[:, 0]
    abs_delta = np.abs(matrix - baseline[:, None])
    max_abs_change = np.max(abs_delta, axis=1)
    mean_abs_change = np.mean(abs_delta, axis=1)
    trajectory_total_variation = np.sum(np.abs(np.diff(matrix, axis=1)), axis=1)

    recoverability = np.array([float(row["recoverability"]) for row in rows])
    delta_persistence = np.array([float(row["delta_persistence"]) for row in rows])
    collapse_pressure = np.maximum(0.0, COLLAPSE_THRESHOLD - recoverability)
    drop_pressure = np.maximum(0.0, -delta_persistence)
    weights = collapse_pressure + drop_pressure
    if float(np.sum(weights)) <= 1e-12:
        weights = np.ones_like(recoverability)
    weights = weights / np.sum(weights)
    haos_loss_score = abs_delta @ weights

    ranked_rows: list[dict[str, float | str]] = []
    for idx, feature_id in enumerate(dataset.feature_ids):
        ranked_rows.append(
            {
                "feature_id": feature_id,
                "haos_loss_score": float(haos_loss_score[idx]),
                "max_abs_change": float(max_abs_change[idx]),
                "mean_abs_change": float(mean_abs_change[idx]),
                "trajectory_total_variation": float(trajectory_total_variation[idx]),
            }
        )
    return ranked_rows


def evaluate_gene_set(
    rankings: list[dict[str, float | str]],
    gene_set: set[str],
    min_average_precision_margin: float = 0.02,
) -> dict[str, object]:
    normalized_ids = [_normalize_gene_id(str(row["feature_id"])) for row in rankings]
    labels = np.array([gene_id in gene_set for gene_id in normalized_ids], dtype=bool)
    matched = int(np.sum(labels))
    total = len(labels)

    if total == 0 or matched == 0:
        return {
            "gene_set_status": "OPEN_NO_MATCHED_GENES",
            "feature_count": total,
            "known_gene_set_size": len(gene_set),
            "matched_gene_count": matched,
            "haos_average_precision": None,
            "fold_change_average_precision": None,
            "average_precision_margin": None,
            "haos_beats_fold_change": False,
            "failure_rule": "Gene-set comparison is open when no supplied known genes match feature identifiers.",
        }

    haos_scores = np.array([float(row["haos_loss_score"]) for row in rankings])
    fold_scores = np.array([float(row["max_abs_change"]) for row in rankings])
    haos_ap = _average_precision(labels, haos_scores)
    fold_ap = _average_precision(labels, fold_scores)
    margin = haos_ap - fold_ap
    haos_beats = bool(margin >= min_average_precision_margin)

    top_values = {}
    for k in (25, 50, 100, 250):
        top_values[f"haos_top_{k}_enrichment"] = _top_k_enrichment(labels, haos_scores, k)
        top_values[f"fold_change_top_{k}_enrichment"] = _top_k_enrichment(labels, fold_scores, k)

    return {
        "gene_set_status": "GENE_SET_EVALUATED",
        "feature_count": total,
        "known_gene_set_size": len(gene_set),
        "matched_gene_count": matched,
        "matched_gene_fraction": float(matched / total),
        "haos_average_precision": float(haos_ap),
        "fold_change_average_precision": float(fold_ap),
        "average_precision_margin": float(margin),
        "haos_beats_fold_change": haos_beats,
        "minimum_average_precision_margin": min_average_precision_margin,
        "failure_rule": "If HAOS average precision does not exceed simple fold-change average precision by the predeclared margin, the biology-specific ranking claim fails.",
        **top_values,
    }


def compute_recoverability_metrics(
    dataset: ExpressionDataset,
    collapse_threshold: float = COLLAPSE_THRESHOLD,
    sustain_steps: int = SUSTAIN_STEPS,
) -> list[dict[str, float | int | bool | str | None]]:
    matrix = dataset.matrix
    identity_retention = compute_trajectory_identity_retention(matrix)
    baseline = matrix[:, 0]
    distances = np.sqrt(np.mean((matrix - baseline[:, None]) ** 2, axis=0))
    max_distance = max(float(np.max(distances)), 1e-9)
    expression_match = np.clip(1.0 - distances / max_distance, 0.0, 1.0)

    step_distances = np.zeros(matrix.shape[1], dtype=float)
    if matrix.shape[1] > 1:
        step_distances[1:] = np.sqrt(np.mean(np.diff(matrix, axis=1) ** 2, axis=0))
    max_step = max(float(np.max(step_distances)), 1e-9)
    step_coherence = np.clip(1.0 - step_distances / max_step, 0.0, 1.0)

    support_retention = np.ones(matrix.shape[1], dtype=float)
    recoverability = np.clip(
        0.72 * expression_match + 0.20 * step_coherence + 0.08 * support_retention,
        0.0,
        1.0,
    )

    delta = np.zeros_like(recoverability)
    delta[1:] = recoverability[1:] - recoverability[:-1]
    k_star_index = find_k_star(recoverability, collapse_threshold, sustain_steps)
    k_star_time = dataset.time_values[k_star_index] if k_star_index is not None else None

    rows: list[dict[str, float | int | bool | str | None]] = []
    min_delta_so_far = 0.0
    for idx, value in enumerate(recoverability):
        min_delta_so_far = min(min_delta_so_far, float(delta[idx]))
        trajectory_identity_failure = bool(
            idx > 0 and identity_retention[idx] < TRAJECTORY_IDENTITY_FAILURE_THRESHOLD
        )
        visible_failure = bool(
            expression_match[idx] < VISIBLE_FAILURE_MATCH_THRESHOLD
            or trajectory_identity_failure
        )
        safety_margin = compute_safety_margin(dataset.time_values[idx], k_star_time, dataset.time_values[-1])
        rows.append(
            {
                "sample_index": idx,
                "sample_label": dataset.sample_labels[idx],
                "time_value": dataset.time_values[idx],
                "recoverability": float(value),
                "delta_persistence": float(delta[idx]),
                "k_star": k_star_index,
                "k_star_time": k_star_time,
                "safety_margin": safety_margin,
                "confidence": abs(min_delta_so_far),
                "expression_match": float(expression_match[idx]),
                "step_coherence": float(step_coherence[idx]),
                "trajectory_identity_retention": float(identity_retention[idx]),
                "trajectory_identity_failure": trajectory_identity_failure,
                "support_retention": float(support_retention[idx]),
                "visible_failure": visible_failure,
            }
        )
    return rows


def find_k_star(values: Iterable[float], threshold: float, sustain_steps: int) -> int | None:
    series = list(values)
    for idx, value in enumerate(series):
        end = idx + sustain_steps + 1
        if end > len(series):
            continue
        if value < threshold and all(next_value < threshold for next_value in series[idx:end]):
            return idx
    return None


def compute_safety_margin(current_time: float, k_star_time: float | None, max_time: float) -> float:
    if k_star_time is None:
        return max(0.0, max_time - current_time)
    return k_star_time - current_time


def first_visible_failure(rows: list[dict[str, object]]) -> float | None:
    for row in rows:
        if bool(row["visible_failure"]):
            return float(row["time_value"])
    return None


def first_visible_failure_index(rows: list[dict[str, object]]) -> int | None:
    for idx, row in enumerate(rows):
        if bool(row["visible_failure"]):
            return idx
    return None


def _first_true_index(values: Iterable[bool]) -> int | None:
    for idx, value in enumerate(values):
        if value:
            return idx
    return None


def _safe_correlation(left: np.ndarray, right: np.ndarray) -> float | None:
    if left.size != right.size or left.size < 2:
        return None
    if float(np.std(left)) <= 1e-12 or float(np.std(right)) <= 1e-12:
        return None
    return float(np.corrcoef(left, right)[0, 1])


def analyze_failure_modes(
    rows: list[dict[str, object]],
    standard_metrics: list[dict[str, object]],
) -> dict[str, object]:
    recoverability = np.array([float(row["recoverability"]) for row in rows])
    expression_match = np.array([float(row["expression_match"]) for row in rows])
    step_coherence = np.array([float(row["step_coherence"]) for row in rows])
    identity_retention = np.array([float(row["trajectory_identity_retention"]) for row in rows])
    visible_flags = [bool(row["visible_failure"]) for row in rows]
    k_star_index = rows[0]["k_star"] if rows else None
    first_visible_index = _first_true_index(visible_flags)

    if (
        k_star_index is not None
        and first_visible_index is not None
        and int(k_star_index) < int(first_visible_index)
    ):
        failure_mode = "OBSERVED_EARLY_DETECTION_WINDOW"
        interpretation = "Observed data have k_star before the visible-failure proxy; controls decide whether this is a probe pass."
    elif first_visible_index == 1:
        failure_mode = "NO_PRE_VISIBLE_POST_BASELINE_SAMPLE"
        interpretation = (
            "The visible-failure proxy fires at the first post-baseline sample, "
            "so this dataset has no measured pre-visible window for early detection."
        )
    elif k_star_index == first_visible_index:
        failure_mode = "COLLAPSE_AND_VISIBLE_FAILURE_COINCIDE"
        interpretation = "k_star and visible failure occur at the same sampled time point."
    elif k_star_index is None:
        failure_mode = "NO_K_STAR"
        interpretation = "recoverability never crossed the sustained collapse threshold."
    else:
        failure_mode = "EARLY_WINDOW_AVAILABLE_BUT_NOT_RESOLVED"
        interpretation = "a pre-visible window exists, but the current proxy did not produce an earlier k_star."

    recoverability_expression_corr = _safe_correlation(recoverability, expression_match)
    recoverability_step_corr = _safe_correlation(recoverability, step_coherence)
    recoverability_identity_corr = _safe_correlation(recoverability, identity_retention)
    first_identity_failure_index = _first_true_index(
        bool(row["trajectory_identity_failure"]) for row in rows
    )

    first_standard_large_change_index = None
    for row in standard_metrics:
        if float(row["fraction_abs_change_ge_1_0"]) >= 0.50:
            first_standard_large_change_index = int(row["sample_index"])
            break

    return {
        "failure_mode": failure_mode,
        "interpretation": interpretation,
        "k_star_index": k_star_index,
        "k_star_time": rows[0]["k_star_time"] if rows else None,
        "first_visible_index": first_visible_index,
        "first_visible_time": rows[first_visible_index]["time_value"] if first_visible_index is not None else None,
        "first_standard_large_change_index": first_standard_large_change_index,
        "first_standard_large_change_time": (
            standard_metrics[first_standard_large_change_index]["time_value"]
            if first_standard_large_change_index is not None
            else None
        ),
        "recoverability_expression_match_correlation": recoverability_expression_corr,
        "recoverability_step_coherence_correlation": recoverability_step_corr,
        "recoverability_identity_retention_correlation": recoverability_identity_corr,
        "first_identity_failure_index": first_identity_failure_index,
        "first_identity_failure_time": (
            rows[first_identity_failure_index]["time_value"]
            if first_identity_failure_index is not None
            else None
        ),
        "dominant_failure_note": "Current recoverability is strongly tied to expression distance from baseline when correlation magnitude is high.",
    }


def threshold_sensitivity_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    recoverability = [float(row["recoverability"]) for row in rows]
    expression_match = [float(row["expression_match"]) for row in rows]
    output_rows: list[dict[str, object]] = []
    for collapse_threshold in COLLAPSE_THRESHOLD_GRID:
        k_star_index = find_k_star(recoverability, collapse_threshold, SUSTAIN_STEPS)
        k_star_time = rows[k_star_index]["time_value"] if k_star_index is not None else None
        for visible_threshold in VISIBLE_MATCH_THRESHOLD_GRID:
            visible_index = _first_true_index(value < visible_threshold for value in expression_match)
            visible_time = rows[visible_index]["time_value"] if visible_index is not None else None
            early_detection = bool(
                k_star_index is not None
                and visible_index is not None
                and k_star_index < visible_index
            )
            output_rows.append(
                {
                    "collapse_threshold": collapse_threshold,
                    "visible_match_threshold": visible_threshold,
                    "k_star_index": k_star_index,
                    "k_star_time": k_star_time,
                    "first_visible_index": visible_index,
                    "first_visible_time": visible_time,
                    "early_detection": early_detection,
                }
            )
    return output_rows


def summarize_bridge(
    dataset: ExpressionDataset,
    rows: list[dict[str, object]],
    standard_metrics: list[dict[str, object]] | None = None,
    label: str = "observed",
) -> dict[str, object]:
    baseline_stable = bool(rows and float(rows[0]["recoverability"]) >= 0.95)
    k_star_index = rows[0]["k_star"] if rows else None
    k_star_time = rows[0]["k_star_time"] if rows else None
    visible_index = first_visible_failure_index(rows)
    visible_time = first_visible_failure(rows)
    early_detection = bool(
        k_star_index is not None
        and visible_index is not None
        and int(k_star_index) < int(visible_index)
    )
    max_fraction_abs_change = 0.0
    if standard_metrics:
        max_fraction_abs_change = max(float(row["fraction_abs_change_ge_1_0"]) for row in standard_metrics)

    return {
        "label": label,
        "source_name": dataset.source_name,
        "feature_count": int(dataset.matrix.shape[0]),
        "sample_count": int(dataset.matrix.shape[1]),
        "baseline_stable": baseline_stable,
        "k_star_index": k_star_index,
        "k_star_time": k_star_time,
        "first_visible_failure_index": visible_index,
        "first_visible_failure_time": visible_time,
        "early_detection": early_detection,
        "max_fraction_abs_change_ge_1_0": max_fraction_abs_change,
        "scope_note": "External expression-pattern proxy only; not biological failure validation.",
    }


def run_single_bridge(dataset: ExpressionDataset, label: str) -> BridgeRun:
    standard_metrics = compute_standard_expression_metrics(dataset)
    rows = compute_recoverability_metrics(dataset)
    summary = summarize_bridge(dataset, rows, standard_metrics, label=label)
    return BridgeRun(label, dataset, rows, summary, standard_metrics)


def run_bridge_with_controls(dataset: ExpressionDataset) -> dict[str, object]:
    observed = run_single_bridge(dataset, "observed")
    controls = [
        run_single_bridge(make_time_shuffle_control(dataset), "time_shuffle_control"),
        run_single_bridge(make_feature_shuffle_control(dataset), "feature_shuffle_control"),
        run_single_bridge(make_row_permutation_control(dataset), "row_permutation_control"),
    ]

    control_early_count = sum(bool(control.summary["early_detection"]) for control in controls)
    observed_early = bool(observed.summary["early_detection"])
    contrast_pass = observed_early and control_early_count == 0
    if contrast_pass:
        bridge_status = "PROBE_PASS_WITH_CONTROLS"
    elif not observed_early:
        bridge_status = "PROBE_FAIL_NO_EARLY_DETECTION"
    else:
        bridge_status = "PROBE_FAIL_CONTROL_MATCH"

    return {
        "bridge_status": bridge_status,
        "contrast_pass": contrast_pass,
        "control_early_detection_count": control_early_count,
        "runs": [observed, *controls],
        "failure_rule": "The bridge fails if controls produce the same early-detection signature as observed data or observed data lacks early detection.",
    }


def write_feature_rankings_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "feature_id",
        "haos_loss_score",
        "max_abs_change",
        "mean_abs_change",
        "trajectory_total_variation",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_results_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_index",
        "sample_label",
        "time_value",
        "recoverability",
        "delta_persistence",
        "k_star",
        "k_star_time",
        "safety_margin",
        "confidence",
        "expression_match",
        "step_coherence",
        "trajectory_identity_retention",
        "trajectory_identity_failure",
        "support_retention",
        "visible_failure",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_standard_metrics_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_index",
        "sample_label",
        "time_value",
        "mean_abs_change",
        "median_abs_change",
        "fraction_abs_change_ge_0_5",
        "fraction_abs_change_ge_1_0",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_metric_decomposition_csv(
    path: Path,
    rows: list[dict[str, object]],
    standard_metrics: list[dict[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_index",
        "sample_label",
        "time_value",
        "recoverability",
        "expression_match",
        "step_coherence",
        "trajectory_identity_retention",
        "trajectory_identity_failure",
        "support_retention",
        "delta_persistence",
        "mean_abs_change",
        "median_abs_change",
        "fraction_abs_change_ge_0_5",
        "fraction_abs_change_ge_1_0",
        "visible_failure",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row, standard_row in zip(rows, standard_metrics):
            writer.writerow(
                {
                    "sample_index": row["sample_index"],
                    "sample_label": row["sample_label"],
                    "time_value": row["time_value"],
                    "recoverability": row["recoverability"],
                    "expression_match": row["expression_match"],
                    "step_coherence": row["step_coherence"],
                    "trajectory_identity_retention": row["trajectory_identity_retention"],
                    "trajectory_identity_failure": row["trajectory_identity_failure"],
                    "support_retention": row["support_retention"],
                    "delta_persistence": row["delta_persistence"],
                    "mean_abs_change": standard_row["mean_abs_change"],
                    "median_abs_change": standard_row["median_abs_change"],
                    "fraction_abs_change_ge_0_5": standard_row["fraction_abs_change_ge_0_5"],
                    "fraction_abs_change_ge_1_0": standard_row["fraction_abs_change_ge_1_0"],
                    "visible_failure": row["visible_failure"],
                }
            )


def write_threshold_sensitivity_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "collapse_threshold",
        "visible_match_threshold",
        "k_star_index",
        "k_star_time",
        "first_visible_index",
        "first_visible_time",
        "early_detection",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_control_summary_csv(path: Path, runs: list[BridgeRun]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "label",
        "source_name",
        "feature_count",
        "sample_count",
        "baseline_stable",
        "k_star_index",
        "k_star_time",
        "first_visible_failure_index",
        "first_visible_failure_time",
        "early_detection",
        "max_fraction_abs_change_ge_1_0",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for run in runs:
            writer.writerow({field: run.summary[field] for field in fieldnames})


def write_status_json(path: Path, status: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(status, handle, indent=2, sort_keys=True)
        handle.write("\n")
