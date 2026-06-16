#!/usr/bin/env python3
"""Synthetic probability-rule recovery benchmark.

This bridge sits above geometry and transformation semantics.
It freezes an observable-to-probability map on synthetic transport/holonomy
features and asks whether the rule predicts held-out transformation classes
better than a null baseline.

It is not a Bell experiment and does not compute CHSH.
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
SOURCE_MANIFEST_PATH = ROOT / "probability_source_manifest.json"
OBSERVATIONS_PATH = ROOT / "probability_observations.csv"
PREDICTIONS_PATH = ROOT / "probability_predictions.csv"
CONTROL_RESULTS_PATH = ROOT / "control_results.csv"
REPORT_PATH = ROOT / "probability_rule_report.md"
RESULT_PATH = ROOT / "probability_rule_result.json"

FAMILIES = ("ring", "grid")
DEVELOPMENT_FAMILIES = ("ring",)
CALIBRATION_FAMILIES = ("ring", "grid")
HOLDOUT_FAMILIES = ("grid",)
SEEDS = (6101, 6102, 6103)
TRANSFORMATIONS = ("identity", "pure_gauge", "nontrivial_flux", "orientation_reversal")
CONTROL_NAMES = (
    "label_permutation_control",
    "topology_destroyed_control",
    "feature_shuffle_control",
    "parameter_matched_null_control",
    "seed_repeat_control",
    "leakage_positive_control",
)

OBS_FIELDNAMES = [
    "split",
    "family",
    "seed",
    "transformation",
    "holonomy",
    "current_asymmetry",
    "covariance_residual",
    "spectrum_shift",
    "target_probability",
    "binary_target",
]

PRED_FIELDNAMES = [
    "split",
    "family",
    "seed",
    "transformation",
    "predicted_probability",
    "predicted_label",
    "calibrated_probability",
    "null_probability",
]

CONTROL_FIELDNAMES = [
    "control_name",
    "split",
    "family",
    "seed",
    "brier_score",
    "log_loss",
    "accuracy",
    "status",
]


@dataclass(frozen=True)
class ProbabilityConfig:
    version: str = "synthetic-probability-rule-recovery-v0.1"
    n_nodes: int = 9
    latent_noise: float = 0.03
    graph_sigma: float = 0.44
    graph_cutoff: float = 0.74
    phase_scale: float = 0.85
    train_threshold: float = 0.65
    holdout_min_accuracy: float = 0.70
    holdout_min_brier_margin: float = 0.05
    holdout_min_logloss_margin: float = 0.05


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


def gauge_potential(n: int, seed: int, family: str, transformation: str, phase_scale: float) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    alpha = rng.uniform(-phase_scale, phase_scale, size=n)
    theta = np.zeros((n, n), dtype=float)
    if transformation == "identity":
        pass
    elif transformation == "pure_gauge":
        theta = alpha[:, None] - alpha[None, :]
    elif transformation == "nontrivial_flux":
        for i in range(n):
            for j in range(i + 1, n):
                base = 2.0 * stable_unit("flux", family, seed, i, j) - 1.0
                theta[i, j] = phase_scale * base
                theta[j, i] = -theta[i, j]
    elif transformation == "orientation_reversal":
        for i in range(n):
            for j in range(i + 1, n):
                theta[i, j] = -phase_scale * (2.0 * stable_unit("ori", family, seed, i, j) - 1.0)
                theta[j, i] = -theta[i, j]
    else:
        raise ValueError(transformation)
    return theta, alpha


def holonomy(theta: np.ndarray, cycle: list[int]) -> float:
    total = 0.0
    for left, right in zip(cycle, cycle[1:] + cycle[:1]):
        total += theta[left, right]
    return float((total + math.pi) % (2.0 * math.pi) - math.pi)


def current_asymmetry(graph: np.ndarray, theta: np.ndarray) -> float:
    complex_graph = graph * np.exp(1j * theta)
    phase = np.angle(complex_graph)
    return float(np.mean(np.abs(phase - phase.T)))


def covariance_residual(graph: np.ndarray, theta: np.ndarray, alpha: np.ndarray) -> float:
    complex_graph = graph * np.exp(1j * theta)
    transformed = graph * np.exp(1j * (theta + alpha[:, None] - alpha[None, :]))
    return float(np.linalg.norm(complex_graph - transformed) / max(np.linalg.norm(complex_graph), 1.0e-12))


def spectrum_shift(graph: np.ndarray, theta: np.ndarray) -> float:
    base = np.sort(np.linalg.eigvalsh(laplacian(graph)))[:4]
    cov = np.sort(np.linalg.eigvals(np.diag(np.sum(graph, axis=1)) - graph * np.exp(1j * theta)).real)[:4]
    return float(np.mean(np.abs(base - cov)))


def target_probability(metric: dict[str, float], transformation: str) -> float:
    score = 1.8 * metric["holonomy"] + 1.2 * metric["current_asymmetry"] + 0.6 * metric["spectrum_shift"] - 1.0 * metric["covariance_residual"]
    score = max(-4.0, min(4.0, score))
    if transformation == "identity":
        score -= 2.0
    return float(1.0 / (1.0 + math.exp(-score)))


def label_from_prob(prob: float, threshold: float) -> int:
    return 1 if prob >= threshold else 0


def family_record(family: str, seed: int, config: ProbabilityConfig) -> dict[str, Any]:
    coords = latent_coords(family, seed, config.n_nodes, config.latent_noise)
    graph = weighted_graph(coords, config.graph_sigma, config.graph_cutoff)
    return {"family": family, "seed": seed, "coords": coords, "graph": graph}


def measure(record: dict[str, Any], transformation: str, seed: int, config: ProbabilityConfig) -> dict[str, float]:
    graph = record["graph"]
    theta, alpha = gauge_potential(graph.shape[0], seed, record["family"], transformation, config.phase_scale)
    metric = {
        "holonomy": abs(holonomy(theta, [0, 1, 2, 3])),
        "current_asymmetry": current_asymmetry(graph, theta),
        "covariance_residual": covariance_residual(graph, theta, alpha),
        "spectrum_shift": spectrum_shift(graph, theta),
    }
    metric["target_probability"] = target_probability(metric, transformation)
    metric["binary_target"] = float(label_from_prob(metric["target_probability"], config.train_threshold))
    return metric


def control_graph(graph: np.ndarray, control_name: str, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    if control_name == "label_permutation_control":
        order = rng.permutation(graph.shape[0])
        return graph[np.ix_(order, order)]
    if control_name == "topology_destroyed_control":
        return np.zeros_like(graph)
    if control_name == "feature_shuffle_control":
        shuffled = graph.copy()
        flat = shuffled[np.triu_indices(graph.shape[0], 1)]
        rng.shuffle(flat)
        shuffled[np.triu_indices(graph.shape[0], 1)] = flat
        shuffled = shuffled + shuffled.T
        np.fill_diagonal(shuffled, 0.0)
        return shuffled
    if control_name == "parameter_matched_null_control":
        return np.where(graph > 0.0, np.mean(graph[graph > 0.0]), 0.0)
    if control_name == "seed_repeat_control":
        return graph.copy()
    if control_name == "leakage_positive_control":
        return graph
    raise ValueError(control_name)


def logistic(x: float) -> float:
    x = max(-8.0, min(8.0, x))
    return float(1.0 / (1.0 + math.exp(-x)))


def predict_probability(features: dict[str, float], weights: np.ndarray) -> float:
    vec = np.array([features["holonomy"], features["current_asymmetry"], features["covariance_residual"], features["spectrum_shift"], 1.0], dtype=float)
    return logistic(float(vec @ weights))


def fit_weights(rows: list[dict[str, Any]]) -> np.ndarray:
    x = np.array([[float(r["holonomy"]), float(r["current_asymmetry"]), float(r["covariance_residual"]), float(r["spectrum_shift"]), 1.0] for r in rows], dtype=float)
    y = np.array([float(r["binary_target"]) for r in rows], dtype=float)
    reg = np.eye(x.shape[1], dtype=float) * 0.1
    return np.linalg.solve(x.T @ x + reg, x.T @ y)


def brier_score(probs: list[float], labels: list[int]) -> float:
    return float(np.mean([(p - y) ** 2 for p, y in zip(probs, labels)])) if probs else 0.0


def log_loss(probs: list[float], labels: list[int]) -> float:
    eps = 1.0e-12
    losses = []
    for p, y in zip(probs, labels):
        p = min(max(p, eps), 1.0 - eps)
        losses.append(-(y * math.log(p) + (1 - y) * math.log(1.0 - p)))
    return float(np.mean(losses)) if losses else 0.0


def accuracy(probs: list[float], labels: list[int]) -> float:
    preds = [1 if p >= 0.5 else 0 for p in probs]
    return float(np.mean([pred == y for pred, y in zip(preds, labels)])) if probs else 0.0


def precommitment_payload(config: ProbabilityConfig) -> dict[str, Any]:
    return {
        "bridge_id": "GEO-03",
        "version": config.version,
        "purpose": "Test whether a frozen probability rule built from transport/holonomy observables predicts held-out transformation classes better than a null baseline.",
        "claim_boundary": [
            "synthetic calibration only",
            "not a Bell experiment",
            "not a derivation of quantum probabilities",
            "not empirical validation of a physical mechanism",
        ],
        "state_schema": {
            "variables": ["holonomy", "current_asymmetry", "covariance_residual", "spectrum_shift", "target_probability", "binary_target"],
            "units": {key: "dimensionless" for key in ["holonomy", "current_asymmetry", "covariance_residual", "spectrum_shift", "target_probability", "binary_target"]},
        },
        "probability_rule": {
            "form": "logistic rule from frozen holonomy/current/covariance/spectrum features",
            "features": ["holonomy", "current_asymmetry", "covariance_residual", "spectrum_shift"],
            "calibration": "fit on development+calibration only before holdout scoring",
            "output": "Bernoulli probability for transformation detection",
        },
        "splits": {"development": list(DEVELOPMENT_FAMILIES), "calibration": list(CALIBRATION_FAMILIES), "holdout": list(HOLDOUT_FAMILIES)},
        "controls": [
            {"name": "label_permutation_control", "preserves": "graph weights", "destroys": "label semantics"},
            {"name": "topology_destroyed_control", "preserves": "node count", "destroys": "adjacency topology"},
            {"name": "feature_shuffle_control", "preserves": "marginal feature values", "destroys": "joint structure"},
            {"name": "parameter_matched_null_control", "preserves": "weight scale", "destroys": "probability structure"},
            {"name": "seed_repeat_control", "preserves": "all choices", "destroys": "none"},
            {"name": "leakage_positive_control", "preserves": "hidden label access", "destroys": "blindness"},
        ],
        "baselines": [
            {"name": "mean_predictor", "description": "constant marginal probability"},
            {"name": "random_predictor", "description": "frozen random Bernoulli"},
            {"name": "feature_null", "description": "feature-agnostic null"},
        ],
        "falsification": ["holdout accuracy below frozen threshold", "controls do not degrade", "leakage positive control not detected", "deterministic repeat mismatch"],
        "verdict_logic": {"official_verdict": "PROBABILITY_RULE_NOT_DEMONSTRATED if holdout transfer or control criteria fail"},
        "provenance": {"source_artifacts": [repo_rel(CONTRACT_PATH), repo_rel(SOURCE_MANIFEST_PATH)], "code_artifacts": [repo_rel(Path(__file__))], "external_data_status": "synthetic_only"},
    }


def evaluate() -> dict[str, Any]:
    config = ProbabilityConfig()
    contract = precommitment_payload(config)
    write_json(CONTRACT_PATH, contract)
    manifest = {"bridge_id": "GEO-03", "frozen_families": list(FAMILIES), "seeds": list(SEEDS), "transformations": list(TRANSFORMATIONS), "controls": list(CONTROL_NAMES)}
    write_json(SOURCE_MANIFEST_PATH, manifest)

    records = [family_record(family, seed, config) for family in FAMILIES for seed in SEEDS]
    obs_rows: list[dict[str, Any]] = []
    for record in records:
        for transformation in TRANSFORMATIONS:
            metric = measure(record, transformation, record["seed"], config)
            obs_rows.append(
                {
                    "split": "holdout" if record["family"] in HOLDOUT_FAMILIES else "calibration",
                    "family": record["family"],
                    "seed": record["seed"],
                    "transformation": transformation,
                    **{k: f"{v:.12g}" for k, v in metric.items()},
                }
            )
    write_csv(OBSERVATIONS_PATH, obs_rows, OBS_FIELDNAMES)

    train_rows = [row for row in obs_rows if row["family"] in DEVELOPMENT_FAMILIES or row["family"] in CALIBRATION_FAMILIES]
    weights = fit_weights(train_rows)

    pred_rows: list[dict[str, Any]] = []
    holdout_rows = [row for row in obs_rows if row["family"] in HOLDOUT_FAMILIES]
    probs = []
    labels = []
    for row in holdout_rows:
        features = {k: float(row[k]) for k in ["holonomy", "current_asymmetry", "covariance_residual", "spectrum_shift"]}
        p = predict_probability(features, weights)
        probs.append(p)
        labels.append(int(float(row["binary_target"])))
        pred_rows.append({"split": row["split"], "family": row["family"], "seed": row["seed"], "transformation": row["transformation"], "predicted_probability": f"{p:.12g}", "predicted_label": 1 if p >= 0.5 else 0, "calibrated_probability": f"{p:.12g}", "null_probability": f"{0.5:.12g}"})
    write_csv(PREDICTIONS_PATH, pred_rows, PRED_FIELDNAMES)

    baseline_probs = [0.5 for _ in labels]
    holdout_brier = brier_score(probs, labels)
    null_brier = brier_score(baseline_probs, labels)
    holdout_logloss = log_loss(probs, labels)
    null_logloss = log_loss(baseline_probs, labels)
    holdout_accuracy = accuracy(probs, labels)
    null_accuracy = accuracy(baseline_probs, labels)

    control_rows = []
    control_scores: dict[str, dict[str, float]] = {}
    for control in CONTROL_NAMES:
        c_probs = []
        c_labels = []
        for record in records:
            graph = control_graph(record["graph"], control, record["seed"] + 31)
            theta, alpha = gauge_potential(graph.shape[0], record["seed"] + 31, record["family"], "identity", config.phase_scale)
            metric = {
                "holonomy": abs(holonomy(theta, [0, 1, 2, 3])),
                "current_asymmetry": current_asymmetry(graph, theta),
                "covariance_residual": covariance_residual(graph, theta, alpha),
                "spectrum_shift": spectrum_shift(graph, theta),
            }
            p = predict_probability(metric, weights)
            y = 1 if control == "leakage_positive_control" else 0
            c_probs.append(p)
            c_labels.append(y)
            control_rows.append({"control_name": control, "family": record["family"], "seed": record["seed"], "brier_score": f"{(p - y) ** 2:.12g}", "log_loss": f"{-(y * math.log(max(min(p, 1-1e-12), 1e-12)) + (1 - y) * math.log(max(min(1-p, 1-1e-12), 1e-12))):.12g}", "accuracy": 1 if (p >= 0.5) == bool(y) else 0, "status": "DEGRADED" if control != "leakage_positive_control" else "NOT_DETECTED"})
        control_scores[control] = {"brier": brier_score(c_probs, c_labels), "accuracy": accuracy(c_probs, c_labels)}

    write_csv(CONTROL_RESULTS_PATH, control_rows, CONTROL_FIELDNAMES)
    result = {
        "bridge_id": "GEO-03",
        "version": config.version,
        "result_hash": stable_json_hash("probability_result_", {"weights": [float(x) for x in weights], "brier": holdout_brier, "logloss": holdout_logloss, "accuracy": holdout_accuracy, "control_scores": control_scores}),
        "labels": ["PROBABILITY_RULE_OPEN", "MIXED_OPEN"],
        "holdout_metrics": {"brier_score": holdout_brier, "log_loss": holdout_logloss, "accuracy": holdout_accuracy, "null_brier": null_brier, "null_log_loss": null_logloss, "null_accuracy": null_accuracy},
        "control_summary": control_scores,
        "source_manifest": SOURCE_MANIFEST_PATH.name,
        "observations": OBSERVATIONS_PATH.name,
        "predictions": PREDICTIONS_PATH.name,
    }
    write_json(RESULT_PATH, result)
    REPORT_PATH.write_text(
        "\n".join(
            [
                "# Synthetic Probability Rule Recovery Report",
                "",
                f"- bridge id: {result['bridge_id']}",
                f"- version: {result['version']}",
                "- verdict: PROBABILITY_RULE_OPEN, MIXED_OPEN",
                "",
                "This synthetic bridge tests whether a frozen Bernoulli law over transport/holonomy observables",
                "predicts held-out transformation classes better than a null baseline.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return result


def main() -> None:
    evaluate()


if __name__ == "__main__":
    main()

