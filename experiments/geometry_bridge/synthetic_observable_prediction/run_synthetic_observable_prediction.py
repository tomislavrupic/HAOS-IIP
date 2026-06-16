#!/usr/bin/env python3
"""Synthetic observable-prediction benchmark.

This bridge turns the frozen probability rule into a prediction about a held-
out observable: whether a transformation produces a nontrivial transport
response.

It is synthetic, deterministic, and not a Bell experiment.
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
SOURCE_MANIFEST_PATH = ROOT / "observable_source_manifest.json"
OBSERVATIONS_PATH = ROOT / "observable_observations.csv"
PREDICTIONS_PATH = ROOT / "observable_predictions.csv"
REPORT_PATH = ROOT / "observable_prediction_report.md"
RESULT_PATH = ROOT / "observable_prediction_result.json"

FAMILIES = ("ring", "grid")
SEEDS = (7101, 7102, 7103)
TRANSFORMATIONS = ("identity", "pure_gauge", "nontrivial_flux", "orientation_reversal")

FIELDNAMES = [
    "split",
    "family",
    "seed",
    "transformation",
    "observable_value",
    "observable_label",
    "predicted_probability",
    "predicted_label",
]


@dataclass(frozen=True)
class ObservableConfig:
    version: str = "synthetic-observable-prediction-v0.1"
    n_nodes: int = 9
    latent_noise: float = 0.03
    graph_sigma: float = 0.44
    graph_cutoff: float = 0.74
    phase_scale: float = 0.85
    observable_threshold: float = 0.0


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


def holonomy_like(graph: np.ndarray, seed: int) -> float:
    n = graph.shape[0]
    cycle = [0, 1, 2, 3]
    total = 0.0
    for left, right in zip(cycle, cycle[1:] + cycle[:1]):
        total += graph[left, right] * (2.0 * stable_unit("obs", seed, left, right) - 1.0)
    return float(abs(total) / max(np.sum(graph), 1.0e-12))


def observable_value(family: str, transformation: str, graph: np.ndarray, seed: int) -> float:
    value = holonomy_like(graph, seed)
    if transformation == "identity":
        value -= 0.08
    elif transformation == "pure_gauge":
        value -= 0.02
    elif transformation == "nontrivial_flux":
        value += 0.11
    elif transformation == "orientation_reversal":
        value += 0.06
    if family == "grid":
        value += 0.015
    return float(value)


def predicted_probability(observable: float) -> float:
    return float(1.0 / (1.0 + math.exp(-10.0 * (observable - 0.2))))


def predict_label(probability: float) -> int:
    return 1 if probability >= 0.5 else 0


def make_record(family: str, seed: int, config: ObservableConfig) -> dict[str, Any]:
    coords = latent_coords(family, seed, config.n_nodes, config.latent_noise)
    graph = weighted_graph(coords, config.graph_sigma, config.graph_cutoff)
    return {"family": family, "seed": seed, "graph": graph}


def precommitment_payload(config: ObservableConfig) -> dict[str, Any]:
    return {
        "bridge_id": "GEO-04",
        "version": config.version,
        "purpose": "Test whether a frozen probability rule can predict a synthetic held-out observable derived from transport/holonomy structure.",
        "claim_boundary": [
            "synthetic calibration only",
            "not a Bell experiment",
            "not a physical prediction claim",
            "not empirical validation of a physical mechanism",
        ],
        "observable": {
            "name": "transformation_class",
            "definition": "which frozen transformation generated the observable response",
            "classes": list(TRANSFORMATIONS),
            "threshold_policy": "class centroids fit on development+calibration only",
        },
        "prediction_rule": "nearest centroid over frozen observable values",
        "provenance": {
            "source_artifacts": [repo_rel(CONTRACT_PATH), repo_rel(SOURCE_MANIFEST_PATH)],
            "code_artifacts": [repo_rel(Path(__file__))],
            "external_data_status": "synthetic_only",
        },
    }


def evaluate() -> dict[str, Any]:
    config = ObservableConfig()
    contract = precommitment_payload(config)
    write_json(CONTRACT_PATH, contract)
    write_json(SOURCE_MANIFEST_PATH, {"bridge_id": "GEO-04", "families": list(FAMILIES), "seeds": list(SEEDS), "transformations": list(TRANSFORMATIONS)})

    obs_rows = []
    pred_rows = []
    calibration_values = []
    all_rows: list[dict[str, Any]] = []
    for family in FAMILIES:
        for seed in SEEDS:
            record = make_record(family, seed, config)
            for transformation in TRANSFORMATIONS:
                obs = observable_value(family, transformation, record["graph"], seed)
                row = {
                    "split": "holdout" if family == "grid" else "calibration",
                    "family": family,
                    "seed": seed,
                    "transformation": transformation,
                    "observable_value": obs,
                }
                if family != "grid":
                    calibration_values.append(obs)
                all_rows.append(row)

    class_means = {}
    for transformation in TRANSFORMATIONS:
        class_rows = [row for row in all_rows if row["split"] == "calibration" and row["transformation"] == transformation]
        class_means[transformation] = float(np.mean([row["observable_value"] for row in class_rows])) if class_rows else 0.0

    y_true: list[str] = []
    y_pred: list[str] = []
    for row in all_rows:
        label = row["transformation"]
        diffs = {name: abs(row["observable_value"] - center) for name, center in class_means.items()}
        pred_label = min(diffs, key=diffs.get)
        prob = predicted_probability(row["observable_value"])
        pred = 1 if pred_label in {"nontrivial_flux", "orientation_reversal"} else 0
        y_true.append(label)
        y_pred.append(pred_label)
        rendered = {
            "split": row["split"],
            "family": row["family"],
            "seed": row["seed"],
            "transformation": row["transformation"],
            "observable_value": f"{row['observable_value']:.12g}",
            "observable_label": label,
            "predicted_probability": f"{prob:.12g}",
            "predicted_label": pred_label,
        }
        obs_rows.append(rendered)
        pred_rows.append(rendered)

    write_csv(OBSERVATIONS_PATH, obs_rows, FIELDNAMES)
    write_csv(PREDICTIONS_PATH, pred_rows, FIELDNAMES)
    holdout_rows = [row for row in obs_rows if row["split"] == "holdout"]
    holdout_y_true = [row["observable_label"] for row in holdout_rows]
    holdout_y_pred = [row["predicted_label"] for row in holdout_rows]
    holdout_probs = [float(row["predicted_probability"]) for row in holdout_rows]
    holdout_obs = [float(row["observable_value"]) for row in holdout_rows]
    accuracy = float(np.mean([a == b for a, b in zip(holdout_y_true, holdout_y_pred)])) if holdout_y_true else 0.0
    null_accuracy = max((holdout_y_true.count(name) / len(holdout_y_true)) for name in TRANSFORMATIONS) if holdout_y_true else 0.0
    observable_rmse = float(np.sqrt(np.mean([(p - o) ** 2 for p, o in zip(holdout_probs, holdout_obs)]))) if holdout_probs else 0.0

    holdout_transform_means = {
        transformation: float(
            np.mean(
                [
                    row["observable_value"]
                    for row in all_rows
                    if row["split"] == "holdout" and row["transformation"] == transformation
                ]
            )
        )
        for transformation in TRANSFORMATIONS
    }
    pairwise_total = 0
    pairwise_correct = 0
    for i, left in enumerate(TRANSFORMATIONS):
        for right in TRANSFORMATIONS[i + 1 :]:
            pairwise_total += 1
            if (holdout_transform_means[left] > holdout_transform_means[right]) == (
                class_means[left] > class_means[right]
            ):
                pairwise_correct += 1
    pairwise_accuracy = pairwise_correct / float(pairwise_total) if pairwise_total else 0.0
    result = {
        "bridge_id": "GEO-04",
        "version": config.version,
        "result_hash": stable_json_hash("observable_result_", {"accuracy": accuracy, "labels": holdout_y_true, "preds": holdout_y_pred, "class_means": class_means}),
        "labels": ["OBSERVABLE_PREDICTION_OPEN", "MIXED_OPEN"],
        "holdout_metrics": {
            "accuracy": accuracy,
            "null_accuracy": null_accuracy,
            "pairwise_accuracy": pairwise_accuracy,
            "pairwise_null_accuracy": 0.5,
            "observable_rmse": observable_rmse,
            "class_means": class_means,
            "holdout_transform_means": holdout_transform_means,
        },
    }
    write_json(RESULT_PATH, result)
    REPORT_PATH.write_text(
        "\n".join(
            [
                "# Synthetic Observable Prediction Report",
                "",
                f"- bridge id: {result['bridge_id']}",
                f"- version: {result['version']}",
                "- verdict: OBSERVABLE_PREDICTION_OPEN, MIXED_OPEN",
                "",
                "This synthetic bridge converts a frozen probability rule into an observable prediction for a held-out transport response label.",
                f"Calibration class means: {json.dumps(class_means, sort_keys=True)}",
                f"Holdout pairwise accuracy: {pairwise_accuracy:.6f}",
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
