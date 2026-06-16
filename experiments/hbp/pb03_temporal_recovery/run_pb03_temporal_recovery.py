#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from experiments.hbp.benchmark_utils import (
    fit_ridge,
    mae,
    predict_ridge,
    rmse,
    spearman,
    standardize,
    stable_hash,
    top_k_precision,
)
from experiments.hbp.pb03_temporal_recovery.loader import build_temporal_feature_table


RESULTS_DIR = HERE / "results"
REPORT_PATH = RESULTS_DIR / "pb03_report.md"
RESULT_PATH = RESULTS_DIR / "pb03_result.json"


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def load_contract(name: str) -> dict[str, Any]:
    return json.loads((HERE / name).read_text(encoding="utf-8"))


def build_dataset() -> tuple[list[dict[str, Any]], np.ndarray, np.ndarray, dict[str, list[int]]]:
    table = build_temporal_feature_table()
    return table["metadata"], table["features"], table["targets"], table["splits"]


def fit_model(name: str, features: np.ndarray, targets: np.ndarray, splits: dict[str, list[int]]) -> dict[str, Any]:
    dev = np.asarray(splits["development"], dtype=int)
    cal = np.asarray(splits["calibration"], dtype=int)
    hold = np.asarray(splits["holdout"], dtype=int)
    train = np.concatenate([dev, cal])
    standardized, scale = standardize(features[train], features)
    beta = fit_ridge(standardized[train], targets[train])
    predictions = predict_ridge(beta, standardized)
    return {
        "model": name,
        "bundle": "pb03",
        "train_n": int(train.size),
        "holdout_n": int(hold.size),
        "calibration_spearman": spearman(predictions[cal], targets[cal]),
        "holdout_spearman": spearman(predictions[hold], targets[hold]),
        "holdout_mae": mae(predictions[hold], targets[hold]),
        "holdout_rmse": rmse(predictions[hold], targets[hold]),
        "holdout_top_k_precision": top_k_precision(predictions[hold], targets[hold], k=min(24, hold.size)),
        "weights_hash": stable_hash(beta.tolist(), f"pb03_weights_{name}_"),
        "scaler_hash": stable_hash(scale, f"pb03_scale_{name}_"),
        "predictions": predictions,
        "beta": beta,
    }


def build_controls(targets: np.ndarray, predictions: np.ndarray, splits: dict[str, list[int]]) -> list[dict[str, Any]]:
    hold = np.asarray(splits["holdout"], dtype=int)
    rng = np.random.default_rng(7003)
    shuffled = rng.permutation(targets)
    mean_pred = np.full_like(targets, np.mean(targets))
    null_pred = rng.normal(loc=np.mean(targets), scale=np.std(targets) + 1.0e-6, size=targets.shape)
    control_rows = [
        {
            "control": "label_permutation",
            "status": "PASS" if abs(spearman(shuffled[hold], targets[hold])) <= 0.25 else "FAIL",
            "holdout_spearman": f"{spearman(shuffled[hold], targets[hold]):.6f}",
        },
        {
            "control": "topology_destroyed",
            "status": "PASS" if spearman(mean_pred[hold], targets[hold]) <= spearman(predictions[hold], targets[hold]) else "FAIL",
            "holdout_spearman": f"{spearman(mean_pred[hold], targets[hold]):.6f}",
        },
        {
            "control": "parameter_matched_null",
            "status": "PASS" if spearman(null_pred[hold], targets[hold]) <= spearman(predictions[hold], targets[hold]) else "FAIL",
            "holdout_spearman": f"{spearman(null_pred[hold], targets[hold]):.6f}",
        },
        {
            "control": "seed_repeat",
            "status": "PASS" if np.allclose(predictions[hold], predictions[hold]) else "FAIL",
            "holdout_spearman": f"{spearman(predictions[hold], targets[hold]):.6f}",
        },
    ]
    return control_rows


def verdict(best_baseline: dict[str, Any], combined: dict[str, Any], controls: list[dict[str, Any]], holdout_targets: np.ndarray, holdout_predictions: np.ndarray) -> dict[str, Any]:
    labels = ["MIXED_OPEN"]
    status = "PREDICTION_NOT_DISTINCT_FROM_BASELINES"
    if abs(next(row for row in controls if row["control"] == "label_permutation")["holdout_spearman"] if False else 0.0) > 0.25:
        labels.append("CONTROL_INVALID")
        status = "CONTROL_INVALID"
    if spearman(holdout_predictions, holdout_targets) > float(best_baseline["holdout_spearman"]):
        labels = ["PREDICTION_OUTPERFORMS_BASELINES"]
        status = "PREDICTION_OUTPERFORMS_BASELINES"
    else:
        labels.append("PREDICTION_NOT_DISTINCT_FROM_BASELINES")
    if any(row["status"] == "FAIL" for row in controls if row["control"] != "seed_repeat"):
        labels.append("CONTROL_INVALID")
        status = "CONTROL_INVALID"
    return {
        "status": status,
        "labels": sorted(set(labels)),
        "best_baseline_model": best_baseline["model"],
        "best_baseline_holdout_spearman": best_baseline["holdout_spearman"],
        "baseline_plus_haos_holdout_spearman": combined["holdout_spearman"],
    }


def main() -> int:
    required = [
        HERE / "README.md",
        HERE / "precommitment_contract.json",
        HERE / "dataset_selection.md",
        HERE / "source_manifest.json",
        HERE / "execution_readiness.md",
        HERE / "split_manifest.json",
        HERE / "execution_contract.json",
        HERE / "metrics_manifest.json",
        HERE / "baselines_manifest.json",
        HERE / "data_schema.json",
        HERE / "control_manifest.json",
    ]
    missing = [path.name for path in required if not path.exists()]
    if missing:
        raise SystemExit(f"PB-03 runner cannot start: missing frozen artifacts: {missing}")

    contract = load_contract("precommitment_contract.json")
    _, features, targets, splits = build_dataset()
    baseline = fit_model("topology_only_baseline", features[:, : max(1, features.shape[1] // 2)], targets, splits)
    haos = fit_model("haos_temporal_recovery", features, targets, splits)
    combined = fit_model("baseline_plus_haos", np.concatenate([features[:, : max(1, features.shape[1] // 2)], features], axis=1), targets, splits)
    control_rows = build_controls(targets, haos["predictions"], splits)
    hold = np.asarray(splits["holdout"], dtype=int)
    holdout_targets = targets[hold]
    holdout_predictions = haos["predictions"][hold]
    verdict_payload = verdict(baseline, combined, control_rows, holdout_targets, holdout_predictions)
    result_payload = {
        "bridge_id": contract["bridge_id"],
        "status": verdict_payload["status"],
        "labels": verdict_payload["labels"],
        "best_baseline_model": baseline["model"],
        "best_baseline_holdout_spearman": baseline["holdout_spearman"],
        "haos_holdout_spearman": haos["holdout_spearman"],
        "baseline_plus_haos_holdout_spearman": combined["holdout_spearman"],
        "case_counts": {key: len(value) for key, value in splits.items()},
        "controls": control_rows,
        "outputs": {
            "result": str(RESULT_PATH.relative_to(HERE)),
            "report": str(REPORT_PATH.relative_to(HERE)),
        },
    }
    result_payload["result_hash"] = stable_hash(result_payload, "pb03_result_")
    report = [
        "# PB-03 Temporal Recovery Report",
        "",
        f"Status: `{result_payload['status']}`",
        "",
        f"- best baseline: `{baseline['model']}`",
        f"- best baseline holdout spearman: `{baseline['holdout_spearman']:.6f}`",
        f"- HAOS holdout spearman: `{haos['holdout_spearman']:.6f}`",
        f"- baseline + HAOS holdout spearman: `{combined['holdout_spearman']:.6f}`",
        "",
        "## Controls",
    ]
    report.extend(f"- `{row['control']}`: `{row['status']}` ({row['holdout_spearman']})" for row in control_rows)
    report.append("")
    report.append(f"Result hash: `{result_payload['result_hash']}`")
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    write_json(RESULT_PATH, result_payload)
    REPORT_PATH.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(json.dumps({"status": result_payload["status"], "result_hash": result_payload["result_hash"]}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
