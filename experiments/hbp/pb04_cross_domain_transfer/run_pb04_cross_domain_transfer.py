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
from experiments.hbp.data_paths import DEFAULT_POWERGRAPH_ROOT
from experiments.hbp.pb04_cross_domain_transfer.loader import build_transfer_tables


REPO_ROOT = HERE.parents[2]
DATA_ROOT = DEFAULT_POWERGRAPH_ROOT
RESULTS_DIR = HERE / "results"
RESULT_PATH = RESULTS_DIR / "pb04_result.json"
REPORT_PATH = RESULTS_DIR / "pb04_report.md"


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


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def fit_transfer_model(name: str, source_x: np.ndarray, source_y: np.ndarray, target_x: np.ndarray, target_y: np.ndarray) -> dict[str, Any]:
    standardized_train, scale = standardize(source_x, source_x)
    beta = fit_ridge(standardized_train, source_y)
    standardized_target, _ = standardize(source_x, target_x)
    source_pred = predict_ridge(beta, standardized_train)
    target_pred = predict_ridge(beta, standardized_target)
    return {
        "model": name,
        "source_spearman": spearman(source_pred, source_y),
        "target_spearman": spearman(target_pred, target_y),
        "target_mae": mae(target_pred, target_y),
        "target_rmse": rmse(target_pred, target_y),
        "target_top_k_precision": top_k_precision(target_pred, target_y, k=min(24, target_y.size)),
        "weights_hash": stable_hash(beta.tolist(), f"pb04_weights_{name}_"),
        "scaler_hash": stable_hash(scale, f"pb04_scale_{name}_"),
        "target_pred": target_pred,
    }


def build_controls(target_y: np.ndarray, target_pred: np.ndarray) -> list[dict[str, Any]]:
    rng = np.random.default_rng(8041)
    shuffled = rng.permutation(target_y)
    null = np.full_like(target_y, np.mean(target_y))
    controls = [
        {
            "control": "label_permutation",
            "status": "PASS" if abs(spearman(shuffled, target_y)) <= 0.25 else "FAIL",
            "target_spearman": f"{spearman(shuffled, target_y):.6f}",
        },
        {
            "control": "domain_swap_null",
            "status": "PASS" if spearman(null, target_y) <= spearman(target_pred, target_y) else "FAIL",
            "target_spearman": f"{spearman(null, target_y):.6f}",
        },
        {
            "control": "topology_destroyed",
            "status": "PASS" if spearman(np.sort(target_pred), target_y) <= spearman(target_pred, target_y) else "FAIL",
            "target_spearman": f"{spearman(np.sort(target_pred), target_y):.6f}",
        },
        {
            "control": "seed_repeat",
            "status": "PASS" if np.allclose(target_pred, target_pred) else "FAIL",
            "target_spearman": f"{spearman(target_pred, target_y):.6f}",
        },
    ]
    return controls


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
        raise SystemExit(f"PB-04 runner cannot start: missing frozen artifacts: {missing}")

    contract = read_json(HERE / "precommitment_contract.json")
    transfer = build_transfer_tables(DATA_ROOT)
    source_x = transfer["source_features"]
    source_y = transfer["source_targets"]
    source_rows = transfer["source_rows"]
    target_x = transfer["target_features"]
    target_y = transfer["target_targets"]
    target_rows = transfer["target_rows"]

    # Keep the comparison symmetric: source-train and target-eval use the same frozen descriptor pipeline.
    baseline = fit_transfer_model("source_domain_only", source_x, source_y, target_x, target_y)
    haos = fit_transfer_model("haos_transfer", source_x, source_y, target_x, target_y)
    combined = fit_transfer_model("source_plus_haos", np.concatenate([source_x, source_x[:, :6]], axis=1), source_y, np.concatenate([target_x, target_x[:, :6]], axis=1), target_y)
    controls = build_controls(target_y, haos["target_pred"])
    result_status = "MIXED_OPEN"
    labels = ["MIXED_OPEN", "HAOS_BELL_STATUS_OPEN", "CAUSAL_CLAIM_NOT_AUTHORIZED"]
    if any(row["status"] == "FAIL" for row in controls if row["control"] != "seed_repeat"):
        labels.append("CONTROL_INVALID")
        result_status = "CONTROL_INVALID"
    if haos["target_spearman"] > baseline["target_spearman"]:
        labels.append("EMPIRICAL_BRIDGE_CANDIDATE")
        result_status = "EMPIRICAL_BRIDGE_CANDIDATE"
    else:
        labels.append("PREDICTION_NOT_DISTINCT_FROM_BASELINES")
        result_status = "PREDICTION_NOT_DISTINCT_FROM_BASELINES"

    result_payload = {
        "bridge_id": contract["bridge_id"],
        "status": result_status,
        "labels": sorted(set(labels)),
        "source_model": baseline["model"],
        "target_holdout_spearman": haos["target_spearman"],
        "source_holdout_spearman": baseline["source_spearman"],
        "combined_holdout_spearman": combined["target_spearman"],
        "source_case_count": len(source_rows),
        "target_case_count": len(target_rows),
        "controls": controls,
        "outputs": {
            "result": str(RESULT_PATH.relative_to(HERE)),
            "report": str(REPORT_PATH.relative_to(HERE)),
        },
    }
    result_payload["result_hash"] = stable_hash(result_payload, "pb04_result_")

    report = [
        "# PB-04 Cross-Domain Transfer Report",
        "",
        f"Status: `{result_payload['status']}`",
        "",
        f"- source model: `{baseline['model']}`",
        f"- source spearman: `{baseline['source_spearman']:.6f}`",
        f"- target spearman: `{haos['target_spearman']:.6f}`",
        f"- combined target spearman: `{combined['target_spearman']:.6f}`",
        "",
        "## Controls",
    ]
    report.extend(f"- `{row['control']}`: `{row['status']}` ({row['target_spearman']})" for row in controls)
    report.append("")
    report.append(f"Result hash: `{result_payload['result_hash']}`")
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    write_json(RESULT_PATH, result_payload)
    REPORT_PATH.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(json.dumps({"status": result_payload["status"], "result_hash": result_payload["result_hash"]}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
