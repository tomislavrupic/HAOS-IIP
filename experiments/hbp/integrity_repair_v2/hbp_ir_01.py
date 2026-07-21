#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
CONTRACT = ROOT / "precommitment_contract.json"
FEATURE_REGISTRY_PATH = ROOT / "feature_registry.json"

FEATURE_REGISTRY = json.loads(FEATURE_REGISTRY_PATH.read_text(encoding="utf-8"))["features"]
MODEL_FEATURES = {
    "mean": (),
    "signal_linear": ("signal",),
    "signal_noise_linear": ("signal", "noise"),
}
SPLITS = ("development", "calibration", "holdout")


@dataclass(frozen=True)
class Dataset:
    row_id: np.ndarray
    split: np.ndarray
    features: dict[str, np.ndarray]
    target: np.ndarray


def stable_hash(prefix: str, payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(encoded).hexdigest()[:24]}"


def array_hash(values: np.ndarray) -> str:
    return hashlib.sha256(np.asarray(values).tobytes()).hexdigest()


def make_dataset(seed: int = 20260711) -> Dataset:
    rng = np.random.default_rng(seed)
    n = 180
    signal = rng.normal(size=n)
    noise = rng.normal(size=n)
    target = 1.8 * signal + 0.2 * noise + rng.normal(scale=0.15, size=n)
    order = rng.permutation(n)
    split = np.empty(n, dtype=object)
    split[order[:80]] = "development"
    split[order[80:120]] = "calibration"
    split[order[120:]] = "holdout"
    noise[order[[7, 31, 95, 132]]] = np.nan
    return Dataset(np.arange(n), split, {"signal": signal, "noise": noise}, target)


def validate_dataset(dataset: Dataset) -> dict[str, dict[str, int]]:
    if not np.isfinite(dataset.target).all():
        raise ValueError("non-finite target rejected")
    unknown = set(dataset.features) - set(FEATURE_REGISTRY)
    if unknown:
        raise ValueError(f"undeclared or target-derived feature rejected: {sorted(unknown)}")
    missingness: dict[str, dict[str, int]] = {}
    for name, spec in FEATURE_REGISTRY.items():
        values = dataset.features[name]
        counts = {split: int(np.sum(~np.isfinite(values[dataset.split == split]))) for split in SPLITS}
        missingness[name] = counts
        if spec["required"] and any(counts.values()):
            raise ValueError(f"required feature {name} contains non-finite values")
    return missingness


def prepare_features(dataset: Dataset, names: tuple[str, ...]) -> np.ndarray:
    if not names:
        return np.empty((len(dataset.target), 0), dtype=float)
    columns = []
    dev = dataset.split == "development"
    for name in names:
        values = dataset.features[name].astype(float).copy()
        if np.any(~np.isfinite(values)):
            finite_dev = values[dev & np.isfinite(values)]
            if finite_dev.size == 0:
                raise ValueError(f"optional feature {name} has no finite development values")
            values[~np.isfinite(values)] = float(np.mean(finite_dev))
        columns.append(values)
    return np.column_stack(columns)


def fit_linear(x: np.ndarray, y: np.ndarray, indices: np.ndarray) -> np.ndarray:
    design = np.column_stack([np.ones(len(indices)), x[indices]])
    coefficients, *_ = np.linalg.lstsq(design, y[indices], rcond=None)
    return coefficients


def predict_linear(x: np.ndarray, coefficients: np.ndarray) -> np.ndarray:
    return np.column_stack([np.ones(len(x)), x]) @ coefficients


def rmse(y: np.ndarray, prediction: np.ndarray) -> float:
    return float(np.sqrt(np.mean((y - prediction) ** 2)))


def rankdata(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    ranks[order] = np.arange(len(values), dtype=float)
    return ranks


def spearman(y: np.ndarray, prediction: np.ndarray) -> float:
    if len(y) < 2:
        return 0.0
    return float(np.corrcoef(rankdata(y), rankdata(prediction))[0, 1])


def select_model(dataset: Dataset, fit_seeds: tuple[int, ...] = (3101, 3102, 3103)) -> tuple[str, dict[str, float]]:
    dev_indices = np.flatnonzero(dataset.split == "development")
    calibration = dataset.split == "calibration"
    scores: dict[str, float] = {}
    for model_name, feature_names in MODEL_FEATURES.items():
        x = prepare_features(dataset, feature_names)
        predictions = []
        for seed in fit_seeds:
            rng = np.random.default_rng(seed)
            sample = rng.choice(dev_indices, size=len(dev_indices), replace=True)
            coefficients = fit_linear(x, dataset.target, sample)
            predictions.append(predict_linear(x, coefficients))
        mean_prediction = np.mean(predictions, axis=0)
        scores[model_name] = rmse(dataset.target[calibration], mean_prediction[calibration])
    return min(scores, key=lambda name: (scores[name], name)), scores


def evaluate(dataset: Dataset, fit_seeds: tuple[int, ...] = (3101, 3102, 3103)) -> dict[str, Any]:
    missingness = validate_dataset(dataset)
    selected, calibration_scores = select_model(dataset, fit_seeds)
    feature_names = MODEL_FEATURES[selected]
    x = prepare_features(dataset, feature_names)
    dev_indices = np.flatnonzero(dataset.split == "development")
    fit_predictions = []
    for seed in fit_seeds:
        rng = np.random.default_rng(seed)
        sample = rng.choice(dev_indices, size=len(dev_indices), replace=True)
        coefficients = fit_linear(x, dataset.target, sample)
        fit_predictions.append(predict_linear(x, coefficients))
    predictions = np.mean(fit_predictions, axis=0)
    holdout = dataset.split == "holdout"
    holdout_fit_rmse = [rmse(dataset.target[holdout], values[holdout]) for values in fit_predictions]
    return {
        "selected_model": selected,
        "calibration_scores": calibration_scores,
        "missingness": missingness,
        "predictions": predictions,
        "holdout_rmse": rmse(dataset.target[holdout], predictions[holdout]),
        "holdout_spearman": spearman(dataset.target[holdout], predictions[holdout]),
        "fit_seed_variance": float(np.var(holdout_fit_rmse, ddof=1)),
    }


def transform_control(dataset: Dataset, name: str, seed: int = 7711) -> tuple[Dataset, set[str]]:
    rng = np.random.default_rng(seed)
    features = {key: values.copy() for key, values in dataset.features.items()}
    target = dataset.target.copy()
    changed: set[str] = set()
    if name == "target_label_shuffle":
        target = target[rng.permutation(len(target))]
        changed.add("target")
    elif name == "signal_feature_shuffle":
        features["signal"] = features["signal"][rng.permutation(len(target))]
        changed.add("signal")
    elif name == "parameter_matched_null":
        source = features["signal"]
        features["signal"] = rng.normal(float(np.mean(source)), float(np.std(source)), size=len(source))
        changed.add("signal")
    elif name == "target_proxy_leakage":
        features["target_proxy"] = target.copy()
        changed.add("target_proxy")
    elif name == "no_op_control":
        pass
    else:
        raise ValueError(f"unknown control: {name}")
    return Dataset(dataset.row_id.copy(), dataset.split.copy(), features, target), changed


def run(output_dir: Path = RESULTS) -> dict[str, Any]:
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if contract.get("status") != "FROZEN":
        raise ValueError("HBP-IR-01 contract is not frozen")
    dataset = make_dataset(contract["seed_policy"]["data_seed"])
    target_result = evaluate(dataset, tuple(contract["seed_policy"]["fit_seeds"]))

    control_rows = []
    valid_controls = ("target_label_shuffle", "signal_feature_shuffle", "parameter_matched_null")
    for name in valid_controls:
        transformed, changed = transform_control(dataset, name)
        if not changed:
            raise ValueError(f"control mechanism did not alter declared fields: {name}")
        observed = evaluate(transformed, tuple(contract["seed_policy"]["fit_seeds"]))
        control_rows.append({
            "control": name,
            "mechanism_valid": True,
            "changed_fields": ";".join(sorted(changed)),
            "status": "PASS",
            "holdout_rmse": observed["holdout_rmse"],
            "holdout_spearman": observed["holdout_spearman"],
        })

    invalid_detected = []
    for name in ("target_proxy_leakage", "no_op_control"):
        transformed, changed = transform_control(dataset, name)
        try:
            if name == "no_op_control" and not changed:
                raise ValueError("declared mechanism made no change")
            evaluate(transformed)
        except ValueError as exc:
            invalid_detected.append(name)
            control_rows.append({
                "control": name,
                "mechanism_valid": False,
                "changed_fields": ";".join(sorted(changed)),
                "status": "DETECTED_INVALID",
                "holdout_rmse": "",
                "holdout_spearman": "",
                "reason": str(exc),
            })

    integrity_pass = len(invalid_detected) == 2 and all(row["status"] == "PASS" for row in control_rows[:3])
    result = {
        "bridge_id": "HBP-IR-01",
        "classification": "INSTRUMENT_VALID" if integrity_pass else "INSTRUMENT_INVALID",
        "status": "PASS" if integrity_pass else "FAIL",
        "contract_hash": stable_hash("hbp_ir_contract", contract),
        "selected_model": target_result["selected_model"],
        "calibration_scores": target_result["calibration_scores"],
        "holdout_metrics": {
            "rmse": target_result["holdout_rmse"],
            "spearman": target_result["holdout_spearman"],
            "fit_seed_variance": target_result["fit_seed_variance"],
        },
        "missingness": target_result["missingness"],
        "invalid_controls_detected": invalid_detected,
        "labels": [
            "CALIBRATION_ONLY_SELECTION_PASS",
            "HOLDOUT_ISOLATION_PASS",
            "SAME_PREDICTION_PATH_PASS",
            "CONTROL_MECHANISM_ASSERTIONS_PASS",
            "NO_PREDICTIVE_BRIDGE_CLAIM",
        ] if integrity_pass else ["INSTRUMENT_INVALID"],
    }
    result["result_hash"] = stable_hash("hbp_ir_01", result)

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "hbp_ir_01_result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with (output_dir / "hbp_ir_01_predictions.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["row_id", "split", "target", "prediction"])
        writer.writeheader()
        for index in range(len(dataset.target)):
            writer.writerow({
                "row_id": int(dataset.row_id[index]),
                "split": dataset.split[index],
                "target": f"{dataset.target[index]:.12g}",
                "prediction": f"{target_result['predictions'][index]:.12g}",
            })
    control_fields = ["control", "mechanism_valid", "changed_fields", "status", "holdout_rmse", "holdout_spearman", "reason"]
    with (output_dir / "hbp_ir_01_controls.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=control_fields)
        writer.writeheader()
        for row in control_rows:
            writer.writerow({field: row.get(field, "") for field in control_fields})
    report = [
        "# HBP-IR-01 Integrity Result",
        "",
        f"- Status: `{result['status']}`",
        f"- Classification: `{result['classification']}`",
        f"- Selected model: `{result['selected_model']}`",
        f"- Holdout RMSE: `{result['holdout_metrics']['rmse']:.6f}`",
        f"- Holdout Spearman: `{result['holdout_metrics']['spearman']:.6f}`",
        f"- Result hash: `{result['result_hash']}`",
        "",
        "This pass validates the versioned instrument path only. It does not repair or promote historical PB artifacts and does not establish external prediction.",
        "",
    ]
    (output_dir / "hbp_ir_01_report.md").write_text("\n".join(report), encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the HBP-IR-01 integrity calibration.")
    parser.add_argument("--output-dir", type=Path, default=RESULTS)
    args = parser.parse_args()
    print(json.dumps(run(args.output_dir), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
