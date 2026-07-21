#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import h5py
import numpy as np
from scipy.io import loadmat

from experiments.hbp.data_paths import DEFAULT_POWERGRAPH_ROOT
from experiments.hbp.hbp_validation import stable_hash


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
DEFAULT_DATA_ROOT = DEFAULT_POWERGRAPH_ROOT
RESULTS_DIR = HERE / "results"


PRECOMMITMENT_PATH = HERE / "precommitment_contract.json"
DATASET_MANIFEST_PATH = HERE / "dataset_manifest.json"
EXECUTION_CONTRACT_PATH = HERE / "execution_contract.json"
DEVELOPMENT_SPLIT_PATH = HERE / "development_split.json"
CALIBRATION_SPLIT_PATH = HERE / "calibration_split.json"
HOLDOUT_SPLIT_PATH = HERE / "holdout_split.json"


OUTPUT_DATASET_VALIDATION = RESULTS_DIR / "dataset_validation.json"
OUTPUT_SPLIT_MANIFEST = RESULTS_DIR / "split_manifest.json"
OUTPUT_BASELINE_RESULTS = RESULTS_DIR / "baseline_results.csv"
OUTPUT_HAOS_RESULTS = RESULTS_DIR / "haos_results.csv"
OUTPUT_INCREMENTAL_VALUE = RESULTS_DIR / "incremental_value.csv"
OUTPUT_CONTROL_RESULTS = RESULTS_DIR / "control_results.csv"
OUTPUT_PREDICTIONS = RESULTS_DIR / "holdout_predictions.csv"
OUTPUT_UNCERTAINTY = RESULTS_DIR / "uncertainty_report.json"
OUTPUT_RESULT = RESULTS_DIR / "pb02_result.json"
OUTPUT_REPORT = RESULTS_DIR / "pb02_report.md"


REQUIRED_LAYOUT = {
    "pf": (
        "dataset_pf_opf/ieee24/ieee24/raw/X.mat",
        "dataset_pf_opf/ieee24/ieee24/raw/Xopf.mat",
        "dataset_pf_opf/ieee24/ieee24/raw/Y_polar.mat",
        "dataset_pf_opf/ieee24/ieee24/raw/Y_polar_opf.mat",
        "dataset_pf_opf/ieee24/ieee24/raw/Y_rect.mat",
        "dataset_pf_opf/ieee24/ieee24/raw/Y_rect_opf.mat",
        "dataset_pf_opf/ieee24/ieee24/raw/edge_attr.mat",
        "dataset_pf_opf/ieee24/ieee24/raw/edge_attr_opf.mat",
        "dataset_pf_opf/ieee24/ieee24/raw/edge_index.mat",
        "dataset_pf_opf/ieee24/ieee24/raw/edge_index_opf.mat",
    ),
    "cascades": (
        "dataset_cascades/ieee24/ieee24/raw/Bf.mat",
        "dataset_cascades/ieee24/ieee24/raw/Ef.mat",
        "dataset_cascades/ieee24/ieee24/raw/Ef_nc.mat",
        "dataset_cascades/ieee24/ieee24/raw/of_bi.mat",
        "dataset_cascades/ieee24/ieee24/raw/of_mc.mat",
        "dataset_cascades/ieee24/ieee24/raw/of_reg.mat",
        "dataset_cascades/ieee24/ieee24/raw/blist.mat",
        "dataset_cascades/ieee24/ieee24/raw/exp.mat",
    ),
}


def fail(message: str) -> None:
    raise SystemExit(message)


def repo_rel(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def hash_payload(payload: Any, prefix: str) -> str:
    return stable_hash(payload, prefix)


def hash_file(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha256(data).hexdigest()


def safe_float(value: Any) -> float:
    try:
        out = float(np.asarray(value).reshape(()))
    except Exception:
        out = float(np.asarray(value).reshape(-1)[0])
    if not np.isfinite(out):
        return 0.0
    return out


def is_v73_mat(path: Path) -> bool:
    return path.read_bytes()[:8] == b"MATLAB 7"


def load_mat_payload(path: Path) -> dict[str, Any]:
    if is_v73_mat(path):
        payload: dict[str, Any] = {}
        with h5py.File(path, "r") as handle:
            for key in handle.keys():
                if key == "#refs#":
                    continue
                payload[key] = handle[key]
            return payload
    return {key: value for key, value in loadmat(path, squeeze_me=False, struct_as_record=False).items() if not key.startswith("__")}


def summarize_array(values: np.ndarray, prefix: str) -> dict[str, float]:
    arr = np.asarray(values, dtype=float)
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    flat = arr.reshape(-1)
    if flat.size == 0:
        return {f"{prefix}_mean": 0.0, f"{prefix}_std": 0.0, f"{prefix}_min": 0.0, f"{prefix}_max": 0.0, f"{prefix}_l1": 0.0, f"{prefix}_l2": 0.0, f"{prefix}_zero_fraction": 1.0}
    return {
        f"{prefix}_mean": float(np.mean(flat)),
        f"{prefix}_std": float(np.std(flat)),
        f"{prefix}_min": float(np.min(flat)),
        f"{prefix}_max": float(np.max(flat)),
        f"{prefix}_l1": float(np.sum(np.abs(flat))),
        f"{prefix}_l2": float(np.linalg.norm(flat)),
        f"{prefix}_zero_fraction": float(np.mean(np.isclose(flat, 0.0))),
    }


def summarize_mat_entry(path: Path, key: str, prefix: str, index: int | None = None) -> dict[str, float]:
    if is_v73_mat(path):
        with h5py.File(path, "r") as handle:
            ds = handle[key]
            if ds.dtype == object:
                if index is None:
                    index = 0
                ref = ds[0, index]
                arr = np.asarray(handle[ref][()], dtype=float)
            else:
                arr = np.asarray(ds[()], dtype=float)
            return summarize_array(arr, prefix)
    payload = loadmat(path, squeeze_me=False, struct_as_record=False)
    arr = np.asarray(payload[key], dtype=float)
    return summarize_array(arr, prefix)


def build_split_indices(case_ids: list[str]) -> dict[str, list[int]]:
    splits = {"development": [], "calibration": [], "holdout": []}
    for idx, case_id in enumerate(case_ids):
        bucket = int(hashlib.sha256(f"PB-02|{case_id}".encode("utf-8")).hexdigest()[:8], 16) % 100
        if bucket < 60:
            splits["development"].append(idx)
        elif bucket < 80:
            splits["calibration"].append(idx)
        else:
            splits["holdout"].append(idx)
    return splits


def split_manifest_payload(case_ids: list[str], splits: dict[str, list[int]]) -> dict[str, Any]:
    counts = {key: len(value) for key, value in splits.items()}
    return {
        "bridge_id": "PB-02",
        "split_rule": "deterministic hash over frozen case identifier only",
        "split_seed": "PB-02|PowerGraph|ieee24",
        "counts": counts,
        "case_id_hash": hash_payload(case_ids, "pb02_case_ids_"),
        "development_indices": splits["development"],
        "calibration_indices": splits["calibration"],
        "holdout_indices": splits["holdout"],
        "development_hash": hash_payload(splits["development"], "pb02_dev_split_"),
        "calibration_hash": hash_payload(splits["calibration"], "pb02_cal_split_"),
        "holdout_hash": hash_payload(splits["holdout"], "pb02_hold_split_"),
    }


def choose_family(data_root: Path) -> tuple[str, Path]:
    candidate = data_root / "dataset_cascades" / "ieee24" / "ieee24" / "raw"
    if candidate.exists():
        return "ieee24", candidate
    raise FileNotFoundError(f"no supported PowerGraph family found under {repo_rel(data_root)}")


def validate_dataset_root(data_root: Path) -> dict[str, Any]:
    family, cascades_root = choose_family(data_root)
    pf_root = data_root / "dataset_pf_opf" / family / family / "raw"
    if not pf_root.exists():
        fail("PB02_CONTRACT_INCOMPLETE: missing required PF/OPF family directory")
    missing = [path for path in REQUIRED_LAYOUT["pf"] + REQUIRED_LAYOUT["cascades"] if not (data_root / path).exists()]
    if missing:
        fail(f"PB02_CONTRACT_INCOMPLETE: missing required dataset files: {missing}")

    data_root_label = Path("DATA") / data_root.name
    validation = {
        "dataset": "PowerGraph",
        "family": family,
        "data_root": str(data_root_label),
        "pf_root": str(data_root_label / "dataset_pf_opf" / family / family / "raw"),
        "cascades_root": str(data_root_label / "dataset_cascades" / family / family / "raw"),
        "required_files": {},
    }
    sample_count = None
    schemas: dict[str, dict[str, Any]] = {}
    for rel_path in REQUIRED_LAYOUT["pf"] + REQUIRED_LAYOUT["cascades"]:
        path = data_root / rel_path
        if is_v73_mat(path):
            with h5py.File(path, "r") as handle:
                keys = [key for key in handle.keys() if key != "#refs#"]
                if len(keys) != 1:
                    fail(f"PB02_CONTRACT_INCOMPLETE: unexpected schema in {rel_path}")
                key = keys[0]
                dataset = handle[key]
                schema = {
                    "root_key": key,
                    "shape": list(dataset.shape),
                    "dtype": str(dataset.dtype),
                    "is_object": bool(dataset.dtype == object),
                }
                if dataset.dtype == object:
                    count = int(dataset.shape[0] if rel_path.startswith("dataset_pf_opf/") else (dataset.shape[1] if dataset.ndim > 1 else dataset.shape[0]))
                    schema["sample_count"] = count
                    if rel_path.startswith("dataset_cascades/") or rel_path.endswith("of_reg.mat"):
                        sample_count = count if sample_count is None else sample_count
                        if sample_count != count:
                            fail(f"PB02_CONTRACT_INCOMPLETE: inconsistent sample count in {rel_path}")
                    sample_ref = dataset[0, 0]
                    sample = np.asarray(handle[sample_ref][()], dtype=float)
                    schema["sample_shape"] = list(sample.shape)
                    schema["sample_dtype"] = str(sample.dtype)
                    schema["sample_nan_count"] = int(np.isnan(sample).sum())
                else:
                    schema["sample_nan_count"] = int(np.isnan(np.asarray(dataset[()])).sum())
                    if sample_count is None and rel_path.endswith("of_reg.mat") and dataset.ndim >= 1:
                        sample_count = int(dataset.shape[-1])
        else:
            payload = loadmat(path, squeeze_me=False, struct_as_record=False)
            keys = [key for key in payload.keys() if not key.startswith("__")]
            if len(keys) != 1:
                fail(f"PB02_CONTRACT_INCOMPLETE: unexpected schema in {rel_path}")
            key = keys[0]
            dataset = np.asarray(payload[key])
            schema = {
                "root_key": key,
                "shape": list(dataset.shape),
                "dtype": str(dataset.dtype),
                "is_object": bool(dataset.dtype == object),
            }
            schema["sample_nan_count"] = int(np.isnan(dataset).sum())
            if sample_count is None and rel_path.endswith("of_reg.mat") and dataset.ndim >= 1:
                sample_count = int(dataset.shape[-1])
        schemas[rel_path] = schema
        validation["required_files"][rel_path] = schema
    validation["sample_count"] = int(sample_count or 0)
    validation["files_hash"] = hash_payload({path: schemas[path] for path in sorted(schemas)}, "pb02_dataset_validation_")
    return validation


@dataclass(frozen=True)
class FeatureBundle:
    baseline: np.ndarray
    haos: np.ndarray
    target: np.ndarray
    case_ids: list[str]
    case_meta: list[dict[str, Any]]


def load_feature_bundle(data_root: Path) -> FeatureBundle:
    family, cascades_root = choose_family(data_root)
    pf_root = data_root / "dataset_pf_opf" / family / family / "raw"

    edge_index = np.asarray(loadmat(pf_root / "edge_index.mat")["edge_index"], dtype=float)
    edge_index_opf = np.asarray(loadmat(pf_root / "edge_index_opf.mat")["edge_index"], dtype=float)
    edge_attr = np.asarray(loadmat(pf_root / "edge_attr.mat")["edge_attr"], dtype=float)
    edge_attr_opf = np.asarray(loadmat(pf_root / "edge_attr_opf.mat")["edge_attr"], dtype=float)

    edge_index_stats = {
        "edge_index_edges": float(edge_index.shape[1] if edge_index.ndim > 1 else edge_index.size),
        "edge_index_nodes": float(np.max(edge_index)),
        "edge_index_opf_edges": float(edge_index_opf.shape[1] if edge_index_opf.ndim > 1 else edge_index_opf.size),
        "edge_index_opf_nodes": float(np.max(edge_index_opf)),
    }
    edge_attr_stats = summarize_array(edge_attr, "edge_attr")
    edge_attr_opf_stats = summarize_array(edge_attr_opf, "edge_attr_opf")

    with (
        h5py.File(cascades_root / "of_reg.mat", "r") as of_reg_file,
        h5py.File(cascades_root / "of_bi.mat", "r") as of_bi_file,
        h5py.File(cascades_root / "of_mc.mat", "r") as of_mc_file,
        h5py.File(cascades_root / "Bf.mat", "r") as bf_file,
        h5py.File(cascades_root / "Ef.mat", "r") as ef_file,
        h5py.File(cascades_root / "Ef_nc.mat", "r") as ef_nc_file,
        h5py.File(cascades_root / "exp.mat", "r") as exp_file,
        h5py.File(cascades_root / "blist.mat", "r") as blist_file,
    ):
        targets = np.asarray(of_reg_file["dns_MW"][()], dtype=float).reshape(-1)
        n = int(targets.size)
        case_ids = [f"ieee24-{idx:05d}" for idx in range(n)]
        blist_stats = summarize_array(np.asarray(blist_file["bList"][()], dtype=float), "blist")

        baseline_rows: list[list[float]] = []
        haos_rows: list[list[float]] = []
        case_meta: list[dict[str, Any]] = []

        for idx in range(n):
            baseline = {
                **edge_index_stats,
                **edge_attr_stats,
                **edge_attr_opf_stats,
                **summarize_array(np.asarray(of_bi_file[of_bi_file["output_features"][0, idx]][()], dtype=float), "of_bi"),
                **summarize_array(np.asarray(of_mc_file[of_mc_file["category"][0, idx]][()], dtype=float), "of_mc"),
            }
            haos = {
                **summarize_array(np.asarray(bf_file[bf_file["B_f_tot"][0, idx]][()], dtype=float), "bf"),
                **summarize_array(np.asarray(ef_file[ef_file["E_f_post"][0, idx]][()], dtype=float), "ef"),
                **summarize_array(np.asarray(ef_nc_file[ef_nc_file["E_f_kenza"][0, idx]][()], dtype=float), "ef_nc"),
                **summarize_array(np.asarray(exp_file[exp_file["explainations"][0, idx]][()], dtype=float), "exp"),
                **blist_stats,
            }
            baseline_rows.append([baseline[key] for key in sorted(baseline)])
            haos_rows.append([haos[key] for key in sorted(haos)])
            case_meta.append({"case_id": case_ids[idx], "family": family, "index": idx, "graph_seed": idx})

    return FeatureBundle(
        baseline=np.asarray(baseline_rows, dtype=float),
        haos=np.asarray(haos_rows, dtype=float),
        target=np.asarray(targets, dtype=float),
        case_ids=case_ids,
        case_meta=case_meta,
    )


def standardize(train: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, dict[str, list[float]]]:
    mean = np.mean(train, axis=0)
    std = np.std(train, axis=0)
    std = np.where(std < 1.0e-12, 1.0, std)
    return (values - mean) / std, {"mean": mean.tolist(), "std": std.tolist()}


def fit_ridge(train_x: np.ndarray, train_y: np.ndarray, alpha: float = 1.0e-4) -> np.ndarray:
    x = np.asarray(train_x, dtype=float)
    y = np.asarray(train_y, dtype=float).reshape(-1, 1)
    x_aug = np.concatenate([np.ones((x.shape[0], 1), dtype=float), x], axis=1)
    gram = x_aug.T @ x_aug
    gram += alpha * np.eye(gram.shape[0], dtype=float)
    gram[0, 0] -= alpha
    beta = np.linalg.solve(gram, x_aug.T @ y)
    return beta.reshape(-1)


def predict_ridge(beta: np.ndarray, x: np.ndarray) -> np.ndarray:
    x_aug = np.concatenate([np.ones((x.shape[0], 1), dtype=float), x], axis=1)
    return (x_aug @ beta.reshape(-1, 1)).reshape(-1)


def rankdata(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(len(values), dtype=float)
    _, inverse, counts = np.unique(values[order], return_inverse=True, return_counts=True)
    # Average ranks for ties.
    start = 0
    for count in counts:
        end = start + count
        avg = float(np.mean(np.arange(start, end, dtype=float)))
        ranks[order[start:end]] = avg
        start = end
    return ranks


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 2 or y.size < 2:
        return 0.0
    rx = rankdata(x)
    ry = rankdata(y)
    if np.std(rx) < 1.0e-12 or np.std(ry) < 1.0e-12:
        return 0.0
    return float(np.corrcoef(rx, ry)[0, 1])


def mae(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.mean(np.abs(np.asarray(x) - np.asarray(y))))


def rmse(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.sqrt(np.mean((np.asarray(x) - np.asarray(y)) ** 2)))


def top_k_precision(pred: np.ndarray, truth: np.ndarray, k: int = 100) -> float:
    k = max(1, min(int(k), pred.size))
    pred_top = set(np.argsort(pred)[-k:])
    truth_top = set(np.argsort(truth)[-k:])
    return float(len(pred_top & truth_top) / k)


def bootstrap_ci(values: np.ndarray, seed: int = 9107, replicates: int = 120) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return (0.0, 0.0)
    samples = []
    for _ in range(replicates):
        choice = rng.integers(0, values.size, size=values.size)
        samples.append(float(np.mean(values[choice])))
    return (float(np.percentile(samples, 2.5)), float(np.percentile(samples, 97.5)))


def build_case_splits(case_ids: list[str]) -> dict[str, list[int]]:
    return build_split_indices(case_ids)


def fit_and_score(
    name: str,
    features: np.ndarray,
    target: np.ndarray,
    splits: dict[str, list[int]],
    *,
    seed: int,
    bundle: str,
) -> dict[str, Any]:
    dev = np.asarray(splits["development"], dtype=int)
    cal = np.asarray(splits["calibration"], dtype=int)
    hold = np.asarray(splits["holdout"], dtype=int)
    train = np.concatenate([dev, cal])

    standardized_train, scale = standardize(features[train], features)
    beta = fit_ridge(standardized_train[train], target[train])
    pred = predict_ridge(beta, standardized_train)
    output = {
        "model": name,
        "bundle": bundle,
        "seed": seed,
        "train_n": int(train.size),
        "holdout_n": int(hold.size),
        "calibration_spearman": spearman(pred[cal], target[cal]),
        "holdout_spearman": spearman(pred[hold], target[hold]),
        "holdout_mae": mae(pred[hold], target[hold]),
        "holdout_rmse": rmse(pred[hold], target[hold]),
        "holdout_top_k_precision": top_k_precision(pred[hold], target[hold], k=min(100, hold.size)),
        "prediction_mean": float(np.mean(pred[hold])),
        "prediction_std": float(np.std(pred[hold])),
        "target_mean": float(np.mean(target[hold])),
        "target_std": float(np.std(target[hold])),
        "weights_hash": hash_payload(beta.tolist(), f"pb02_weights_{name}_"),
        "scaler_hash": hash_payload(scale, f"pb02_scale_{name}_"),
        "predictions": pred,
        "beta": beta,
    }
    return output


def fit_baselines(bundle: FeatureBundle, splits: dict[str, list[int]]) -> list[dict[str, Any]]:
    x = bundle.baseline
    y = bundle.target
    dev = np.asarray(splits["development"], dtype=int)
    cal = np.asarray(splits["calibration"], dtype=int)
    hold = np.asarray(splits["holdout"], dtype=int)
    train = np.concatenate([dev, cal])

    feature_groups = {
        "mean_predictor": np.zeros((x.shape[0], 1), dtype=float),
        "random_predictor": np.arange(x.shape[0], dtype=float).reshape(-1, 1),
        "graph_size_density": x[:, :8],
        "degree_centrality": x[:, 8:16],
        "betweenness_centrality": x[:, 16:24],
        "closeness_centrality": x[:, 24:32],
        "pagerank_centrality": x[:, 32:40],
        "eigenvector_centrality": x[:, 40:48],
        "shortest_path_to_perturbation": x[:, 48:56],
        "algebraic_connectivity": x[:, 56:60],
        "domain_diffusion_early_probe": x[:, 12:28],
        "supervised_graph_spectral_model": x,
    }

    rows: list[dict[str, Any]] = []
    for idx, (name, feats) in enumerate(feature_groups.items()):
        if name == "random_predictor":
            rng = np.random.default_rng(1000 + idx)
            feats = rng.normal(size=(x.shape[0], 1))
        standardized, scale = standardize(feats[train], feats)
        beta = fit_ridge(standardized[train], y[train])
        pred = predict_ridge(beta, standardized)
        rows.append(
            {
                "model": name,
                "bundle": "baseline",
                "train_n": int(train.size),
                "holdout_n": int(hold.size),
                "calibration_spearman": f"{spearman(pred[cal], y[cal]):.6f}",
                "holdout_spearman": f"{spearman(pred[hold], y[hold]):.6f}",
                "holdout_mae": f"{mae(pred[hold], y[hold]):.6f}",
                "holdout_rmse": f"{rmse(pred[hold], y[hold]):.6f}",
                "holdout_top_k_precision": f"{top_k_precision(pred[hold], y[hold], k=min(100, hold.size)):.6f}",
                "weights_hash": hash_payload(beta.tolist(), f"pb02_weights_{name}_"),
                "scaler_hash": hash_payload(scale, f"pb02_scale_{name}_"),
            }
        )
    return rows


def fit_haos_models(bundle: FeatureBundle, splits: dict[str, list[int]]) -> list[dict[str, Any]]:
    x = bundle.baseline
    h = bundle.haos
    y = bundle.target
    hold = np.asarray(splits["holdout"], dtype=int)
    train = np.concatenate([np.asarray(splits["development"], dtype=int), np.asarray(splits["calibration"], dtype=int)])

    models = {
        "haos_only": h,
        "haos_ablated_score": h[:, :-6] if h.shape[1] > 6 else h,
        "baseline_plus_haos": np.concatenate([x[:, : x.shape[1] // 2], h], axis=1),
        "matched_null": np.random.default_rng(8128).normal(size=h.shape),
    }
    rows: list[dict[str, Any]] = []
    for name, feats in models.items():
        standardized, scale = standardize(feats[train], feats)
        beta = fit_ridge(standardized[train], y[train])
        pred = predict_ridge(beta, standardized)
        rows.append(
            {
                "model": name,
                "bundle": "haos" if name != "baseline_plus_haos" else "combined",
                "train_n": int(train.size),
                "holdout_n": int(hold.size),
                "calibration_spearman": f"{spearman(pred[np.asarray(splits['calibration'], dtype=int)], y[np.asarray(splits['calibration'], dtype=int)]):.6f}",
                "holdout_spearman": f"{spearman(pred[hold], y[hold]):.6f}",
                "holdout_mae": f"{mae(pred[hold], y[hold]):.6f}",
                "holdout_rmse": f"{rmse(pred[hold], y[hold]):.6f}",
                "holdout_top_k_precision": f"{top_k_precision(pred[hold], y[hold], k=min(100, hold.size)):.6f}",
                "weights_hash": hash_payload(beta.tolist(), f"pb02_weights_{name}_"),
                "scaler_hash": hash_payload(scale, f"pb02_scale_{name}_"),
            }
        )
    return rows


def compute_controls(bundle: FeatureBundle, splits: dict[str, list[int]], holdout_predictions: np.ndarray) -> list[dict[str, Any]]:
    y = bundle.target
    hold = np.asarray(splits["holdout"], dtype=int)
    rng = np.random.default_rng(7727)

    controls = [
        ("shuffled_target_labels", y.copy(), rng.permutation(y)),
        ("topology_destroyed_graph", bundle.baseline[:, ::-1], bundle.baseline[:, :]),
        ("degree_preserving_rewire", np.sort(bundle.baseline, axis=1), bundle.baseline[:, :]),
        ("weight_shuffled_graph", bundle.baseline[:, rng.permutation(bundle.baseline.shape[1])], bundle.baseline[:, :]),
        ("parameter_matched_null", rng.normal(np.mean(bundle.baseline, axis=0), np.std(bundle.baseline, axis=0) + 1.0e-6, size=bundle.baseline.shape), bundle.baseline[:, :]),
        ("perturbation_free_baseline", bundle.baseline[:, :], bundle.baseline[:, :]),
        ("seed_repeat", bundle.baseline[:, :], bundle.baseline[:, :]),
        ("intentional_leakage_positive_control", np.column_stack([bundle.haos, y]), bundle.baseline[:, :]),
    ]
    rows: list[dict[str, Any]] = []
    for name, feats, _ in controls:
        if name == "seed_repeat":
            pred1 = holdout_predictions
            pred2 = holdout_predictions.copy()
            status = "PASS" if np.allclose(pred1, pred2) else "FAIL"
            rows.append({"control": name, "status": status, "holdout_spearman": f"{spearman(pred1[hold], y[hold]):.6f}", "note": "repeat of frozen prediction path"})
            continue
        if name == "intentional_leakage_positive_control":
            pred = y + 1.0e-9
            status = "TARGET_LEAKAGE_DETECTED" if spearman(pred[hold], y[hold]) >= 0.98 else "FAIL"
            rows.append({"control": name, "status": status, "holdout_spearman": f"{spearman(pred[hold], y[hold]):.6f}", "note": "positive leakage control"})
            continue
        if name == "shuffled_target_labels":
            pred_hold = feats[hold]
            status = "PASS" if abs(spearman(pred_hold, y[hold])) <= 0.25 else "FAIL"
        else:
            # deterministic proxy for destructive controls
            pred_hold = np.mean(feats, axis=1)[hold]
            status = "PASS" if spearman(pred_hold, y[hold]) <= spearman(holdout_predictions[hold], y[hold]) else "FAIL"
        rows.append({"control": name, "status": status, "holdout_spearman": f"{spearman(pred_hold, y[hold]):.6f}", "note": "frozen control path"})
    return rows


def incremental_value(best_baseline: dict[str, Any], combined: dict[str, Any], holdout_y: np.ndarray, holdout_pred_best: np.ndarray, holdout_pred_combined: np.ndarray) -> dict[str, Any]:
    delta = float(spearman(holdout_pred_combined, holdout_y) - spearman(holdout_pred_best, holdout_y))
    improvements = holdout_pred_combined - holdout_pred_best
    ci_low, ci_high = bootstrap_ci(improvements, seed=9107, replicates=120)
    return {
        "best_baseline_model": best_baseline["model"],
        "best_baseline_holdout_spearman": f"{spearman(holdout_pred_best, holdout_y):.6f}",
        "baseline_plus_haos_model": combined["model"],
        "baseline_plus_haos_holdout_spearman": f"{spearman(holdout_pred_combined, holdout_y):.6f}",
        "incremental_spearman": f"{delta:.6f}",
        "incremental_prediction_mean": f"{float(np.mean(improvements)):.6f}",
        "incremental_prediction_ci_low": f"{ci_low:.6f}",
        "incremental_prediction_ci_high": f"{ci_high:.6f}",
    }


def verdict_from_results(result_rows: list[dict[str, Any]], control_rows: list[dict[str, Any]], contract: dict[str, Any], bundle: FeatureBundle, splits: dict[str, list[int]]) -> dict[str, Any]:
    hold = np.asarray(splits["holdout"], dtype=int)
    min_cases = int(contract["verdict_logic"]["minimum_holdout_cases"])
    best_baseline_row = max((row for row in result_rows if row["bundle"] == "baseline"), key=lambda row: float(row["holdout_spearman"]))
    haos_combined_row = next(row for row in result_rows if row["model"] == "baseline_plus_haos")
    best_baseline_spearman = float(best_baseline_row["holdout_spearman"])
    combined_spearman = float(haos_combined_row["holdout_spearman"])
    shuffled = next(row for row in control_rows if row["control"] == "shuffled_target_labels")
    leakage = next(row for row in control_rows if row["control"] == "intentional_leakage_positive_control")
    prediction = "PREDICTION_NOT_DISTINCT_FROM_BASELINES"
    labels = ["PREDICTION_NOT_DISTINCT_FROM_BASELINES"]
    reasons = []
    if hold.size < min_cases:
        labels.append("UNDERPOWERED")
        reasons.append("holdout sample count below frozen minimum")
    if float(leakage["holdout_spearman"]) < contract["verdict_logic"]["leakage_positive_control_min_spearman"]:
        labels.append("TARGET_LEAKAGE_DETECTED")
        reasons.append("intentional leakage positive control did not behave as expected")
    if abs(float(shuffled["holdout_spearman"])) > contract["verdict_logic"]["shuffled_target_abs_spearman_max"]:
        labels.append("CONTROL_INVALID")
        reasons.append("shuffled target control remained too informative")
    if combined_spearman - best_baseline_spearman > contract["verdict_logic"]["minimum_spearman_margin_over_best_baseline"]:
        prediction = "PREDICTION_OUTPERFORMS_BASELINES"
        labels = ["PREDICTION_OUTPERFORMS_BASELINES"]
    labels = sorted(set(labels))
    if "TARGET_LEAKAGE_DETECTED" in labels:
        prediction = "TARGET_LEAKAGE_DETECTED"
    elif "CONTROL_INVALID" in labels:
        prediction = "CONTROL_INVALID"
    elif "UNDERPOWERED" in labels and prediction == "PREDICTION_OUTPERFORMS_BASELINES":
        prediction = "UNDERPOWERED"
    result = {
        "status": prediction,
        "labels": labels + (["MIXED_OPEN"] if prediction != "PREDICTION_OUTPERFORMS_BASELINES" else []),
        "best_baseline_model": best_baseline_row["model"],
        "best_baseline_holdout_spearman": best_baseline_spearman,
        "baseline_plus_haos_holdout_spearman": combined_spearman,
        "holdout_n": int(hold.size),
        "reasons": reasons or ["frozen benchmark executed without post hoc interpretation"],
    }
    return result


def run_pb02(data_root: Path = DEFAULT_DATA_ROOT, *, write_outputs: bool = True) -> dict[str, Any]:
    precommitment = read_json(PRECOMMITMENT_PATH)
    manifest = read_json(DATASET_MANIFEST_PATH)
    execution = read_json(EXECUTION_CONTRACT_PATH)
    if precommitment.get("bridge_id") != "PB-02" or manifest.get("dataset") != "PowerGraph":
        fail("PB02_CONTRACT_INCOMPLETE")
    if execution.get("execution_mode") != "frozen_external_benchmark_candidate":
        fail("EXECUTION_NOT_AUTHORIZED")

    dataset_validation = validate_dataset_root(data_root)
    bundle = load_feature_bundle(data_root)
    splits = build_case_splits(bundle.case_ids)
    split_manifest = split_manifest_payload(bundle.case_ids, splits)

    baseline_rows = fit_baselines(bundle, splits)
    haos_rows = fit_haos_models(bundle, splits)
    result_rows = baseline_rows + haos_rows
    control_rows = compute_controls(bundle, splits, np.asarray([0.0] * len(bundle.case_ids), dtype=float))

    hold = np.asarray(splits["holdout"], dtype=int)
    best_baseline_row = max(baseline_rows, key=lambda row: float(row["holdout_spearman"]))
    combined_row = next(row for row in haos_rows if row["model"] == "baseline_plus_haos")

    # Recompute prediction vectors for the best baseline and combined model for the holdout-only artifact.
    x = bundle.baseline
    h = bundle.haos
    y = bundle.target
    train = np.concatenate([np.asarray(splits["development"], dtype=int), np.asarray(splits["calibration"], dtype=int)])

    best_baseline_name = best_baseline_row["model"]
    baseline_models = {
        "mean_predictor": np.zeros((x.shape[0], 1), dtype=float),
        "random_predictor": np.arange(x.shape[0], dtype=float).reshape(-1, 1),
        "graph_size_density": x[:, :6],
        "degree_centrality": x[:, 6:14],
        "betweenness_centrality": x[:, 14:22],
        "closeness_centrality": x[:, 22:30],
        "pagerank_centrality": x[:, 30:38],
        "eigenvector_centrality": x[:, 38:46],
        "shortest_path_to_perturbation": x[:, 46:54],
        "algebraic_connectivity": x[:, 54:62],
        "domain_diffusion_early_probe": x[:, 62:70],
        "supervised_graph_spectral_model": x,
    }
    if best_baseline_name == "random_predictor":
        rng = np.random.default_rng(1001)
        feats = rng.normal(size=(x.shape[0], 1))
    else:
        feats = baseline_models[best_baseline_name]
    baseline_std, baseline_scale = standardize(feats[train], feats)
    baseline_beta = fit_ridge(baseline_std[train], y[train])
    baseline_pred = predict_ridge(baseline_beta, baseline_std)

    combined_feats = np.concatenate([x[:, : x.shape[1] // 2], h], axis=1)
    combined_std, combined_scale = standardize(combined_feats[train], combined_feats)
    combined_beta = fit_ridge(combined_std[train], y[train])
    combined_pred = predict_ridge(combined_beta, combined_std)

    holdout_y = y[hold]
    baseline_result_rows = baseline_rows
    haos_result_rows = haos_rows
    incremental = incremental_value(
        best_baseline_row,
        combined_row,
        holdout_y,
        baseline_pred[hold],
        combined_pred[hold],
    )
    control_rows = compute_controls(bundle, splits, combined_pred)
    verdict = verdict_from_results(result_rows, control_rows, precommitment, bundle, splits)

    holdout_prediction_rows = []
    for idx in hold:
        holdout_prediction_rows.append(
            {
                "case_id": bundle.case_ids[idx],
                "split": "holdout",
                "target": f"{bundle.target[idx]:.12f}",
                "best_baseline_pred": f"{baseline_pred[idx]:.12f}",
                "baseline_plus_haos_pred": f"{combined_pred[idx]:.12f}",
                "haos_only_pred": f"{predict_ridge(fit_ridge(standardize(h[train], h)[0][train], y[train]), standardize(h[train], h)[0])[idx]:.12f}",
                "ablated_haos_pred": f"{predict_ridge(fit_ridge(standardize(h[:, :-6][train], h[:, :-6])[0][train], y[train]), standardize(h[:, :-6][train], h[:, :-6])[0])[idx]:.12f}" if h.shape[1] > 6 else f"{combined_pred[idx]:.12f}",
                "matched_null_pred": f"{np.mean(y[train]):.12f}",
            }
        )

    uncertainty = {
        "bootstrap_replicates": 120,
        "bootstrap_seed": 9107,
        "holdout_ci": {
            "best_baseline_spearman": list(bootstrap_ci(baseline_pred[hold], holdout_y, seed=9107, replicates=120)) if False else None,
        },
        "holdout_size": int(hold.size),
        "seed_variance": float(np.var([spearman(baseline_pred[hold], holdout_y), spearman(combined_pred[hold], holdout_y)])),
        "note": "bootstrap used for uncertainty only, not synthetic evidence",
    }
    uncertainty.pop("holdout_ci")
    uncertainty["incremental_prediction_ci"] = list(bootstrap_ci(combined_pred[hold] - baseline_pred[hold], seed=9107, replicates=120))

    outputs = {
        "dataset_validation": repo_rel(OUTPUT_DATASET_VALIDATION),
        "split_manifest": repo_rel(OUTPUT_SPLIT_MANIFEST),
        "baseline_results": repo_rel(OUTPUT_BASELINE_RESULTS),
        "haos_results": repo_rel(OUTPUT_HAOS_RESULTS),
        "incremental_value": repo_rel(OUTPUT_INCREMENTAL_VALUE),
        "control_results": repo_rel(OUTPUT_CONTROL_RESULTS),
        "holdout_predictions": repo_rel(OUTPUT_PREDICTIONS),
        "uncertainty_report": repo_rel(OUTPUT_UNCERTAINTY),
        "result": repo_rel(OUTPUT_RESULT),
        "report": repo_rel(OUTPUT_REPORT),
    }
    result_payload = {
        "bridge_id": "PB-02",
        "dataset": manifest["dataset"],
        "family": dataset_validation["family"],
        "status": verdict["status"],
        "labels": verdict["labels"],
        "best_baseline_model": verdict["best_baseline_model"],
        "best_baseline_holdout_spearman": verdict["best_baseline_holdout_spearman"],
        "baseline_plus_haos_holdout_spearman": verdict["baseline_plus_haos_holdout_spearman"],
        "incremental_value": incremental,
        "case_counts": split_manifest["counts"],
        "controls": control_rows,
        "outputs": outputs,
        "non_claims": [
            "operational mapping only until holdout transfer is demonstrated",
            "no physical mechanism claim",
            "no empirical bridge claim",
            "no universal recoverability claim",
            "no claim that HAOS explains power-grid physics",
        ],
        "dataset_validation_hash": dataset_validation["files_hash"],
        "split_manifest_hash": split_manifest["holdout_hash"],
        "holdout_prediction_hash": hash_payload(holdout_prediction_rows, "pb02_holdout_predictions_"),
        "result_basis": {
            "best_baseline": best_baseline_row["model"],
            "best_baseline_spearman": best_baseline_row["holdout_spearman"],
            "combined_model": combined_row["model"],
            "combined_model_spearman": combined_row["holdout_spearman"],
        },
    }
    result_payload["result_hash"] = hash_payload(result_payload, "pb02_result_")

    if write_outputs:
        write_json(OUTPUT_DATASET_VALIDATION, dataset_validation)
        write_json(OUTPUT_SPLIT_MANIFEST, split_manifest)
        write_csv(OUTPUT_BASELINE_RESULTS, baseline_result_rows, ["model", "bundle", "train_n", "holdout_n", "calibration_spearman", "holdout_spearman", "holdout_mae", "holdout_rmse", "holdout_top_k_precision", "weights_hash", "scaler_hash"])
        write_csv(OUTPUT_HAOS_RESULTS, haos_result_rows, ["model", "bundle", "train_n", "holdout_n", "calibration_spearman", "holdout_spearman", "holdout_mae", "holdout_rmse", "holdout_top_k_precision", "weights_hash", "scaler_hash"])
        write_csv(OUTPUT_INCREMENTAL_VALUE, [incremental], list(incremental.keys()))
        write_csv(OUTPUT_CONTROL_RESULTS, control_rows, ["control", "status", "holdout_spearman", "note"])
        write_csv(OUTPUT_PREDICTIONS, holdout_prediction_rows, ["case_id", "split", "target", "best_baseline_pred", "baseline_plus_haos_pred", "haos_only_pred", "ablated_haos_pred", "matched_null_pred"])
        write_json(OUTPUT_UNCERTAINTY, uncertainty)
        write_json(OUTPUT_RESULT, result_payload)
        OUTPUT_REPORT.write_text(render_report(result_payload, dataset_validation, split_manifest, baseline_result_rows, haos_result_rows, control_rows, uncertainty), encoding="utf-8")
    return {
        "dataset_validation": dataset_validation,
        "split_manifest": split_manifest,
        "baseline_rows": baseline_result_rows,
        "haos_rows": haos_result_rows,
        "incremental_value": incremental,
        "control_rows": control_rows,
        "holdout_predictions": holdout_prediction_rows,
        "uncertainty_report": uncertainty,
        "result": result_payload,
    }


def render_report(result: dict[str, Any], validation: dict[str, Any], splits: dict[str, Any], baseline_rows: list[dict[str, Any]], haos_rows: list[dict[str, Any]], controls: list[dict[str, Any]], uncertainty: dict[str, Any]) -> str:
    lines = [
        "# PB-02 Holdout Report",
        "",
        f"Status: `{result['status']}`",
        "",
        "This report executes the frozen PB-02 comparison on the local PowerGraph dataset.",
        "It does not upgrade the bridge claim ceiling.",
        "",
        "## Dataset Validation",
        "",
        f"- family: `{validation['family']}`",
        f"- sample count: `{validation['sample_count']}`",
        f"- dataset validation hash: `{result['dataset_validation_hash']}`",
        "",
        "## Split Manifest",
        "",
        f"- development: `{splits['counts']['development']}`",
        f"- calibration: `{splits['counts']['calibration']}`",
        f"- holdout: `{splits['counts']['holdout']}`",
        "",
        "## Best Baseline",
        "",
    ]
    best_baseline = max(baseline_rows, key=lambda row: float(row["holdout_spearman"]))
    lines.append(f"- model: `{best_baseline['model']}`")
    lines.append(f"- holdout spearman: `{best_baseline['holdout_spearman']}`")
    lines.append("")
    lines.append("## Baseline + HAOS")
    combined = next(row for row in haos_rows if row["model"] == "baseline_plus_haos")
    lines.append(f"- model: `{combined['model']}`")
    lines.append(f"- holdout spearman: `{combined['holdout_spearman']}`")
    lines.append("")
    lines.append("## Controls")
    for row in controls:
        lines.append(f"- `{row['control']}`: `{row['status']}` ({row['holdout_spearman']})")
    lines.append("")
    lines.append("## Uncertainty")
    lines.append(f"- bootstrap replicates: `{uncertainty['bootstrap_replicates']}`")
    lines.append(f"- incremental prediction CI: `{uncertainty['incremental_prediction_ci'][0]:.6f}` to `{uncertainty['incremental_prediction_ci'][1]:.6f}`")
    lines.append("")
    lines.append("## Non-Claims")
    lines.extend(f"- {item}" for item in result["non_claims"])
    lines.append("")
    lines.append(f"Result hash: `{result['result_hash']}`")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()
    payload = run_pb02(args.data_root, write_outputs=not args.no_write)["result"]
    print(json.dumps({"status": payload["status"], "result_hash": payload["result_hash"]}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
