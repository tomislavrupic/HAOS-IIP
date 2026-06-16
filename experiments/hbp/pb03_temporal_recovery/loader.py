from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import numpy as np
from scipy.io import loadmat

from experiments.hbp.benchmark_utils import build_feature_vector, summarize_array


HERE = Path(__file__).resolve().parent
DEFAULT_DATA_ROOT = Path("/Volumes/Samsung T5/2026/HAOS/HAOS DOCS/DATA/Powergraph")


@dataclass(frozen=True)
class PowerGraphCascadeSource:
    data_root: Path
    family: str
    cascade_root: Path
    file_paths: dict[str, Path]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def choose_powergraph_family(data_root: Path) -> PowerGraphCascadeSource:
    manifest = read_json(HERE / "source_manifest.json")
    selected_cases = list(manifest.get("selected_cases", []))
    for family in selected_cases:
        cascade_root = data_root / "dataset_cascades" / family / family / "raw"
        if cascade_root.exists():
            file_paths = {path.name: path for path in sorted(cascade_root.glob("*.mat"))}
            return PowerGraphCascadeSource(
                data_root=data_root,
                family=family,
                cascade_root=cascade_root,
                file_paths=file_paths,
            )
    raise FileNotFoundError(f"no frozen PB-03 family found under {data_root}")


def mat_shape(path: Path) -> list[int]:
    if not path.exists():
        return []
    with h5py.File(path, "r") as handle:
        keys = [key for key in handle.keys() if key != "#refs#"]
        if not keys:
            return []
        dataset = handle[keys[0]]
        return list(dataset.shape)


def build_loader_summary(data_root: Path = DEFAULT_DATA_ROOT) -> dict[str, Any]:
    source = choose_powergraph_family(data_root)
    file_summaries = {}
    for name, path in source.file_paths.items():
        file_summaries[name] = {
            "exists": path.exists(),
            "size_bytes": path.stat().st_size if path.exists() else 0,
            "shape": mat_shape(path),
        }
    return {
        "bridge_id": "PB-03-TEMPORAL-RECOVERY",
        "selection_status": "frozen",
        "family": source.family,
        "data_root": str(source.data_root),
        "cascade_root": str(source.cascade_root),
        "file_count": len(source.file_paths),
        "file_summaries": file_summaries,
    }


def _read_object_entry(handle: h5py.File, key: str, index: int) -> np.ndarray:
    dataset = handle[key]
    ref = dataset[0, index]
    return np.asarray(handle[ref][()], dtype=float)


def _safe_float(value: Any) -> float:
    try:
        out = float(np.asarray(value).reshape(()))
    except Exception:
        out = float(np.asarray(value).reshape(-1)[0])
    if not np.isfinite(out):
        return 0.0
    return out


def load_case_table(data_root: Path = DEFAULT_DATA_ROOT, *, sample_size: int = 96) -> dict[str, Any]:
    source = choose_powergraph_family(data_root)
    pf_root = data_root / "dataset_pf_opf" / source.family / source.family / "raw"
    with h5py.File(source.cascade_root / "of_reg.mat", "r") as of_reg_file:
        targets = np.asarray(of_reg_file["dns_MW"][()], dtype=float).reshape(-1)
    n_cases = int(targets.size)
    indices = list(range(n_cases))
    ranked = sorted(indices, key=lambda idx: (idx * 2654435761) % 2_147_483_647)
    selected = ranked[:sample_size]

    case_rows: list[dict[str, Any]] = []
    with (
        h5py.File(source.cascade_root / "of_bi.mat", "r") as of_bi_file,
        h5py.File(source.cascade_root / "of_mc.mat", "r") as of_mc_file,
        h5py.File(source.cascade_root / "Bf.mat", "r") as bf_file,
        h5py.File(source.cascade_root / "Ef.mat", "r") as ef_file,
        h5py.File(source.cascade_root / "Ef_nc.mat", "r") as ef_nc_file,
        h5py.File(source.cascade_root / "exp.mat", "r") as exp_file,
        h5py.File(source.cascade_root / "blist.mat", "r") as blist_file,
    ):
        blist = np.asarray(blist_file["bList"][()], dtype=float)
        for idx in selected:
            arrays = {
                "of_bi": _read_object_entry(of_bi_file, "output_features", idx),
                "of_mc": _read_object_entry(of_mc_file, "category", idx),
                "bf": _read_object_entry(bf_file, "B_f_tot", idx),
                "ef": _read_object_entry(ef_file, "E_f_post", idx),
                "ef_nc": _read_object_entry(ef_nc_file, "E_f_kenza", idx),
                "exp": _read_object_entry(exp_file, "explainations", idx),
                "blist": blist,
            }
            case_rows.append(
                {
                    "case_id": f"{source.family}-{idx:05d}",
                    "index": idx,
                    "target": _safe_float(targets[idx]),
                    "arrays": arrays,
                }
            )

    split_order = sorted(
        range(len(case_rows)),
        key=lambda pos: hashlib.sha256(case_rows[pos]["case_id"].encode("utf-8")).hexdigest(),
    )
    dev_cut = sample_size // 2
    cal_cut = dev_cut + sample_size // 4
    splits = {
        "development": sorted(split_order[:dev_cut]),
        "calibration": sorted(split_order[dev_cut:cal_cut]),
        "holdout": sorted(split_order[cal_cut:sample_size]),
    }

    return {
        "bridge_id": "PB-03-TEMPORAL-RECOVERY",
        "family": source.family,
        "case_rows": case_rows,
        "splits": splits,
        "sample_size": sample_size,
        "target_name": "dns_MW",
        "pf_root": str(pf_root),
    }


def build_temporal_feature_table(data_root: Path = DEFAULT_DATA_ROOT, *, sample_size: int = 96) -> dict[str, Any]:
    table = load_case_table(data_root, sample_size=sample_size)
    rows = table["case_rows"]
    features = []
    metadata = []

    def temporal_signature(values: np.ndarray) -> list[float]:
        arr = np.asarray(values, dtype=float).reshape(-1)
        if arr.size == 0:
            return [0.0, 0.0, 0.0, 0.0, 0.0]
        split = max(1, arr.size // 2)
        head = arr[:split]
        tail = arr[-split:]
        delta = tail.mean() - head.mean()
        slope = float((arr[-1] - arr[0]) / max(1, arr.size - 1))
        center = float(np.mean(arr))
        history = float(np.mean(np.diff(arr))) if arr.size > 1 else 0.0
        spread = float(np.std(arr))
        return [center, delta, slope, history, spread]

    for row in rows:
        arrays = row["arrays"]
        pre = np.asarray(arrays["of_bi"], dtype=float).reshape(-1)
        post = np.asarray(arrays["of_mc"], dtype=float).reshape(-1)
        flow = np.asarray(arrays["bf"], dtype=float).reshape(-1)
        response = np.asarray(arrays["ef"], dtype=float).reshape(-1)
        no_cascade = np.asarray(arrays["ef_nc"], dtype=float).reshape(-1)
        explanation = np.asarray(arrays["exp"], dtype=float).reshape(-1)
        branch = np.asarray(arrays["blist"], dtype=float).reshape(-1)
        ordered = [
            pre,
            post,
            flow,
            response,
            no_cascade,
            explanation,
            branch,
        ]
        prefixes = [
            "pre_damage",
            "post_damage",
            "flow_balance",
            "flow_response",
            "no_cascade_response",
            "explanation",
            "branch_list",
        ]
        features.append(build_feature_vector(ordered, prefixes))
        center_pre, delta_pre, slope_pre, history_pre, spread_pre = temporal_signature(pre)
        center_post, delta_post, slope_post, history_post, spread_post = temporal_signature(post)
        center_flow, delta_flow, slope_flow, history_flow, spread_flow = temporal_signature(flow)
        center_response, delta_response, slope_response, history_response, spread_response = temporal_signature(response)
        center_no_cascade, delta_no_cascade, slope_no_cascade, history_no_cascade, spread_no_cascade = temporal_signature(no_cascade)
        center_explanation, delta_explanation, slope_explanation, history_explanation, spread_explanation = temporal_signature(explanation)
        center_branch, delta_branch, slope_branch, history_branch, spread_branch = temporal_signature(branch)
        features[-1] = np.concatenate(
            [
                features[-1],
                np.asarray(
                    [
                        center_pre,
                        delta_pre,
                        slope_pre,
                        history_pre,
                        spread_pre,
                        center_post,
                        delta_post,
                        slope_post,
                        history_post,
                        spread_post,
                        center_flow,
                        delta_flow,
                        slope_flow,
                        history_flow,
                        spread_flow,
                        center_response,
                        delta_response,
                        slope_response,
                        history_response,
                        spread_response,
                        center_no_cascade,
                        delta_no_cascade,
                        slope_no_cascade,
                        history_no_cascade,
                        spread_no_cascade,
                        center_explanation,
                        delta_explanation,
                        slope_explanation,
                        history_explanation,
                        spread_explanation,
                        center_branch,
                        delta_branch,
                        slope_branch,
                        history_branch,
                        spread_branch,
                        float(np.mean(flow) - np.mean(pre)),
                        float(np.mean(response) - np.mean(post)),
                        float(np.mean(no_cascade) - np.mean(flow)),
                        float(np.mean(explanation) - np.mean(branch)),
                    ],
                    dtype=float,
                ),
            ]
        )
        metadata.append(
            {
                "case_id": row["case_id"],
                "target": row["target"],
                "damage_intensity": float(np.mean(np.abs(flow))),
                "recovery_proxy": float(np.mean(response)),
                "branch_persistence": float(np.mean(explanation)),
                "temporal_slope": float(np.mean(np.diff(flow))) if flow.size > 1 else 0.0,
            }
        )
    return {
        "bridge_id": table["bridge_id"],
        "family": table["family"],
        "features": np.asarray(features, dtype=float),
        "targets": np.asarray([row["target"] for row in rows], dtype=float),
        "metadata": metadata,
        "splits": table["splits"],
        "sample_size": table["sample_size"],
        "target_name": table["target_name"],
    }
