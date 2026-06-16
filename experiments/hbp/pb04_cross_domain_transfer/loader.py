from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from experiments.biology.gene_network_demo.gene_network_model import (
    apply_edge_weakening,
    build_default_network,
    final_window_mean,
    normalized_distance,
    simulate_network,
)
from experiments.hbp.benchmark_utils import build_feature_vector


HERE = Path(__file__).resolve().parent
DEFAULT_DATA_ROOT = Path("/Volumes/Samsung T5/2026/HAOS/HAOS DOCS/DATA/Powergraph")


@dataclass(frozen=True)
class CrossDomainSource:
    source_root: Path
    target_root: Path
    source_file_count: int
    target_file_count: int


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def build_loader_summary(data_root: Path = DEFAULT_DATA_ROOT) -> dict[str, Any]:
    source_manifest = read_json(HERE / "source_manifest.json")
    source_path = Path("/Volumes/Samsung T5/2026/HAOS/HAOS DOCS/HAOS-IIP/experiments/biology/gene_network_demo")
    target_root = data_root / "dataset_pf_opf" / "ieee24" / "ieee24" / "raw"
    source_files = list(source_path.glob("*"))
    target_files = list(target_root.glob("*")) if target_root.exists() else []
    return {
        "bridge_id": "PB-04-CROSS-DOMAIN-TRANSFER",
        "selection_status": "frozen",
        "source_domain": source_manifest["source_domain_candidate"]["path"],
        "target_domain": source_manifest["target_domain_candidate"]["path"],
        "source_file_count": len(source_files),
        "target_file_count": len(target_files),
        "source_exists": source_path.exists(),
        "target_exists": target_root.exists(),
    }


def build_transfer_tables(data_root: Path = DEFAULT_DATA_ROOT) -> dict[str, Any]:
    """Build cross-domain recovery descriptors with a matched semantic vocabulary.

    The source side is the deterministic biology proxy. The target side is the
    frozen PowerGraph holdout bundle, but only its recovery-oriented descriptors
    are used here to avoid handing the transfer model direct topology cues.
    """
    source_network = build_default_network()
    source_baseline = simulate_network(source_network)
    source_rows: list[dict[str, Any]] = []
    source_features: list[np.ndarray] = []
    source_targets: list[float] = []
    for idx, level in enumerate(np.linspace(0.0, 1.0, 31)):
        perturbed = apply_edge_weakening(source_network, float(level))
        perturbed_trajectory = simulate_network(perturbed)
        response = final_window_mean(perturbed_trajectory, 30)
        baseline = final_window_mean(source_baseline, 30)
        delta = response - baseline
        feature_arrays = [
            response,
            baseline,
            delta,
            np.abs(delta),
            np.asarray([float(np.mean(delta))], dtype=float),
            np.asarray([float(np.std(delta))], dtype=float),
        ]
        source_features.append(build_feature_vector(feature_arrays, ["response", "baseline", "delta", "abs_delta", "delta_mean", "delta_std"]))
        source_targets.append(float(np.clip(1.0 - normalized_distance(response, baseline), 0.0, 1.0)))
        source_rows.append(
            {
                "case_id": f"biology-{idx:03d}",
                "perturbation_level": float(level),
                "recovery_target": source_targets[-1],
            }
        )

    from experiments.hbp.pb02_external_network_recovery.pb02_holdout import load_feature_bundle

    bundle = load_feature_bundle(data_root)
    target_features: list[np.ndarray] = []
    target_targets: list[float] = []
    target_rows: list[dict[str, Any]] = []
    for idx, case_id in enumerate(bundle.case_ids):
        baseline_row = np.asarray(bundle.baseline[idx], dtype=float)
        haos_row = np.asarray(bundle.haos[idx], dtype=float)
        shared = min(baseline_row.size, haos_row.size)
        delta = haos_row[:shared] - baseline_row[:shared]
        feature_arrays = [
            haos_row,
            baseline_row,
            delta,
            np.abs(delta),
            np.asarray([float(np.mean(delta))], dtype=float),
            np.asarray([float(np.std(delta))], dtype=float),
        ]
        target_features.append(build_feature_vector(feature_arrays, ["haos", "baseline", "delta", "abs_delta", "delta_mean", "delta_std"]))
        target_targets.append(float(bundle.target[idx]))
        target_rows.append(
            {
                "case_id": case_id,
                "target": float(bundle.target[idx]),
            }
        )

    return {
        "bridge_id": "PB-04-CROSS-DOMAIN-TRANSFER",
        "source_features": np.asarray(source_features, dtype=float),
        "source_targets": np.asarray(source_targets, dtype=float),
        "source_rows": source_rows,
        "target_features": np.asarray(target_features, dtype=float),
        "target_targets": np.asarray(target_targets, dtype=float),
        "target_rows": target_rows,
        "source_case_count": len(source_rows),
        "target_case_count": len(target_rows),
    }
