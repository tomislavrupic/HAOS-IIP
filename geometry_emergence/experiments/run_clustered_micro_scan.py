#!/usr/bin/env python3

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
MODULE_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = MODULE_ROOT / "experiments" / "clustered_micro_scan_config.json"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from geometry_emergence.metrics import (
    compute_neighborhood_persistence,
    compute_path_distortion,
    compute_recoverability,
)
from geometry_emergence.operators import build_clustered_graph, build_transport_operator


EPSILON = 1.0e-12
PRIMARY_METRICS = (
    "distortion_score",
    "neighborhood_retention",
    "recoverability_score",
)


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)
        handle.write("\n")


def round_float(value: float, digits: int = 6) -> float:
    return round(float(value), digits)


def rounded_widths(start: float, stop: float, step: float) -> list[float]:
    count = int(np.floor((stop - start) / step + 0.5)) + 1
    return [round_float(start + index * step, 6) for index in range(count)]


def load_config() -> dict[str, Any]:
    config = read_json(CONFIG_PATH)
    base_config = read_json(MODULE_ROOT / str(config["base_config_path"]))
    base_config.update(config)
    return base_config


def build_interaction_graph(
    config: dict[str, Any],
    kernel_width: float,
    seed: int,
):
    return build_clustered_graph(
        n_nodes=int(config["n_nodes"]),
        kernel_width=float(kernel_width),
        locality_radius=float(config["locality_radius_factor"]) * float(kernel_width),
        seed=int(seed),
        embedding_dim=int(config["embedding_dim"]),
        n_clusters=int(config.get("cluster_count", 3)),
        cluster_spread=float(config.get("cluster_spread", 0.08)),
    )


def mean_local_clustering_coefficient(adjacency: np.ndarray) -> float:
    coefficients: list[float] = []
    for node in range(adjacency.shape[0]):
        neighbors = np.flatnonzero(adjacency[node])
        degree = int(neighbors.size)
        if degree < 2:
            coefficients.append(0.0)
            continue
        induced = adjacency[np.ix_(neighbors, neighbors)]
        edge_count = float(np.sum(induced) / 2.0)
        possible_edges = degree * (degree - 1) / 2.0
        coefficients.append(edge_count / max(possible_edges, 1.0))
    return float(np.mean(coefficients))


def structural_logs(interaction_graph) -> dict[str, float]:
    components = interaction_graph.connected_components()
    component_sizes = np.array([len(component) for component in components], dtype=float)
    adjacency = interaction_graph.affinity > 0.0
    return {
        "mean_cluster_size": float(np.mean(component_sizes)),
        "largest_connected_component_fraction": float(np.max(component_sizes) / interaction_graph.n_nodes),
        "clustering_coefficient": mean_local_clustering_coefficient(adjacency),
    }


def evaluate_row(
    config: dict[str, Any],
    kernel_width: float,
    seed: int,
) -> dict[str, Any]:
    interaction_graph = build_interaction_graph(
        config=config,
        kernel_width=kernel_width,
        seed=seed,
    )
    transport_operator = build_transport_operator(
        interaction_graph=interaction_graph,
        self_weight=float(config.get("transport_self_weight", 0.05)),
    )

    distortion = compute_path_distortion(
        interaction_graph=interaction_graph,
        transport_operator=transport_operator,
        max_steps=int(config.get("transport_steps", 12)),
        threshold=float(config.get("transport_threshold", 0.01)),
    )
    persistence = compute_neighborhood_persistence(
        transport_operator=transport_operator,
        k_nearest=int(config.get("k_nearest", 6)),
        steps=int(config.get("transport_steps", 12)),
        support_threshold=float(
            config.get(
                "neighborhood_support_threshold",
                config.get("transport_threshold", 0.01),
            )
        ),
    )
    recoverability = compute_recoverability(
        transport_operator=transport_operator,
        reachable_fraction=float(distortion["reachable_fraction"]),
        seed=int(seed),
        burn_in_steps=int(config.get("burn_in_steps", 6)),
        recovery_steps=int(config.get("recovery_steps", 10)),
        perturbation_scale=float(config.get("perturbation_scale", 0.15)),
    )

    row = {
        "seed": int(seed),
        "kernel_width": float(kernel_width),
        "distortion_score": float(distortion["distortion_score"]),
        "neighborhood_retention": float(persistence["neighborhood_retention"]),
        "recoverability_score": float(recoverability["recoverability_score"]),
        "fraction_participation_score": float(persistence["support_coverage"]),
    }
    row.update(structural_logs(interaction_graph))
    return row


def _normalize_progress(
    values: np.ndarray,
    improving_when_lower: bool,
) -> np.ndarray:
    start = float(values[0])
    end = float(np.min(values)) if improving_when_lower else float(np.max(values))
    scale = abs(end - start)
    if scale <= EPSILON:
        return np.zeros_like(values, dtype=float)
    if improving_when_lower:
        progress = (start - values) / scale
    else:
        progress = (values - start) / scale
    return np.clip(progress, 0.0, 1.0)


def detect_coupled_onset(
    rows: list[dict[str, Any]],
    coupled_threshold: float,
) -> dict[str, Any]:
    ordered = sorted(rows, key=lambda item: float(item["kernel_width"]))
    distortion = np.array([row["distortion_score"] for row in ordered], dtype=float)
    persistence = np.array([row["neighborhood_retention"] for row in ordered], dtype=float)
    recoverability = np.array([row["recoverability_score"] for row in ordered], dtype=float)

    distortion_progress = _normalize_progress(distortion, improving_when_lower=True)
    persistence_progress = _normalize_progress(persistence, improving_when_lower=False)
    recoverability_progress = _normalize_progress(recoverability, improving_when_lower=False)
    coupled_progress = np.minimum.reduce(
        [distortion_progress, persistence_progress, recoverability_progress]
    )
    peak_index = int(np.argmax(coupled_progress))
    peak_row = ordered[peak_index]

    active = coupled_progress >= float(coupled_threshold)
    if not np.any(active):
        return {
            "onset_detected": False,
            "first_crossing": None,
            "reference_row": peak_row,
        }

    onset_index = int(np.flatnonzero(active)[0])
    onset_row = ordered[onset_index]
    return {
        "onset_detected": True,
        "first_crossing": {
            "kernel_width": round_float(onset_row["kernel_width"]),
            "metric_triplet": {
                "distortion_score": round_float(onset_row["distortion_score"]),
                "neighborhood_retention": round_float(onset_row["neighborhood_retention"]),
                "recoverability_score": round_float(onset_row["recoverability_score"]),
            },
            "fraction_participation_score": round_float(
                onset_row["fraction_participation_score"]
            ),
            "structural_logs": {
                "mean_cluster_size": round_float(onset_row["mean_cluster_size"]),
                "largest_connected_component_fraction": round_float(
                    onset_row["largest_connected_component_fraction"]
                ),
                "clustering_coefficient": round_float(
                    onset_row["clustering_coefficient"]
                ),
            },
        },
        "reference_row": onset_row,
    }


def summarize(seed_results: dict[int, list[dict[str, Any]]], config: dict[str, Any]) -> dict[str, Any]:
    seed_summaries = []
    onset_values = []
    missing_onset_seeds = []

    for seed in [int(value) for value in config["seed_panel"]]:
        onset = detect_coupled_onset(
            seed_results[seed],
            coupled_threshold=float(config["coupled_threshold"]),
        )
        if onset["onset_detected"]:
            onset_values.append(float(onset["first_crossing"]["kernel_width"]))
        else:
            missing_onset_seeds.append(int(seed))

        seed_summaries.append(
            {
                "seed": int(seed),
                "onset_detected": bool(onset["onset_detected"]),
                "first_crossing": onset["first_crossing"],
                "structural_reference": {
                    "kernel_width": round_float(onset["reference_row"]["kernel_width"]),
                    "mean_cluster_size": round_float(
                        onset["reference_row"]["mean_cluster_size"]
                    ),
                    "largest_connected_component_fraction": round_float(
                        onset["reference_row"]["largest_connected_component_fraction"]
                    ),
                    "clustering_coefficient": round_float(
                        onset["reference_row"]["clustering_coefficient"]
                    ),
                },
            }
        )

    onset_mean = None
    onset_std = None
    onset_interval = None
    if onset_values:
        onset_mean = round_float(float(np.mean(onset_values)))
        onset_std = round_float(float(np.std(onset_values)))
        onset_interval = [
            round_float(float(np.min(onset_values))),
            round_float(float(np.max(onset_values))),
        ]

    return {
        "graph_kind": str(config["graph_kind"]),
        "kernel_width_window": {
            "start": round_float(float(config["kernel_width_start"])),
            "stop": round_float(float(config["kernel_width_stop"])),
            "step": round_float(float(config["kernel_width_step"])),
        },
        "coupled_threshold": float(config["coupled_threshold"]),
        "seeds": [int(value) for value in config["seed_panel"]],
        "onset_detected_count": int(len(onset_values)),
        "onset_missing_seeds": missing_onset_seeds,
        "onset_values": [round_float(value) for value in onset_values],
        "onset_mean": onset_mean,
        "onset_std": onset_std,
        "onset_interval": onset_interval,
        "seed_summaries": seed_summaries,
    }


def main() -> None:
    config = load_config()
    kernel_widths = rounded_widths(
        start=float(config["kernel_width_start"]),
        stop=float(config["kernel_width_stop"]),
        step=float(config["kernel_width_step"]),
    )
    seed_results = {
        int(seed): [
            evaluate_row(
                config=config,
                kernel_width=kernel_width,
                seed=int(seed),
            )
            for kernel_width in kernel_widths
        ]
        for seed in [int(value) for value in config["seed_panel"]]
    }

    summary = summarize(seed_results=seed_results, config=config)
    output_path = MODULE_ROOT / str(config["summary_output_path"])
    write_json(output_path, summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
