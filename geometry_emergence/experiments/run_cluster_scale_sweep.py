#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path
from typing import Any, Optional

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
MODULE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG_PATH = MODULE_ROOT / "configs" / "cluster_scale_sweep_config.json"
FULL_OUTPUT_PATH = MODULE_ROOT / "diagnostics" / "cluster_scale_sweep_full.json"
SUMMARY_OUTPUT_PATH = MODULE_ROOT / "diagnostics" / "cluster_scale_sweep_summary.json"
PLOT_OUTPUT_PATH = MODULE_ROOT / "diagnostics" / "cluster_scale_transition.svg"
PLOT_PDF_OUTPUT_PATH = MODULE_ROOT / "diagnostics" / "cluster_scale_transition.pdf"
PLOT_PNG_OUTPUT_PATH = MODULE_ROOT / "diagnostics" / "cluster_scale_transition.png"
CANONICAL_SUMMARY_PATH = MODULE_ROOT / "diagnostics" / "canonical_transition_summary.json"

os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from geometry_emergence.figure_protocol import build_cluster_scale_phase_diagram, save_figure_bundle
from geometry_emergence.metrics import (
    compute_effective_dimension,
    compute_flow_bending_index,
    compute_neighborhood_persistence,
    compute_path_distortion,
    compute_recoverability,
)
from geometry_emergence.operators import InteractionGraph, KernelField, build_transport_operator


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


def resolve_config_path(config_arg: Optional[str]) -> Path:
    if config_arg is None:
        return DEFAULT_CONFIG_PATH

    candidate = Path(config_arg)
    candidates = [
        candidate,
        Path.cwd() / candidate,
        MODULE_ROOT / candidate,
        MODULE_ROOT / "configs" / candidate.name,
    ]
    for path in candidates:
        if path.exists():
            return path.resolve()
    raise FileNotFoundError(f"unable to resolve config path: {config_arg}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        default=None,
        help="Optional config path. Supports paths relative to the current directory or geometry_emergence root.",
    )
    return parser.parse_args()


def load_config(config_path: Path) -> dict[str, Any]:
    config = read_json(config_path)
    defaults = {
        "embedding_dim": 2,
        "cluster_count": 3,
        "base_cluster_spread": 0.075,
        "base_center_radius": 0.18,
        "transport_self_weight": 0.05,
        "locality_radius_factor": 3.0,
        "transport_threshold": 0.01,
        "neighborhood_support_threshold": 0.01,
        "k_nearest": 6,
        "burn_in_steps": 6,
        "recovery_steps": 10,
        "perturbation_scale": 0.15,
        "dimension_horizon": 6,
        "dimension_local_k": 8,
        "structural_radius_factor": 1.6,
    }
    for key, value in defaults.items():
        config.setdefault(key, value)
    return config


def reflect_unit_box(positions: np.ndarray) -> np.ndarray:
    reflected = np.asarray(positions, dtype=float).copy()
    for _ in range(4):
        reflected = np.where(reflected < 0.0, -reflected, reflected)
        reflected = np.where(reflected > 1.0, 2.0 - reflected, reflected)
    return np.clip(reflected, 0.0, 1.0)


def build_scaled_clustered_graph(
    graph_size: int,
    kernel_width: float,
    locality_radius: float,
    seed: int,
    cluster_scale: float,
    embedding_dim: int,
    cluster_count: int,
    base_cluster_spread: float,
    base_center_radius: float,
) -> InteractionGraph:
    if embedding_dim != 2:
        raise ValueError("cluster-scale sweep currently supports embedding_dim=2 only")

    rng = np.random.default_rng(seed)
    local_spread = float(base_cluster_spread) * float(cluster_scale)

    # Preserve the overall second moment of the point cloud across scales so
    # graph size and global spatial density remain comparable.
    target_second_moment = float(base_center_radius) ** 2 + 2.0 * float(base_cluster_spread) ** 2
    center_radius_sq = max(target_second_moment - 2.0 * local_spread ** 2, 0.0)
    center_radius = math.sqrt(center_radius_sq)

    center_angles = np.linspace(0.0, 2.0 * math.pi, int(cluster_count), endpoint=False)
    centers = np.column_stack(
        (
            0.5 + center_radius * np.cos(center_angles),
            0.5 + center_radius * np.sin(center_angles),
        )
    )

    counts = np.full(int(cluster_count), int(graph_size) // int(cluster_count), dtype=int)
    counts[: int(graph_size) % int(cluster_count)] += 1
    assignments = np.repeat(np.arange(int(cluster_count)), counts)
    rng.shuffle(assignments)

    positions = centers[assignments] + rng.normal(
        loc=0.0,
        scale=local_spread,
        size=(int(graph_size), int(embedding_dim)),
    )
    positions = reflect_unit_box(positions)

    kernel_field = KernelField(
        kernel_width=float(kernel_width),
        locality_radius=float(locality_radius),
    )
    affinity, distances = kernel_field.build_affinity(positions)
    return InteractionGraph(
        positions=positions,
        affinity=np.asarray(affinity, dtype=float),
        distances=np.asarray(distances, dtype=float),
        kernel_width=float(kernel_width),
        locality_radius=float(locality_radius),
    )


def connected_components_from_adjacency(adjacency: np.ndarray) -> list[list[int]]:
    seen = np.zeros(adjacency.shape[0], dtype=bool)
    components: list[list[int]] = []
    for start in range(adjacency.shape[0]):
        if seen[start]:
            continue
        stack = [start]
        seen[start] = True
        component: list[int] = []
        while stack:
            node = stack.pop()
            component.append(int(node))
            neighbors = np.flatnonzero(adjacency[node])
            for neighbor in neighbors.tolist():
                if not seen[neighbor]:
                    seen[neighbor] = True
                    stack.append(int(neighbor))
        components.append(sorted(component))
    return components


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


def structural_cluster_logs(
    interaction_graph: InteractionGraph,
    base_cluster_spread: float,
    structural_radius_factor: float,
) -> dict[str, float]:
    structural_radius = float(structural_radius_factor) * float(base_cluster_spread)
    structural_adjacency = (
        (interaction_graph.distances <= structural_radius)
        & ~np.eye(interaction_graph.n_nodes, dtype=bool)
    )
    components = connected_components_from_adjacency(structural_adjacency)
    component_sizes = np.array([len(component) for component in components], dtype=float)
    return {
        "mean_cluster_size": float(np.mean(component_sizes)),
        "largest_component_fraction": float(np.max(component_sizes) / interaction_graph.n_nodes),
    }


def _normalize_progress(values: np.ndarray, improving_when_lower: bool) -> np.ndarray:
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


def compute_coupled_progress(rows: list[dict[str, Any]]) -> dict[str, list[float]]:
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
    return {
        "distortion_progress": [round_float(value) for value in distortion_progress.tolist()],
        "persistence_progress": [round_float(value) for value in persistence_progress.tolist()],
        "recoverability_progress": [round_float(value) for value in recoverability_progress.tolist()],
        "coupled_diagnostic_score": [round_float(value) for value in coupled_progress.tolist()],
    }


def evaluate_run(
    config: dict[str, Any],
    cluster_scale: float,
    kernel_width: float,
    seed: int,
) -> dict[str, Any]:
    interaction_graph = build_scaled_clustered_graph(
        graph_size=int(config["graph_size"]),
        kernel_width=float(kernel_width),
        locality_radius=float(config["locality_radius_factor"]) * float(kernel_width),
        seed=int(seed),
        cluster_scale=float(cluster_scale),
        embedding_dim=int(config["embedding_dim"]),
        cluster_count=int(config["cluster_count"]),
        base_cluster_spread=float(config["base_cluster_spread"]),
        base_center_radius=float(config["base_center_radius"]),
    )
    transport_operator = build_transport_operator(
        interaction_graph=interaction_graph,
        self_weight=float(config["transport_self_weight"]),
    )

    distortion = compute_path_distortion(
        interaction_graph=interaction_graph,
        transport_operator=transport_operator,
        max_steps=int(config["transport_steps"]),
        threshold=float(config["transport_threshold"]),
    )
    persistence = compute_neighborhood_persistence(
        transport_operator=transport_operator,
        k_nearest=int(config["k_nearest"]),
        steps=int(config["transport_steps"]),
        support_threshold=float(config["neighborhood_support_threshold"]),
    )
    recoverability = compute_recoverability(
        transport_operator=transport_operator,
        reachable_fraction=float(distortion["reachable_fraction"]),
        seed=int(seed),
        burn_in_steps=int(config["burn_in_steps"]),
        recovery_steps=int(config["recovery_steps"]),
        perturbation_scale=float(config["perturbation_scale"]),
    )
    dimensionality = compute_effective_dimension(
        transport_operator=transport_operator,
        signature_horizon=int(config["dimension_horizon"]),
        local_k=int(config["dimension_local_k"]),
    )
    flow = compute_flow_bending_index(
        transport_operator=transport_operator,
        steps=int(config["transport_steps"]),
    )

    row = {
        "kernel_width": float(kernel_width),
        "seed": int(seed),
        "scale": float(cluster_scale),
        "distortion_score": float(distortion["distortion_score"]),
        "neighborhood_retention": float(persistence["neighborhood_retention"]),
        "recoverability_score": float(recoverability["recoverability_score"]),
        "fraction_participation_score": float(persistence["support_coverage"]),
        "clustering_coefficient_mean": float(
            mean_local_clustering_coefficient(interaction_graph.affinity > 0.0)
        ),
        "effective_dimension": float(dimensionality["effective_dimension"]),
        "flow_bending_index": float(flow["flow_bending_index"]),
    }
    row.update(
        structural_cluster_logs(
            interaction_graph=interaction_graph,
            base_cluster_spread=float(config["base_cluster_spread"]),
            structural_radius_factor=float(config["structural_radius_factor"]),
        )
    )
    return row


def detect_onset(rows: list[dict[str, Any]], coupled_threshold: float) -> dict[str, Any]:
    ordered = sorted(rows, key=lambda item: float(item["kernel_width"]))
    progress = compute_coupled_progress(ordered)
    coupled = np.array(progress["coupled_diagnostic_score"], dtype=float)
    crossing = np.flatnonzero(coupled >= float(coupled_threshold))
    onset_index = None if crossing.size == 0 else int(crossing[0])
    onset_row = None if onset_index is None else ordered[onset_index]
    return {
        "onset_kernel_width": None if onset_row is None else float(onset_row["kernel_width"]),
        "onset_detected": onset_row is not None,
        "onset_index": onset_index,
        "onset_row": None
        if onset_row is None
        else {
            "kernel_width": round_float(onset_row["kernel_width"]),
            "path_distortion_metric": round_float(1.0 - onset_row["distortion_score"]),
            "neighborhood_persistence": round_float(onset_row["neighborhood_retention"]),
            "recoverability_coupling": round_float(onset_row["recoverability_score"]),
            "fraction_participation_score": round_float(onset_row["fraction_participation_score"]),
            "clustering_coefficient_mean": round_float(onset_row["clustering_coefficient_mean"]),
            "effective_dimension": round_float(onset_row["effective_dimension"]),
            "mean_cluster_size": round_float(onset_row["mean_cluster_size"]),
            "largest_component_fraction": round_float(onset_row["largest_component_fraction"]),
            "flow_bending_index": round_float(onset_row["flow_bending_index"]),
        },
        "coupled_progress": progress,
    }


def compute_rank(values: list[float]) -> np.ndarray:
    order = np.argsort(np.asarray(values, dtype=float), kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    ranks[order] = np.arange(len(values), dtype=float)
    return ranks


def onset_scale_monotonicity(per_scale_summary: list[dict[str, Any]]) -> Optional[float]:
    valid = [
        (float(item["scale"]), float(item["mean_onset"]))
        for item in per_scale_summary
        if item["mean_onset"] is not None
    ]
    if len(valid) < 2:
        return None
    scales = [item[0] for item in valid]
    onsets = [item[1] for item in valid]
    scale_ranks = compute_rank(scales)
    onset_ranks = compute_rank(onsets)
    if np.std(scale_ranks) <= EPSILON or np.std(onset_ranks) <= EPSILON:
        return None
    return float(np.corrcoef(scale_ranks, onset_ranks)[0, 1])


def canonical_clustered_variance() -> Optional[float]:
    if not CANONICAL_SUMMARY_PATH.exists():
        return None
    payload = read_json(CANONICAL_SUMMARY_PATH)
    panel = payload.get("graph_panels", {}).get("clustered", {})
    onset_std = panel.get("onset_std")
    if onset_std is None:
        return None
    return float(onset_std) ** 2


def variance_reduction_vs_canonical(per_scale_summary: list[dict[str, Any]]) -> Optional[float]:
    canonical_variance = canonical_clustered_variance()
    if canonical_variance is None or canonical_variance <= EPSILON:
        return None

    scale_variances = [
        float(item["std_onset"]) ** 2
        for item in per_scale_summary
        if item["std_onset"] is not None
    ]
    if not scale_variances:
        return None
    mean_scale_variance = float(np.mean(scale_variances))
    return float(1.0 - mean_scale_variance / canonical_variance)


def build_scale_summary(scale: float, seed_runs: list[dict[str, Any]], seed_count: int) -> dict[str, Any]:
    onset_values = [
        float(seed_run["onset_kernel_width"])
        for seed_run in seed_runs
        if seed_run["onset_kernel_width"] is not None
    ]
    detected_fraction = float(len(onset_values) / max(seed_count, 1))
    if onset_values:
        mean_onset = float(np.mean(onset_values))
        std_onset = float(np.std(onset_values))
        min_onset = float(np.min(onset_values))
        max_onset = float(np.max(onset_values))
    else:
        mean_onset = None
        std_onset = None
        min_onset = None
        max_onset = None

    return {
        "scale": round_float(scale),
        "detected_fraction": round_float(detected_fraction),
        "mean_onset": None if mean_onset is None else round_float(mean_onset),
        "std_onset": None if std_onset is None else round_float(std_onset),
        "min_onset": None if min_onset is None else round_float(min_onset),
        "max_onset": None if max_onset is None else round_float(max_onset),
    }


def plot_transition(full_payload: dict[str, Any], summary_payload: dict[str, Any]) -> None:
    figure = build_cluster_scale_phase_diagram(
        full_payload=full_payload,
        summary_payload=summary_payload,
    )
    save_figure_bundle(
        figure,
        svg_path=PLOT_OUTPUT_PATH,
        pdf_path=PLOT_PDF_OUTPUT_PATH,
        png_path=PLOT_PNG_OUTPUT_PATH,
    )
    figure.clf()


def run_cluster_scale_sweep(config: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    kernel_widths = rounded_widths(
        start=float(config["kernel_width_range"][0]),
        stop=float(config["kernel_width_range"][1]),
        step=float(config["kernel_step"]),
    )

    scale_runs: list[dict[str, Any]] = []
    per_scale_summary: list[dict[str, Any]] = []
    for scale in [float(value) for value in config["cluster_scales"]]:
        seed_runs = []
        for seed in [int(value) for value in config["seeds"]]:
            rows = [
                evaluate_run(
                    config=config,
                    cluster_scale=scale,
                    kernel_width=kernel_width,
                    seed=seed,
                )
                for kernel_width in kernel_widths
            ]
            onset = detect_onset(rows=rows, coupled_threshold=float(config["coupled_threshold"]))
            seed_runs.append(
                {
                    "seed": int(seed),
                    "rows": [
                        {
                            key: (
                                round_float(value)
                                if isinstance(value, float)
                                else value
                            )
                            for key, value in row.items()
                        }
                        for row in rows
                    ],
                    "onset_kernel_width": None
                    if onset["onset_kernel_width"] is None
                    else round_float(onset["onset_kernel_width"]),
                    "onset_detection": onset,
                }
            )

        scale_runs.append(
            {
                "scale": round_float(scale),
                "seed_runs": seed_runs,
            }
        )
        per_scale_summary.append(
            build_scale_summary(
                scale=scale,
                seed_runs=seed_runs,
                seed_count=len(config["seeds"]),
            )
        )

    full_payload = {
        "statement": "Deterministic cluster-scale sweep for geometry emergence transport diagnostics.",
        "config": config,
        "kernel_widths": kernel_widths,
        "scale_runs": scale_runs,
    }

    summary_payload = {
        "statement": "Deterministic cluster-scale sweep summary for transport-geometry diagnostics.",
        "config": {
            "graph_size": int(config["graph_size"]),
            "seeds": [int(value) for value in config["seeds"]],
            "cluster_scales": [round_float(value) for value in config["cluster_scales"]],
            "kernel_width_range": [round_float(value) for value in config["kernel_width_range"]],
            "kernel_step": round_float(config["kernel_step"]),
            "transport_steps": int(config["transport_steps"]),
            "coupled_threshold": float(config["coupled_threshold"]),
        },
        "per_scale": per_scale_summary,
        "global": {
            "onset_scale_monotonicity_score": None,
            "variance_reduction_vs_previous_canonical": None,
        },
    }
    monotonicity = onset_scale_monotonicity(per_scale_summary)
    variance_reduction = variance_reduction_vs_canonical(per_scale_summary)
    if monotonicity is not None:
        summary_payload["global"]["onset_scale_monotonicity_score"] = round_float(monotonicity)
    if variance_reduction is not None:
        summary_payload["global"]["variance_reduction_vs_previous_canonical"] = round_float(
            variance_reduction
        )

    plot_transition(full_payload=full_payload, summary_payload=summary_payload)
    return full_payload, summary_payload


def main() -> None:
    args = parse_args()
    config_path = resolve_config_path(args.config)
    config = load_config(config_path)
    full_payload, summary_payload = run_cluster_scale_sweep(config)
    write_json(FULL_OUTPUT_PATH, full_payload)
    write_json(SUMMARY_OUTPUT_PATH, summary_payload)
    print(json.dumps(summary_payload, indent=2))


if __name__ == "__main__":
    main()
