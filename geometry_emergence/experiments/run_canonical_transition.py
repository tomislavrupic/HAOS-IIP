#!/usr/bin/env python3

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
MODULE_ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = MODULE_ROOT / "experiments" / "canonical_transition_config.json"

os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from geometry_emergence.metrics import (
    compute_neighborhood_persistence,
    compute_path_distortion,
    compute_recoverability,
)
from geometry_emergence.operators import (
    build_clustered_graph,
    build_random_interaction_graph,
    build_transport_operator,
)


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


def build_interaction_graph(
    config: dict[str, Any],
    graph_kind: str,
    kernel_width: float,
    seed: int,
):
    locality_radius = float(config["locality_radius_factor"]) * float(kernel_width)
    common_kwargs = {
        "n_nodes": int(config["n_nodes"]),
        "kernel_width": float(kernel_width),
        "locality_radius": locality_radius,
        "seed": int(seed),
        "embedding_dim": int(config["embedding_dim"]),
    }
    if graph_kind == "clustered":
        return build_clustered_graph(
            **common_kwargs,
            n_clusters=int(config.get("cluster_count", 3)),
            cluster_spread=float(config.get("cluster_spread", 0.08)),
        )
    if graph_kind == "random":
        return build_random_interaction_graph(**common_kwargs)
    raise ValueError(f"unsupported graph_kind: {graph_kind}")


def evaluate_primary_metrics(
    config: dict[str, Any],
    graph_kind: str,
    kernel_width: float,
    seed: int,
) -> dict[str, Any]:
    interaction_graph = build_interaction_graph(
        config=config,
        graph_kind=graph_kind,
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

    return {
        "kernel_width": float(kernel_width),
        "graph_kind": graph_kind,
        "seed": int(seed),
        "distortion_score": float(distortion["distortion_score"]),
        "neighborhood_retention": float(persistence["neighborhood_retention"]),
        "recoverability_score": float(recoverability["recoverability_score"]),
        "reachable_fraction": float(distortion["reachable_fraction"]),
    }


def run_seed_sweep(
    config: dict[str, Any],
    graph_kind: str,
    seed: int,
    kernel_widths: list[float],
) -> list[dict[str, Any]]:
    return [
        evaluate_primary_metrics(
            config=config,
            graph_kind=graph_kind,
            kernel_width=kernel_width,
            seed=seed,
        )
        for kernel_width in kernel_widths
    ]


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
    results: list[dict[str, Any]],
    coupled_threshold: float,
) -> dict[str, Any]:
    ordered = sorted(results, key=lambda item: float(item["kernel_width"]))
    kernel_widths = np.array([row["kernel_width"] for row in ordered], dtype=float)
    distortion = np.array([row["distortion_score"] for row in ordered], dtype=float)
    persistence = np.array([row["neighborhood_retention"] for row in ordered], dtype=float)
    recoverability = np.array([row["recoverability_score"] for row in ordered], dtype=float)

    distortion_progress = _normalize_progress(distortion, improving_when_lower=True)
    persistence_progress = _normalize_progress(persistence, improving_when_lower=False)
    recoverability_progress = _normalize_progress(recoverability, improving_when_lower=False)
    coupled_progress = np.minimum.reduce(
        [distortion_progress, persistence_progress, recoverability_progress]
    )

    active = coupled_progress >= float(coupled_threshold)
    onset_kernel_width = None
    onset_index = None
    if np.any(active):
        onset_index = int(np.flatnonzero(active)[0])
        onset_kernel_width = float(kernel_widths[onset_index])

    return {
        "kernel_widths": [round_float(value) for value in kernel_widths.tolist()],
        "coupled_progress": [round_float(value) for value in coupled_progress.tolist()],
        "distortion_progress": [round_float(value) for value in distortion_progress.tolist()],
        "persistence_progress": [round_float(value) for value in persistence_progress.tolist()],
        "recoverability_progress": [round_float(value) for value in recoverability_progress.tolist()],
        "coupled_threshold": float(coupled_threshold),
        "onset_kernel_width": onset_kernel_width,
        "onset_index": onset_index,
    }


def aggregate_seed_panel(
    graph_kind: str,
    seed_runs: dict[int, list[dict[str, Any]]],
    coupled_threshold: float,
) -> dict[str, Any]:
    seeds = sorted(seed_runs)
    kernel_widths = [float(row["kernel_width"]) for row in seed_runs[seeds[0]]]

    per_seed_onsets: list[float] = []
    missing_onset_seeds: list[int] = []
    seed_curves: list[dict[str, Any]] = []
    metric_arrays = {
        metric: np.array(
            [[row[metric] for row in seed_runs[seed]] for seed in seeds],
            dtype=float,
        )
        for metric in PRIMARY_METRICS
    }

    coupled_arrays = []
    for seed in seeds:
        onset = detect_coupled_onset(seed_runs[seed], coupled_threshold=coupled_threshold)
        if onset["onset_kernel_width"] is not None:
            per_seed_onsets.append(float(onset["onset_kernel_width"]))
        else:
            missing_onset_seeds.append(int(seed))
        coupled_arrays.append(np.array(onset["coupled_progress"], dtype=float))
        seed_curves.append(
            {
                "seed": int(seed),
                "onset_kernel_width": onset["onset_kernel_width"],
                "coupled_progress": onset["coupled_progress"],
            }
        )

    coupled_array = np.vstack(coupled_arrays)
    onset_mean = None
    onset_std = None
    onset_interval = None
    if per_seed_onsets:
        onset_mean = float(np.mean(per_seed_onsets))
        onset_std = float(np.std(per_seed_onsets))
        onset_interval = [
            float(np.min(per_seed_onsets)),
            float(np.max(per_seed_onsets)),
        ]

    return {
        "graph_kind": graph_kind,
        "seeds": [int(seed) for seed in seeds],
        "onset_values": [round_float(value) for value in per_seed_onsets],
        "onset_detected_count": int(len(per_seed_onsets)),
        "onset_missing_seeds": missing_onset_seeds,
        "onset_mean": None if onset_mean is None else round_float(onset_mean),
        "onset_std": None if onset_std is None else round_float(onset_std),
        "onset_interval": None
        if onset_interval is None
        else [round_float(value) for value in onset_interval],
        "kernel_widths": [round_float(value) for value in kernel_widths],
        "mean_curves": {
            metric: [round_float(value) for value in np.mean(values, axis=0).tolist()]
            for metric, values in metric_arrays.items()
        },
        "std_curves": {
            metric: [round_float(value) for value in np.std(values, axis=0).tolist()]
            for metric, values in metric_arrays.items()
        },
        "mean_coupled_progress": [
            round_float(value) for value in np.mean(coupled_array, axis=0).tolist()
        ],
        "seed_curves": seed_curves,
    }


def plot_primary_metrics(
    summary: dict[str, Any],
    output_path: Path,
) -> None:
    colors = {
        "clustered": "#1f4e79",
        "random": "#c6651a",
    }
    labels = {
        "distortion_score": "Distortion score",
        "neighborhood_retention": "Neighborhood persistence",
        "recoverability_score": "Recoverability coupling",
    }

    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
    graph_panels = summary["graph_panels"]

    for axis, metric in zip(axes, PRIMARY_METRICS):
        for graph_kind in summary["graph_kinds"]:
            panel = graph_panels[graph_kind]
            kernel_widths = np.array(panel["kernel_widths"], dtype=float)
            mean_values = np.array(panel["mean_curves"][metric], dtype=float)
            std_values = np.array(panel["std_curves"][metric], dtype=float)
            color = colors[graph_kind]

            axis.plot(
                kernel_widths,
                mean_values,
                color=color,
                linewidth=2.2,
                label=graph_kind.capitalize(),
            )
            axis.fill_between(
                kernel_widths,
                mean_values - std_values,
                mean_values + std_values,
                color=color,
                alpha=0.14,
            )

            onset_interval = panel["onset_interval"]
            if onset_interval is not None:
                axis.axvspan(
                    float(onset_interval[0]),
                    float(onset_interval[1]),
                    color=color,
                    alpha=0.06,
                )
                axis.axvline(
                    float(panel["onset_mean"]),
                    color=color,
                    linestyle="--",
                    linewidth=1.2,
                )

        axis.set_ylabel(labels[metric])
        axis.grid(alpha=0.25)

    axes[0].set_title("Canonical Geometry Transition Panel")
    axes[-1].set_xlabel("Kernel width")
    axes[0].legend(frameon=False, ncol=2)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    config = read_json(CONFIG_PATH)

    coarse_results = run_seed_sweep(
        config=config,
        graph_kind=str(config["coarse_reference_graph_kind"]),
        seed=int(config["coarse_reference_seed"]),
        kernel_widths=[float(value) for value in config["coarse_kernel_widths"]],
    )
    coarse_transition = detect_coupled_onset(
        coarse_results,
        coupled_threshold=float(config["coupled_threshold"]),
    )
    coarse_onset = coarse_transition["onset_kernel_width"]
    if coarse_onset is None:
        coarse_onset = float(config["coarse_kernel_widths"][0])

    refine_start = max(
        float(min(config["coarse_kernel_widths"])),
        float(coarse_onset) - float(config["refine_half_width"]),
    )
    refine_stop = min(
        float(max(config["coarse_kernel_widths"])),
        float(coarse_onset) + float(config["refine_half_width"]),
    )
    refined_kernel_widths = rounded_widths(
        start=refine_start,
        stop=refine_stop,
        step=float(config["refine_step"]),
    )

    graph_panels: dict[str, Any] = {}
    for graph_kind in [str(value) for value in config["graph_kinds"]]:
        seed_runs = {
            int(seed): run_seed_sweep(
                config=config,
                graph_kind=graph_kind,
                seed=int(seed),
                kernel_widths=refined_kernel_widths,
            )
            for seed in [int(value) for value in config["seed_panel"]]
        }
        graph_panels[graph_kind] = aggregate_seed_panel(
            graph_kind=graph_kind,
            seed_runs=seed_runs,
            coupled_threshold=float(config["coupled_threshold"]),
        )

    summary = {
        "statement": str(config["statement"]),
        "primary_metrics": list(PRIMARY_METRICS),
        "graph_kinds": [str(value) for value in config["graph_kinds"]],
        "coarse_reference": {
            "graph_kind": str(config["coarse_reference_graph_kind"]),
            "seed": int(config["coarse_reference_seed"]),
            "coarse_kernel_widths": [round_float(value) for value in config["coarse_kernel_widths"]],
            "coarse_onset": None if coarse_transition["onset_kernel_width"] is None else round_float(coarse_transition["onset_kernel_width"]),
            "refined_kernel_widths": refined_kernel_widths,
        },
        "graph_panels": graph_panels,
    }

    summary_output_path = MODULE_ROOT / str(config["summary_output_path"])
    plot_output_path = MODULE_ROOT / str(config["plot_output_path"])

    write_json(summary_output_path, summary)
    plot_primary_metrics(summary=summary, output_path=plot_output_path)

    print(
        json.dumps(
            {
                "coarse_onset": summary["coarse_reference"]["coarse_onset"],
                "refined_window": [
                    summary["coarse_reference"]["refined_kernel_widths"][0],
                    summary["coarse_reference"]["refined_kernel_widths"][-1],
                ],
                "clustered_onset_mean": graph_panels["clustered"]["onset_mean"],
                "random_onset_mean": graph_panels["random"]["onset_mean"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
