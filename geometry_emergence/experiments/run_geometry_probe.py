#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Optional


ROOT = Path(__file__).resolve().parents[2]
MODULE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG_PATH = MODULE_ROOT / "experiments" / "default_geometry_probe_config.json"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from geometry_emergence.metrics import (
    compute_effective_dimension,
    compute_flow_bending_index,
    compute_neighborhood_persistence,
    compute_path_distortion,
    compute_recoverability,
)
from geometry_emergence.operators import (
    build_clustered_graph,
    build_random_interaction_graph,
    build_transport_operator,
)


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)
        handle.write("\n")


def load_config(config_path: Optional[Path]) -> dict[str, Any]:
    config = read_json(DEFAULT_CONFIG_PATH)
    if config_path is None:
        return config
    override = read_json(config_path)
    config.update(override)
    return config


def _stable_geometry_score(result: dict[str, Any]) -> float:
    distortion_term = max(0.0, 1.0 - float(result["distortion_score"]))
    persistence_term = max(0.0, float(result["neighborhood_retention"]))
    recoverability_term = max(0.0, float(result["recoverability_score"]))
    return float((distortion_term * persistence_term * recoverability_term) ** (1.0 / 3.0))


def _detect_transition(results: list[dict[str, Any]]) -> dict[str, Any]:
    if not results:
        return {
            "transition_detected": False,
            "transition_kernel_width": None,
            "best_kernel_width": None,
        }

    scored = []
    for result in results:
        score = _stable_geometry_score(result)
        enriched = dict(result)
        enriched["stable_geometry_score"] = score
        scored.append(enriched)

    best_result = max(scored, key=lambda item: float(item["stable_geometry_score"]))
    score_floor = 0.9 * float(best_result["stable_geometry_score"])

    onset = None
    for result in scored:
        if (
            float(result["stable_geometry_score"]) >= score_floor
            and float(result["reachable_fraction"]) >= 0.9
        ):
            onset = result
            break

    return {
        "transition_detected": onset is not None,
        "transition_kernel_width": None if onset is None else float(onset["kernel_width"]),
        "best_kernel_width": float(best_result["kernel_width"]),
        "best_stable_geometry_score": float(best_result["stable_geometry_score"]),
    }


def _build_graph(config: dict[str, Any], kernel_width: float):
    locality_radius = float(config["locality_radius_factor"]) * float(kernel_width)
    common_kwargs = {
        "n_nodes": int(config["n_nodes"]),
        "kernel_width": float(kernel_width),
        "locality_radius": locality_radius,
        "seed": int(config["seed"]),
        "embedding_dim": int(config["embedding_dim"]),
    }

    graph_kind = str(config.get("graph_kind", "clustered"))
    if graph_kind == "random":
        return build_random_interaction_graph(**common_kwargs)
    if graph_kind == "clustered":
        return build_clustered_graph(
            **common_kwargs,
            n_clusters=int(config.get("cluster_count", 3)),
            cluster_spread=float(config.get("cluster_spread", 0.08)),
        )
    raise ValueError(f"unsupported graph_kind: {graph_kind}")


def run_geometry_probe(config: dict[str, Any]) -> dict[str, Any]:
    results: list[dict[str, Any]] = []

    for kernel_width in [float(value) for value in config["kernel_widths"]]:
        interaction_graph = _build_graph(config, kernel_width)
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
        dimensionality = compute_effective_dimension(
            transport_operator=transport_operator,
            signature_horizon=int(config.get("signature_horizon", 6)),
            local_k=int(config.get("dimension_local_k", 8)),
        )
        flow = compute_flow_bending_index(
            transport_operator=transport_operator,
            steps=int(config.get("transport_steps", 12)),
        )
        recoverability = compute_recoverability(
            transport_operator=transport_operator,
            reachable_fraction=float(distortion["reachable_fraction"]),
            seed=int(config.get("seed", 17)),
            burn_in_steps=int(config.get("burn_in_steps", 6)),
            recovery_steps=int(config.get("recovery_steps", 10)),
            perturbation_scale=float(config.get("perturbation_scale", 0.15)),
        )

        result = {
            "kernel_width": float(kernel_width),
            "locality_radius": float(interaction_graph.locality_radius),
            "component_count": int(interaction_graph.connected_component_count()),
            "mean_degree": float(interaction_graph.mean_degree()),
            "edge_density": float(interaction_graph.edge_density()),
            "distortion_score": float(distortion["distortion_score"]),
            "reachable_fraction": float(distortion["reachable_fraction"]),
            "neighborhood_retention": float(persistence["neighborhood_retention"]),
            "step_retention": [float(value) for value in persistence["step_retention"]],
            "local_pca_dimension": float(dimensionality["local_pca_dimension"]),
            "correlation_dimension": float(dimensionality["correlation_dimension"]),
            "effective_dimension": float(dimensionality["effective_dimension"]),
            "local_dimension_trace": [float(value) for value in dimensionality["local_dimension_trace"]],
            "flow_bending_index": float(flow["flow_bending_index"]),
            "step_flow_bending": [float(value) for value in flow["step_flow_bending"]],
            "recoverability_score": float(recoverability["recoverability_score"]),
            "gap_decay_score": float(recoverability["gap_decay_score"]),
            "recoverability_trace": [float(value) for value in recoverability["recoverability_trace"]],
        }
        result["stable_geometry_score"] = _stable_geometry_score(result)
        results.append(result)

    transition = _detect_transition(results)
    output_path = MODULE_ROOT / str(config.get("output_path", "diagnostics/geometry_probe_results.json"))

    payload = {
        "statement": "This module explores when interaction systems acquire geometry-like structure.",
        "module": "geometry_emergence",
        "config": config,
        "results": results,
        "regime_transition": transition,
    }
    write_json(output_path, payload)
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Optional JSON config override.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    payload = run_geometry_probe(config)
    print(json.dumps(payload["regime_transition"], indent=2))


if __name__ == "__main__":
    main()
