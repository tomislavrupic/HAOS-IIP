#!/usr/bin/env python3

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

PLOTS = REPO_ROOT / "plots"
PLOTS.mkdir(exist_ok=True)
MPLCONFIG = REPO_ROOT / ".mplconfig"
MPLCONFIG.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from numerics.simulations.hodge_modes import run_hodge_test
from numerics.simulations.laplacian_modes import load_config as load_base_config
from numerics.simulations.laplacian_modes import run_laplacian_test

DEFAULT_SWEEP: dict[str, Any] = {
    "epsilon_values": [0.12, 0.2, 0.28],
    "node_values": [216, 343],
    "substrates": ["random_geometric", "cubic_lattice"],
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = load_base_config(config_path=config_path)
    sweep_cfg = DEFAULT_SWEEP.copy()
    existing = merged.get("sweep", {})
    if isinstance(existing, dict):
        sweep_cfg.update(existing)
    merged["sweep"] = sweep_cfg
    if config is not None:
        merged.update(config)
        if isinstance(config.get("sweep"), dict):
            merged["sweep"].update(config["sweep"])
    return merged


def lattice_side_from_nodes(nodes: int) -> int | None:
    side = round(nodes ** (1.0 / 3.0))
    return side if side**3 == nodes else None


def make_summary_plots(laplacian_records: list[dict[str, Any]], hodge_records: list[dict[str, Any]], name: str) -> list[str]:
    plot_paths: list[str] = []

    if laplacian_records:
        gap_path = PLOTS / f"{name}_laplacian_gap.png"
        fig, ax = plt.subplots(figsize=(8, 4))
        grouped: dict[tuple[str, int], list[dict[str, Any]]] = {}
        for record in laplacian_records:
            key = (str(record["substrate"]), int(record["nodes_requested"]))
            grouped.setdefault(key, []).append(record)
        for (substrate, nodes), group in sorted(grouped.items()):
            group.sort(key=lambda item: item["epsilon"])
            ax.plot(
                [item["epsilon"] for item in group],
                [item["spectral_gap"] for item in group],
                marker="o",
                label=f"{substrate}, N={nodes}",
            )
        ax.set_xlabel("epsilon")
        ax.set_ylabel("lambda_1")
        ax.set_title("Laplacian spectral gap sweep")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
        fig.savefig(gap_path, dpi=180, bbox_inches="tight")
        plt.close(fig)
        plot_paths.append(str(gap_path.relative_to(REPO_ROOT)))

        ipr_path = PLOTS / f"{name}_laplacian_ipr.png"
        fig, ax = plt.subplots(figsize=(8, 4))
        grouped = {}
        for record in laplacian_records:
            key = (str(record["substrate"]), int(record["nodes_requested"]))
            grouped.setdefault(key, []).append(record)
        for (substrate, nodes), group in sorted(grouped.items()):
            group.sort(key=lambda item: item["epsilon"])
            ax.plot(
                [item["epsilon"] for item in group],
                [item["mean_ipr_first3"] for item in group],
                marker="o",
                label=f"{substrate}, N={nodes}",
            )
        ax.set_xlabel("epsilon")
        ax.set_ylabel("mean IPR (first 3 modes)")
        ax.set_title("Laplacian localization sweep")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
        fig.savefig(ipr_path, dpi=180, bbox_inches="tight")
        plt.close(fig)
        plot_paths.append(str(ipr_path.relative_to(REPO_ROOT)))

    if hodge_records:
        hodge_path = PLOTS / f"{name}_hodge_coexact.png"
        fig, ax = plt.subplots(figsize=(8, 4))
        grouped_hodge: dict[int, list[dict[str, Any]]] = {}
        for record in hodge_records:
            grouped_hodge.setdefault(int(record["nodes_requested"]), []).append(record)
        for nodes, group in sorted(grouped_hodge.items()):
            group.sort(key=lambda item: item["epsilon"])
            ax.plot(
                [item["epsilon"] for item in group],
                [item["lowest_mode_coexact_fraction"] for item in group],
                marker="o",
                label=f"cubic_lattice, N={nodes}",
            )
        ax.set_xlabel("epsilon")
        ax.set_ylabel("coexact fraction of first nonzero L1 mode")
        ax.set_title("Hodge transverse content sweep")
        ax.set_ylim(0.0, 1.05)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
        fig.savefig(hodge_path, dpi=180, bbox_inches="tight")
        plt.close(fig)
        plot_paths.append(str(hodge_path.relative_to(REPO_ROOT)))

    return plot_paths


def run_parameter_sweep(config: dict[str, Any] | None = None, plot_name: str = "parameter_sweep") -> dict[str, Any]:
    cfg = load_config(config)
    sweep_cfg = cfg["sweep"]
    epsilon_values = [float(value) for value in sweep_cfg.get("epsilon_values", DEFAULT_SWEEP["epsilon_values"])]
    node_values = [int(value) for value in sweep_cfg.get("node_values", DEFAULT_SWEEP["node_values"])]
    substrates = [str(value) for value in sweep_cfg.get("substrates", DEFAULT_SWEEP["substrates"])]

    laplacian_records: list[dict[str, Any]] = []
    hodge_records: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []

    for substrate in substrates:
        for nodes in node_values:
            for epsilon in epsilon_values:
                lap_cfg = dict(cfg)
                lap_cfg.update({"substrate": substrate, "nodes": nodes, "epsilon": epsilon})
                if substrate == "cubic_lattice":
                    side = lattice_side_from_nodes(nodes)
                    if side is None:
                        skipped.append({"branch": "laplacian", "substrate": substrate, "nodes": nodes, "epsilon": epsilon, "reason": "node count is not a perfect cube"})
                        continue
                    lap_cfg["lattice_side"] = side
                result = run_laplacian_test(lap_cfg, with_plots=False)
                laplacian_records.append(
                    {
                        "substrate": substrate,
                        "nodes_requested": nodes,
                        "nodes_effective": int(result["config"]["nodes"]),
                        "epsilon": epsilon,
                        "spectral_gap": float(result["spectrum"]["spectral_gap"]),
                        "mean_ipr_first3": float(result["spectrum"]["mean_ipr_first3"]),
                        "conclusion": result["conclusion"],
                    }
                )

                if substrate != "cubic_lattice":
                    continue
                side = lattice_side_from_nodes(nodes)
                if side is None:
                    continue
                hodge_cfg = dict(cfg)
                hodge_cfg.update({"epsilon": epsilon, "hodge_substrate": "cubic_lattice", "hodge_lattice_side": side})
                hodge_result = run_hodge_test(hodge_cfg, with_plots=False)
                first_mode = hodge_result["spectrum"]["mode_records"][0]
                hodge_records.append(
                    {
                        "substrate": "cubic_lattice",
                        "nodes_requested": nodes,
                        "nodes_effective": int(hodge_result["config"]["nodes"]),
                        "epsilon": epsilon,
                        "lowest_mode_label": str(first_mode["label"]),
                        "lowest_mode_exact_fraction": float(first_mode["exact_fraction"]),
                        "lowest_mode_coexact_fraction": float(first_mode["coexact_fraction"]),
                        "conclusion": hodge_result["conclusion"],
                    }
                )

    plot_paths = make_summary_plots(laplacian_records, hodge_records, plot_name)

    cubic_records = [record for record in laplacian_records if record["substrate"] == "cubic_lattice"]
    random_records = [record for record in laplacian_records if record["substrate"] == "random_geometric"]
    hodge_transverse = sum(1 for record in hodge_records if record["lowest_mode_label"] == "transverse-like")

    if cubic_records and random_records:
        avg_cubic_ipr = float(np.mean([record["mean_ipr_first3"] for record in cubic_records]))
        avg_random_ipr = float(np.mean([record["mean_ipr_first3"] for record in random_records]))
        if avg_cubic_ipr < avg_random_ipr and hodge_transverse > 0:
            observation = "cubic substrates stay smoother than random graphs, and some L1 low modes become transverse-like"
            conclusion = "the sweep supports a geometry branch on regular substrates and a separate candidate edge-vector branch"
        elif avg_cubic_ipr < avg_random_ipr:
            observation = "cubic substrates stay smoother than random graphs across the sweep"
            conclusion = "substrate regularity matters strongly for clean low-mode geometry, while the vector branch remains only partially separated"
        else:
            observation = "low-mode localization remains sensitive to both epsilon and substrate choice"
            conclusion = "the sweep does not yet isolate a universally clean continuum regime"
    else:
        observation = "parameter sweep completed on the requested grid"
        conclusion = "the sweep provides scaling data, but cross-substrate comparison is incomplete"

    return {
        "experiment": "parameter_sweep",
        "config": {
            "epsilon_values": epsilon_values,
            "node_values": node_values,
            "substrates": substrates,
        },
        "laplacian_records": laplacian_records,
        "hodge_records": hodge_records,
        "skipped": skipped,
        "plots": plot_paths,
        "observation": observation,
        "conclusion": conclusion,
    }


def main() -> None:
    result = run_parameter_sweep()
    print("grid =", result["config"])
    print("laplacian cases =", len(result["laplacian_records"]))
    print("hodge cases =", len(result["hodge_records"]))
    print("plots =", result["plots"])
    print("observation =", result["observation"])
    print("conclusion =", result["conclusion"])


if __name__ == "__main__":
    main()
