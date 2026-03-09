#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import os
import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

PLOTS = REPO_ROOT / "plots"
PLOTS.mkdir(exist_ok=True)
RESULTS = REPO_ROOT / "data"
RESULTS.mkdir(exist_ok=True)
MPLCONFIG = REPO_ROOT / ".mplconfig"
MPLCONFIG.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from numerics.simulations.hodge_modes import build_cubic_lattice_complex
from numerics.simulations.periodic_twisted_l1 import analyze_complex, build_periodic_twisted_complex, plot_edge_mode

DEFAULT_CONFIG: dict[str, Any] = {
    "epsilon": 0.2,
    "periodic_twisted_l1_flux_scan": {
        "sizes": [4, 5],
        "flux_quanta": [0, 1, 2, 3, 4],
        "open_compare_sizes": [4, 5],
        "low_modes": 4,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged["periodic_twisted_l1_flux_scan"] = dict(DEFAULT_CONFIG["periodic_twisted_l1_flux_scan"])
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != "periodic_twisted_l1_flux_scan"})
        if isinstance(on_disk.get("periodic_twisted_l1_flux_scan"), dict):
            merged["periodic_twisted_l1_flux_scan"].update(on_disk["periodic_twisted_l1_flux_scan"])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != "periodic_twisted_l1_flux_scan"})
        if isinstance(config.get("periodic_twisted_l1_flux_scan"), dict):
            merged["periodic_twisted_l1_flux_scan"].update(config["periodic_twisted_l1_flux_scan"])
    return merged


def build_open_case(n_side: int, epsilon: float) -> dict[str, Any]:
    complex_data = build_cubic_lattice_complex(n_side=n_side, epsilon=epsilon)
    complex_data["flux_quanta"] = 0
    complex_data["plaquette_angle"] = 0.0
    complex_data["face_holonomies"] = np.ones(len(complex_data["face_weights"]), dtype=complex)
    return analyze_complex(complex_data, f"open_n{n_side}")


def select_low_coexact_branch(case: dict[str, Any], low_modes: int) -> dict[str, Any]:
    candidates = list(case["full_modes"][:low_modes])
    coexact_candidates = [rec for rec in candidates if float(rec["coexact_fraction"]) >= 0.8]
    selected_pool = coexact_candidates if coexact_candidates else candidates
    selected = min(
        selected_pool,
        key=lambda rec: (
            float(rec["eigenvalue"]),
            float(rec["divergence_norm"]),
            -float(rec["coexact_fraction"]),
        ),
    )
    selected_index = int(selected["mode_index"])
    return {
        "mode_index": selected_index,
        "eigenvalue": float(selected["eigenvalue"]),
        "exact_fraction": float(selected["exact_fraction"]),
        "coexact_fraction": float(selected["coexact_fraction"]),
        "divergence_norm": float(selected["divergence_norm"]),
        "curl_norm": float(selected["curl_norm"]),
        "ipr": float(selected["ipr"]),
        "support_pattern": str(selected["support_pattern"]),
    }


def sanitize_case(case: dict[str, Any]) -> dict[str, Any]:
    return {
        "label": case["label"],
        "config": case["config"],
        "face_holonomy_phase": float(case["face_holonomy_phase"]),
        "full_spectrum": [float(value) for value in case["full_spectrum"]],
        "projected_spectrum": [float(value) for value in case["projected_spectrum"]],
        "full_modes": case["full_modes"],
        "projected_modes": case["projected_modes"],
    }


def make_plots(result: dict[str, Any], plot_name: str) -> list[str]:
    plot_paths: list[str] = []
    periodic_cases = result["periodic_cases"]
    open_cases = result["open_cases"]
    low_modes = int(result["config"]["low_modes"])
    sizes = [int(value) for value in result["config"]["sizes"]]

    spectral_path = PLOTS / f"{plot_name}_spectral_flow.png"
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharex=False)
    for n_side in sizes:
        labels = sorted(
            [label for label in periodic_cases if periodic_cases[label]["config"]["n_side"] == n_side],
            key=lambda label: periodic_cases[label]["config"]["flux_quanta"],
        )
        flux_values = [periodic_cases[label]["config"]["plaquette_angle"] for label in labels]
        for mode_index in range(low_modes):
            axes[0].plot(
                flux_values,
                [periodic_cases[label]["full_spectrum"][mode_index] for label in labels],
                marker="o",
                alpha=0.75,
                label=f"N={n_side**3}, full {mode_index}" if mode_index == 0 else None,
            )
            axes[1].plot(
                flux_values,
                [periodic_cases[label]["projected_spectrum"][mode_index] for label in labels],
                marker="o",
                alpha=0.75,
                label=f"N={n_side**3}, proj {mode_index}" if mode_index == 0 else None,
            )
    axes[0].set_title("Full L1 spectral flow")
    axes[1].set_title("Coexact-projected spectral flow")
    for ax in axes:
        ax.set_xlabel("plaquette angle")
        ax.set_ylabel("eigenvalue")
        ax.grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(spectral_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(spectral_path.relative_to(REPO_ROOT)))

    fraction_path = PLOTS / f"{plot_name}_fractions.png"
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharex=True, sharey=True)
    for ax, metric in zip(axes, ["exact_fraction", "coexact_fraction"]):
        for n_side in sizes:
            labels = sorted(
                [label for label in periodic_cases if periodic_cases[label]["config"]["n_side"] == n_side],
                key=lambda label: periodic_cases[label]["config"]["flux_quanta"],
            )
            flux_values = [periodic_cases[label]["config"]["plaquette_angle"] for label in labels]
            values = [periodic_cases[label]["selected_branch"][metric] for label in labels]
            ax.plot(flux_values, values, marker="o", label=f"N={n_side**3}")
        ax.set_title(metric.replace("_", " "))
        ax.set_xlabel("plaquette angle")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("branch fraction")
    axes[1].legend(fontsize=8)
    fig.savefig(fraction_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(fraction_path.relative_to(REPO_ROOT)))

    divergence_path = PLOTS / f"{plot_name}_divergence.png"
    fig, ax = plt.subplots(figsize=(8, 4))
    for n_side in sizes:
        labels = sorted(
            [label for label in periodic_cases if periodic_cases[label]["config"]["n_side"] == n_side],
            key=lambda label: periodic_cases[label]["config"]["flux_quanta"],
        )
        flux_values = [periodic_cases[label]["config"]["plaquette_angle"] for label in labels]
        values = [periodic_cases[label]["selected_branch"]["divergence_norm"] for label in labels]
        ax.plot(flux_values, values, marker="o", label=f"N={n_side**3}")
    ax.set_xlabel("plaquette angle")
    ax.set_ylabel(r"$||d_0^* a||$")
    ax.set_title("Selected coexact branch divergence norm")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(divergence_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(divergence_path.relative_to(REPO_ROOT)))

    curl_path = PLOTS / f"{plot_name}_curl.png"
    fig, ax = plt.subplots(figsize=(8, 4))
    for n_side in sizes:
        labels = sorted(
            [label for label in periodic_cases if periodic_cases[label]["config"]["n_side"] == n_side],
            key=lambda label: periodic_cases[label]["config"]["flux_quanta"],
        )
        flux_values = [periodic_cases[label]["config"]["plaquette_angle"] for label in labels]
        values = [periodic_cases[label]["selected_branch"]["curl_norm"] for label in labels]
        ax.plot(flux_values, values, marker="o", label=f"N={n_side**3}")
    ax.set_xlabel("plaquette angle")
    ax.set_ylabel(r"$||d_1 a||$")
    ax.set_title("Selected coexact branch curl norm")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(curl_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(curl_path.relative_to(REPO_ROOT)))

    persistence_path = PLOTS / f"{plot_name}_persistence.png"
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharex=True)
    for n_side in sizes:
        labels = sorted(
            [label for label in periodic_cases if periodic_cases[label]["config"]["n_side"] == n_side],
            key=lambda label: periodic_cases[label]["config"]["flux_quanta"],
        )
        flux_values = [periodic_cases[label]["config"]["plaquette_angle"] for label in labels]
        axes[0].plot(
            flux_values,
            [periodic_cases[label]["selected_branch"]["eigenvalue"] for label in labels],
            marker="o",
            label=f"N={n_side**3}",
        )
        axes[1].plot(
            flux_values,
            [periodic_cases[label]["selected_branch"]["coexact_fraction"] for label in labels],
            marker="o",
            label=f"N={n_side**3}",
        )
    axes[0].set_title("Selected branch eigenvalue persistence")
    axes[1].set_title("Selected branch coexact persistence")
    for ax in axes:
        ax.set_xlabel("plaquette angle")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("eigenvalue")
    axes[1].set_ylabel("coexact fraction")
    fig.savefig(persistence_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(persistence_path.relative_to(REPO_ROOT)))

    baseline_path = PLOTS / f"{plot_name}_open_vs_periodic.png"
    fig, axes = plt.subplots(1, len(sizes), figsize=(7 * len(sizes), 4), sharey=True)
    if len(sizes) == 1:
        axes = [axes]
    for ax, n_side in zip(axes, sizes):
        open_label = f"open_n{n_side}"
        periodic_label = f"periodic_n{n_side}_flux0"
        ax.plot(range(low_modes), open_cases[open_label]["full_spectrum"][:low_modes], marker="o", label="open")
        ax.plot(range(low_modes), periodic_cases[periodic_label]["full_spectrum"][:low_modes], marker="o", label="periodic flux=0")
        ax.set_title(f"Open vs periodic, N={n_side**3}")
        ax.set_xlabel("mode index")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("eigenvalue")
    fig.savefig(baseline_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(baseline_path.relative_to(REPO_ROOT)))

    mode_path = PLOTS / f"{plot_name}_modes.png"
    representative_fluxes = []
    available_flux = sorted(set(int(periodic_cases[label]["config"]["flux_quanta"]) for label in periodic_cases if periodic_cases[label]["config"]["n_side"] == max(sizes)))
    if available_flux:
        representative_fluxes = [available_flux[0], available_flux[len(available_flux) // 2], available_flux[-1]]
    fig = plt.figure(figsize=(14, 8))
    panel_index = 1
    for n_side in sizes:
        for flux in representative_fluxes:
            label = f"periodic_n{n_side}_flux{flux}"
            case = periodic_cases[label]
            vector = case["full_vectors"][case["selected_branch"]["mode_index"]]
            ax = fig.add_subplot(len(sizes), len(representative_fluxes), panel_index, projection="3d")
            plot_edge_mode(
                ax,
                case["midpoints"],
                case["directions"],
                vector,
                f"N={n_side**3}, flux={flux}, coexact={case['selected_branch']['coexact_fraction']:.2f}",
            )
            panel_index += 1
    fig.savefig(mode_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(mode_path.relative_to(REPO_ROOT)))

    return plot_paths


def summarize_verdict(periodic_cases: dict[str, Any], sizes: list[int], flux_quanta: list[int]) -> tuple[str, str, dict[str, Any]]:
    branch_continuous = True
    branch_mostly_coexact = True
    robust_across_sizes = True
    still_topological = True
    low_branch_records: dict[str, list[dict[str, float]]] = {}

    for n_side in sizes:
        labels = sorted(
            [label for label in periodic_cases if periodic_cases[label]["config"]["n_side"] == n_side],
            key=lambda label: periodic_cases[label]["config"]["flux_quanta"],
        )
        records = [periodic_cases[label]["selected_branch"] for label in labels]
        low_branch_records[str(n_side)] = [
            {
                "flux_quanta": int(periodic_cases[label]["config"]["flux_quanta"]),
                "plaquette_angle": float(periodic_cases[label]["config"]["plaquette_angle"]),
                "eigenvalue": float(records[idx]["eigenvalue"]),
                "coexact_fraction": float(records[idx]["coexact_fraction"]),
            }
            for idx, label in enumerate(labels)
        ]
        if any(record["coexact_fraction"] < 0.8 for record in records):
            branch_mostly_coexact = False
            robust_across_sizes = False
        support_patterns = [str(record["support_pattern"]) for record in records if int(record["mode_index"]) < 2]
        if any("localized" in pattern for pattern in support_patterns):
            still_topological = False
        eigenvalues = [float(record["eigenvalue"]) for record in records]
        if len(eigenvalues) >= 3:
            total_range = max(max(eigenvalues) - min(eigenvalues), 1e-9)
            second_diff = max(abs(eigenvalues[idx + 1] - 2.0 * eigenvalues[idx] + eigenvalues[idx - 1]) for idx in range(1, len(eigenvalues) - 1))
            if second_diff > 0.8 * total_range:
                branch_continuous = False
        zero_case = next(periodic_cases[label] for label in labels if periodic_cases[label]["config"]["flux_quanta"] == 0)
        max_flux_case = next(periodic_cases[label] for label in labels if periodic_cases[label]["config"]["flux_quanta"] == max(flux_quanta))
        if abs(max_flux_case["selected_branch"]["eigenvalue"] - zero_case["selected_branch"]["eigenvalue"]) < 0.02:
            branch_continuous = False
            robust_across_sizes = False

    if branch_continuous and branch_mostly_coexact and robust_across_sizes:
        observation = "the low mostly coexact L1 branch moves smoothly under flux and remains stable across the scanned sizes"
    elif branch_mostly_coexact and robust_across_sizes:
        observation = "the low L1 branch stays mostly coexact across the scan, but its ordering under flux is not yet a single clean smooth band"
    else:
        observation = "the low periodic L1 branch changes under flux, but its low mostly coexact character is not stable across the full scan"

    if still_topological:
        conclusion = "the scanned branch is robust and mostly coexact, but it still looks dominated by harmonic or topological edge structure rather than a clean propagating vector band"
    else:
        conclusion = "the scan shows some departure from pure torus-cycle behavior, but no isolated propagating Maxwell-like band is established yet"

    verdict = {
        "moves_continuously_under_flux": branch_continuous,
        "stays_mostly_coexact": branch_mostly_coexact,
        "robust_across_sizes": robust_across_sizes,
        "still_harmonic_or_topological": still_topological,
        "low_branch_records": low_branch_records,
    }
    return observation, conclusion, verdict


def save_results(result: dict[str, Any], plot_paths: list[str], plot_name: str) -> tuple[Path, list[str], str]:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output = dict(result)
    output["plots"] = plot_paths
    payload = json.dumps(output, indent=2)
    stamped = RESULTS / f"{timestamp}_{plot_name}.json"
    latest = RESULTS / f"{plot_name}_latest.json"
    stamped.write_text(payload, encoding="utf-8")
    latest.write_text(payload, encoding="utf-8")

    stamped_plots: list[str] = []
    for rel_path in plot_paths:
        src = REPO_ROOT / rel_path
        dst = PLOTS / f"{timestamp}_{plot_name}_{src.name}"
        shutil.copy2(src, dst)
        stamped_plots.append(str(dst.relative_to(REPO_ROOT)))
    return stamped, stamped_plots, timestamp


def run_periodic_twisted_l1_flux_scan(
    config: dict[str, Any] | None = None,
    config_path: Path | None = None,
    plot_name: str = "periodic_twisted_l1_flux_scan",
) -> tuple[dict[str, Any], Path, list[str], str]:
    cfg = load_config(config, config_path=config_path)
    epsilon = float(cfg.get("epsilon", 0.2))
    scan_cfg = cfg["periodic_twisted_l1_flux_scan"]
    sizes = [int(value) for value in scan_cfg.get("sizes", [4, 5])]
    flux_quanta = [int(value) for value in scan_cfg.get("flux_quanta", [0, 1, 2, 3, 4])]
    open_compare_sizes = [int(value) for value in scan_cfg.get("open_compare_sizes", sizes)]
    low_modes = int(scan_cfg.get("low_modes", 4))

    open_cases: dict[str, Any] = {}
    periodic_cases: dict[str, Any] = {}

    for n_side in open_compare_sizes:
        open_label = f"open_n{n_side}"
        open_cases[open_label] = build_open_case(n_side=n_side, epsilon=epsilon)

    for n_side in sizes:
        for flux in flux_quanta:
            label = f"periodic_n{n_side}_flux{flux}"
            periodic_cases[label] = build_periodic_twisted_complex(n_side=n_side, epsilon=epsilon, flux_quanta=flux)
            periodic_cases[label] = analyze_complex(periodic_cases[label], label)
            periodic_cases[label]["selected_branch"] = select_low_coexact_branch(periodic_cases[label], low_modes=low_modes)

    plot_ready = {
        "config": {
            "epsilon": epsilon,
            "sizes": sizes,
            "flux_quanta": flux_quanta,
            "open_compare_sizes": open_compare_sizes,
            "low_modes": low_modes,
        },
        "open_cases": open_cases,
        "periodic_cases": periodic_cases,
    }
    plot_paths = make_plots(plot_ready, plot_name=plot_name)
    observation, conclusion, verdict = summarize_verdict(periodic_cases, sizes, flux_quanta)
    result = {
        "experiment": "periodic_twisted_l1_flux_scan",
        "config": plot_ready["config"],
        "open_cases": {label: sanitize_case(case) for label, case in open_cases.items()},
        "periodic_cases": {label: sanitize_case(case) | {"selected_branch": case["selected_branch"]} for label, case in periodic_cases.items()},
        "verdict": verdict,
        "observation": observation,
        "conclusion": conclusion,
    }
    result_path, stamped_plots, timestamp = save_results(result, plot_paths, plot_name=plot_name)
    return result, result_path, stamped_plots, timestamp


def main() -> None:
    result, result_path, stamped_plots, _ = run_periodic_twisted_l1_flux_scan()
    print(f"Saved: {result_path.relative_to(REPO_ROOT)}")
    print(f"Plots: {', '.join(stamped_plots)}")
    print(f"Observation: {result['observation']}")
    print(f"Conclusion: {result['conclusion']}")
    print(f"Verdict: {result['verdict']}")


if __name__ == "__main__":
    main()
