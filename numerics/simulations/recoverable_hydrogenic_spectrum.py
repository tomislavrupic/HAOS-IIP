#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import os
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "data"
PLOTS = REPO_ROOT / "plots"
EXPERIMENT_LOG = REPO_ROOT / "experiments" / "EXPERIMENT_LOG.md"
for path in (DATA, PLOTS):
    path.mkdir(exist_ok=True)

MPLCONFIG = REPO_ROOT / ".mplconfig"
MPLCONFIG.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_CONFIG: dict[str, Any] = {
    "D": 1.0,
    "alpha": 1.0,
    "r_min": 1.0e-3,
    "r_max": 80.0,
    "n_grid": 1400,
    "l_values": [0, 1, 2],
    "states_per_l": 4,
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        raw = json.loads(path.read_text())
        if isinstance(raw.get("recoverable_hydrogenic"), dict):
            merged.update(raw["recoverable_hydrogenic"])
    if config is not None:
        merged.update(config)
    return merged


def analytic_energy(n: int, D: float, alpha: float) -> float:
    return -(alpha * alpha) / (4.0 * D * n * n)


def build_radial_hamiltonian(r: np.ndarray, ell: int, D: float, alpha: float) -> tuple[np.ndarray, np.ndarray, float]:
    dr = float(r[1] - r[0])
    rin = r[1:-1]
    centrifugal = D * ell * (ell + 1.0) / (rin * rin)
    potential = -alpha / rin
    main = (2.0 * D / (dr * dr)) + centrifugal + potential
    off = np.full(rin.size - 1, -D / (dr * dr), dtype=float)
    H = np.diag(main) + np.diag(off, 1) + np.diag(off, -1)
    return rin, H, dr


def normalize_radial_mode(chi: np.ndarray, dr: float) -> np.ndarray:
    norm = math.sqrt(float(np.sum(np.abs(chi) ** 2) * dr))
    if norm == 0.0:
        return chi
    return chi / norm


def solve_l_channel(
    r: np.ndarray,
    ell: int,
    D: float,
    alpha: float,
    states_per_l: int,
) -> dict[str, Any]:
    rin, H, dr = build_radial_hamiltonian(r, ell=ell, D=D, alpha=alpha)
    evals, evecs = np.linalg.eigh(H)
    negative_idx = [idx for idx, val in enumerate(evals) if val < 0.0]
    keep = negative_idx[:states_per_l]

    bound_states: list[dict[str, Any]] = []
    for radial_index, idx in enumerate(keep):
        chi = normalize_radial_mode(np.asarray(evecs[:, idx], dtype=float), dr)
        principal_n = radial_index + ell + 1
        theory = analytic_energy(principal_n, D=D, alpha=alpha)
        density = np.abs(chi) ** 2
        bound_states.append(
            {
                "ell": ell,
                "radial_index": radial_index,
                "principal_n": principal_n,
                "energy": float(evals[idx]),
                "analytic_energy": theory,
                "relative_error": abs(float(evals[idx]) - theory) / max(abs(theory), 1.0e-12),
                "r": rin.tolist(),
                "chi": chi.tolist(),
                "density_peak_r": float(rin[int(np.argmax(density))]),
            }
        )

    return {
        "ell": ell,
        "bound_states": bound_states,
        "num_negative_states": len(negative_idx),
    }


def summarize_scaling(channels: list[dict[str, Any]], D: float, alpha: float) -> dict[str, Any]:
    s_wave = next(channel for channel in channels if channel["ell"] == 0)
    n_values = np.array([state["principal_n"] for state in s_wave["bound_states"]], dtype=float)
    energies = np.array([state["energy"] for state in s_wave["bound_states"]], dtype=float)
    predicted = np.array([analytic_energy(int(n), D=D, alpha=alpha) for n in n_values], dtype=float)
    n2_abs_e = np.abs(energies) * (n_values**2)
    return {
        "s_wave_n_values": n_values.tolist(),
        "s_wave_energies": energies.tolist(),
        "predicted_energies": predicted.tolist(),
        "n2_abs_energy": n2_abs_e.tolist(),
        "n2_abs_energy_mean": float(np.mean(n2_abs_e)),
        "n2_abs_energy_std": float(np.std(n2_abs_e)),
    }


def build_observation(channels: list[dict[str, Any]], scaling: dict[str, Any]) -> tuple[str, str]:
    rel_errors = [
        state["relative_error"]
        for channel in channels
        for state in channel["bound_states"]
        if state["principal_n"] <= 4
    ]
    max_error = max(rel_errors) if rel_errors else float("inf")
    std = float(scaling["n2_abs_energy_std"])
    mean = max(float(scaling["n2_abs_energy_mean"]), 1.0e-12)

    if max_error < 0.03 and std / mean < 0.03:
        observation = "bound radial modes form a discrete ladder with near-constant n^2 |E_n|"
        conclusion = "discreteness emerges from quadratic dispersion, inverse-distance geometry, and normalizability without a separate quantization postulate"
    else:
        observation = "bound radial modes are discrete, but numerical scaling is only approximate in the present grid range"
        conclusion = "the substrate supports bound coherent modes, though tighter numerics are needed before claiming clean hydrogenic scaling"
    return observation, conclusion


def make_plots(channels: list[dict[str, Any]], scaling: dict[str, Any], stamp: str) -> list[str]:
    plot_paths: list[str] = []

    scaling_path = PLOTS / f"{stamp}_recoverable_hydrogenic_energy_scaling.png"
    n_vals = np.array(scaling["s_wave_n_values"], dtype=float)
    numerical = np.array(scaling["s_wave_energies"], dtype=float)
    predicted = np.array(scaling["predicted_energies"], dtype=float)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(n_vals, numerical, "o-", label="numerical s-wave")
    ax.plot(n_vals, predicted, "--", label="analytic -1/n^2")
    ax.set_xlabel("principal index n")
    ax.set_ylabel("energy")
    ax.set_title("Recoverable hydrogenic scaling")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(scaling_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(scaling_path.relative_to(REPO_ROOT)))

    invariant_path = PLOTS / f"{stamp}_recoverable_hydrogenic_n2_energy.png"
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(n_vals, np.array(scaling["n2_abs_energy"], dtype=float), "o-")
    ax.axhline(float(scaling["n2_abs_energy_mean"]), color="k", linestyle="--", label="mean")
    ax.set_xlabel("principal index n")
    ax.set_ylabel(r"$n^2 |E_n|$")
    ax.set_title("Scaling invariant for s-wave states")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(invariant_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(invariant_path.relative_to(REPO_ROOT)))

    modes_path = PLOTS / f"{stamp}_recoverable_hydrogenic_modes.png"
    s_wave = next(channel for channel in channels if channel["ell"] == 0)
    fig, axes = plt.subplots(2, 2, figsize=(9, 6), sharex=True)
    axes = axes.ravel()
    for ax, state in zip(axes, s_wave["bound_states"][:4]):
        r = np.array(state["r"], dtype=float)
        density = np.abs(np.array(state["chi"], dtype=float)) ** 2
        ax.plot(r, density, linewidth=1.6)
        ax.set_title(f"n={state['principal_n']}, ell=0")
        ax.grid(alpha=0.25)
    for ax in axes:
        ax.set_xlim(0.0, 40.0)
        ax.set_xlabel("r")
        ax.set_ylabel(r"$|\chi(r)|^2$")
    fig.suptitle("First recoverable radial bound modes")
    fig.tight_layout()
    fig.savefig(modes_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(modes_path.relative_to(REPO_ROOT)))

    return plot_paths


def save_results(result: dict[str, Any], stamp: str) -> tuple[str, str]:
    stamped = DATA / f"{stamp}_recoverable_hydrogenic_spectrum.json"
    latest = DATA / "recoverable_hydrogenic_spectrum_latest.json"
    stamped.write_text(json.dumps(result, indent=2))
    latest.write_text(json.dumps(result, indent=2))
    return str(stamped.relative_to(REPO_ROOT)), str(latest.relative_to(REPO_ROOT))


def append_log(stamp: str, result_path: str, plot_paths: list[str], config: dict[str, Any], observation: str, conclusion: str) -> None:
    with EXPERIMENT_LOG.open("a") as handle:
        handle.write("\n## Recoverable hydrogenic spectrum\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            "- Config: "
            f"D={config['D']}, alpha={config['alpha']}, r_min={config['r_min']}, "
            f"r_max={config['r_max']}, n_grid={config['n_grid']}, "
            f"l_values={config['l_values']}, states_per_l={config['states_per_l']}\n"
        )
        handle.write(f"- Results: `{result_path}`\n")
        handle.write("- Plots: " + ", ".join(f"`{path}`" for path in plot_paths) + "\n")
        handle.write(f"- Observation: {observation}\n")
        handle.write(f"- Conclusion: {conclusion}\n")


def run_recoverable_hydrogenic_spectrum(config: dict[str, Any] | None = None) -> dict[str, Any]:
    cfg = load_config(config)
    D = float(cfg["D"])
    alpha = float(cfg["alpha"])
    r = np.linspace(float(cfg["r_min"]), float(cfg["r_max"]), int(cfg["n_grid"]))
    states_per_l = int(cfg["states_per_l"])
    l_values = [int(val) for val in cfg["l_values"]]

    channels = [
        solve_l_channel(r=r, ell=ell, D=D, alpha=alpha, states_per_l=states_per_l)
        for ell in l_values
    ]
    scaling = summarize_scaling(channels, D=D, alpha=alpha)
    observation, conclusion = build_observation(channels, scaling)

    return {
        "experiment": "recoverable_hydrogenic_spectrum",
        "config": {
            "D": D,
            "alpha": alpha,
            "r_min": float(cfg["r_min"]),
            "r_max": float(cfg["r_max"]),
            "n_grid": int(cfg["n_grid"]),
            "l_values": l_values,
            "states_per_l": states_per_l,
        },
        "derivation": {
            "time_dependent_equation": "d_t psi = i (D Laplacian - V) psi",
            "time_independent_equation": "-D Laplacian u + V u = E u",
            "radial_equation_for_chi": "-D chi'' + [D ell(ell+1)/r^2 - alpha/r] chi = E chi",
            "analytic_scaling": "E_n = -alpha^2 / (4 D n^2)",
        },
        "channels": channels,
        "scaling": scaling,
        "observation": observation,
        "conclusion": conclusion,
    }


def main() -> None:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    result = run_recoverable_hydrogenic_spectrum()
    plot_paths = make_plots(result["channels"], result["scaling"], stamp=stamp)
    result["plots"] = plot_paths
    stamped_path, latest_path = save_results(result, stamp=stamp)
    append_log(
        stamp=stamp,
        result_path=stamped_path,
        plot_paths=plot_paths,
        config=result["config"],
        observation=result["observation"],
        conclusion=result["conclusion"],
    )
    print("results =", stamped_path)
    print("latest =", latest_path)
    print("plots =", plot_paths)
    print("s-wave energies =", [round(v, 6) for v in result["scaling"]["s_wave_energies"]])
    print("observation =", result["observation"])
    print("conclusion =", result["conclusion"])


if __name__ == "__main__":
    main()
