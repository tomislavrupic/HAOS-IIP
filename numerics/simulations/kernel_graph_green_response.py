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
    "n_sides": [9, 11, 13],
    "epsilon_coeffs": [0.5, 1.0],
    "cutoff_factor": 2.5,
    "fit_r_max": 0.55,
    "random_seed": 42,
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        raw = json.loads(path.read_text())
        if isinstance(raw.get("kernel_graph_green_response"), dict):
            merged.update(raw["kernel_graph_green_response"])
    if config is not None:
        merged.update(config)
    return merged


def build_cubic_points(n_side: int) -> tuple[np.ndarray, np.ndarray, int]:
    grid = np.linspace(0.0, 1.0, n_side)
    X, Y, Z = np.meshgrid(grid, grid, grid, indexing="ij")
    points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    node_index = np.arange(n_side**3).reshape((n_side, n_side, n_side))
    center = n_side // 2
    source_idx = int(node_index[center, center, center])
    return points, node_index, source_idx


def build_kernel_graph(points: np.ndarray, epsilon_k: float, cutoff_factor: float) -> tuple[np.ndarray, np.ndarray]:
    diff = points[:, None, :] - points[None, :, :]
    d2 = np.sum(diff * diff, axis=2)
    cutoff = cutoff_factor * math.sqrt(epsilon_k)
    mask = (d2 > 0.0) & (d2 <= cutoff * cutoff)
    W = np.zeros_like(d2)
    W[mask] = np.exp(-d2[mask] / (2.0 * epsilon_k))
    degrees = np.sum(W, axis=1)
    L = np.diag(degrees) - W
    return W, L


def solve_mean_zero_poisson(L: np.ndarray, source_idx: int) -> np.ndarray:
    n = L.shape[0]
    s = np.full(n, -1.0 / n, dtype=float)
    s[source_idx] += 1.0
    J = np.ones((n, n), dtype=float) / n
    phi = np.linalg.solve(L + J, s)
    return phi - np.mean(phi)


def radial_profile(points: np.ndarray, values: np.ndarray, source_idx: int, spacing: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    center = points[source_idx]
    r = np.linalg.norm(points - center, axis=1)
    rounded = np.round(r / max(spacing, 1.0e-12), 8) * spacing
    unique = np.unique(rounded)
    centers: list[float] = []
    means: list[float] = []
    counts: list[int] = []
    for radius in unique:
        mask = np.isclose(rounded, radius)
        if np.sum(mask) == 0:
            continue
        centers.append(float(np.mean(r[mask])))
        means.append(float(np.mean(values[mask])))
        counts.append(int(np.sum(mask)))
    order = np.argsort(centers)
    return np.asarray(centers)[order], np.asarray(means)[order], np.asarray(counts)[order]


def fit_inverse_r(r: np.ndarray, phi: np.ndarray, r_min: float, r_max: float) -> dict[str, float]:
    mask = (r >= r_min) & (r <= r_max)
    if np.sum(mask) < 3:
        raise ValueError(f"insufficient fit samples in window [{r_min}, {r_max}]")
    X = np.column_stack([1.0 / r[mask], np.ones(np.sum(mask), dtype=float)])
    coeffs, *_ = np.linalg.lstsq(X, phi[mask], rcond=None)
    pred = X @ coeffs
    ss_res = float(np.sum((phi[mask] - pred) ** 2))
    ss_tot = float(np.sum((phi[mask] - np.mean(phi[mask])) ** 2))
    r2 = 1.0 - ss_res / max(ss_tot, 1.0e-18)
    residual = phi[mask] - float(coeffs[1])
    power = fit_power(r[mask], np.abs(residual))
    return {
        "B_over_r": float(coeffs[0]),
        "A_offset": float(coeffs[1]),
        "r2": float(r2),
        "residual_slope": float(power["slope"]),
    }


def fit_power(r: np.ndarray, y: np.ndarray) -> dict[str, float]:
    mask = (r > 0.0) & (y > 0.0)
    if np.sum(mask) < 3:
        raise ValueError("insufficient positive samples for power-law fit")
    x = np.log(r[mask])
    z = np.log(y[mask])
    slope, intercept = np.polyfit(x, z, 1)
    pred = slope * x + intercept
    ss_res = float(np.sum((z - pred) ** 2))
    ss_tot = float(np.sum((z - np.mean(z)) ** 2))
    r2 = 1.0 - ss_res / max(ss_tot, 1.0e-18)
    return {"slope": float(slope), "intercept": float(intercept), "r2": float(r2)}


def smooth_profile(values: np.ndarray) -> np.ndarray:
    if values.size < 5:
        return values
    kernel = np.array([1.0, 2.0, 3.0, 2.0, 1.0], dtype=float)
    kernel /= np.sum(kernel)
    padded = np.pad(values, (2, 2), mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def run_case(n_side: int, epsilon_coeff: float, cutoff_factor: float, fit_r_max: float) -> dict[str, Any]:
    points, _, source_idx = build_cubic_points(n_side)
    h = 1.0 / (n_side - 1)
    epsilon_k = epsilon_coeff * h * h
    _, L = build_kernel_graph(points, epsilon_k=epsilon_k, cutoff_factor=cutoff_factor)
    phi = solve_mean_zero_poisson(L, source_idx=source_idx)
    r_nodes = np.linalg.norm(points - points[source_idx], axis=1)
    r, phi_profile, counts = radial_profile(points, phi, source_idx=source_idx, spacing=h)
    fit_r_min = max(2.0 * h, 1.5 * math.sqrt(epsilon_k))
    fit = fit_inverse_r(r_nodes[r_nodes > 0.0], phi[r_nodes > 0.0], r_min=fit_r_min, r_max=fit_r_max)
    grad = np.abs(np.gradient(smooth_profile(phi_profile), r))
    force_mask = (r >= fit_r_min) & (r <= fit_r_max)
    force_fit = fit_power(r[force_mask], grad[force_mask])

    return {
        "n_side": n_side,
        "nodes": int(points.shape[0]),
        "lattice_spacing": h,
        "epsilon_coeff": epsilon_coeff,
        "epsilon_k": epsilon_k,
        "fit_r_min": fit_r_min,
        "fit_r_max": fit_r_max,
        "fit": fit,
        "force_fit": force_fit,
        "radial_r": r.tolist(),
        "radial_phi": phi_profile.tolist(),
        "radial_counts": counts.tolist(),
    }


def summarize_observation(cases: list[dict[str, Any]]) -> tuple[str, str]:
    residual_slopes = np.array([case["fit"]["residual_slope"] for case in cases], dtype=float)
    fit_r2 = np.array([case["fit"]["r2"] for case in cases], dtype=float)
    if np.max(np.abs(residual_slopes + 1.0)) < 0.18 and np.min(fit_r2) > 0.98:
        observation = "kernel-graph Green response approaches A + B/r on the cubic substrate scan"
        conclusion = "the weighted interaction kernel induces a graph Laplacian whose far field is consistent with inverse-distance geometry; direct shell-derivative force estimates remain noisier than the field fit"
    else:
        observation = "kernel-graph Green response shows partial inverse-distance structure, but convergence is not yet clean across the scan"
        conclusion = "the graph operator is moving toward the expected Green law, but stronger continuum control is still needed"
    return observation, conclusion


def make_plots(cases: list[dict[str, Any]], stamp: str) -> list[str]:
    plot_paths: list[str] = []

    profile_path = PLOTS / f"{stamp}_kernel_graph_green_profiles.png"
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=False)
    for eps_coeff in sorted({case["epsilon_coeff"] for case in cases}):
        subset = [case for case in cases if case["epsilon_coeff"] == eps_coeff]
        subset = sorted(subset, key=lambda item: item["n_side"])
        rep = subset[-1]
        r = np.array(rep["radial_r"], dtype=float)
        phi = np.array(rep["radial_phi"], dtype=float)
        A = float(rep["fit"]["A_offset"])
        axes[0].loglog(r[1:], np.abs(phi[1:] - A), marker="o", label=f"eps_coeff={eps_coeff}, n={rep['n_side']}")
        fit_curve = np.abs(float(rep["fit"]["B_over_r"]) / np.maximum(r[1:], 1.0e-12))
        axes[0].loglog(r[1:], fit_curve, "--")

        grad = np.abs(np.gradient(phi, r))
        axes[1].loglog(r[1:], grad[1:], marker="o", label=f"eps_coeff={eps_coeff}, n={rep['n_side']}")
        force_curve = np.abs(float(rep["fit"]["B_over_r"])) / np.maximum(r[1:], 1.0e-12) ** 2
        axes[1].loglog(r[1:], force_curve, "--")

    axes[0].set_title("Field profile after offset subtraction")
    axes[0].set_xlabel("r")
    axes[0].set_ylabel(r"$|\phi(r)-A|$")
    axes[0].grid(alpha=0.25)
    axes[0].legend()

    axes[1].set_title("Radial gradient profile")
    axes[1].set_xlabel("r")
    axes[1].set_ylabel(r"$|d\phi/dr|$")
    axes[1].grid(alpha=0.25)
    axes[1].legend()

    fig.savefig(profile_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(profile_path.relative_to(REPO_ROOT)))

    exponent_path = PLOTS / f"{stamp}_kernel_graph_green_exponents.png"
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharex=True)
    for eps_coeff in sorted({case["epsilon_coeff"] for case in cases}):
        subset = sorted([case for case in cases if case["epsilon_coeff"] == eps_coeff], key=lambda item: item["n_side"])
        n_vals = [case["n_side"] for case in subset]
        field_slopes = [case["fit"]["residual_slope"] for case in subset]
        force_slopes = [case["force_fit"]["slope"] for case in subset]
        axes[0].plot(n_vals, field_slopes, marker="o", label=f"eps_coeff={eps_coeff}")
        axes[1].plot(n_vals, force_slopes, marker="o", label=f"eps_coeff={eps_coeff}")

    axes[0].axhline(-1.0, color="k", linestyle="--")
    axes[0].set_title("Field exponent convergence")
    axes[0].set_xlabel("cubic side n")
    axes[0].set_ylabel("slope of |phi-A| vs r")
    axes[0].grid(alpha=0.25)
    axes[0].legend()

    axes[1].axhline(-2.0, color="k", linestyle="--")
    axes[1].set_title("Force exponent convergence")
    axes[1].set_xlabel("cubic side n")
    axes[1].set_ylabel("slope of |dphi/dr| vs r")
    axes[1].grid(alpha=0.25)
    axes[1].legend()

    fig.savefig(exponent_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(exponent_path.relative_to(REPO_ROOT)))

    fit_path = PLOTS / f"{stamp}_kernel_graph_green_fit_quality.png"
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for eps_coeff in sorted({case["epsilon_coeff"] for case in cases}):
        subset = sorted([case for case in cases if case["epsilon_coeff"] == eps_coeff], key=lambda item: item["n_side"])
        n_vals = [case["n_side"] for case in subset]
        fit_r2 = [case["fit"]["r2"] for case in subset]
        ax.plot(n_vals, fit_r2, marker="o", label=f"eps_coeff={eps_coeff}")
    ax.set_xlabel("cubic side n")
    ax.set_ylabel(r"$R^2$ for A + B/r fit")
    ax.set_title("Graph Green fit quality")
    ax.set_ylim(0.9, 1.001)
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(fit_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(fit_path.relative_to(REPO_ROOT)))

    return plot_paths


def append_log(result_path: str, plot_paths: list[str], config: dict[str, Any], observation: str, conclusion: str) -> None:
    with EXPERIMENT_LOG.open("a") as handle:
        handle.write("\n## Kernel graph Green response\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            "- Config: "
            f"n_sides={config['n_sides']}, epsilon_coeffs={config['epsilon_coeffs']}, "
            f"cutoff_factor={config['cutoff_factor']}, fit_r_max={config['fit_r_max']}\n"
        )
        handle.write(f"- Results: `{result_path}`\n")
        handle.write("- Plots: " + ", ".join(f"`{path}`" for path in plot_paths) + "\n")
        handle.write(f"- Observation: {observation}\n")
        handle.write(f"- Conclusion: {conclusion}\n")


def save_results(result: dict[str, Any], stamp: str) -> tuple[str, str]:
    stamped = DATA / f"{stamp}_kernel_graph_green_response.json"
    latest = DATA / "kernel_graph_green_response_latest.json"
    stamped.write_text(json.dumps(result, indent=2))
    latest.write_text(json.dumps(result, indent=2))
    return str(stamped.relative_to(REPO_ROOT)), str(latest.relative_to(REPO_ROOT))


def run_kernel_graph_green_response(config: dict[str, Any] | None = None) -> dict[str, Any]:
    cfg = load_config(config)
    cases: list[dict[str, Any]] = []
    for eps_coeff in [float(v) for v in cfg["epsilon_coeffs"]]:
        for n_side in [int(v) for v in cfg["n_sides"]]:
            cases.append(
                run_case(
                    n_side=n_side,
                    epsilon_coeff=eps_coeff,
                    cutoff_factor=float(cfg["cutoff_factor"]),
                    fit_r_max=float(cfg["fit_r_max"]),
                )
            )

    observation, conclusion = summarize_observation(cases)
    return {
        "experiment": "kernel_graph_green_response",
        "config": cfg,
        "derivation": {
            "graph_operator": "L = D - W with W_ij = exp(-|x_i-x_j|^2 / (2 epsilon_k)) on the cutoff graph",
            "poisson_problem": "L phi = s in the mean-zero subspace",
            "source_handling": "s = delta_source - (1/N) 1",
            "far_field_target": "phi(r) = A + B/r, |d phi / dr| ~ 1/r^2 if the graph Laplacian converges to the 3D continuum Laplacian",
        },
        "cases": cases,
        "observation": observation,
        "conclusion": conclusion,
    }


def main() -> None:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    result = run_kernel_graph_green_response()
    plots = make_plots(result["cases"], stamp=stamp)
    result["plots"] = plots
    result_path, latest_path = save_results(result, stamp=stamp)
    append_log(
        result_path=result_path,
        plot_paths=plots,
        config=result["config"],
        observation=result["observation"],
        conclusion=result["conclusion"],
    )
    print("results =", result_path)
    print("latest =", latest_path)
    print("plots =", plots)
    for case in result["cases"]:
        print(
            "case",
            f"n={case['n_side']}",
            f"eps_coeff={case['epsilon_coeff']}",
            f"field_slope={case['fit']['residual_slope']:.4f}",
            f"force_slope={case['force_fit']['slope']:.4f}",
            f"fit_R2={case['fit']['r2']:.4f}",
        )
    print("observation =", result["observation"])
    print("conclusion =", result["conclusion"])


if __name__ == "__main__":
    main()
