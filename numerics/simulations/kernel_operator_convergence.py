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
    "n_sides": [9, 11, 13, 17, 21],
    "epsilon_coeffs": [0.5, 1.0, 2.0],
    "cutoff_factor": 2.5,
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        raw = json.loads(path.read_text())
        if isinstance(raw.get("kernel_operator_convergence"), dict):
            merged.update(raw["kernel_operator_convergence"])
    if config is not None:
        merged.update(config)
    return merged


def build_grid(n_side: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    grid = np.linspace(0.0, 1.0, n_side)
    x, y, z = np.meshgrid(grid, grid, grid, indexing="ij")
    h = float(grid[1] - grid[0])
    return x, y, z, h


def analytic_functions(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> dict[str, dict[str, np.ndarray]]:
    r2 = x * x + y * y + z * z
    f1 = r2
    lap1 = np.full_like(f1, 6.0)

    f2 = np.sin(math.pi * x) * np.sin(math.pi * y) * np.sin(math.pi * z)
    lap2 = -3.0 * math.pi * math.pi * f2

    f3 = np.exp(-r2)
    lap3 = (4.0 * r2 - 6.0) * f3

    return {
        "f1_quadratic": {"values": f1, "laplacian": lap1},
        "f2_trigonometric": {"values": f2, "laplacian": lap2},
        "f3_gaussian": {"values": f3, "laplacian": lap3},
    }


def build_offsets(epsilon_coeff: float, cutoff_factor: float) -> list[tuple[int, int, int, float]]:
    cutoff_sq = cutoff_factor * cutoff_factor * epsilon_coeff
    max_offset = int(math.floor(cutoff_factor * math.sqrt(epsilon_coeff)))
    offsets: list[tuple[int, int, int, float]] = []
    for dx in range(-max_offset, max_offset + 1):
        for dy in range(-max_offset, max_offset + 1):
            for dz in range(-max_offset, max_offset + 1):
                if dx == 0 and dy == 0 and dz == 0:
                    continue
                offset_sq = dx * dx + dy * dy + dz * dz
                if offset_sq <= cutoff_sq + 1.0e-12:
                    weight = math.exp(-offset_sq / epsilon_coeff)
                    offsets.append((dx, dy, dz, weight))
    return offsets


def apply_graph_laplacian(values: np.ndarray, offsets: list[tuple[int, int, int, float]]) -> np.ndarray:
    out = np.zeros_like(values)
    n = values.shape[0]
    for dx, dy, dz, weight in offsets:
        i0 = max(0, -dx)
        i1 = min(n, n - dx)
        j0 = max(0, -dy)
        j1 = min(n, n - dy)
        k0 = max(0, -dz)
        k1 = min(n, n - dz)

        src = (slice(i0, i1), slice(j0, j1), slice(k0, k1))
        nbr = (slice(i0 + dx, i1 + dx), slice(j0 + dy, j1 + dy), slice(k0 + dz, k1 + dz))
        out[src] += weight * (values[src] - values[nbr])
    return out


def induced_operator_scale(h: float, offsets: list[tuple[int, int, int, float]]) -> float:
    moment_x = sum(weight * dx * dx for dx, _, _, weight in offsets)
    moment_y = sum(weight * dy * dy for _, dy, _, weight in offsets)
    moment_z = sum(weight * dz * dz for _, _, dz, weight in offsets)
    mu = (moment_x + moment_y + moment_z) / 3.0
    return -2.0 / (mu * h * h)


def interior_margin(offsets: list[tuple[int, int, int, float]]) -> int:
    return max(max(abs(dx), abs(dy), abs(dz)) for dx, dy, dz, _ in offsets)


def restricted_l2_error(discrete: np.ndarray, analytic: np.ndarray, margin: int, h: float) -> float:
    interior = (
        slice(margin, discrete.shape[0] - margin),
        slice(margin, discrete.shape[1] - margin),
        slice(margin, discrete.shape[2] - margin),
    )
    diff = discrete[interior] - analytic[interior]
    return float(math.sqrt(np.sum(diff * diff) * (h**3)))


def profile_line(values: np.ndarray) -> np.ndarray:
    center_y = values.shape[1] // 2
    center_z = values.shape[2] // 2
    return np.asarray(values[:, center_y, center_z], dtype=float)


def fit_error_scaling(h_values: np.ndarray, errors: np.ndarray) -> dict[str, float]:
    x = np.log(h_values)
    y = np.log(errors)
    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * x + intercept
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / max(ss_tot, 1.0e-18)
    return {"order_p": float(slope), "intercept": float(intercept), "r2": float(r2)}


def run_case(n_side: int, epsilon_coeff: float, cutoff_factor: float) -> dict[str, Any]:
    x, y, z, h = build_grid(n_side)
    funcs = analytic_functions(x, y, z)
    offsets = build_offsets(epsilon_coeff=epsilon_coeff, cutoff_factor=cutoff_factor)
    scale = induced_operator_scale(h=h, offsets=offsets)
    margin = interior_margin(offsets)

    results: dict[str, Any] = {}
    for name, payload in funcs.items():
        values = payload["values"]
        analytic_lap = payload["laplacian"]
        raw = apply_graph_laplacian(values, offsets)
        induced = scale * raw
        error = restricted_l2_error(induced, analytic_lap, margin=margin, h=h)
        results[name] = {
            "error_l2": error,
            "line_discrete": profile_line(induced).tolist(),
            "line_analytic": profile_line(analytic_lap).tolist(),
        }

    return {
        "n_side": n_side,
        "nodes": int(n_side**3),
        "h": h,
        "epsilon_coeff": epsilon_coeff,
        "epsilon_k": epsilon_coeff * h * h,
        "scale_factor": scale,
        "margin": margin,
        "offset_count": len(offsets),
        "functions": results,
    }


def summarize_cases(cases: list[dict[str, Any]]) -> tuple[str, str, dict[str, dict[str, Any]]]:
    grouped: dict[str, dict[str, Any]] = {}
    for func_name in cases[0]["functions"].keys():
        grouped[func_name] = {}
        for eps in sorted({case["epsilon_coeff"] for case in cases}):
            subset = sorted([case for case in cases if case["epsilon_coeff"] == eps], key=lambda item: item["h"])
            h_values = np.array([case["h"] for case in subset], dtype=float)
            errors = np.array([case["functions"][func_name]["error_l2"] for case in subset], dtype=float)
            grouped[func_name][str(eps)] = {
                "h_values": h_values.tolist(),
                "errors": errors.tolist(),
                "fit": fit_error_scaling(h_values, errors),
            }

    quadratic_errors = [
        case["functions"]["f1_quadratic"]["error_l2"]
        for case in cases
    ]
    smooth_orders = [
        grouped[func][eps]["fit"]["order_p"]
        for func in ("f2_trigonometric", "f3_gaussian")
        for eps in grouped[func]
    ]
    local_orders = [
        grouped[func][eps]["fit"]["order_p"]
        for func in ("f2_trigonometric", "f3_gaussian")
        for eps in ("0.5", "1.0")
        if eps in grouped[func]
    ]

    if max(quadratic_errors) < 1.0e-10 and min(local_orders) > 1.0:
        observation = "the kernel-induced operator reproduces the quadratic test exactly and converges toward the continuum Laplacian for smooth test functions"
        conclusion = "after discrete second-moment normalization, the graph operator shows clear Laplacian convergence on the cubic scan; broader kernels remain less accurate at finite resolution"
    elif float(np.min(smooth_orders)) > 0.5:
        observation = "the kernel-induced operator shows positive-order convergence toward the continuum Laplacian across the cubic scan"
        conclusion = "the interaction kernel induces a Laplacian-type operator, though convergence quality depends strongly on the kernel width coefficient"
    else:
        observation = "the kernel-induced operator shows only weak or inconsistent convergence across the current scan"
        conclusion = "the present normalization or stencil regime is not yet giving a clean continuum Laplacian limit"
    return observation, conclusion, grouped


def make_plots(cases: list[dict[str, Any]], grouped: dict[str, dict[str, Any]], stamp: str) -> list[str]:
    plot_paths: list[str] = []

    error_path = PLOTS / f"{stamp}_operator_error_vs_h.png"
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharex=True)
    for ax, func_name in zip(axes, grouped.keys()):
        for eps, payload in grouped[func_name].items():
            h_values = np.array(payload["h_values"], dtype=float)
            errors = np.array(payload["errors"], dtype=float)
            order = payload["fit"]["order_p"]
            ax.loglog(h_values, errors, marker="o", label=f"c_eps={eps}, p={order:.2f}")
        ax.set_title(func_name.replace("_", " "))
        ax.set_xlabel("h")
        ax.set_ylabel(r"$\| \widehat{L}_h f - \Delta f \|_{L^2}$")
        ax.grid(alpha=0.25)
        ax.legend()
    fig.savefig(error_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(error_path.relative_to(REPO_ROOT)))

    profile_path = PLOTS / f"{stamp}_operator_profiles.png"
    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    finest = max(cases, key=lambda item: item["n_side"])
    xline = np.linspace(0.0, 1.0, finest["n_side"])
    for ax, func_name in zip(axes, finest["functions"].keys()):
        analytic = np.array(finest["functions"][func_name]["line_analytic"], dtype=float)
        ax.plot(xline, analytic, color="k", linewidth=2.0, label="analytic")
        for eps in sorted({case["epsilon_coeff"] for case in cases}):
            case = next(item for item in cases if item["n_side"] == finest["n_side"] and item["epsilon_coeff"] == eps)
            discrete = np.array(case["functions"][func_name]["line_discrete"], dtype=float)
            ax.plot(xline, discrete, label=f"c_eps={eps}")
        ax.set_title(func_name.replace("_", " "))
        ax.set_ylabel(r"$\widehat{L}_h f$ / $\Delta f$")
        ax.grid(alpha=0.25)
        ax.legend()
    axes[-1].set_xlabel("central x-line")
    fig.savefig(profile_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(profile_path.relative_to(REPO_ROOT)))

    return plot_paths


def save_results(result: dict[str, Any], stamp: str) -> tuple[str, str]:
    stamped = DATA / f"{stamp}_kernel_operator_convergence.json"
    latest = DATA / "kernel_operator_convergence_latest.json"
    stamped.write_text(json.dumps(result, indent=2))
    latest.write_text(json.dumps(result, indent=2))
    return str(stamped.relative_to(REPO_ROOT)), str(latest.relative_to(REPO_ROOT))


def append_log(result_path: str, plot_paths: list[str], config: dict[str, Any], observation: str, conclusion: str) -> None:
    with EXPERIMENT_LOG.open("a") as handle:
        handle.write("\n## Kernel Laplacian convergence\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            "- Config: "
            f"n_sides={config['n_sides']}, epsilon_coeffs={config['epsilon_coeffs']}, cutoff_factor={config['cutoff_factor']}\n"
        )
        handle.write(f"- Results: `{result_path}`\n")
        handle.write("- Plots: " + ", ".join(f"`{path}`" for path in plot_paths) + "\n")
        handle.write(f"- Observation: {observation}\n")
        handle.write(f"- Conclusion: {conclusion}\n")


def run_kernel_operator_convergence(config: dict[str, Any] | None = None) -> dict[str, Any]:
    cfg = load_config(config)
    cases: list[dict[str, Any]] = []
    for eps in [float(val) for val in cfg["epsilon_coeffs"]]:
        for n_side in [int(val) for val in cfg["n_sides"]]:
            cases.append(run_case(n_side=n_side, epsilon_coeff=eps, cutoff_factor=float(cfg["cutoff_factor"])))

    observation, conclusion, grouped = summarize_cases(cases)
    return {
        "experiment": "kernel_operator_convergence",
        "config": cfg,
        "derivation": {
            "graph_laplacian": "L = D - W with W_ij = exp(-|x_i-x_j|^2 / epsilon_k) on the cutoff graph",
            "induced_operator": "Lhat_h = -(2 / (mu_2 h^2)) L, with mu_2 the discrete second moment of the kernel stencil",
            "target": "Lhat_h f -> Delta f on interior nodes as h -> 0",
        },
        "cases": cases,
        "grouped": grouped,
        "observation": observation,
        "conclusion": conclusion,
    }


def main() -> None:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    result = run_kernel_operator_convergence()
    plots = make_plots(result["cases"], result["grouped"], stamp=stamp)
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
    for func_name, by_eps in result["grouped"].items():
        for eps, payload in by_eps.items():
            print(
                func_name,
                f"c_eps={eps}",
                f"order={payload['fit']['order_p']:.4f}",
                f"R2={payload['fit']['r2']:.4f}",
            )
    print("observation =", result["observation"])
    print("conclusion =", result["conclusion"])


if __name__ == "__main__":
    main()
