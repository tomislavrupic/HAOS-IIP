#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import os
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
from scipy.sparse import coo_matrix, csc_matrix
from scipy.sparse.linalg import eigsh, expm_multiply, spsolve

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "data"
PLOTS = REPO_ROOT / "plots"
EXPERIMENT_LOG = REPO_ROOT / "experiments" / "EXPERIMENT_LOG.md"
NOTE_PATH = REPO_ROOT / "experiments" / "operators" / "Scalar_Kernel_Graph_Geometry_Bridge_v1.md"
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
    "epsilon_coeff": 0.5,
    "cutoff_factor": 2.5,
    "green_fit_r_max": 0.55,
    "heat_time_max": 0.03,
    "heat_time_steps": 61,
    "msd_fit_window": [0.004, 0.02],
    "arrival_threshold_fraction": 0.5,
    "arrival_scan_shells": 9,
    "low_mode_count": 8,
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        raw = json.loads(path.read_text())
        if isinstance(raw.get("scalar_kernel_graph_geometry_bridge"), dict):
            merged.update(raw["scalar_kernel_graph_geometry_bridge"])
    if config is not None:
        merged.update(config)
    return merged


def build_cubic_points(n_side: int) -> tuple[np.ndarray, np.ndarray, int, float]:
    grid = np.linspace(0.0, 1.0, n_side)
    x, y, z = np.meshgrid(grid, grid, grid, indexing="ij")
    points = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
    node_index = np.arange(n_side**3).reshape((n_side, n_side, n_side))
    source_idx = int(node_index[n_side // 2, n_side // 2, n_side // 2])
    h = float(grid[1] - grid[0])
    return points, node_index, source_idx, h


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
                    offsets.append((dx, dy, dz, math.exp(-offset_sq / epsilon_coeff)))
    return offsets


def build_sparse_laplacian(n_side: int, node_index: np.ndarray, offsets: list[tuple[int, int, int, float]]) -> coo_matrix:
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    degrees = np.zeros(n_side**3, dtype=float)

    for dx, dy, dz, weight in offsets:
        i0 = max(0, -dx)
        i1 = min(n_side, n_side - dx)
        j0 = max(0, -dy)
        j1 = min(n_side, n_side - dy)
        k0 = max(0, -dz)
        k1 = min(n_side, n_side - dz)

        src = node_index[i0:i1, j0:j1, k0:k1].ravel()
        dst = node_index[i0 + dx:i1 + dx, j0 + dy:j1 + dy, k0 + dz:k1 + dz].ravel()

        rows.extend(src.tolist())
        cols.extend(dst.tolist())
        data.extend(((-weight) * np.ones(src.size, dtype=float)).tolist())
        degrees[src] += weight

    diagonal = np.arange(n_side**3, dtype=int)
    rows.extend(diagonal.tolist())
    cols.extend(diagonal.tolist())
    data.extend(degrees.tolist())
    return coo_matrix((data, (rows, cols)), shape=(n_side**3, n_side**3)).tocsr()


def induced_operator_scale(h: float, offsets: list[tuple[int, int, int, float]]) -> float:
    moment_x = sum(weight * dx * dx for dx, _, _, weight in offsets)
    moment_y = sum(weight * dy * dy for _, dy, _, weight in offsets)
    moment_z = sum(weight * dz * dz for _, _, dz, weight in offsets)
    mu = (moment_x + moment_y + moment_z) / 3.0
    return -2.0 / (mu * h * h)


def build_scaled_operator(n_side: int, epsilon_coeff: float, cutoff_factor: float) -> tuple[np.ndarray, int, float, float, coo_matrix, coo_matrix]:
    points, node_index, source_idx, h = build_cubic_points(n_side)
    offsets = build_offsets(epsilon_coeff=epsilon_coeff, cutoff_factor=cutoff_factor)
    laplacian = build_sparse_laplacian(n_side=n_side, node_index=node_index, offsets=offsets)
    scale = induced_operator_scale(h=h, offsets=offsets)
    delta = laplacian * scale
    epsilon_k = epsilon_coeff * h * h
    return points, source_idx, h, epsilon_k, laplacian, delta


def fit_power(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    mask = (x > 0.0) & (y > 0.0)
    if np.sum(mask) < 3:
        raise ValueError("insufficient positive samples for power-law fit")
    log_x = np.log(x[mask])
    log_y = np.log(y[mask])
    slope, intercept = np.polyfit(log_x, log_y, 1)
    pred = slope * log_x + intercept
    ss_res = float(np.sum((log_y - pred) ** 2))
    ss_tot = float(np.sum((log_y - np.mean(log_y)) ** 2))
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "r2": float(1.0 - ss_res / max(ss_tot, 1.0e-18)),
    }


def solve_green_field(delta: coo_matrix, source_idx: int) -> np.ndarray:
    n = delta.shape[0]
    source = np.full(n, -1.0 / n, dtype=float)
    source[source_idx] += 1.0
    gauge = csc_matrix(np.full((n, n), 1.0 / n, dtype=float))
    phi = spsolve((-delta + gauge).tocsc(), source)
    return np.asarray(phi, dtype=float) - float(np.mean(phi))


def green_metrics(
    points: np.ndarray,
    source_idx: int,
    h: float,
    epsilon_k: float,
    phi: np.ndarray,
    fit_r_max: float,
) -> dict[str, Any]:
    radii = np.linalg.norm(points - points[source_idx], axis=1)
    fit_r_min = max(2.0 * h, 1.5 * math.sqrt(epsilon_k))
    mask = (radii >= fit_r_min) & (radii <= fit_r_max)
    x = np.column_stack([1.0 / radii[mask], np.ones(np.sum(mask), dtype=float)])
    coeffs, *_ = np.linalg.lstsq(x, phi[mask], rcond=None)
    pred = x @ coeffs
    ss_res = float(np.sum((phi[mask] - pred) ** 2))
    ss_tot = float(np.sum((phi[mask] - np.mean(phi[mask])) ** 2))
    fit = fit_power(radii[mask], np.abs(phi[mask] - float(coeffs[1])))
    shell_r, shell_phi, _ = radial_profile(points, phi, source_idx=source_idx, spacing=h)
    return {
        "fit_r_min": float(fit_r_min),
        "fit_r_max": float(fit_r_max),
        "A_offset": float(coeffs[1]),
        "B_over_r": float(coeffs[0]),
        "fit_r2": float(1.0 - ss_res / max(ss_tot, 1.0e-18)),
        "residual_slope": float(fit["slope"]),
        "effective_dimension_hint": float(2.0 - fit["slope"]),
        "radial_r": shell_r.tolist(),
        "radial_phi": shell_phi.tolist(),
    }


def radial_profile(
    points: np.ndarray,
    values: np.ndarray,
    source_idx: int,
    spacing: float,
) -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    center = points[source_idx]
    radii = np.linalg.norm(points - center, axis=1)
    rounded = np.round(radii / max(spacing, 1.0e-12), 8) * spacing
    unique = np.unique(rounded)
    shell_radii: list[float] = []
    shell_values: list[float] = []
    shell_indices: list[np.ndarray] = []
    for radius in unique:
        mask = np.isclose(rounded, radius)
        if np.sum(mask) == 0:
            continue
        shell_radii.append(float(np.mean(radii[mask])))
        shell_values.append(float(np.mean(values[mask])))
        shell_indices.append(np.flatnonzero(mask))
    order = np.argsort(shell_radii)
    return (
        np.asarray(shell_radii, dtype=float)[order],
        np.asarray(shell_values, dtype=float)[order],
        [shell_indices[int(i)] for i in order],
    )


def heat_trajectory(delta: coo_matrix, source_idx: int, time_max: float, time_steps: int) -> tuple[np.ndarray, np.ndarray]:
    initial = np.zeros(delta.shape[0], dtype=float)
    initial[source_idx] = 1.0
    times = np.linspace(0.0, time_max, time_steps)
    trajectory = expm_multiply(delta, initial, start=times[0], stop=times[-1], num=len(times), endpoint=True)
    return times, np.asarray(trajectory, dtype=float)


def heat_metrics(
    points: np.ndarray,
    source_idx: int,
    times: np.ndarray,
    trajectory: np.ndarray,
    fit_window: list[float],
) -> dict[str, Any]:
    radii_sq = np.sum((points - points[source_idx]) ** 2, axis=1)
    msd = trajectory @ radii_sq
    mask = (times >= float(fit_window[0])) & (times <= float(fit_window[1]))
    slope, intercept = np.polyfit(times[mask], msd[mask], 1)
    pred = slope * times[mask] + intercept
    ss_res = float(np.sum((msd[mask] - pred) ** 2))
    ss_tot = float(np.sum((msd[mask] - np.mean(msd[mask])) ** 2))
    return {
        "msd_slope": float(slope),
        "msd_intercept": float(intercept),
        "fit_r2": float(1.0 - ss_res / max(ss_tot, 1.0e-18)),
        "effective_diffusivity": float(slope / 6.0),
        "times": times.tolist(),
        "mean_squared_radius": msd.tolist(),
    }


def shell_arrival_metrics(
    shell_radii: np.ndarray,
    shell_indices: list[np.ndarray],
    times: np.ndarray,
    trajectory: np.ndarray,
    threshold_fraction: float,
    scan_shells: int,
) -> dict[str, Any]:
    rows: list[dict[str, float]] = []
    max_shell = min(scan_shells + 1, len(shell_indices))
    for shell_idx in range(1, max_shell):
        series = np.mean(trajectory[:, shell_indices[shell_idx]], axis=1)
        peak = float(np.max(series))
        threshold = threshold_fraction * peak
        hits = np.flatnonzero(series >= threshold)
        if hits.size == 0:
            continue
        rows.append(
            {
                "radius": float(shell_radii[shell_idx]),
                "radius_squared": float(shell_radii[shell_idx] ** 2),
                "arrival_time": float(times[int(hits[0])]),
                "peak_shell_mean": peak,
            }
        )

    radius_sq = np.asarray([row["radius_squared"] for row in rows], dtype=float)
    arrivals = np.asarray([row["arrival_time"] for row in rows], dtype=float)
    slope, intercept = np.polyfit(radius_sq, arrivals, 1)
    pred = slope * radius_sq + intercept
    ss_res = float(np.sum((arrivals - pred) ** 2))
    ss_tot = float(np.sum((arrivals - np.mean(arrivals)) ** 2))
    return {
        "threshold_fraction": float(threshold_fraction),
        "fit_slope": float(slope),
        "fit_intercept": float(intercept),
        "fit_r2": float(1.0 - ss_res / max(ss_tot, 1.0e-18)),
        "rows": rows,
    }


def low_mode_metrics(points: np.ndarray, delta: coo_matrix, low_mode_count: int) -> dict[str, Any]:
    eigenvalues, eigenvectors = eigsh((-delta).astype(float), k=low_mode_count, which="SM")
    order = np.argsort(eigenvalues)
    eigenvalues = np.asarray(eigenvalues[order], dtype=float)
    eigenvectors = np.asarray(eigenvectors[:, order], dtype=float)

    first_shell_vals = eigenvalues[1:4]
    coordinate_basis, _ = np.linalg.qr(points - np.mean(points, axis=0))
    mode_basis, _ = np.linalg.qr(eigenvectors[:, 1:4])
    overlap = coordinate_basis.T @ mode_basis
    principal_cosines = np.linalg.svd(overlap, compute_uv=False)
    return {
        "eigenvalues": eigenvalues.tolist(),
        "first_shell_eigenvalues": first_shell_vals.tolist(),
        "triplet_relative_spread": float((np.max(first_shell_vals) - np.min(first_shell_vals)) / np.mean(first_shell_vals)),
        "principal_cosines": principal_cosines.tolist(),
        "subspace_affinity": float(np.mean(principal_cosines**2)),
        "coordinate_overlap_matrix": np.abs(overlap).tolist(),
    }


def run_case(n_side: int, epsilon_coeff: float, cutoff_factor: float, config: dict[str, Any]) -> dict[str, Any]:
    points, source_idx, h, epsilon_k, _, delta = build_scaled_operator(
        n_side=n_side,
        epsilon_coeff=epsilon_coeff,
        cutoff_factor=cutoff_factor,
    )
    phi = solve_green_field(delta=delta, source_idx=source_idx)
    green = green_metrics(
        points=points,
        source_idx=source_idx,
        h=h,
        epsilon_k=epsilon_k,
        phi=phi,
        fit_r_max=float(config["green_fit_r_max"]),
    )
    times, trajectory = heat_trajectory(
        delta=delta,
        source_idx=source_idx,
        time_max=float(config["heat_time_max"]),
        time_steps=int(config["heat_time_steps"]),
    )
    shell_radii, _, shell_indices = radial_profile(points, phi, source_idx=source_idx, spacing=h)
    heat = heat_metrics(
        points=points,
        source_idx=source_idx,
        times=times,
        trajectory=trajectory,
        fit_window=[float(v) for v in config["msd_fit_window"]],
    )
    shell_arrival = shell_arrival_metrics(
        shell_radii=shell_radii,
        shell_indices=shell_indices,
        times=times,
        trajectory=trajectory,
        threshold_fraction=float(config["arrival_threshold_fraction"]),
        scan_shells=int(config["arrival_scan_shells"]),
    )
    low_mode = low_mode_metrics(points=points, delta=delta, low_mode_count=int(config["low_mode_count"]))
    return {
        "n_side": int(n_side),
        "nodes": int(points.shape[0]),
        "lattice_spacing": float(h),
        "epsilon_coeff": float(epsilon_coeff),
        "epsilon_k": float(epsilon_k),
        "green": green,
        "heat": heat,
        "shell_arrival": shell_arrival,
        "low_mode": low_mode,
    }


def summarize_observation(cases: list[dict[str, Any]]) -> tuple[str, str]:
    green_fit_r2 = np.asarray([case["green"]["fit_r2"] for case in cases], dtype=float)
    green_dim = np.asarray([case["green"]["effective_dimension_hint"] for case in cases], dtype=float)
    heat_fit_r2 = np.asarray([case["heat"]["fit_r2"] for case in cases], dtype=float)
    arrival_fit_r2 = np.asarray([case["shell_arrival"]["fit_r2"] for case in cases], dtype=float)
    low_cos = np.asarray([min(case["low_mode"]["principal_cosines"]) for case in cases], dtype=float)

    if (
        np.min(green_fit_r2) > 0.97
        and np.max(np.abs(green_dim - 3.0)) < 0.06
        and np.min(heat_fit_r2) > 0.995
        and np.min(arrival_fit_r2) > 0.97
        and np.min(low_cos) > 0.99
    ):
        observation = (
            "the shared-family scalar kernel-graph carrier now gives one coherent 3D-like geometry read across all four channels: "
            "inverse-distance Green response, linear diffusion on shell second moments, shell-arrival times linear in radius squared, "
            "and a first low-mode triplet that matches the Euclidean coordinate subspace"
        )
        conclusion = (
            "the previous geometry-closure preflight mismatch is removed on the tested cubic scalar carrier: one bounded metric-like reconstruction "
            "now organizes Green response, heat behavior, shell-arrival structure, and low-mode organization on the same operator/substrate family. "
            "This is a carrier-level bridge, not yet a universality or full geometry-closure claim beyond the tested cubic local-kernel regime."
        )
    else:
        observation = (
            "the scalar kernel-graph carrier now supports a single shared-family geometry test, but at least one of the four channels is still not yet "
            "cleanly organized by the same bounded metric reconstruction"
        )
        conclusion = (
            "the category mistake is fixed, but the first shared-family bridge does not yet support a positive geometry-style read across all selected channels"
        )
    return observation, conclusion


def make_plots(cases: list[dict[str, Any]], stamp: str) -> list[str]:
    plot_paths: list[str] = []

    summary_path = PLOTS / f"{stamp}_scalar_kernel_graph_geometry_bridge_summary.png"
    n_vals = [case["n_side"] for case in cases]
    green_dim = [case["green"]["effective_dimension_hint"] for case in cases]
    green_r2 = [case["green"]["fit_r2"] for case in cases]
    heat_diff = [case["heat"]["effective_diffusivity"] for case in cases]
    arrival_r2 = [case["shell_arrival"]["fit_r2"] for case in cases]
    low_cos = [min(case["low_mode"]["principal_cosines"]) for case in cases]
    low_affinity = [case["low_mode"]["subspace_affinity"] for case in cases]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    axes[0, 0].plot(n_vals, green_dim, marker="o", label="dimension hint")
    axes[0, 0].axhline(3.0, color="k", linestyle="--", linewidth=1.0)
    axes[0, 0].set_title("Green field geometry hint")
    axes[0, 0].set_xlabel("cubic side n")
    axes[0, 0].set_ylabel("effective dimension")
    axes[0, 0].grid(alpha=0.25)
    axes[0, 0].legend()

    axes[0, 1].plot(n_vals, heat_diff, marker="o", label="effective diffusivity")
    axes[0, 1].axhline(1.0, color="k", linestyle="--", linewidth=1.0)
    axes[0, 1].set_title("Heat diffusion closure")
    axes[0, 1].set_xlabel("cubic side n")
    axes[0, 1].set_ylabel("D_eff")
    axes[0, 1].grid(alpha=0.25)
    axes[0, 1].legend()

    axes[1, 0].plot(n_vals, green_r2, marker="o", label="Green fit $R^2$")
    axes[1, 0].plot(n_vals, arrival_r2, marker="s", label="arrival fit $R^2$")
    axes[1, 0].set_title("Metric-fit quality")
    axes[1, 0].set_xlabel("cubic side n")
    axes[1, 0].set_ylabel("$R^2$")
    axes[1, 0].set_ylim(0.95, 1.001)
    axes[1, 0].grid(alpha=0.25)
    axes[1, 0].legend()

    axes[1, 1].plot(n_vals, low_cos, marker="o", label="min principal cosine")
    axes[1, 1].plot(n_vals, low_affinity, marker="s", label="subspace affinity")
    axes[1, 1].set_title("First-shell low-mode organization")
    axes[1, 1].set_xlabel("cubic side n")
    axes[1, 1].set_ylabel("subspace alignment")
    axes[1, 1].set_ylim(0.95, 1.001)
    axes[1, 1].grid(alpha=0.25)
    axes[1, 1].legend()

    fig.savefig(summary_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(summary_path.relative_to(REPO_ROOT)))

    rep = cases[-1]
    detail_path = PLOTS / f"{stamp}_scalar_kernel_graph_geometry_bridge_detail.png"
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    radial_r = np.asarray(rep["green"]["radial_r"], dtype=float)
    radial_phi = np.asarray(rep["green"]["radial_phi"], dtype=float)
    offset = float(rep["green"]["A_offset"])
    coeff = float(rep["green"]["B_over_r"])
    axes[0, 0].loglog(radial_r[1:], np.abs(radial_phi[1:] - offset), marker="o", label="shell mean")
    axes[0, 0].loglog(radial_r[1:], np.abs(coeff) / np.maximum(radial_r[1:], 1.0e-12), "--", label="|B|/r")
    axes[0, 0].set_title(f"Green profile (n={rep['n_side']})")
    axes[0, 0].set_xlabel("r")
    axes[0, 0].set_ylabel(r"$|\phi(r)-A|$")
    axes[0, 0].grid(alpha=0.25)
    axes[0, 0].legend()

    times = np.asarray(rep["heat"]["times"], dtype=float)
    msd = np.asarray(rep["heat"]["mean_squared_radius"], dtype=float)
    slope = float(rep["heat"]["msd_slope"])
    intercept = float(rep["heat"]["msd_intercept"])
    axes[0, 1].plot(times, msd, label="shell MSD")
    axes[0, 1].plot(times, slope * times + intercept, "--", label="linear fit")
    axes[0, 1].set_title(f"Heat shell second moment (n={rep['n_side']})")
    axes[0, 1].set_xlabel("t")
    axes[0, 1].set_ylabel(r"$\langle r^2 \rangle$")
    axes[0, 1].grid(alpha=0.25)
    axes[0, 1].legend()

    shell_rows = rep["shell_arrival"]["rows"]
    radius_sq = np.asarray([row["radius_squared"] for row in shell_rows], dtype=float)
    arrivals = np.asarray([row["arrival_time"] for row in shell_rows], dtype=float)
    arr_slope = float(rep["shell_arrival"]["fit_slope"])
    arr_intercept = float(rep["shell_arrival"]["fit_intercept"])
    axes[1, 0].plot(radius_sq, arrivals, "o", label="shell arrivals")
    axes[1, 0].plot(radius_sq, arr_slope * radius_sq + arr_intercept, "--", label="linear fit")
    axes[1, 0].set_title(f"Shell-arrival geometry (n={rep['n_side']})")
    axes[1, 0].set_xlabel(r"$r^2$")
    axes[1, 0].set_ylabel("arrival time")
    axes[1, 0].grid(alpha=0.25)
    axes[1, 0].legend()

    overlap = np.asarray(rep["low_mode"]["coordinate_overlap_matrix"], dtype=float)
    im = axes[1, 1].imshow(overlap, cmap="viridis", vmin=0.0, vmax=1.0)
    axes[1, 1].set_title(f"Coordinate / first-shell overlap (n={rep['n_side']})")
    axes[1, 1].set_xticks(range(3), labels=["mode 1", "mode 2", "mode 3"])
    axes[1, 1].set_yticks(range(3), labels=["x", "y", "z"])
    fig.colorbar(im, ax=axes[1, 1], fraction=0.046, pad=0.04)

    fig.savefig(detail_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(detail_path.relative_to(REPO_ROOT)))
    return plot_paths


def save_results(result: dict[str, Any], stamp: str) -> tuple[str, str]:
    stamped = DATA / f"{stamp}_scalar_kernel_graph_geometry_bridge.json"
    latest = DATA / "scalar_kernel_graph_geometry_bridge_latest.json"
    stamped.write_text(json.dumps(result, indent=2))
    latest.write_text(json.dumps(result, indent=2))
    return str(stamped.relative_to(REPO_ROOT)), str(latest.relative_to(REPO_ROOT))


def write_note(result: dict[str, Any], result_path: str, plot_paths: list[str]) -> None:
    rows = "\n".join(
        "| "
        + " | ".join(
            [
                f"`{case['n_side']}`",
                f"`{case['green']['effective_dimension_hint']:.4f}`",
                f"`{case['green']['fit_r2']:.4f}`",
                f"`{case['heat']['effective_diffusivity']:.4f}`",
                f"`{case['heat']['fit_r2']:.4f}`",
                f"`{case['shell_arrival']['fit_slope']:.4f}`",
                f"`{case['shell_arrival']['fit_r2']:.4f}`",
                f"`{min(case['low_mode']['principal_cosines']):.4f}`",
                f"`{case['low_mode']['subspace_affinity']:.4f}`",
            ]
        )
        + " |"
        for case in result["cases"]
    )
    note = f"""# Scalar Kernel Graph Geometry Bridge v1

## Purpose

Close the preflight blocker for geometry closure by rebuilding all four geometry-facing observables on **one common scalar kernel-graph family** instead of mixing the old 3D Green-response line with the frozen 2D branch-local torus line.

The question here is deliberately bounded:

> on the same 3D cubic scalar kernel graph, under the same Euclidean shell reconstruction, do Green response, heat behavior, shell-arrival structure, and low-mode organization point to one coherent metric-like geometry read?

## Shared carrier

Use the strongest current scalar operator family:

- cubic 3D grids with `n = 9, 11, 13`
- local kernel regime `c_epsilon = 0.5`
- cutoff factor `2.5`
- kernel weights

$$
w_{{ij}} = \\exp\\left(-\\frac{{|x_i-x_j|^2}}{{\\epsilon_k}}\\right),
\\qquad
\\epsilon_k = c_\\epsilon h^2
$$

- induced scalar operator

$$
\\widehat{{L}}_h = -\\frac{{2}}{{\\mu_2 h^2}} (D-W)
$$

This is the same local scalar carrier already used for the strongest CP1 operator-control result.

## One common metric-like coarse reconstruction

All four channels use the same Euclidean coarse geometry centered at the source node:

- shell radius `r`
- shell second moment `r^2`
- coordinate subspace `span(x, y, z)`

No branch-local torus artifacts enter this bridge.

## Measured channels

1. **Green response**
   Fit the mean-zero Poisson field from `-\\widehat{{L}}_h \\phi = \\delta - N^{{-1}} 1` against `A + B/r`.

2. **Heat behavior**
   Evolve `u_t = \\widehat{{L}}_h u` from a point source and test whether the shell second moment `\\langle r^2 \\rangle(t)` is linear in `t`.

3. **Shell-arrival structure**
   On the same heat trajectory, define shell arrival as the first half-peak crossing of the shell-mean signal and fit arrival time against `r^2`.

4. **Low-mode organization**
   Compare the first positive eigen-triplet subspace of `-\\widehat{{L}}_h` against the Euclidean coordinate subspace `span(x, y, z)` using principal angles. This avoids fake instability from arbitrary rotations inside the degenerate first shell.

## Artifacts

- script: `numerics/simulations/scalar_kernel_graph_geometry_bridge.py`
- config: `config.json -> scalar_kernel_graph_geometry_bridge`
- results: `{result_path}`
- latest: `data/scalar_kernel_graph_geometry_bridge_latest.json`
- plots:
  - `{plot_paths[0]}`
  - `{plot_paths[1]}`

## Results

| `n` | Green dim hint | Green fit `R^2` | `D_eff` | Heat fit `R^2` | arrival slope | arrival fit `R^2` | min principal cosine | subspace affinity |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
{rows}

## Interpretation

- observation: {result['observation']}
- conclusion: {result['conclusion']}

## Current verdict

This bridge fixes the **category mistake** identified in the geometry-closure preflight note: all four channels now live on the same scalar operator/substrate family.

The positive read should stay bounded:

- this supports a common 3D-like metric organization on the tested local cubic scalar carrier;
- it does **not** yet claim universality across other kernel families, point-set disorder families, or full geometry emergence beyond the tested bridge regime.
"""
    NOTE_PATH.write_text(note, encoding="utf-8")


def append_log(result_path: str, plot_paths: list[str], config: dict[str, Any], observation: str, conclusion: str) -> None:
    with EXPERIMENT_LOG.open("a") as handle:
        handle.write("\n## scalar kernel graph geometry bridge\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            "- Config: "
            f"n_sides={config['n_sides']}, epsilon_coeff={config['epsilon_coeff']}, "
            f"cutoff_factor={config['cutoff_factor']}, heat_time_max={config['heat_time_max']}, "
            f"arrival_threshold_fraction={config['arrival_threshold_fraction']}\n"
        )
        handle.write(f"- Results: `{result_path}`\n")
        handle.write("- Plots: " + ", ".join(f"`{path}`" for path in plot_paths) + "\n")
        handle.write(f"- Observation: {observation}\n")
        handle.write(f"- Conclusion: {conclusion}\n")


def run_scalar_kernel_graph_geometry_bridge(config: dict[str, Any] | None = None) -> dict[str, Any]:
    cfg = load_config(config)
    cases = [
        run_case(
            n_side=int(n_side),
            epsilon_coeff=float(cfg["epsilon_coeff"]),
            cutoff_factor=float(cfg["cutoff_factor"]),
            config=cfg,
        )
        for n_side in cfg["n_sides"]
    ]
    observation, conclusion = summarize_observation(cases)
    return {
        "experiment": "scalar_kernel_graph_geometry_bridge",
        "config": cfg,
        "derivation": {
            "carrier": "one local scalar kernel-graph family on cubic 3D grids",
            "metric_like_reconstruction": "Euclidean shell radius and coordinate subspace on the same cubic carrier",
            "green_target": "A + B/r",
            "heat_target": "<r^2>(t) linear in t",
            "shell_arrival_target": "arrival time linear in r^2",
            "low_mode_target": "first positive eigen-triplet spans the Euclidean coordinate shell",
        },
        "cases": cases,
        "observation": observation,
        "conclusion": conclusion,
    }


def main() -> None:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    result = run_scalar_kernel_graph_geometry_bridge()
    plots = make_plots(result["cases"], stamp=stamp)
    result["plots"] = plots
    result_path, latest_path = save_results(result, stamp=stamp)
    write_note(result, result_path=result_path, plot_paths=plots)
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
            f"green_dim={case['green']['effective_dimension_hint']:.4f}",
            f"green_R2={case['green']['fit_r2']:.4f}",
            f"D_eff={case['heat']['effective_diffusivity']:.4f}",
            f"arrival_R2={case['shell_arrival']['fit_r2']:.4f}",
            f"low_cos={min(case['low_mode']['principal_cosines']):.4f}",
        )
    print("observation =", result["observation"])
    print("conclusion =", result["conclusion"])


if __name__ == "__main__":
    main()
