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
from scipy.spatial import cKDTree

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "data"
PLOTS = REPO_ROOT / "plots"
EXPERIMENT_LOG = REPO_ROOT / "experiments" / "EXPERIMENT_LOG.md"
NOTE_PATH = REPO_ROOT / "experiments" / "operators" / "Scalar_Kernel_Graph_Geometry_Robustness_v1.md"
for path in (DATA, PLOTS):
    path.mkdir(exist_ok=True)

MPLCONFIG = REPO_ROOT / ".mplconfig"
MPLCONFIG.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_CONFIG: dict[str, Any] = {
    "epsilon_coeff": 0.5,
    "cutoff_factor": 2.5,
    "green_fit_r_max": 0.55,
    "heat_time_max": 0.03,
    "heat_time_steps": 61,
    "msd_fit_window": [0.004, 0.02],
    "arrival_threshold_fraction": 0.5,
    "arrival_scan_shells": 9,
    "low_mode_count": 8,
    "disorder_n_sides": [11, 13],
    "jitter_fractions": [0.0, 0.02, 0.04],
    "jitter_seed": 42,
    "kernel_n_sides": [11, 13],
    "kernel_families": ["gaussian_local", "gaussian_half", "inverse_quadratic"],
    "thresholds": {
        "min_green_fit_r2": 0.97,
        "max_green_dimension_deviation": 0.06,
        "min_heat_fit_r2": 0.995,
        "min_arrival_fit_r2": 0.97,
        "min_low_mode_cosine": 0.99,
        "max_metric_anisotropy": 0.02,
        "max_metric_offdiag_ratio": 0.02,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = json.loads(json.dumps(DEFAULT_CONFIG))
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        raw = json.loads(path.read_text())
        if isinstance(raw.get("scalar_kernel_graph_geometry_robustness"), dict):
            merged.update({k: v for k, v in raw["scalar_kernel_graph_geometry_robustness"].items() if k != "thresholds"})
            if isinstance(raw["scalar_kernel_graph_geometry_robustness"].get("thresholds"), dict):
                merged["thresholds"].update(raw["scalar_kernel_graph_geometry_robustness"]["thresholds"])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != "thresholds"})
        if isinstance(config.get("thresholds"), dict):
            merged["thresholds"].update(config["thresholds"])
    return merged


def build_cubic_points(n_side: int) -> tuple[np.ndarray, int, float]:
    grid = np.linspace(0.0, 1.0, n_side)
    x, y, z = np.meshgrid(grid, grid, grid, indexing="ij")
    points = np.column_stack([x.ravel(), y.ravel(), z.ravel()]).astype(float)
    source_idx = int(((n_side // 2) * n_side + (n_side // 2)) * n_side + (n_side // 2))
    h = float(grid[1] - grid[0])
    return points, source_idx, h


def build_jittered_points(n_side: int, jitter_fraction: float, seed: int) -> tuple[np.ndarray, int, float]:
    points, source_idx, h = build_cubic_points(n_side)
    if jitter_fraction <= 0.0:
        return points, source_idx, h
    grid_ids = np.indices((n_side, n_side, n_side)).reshape(3, -1).T
    boundary_mask = np.any((grid_ids == 0) | (grid_ids == n_side - 1), axis=1)
    rng = np.random.default_rng(seed)
    jitter = jitter_fraction * h * rng.standard_normal(points.shape)
    jitter[boundary_mask] = 0.0
    points = np.clip(points + jitter, 0.0, 1.0)
    return points, source_idx, h


def kernel_weight(q: float, family: str) -> float:
    if family == "gaussian_local":
        return math.exp(-q)
    if family == "gaussian_half":
        return math.exp(-0.5 * q)
    if family == "inverse_quadratic":
        return 1.0 / (1.0 + q) ** 2
    raise ValueError(f"unknown kernel family: {family}")


def build_point_cloud_operator(
    points: np.ndarray,
    h: float,
    epsilon_coeff: float,
    cutoff_factor: float,
    kernel_family: str,
) -> tuple[coo_matrix, float]:
    epsilon_k = epsilon_coeff * h * h
    cutoff = cutoff_factor * math.sqrt(epsilon_k)
    tree = cKDTree(points)
    neighborhoods = tree.query_ball_point(points, cutoff * (1.0 + 1.0e-12))

    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    degrees = np.zeros(points.shape[0], dtype=float)
    local_moments: list[float] = []

    for index, raw_neighbors in enumerate(neighborhoods):
        local_second_moment = 0.0
        for neighbor in raw_neighbors:
            if neighbor == index:
                continue
            delta = points[neighbor] - points[index]
            d2 = float(np.dot(delta, delta))
            q = d2 / max(epsilon_k, 1.0e-18)
            weight = kernel_weight(q, kernel_family)
            rows.append(index)
            cols.append(neighbor)
            data.append(-weight)
            degrees[index] += weight
            local_second_moment += weight * float(np.dot(delta, delta)) / 3.0
        local_moments.append(local_second_moment)

    diagonal = np.arange(points.shape[0], dtype=int)
    rows.extend(diagonal.tolist())
    cols.extend(diagonal.tolist())
    data.extend(degrees.tolist())
    laplacian = coo_matrix((data, (rows, cols)), shape=(points.shape[0], points.shape[0])).tocsr()
    mean_second_moment = float(np.mean(local_moments))
    scale = -2.0 / max(mean_second_moment, 1.0e-18)
    return laplacian * scale, epsilon_k


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


def green_metrics(points: np.ndarray, source_idx: int, h: float, epsilon_k: float, phi: np.ndarray, fit_r_max: float) -> dict[str, Any]:
    radii = np.linalg.norm(points - points[source_idx], axis=1)
    fit_r_min = max(2.0 * h, 1.5 * math.sqrt(epsilon_k))
    mask = (radii >= fit_r_min) & (radii <= fit_r_max)
    x = np.column_stack([1.0 / radii[mask], np.ones(np.sum(mask), dtype=float)])
    coeffs, *_ = np.linalg.lstsq(x, phi[mask], rcond=None)
    pred = x @ coeffs
    ss_res = float(np.sum((phi[mask] - pred) ** 2))
    ss_tot = float(np.sum((phi[mask] - np.mean(phi[mask])) ** 2))
    fit = fit_power(radii[mask], np.abs(phi[mask] - float(coeffs[1])))
    return {
        "fit_r2": float(1.0 - ss_res / max(ss_tot, 1.0e-18)),
        "residual_slope": float(fit["slope"]),
        "effective_dimension_hint": float(2.0 - fit["slope"]),
        "A_offset": float(coeffs[1]),
        "B_over_r": float(coeffs[0]),
    }


def shell_assignment(points: np.ndarray, source_idx: int, h: float, mode: str) -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    radii = np.linalg.norm(points - points[source_idx], axis=1)
    if mode == "rounded":
        labels = np.round(radii / max(h, 1.0e-12), 8) * h
        shell_radii = np.unique(labels)
        shell_ids = np.argmin(np.abs(labels[:, None] - shell_radii[None, :]), axis=1)
    elif mode == "ideal_clean_reference":
        n_side = int(round(1.0 / h)) + 1
        clean_points, clean_source_idx, _ = build_cubic_points(n_side)
        clean_radii = np.linalg.norm(clean_points - clean_points[clean_source_idx], axis=1)
        shell_radii = np.unique(np.round(clean_radii / max(h, 1.0e-12), 8) * h)
        shell_ids = np.argmin(np.abs(radii[:, None] - shell_radii[None, :]), axis=1)
    else:
        raise ValueError(f"unknown shell mode: {mode}")
    shell_indices = [np.flatnonzero(shell_ids == index) for index in range(len(shell_radii))]
    shell_means = np.asarray([float(np.mean(radii[idx])) for idx in shell_indices], dtype=float)
    return shell_means, radii, shell_indices


def heat_trajectory(delta: coo_matrix, source_idx: int, time_max: float, time_steps: int) -> tuple[np.ndarray, np.ndarray]:
    initial = np.zeros(delta.shape[0], dtype=float)
    initial[source_idx] = 1.0
    times = np.linspace(0.0, time_max, time_steps)
    trajectory = expm_multiply(delta, initial, start=times[0], stop=times[-1], num=len(times), endpoint=True)
    return times, np.asarray(trajectory, dtype=float)


def heat_metrics(points: np.ndarray, source_idx: int, times: np.ndarray, trajectory: np.ndarray, fit_window: list[float]) -> dict[str, Any]:
    radii_sq = np.sum((points - points[source_idx]) ** 2, axis=1)
    msd = trajectory @ radii_sq
    mask = (times >= float(fit_window[0])) & (times <= float(fit_window[1]))
    slope, intercept = np.polyfit(times[mask], msd[mask], 1)
    pred = slope * times[mask] + intercept
    ss_res = float(np.sum((msd[mask] - pred) ** 2))
    ss_tot = float(np.sum((msd[mask] - np.mean(msd[mask])) ** 2))
    return {
        "msd_slope": float(slope),
        "effective_diffusivity": float(slope / 6.0),
        "fit_r2": float(1.0 - ss_res / max(ss_tot, 1.0e-18)),
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
    for shell_idx in range(1, min(scan_shells + 1, len(shell_indices))):
        if shell_indices[shell_idx].size == 0:
            continue
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
    coordinate_basis, _ = np.linalg.qr(points - np.mean(points, axis=0))
    mode_basis, _ = np.linalg.qr(eigenvectors[:, 1:4])
    overlap = coordinate_basis.T @ mode_basis
    principal_cosines = np.linalg.svd(overlap, compute_uv=False)
    return {
        "eigenvalues": eigenvalues.tolist(),
        "principal_cosines": principal_cosines.tolist(),
        "subspace_affinity": float(np.mean(principal_cosines**2)),
    }


def metric_response_tensor(points: np.ndarray, delta: coo_matrix) -> dict[str, Any]:
    coordinate_basis, _ = np.linalg.qr(points - np.mean(points, axis=0))
    tensor = np.asarray(coordinate_basis.T @ ((-delta) @ coordinate_basis), dtype=float)
    diagonal = np.diag(tensor)
    offdiag = tensor - np.diag(diagonal)
    mean_diag = max(float(np.mean(diagonal)), 1.0e-18)
    return {
        "tensor": tensor.tolist(),
        "diagonal": diagonal.tolist(),
        "anisotropy": float((np.max(diagonal) - np.min(diagonal)) / mean_diag),
        "offdiag_ratio": float(np.linalg.norm(offdiag) / mean_diag),
    }


def evaluate_case(
    label: str,
    points: np.ndarray,
    source_idx: int,
    h: float,
    epsilon_coeff: float,
    cutoff_factor: float,
    kernel_family: str,
    shell_mode: str,
    config: dict[str, Any],
) -> dict[str, Any]:
    delta, epsilon_k = build_point_cloud_operator(
        points=points,
        h=h,
        epsilon_coeff=epsilon_coeff,
        cutoff_factor=cutoff_factor,
        kernel_family=kernel_family,
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
    shell_radii, _, shell_indices = shell_assignment(points=points, source_idx=source_idx, h=h, mode=shell_mode)
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
    metric_response = metric_response_tensor(points=points, delta=delta)
    return {
        "label": label,
        "points_count": int(points.shape[0]),
        "lattice_spacing": float(h),
        "epsilon_coeff": float(epsilon_coeff),
        "kernel_family": kernel_family,
        "shell_mode": shell_mode,
        "green": green,
        "heat": heat,
        "shell_arrival": shell_arrival,
        "low_mode": low_mode,
        "metric_response": metric_response,
    }


def case_passes(case: dict[str, Any], thresholds: dict[str, float]) -> bool:
    return bool(
        case["green"]["fit_r2"] >= float(thresholds["min_green_fit_r2"])
        and abs(case["green"]["effective_dimension_hint"] - 3.0) <= float(thresholds["max_green_dimension_deviation"])
        and case["heat"]["fit_r2"] >= float(thresholds["min_heat_fit_r2"])
        and case["shell_arrival"]["fit_r2"] >= float(thresholds["min_arrival_fit_r2"])
        and min(case["low_mode"]["principal_cosines"]) >= float(thresholds["min_low_mode_cosine"])
        and case["metric_response"]["anisotropy"] <= float(thresholds["max_metric_anisotropy"])
        and case["metric_response"]["offdiag_ratio"] <= float(thresholds["max_metric_offdiag_ratio"])
    )


def make_result(config: dict[str, Any] | None = None) -> dict[str, Any]:
    cfg = load_config(config)
    thresholds = cfg["thresholds"]

    disorder_cases: list[dict[str, Any]] = []
    for n_side in [int(v) for v in cfg["disorder_n_sides"]]:
        for jitter_fraction in [float(v) for v in cfg["jitter_fractions"]]:
            points, source_idx, h = build_jittered_points(
                n_side=n_side,
                jitter_fraction=jitter_fraction,
                seed=int(cfg["jitter_seed"]),
            )
            case = evaluate_case(
                label=f"disorder|n={n_side}|jitter={jitter_fraction:.3f}",
                points=points,
                source_idx=source_idx,
                h=h,
                epsilon_coeff=float(cfg["epsilon_coeff"]),
                cutoff_factor=float(cfg["cutoff_factor"]),
                kernel_family="gaussian_local",
                shell_mode="rounded" if jitter_fraction == 0.0 else "ideal_clean_reference",
                config=cfg,
            )
            case["jitter_fraction"] = float(jitter_fraction)
            case["pass"] = case_passes(case, thresholds)
            disorder_cases.append(case)

    kernel_cases: list[dict[str, Any]] = []
    for n_side in [int(v) for v in cfg["kernel_n_sides"]]:
        points, source_idx, h = build_cubic_points(n_side)
        for family in [str(v) for v in cfg["kernel_families"]]:
            case = evaluate_case(
                label=f"kernel|n={n_side}|family={family}",
                points=points,
                source_idx=source_idx,
                h=h,
                epsilon_coeff=float(cfg["epsilon_coeff"]),
                cutoff_factor=float(cfg["cutoff_factor"]),
                kernel_family=family,
                shell_mode="rounded",
                config=cfg,
            )
            case["pass"] = case_passes(case, thresholds)
            kernel_cases.append(case)

    passing_disorder = sum(1 for case in disorder_cases if case["pass"])
    passing_kernels = sum(1 for case in kernel_cases if case["pass"])
    weakest_disorder = min(disorder_cases, key=lambda item: item["shell_arrival"]["fit_r2"])
    weakest_kernel = min(kernel_cases, key=lambda item: item["shell_arrival"]["fit_r2"])

    if passing_disorder == len(disorder_cases) and passing_kernels == len(kernel_cases):
        observation = (
            "the scalar kernel-graph geometry bridge survives the first bounded robustness pass: mild point-set disorder, bounded kernel-family variation, "
            "and coordinate-response extraction all stay compatible with the same 3D-like metric read on the tested carrier"
        )
        conclusion = (
            f"all {len(disorder_cases)} disorder cases and all {len(kernel_cases)} kernel-family cases pass the shared CP4 thresholds; "
            f"the weakest disorder case is `{weakest_disorder['label']}` with shell-arrival fit R^2 = {weakest_disorder['shell_arrival']['fit_r2']:.4f}, "
            f"and the weakest kernel case is `{weakest_kernel['label']}` with shell-arrival fit R^2 = {weakest_kernel['shell_arrival']['fit_r2']:.4f}; "
            "this means the scalar-carrier geometry closure is now robust across the tested mild disorder and bounded kernel-family window, while remaining a bounded carrier-level statement rather than full universality"
        )
    else:
        failing_disorder = [case["label"] for case in disorder_cases if not case["pass"]]
        failing_kernels = [case["label"] for case in kernel_cases if not case["pass"]]
        observation = (
            "the scalar kernel-graph geometry bridge remains coherent at baseline, but the first robustness pass exposes cases where the same geometry read does not survive unchanged"
        )
        conclusion = (
            f"open robustness cases remain: disorder={failing_disorder}, kernels={failing_kernels}; "
            "the correct boundary therefore stays at the original bounded scalar-carrier closure without a positive robustness extension"
        )

    return {
        "experiment": "scalar_kernel_graph_geometry_robustness",
        "config": cfg,
        "disorder_cases": disorder_cases,
        "kernel_cases": kernel_cases,
        "observation": observation,
        "conclusion": conclusion,
    }


def make_disorder_plot(cases: list[dict[str, Any]], stamp: str) -> str:
    plot_path = PLOTS / f"{stamp}_scalar_kernel_graph_geometry_robustness_disorder.png"
    n_values = sorted({int(case["label"].split("|")[1].split("=")[1]) for case in cases})
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    metrics = [
        ("green.fit_r2", "Green fit $R^2$", axes[0, 0]),
        ("shell_arrival.fit_r2", "Shell-arrival fit $R^2$", axes[0, 1]),
        ("low_mode.min_cos", "Min low-mode cosine", axes[1, 0]),
        ("metric_response.anisotropy", "Metric anisotropy", axes[1, 1]),
    ]
    for n_side in n_values:
        subset = [case for case in cases if int(case["label"].split("|")[1].split("=")[1]) == n_side]
        subset = sorted(subset, key=lambda item: float(item["jitter_fraction"]))
        x = [float(case["jitter_fraction"]) for case in subset]
        values = {
            "green.fit_r2": [case["green"]["fit_r2"] for case in subset],
            "shell_arrival.fit_r2": [case["shell_arrival"]["fit_r2"] for case in subset],
            "low_mode.min_cos": [min(case["low_mode"]["principal_cosines"]) for case in subset],
            "metric_response.anisotropy": [case["metric_response"]["anisotropy"] for case in subset],
        }
        for key, title, axis in metrics:
            axis.plot(x, values[key], marker="o", label=f"n={n_side}")
            axis.set_title(title)
            axis.set_xlabel("jitter fraction")
            axis.grid(alpha=0.25)
    axes[0, 0].legend()
    axes[0, 1].legend()
    axes[1, 0].legend()
    axes[1, 1].legend()
    fig.savefig(plot_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return str(plot_path.relative_to(REPO_ROOT))


def make_kernel_plot(cases: list[dict[str, Any]], stamp: str) -> str:
    plot_path = PLOTS / f"{stamp}_scalar_kernel_graph_geometry_robustness_kernels.png"
    families = [str(case["kernel_family"]) for case in cases if int(case["label"].split("|")[1].split("=")[1]) == min(int(item["label"].split("|")[1].split("=")[1]) for item in cases)]
    family_order = [str(v) for v in DEFAULT_CONFIG["kernel_families"]]
    n_values = sorted({int(case["label"].split("|")[1].split("=")[1]) for case in cases})
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    metrics = [
        ("green.fit_r2", "Green fit $R^2$", axes[0, 0]),
        ("shell_arrival.fit_r2", "Shell-arrival fit $R^2$", axes[0, 1]),
        ("low_mode.min_cos", "Min low-mode cosine", axes[1, 0]),
        ("metric_response.offdiag", "Metric offdiag ratio", axes[1, 1]),
    ]
    for n_side in n_values:
        subset = [case for case in cases if int(case["label"].split("|")[1].split("=")[1]) == n_side]
        subset = sorted(subset, key=lambda item: family_order.index(str(item["kernel_family"])))
        x = range(len(subset))
        values = {
            "green.fit_r2": [case["green"]["fit_r2"] for case in subset],
            "shell_arrival.fit_r2": [case["shell_arrival"]["fit_r2"] for case in subset],
            "low_mode.min_cos": [min(case["low_mode"]["principal_cosines"]) for case in subset],
            "metric_response.offdiag": [case["metric_response"]["offdiag_ratio"] for case in subset],
        }
        labels = [str(case["kernel_family"]) for case in subset]
        for key, title, axis in metrics:
            axis.plot(x, values[key], marker="o", label=f"n={n_side}")
            axis.set_title(title)
            axis.set_xticks(list(x), labels=labels, rotation=20)
            axis.grid(alpha=0.25)
    axes[0, 0].legend()
    axes[0, 1].legend()
    axes[1, 0].legend()
    axes[1, 1].legend()
    fig.savefig(plot_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return str(plot_path.relative_to(REPO_ROOT))


def make_metric_tensor_plot(disorder_cases: list[dict[str, Any]], kernel_cases: list[dict[str, Any]], stamp: str) -> str:
    plot_path = PLOTS / f"{stamp}_scalar_kernel_graph_geometry_robustness_metric_tensor.png"
    selected = [
        min(disorder_cases, key=lambda item: abs(float(item["jitter_fraction"]) - 0.04)),
        [case for case in kernel_cases if str(case["kernel_family"]) == "gaussian_local"][-1],
        [case for case in kernel_cases if str(case["kernel_family"]) == "gaussian_half"][-1],
        [case for case in kernel_cases if str(case["kernel_family"]) == "inverse_quadratic"][-1],
    ]
    fig, axes = plt.subplots(2, 2, figsize=(9, 8))
    for axis, case in zip(axes.ravel(), selected):
        tensor = np.asarray(case["metric_response"]["tensor"], dtype=float)
        mean_diag = np.mean(np.diag(tensor))
        im = axis.imshow(tensor / max(mean_diag, 1.0e-18), cmap="viridis", vmin=0.0, vmax=1.1)
        axis.set_title(case["label"].replace("|", "\n"))
        axis.set_xticks(range(3), labels=["x", "y", "z"])
        axis.set_yticks(range(3), labels=["x", "y", "z"])
        fig.colorbar(im, ax=axis, fraction=0.046, pad=0.04)
    fig.savefig(plot_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return str(plot_path.relative_to(REPO_ROOT))


def save_results(result: dict[str, Any], stamp: str) -> tuple[str, str]:
    stamped = DATA / f"{stamp}_scalar_kernel_graph_geometry_robustness.json"
    latest = DATA / "scalar_kernel_graph_geometry_robustness_latest.json"
    stamped.write_text(json.dumps(result, indent=2))
    latest.write_text(json.dumps(result, indent=2))
    return str(stamped.relative_to(REPO_ROOT)), str(latest.relative_to(REPO_ROOT))


def write_note(result: dict[str, Any], result_path: str, plot_paths: list[str]) -> None:
    disorder_rows = "\n".join(
        f"| `{case['label']}` | `{case['green']['fit_r2']:.4f}` | `{case['green']['effective_dimension_hint']:.4f}` | `{case['heat']['fit_r2']:.4f}` | `{case['shell_arrival']['fit_r2']:.4f}` | `{min(case['low_mode']['principal_cosines']):.4f}` | `{case['metric_response']['anisotropy']:.4f}` | `{case['metric_response']['offdiag_ratio']:.4f}` | {'PASS' if case['pass'] else 'OPEN'} |"
        for case in result["disorder_cases"]
    )
    kernel_rows = "\n".join(
        f"| `{case['label']}` | `{case['green']['fit_r2']:.4f}` | `{case['green']['effective_dimension_hint']:.4f}` | `{case['heat']['fit_r2']:.4f}` | `{case['shell_arrival']['fit_r2']:.4f}` | `{min(case['low_mode']['principal_cosines']):.4f}` | `{case['metric_response']['anisotropy']:.4f}` | `{case['metric_response']['offdiag_ratio']:.4f}` | {'PASS' if case['pass'] else 'OPEN'} |"
        for case in result["kernel_cases"]
    )
    note = f"""# Scalar Kernel Graph Geometry Robustness v1

## Purpose

Run the next honest move after the bounded scalar-carrier `CP4` bridge:

1. mild point-set disorder on the same carrier
2. bounded kernel-family variation on the same carrier
3. response / metric extraction on that same coarse reconstruction

This note does **not** introduce a new geometry story. It tests whether the existing scalar-carrier geometry bridge survives these bounded stresses.

## Construction

The carrier remains the same:

- scalar kernel graph on `3D` cubic point clouds
- local regime `c_epsilon = 0.5`
- same Green / heat / shell-arrival / first-shell low-mode observables
- same Euclidean reconstruction centered at the source

Two bounded extensions are tested.

### 1. Disorder pass

- `n = {result['config']['disorder_n_sides']}`
- jitter fractions `= {result['config']['jitter_fractions']}`
- kernel family fixed to `gaussian_local`

For jittered cases, shell-arrival bins are assigned to the nearest clean-grid shell at the same `n`, so the coarse reconstruction itself stays fixed while the point cloud is perturbed.

### 2. Kernel-family pass

- `n = {result['config']['kernel_n_sides']}`
- families `= {result['config']['kernel_families']}`

Families tested:

- `gaussian_local`: `exp(-|x_i-x_j|^2 / epsilon_k)`
- `gaussian_half`: `exp(-|x_i-x_j|^2 / (2 epsilon_k))`
- `inverse_quadratic`: `(1 + |x_i-x_j|^2 / epsilon_k)^(-2)`

### 3. Response / metric extraction

On every case, extract the coordinate-response stiffness tensor on the orthonormalized coordinate basis `span(x,y,z)`:

$$
K_{{coord}} = Q^T (-\\widehat{{L}}_h) Q.
$$

This is the bounded metric-response read for the same carrier. In the isotropic case it should stay close to a scalar multiple of the identity.

## Artifacts

- script: `numerics/simulations/scalar_kernel_graph_geometry_robustness.py`
- config: `config.json -> scalar_kernel_graph_geometry_robustness`
- results: `{result_path}`
- latest: `data/scalar_kernel_graph_geometry_robustness_latest.json`
- plots:
  - `{plot_paths[0]}`
  - `{plot_paths[1]}`
  - `{plot_paths[2]}`

## Disorder pass

| case | Green fit `R^2` | Green dim hint | Heat fit `R^2` | shell-arrival fit `R^2` | min low-mode cosine | metric anisotropy | metric offdiag ratio | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
{disorder_rows}

## Kernel-family pass

| case | Green fit `R^2` | Green dim hint | Heat fit `R^2` | shell-arrival fit `R^2` | min low-mode cosine | metric anisotropy | metric offdiag ratio | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
{kernel_rows}

## Interpretation

- observation: {result['observation']}
- conclusion: {result['conclusion']}

## Current boundary

If this note is read positively, the correct bounded statement is:

> the scalar-carrier geometry closure now survives the first mild disorder and bounded kernel-family robustness pass, while the same coordinate-response tensor remains close to isotropic on the tested window.

The remaining open step is still stronger than this:

- broader disorder families
- broader kernel families
- stronger response / metric extraction beyond the coordinate stiffness tensor
"""
    NOTE_PATH.write_text(note, encoding="utf-8")


def append_log(result_path: str, plot_paths: list[str], config: dict[str, Any], observation: str, conclusion: str) -> None:
    with EXPERIMENT_LOG.open("a") as handle:
        handle.write("\n## scalar kernel graph geometry robustness\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            "- Config: "
            f"disorder_n_sides={config['disorder_n_sides']}, jitter_fractions={config['jitter_fractions']}, "
            f"kernel_n_sides={config['kernel_n_sides']}, kernel_families={config['kernel_families']}\n"
        )
        handle.write(f"- Results: `{result_path}`\n")
        handle.write("- Plots: " + ", ".join(f"`{path}`" for path in plot_paths) + "\n")
        handle.write(f"- Observation: {observation}\n")
        handle.write(f"- Conclusion: {conclusion}\n")


def main() -> None:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    result = make_result()
    plots = [
        make_disorder_plot(result["disorder_cases"], stamp=stamp),
        make_kernel_plot(result["kernel_cases"], stamp=stamp),
        make_metric_tensor_plot(result["disorder_cases"], result["kernel_cases"], stamp=stamp),
    ]
    result["plots"] = plots
    result_path, latest_path = save_results(result, stamp=stamp)
    write_note(result, result_path=result_path, plot_paths=plots)
    append_log(result_path, plots, result["config"], result["observation"], result["conclusion"])
    print("results =", result_path)
    print("latest =", latest_path)
    print("plots =", plots)
    print("observation =", result["observation"])
    print("conclusion =", result["conclusion"])


if __name__ == "__main__":
    main()
