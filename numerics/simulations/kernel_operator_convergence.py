#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import os
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

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
    "point_set_n_sides": [11, 13, 15, 17],
    "point_set_epsilon_coeffs": [0.5],
    "point_set_jitter_fraction": 0.04,
    "point_set_seed": 42,
    "point_set_boundary_buffer_factor": 1.05,
    "boundary_n_sides": [9, 11, 13, 17, 21],
    "boundary_epsilon_coeffs": [0.5, 1.0, 2.0],
}

QUADRATIC_LAPLACIAN = 6.0
POINT_SET_BASIS_SIZE = 9
LAPLACIAN_SELECTOR = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0], dtype=float)


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
    lap1 = np.full_like(f1, QUADRATIC_LAPLACIAN)

    f2 = np.sin(math.pi * x) * np.sin(math.pi * y) * np.sin(math.pi * z)
    lap2 = -3.0 * math.pi * math.pi * f2

    f3 = np.exp(-r2)
    lap3 = (4.0 * r2 - 6.0) * f3

    return {
        "f1_quadratic": {"values": f1, "laplacian": lap1},
        "f2_trigonometric": {"values": f2, "laplacian": lap2},
        "f3_gaussian": {"values": f3, "laplacian": lap3},
    }


def boundary_functions(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> dict[str, dict[str, np.ndarray | str]]:
    f_dirichlet = np.sin(math.pi * x) * np.sin(math.pi * y) * np.sin(math.pi * z)
    lap_dirichlet = -3.0 * math.pi * math.pi * f_dirichlet

    f_neumann = np.cos(math.pi * x) * np.cos(math.pi * y) * np.cos(math.pi * z)
    lap_neumann = -3.0 * math.pi * math.pi * f_neumann

    return {
        "dirichlet_trigonometric": {
            "values": f_dirichlet,
            "laplacian": lap_dirichlet,
            "boundary_mode": "dirichlet",
        },
        "neumann_trigonometric": {
            "values": f_neumann,
            "laplacian": lap_neumann,
            "boundary_mode": "neumann",
        },
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


def reflect_index(index: int, n: int, boundary_mode: str) -> tuple[int, float]:
    reflected = index
    sign = 1.0
    while reflected < 0 or reflected >= n:
        if reflected < 0:
            reflected = -reflected
        else:
            reflected = 2 * (n - 1) - reflected
        if boundary_mode == "dirichlet":
            sign *= -1.0
    return reflected, sign


def apply_graph_laplacian_reflected(
    values: np.ndarray,
    offsets: list[tuple[int, int, int, float]],
    boundary_mode: str,
) -> np.ndarray:
    out = np.zeros_like(values)
    n = values.shape[0]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                center = float(values[i, j, k])
                accum = 0.0
                for dx, dy, dz, weight in offsets:
                    ni, sx = reflect_index(i + dx, n, boundary_mode)
                    nj, sy = reflect_index(j + dy, n, boundary_mode)
                    nk, sz = reflect_index(k + dz, n, boundary_mode)
                    accum += weight * (center - float(sx * sy * sz * values[ni, nj, nk]))
                out[i, j, k] = accum
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


def full_domain_l2_error(discrete: np.ndarray, analytic: np.ndarray, h: float) -> float:
    diff = discrete - analytic
    return float(math.sqrt(np.sum(diff * diff) * (h**3)))


def profile_line(values: np.ndarray) -> np.ndarray:
    center_y = values.shape[1] // 2
    center_z = values.shape[2] // 2
    return np.asarray(values[:, center_y, center_z], dtype=float)


def fit_error_scaling(h_values: np.ndarray, errors: np.ndarray) -> dict[str, float]:
    safe_errors = np.maximum(errors, 1.0e-18)
    x = np.log(h_values)
    y = np.log(safe_errors)
    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * x + intercept
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / max(ss_tot, 1.0e-18)
    return {"order_p": float(slope), "intercept": float(intercept), "r2": float(r2)}


def run_baseline_case(n_side: int, epsilon_coeff: float, cutoff_factor: float) -> dict[str, Any]:
    x, y, z, h = build_grid(n_side)
    funcs = analytic_functions(x, y, z)
    offsets = build_offsets(epsilon_coeff=epsilon_coeff, cutoff_factor=cutoff_factor)
    scale = induced_operator_scale(h=h, offsets=offsets)
    margin = interior_margin(offsets)

    results: dict[str, Any] = {}
    for name, payload in funcs.items():
        values = np.asarray(payload["values"], dtype=float)
        analytic_lap = np.asarray(payload["laplacian"], dtype=float)
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


def run_boundary_case(n_side: int, epsilon_coeff: float, cutoff_factor: float) -> dict[str, Any]:
    x, y, z, h = build_grid(n_side)
    funcs = boundary_functions(x, y, z)
    offsets = build_offsets(epsilon_coeff=epsilon_coeff, cutoff_factor=cutoff_factor)
    scale = induced_operator_scale(h=h, offsets=offsets)

    results: dict[str, Any] = {}
    for name, payload in funcs.items():
        values = np.asarray(payload["values"], dtype=float)
        analytic_lap = np.asarray(payload["laplacian"], dtype=float)
        boundary_mode = str(payload["boundary_mode"])
        raw = apply_graph_laplacian_reflected(values, offsets, boundary_mode=boundary_mode)
        induced = scale * raw
        error = full_domain_l2_error(induced, analytic_lap, h=h)
        results[name] = {
            "error_l2": error,
            "boundary_mode": boundary_mode,
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
        "offset_count": len(offsets),
        "functions": results,
    }


def build_jittered_points(n_side: int, jitter_fraction: float, seed: int) -> tuple[np.ndarray, float]:
    x, y, z, h = build_grid(n_side)
    points = np.column_stack([x.ravel(), y.ravel(), z.ravel()]).astype(float)
    grid_ids = np.indices((n_side, n_side, n_side)).reshape(3, -1).T
    boundary_mask = np.any((grid_ids == 0) | (grid_ids == n_side - 1), axis=1)
    rng = np.random.default_rng(seed)
    jitter = jitter_fraction * h * rng.standard_normal(points.shape)
    jitter[boundary_mask] = 0.0
    points = np.clip(points + jitter, 0.0, 1.0)
    return points, h


def point_set_basis(deltas: np.ndarray) -> np.ndarray:
    dx = deltas[:, 0]
    dy = deltas[:, 1]
    dz = deltas[:, 2]
    return np.column_stack(
        [
            dx,
            dy,
            dz,
            0.5 * dx * dx,
            0.5 * dy * dy,
            0.5 * dz * dz,
            dx * dy,
            dx * dz,
            dy * dz,
        ]
    )


def build_point_set_operator(points: np.ndarray, epsilon_k: float, cutoff_factor: float) -> dict[str, Any]:
    cutoff = cutoff_factor * math.sqrt(epsilon_k)
    tree = cKDTree(points)
    neighborhoods = tree.query_ball_point(points, cutoff * (1.0 + 1.0e-12))

    stencils: list[dict[str, Any]] = []
    neighbor_counts: list[int] = []
    condition_numbers: list[float] = []
    valid_count = 0

    for index, raw_neighbors in enumerate(neighborhoods):
        neighbors = np.asarray([value for value in raw_neighbors if value != index], dtype=int)
        neighbor_counts.append(int(neighbors.size))
        if neighbors.size < POINT_SET_BASIS_SIZE:
            condition_numbers.append(float("inf"))
            stencils.append({"valid": False, "neighbors": neighbors, "lap_weights": np.zeros(0, dtype=float)})
            continue

        deltas = points[neighbors] - points[index]
        d2 = np.sum(deltas * deltas, axis=1)
        weights = np.exp(-d2 / max(epsilon_k, 1.0e-18))
        basis = point_set_basis(deltas)
        weighted_basis = np.sqrt(weights)[:, None] * basis
        singular_values = np.linalg.svd(weighted_basis, compute_uv=False)
        if singular_values[-1] <= 1.0e-14:
            condition = float("inf")
            rank = 0
        else:
            condition = float(singular_values[0] / singular_values[-1])
            rank = int(np.sum(singular_values > singular_values[0] * 1.0e-12))
        condition_numbers.append(condition)
        if rank < POINT_SET_BASIS_SIZE:
            stencils.append({"valid": False, "neighbors": neighbors, "lap_weights": np.zeros(0, dtype=float)})
            continue

        normal = basis.T @ (weights[:, None] * basis)
        rhs = basis.T * weights
        try:
            projection = np.linalg.solve(normal, rhs)
        except np.linalg.LinAlgError:
            reg_scale = max(float(np.trace(normal)) / POINT_SET_BASIS_SIZE, 1.0e-14)
            projection = np.linalg.solve(normal + 1.0e-12 * reg_scale * np.eye(POINT_SET_BASIS_SIZE), rhs)
        lap_weights = np.asarray(LAPLACIAN_SELECTOR @ projection, dtype=float)
        stencils.append({"valid": True, "neighbors": neighbors, "lap_weights": lap_weights})
        valid_count += 1

    finite_conditions = [value for value in condition_numbers if math.isfinite(value)]
    return {
        "cutoff_radius": cutoff,
        "stencils": stencils,
        "average_neighbor_count": float(np.mean(neighbor_counts)) if neighbor_counts else 0.0,
        "min_neighbor_count": int(min(neighbor_counts)) if neighbor_counts else 0,
        "median_condition_number": float(np.median(finite_conditions)) if finite_conditions else float("inf"),
        "max_condition_number": float(max(finite_conditions)) if finite_conditions else float("inf"),
        "valid_stencil_fraction": float(valid_count / max(len(stencils), 1)),
    }


def apply_point_set_operator(values: np.ndarray, operator: dict[str, Any]) -> np.ndarray:
    output = np.full(values.shape[0], np.nan, dtype=float)
    for index, stencil in enumerate(operator["stencils"]):
        if not bool(stencil["valid"]):
            continue
        neighbors = np.asarray(stencil["neighbors"], dtype=int)
        if neighbors.size == 0:
            continue
        local_diff = values[neighbors] - values[index]
        output[index] = float(np.dot(np.asarray(stencil["lap_weights"], dtype=float), local_diff))
    return output


def point_set_l2_error(discrete: np.ndarray, analytic: np.ndarray, valid_mask: np.ndarray, h: float) -> float:
    diff = discrete[valid_mask] - analytic[valid_mask]
    return float(math.sqrt(np.sum(diff * diff) * (h**3)))


def run_point_set_case(
    n_side: int,
    epsilon_coeff: float,
    cutoff_factor: float,
    jitter_fraction: float,
    seed: int,
    boundary_buffer_factor: float,
) -> dict[str, Any]:
    points, h = build_jittered_points(n_side=n_side, jitter_fraction=jitter_fraction, seed=seed)
    epsilon_k = epsilon_coeff * h * h
    operator = build_point_set_operator(points=points, epsilon_k=epsilon_k, cutoff_factor=cutoff_factor)
    funcs = analytic_functions(points[:, 0], points[:, 1], points[:, 2])
    boundary_distance = np.min(
        np.column_stack(
            [
                points[:, 0],
                1.0 - points[:, 0],
                points[:, 1],
                1.0 - points[:, 1],
                points[:, 2],
                1.0 - points[:, 2],
            ]
        ),
        axis=1,
    )
    stencil_valid = np.array([bool(item["valid"]) for item in operator["stencils"]], dtype=bool)
    valid_mask = stencil_valid & (boundary_distance >= boundary_buffer_factor * float(operator["cutoff_radius"]))

    results: dict[str, Any] = {}
    for name, payload in funcs.items():
        values = np.asarray(payload["values"], dtype=float).reshape(-1)
        analytic_lap = np.asarray(payload["laplacian"], dtype=float).reshape(-1)
        induced = apply_point_set_operator(values=values, operator=operator)
        error = point_set_l2_error(induced, analytic_lap, valid_mask=valid_mask, h=h)
        results[name] = {
            "error_l2": error,
            "valid_nodes": int(np.sum(valid_mask)),
            "total_nodes": int(values.size),
        }

    return {
        "n_side": n_side,
        "nodes": int(points.shape[0]),
        "h": h,
        "epsilon_coeff": epsilon_coeff,
        "epsilon_k": epsilon_k,
        "jitter_fraction": jitter_fraction,
        "seed": seed,
        "valid_nodes": int(np.sum(valid_mask)),
        "operator": {
            "cutoff_radius": float(operator["cutoff_radius"]),
            "average_neighbor_count": float(operator["average_neighbor_count"]),
            "min_neighbor_count": int(operator["min_neighbor_count"]),
            "median_condition_number": float(operator["median_condition_number"]),
            "max_condition_number": float(operator["max_condition_number"]),
            "valid_stencil_fraction": float(operator["valid_stencil_fraction"]),
        },
        "functions": results,
    }


def group_cases(cases: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    grouped: dict[str, dict[str, Any]] = {}
    for func_name in cases[0]["functions"].keys():
        grouped[func_name] = {}
        for eps in sorted({float(case["epsilon_coeff"]) for case in cases}):
            subset = sorted([case for case in cases if float(case["epsilon_coeff"]) == eps], key=lambda item: float(item["h"]))
            h_values = np.array([float(case["h"]) for case in subset], dtype=float)
            errors = np.array([float(case["functions"][func_name]["error_l2"]) for case in subset], dtype=float)
            grouped[func_name][str(eps)] = {
                "h_values": h_values.tolist(),
                "errors": errors.tolist(),
                "fit": fit_error_scaling(h_values, errors),
            }
    return grouped


def summarize_baseline(cases: list[dict[str, Any]], grouped: dict[str, dict[str, Any]]) -> tuple[str, str]:
    quadratic_errors = [float(case["functions"]["f1_quadratic"]["error_l2"]) for case in cases]
    local_orders = [
        grouped[func_name][eps]["fit"]["order_p"]
        for func_name in ("f2_trigonometric", "f3_gaussian")
        for eps in ("0.5", "1.0")
        if eps in grouped[func_name]
    ]
    if max(quadratic_errors) < 1.0e-10 and min(local_orders) > 1.0:
        observation = "the kernel-induced graph operator reproduces the quadratic test exactly and converges toward the continuum Laplacian on the cubic interior scan"
        conclusion = "after discrete second-moment normalization, the cubic baseline remains a clean scalar continuum control; broader kernels remain less accurate at finite resolution"
    else:
        observation = "the cubic baseline shows only partial or inconsistent scalar convergence"
        conclusion = "the current stencil or normalization no longer gives a clean continuum control on the interior cubic scan"
    return observation, conclusion


def summarize_point_set(cases: list[dict[str, Any]], grouped: dict[str, dict[str, Any]]) -> tuple[str, str]:
    quadratic_errors = [float(case["functions"]["f1_quadratic"]["error_l2"]) for case in cases]
    smooth_orders = [
        grouped[func_name][eps]["fit"]["order_p"]
        for func_name in ("f2_trigonometric", "f3_gaussian")
        for eps in grouped[func_name]
    ]
    valid_fractions = [float(case["operator"]["valid_stencil_fraction"]) for case in cases]
    if max(quadratic_errors) < 5.0e-10 and min(smooth_orders) > 0.9 and min(valid_fractions) > 0.95:
        observation = "the kernel-weighted local operator remains convergent on weakly perturbed 3D point sets"
        conclusion = "scalar continuum control survives mild geometric disorder in the local-kernel regime once the operator is recovered from the same kernel weights by local quadratic reproduction"
    else:
        observation = "the weakly perturbed point-set scan only partially preserves scalar convergence"
        conclusion = "mild geometric disorder still distorts the current operator recovery enough that CP1 is not yet closed on point-set controls"
    return observation, conclusion


def summarize_boundary(cases: list[dict[str, Any]], grouped: dict[str, dict[str, Any]]) -> tuple[str, str]:
    dirichlet_orders = [
        grouped["dirichlet_trigonometric"][eps]["fit"]["order_p"]
        for eps in grouped["dirichlet_trigonometric"]
    ]
    neumann_orders = [
        grouped["neumann_trigonometric"][eps]["fit"]["order_p"]
        for eps in grouped["neumann_trigonometric"]
    ]
    if min(dirichlet_orders) > 1.0 and min(neumann_orders) > 1.0:
        observation = "reflected boundary controls retain convergence on both Dirichlet and Neumann-compatible test families"
        conclusion = "the same local kernel regime remains stable when the cubic scan is extended from interior-only error to explicit homogeneous boundary-control families"
    else:
        observation = "the reflected boundary controls are only partially convergent across the current scan"
        conclusion = "boundary handling is still too resolution-sensitive to count as a clean CP1 boundary-control pass"
    return observation, conclusion


def summarize_overall(
    baseline: dict[str, Any],
    point_set: dict[str, Any],
    boundary: dict[str, Any],
) -> tuple[str, str]:
    success = all(
        [
            "clean scalar continuum control" in baseline["conclusion"],
            "scalar continuum control survives mild geometric disorder" in point_set["conclusion"],
            "same local kernel regime remains stable" in boundary["conclusion"],
        ]
    )
    if success:
        observation = "the scalar operator remains convergent on the cubic baseline, weakly perturbed point sets, and reflected Dirichlet or Neumann boundary controls"
        conclusion = "CP1 now extends beyond the clean interior cubic scan: the local-kernel scalar branch has bounded continuum control across mild geometric disorder and explicit boundary families, while broader kernels remain the finite-resolution weak point"
    else:
        observation = "the scalar operator retains some continuum structure, but at least one CP1 control family is still unstable"
        conclusion = "CP1 is only partially closed: the cubic baseline remains informative, but point-set or boundary-control robustness still needs more work"
    return observation, conclusion


def make_baseline_plots(cases: list[dict[str, Any]], grouped: dict[str, dict[str, Any]], stamp: str) -> list[str]:
    plot_paths: list[str] = []

    error_path = PLOTS / f"{stamp}_operator_error_vs_h.png"
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharex=True)
    for ax, func_name in zip(axes, grouped.keys()):
        for eps, payload in grouped[func_name].items():
            h_values = np.array(payload["h_values"], dtype=float)
            errors = np.array(payload["errors"], dtype=float)
            order = float(payload["fit"]["order_p"])
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
    finest = max(cases, key=lambda item: int(item["n_side"]))
    xline = np.linspace(0.0, 1.0, int(finest["n_side"]))
    for ax, func_name in zip(axes, finest["functions"].keys()):
        analytic = np.array(finest["functions"][func_name]["line_analytic"], dtype=float)
        ax.plot(xline, analytic, color="k", linewidth=2.0, label="analytic")
        for eps in sorted({float(case["epsilon_coeff"]) for case in cases}):
            case = next(item for item in cases if int(item["n_side"]) == int(finest["n_side"]) and float(item["epsilon_coeff"]) == eps)
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


def make_point_set_plot(grouped: dict[str, dict[str, Any]], stamp: str) -> str:
    plot_path = PLOTS / f"{stamp}_point_set_operator_error_vs_h.png"
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharex=True)
    for ax, func_name in zip(axes, grouped.keys()):
        for eps, payload in grouped[func_name].items():
            h_values = np.array(payload["h_values"], dtype=float)
            errors = np.array(payload["errors"], dtype=float)
            order = float(payload["fit"]["order_p"])
            ax.loglog(h_values, errors, marker="o", label=f"c_eps={eps}, p={order:.2f}")
        ax.set_title(f"{func_name.replace('_', ' ')} (jittered)")
        ax.set_xlabel("h")
        ax.set_ylabel(r"$\| \widehat{L}_h f - \Delta f \|_{L^2}$")
        ax.grid(alpha=0.25)
        ax.legend()
    fig.savefig(plot_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return str(plot_path.relative_to(REPO_ROOT))


def make_boundary_plots(cases: list[dict[str, Any]], grouped: dict[str, dict[str, Any]], stamp: str) -> list[str]:
    plot_paths: list[str] = []

    error_path = PLOTS / f"{stamp}_boundary_condition_error_vs_h.png"
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharex=True)
    for ax, func_name in zip(axes, ("dirichlet_trigonometric", "neumann_trigonometric")):
        for eps, payload in grouped[func_name].items():
            h_values = np.array(payload["h_values"], dtype=float)
            errors = np.array(payload["errors"], dtype=float)
            order = float(payload["fit"]["order_p"])
            ax.loglog(h_values, errors, marker="o", label=f"c_eps={eps}, p={order:.2f}")
        ax.set_title(func_name.replace("_", " "))
        ax.set_xlabel("h")
        ax.set_ylabel(r"$\| \widehat{L}_h f - \Delta f \|_{L^2}$")
        ax.grid(alpha=0.25)
        ax.legend()
    fig.savefig(error_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(error_path.relative_to(REPO_ROOT)))

    profile_path = PLOTS / f"{stamp}_boundary_condition_profiles.png"
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    finest = max(cases, key=lambda item: int(item["n_side"]))
    xline = np.linspace(0.0, 1.0, int(finest["n_side"]))
    for ax, func_name in zip(axes, ("dirichlet_trigonometric", "neumann_trigonometric")):
        analytic = np.array(finest["functions"][func_name]["line_analytic"], dtype=float)
        ax.plot(xline, analytic, color="k", linewidth=2.0, label="analytic")
        for eps in sorted({float(case["epsilon_coeff"]) for case in cases}):
            case = next(item for item in cases if int(item["n_side"]) == int(finest["n_side"]) and float(item["epsilon_coeff"]) == eps)
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
    payload = json.dumps(result, indent=2)
    stamped.write_text(payload, encoding="utf-8")
    latest.write_text(payload, encoding="utf-8")
    return str(stamped.relative_to(REPO_ROOT)), str(latest.relative_to(REPO_ROOT))


def append_log(result_path: str, plot_paths: list[str], config: dict[str, Any], observation: str, conclusion: str) -> None:
    with EXPERIMENT_LOG.open("a", encoding="utf-8") as handle:
        handle.write("\n## Kernel Laplacian convergence\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            "- Config: "
            f"baseline_n_sides={config['n_sides']}, "
            f"baseline_eps={config['epsilon_coeffs']}, "
            f"point_set_n_sides={config['point_set_n_sides']}, "
            f"point_set_eps={config['point_set_epsilon_coeffs']}, "
            f"jitter_fraction={config['point_set_jitter_fraction']}, "
            f"boundary_n_sides={config['boundary_n_sides']}, "
            f"boundary_eps={config['boundary_epsilon_coeffs']}, "
            f"cutoff_factor={config['cutoff_factor']}\n"
        )
        handle.write(f"- Results: `{result_path}`\n")
        handle.write("- Plots: " + ", ".join(f"`{path}`" for path in plot_paths) + "\n")
        handle.write(f"- Observation: {observation}\n")
        handle.write(f"- Conclusion: {conclusion}\n")


def run_kernel_operator_convergence(
    config: dict[str, Any] | None = None,
    config_path: Path | None = None,
) -> dict[str, Any]:
    cfg = load_config(config=config, config_path=config_path)

    baseline_cases = [
        run_baseline_case(n_side=int(n_side), epsilon_coeff=float(eps), cutoff_factor=float(cfg["cutoff_factor"]))
        for eps in cfg["epsilon_coeffs"]
        for n_side in cfg["n_sides"]
    ]
    baseline_grouped = group_cases(baseline_cases)
    baseline_observation, baseline_conclusion = summarize_baseline(baseline_cases, baseline_grouped)

    point_set_cases = [
        run_point_set_case(
            n_side=int(n_side),
            epsilon_coeff=float(eps),
            cutoff_factor=float(cfg["cutoff_factor"]),
            jitter_fraction=float(cfg["point_set_jitter_fraction"]),
            seed=int(cfg["point_set_seed"]),
            boundary_buffer_factor=float(cfg["point_set_boundary_buffer_factor"]),
        )
        for eps in cfg["point_set_epsilon_coeffs"]
        for n_side in cfg["point_set_n_sides"]
    ]
    point_set_grouped = group_cases(point_set_cases)
    point_set_observation, point_set_conclusion = summarize_point_set(point_set_cases, point_set_grouped)

    boundary_cases = [
        run_boundary_case(n_side=int(n_side), epsilon_coeff=float(eps), cutoff_factor=float(cfg["cutoff_factor"]))
        for eps in cfg["boundary_epsilon_coeffs"]
        for n_side in cfg["boundary_n_sides"]
    ]
    boundary_grouped = group_cases(boundary_cases)
    boundary_observation, boundary_conclusion = summarize_boundary(boundary_cases, boundary_grouped)

    overall_observation, overall_conclusion = summarize_overall(
        baseline={"observation": baseline_observation, "conclusion": baseline_conclusion},
        point_set={"observation": point_set_observation, "conclusion": point_set_conclusion},
        boundary={"observation": boundary_observation, "conclusion": boundary_conclusion},
    )

    return {
        "experiment": "kernel_operator_convergence",
        "config": cfg,
        "baseline_cubic_interior": {
            "derivation": {
                "graph_laplacian": "L = D - W with W_ij = exp(-|x_i-x_j|^2 / epsilon_k) on the cutoff stencil",
                "induced_operator": "Lhat_h = -(2 / (mu_2 h^2)) L, with mu_2 the discrete second moment of the kernel stencil",
                "target": "Lhat_h f -> Delta f on interior nodes as h -> 0",
            },
            "cases": baseline_cases,
            "grouped": baseline_grouped,
            "observation": baseline_observation,
            "conclusion": baseline_conclusion,
        },
        "point_set_control": {
            "derivation": {
                "local_operator": "fit a quadratic Taylor polynomial around each point with the same kernel weights and read off the Laplacian from the quadratic coefficients",
                "target": "the kernel-weighted local recovery remains convergent on weakly jittered point sets after excluding a boundary buffer matched to the kernel radius",
            },
            "cases": point_set_cases,
            "grouped": point_set_grouped,
            "observation": point_set_observation,
            "conclusion": point_set_conclusion,
        },
        "boundary_controls": {
            "derivation": {
                "boundary_extension": "extend the same stencil by homogeneous reflected ghost values; odd parity for Dirichlet and even parity for Neumann",
                "target": "the induced operator remains convergent on full-domain boundary-compatible test families",
            },
            "cases": boundary_cases,
            "grouped": boundary_grouped,
            "observation": boundary_observation,
            "conclusion": boundary_conclusion,
        },
        "observation": overall_observation,
        "conclusion": overall_conclusion,
    }


def main() -> None:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    result = run_kernel_operator_convergence()
    plots = make_baseline_plots(result["baseline_cubic_interior"]["cases"], result["baseline_cubic_interior"]["grouped"], stamp=stamp)
    plots.append(make_point_set_plot(result["point_set_control"]["grouped"], stamp=stamp))
    plots.extend(make_boundary_plots(result["boundary_controls"]["cases"], result["boundary_controls"]["grouped"], stamp=stamp))
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
    for section_name in ("baseline_cubic_interior", "point_set_control", "boundary_controls"):
        print(f"{section_name} observation =", result[section_name]["observation"])
        print(f"{section_name} conclusion =", result[section_name]["conclusion"])
        for func_name, by_eps in result[section_name]["grouped"].items():
            for eps, payload in by_eps.items():
                print(
                    section_name,
                    func_name,
                    f"c_eps={eps}",
                    f"order={payload['fit']['order_p']:.4f}",
                    f"R2={payload['fit']['r2']:.4f}",
                )
    print("observation =", result["observation"])
    print("conclusion =", result["conclusion"])


if __name__ == "__main__":
    main()
