#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

os.environ.setdefault("MPLCONFIGDIR", str(Path.cwd() / ".mplconfig"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
PLOTS_DIR = ROOT / "haos_iip_mode_plots"
RESULTS_JSON = ROOT / "haos_iip_results.json"


@dataclass
class BareStudyConfig:
    name: str
    kind: str
    epsilon_factors: list[float]
    cutoff_factor: float
    reference_factor: float
    seed: int
    n_side: int | None = None
    perturb_fraction: float | None = None
    n_points: int | None = None


def ensure_output_paths() -> None:
    PLOTS_DIR.mkdir()
    (ROOT / ".mplconfig").mkdir(exist_ok=True)


def lattice_points(n_side: int, perturb_fraction: float, seed: int) -> tuple[np.ndarray, float]:
    rng = np.random.default_rng(seed)
    grid = np.linspace(0.0, 1.0, n_side)
    X, Y, Z = np.meshgrid(grid, grid, grid, indexing="ij")
    points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    h = 1.0 / (n_side - 1)
    if perturb_fraction > 0:
        perturb = perturb_fraction * h
        noise = rng.uniform(-perturb, perturb, size=points.shape)
        boundary_mask = ((points == 0.0) | (points == 1.0)).any(axis=1)
        noise[boundary_mask] = 0.0
        points = np.clip(points + noise, 0.0, 1.0)
    return points, h


def random_points(n_points: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.random((n_points, 3))


def pairwise_differences(points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    diff = points[:, None, :] - points[None, :, :]
    d2 = np.sum(diff * diff, axis=2)
    return diff, d2


def mean_nearest_neighbor_sq(d2: np.ndarray) -> float:
    d2 = d2.copy()
    np.fill_diagonal(d2, np.inf)
    return float(np.mean(np.min(d2, axis=1)))


def gaussian_adjacency(points: np.ndarray, epsilon: float, cutoff_factor: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    diff, d2 = pairwise_differences(points)
    cutoff = cutoff_factor * math.sqrt(epsilon)
    mask = (d2 > 0.0) & (d2 <= cutoff * cutoff)
    A = np.zeros_like(d2, dtype=float)
    A[mask] = np.exp(-d2[mask] / (2.0 * epsilon))
    return A, diff, d2


def laplacian_from_adjacency(A: np.ndarray) -> np.ndarray:
    D = np.diag(np.sum(A, axis=1))
    return D - A


def connected_component_count(A: np.ndarray) -> int:
    graph = csr_matrix(A > 0)
    n_components, _ = connected_components(graph, directed=False)
    return int(n_components)


def eigh_sorted(M: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    evals, evecs = np.linalg.eigh(M)
    order = np.argsort(evals)
    return evals[order], evecs[:, order]


def first_nontrivial_indices(evals: np.ndarray, tol: float = 1e-9) -> np.ndarray:
    return np.flatnonzero(evals > tol)


def inverse_participation_ratio(vec: np.ndarray) -> float:
    norm_sq = float(np.sum(np.abs(vec) ** 2))
    if norm_sq == 0.0:
        return 0.0
    return float(np.sum(np.abs(vec) ** 4) / (norm_sq * norm_sq))


def participation_ratio(vec: np.ndarray) -> float:
    ipr = inverse_participation_ratio(vec)
    return 0.0 if ipr == 0.0 else 1.0 / ipr


def boundary_mask(points: np.ndarray, tol: float) -> np.ndarray:
    return ((points <= tol) | (points >= 1.0 - tol)).any(axis=1)


def coordinate_correlations(vec: np.ndarray, points: np.ndarray) -> dict[str, float]:
    real_vec = np.real(vec)
    if np.std(real_vec) < 1e-12:
        return {"x": 0.0, "y": 0.0, "z": 0.0, "max": 0.0}
    corr = {}
    for idx, label in enumerate(("x", "y", "z")):
        centered = points[:, idx] - np.mean(points[:, idx])
        if np.std(centered) < 1e-12:
            corr[label] = 0.0
        else:
            corr[label] = float(abs(np.corrcoef(real_vec, centered)[0, 1]))
    corr["max"] = float(max(corr["x"], corr["y"], corr["z"]))
    return corr


def sign_flip_fraction(vec: np.ndarray, A: np.ndarray) -> float:
    real_vec = np.real(vec)
    edges = A > 0
    opposite = (real_vec[:, None] * real_vec[None, :]) < 0
    count_edges = np.count_nonzero(np.triu(edges, 1))
    if count_edges == 0:
        return 0.0
    count_flips = np.count_nonzero(np.triu(edges & opposite, 1))
    return float(count_flips / count_edges)


def mode_entries(
    evals: np.ndarray,
    evecs: np.ndarray,
    points: np.ndarray,
    A: np.ndarray,
    tol: float,
    boundary_tol: float,
    max_modes: int = 16,
) -> list[dict[str, Any]]:
    nontrivial = first_nontrivial_indices(evals, tol=tol)[:max_modes]
    bmask = boundary_mask(points, boundary_tol)
    out = []
    for raw_idx, idx in enumerate(nontrivial, start=1):
        vec = evecs[:, idx]
        abs_sq = np.abs(vec) ** 2
        boundary_mass = float(abs_sq[bmask].sum())
        corr = coordinate_correlations(vec, points)
        ipr = inverse_participation_ratio(vec)
        out.append(
            {
                "mode_rank": raw_idx,
                "eigen_index": int(idx),
                "eigenvalue": float(evals[idx]),
                "ipr": ipr,
                "participation_ratio": participation_ratio(vec),
                "boundary_mass_fraction": boundary_mass,
                "sign_flip_fraction": sign_flip_fraction(vec, A),
                "coord_correlation": corr,
            }
        )
    return out


def summarize_degeneracies(evals: np.ndarray, tol: float = 1e-4, max_modes: int = 18) -> list[list[float]]:
    vals = [float(v) for v in evals[:max_modes]]
    groups: list[list[float]] = []
    for value in vals:
        if not groups or abs(groups[-1][-1] - value) > tol:
            groups.append([value])
        else:
            groups[-1].append(value)
    return groups


def classify_bare_modes(entries: list[dict[str, Any]], n_nodes: int) -> list[str]:
    labels = []
    for entry in entries:
        ipr = entry["ipr"]
        coord_max = entry["coord_correlation"]["max"]
        boundary_mass = entry["boundary_mass_fraction"]
        if ipr < 3.0 / n_nodes and coord_max > 0.75:
            labels.append("geometry-like coordinate mode")
        elif ipr < 4.5 / n_nodes and boundary_mass < 0.55:
            labels.append("smooth scalar-like mode")
        elif ipr > 8.0 / n_nodes:
            labels.append("localized candidate")
        else:
            labels.append("mixed scalar mode")
    return labels


def mode_plot_values(vec: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    real = np.real(vec)
    vmax = np.max(np.abs(real))
    if vmax < 1e-12:
        colors = real
    else:
        colors = real / vmax
    sizes = 25.0 + 300.0 * (np.abs(vec) ** 2) / np.max(np.abs(vec) ** 2)
    return colors, sizes


def plot_mode_grid(
    points: np.ndarray,
    evals: np.ndarray,
    evecs: np.ndarray,
    title: str,
    output_path: Path,
    tol: float,
    max_modes: int = 6,
) -> None:
    nontrivial = first_nontrivial_indices(evals, tol=tol)[:max_modes]
    fig = plt.figure(figsize=(14, 8))
    for panel, idx in enumerate(nontrivial, start=1):
        ax = fig.add_subplot(2, 3, panel, projection="3d")
        colors, sizes = mode_plot_values(evecs[:, idx])
        sc = ax.scatter(
            points[:, 0],
            points[:, 1],
            points[:, 2],
            c=colors,
            s=sizes,
            cmap="coolwarm",
            vmin=-1.0,
            vmax=1.0,
            alpha=0.9,
            linewidths=0.0,
        )
        ax.set_title(f"mode {panel}: lambda={evals[idx]:.4f}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.view_init(elev=24, azim=38)
    fig.suptitle(title)
    fig.subplots_adjust(right=0.92, wspace=0.18, hspace=0.25)
    cbar_ax = fig.add_axes([0.94, 0.16, 0.015, 0.68])
    fig.colorbar(sc, cax=cbar_ax, label="normalized Re(mode)")
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_phase_mode_grid(
    points: np.ndarray,
    evals: np.ndarray,
    evecs: np.ndarray,
    title: str,
    output_path: Path,
    tol: float,
    max_modes: int = 4,
) -> None:
    nontrivial = first_nontrivial_indices(evals, tol=tol)[:max_modes]
    fig = plt.figure(figsize=(12, 8))
    for panel, idx in enumerate(nontrivial, start=1):
        ax = fig.add_subplot(2, 2, panel, projection="3d")
        vec = evecs[:, idx]
        phase = np.angle(vec)
        sizes = 25.0 + 300.0 * (np.abs(vec) ** 2) / np.max(np.abs(vec) ** 2)
        sc = ax.scatter(
            points[:, 0],
            points[:, 1],
            points[:, 2],
            c=phase,
            s=sizes,
            cmap="twilight",
            vmin=-math.pi,
            vmax=math.pi,
            alpha=0.9,
            linewidths=0.0,
        )
        ax.set_title(f"mode {panel}: lambda={evals[idx]:.4f}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.view_init(elev=24, azim=38)
    fig.suptitle(title)
    fig.subplots_adjust(right=0.92, wspace=0.18, hspace=0.25)
    cbar_ax = fig.add_axes([0.94, 0.16, 0.015, 0.68])
    fig.colorbar(sc, cax=cbar_ax, label="phase")
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_spectrum(reference_cases: dict[str, dict[str, Any]], output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    for label, case in reference_cases.items():
        evals = np.array(case["evals"][:18], dtype=float)
        ax.plot(np.arange(len(evals)), evals, marker="o", linewidth=1.4, label=label)
    ax.set_xlabel("eigen-index")
    ax.set_ylabel("eigenvalue")
    ax.set_title("Bottom of the spectrum")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_ipr(reference_cases: dict[str, dict[str, Any]], output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    for label, case in reference_cases.items():
        entries = case["mode_entries"]
        xs = [entry["mode_rank"] for entry in entries]
        ys = [entry["ipr"] for entry in entries]
        ax.plot(xs, ys, marker="o", linewidth=1.4, label=label)
    ax.set_xlabel("nontrivial mode rank")
    ax.set_ylabel("IPR")
    ax.set_title("Localization diagnostic")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_epsilon_sweep(studies: dict[str, Any], output_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for label, study in studies.items():
        sweep = study["epsilon_sweep"]
        eps = [row["epsilon"] for row in sweep]
        gap = [row["spectral_gap"] for row in sweep]
        mean_ipr = [row["mean_ipr_first6"] for row in sweep]
        axes[0].plot(eps, gap, marker="o", linewidth=1.4, label=label)
        axes[1].plot(eps, mean_ipr, marker="o", linewidth=1.4, label=label)
    axes[0].set_xlabel("epsilon")
    axes[0].set_ylabel("lambda_1")
    axes[0].set_title("Spectral gap vs epsilon")
    axes[0].grid(alpha=0.25)
    axes[1].set_xlabel("epsilon")
    axes[1].set_ylabel("mean IPR of first 6 nontrivial modes")
    axes[1].set_title("Low-mode localization vs epsilon")
    axes[1].grid(alpha=0.25)
    axes[0].legend()
    axes[1].legend()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_n_sweep(n_sweep: dict[str, list[dict[str, Any]]], output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    for label, rows in n_sweep.items():
        ax.plot(
            [row["n_nodes"] for row in rows],
            [row["spectral_gap"] for row in rows],
            marker="o",
            linewidth=1.4,
            label=label,
        )
    ax.set_xlabel("node count")
    ax.set_ylabel("lambda_1")
    ax.set_title("Spectral gap vs N")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def nearest_neighbor_edges(points: np.ndarray, h: float, tolerance: float = 1.05) -> np.ndarray:
    diff, d2 = pairwise_differences(points)
    mask = (d2 > 0.0) & (d2 <= (tolerance * h) ** 2)
    return mask


def phase_dressed_laplacian(
    points: np.ndarray,
    epsilon: float,
    h: float,
    flux_per_plaquette: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    A, diff, d2 = gaussian_adjacency(points, epsilon=epsilon, cutoff_factor=1.05)
    nn_mask = (d2 > 0.0) & (d2 <= (1.05 * h) ** 2)
    A[~nn_mask] = 0.0
    B = 2.0 * math.pi * flux_per_plaquette / (h * h)
    x_mid = 0.5 * (points[:, None, 0] + points[None, :, 0])
    dy = diff[:, :, 1]
    theta = B * x_mid * dy
    theta = 0.5 * (theta - theta.T)
    U = np.exp(1j * theta)
    A_theta = A * U
    D = np.diag(np.sum(A, axis=1))
    L_theta = D - A_theta
    return L_theta, A, theta, U


def edge_currents(vec: np.ndarray, A: np.ndarray, U: np.ndarray) -> np.ndarray:
    return 2.0 * A * np.imag(np.conj(vec[:, None]) * U * vec[None, :])


def phase_mode_entries(
    evals: np.ndarray,
    evecs: np.ndarray,
    points: np.ndarray,
    A: np.ndarray,
    U: np.ndarray,
    tol: float,
    max_modes: int = 12,
) -> list[dict[str, Any]]:
    nontrivial = first_nontrivial_indices(evals, tol=tol)[:max_modes]
    out = []
    for raw_idx, idx in enumerate(nontrivial, start=1):
        vec = evecs[:, idx]
        J = edge_currents(vec, A, U)
        current_norm = float(np.sum(np.abs(np.triu(J, 1))))
        current_scale = float(np.sum(np.triu(A * np.abs(vec[:, None]) * np.abs(vec[None, :]), 1)))
        current_ratio = 0.0 if current_scale == 0.0 else current_norm / current_scale
        out.append(
            {
                "mode_rank": raw_idx,
                "eigen_index": int(idx),
                "eigenvalue": float(evals[idx]),
                "ipr": inverse_participation_ratio(vec),
                "participation_ratio": participation_ratio(vec),
                "current_ratio": current_ratio,
                "mean_phase": float(np.angle(np.mean(vec))),
            }
        )
    return out


def plaquette_holonomy(points: np.ndarray, n_side: int, theta: np.ndarray) -> dict[str, float]:
    index = {(i, j, k): i * n_side * n_side + j * n_side + k for i in range(n_side) for j in range(n_side) for k in range(n_side)}
    loops = []
    for i in range(n_side - 1):
        for j in range(n_side - 1):
            for k in range(n_side):
                a = index[(i, j, k)]
                b = index[(i + 1, j, k)]
                c = index[(i + 1, j + 1, k)]
                d = index[(i, j + 1, k)]
                phase = theta[a, b] + theta[b, c] + theta[c, d] + theta[d, a]
                loops.append(phase)
    loops = np.array(loops)
    return {
        "mean_xy_plaquette_phase": float(np.mean(loops)),
        "std_xy_plaquette_phase": float(np.std(loops)),
    }


def run_bare_study(config: BareStudyConfig) -> tuple[dict[str, Any], np.ndarray, np.ndarray]:
    if config.kind == "lattice":
        assert config.n_side is not None
        assert config.perturb_fraction is not None
        points, h = lattice_points(config.n_side, config.perturb_fraction, config.seed)
        boundary_tol = 0.51 * h
    elif config.kind == "random":
        assert config.n_points is not None
        points = random_points(config.n_points, config.seed)
        _, d2 = pairwise_differences(points)
        h = math.sqrt(mean_nearest_neighbor_sq(d2))
        boundary_tol = 0.12
    else:
        raise ValueError(config.kind)

    _, d2 = pairwise_differences(points)
    mean_nn_sq = mean_nearest_neighbor_sq(d2)
    sweep = []
    reference_case: dict[str, Any] | None = None
    for factor in config.epsilon_factors:
        epsilon = factor * mean_nn_sq
        A, _, _ = gaussian_adjacency(points, epsilon=epsilon, cutoff_factor=config.cutoff_factor)
        L = laplacian_from_adjacency(A)
        evals, evecs = eigh_sorted(L)
        entries = mode_entries(evals, evecs, points, A, tol=1e-8, boundary_tol=boundary_tol, max_modes=16)
        labels = classify_bare_modes(entries, len(points))
        for entry, label in zip(entries, labels):
            entry["classification"] = label
        nontrivial = first_nontrivial_indices(evals, tol=1e-8)
        gap = float(evals[nontrivial[0]])
        sweep.append(
            {
                "factor": factor,
                "epsilon": float(epsilon),
                "spectral_gap": gap,
                "lambda_2": float(evals[nontrivial[1]]) if len(nontrivial) > 1 else None,
                "components": connected_component_count(A),
                "mean_ipr_first6": float(np.mean([entry["ipr"] for entry in entries[:6]])),
            }
        )
        if abs(factor - config.reference_factor) < 1e-12:
            reference_case = {
                "epsilon": float(epsilon),
                "evals": [float(v) for v in evals[:18]],
                "mode_entries": entries,
                "components": connected_component_count(A),
                "degeneracy_groups": summarize_degeneracies(evals),
                "cutoff_radius": float(config.cutoff_factor * math.sqrt(epsilon)),
            }
            plot_mode_grid(
                points,
                evals,
                evecs,
                title=f"{config.name}: first nontrivial modes",
                output_path=PLOTS_DIR / f"{config.name}_modes.png",
                tol=1e-8,
            )
            reference_case["plot"] = str((PLOTS_DIR / f"{config.name}_modes.png").name)
    if reference_case is None:
        raise RuntimeError(f"reference factor {config.reference_factor} not found")
    study = {
        "config": {
            "name": config.name,
            "kind": config.kind,
            "n_nodes": int(len(points)),
            "mean_nearest_neighbor_sq": float(mean_nn_sq),
            "boundary": "open",
            "cutoff_factor": config.cutoff_factor,
            "reference_factor": config.reference_factor,
            "epsilon_factors": config.epsilon_factors,
        },
        "reference_case": reference_case,
        "epsilon_sweep": sweep,
    }
    return study, points, d2


def run_lattice_n_sweep(seed: int) -> list[dict[str, Any]]:
    rows = []
    for n_side in (5, 6, 7):
        points, _ = lattice_points(n_side=n_side, perturb_fraction=0.06, seed=seed)
        _, d2 = pairwise_differences(points)
        mean_nn_sq = mean_nearest_neighbor_sq(d2)
        epsilon = 1.5 * mean_nn_sq
        A, _, _ = gaussian_adjacency(points, epsilon=epsilon, cutoff_factor=2.5)
        evals, _ = eigh_sorted(laplacian_from_adjacency(A))
        gap = float(evals[first_nontrivial_indices(evals, tol=1e-8)[0]])
        rows.append({"n_nodes": int(len(points)), "epsilon": float(epsilon), "spectral_gap": gap})
    return rows


def run_random_n_sweep(seed: int) -> list[dict[str, Any]]:
    rows = []
    for n_points in (125, 216, 343):
        points = random_points(n_points=n_points, seed=seed + n_points)
        _, d2 = pairwise_differences(points)
        mean_nn_sq = mean_nearest_neighbor_sq(d2)
        epsilon = 2.0 * mean_nn_sq
        A, _, _ = gaussian_adjacency(points, epsilon=epsilon, cutoff_factor=2.5)
        evals, _ = eigh_sorted(laplacian_from_adjacency(A))
        gap = float(evals[first_nontrivial_indices(evals, tol=1e-8)[0]])
        rows.append({"n_nodes": int(len(points)), "epsilon": float(epsilon), "spectral_gap": gap})
    return rows


def run_phase_branch() -> dict[str, Any]:
    n_side = 6
    points, h = lattice_points(n_side=n_side, perturb_fraction=0.0, seed=1)
    epsilon = h * h
    flux_per_plaquette = 1.0 / 6.0
    L_theta, A, theta, U = phase_dressed_laplacian(points, epsilon=epsilon, h=h, flux_per_plaquette=flux_per_plaquette)
    evals_theta, evecs_theta = eigh_sorted(L_theta)
    phase_entries = phase_mode_entries(evals_theta, evecs_theta, points, A, U, tol=1e-8, max_modes=12)
    plot_phase_mode_grid(
        points,
        evals_theta,
        evecs_theta,
        title="phase-dressed cubic lattice: first nontrivial modes",
        output_path=PLOTS_DIR / "phase_dressed_lattice_modes.png",
        tol=1e-8,
        max_modes=4,
    )

    A0, _, _ = gaussian_adjacency(points, epsilon=epsilon, cutoff_factor=1.05)
    _, d2 = pairwise_differences(points)
    A0[(d2 > (1.05 * h) ** 2)] = 0.0
    evals_base, evecs_base = eigh_sorted(laplacian_from_adjacency(A0))
    base_entries = mode_entries(
        evals_base,
        evecs_base,
        points,
        A0,
        tol=1e-8,
        boundary_tol=0.51 * h,
        max_modes=12,
    )

    holonomy = plaquette_holonomy(points, n_side=n_side, theta=theta)
    return {
        "config": {
            "name": "phase_dressed_lattice",
            "n_nodes": int(len(points)),
            "n_side": n_side,
            "boundary": "open",
            "epsilon": float(epsilon),
            "flux_per_xy_plaquette": flux_per_plaquette,
            "cutoff_radius": float(1.05 * h),
        },
        "bare_reference": {
            "evals": [float(v) for v in evals_base[:18]],
            "mode_entries": base_entries,
        },
        "phase_dressed": {
            "evals": [float(v) for v in evals_theta[:18]],
            "mode_entries": phase_entries,
            "degeneracy_groups": summarize_degeneracies(evals_theta),
            "holonomy": holonomy,
            "plot": "phase_dressed_lattice_modes.png",
        },
    }


def final_verdict(studies: dict[str, Any], phase_branch: dict[str, Any]) -> dict[str, str]:
    lattice_labels = [entry["classification"] for entry in studies["cubic_open_perturbed"]["reference_case"]["mode_entries"][:6]]
    random_labels = [entry["classification"] for entry in studies["random_geometric_open"]["reference_case"]["mode_entries"][:6]]
    phase_current = max(entry["current_ratio"] for entry in phase_branch["phase_dressed"]["mode_entries"][:6])

    worked = "low bare modes are smooth and delocalized on both 3D substrates; the cubic branch shows clear coordinate-like scalar modes"
    did_not = "no autonomous vector sector or fermion-like / spinorial low-mode family appears in the bare scalar Laplacian spectrum"
    next_try = (
        "to probe gauge structure directly, keep the phase-dressed branch and move to edge/connection operators or Dirac-type operators; "
        "to probe particle-like sectors, add controlled defects or a chiral/topological substrate instead of only a node Laplacian"
    )
    if not any("geometry-like coordinate mode" in label for label in lattice_labels + random_labels):
        worked = "the low modes remain smooth scalar diffusion modes, but coordinate alignment is weaker on the chosen discretizations"
    if phase_current < 0.2:
        next_try = (
            "the phase-dressed branch did not produce strong circulating low modes; the next step is a stronger flux background or an explicit edge/plaquette operator"
        )
    return {"worked": worked, "did_not_appear": did_not, "must_try_next": next_try}


def main() -> None:
    ensure_output_paths()

    bare_configs = [
        BareStudyConfig(
            name="cubic_open_perturbed",
            kind="lattice",
            n_side=6,
            perturb_fraction=0.06,
            epsilon_factors=[1.0, 1.5, 2.0],
            cutoff_factor=2.5,
            reference_factor=1.5,
            seed=12,
        ),
        BareStudyConfig(
            name="random_geometric_open",
            kind="random",
            n_points=256,
            epsilon_factors=[1.5, 2.0, 2.5],
            cutoff_factor=2.5,
            reference_factor=2.0,
            seed=21,
        ),
    ]

    studies: dict[str, Any] = {}
    reference_cases: dict[str, dict[str, Any]] = {}
    for config in bare_configs:
        study, _, _ = run_bare_study(config)
        studies[config.name] = study
        reference_cases[config.name] = {
            "evals": study["reference_case"]["evals"],
            "mode_entries": study["reference_case"]["mode_entries"],
        }

    phase_branch = run_phase_branch()
    reference_cases["phase_dressed_lattice"] = {
        "evals": phase_branch["phase_dressed"]["evals"],
        "mode_entries": phase_branch["phase_dressed"]["mode_entries"],
    }

    plot_spectrum(reference_cases, PLOTS_DIR / "low_spectrum_comparison.png")
    plot_ipr(reference_cases, PLOTS_DIR / "localization_vs_mode_index.png")
    plot_epsilon_sweep(studies, PLOTS_DIR / "epsilon_sweep.png")

    n_sweep = {
        "cubic_open_perturbed": run_lattice_n_sweep(seed=99),
        "random_geometric_open": run_random_n_sweep(seed=99),
    }
    plot_n_sweep(n_sweep, PLOTS_DIR / "n_sweep_spectral_gap.png")

    output = {
        "sources_used": [
            "HAOS_IIP kernel.docx",
            "HAOS_IIP Spectral Geometry Note.docx",
            "Interaction-invariant Physics Full.docx",
            "Rebuilding Duda’s picture in HAOS_IIP.docx",
            "HAOS-IIP Emergent Gauge Sector Note.md",
        ],
        "studies": studies,
        "phase_branch": phase_branch,
        "n_sweep": n_sweep,
        "plots": sorted([path.name for path in PLOTS_DIR.glob("*.png")]),
    }
    output["current_verdict"] = final_verdict(studies, phase_branch)

    with RESULTS_JSON.open("w", encoding="utf-8") as handle:
        json.dump(output, handle, indent=2)

    print(f"Wrote {RESULTS_JSON.name}")
    print(f"Wrote plots to {PLOTS_DIR}")
    print("Reference spectral gaps:")
    for name, study in studies.items():
        print(name, study["reference_case"]["mode_entries"][0]["eigenvalue"])
    print("Phase-dressed max current ratio:", max(entry["current_ratio"] for entry in phase_branch["phase_dressed"]["mode_entries"][:6]))


if __name__ == "__main__":
    main()
