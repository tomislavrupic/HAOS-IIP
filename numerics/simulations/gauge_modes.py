#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import os
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
PLOTS = REPO_ROOT / "plots"
PLOTS.mkdir(exist_ok=True)
MPLCONFIG = REPO_ROOT / ".mplconfig"
MPLCONFIG.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_CONFIG: dict[str, Any] = {
    "nodes": 500,
    "epsilon": 0.2,
    "random_seed": 42,
    "substrate": "random_geometric",
    "cutoff_factor": 2.5,
    "gauge_substrate": "cubic_lattice",
    "gauge_lattice_side": 6,
    "flux_per_plaquette": 0.2,
    "nearest_neighbor_factor": 1.05,
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        merged.update(json.loads(path.read_text()))
    if config is not None:
        merged.update(config)
    return merged


def build_points(cfg: dict[str, Any]) -> tuple[np.ndarray, str, float]:
    substrate = str(cfg.get("gauge_substrate") or cfg.get("substrate", "cubic_lattice"))
    seed = int(cfg.get("random_seed", 42))
    rng = np.random.default_rng(seed)

    if substrate == "cubic_lattice":
        n_side = int(cfg.get("gauge_lattice_side", cfg.get("lattice_side", 6)))
        grid = np.linspace(0.0, 1.0, n_side)
        X, Y, Z = np.meshgrid(grid, grid, grid, indexing="ij")
        points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
        h = 1.0 / (n_side - 1)
        return points, substrate, h

    if substrate == "random_geometric":
        n = int(cfg.get("nodes", 500))
        points = rng.random((n, 3))
        diff = points[:, None, :] - points[None, :, :]
        d2 = np.sum(diff * diff, axis=2)
        np.fill_diagonal(d2, np.inf)
        h = float(np.sqrt(np.mean(np.min(d2, axis=1))))
        return points, substrate, h

    raise ValueError(f"Unsupported gauge substrate: {substrate}")


def edge_currents(vec: np.ndarray, A: np.ndarray, U: np.ndarray) -> np.ndarray:
    return 2.0 * A * np.imag(np.conj(vec[:, None]) * U * vec[None, :])


def build_covariant_operator(points: np.ndarray, epsilon: float, h: float, cfg: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    diff = points[:, None, :] - points[None, :, :]
    d2 = np.sum(diff * diff, axis=2)
    substrate = str(cfg.get("gauge_substrate") or cfg.get("substrate", "cubic_lattice"))
    if substrate == "cubic_lattice":
        mask = (d2 > 0.0) & (d2 <= (float(cfg.get("nearest_neighbor_factor", 1.05)) * h) ** 2)
    else:
        cutoff = float(cfg.get("cutoff_factor", 2.5)) * math.sqrt(epsilon)
        mask = (d2 > 0.0) & (d2 <= cutoff * cutoff)

    A = np.zeros_like(d2)
    A[mask] = np.exp(-d2[mask] / (2.0 * epsilon))

    flux_per_plaquette = float(cfg.get("flux_per_plaquette", 0.2))
    Bfield = 2.0 * math.pi * flux_per_plaquette / max(h * h, 1e-12)
    x_mid = 0.5 * (points[:, None, 0] + points[None, :, 0])
    theta = Bfield * x_mid * diff[:, :, 1]
    theta = 0.5 * (theta - theta.T)
    U = np.exp(1j * theta)

    D = np.diag(np.sum(A, axis=1))
    L_theta = D - A * U
    return L_theta, A, U


def summarize_gauge(evals: np.ndarray, evecs: np.ndarray, A: np.ndarray, U: np.ndarray) -> tuple[str, str, float]:
    ratios = []
    for idx in range(1, min(5, len(evals))):
        vec = evecs[:, idx]
        J = edge_currents(vec, A, U)
        current_norm = float(np.sum(np.abs(np.triu(J, 1))))
        scale = float(np.sum(np.triu(A * np.abs(vec[:, None]) * np.abs(vec[None, :]), 1)))
        ratios.append(0.0 if scale == 0.0 else current_norm / scale)
    max_ratio = float(max(ratios) if ratios else 0.0)
    if max_ratio > 0.2:
        observation = "phase dressing produces circulating low modes"
        conclusion = "connection-sensitive gauge branch present, but still scalar in background transport"
    else:
        observation = "phase dressing weakly perturbs the low spectrum"
        conclusion = "no strong gauge-like circulation signal in this setup"
    return observation, conclusion, max_ratio


def make_plots(points: np.ndarray, evals: np.ndarray, evecs: np.ndarray, name: str) -> list[str]:
    plot_paths: list[str] = []

    spectrum_path = PLOTS / f"{name}_spectrum.png"
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(range(min(12, len(evals))), [v.real for v in evals[:12]], marker="o")
    ax.set_xlabel("eigen-index")
    ax.set_ylabel("eigenvalue")
    ax.set_title(f"{name} spectrum")
    ax.grid(alpha=0.25)
    fig.savefig(spectrum_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(spectrum_path.relative_to(REPO_ROOT)))

    mode = evecs[:, 1] if evecs.shape[1] > 1 else evecs[:, 0]
    phase = np.angle(mode)
    sizes = 30 + 280 * (np.abs(mode) ** 2) / (np.max(np.abs(mode) ** 2) or 1.0)
    phase_path = PLOTS / f"{name}_phase_mode.png"
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111, projection="3d")
    sc = ax.scatter(points[:, 0], points[:, 1], points[:, 2], c=phase, cmap="twilight", s=sizes)
    ax.set_title(f"{name} first nontrivial phase mode")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    fig.colorbar(sc, ax=ax, shrink=0.7, label="phase")
    fig.savefig(phase_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(phase_path.relative_to(REPO_ROOT)))

    return plot_paths


def run_gauge_test(
    config: dict[str, Any] | None = None,
    with_plots: bool = True,
    plot_name: str = "gauge_modes",
) -> dict[str, Any]:
    cfg = load_config(config)
    points, substrate, h = build_points(cfg)
    epsilon = float(cfg["epsilon"])
    L_theta, A, U = build_covariant_operator(points, epsilon=epsilon, h=h, cfg=cfg)
    evals, evecs = np.linalg.eigh(L_theta)
    order = np.argsort(evals.real)
    evals = np.asarray(evals[order])
    evecs = np.asarray(evecs[:, order])
    observation, conclusion, max_current_ratio = summarize_gauge(evals, evecs, A, U)
    plot_paths = make_plots(points, evals, evecs, name=plot_name) if with_plots else []

    return {
        "experiment": "gauge_sector_test",
        "config": {
            "nodes": int(len(points)),
            "epsilon": epsilon,
            "random_seed": int(cfg["random_seed"]),
            "substrate": substrate,
            "flux_per_plaquette": float(cfg.get("flux_per_plaquette", 0.2)),
        },
        "spectrum": {
            "first_eigenvalues": [float(v.real) for v in evals[:12]],
            "lowest_mode_current_ratio": max_current_ratio,
        },
        "plots": plot_paths,
        "observation": observation,
        "conclusion": conclusion,
    }


def main() -> None:
    result = run_gauge_test()
    print("epsilon_k =", result["config"]["epsilon"])
    print("first eigenvalues =", [round(v, 6) for v in result["spectrum"]["first_eigenvalues"][:8]])
    print("plots =", result["plots"])
    print("observation =", result["observation"])
    print("conclusion =", result["conclusion"])


if __name__ == "__main__":
    main()
