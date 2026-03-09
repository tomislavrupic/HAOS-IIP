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
    "lattice_side": 8,
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        merged.update(json.loads(path.read_text()))
    if config is not None:
        merged.update(config)
    return merged


def build_points(config: dict[str, Any]) -> tuple[np.ndarray, str]:
    substrate = str(config.get("substrate", "random_geometric"))
    seed = int(config.get("random_seed", 42))
    rng = np.random.default_rng(seed)

    if substrate == "random_geometric":
        n = int(config.get("nodes", 500))
        return rng.random((n, 3)), substrate

    if substrate == "cubic_lattice":
        n_side = int(config.get("lattice_side", round(int(config.get("nodes", 512)) ** (1.0 / 3.0))))
        grid = np.linspace(0.0, 1.0, n_side)
        X, Y, Z = np.meshgrid(grid, grid, grid, indexing="ij")
        points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
        return points, substrate

    raise ValueError(f"Unsupported substrate: {substrate}")


def inverse_participation_ratio(vec: np.ndarray) -> float:
    abs_sq = np.abs(vec) ** 2
    norm_sq = float(np.sum(abs_sq))
    if norm_sq == 0.0:
        return 0.0
    return float(np.sum(abs_sq * abs_sq) / (norm_sq * norm_sq))


def build_laplacian(points: np.ndarray, epsilon: float, cutoff_factor: float) -> tuple[np.ndarray, np.ndarray]:
    diff = points[:, None, :] - points[None, :, :]
    d2 = np.sum(diff * diff, axis=2)
    cutoff = cutoff_factor * math.sqrt(epsilon)
    mask = (d2 > 0.0) & (d2 <= cutoff * cutoff)
    A = np.zeros_like(d2)
    A[mask] = np.exp(-d2[mask] / (2.0 * epsilon))
    D = np.diag(np.sum(A, axis=1))
    return D - A, A


def summarize_observation(evals: np.ndarray, evecs: np.ndarray) -> tuple[str, str, float]:
    nontrivial = evals[1:4]
    mean_ipr = float(np.mean([inverse_participation_ratio(evecs[:, i]) for i in range(1, min(4, evecs.shape[1]))]))
    spread = float((np.max(nontrivial) - np.min(nontrivial)) / max(np.mean(nontrivial), 1e-12)) if len(nontrivial) else 0.0
    if spread < 0.08 and mean_ipr < 0.05:
        observation = "three near-degenerate smooth low modes"
        conclusion = "geometry sector consistent"
    elif mean_ipr > 0.08:
        observation = "low spectrum contains localized modes"
        conclusion = "geometry sector mixed with substrate or boundary artifacts"
    else:
        observation = "low spectrum shows smooth scalar-like modes without strong degeneracy"
        conclusion = "scalar geometry branch present but not especially symmetric"
    return observation, conclusion, mean_ipr


def make_plots(points: np.ndarray, evals: np.ndarray, evecs: np.ndarray, name: str) -> list[str]:
    plot_paths: list[str] = []

    spectrum_path = PLOTS / f"{name}_spectrum.png"
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(range(min(12, len(evals))), evals[:12], marker="o")
    ax.set_xlabel("eigen-index")
    ax.set_ylabel("eigenvalue")
    ax.set_title(f"{name} spectrum")
    ax.grid(alpha=0.25)
    fig.savefig(spectrum_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(spectrum_path.relative_to(REPO_ROOT)))

    mode_path = PLOTS / f"{name}_mode1.png"
    vec = evecs[:, 1] if evecs.shape[1] > 1 else evecs[:, 0]
    colors = np.real(vec)
    vmax = np.max(np.abs(colors)) or 1.0
    sizes = 30 + 280 * (np.abs(vec) ** 2) / (np.max(np.abs(vec) ** 2) or 1.0)
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111, projection="3d")
    sc = ax.scatter(points[:, 0], points[:, 1], points[:, 2], c=colors / vmax, cmap="coolwarm", s=sizes, vmin=-1.0, vmax=1.0)
    ax.set_title(f"{name} first nontrivial mode")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    fig.colorbar(sc, ax=ax, shrink=0.7, label="normalized Re(mode)")
    fig.savefig(mode_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(mode_path.relative_to(REPO_ROOT)))

    return plot_paths


def run_laplacian_test(config: dict[str, Any] | None = None) -> dict[str, Any]:
    cfg = load_config(config)
    epsilon = float(cfg["epsilon"])
    cutoff_factor = float(cfg.get("cutoff_factor", 2.5))
    points, substrate = build_points(cfg)
    L, A = build_laplacian(points, epsilon=epsilon, cutoff_factor=cutoff_factor)
    evals, evecs = np.linalg.eigh(L)
    order = np.argsort(evals)
    evals = np.asarray(evals[order], dtype=float)
    evecs = np.asarray(evecs[:, order])
    observation, conclusion, mean_ipr = summarize_observation(evals, evecs)
    plot_paths = make_plots(points, evals, evecs, name="laplacian_modes")

    return {
        "experiment": "laplacian_geometry_test",
        "config": {
            "nodes": int(len(points)),
            "epsilon": epsilon,
            "random_seed": int(cfg["random_seed"]),
            "substrate": substrate,
            "cutoff_factor": cutoff_factor,
        },
        "spectrum": {
            "first_eigenvalues": [float(v) for v in evals[:12]],
            "spectral_gap": float(evals[1]) if len(evals) > 1 else None,
            "mean_ipr_first3": mean_ipr,
            "degree_mean": float(np.mean(np.sum(A > 0.0, axis=1))),
        },
        "plots": plot_paths,
        "observation": observation,
        "conclusion": conclusion,
    }


def main() -> None:
    result = run_laplacian_test()
    print("epsilon_k =", result["config"]["epsilon"])
    print("first eigenvalues =", [round(v, 6) for v in result["spectrum"]["first_eigenvalues"][:8]])
    print("plots =", result["plots"])
    print("observation =", result["observation"])
    print("conclusion =", result["conclusion"])


if __name__ == "__main__":
    main()
