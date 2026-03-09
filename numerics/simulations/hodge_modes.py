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
from matplotlib import cm

DEFAULT_CONFIG: dict[str, Any] = {
    "epsilon": 0.2,
    "random_seed": 42,
    "hodge_substrate": "cubic_lattice",
    "hodge_lattice_side": 6,
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        merged.update(json.loads(path.read_text()))
    if config is not None:
        merged.update(config)
    return merged


def build_cubic_lattice_complex(n_side: int, epsilon: float) -> dict[str, np.ndarray | list[tuple[int, int]] | float | int]:
    grid = np.linspace(0.0, 1.0, n_side)
    X, Y, Z = np.meshgrid(grid, grid, grid, indexing="ij")
    points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    node_index = np.arange(n_side**3).reshape((n_side, n_side, n_side))

    edges: list[tuple[int, int]] = []
    directions: list[np.ndarray] = []
    midpoints: list[np.ndarray] = []
    edge_weights: list[float] = []
    edge_map: dict[tuple[str, int, int, int], int] = {}

    def add_edge(axis: str, start: tuple[int, int, int], stop: tuple[int, int, int]) -> None:
        u = int(node_index[start])
        v = int(node_index[stop])
        idx = len(edges)
        edge_map[(axis, start[0], start[1], start[2])] = idx
        edges.append((u, v))
        vec = points[v] - points[u]
        directions.append(vec / np.linalg.norm(vec))
        midpoints.append(0.5 * (points[u] + points[v]))
        d2 = float(np.dot(vec, vec))
        edge_weights.append(math.exp(-d2 / (2.0 * epsilon)))

    for i in range(n_side - 1):
        for j in range(n_side):
            for k in range(n_side):
                add_edge("x", (i, j, k), (i + 1, j, k))
    for i in range(n_side):
        for j in range(n_side - 1):
            for k in range(n_side):
                add_edge("y", (i, j, k), (i, j + 1, k))
    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side - 1):
                add_edge("z", (i, j, k), (i, j, k + 1))

    n_nodes = len(points)
    n_edges = len(edges)
    B0 = np.zeros((n_edges, n_nodes), dtype=float)
    for edge_idx, (u, v) in enumerate(edges):
        B0[edge_idx, u] = -1.0
        B0[edge_idx, v] = 1.0

    faces: list[list[tuple[int, int]]] = []
    face_weights: list[float] = []

    def add_face(boundary: list[tuple[str, int, int, int, int]]) -> None:
        oriented_edges: list[tuple[int, int]] = []
        boundary_weights: list[float] = []
        for axis, i, j, k, sign in boundary:
            edge_idx = edge_map[(axis, i, j, k)]
            oriented_edges.append((edge_idx, sign))
            boundary_weights.append(edge_weights[edge_idx])
        faces.append(oriented_edges)
        face_weights.append(float(np.mean(boundary_weights)))

    for i in range(n_side - 1):
        for j in range(n_side - 1):
            for k in range(n_side):
                add_face(
                    [
                        ("x", i, j, k, +1),
                        ("y", i + 1, j, k, +1),
                        ("x", i, j + 1, k, -1),
                        ("y", i, j, k, -1),
                    ]
                )
    for i in range(n_side - 1):
        for j in range(n_side):
            for k in range(n_side - 1):
                add_face(
                    [
                        ("x", i, j, k, +1),
                        ("z", i + 1, j, k, +1),
                        ("x", i, j, k + 1, -1),
                        ("z", i, j, k, -1),
                    ]
                )
    for i in range(n_side):
        for j in range(n_side - 1):
            for k in range(n_side - 1):
                add_face(
                    [
                        ("y", i, j, k, +1),
                        ("z", i, j + 1, k, +1),
                        ("y", i, j, k + 1, -1),
                        ("z", i, j, k, -1),
                    ]
                )

    n_faces = len(faces)
    C = np.zeros((n_faces, n_edges), dtype=float)
    for face_idx, boundary in enumerate(faces):
        for edge_idx, sign in boundary:
            C[face_idx, edge_idx] = float(sign)

    edge_weights_arr = np.asarray(edge_weights, dtype=float)
    face_weights_arr = np.asarray(face_weights, dtype=float)
    d0 = np.diag(np.sqrt(edge_weights_arr)) @ B0
    d1 = np.diag(np.sqrt(face_weights_arr)) @ C
    lower = d0 @ d0.T
    upper = d1.T @ d1
    L1 = lower + upper

    return {
        "points": points,
        "edges": edges,
        "directions": np.asarray(directions, dtype=float),
        "midpoints": np.asarray(midpoints, dtype=float),
        "edge_weights": edge_weights_arr,
        "face_weights": face_weights_arr,
        "B0": B0,
        "C": C,
        "d0": d0,
        "d1": d1,
        "lower": lower,
        "upper": upper,
        "L1": L1,
        "n_side": n_side,
    }


def classify_mode(vec: np.ndarray, eigenvalue: float, lower: np.ndarray, upper: np.ndarray, d0: np.ndarray, d1: np.ndarray) -> dict[str, float | str]:
    lam = max(float(eigenvalue), 1e-12)
    lower_energy = float(vec.T @ lower @ vec)
    upper_energy = float(vec.T @ upper @ vec)
    divergence_norm = float(np.linalg.norm(d0.T @ vec))
    curl_norm = float(np.linalg.norm(d1 @ vec))
    exact_fraction = lower_energy / lam
    coexact_fraction = upper_energy / lam

    if coexact_fraction > 0.65 and divergence_norm < curl_norm:
        label = "transverse-like"
    elif exact_fraction > 0.65 and curl_norm < divergence_norm:
        label = "gradient-like"
    else:
        label = "mixed"

    return {
        "label": label,
        "exact_fraction": exact_fraction,
        "coexact_fraction": coexact_fraction,
        "divergence_norm": divergence_norm,
        "curl_norm": curl_norm,
    }


def summarize_results(mode_records: list[dict[str, float | str]]) -> tuple[str, str]:
    transverse = sum(1 for rec in mode_records if rec["label"] == "transverse-like")
    gradient = sum(1 for rec in mode_records if rec["label"] == "gradient-like")
    if transverse >= 2:
        observation = "low edge spectrum contains coexact transverse-like modes"
        conclusion = "weighted L1 branch shows non-scalar candidate vector structure"
    elif gradient >= 2 and transverse == 0:
        observation = "low edge spectrum is dominated by exact gradient-like modes"
        conclusion = "edge branch is present but not yet separated into a clear vector sector"
    else:
        observation = "low edge spectrum mixes exact and coexact content"
        conclusion = "weighted L1 branch is active, but mode-family separation remains partial"
    return observation, conclusion


def make_plots(
    points: np.ndarray,
    midpoints: np.ndarray,
    directions: np.ndarray,
    evals: np.ndarray,
    evecs: np.ndarray,
    mode_records: list[dict[str, float | str]],
    name: str,
) -> list[str]:
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

    split_path = PLOTS / f"{name}_sector_split.png"
    fig, ax = plt.subplots(figsize=(8, 4))
    indices = np.arange(1, len(mode_records) + 1)
    exact = [float(rec["exact_fraction"]) for rec in mode_records]
    coexact = [float(rec["coexact_fraction"]) for rec in mode_records]
    ax.bar(indices - 0.18, exact, width=0.36, label="exact")
    ax.bar(indices + 0.18, coexact, width=0.36, label="coexact")
    ax.set_xlabel("mode index")
    ax.set_ylabel("fraction of eigenvalue")
    ax.set_title(f"{name} exact vs coexact split")
    ax.set_ylim(0.0, 1.1)
    ax.legend()
    ax.grid(alpha=0.2, axis="y")
    fig.savefig(split_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(split_path.relative_to(REPO_ROOT)))

    preferred = next((idx for idx, rec in enumerate(mode_records, start=1) if rec["label"] == "transverse-like"), 1)
    mode = evecs[:, preferred]
    mag = np.abs(mode)
    threshold = np.quantile(mag, 0.8) if len(mag) > 10 else 0.0
    select = mag >= threshold
    if not np.any(select):
        select = np.ones_like(mag, dtype=bool)
    max_mag = float(np.max(mag[select])) or 1.0
    scale = 0.8 / max(np.max(points[:, 0]) - np.min(points[:, 0]), 1.0)
    vectors = directions[select] * (mode[select] / max_mag)[:, None] * scale
    colors = cm.coolwarm(0.5 + 0.5 * np.clip(mode[select] / max_mag, -1.0, 1.0))

    mode_path = PLOTS / f"{name}_mode.png"
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(midpoints[select, 0], midpoints[select, 1], midpoints[select, 2], c=mode[select] / max_mag, cmap="coolwarm", s=20 + 180 * mag[select] / max_mag, vmin=-1.0, vmax=1.0)
    ax.quiver(
        midpoints[select, 0],
        midpoints[select, 1],
        midpoints[select, 2],
        vectors[:, 0],
        vectors[:, 1],
        vectors[:, 2],
        colors=colors,
        linewidth=0.8,
        arrow_length_ratio=0.25,
    )
    ax.set_title(f"{name} mode {preferred} ({mode_records[preferred - 1]['label']})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    fig.savefig(mode_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(mode_path.relative_to(REPO_ROOT)))

    return plot_paths


def run_hodge_test(
    config: dict[str, Any] | None = None,
    with_plots: bool = True,
    plot_name: str = "hodge_modes",
) -> dict[str, Any]:
    cfg = load_config(config)
    substrate = str(cfg.get("hodge_substrate", "cubic_lattice"))
    if substrate != "cubic_lattice":
        raise ValueError(f"Unsupported Hodge substrate: {substrate}")

    epsilon = float(cfg["epsilon"])
    n_side = int(cfg.get("hodge_lattice_side", cfg.get("gauge_lattice_side", 6)))
    complex_data = build_cubic_lattice_complex(n_side=n_side, epsilon=epsilon)
    L1 = np.asarray(complex_data["L1"], dtype=float)
    evals, evecs = np.linalg.eigh(L1)
    order = np.argsort(evals)
    evals = np.asarray(evals[order], dtype=float)
    evecs = np.asarray(evecs[:, order], dtype=float)

    mode_records = [
        {
            "mode_index": idx,
            "eigenvalue": float(evals[idx]),
            **classify_mode(
                evecs[:, idx],
                float(evals[idx]),
                np.asarray(complex_data["lower"], dtype=float),
                np.asarray(complex_data["upper"], dtype=float),
                np.asarray(complex_data["d0"], dtype=float),
                np.asarray(complex_data["d1"], dtype=float),
            ),
        }
        for idx in range(1, min(7, len(evals)))
    ]
    observation, conclusion = summarize_results(mode_records)
    plot_paths = (
        make_plots(
            np.asarray(complex_data["points"], dtype=float),
            np.asarray(complex_data["midpoints"], dtype=float),
            np.asarray(complex_data["directions"], dtype=float),
            evals,
            evecs,
            mode_records,
            plot_name,
        )
        if with_plots
        else []
    )

    return {
        "experiment": "hodge_l1_test",
        "config": {
            "nodes": int(len(np.asarray(complex_data["points"]))),
            "edges": int(len(np.asarray(complex_data["midpoints"]))),
            "faces": int(len(np.asarray(complex_data["face_weights"]))),
            "epsilon": epsilon,
            "random_seed": int(cfg.get("random_seed", 42)),
            "substrate": substrate,
            "hodge_lattice_side": n_side,
        },
        "spectrum": {
            "first_eigenvalues": [float(v) for v in evals[:12]],
            "mode_records": mode_records,
        },
        "plots": plot_paths,
        "observation": observation,
        "conclusion": conclusion,
    }


def main() -> None:
    result = run_hodge_test()
    print("epsilon_k =", result["config"]["epsilon"])
    print("first eigenvalues =", [round(v, 6) for v in result["spectrum"]["first_eigenvalues"][:8]])
    print("mode labels =", [rec["label"] for rec in result["spectrum"]["mode_records"]])
    print("plots =", result["plots"])
    print("observation =", result["observation"])
    print("conclusion =", result["conclusion"])


if __name__ == "__main__":
    main()
