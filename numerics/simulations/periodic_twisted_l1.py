#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import os
import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

PLOTS = REPO_ROOT / "plots"
PLOTS.mkdir(exist_ok=True)
RESULTS = REPO_ROOT / "data"
RESULTS.mkdir(exist_ok=True)
MPLCONFIG = REPO_ROOT / ".mplconfig"
MPLCONFIG.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

from numerics.simulations.hodge_modes import build_cubic_lattice_complex

DEFAULT_CONFIG: dict[str, Any] = {
    "epsilon": 0.2,
    "periodic_twisted_l1": {
        "sizes": [4, 5],
        "flux_quanta": [0, 1],
        "open_compare_side": 5,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged["periodic_twisted_l1"] = dict(DEFAULT_CONFIG["periodic_twisted_l1"])
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != "periodic_twisted_l1"})
        if isinstance(on_disk.get("periodic_twisted_l1"), dict):
            merged["periodic_twisted_l1"].update(on_disk["periodic_twisted_l1"])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != "periodic_twisted_l1"})
        if isinstance(config.get("periodic_twisted_l1"), dict):
            merged["periodic_twisted_l1"].update(config["periodic_twisted_l1"])
    return merged


def inverse_participation_ratio(vec: np.ndarray) -> float:
    abs_sq = np.abs(vec) ** 2
    norm_sq = float(np.sum(abs_sq))
    if norm_sq == 0.0:
        return 0.0
    return float(np.sum(abs_sq * abs_sq) / (norm_sq * norm_sq))


def build_periodic_twisted_complex(n_side: int, epsilon: float, flux_quanta: int) -> dict[str, Any]:
    node_index = np.arange(n_side**3).reshape((n_side, n_side, n_side))
    points = np.array(
        [[i / n_side, j / n_side, k / n_side] for i in range(n_side) for j in range(n_side) for k in range(n_side)],
        dtype=float,
    )

    h = 1.0 / n_side
    edge_weight = math.exp(-(h * h) / (2.0 * epsilon))
    edge_weights: list[float] = []
    edges: list[tuple[int, int]] = []
    edge_map: dict[tuple[str, int, int, int], int] = {}
    directions: list[np.ndarray] = []
    midpoints: list[np.ndarray] = []
    axes: list[str] = []
    link_phases: list[complex] = []

    phi = 2.0 * math.pi * flux_quanta / max(n_side * n_side, 1)

    def ux(i: int, j: int, k: int) -> complex:
        return np.exp(-1j * phi * n_side * j) if i == n_side - 1 else 1.0 + 0.0j

    def uy(i: int, j: int, k: int) -> complex:
        return np.exp(1j * phi * i)

    def uz(i: int, j: int, k: int) -> complex:
        return 1.0 + 0.0j

    phase_fn = {"x": ux, "y": uy, "z": uz}
    basis = {
        "x": np.array([1.0, 0.0, 0.0]),
        "y": np.array([0.0, 1.0, 0.0]),
        "z": np.array([0.0, 0.0, 1.0]),
    }

    def add_edge(axis: str, i: int, j: int, k: int, ni: int, nj: int, nk: int) -> None:
        u = int(node_index[i, j, k])
        v = int(node_index[ni, nj, nk])
        edge_idx = len(edges)
        edge_map[(axis, i, j, k)] = edge_idx
        edges.append((u, v))
        axes.append(axis)
        directions.append(basis[axis])
        midpoints.append(
            np.array(
                [
                    ((i + 0.5) / n_side) % 1.0 if axis == "x" else i / n_side,
                    ((j + 0.5) / n_side) % 1.0 if axis == "y" else j / n_side,
                    ((k + 0.5) / n_side) % 1.0 if axis == "z" else k / n_side,
                ],
                dtype=float,
            )
        )
        edge_weights.append(edge_weight)
        link_phases.append(phase_fn[axis](i, j, k))

    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side):
                add_edge("x", i, j, k, (i + 1) % n_side, j, k)
                add_edge("y", i, j, k, i, (j + 1) % n_side, k)
                add_edge("z", i, j, k, i, j, (k + 1) % n_side)

    n_nodes = len(points)
    n_edges = len(edges)
    edge_weights_arr = np.asarray(edge_weights, dtype=float)
    link_phase_arr = np.asarray(link_phases, dtype=complex)
    d0 = np.zeros((n_edges, n_nodes), dtype=complex)
    for edge_idx, (u, v) in enumerate(edges):
        weight = math.sqrt(edge_weights_arr[edge_idx])
        d0[edge_idx, u] = -weight
        d0[edge_idx, v] = weight * link_phase_arr[edge_idx]

    face_boundaries: list[list[tuple[str, int, int, int, int, complex]]] = []
    face_holonomies: list[complex] = []
    face_types: list[str] = []

    def add_face(face_type: str, boundary: list[tuple[str, int, int, int, int, complex]], holonomy: complex) -> None:
        face_boundaries.append(boundary)
        face_holonomies.append(holonomy)
        face_types.append(face_type)

    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side):
                add_face(
                    "xy",
                    [
                        ("x", i, j, k, +1, 1.0 + 0.0j),
                        ("y", (i + 1) % n_side, j, k, +1, np.conj(ux(i, j, k))),
                        ("x", i, (j + 1) % n_side, k, -1, np.conj(uy(i, j, k))),
                        ("y", i, j, k, -1, 1.0 + 0.0j),
                    ],
                    ux(i, j, k)
                    * uy((i + 1) % n_side, j, k)
                    * np.conj(ux(i, (j + 1) % n_side, k))
                    * np.conj(uy(i, j, k)),
                )
                add_face(
                    "xz",
                    [
                        ("x", i, j, k, +1, 1.0 + 0.0j),
                        ("z", (i + 1) % n_side, j, k, +1, np.conj(ux(i, j, k))),
                        ("x", i, j, (k + 1) % n_side, -1, 1.0 + 0.0j),
                        ("z", i, j, k, -1, 1.0 + 0.0j),
                    ],
                    1.0 + 0.0j,
                )
                add_face(
                    "yz",
                    [
                        ("y", i, j, k, +1, 1.0 + 0.0j),
                        ("z", i, (j + 1) % n_side, k, +1, np.conj(uy(i, j, k))),
                        ("y", i, j, (k + 1) % n_side, -1, 1.0 + 0.0j),
                        ("z", i, j, k, -1, 1.0 + 0.0j),
                    ],
                    1.0 + 0.0j,
                )

    face_weights = np.full(len(face_boundaries), edge_weight, dtype=float)
    d1 = np.zeros((len(face_boundaries), n_edges), dtype=complex)
    for face_idx, boundary in enumerate(face_boundaries):
        weight = math.sqrt(face_weights[face_idx])
        for axis, i, j, k, sign, phase in boundary:
            d1[face_idx, edge_map[(axis, i, j, k)]] = weight * sign * phase

    lower = d0 @ d0.conj().T
    upper = d1.conj().T @ d1
    L1 = lower + upper

    return {
        "points": points,
        "edges": edges,
        "midpoints": np.asarray(midpoints, dtype=float),
        "directions": np.asarray(directions, dtype=float),
        "axes": axes,
        "edge_weights": edge_weights_arr,
        "link_phases": link_phase_arr,
        "face_weights": face_weights,
        "face_types": face_types,
        "face_holonomies": np.asarray(face_holonomies, dtype=complex),
        "d0": d0,
        "d1": d1,
        "lower": lower,
        "upper": upper,
        "L1": L1,
        "n_side": n_side,
        "flux_quanta": flux_quanta,
        "plaquette_angle": phi,
    }


def orthonormal_image_basis(matrix: np.ndarray, tol: float = 1e-9) -> np.ndarray:
    u, singular_values, _ = np.linalg.svd(matrix, full_matrices=False)
    rank = int(np.sum(singular_values > tol))
    return u[:, :rank]


def nullspace_basis(matrix: np.ndarray, tol: float = 1e-9) -> np.ndarray:
    _, singular_values, vh = np.linalg.svd(matrix, full_matrices=True)
    rank = int(np.sum(singular_values > tol))
    return vh[rank:].conj().T


def analyze_mode(
    vec: np.ndarray,
    eigenvalue: float,
    exact_projector: np.ndarray,
    divfree_projector: np.ndarray,
    d0: np.ndarray,
    d1: np.ndarray,
    axes: list[str],
) -> dict[str, Any]:
    norm = float(np.real(np.vdot(vec, vec))) or 1.0
    exact_fraction = float(np.real(np.vdot(vec, exact_projector @ vec)) / norm)
    coexact_fraction = float(np.real(np.vdot(vec, divfree_projector @ vec)) / norm)
    divergence_norm = float(np.linalg.norm(d0.conj().T @ vec))
    curl_norm = float(np.linalg.norm(d1 @ vec))
    ipr = inverse_participation_ratio(vec)
    axis_energy = {
        axis: float(np.sum(np.abs(vec[[idx for idx, value in enumerate(axes) if value == axis]]) ** 2))
        for axis in ("x", "y", "z")
    }
    dominant_axis = max(axis_energy, key=axis_energy.get)

    if coexact_fraction > 0.95 and divergence_norm < 1e-8 and curl_norm < 1e-8:
        support = f"delocalized harmonic {dominant_axis}-cycle"
    elif coexact_fraction > 0.8 and curl_norm >= divergence_norm:
        support = f"divergence-free {dominant_axis}-biased circulation"
    elif exact_fraction > 0.7:
        support = f"delocalized exact {dominant_axis}-biased gradient"
    elif ipr > 0.06:
        support = f"localized {dominant_axis}-biased edge concentration"
    else:
        support = f"mixed {dominant_axis}-biased support"

    return {
        "eigenvalue": float(np.real(eigenvalue)),
        "exact_fraction": exact_fraction,
        "coexact_fraction": coexact_fraction,
        "divergence_norm": divergence_norm,
        "curl_norm": curl_norm,
        "ipr": ipr,
        "support_pattern": support,
    }


def infer_axes(directions: np.ndarray) -> list[str]:
    axes: list[str] = []
    for direction in directions:
        direction = np.asarray(direction, dtype=float)
        axis_index = int(np.argmax(np.abs(direction)))
        axes.append(["x", "y", "z"][axis_index])
    return axes


def analyze_complex(complex_data: dict[str, Any], label: str) -> dict[str, Any]:
    L1 = np.asarray(complex_data["L1"], dtype=complex)
    d0 = np.asarray(complex_data["d0"], dtype=complex)
    d1 = np.asarray(complex_data["d1"], dtype=complex)
    axes = list(complex_data.get("axes", infer_axes(np.asarray(complex_data["directions"], dtype=float))))

    evals, evecs = np.linalg.eigh(L1)
    order = np.argsort(evals.real)
    evals = np.asarray(evals[order], dtype=complex)
    evecs = np.asarray(evecs[:, order], dtype=complex)

    exact_basis = orthonormal_image_basis(d0)
    divfree_basis = nullspace_basis(d0.conj().T)
    exact_projector = exact_basis @ exact_basis.conj().T if exact_basis.size else np.zeros_like(L1)
    divfree_projector = divfree_basis @ divfree_basis.conj().T if divfree_basis.size else np.zeros_like(L1)

    projected_eval_list: list[float] = []
    projected_modes: list[dict[str, Any]] = []
    projected_vectors: list[np.ndarray] = []
    if divfree_basis.size:
        projected = divfree_basis.conj().T @ L1 @ divfree_basis
        projected_evals, projected_evecs = np.linalg.eigh(projected)
        projected_order = np.argsort(projected_evals.real)
        projected_evals = np.asarray(projected_evals[projected_order], dtype=complex)
        projected_evecs = np.asarray(projected_evecs[:, projected_order], dtype=complex)
        projected_eval_list = [float(np.real(value)) for value in projected_evals[:12]]
        for idx in range(min(8, projected_evecs.shape[1])):
            lifted = divfree_basis @ projected_evecs[:, idx]
            projected_vectors.append(lifted)
            projected_modes.append(
                {
                    "mode_index": idx,
                    **analyze_mode(lifted, float(np.real(projected_evals[idx])), exact_projector, divfree_projector, d0, d1, axes),
                }
            )

    full_modes = [
        {
            "mode_index": idx,
            **analyze_mode(evecs[:, idx], float(np.real(evals[idx])), exact_projector, divfree_projector, d0, d1, axes),
        }
        for idx in range(min(8, evecs.shape[1]))
    ]

    return {
        "label": label,
        "config": {
            "nodes": int(len(complex_data["points"])),
            "edges": int(len(complex_data["edges"])),
            "faces": int(len(complex_data["face_weights"])),
            "n_side": int(complex_data["n_side"]),
            "flux_quanta": int(complex_data.get("flux_quanta", 0)),
            "plaquette_angle": float(complex_data.get("plaquette_angle", 0.0)),
        },
        "full_spectrum": [float(np.real(value)) for value in evals[:12]],
        "projected_spectrum": projected_eval_list,
        "full_modes": full_modes,
        "projected_modes": projected_modes,
        "points": np.asarray(complex_data["points"], dtype=float),
        "midpoints": np.asarray(complex_data["midpoints"], dtype=float),
        "directions": np.asarray(complex_data["directions"], dtype=float),
        "full_vectors": [evecs[:, idx] for idx in range(min(8, evecs.shape[1]))],
        "projected_vectors": projected_vectors,
        "face_holonomy_phase": float(np.angle(np.mean(np.asarray(complex_data.get("face_holonomies", np.array([1.0])))[0:1]))),
    }


def choose_representative_vector(case: dict[str, Any], prefer_projected: bool = False) -> np.ndarray:
    vectors_key = "projected_vectors" if prefer_projected and case["projected_vectors"] else "full_vectors"
    modes_key = "projected_modes" if prefer_projected and case["projected_modes"] else "full_modes"
    vectors = case[vectors_key]
    modes = case[modes_key]
    if not vectors:
        return case["full_vectors"][0]
    preferred_idx = 0
    for idx, record in enumerate(modes):
        if record["coexact_fraction"] > 0.8:
            preferred_idx = idx
            break
    return np.asarray(vectors[preferred_idx], dtype=complex)


def plot_edge_mode(ax: Any, midpoints: np.ndarray, directions: np.ndarray, mode: np.ndarray, title: str) -> None:
    magnitudes = np.abs(mode)
    threshold = np.quantile(magnitudes, 0.82) if len(magnitudes) > 10 else 0.0
    select = magnitudes >= threshold
    if not np.any(select):
        select = np.ones_like(magnitudes, dtype=bool)
    max_mag = float(np.max(magnitudes[select])) or 1.0
    scale = 0.18
    vectors = directions[select] * (np.real(mode[select]) / max_mag)[:, None] * scale
    phases = np.angle(mode[select])
    colors = cm.twilight((phases + math.pi) / (2.0 * math.pi))
    ax.scatter(
        midpoints[select, 0],
        midpoints[select, 1],
        midpoints[select, 2],
        c=phases,
        cmap="twilight",
        s=20 + 180 * magnitudes[select] / max_mag,
        vmin=-math.pi,
        vmax=math.pi,
    )
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
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")


def make_plots(cases: dict[str, Any], open_label: str, plot_name: str) -> list[str]:
    plot_paths: list[str] = []
    periodic_labels = [label for label in cases if label.startswith("periodic")]
    periodic_labels.sort()
    largest_periodic = max(
        periodic_labels,
        key=lambda label: (cases[label]["config"]["n_side"], cases[label]["config"]["flux_quanta"]),
    )
    trivial_label = next(label for label in periodic_labels if cases[label]["config"]["flux_quanta"] == 0 and cases[label]["config"]["n_side"] == cases[largest_periodic]["config"]["n_side"])
    flux_label = next(label for label in periodic_labels if cases[label]["config"]["flux_quanta"] > 0 and cases[label]["config"]["n_side"] == cases[largest_periodic]["config"]["n_side"])
    flux_plot = PLOTS / f"{plot_name}_spectral_flow.png"
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharex=True)
    for n_side in sorted({cases[label]["config"]["n_side"] for label in periodic_labels}):
        labels = [label for label in periodic_labels if cases[label]["config"]["n_side"] == n_side]
        labels.sort(key=lambda label: cases[label]["config"]["flux_quanta"])
        flux_values = [cases[label]["config"]["plaquette_angle"] for label in labels]
        for mode_index in range(4):
            axes[0].plot(
                flux_values,
                [cases[label]["full_spectrum"][mode_index] for label in labels],
                marker="o",
                label=f"N={n_side**3}, mode {mode_index}" if mode_index == 0 else None,
                alpha=0.8,
            )
            axes[1].plot(
                flux_values,
                [cases[label]["projected_spectrum"][mode_index] for label in labels],
                marker="o",
                label=f"N={n_side**3}, proj {mode_index}" if mode_index == 0 else None,
                alpha=0.8,
            )
    axes[0].set_title("Full L1 spectral flow")
    axes[1].set_title("Coexact-projected spectral flow")
    for ax in axes:
        ax.set_xlabel("plaquette angle")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("eigenvalue")
    axes[1].legend(fontsize=8)
    fig.savefig(flux_plot, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(flux_plot.relative_to(REPO_ROOT)))

    compare_plot = PLOTS / f"{plot_name}_open_vs_periodic.png"
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(range(8), cases[open_label]["full_spectrum"][:8], marker="o", label=f"open, N={cases[open_label]['config']['nodes']}")
    ax.plot(range(8), cases[trivial_label]["full_spectrum"][:8], marker="o", label=f"periodic flux=0, N={cases[trivial_label]['config']['nodes']}")
    ax.plot(range(8), cases[flux_label]["full_spectrum"][:8], marker="o", label=f"periodic flux>0, N={cases[flux_label]['config']['nodes']}")
    ax.set_xlabel("mode index")
    ax.set_ylabel("eigenvalue")
    ax.set_title("Open vs periodic low L1 spectrum")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(compare_plot, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(compare_plot.relative_to(REPO_ROOT)))

    fraction_plot = PLOTS / f"{plot_name}_fractions.png"
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    for ax, label in zip(axes, [trivial_label, flux_label]):
        records = cases[label]["full_modes"][:6]
        indices = np.arange(len(records))
        exact = [record["exact_fraction"] for record in records]
        coexact = [record["coexact_fraction"] for record in records]
        ax.bar(indices - 0.18, exact, width=0.36, label="exact")
        ax.bar(indices + 0.18, coexact, width=0.36, label="coexact/div-free")
        ax.set_title(f"N={cases[label]['config']['nodes']}, flux quanta={cases[label]['config']['flux_quanta']}")
        ax.set_xlabel("mode index")
        ax.grid(alpha=0.2, axis="y")
    axes[0].set_ylabel("fraction")
    axes[0].legend(fontsize=8)
    fig.savefig(fraction_plot, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(fraction_plot.relative_to(REPO_ROOT)))

    divergence_plot = PLOTS / f"{plot_name}_divergence.png"
    fig, ax = plt.subplots(figsize=(8, 4))
    for label in [open_label, trivial_label, flux_label]:
        records = cases[label]["full_modes"][:6]
        ax.plot(range(len(records)), [record["divergence_norm"] for record in records], marker="o", label=label)
    ax.set_xlabel("mode index")
    ax.set_ylabel(r"$||d_0^* a||$")
    ax.set_title("Divergence norm by mode")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(divergence_plot, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(divergence_plot.relative_to(REPO_ROOT)))

    curl_plot = PLOTS / f"{plot_name}_curl.png"
    fig, ax = plt.subplots(figsize=(8, 4))
    for label in [open_label, trivial_label, flux_label]:
        records = cases[label]["full_modes"][:6]
        ax.plot(range(len(records)), [record["curl_norm"] for record in records], marker="o", label=label)
    ax.set_xlabel("mode index")
    ax.set_ylabel(r"$||d_1 a||$")
    ax.set_title("Curl norm by mode")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(curl_plot, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(curl_plot.relative_to(REPO_ROOT)))

    mode_plot = PLOTS / f"{plot_name}_modes.png"
    fig = plt.figure(figsize=(14, 10))
    panels = [
        (open_label, False, "Open box low mode"),
        (trivial_label, False, "Periodic flux=0 low mode"),
        (flux_label, False, "Periodic flux>0 low mode"),
        (flux_label, True, "Periodic flux>0 projected mode"),
    ]
    for panel_index, (label, projected, title) in enumerate(panels, start=1):
        ax = fig.add_subplot(2, 2, panel_index, projection="3d")
        vector = choose_representative_vector(cases[label], prefer_projected=projected)
        plot_edge_mode(ax, cases[label]["midpoints"], cases[label]["directions"], vector, title)
    fig.savefig(mode_plot, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(mode_plot.relative_to(REPO_ROOT)))

    return plot_paths


def summarize_verdict(cases: dict[str, Any], sizes: list[int], flux_quanta: list[int], open_label: str) -> tuple[str, str, dict[str, bool]]:
    periodic_cases = [case for label, case in cases.items() if label != open_label]
    open_case = cases[open_label]
    open_exact_dominated = all(record["exact_fraction"] > 0.6 for record in open_case["full_modes"][:3])
    coexact_low_family = True
    for n_side in sizes:
        matching = [case for case in periodic_cases if case["config"]["n_side"] == n_side]
        if not matching:
            coexact_low_family = False
            continue
        for case in matching:
            low_records = case["full_modes"][:2]
            if not all(record["coexact_fraction"] > 0.75 for record in low_records):
                coexact_low_family = False
    flux_sensitive = False
    positive_fluxes = [q for q in flux_quanta if q > 0]
    if positive_fluxes:
        target_flux = positive_fluxes[0]
        for n_side in sizes:
            zero_case = next(case for case in periodic_cases if case["config"]["n_side"] == n_side and case["config"]["flux_quanta"] == 0)
            flux_case = next(case for case in periodic_cases if case["config"]["n_side"] == n_side and case["config"]["flux_quanta"] == target_flux)
            shift = abs(flux_case["full_spectrum"][0] - zero_case["full_spectrum"][0])
            if shift > 0.01:
                flux_sensitive = True
    persists_across_sizes = True
    for n_side in sizes:
        matching = [case for case in periodic_cases if case["config"]["n_side"] == n_side and case["config"]["flux_quanta"] > 0]
        if not matching or matching[0]["full_modes"][0]["coexact_fraction"] <= 0.75:
            persists_across_sizes = False
    vector_like = coexact_low_family and flux_sensitive and persists_across_sizes and open_exact_dominated

    if vector_like:
        observation = "periodic topology produces a low divergence-free family absent in the open box, and nontrivial flux shifts it cleanly"
        conclusion = "the periodic/twisted L1 branch shows a robust non-scalar flux-sensitive edge sector, but it is not yet a Maxwell-like propagating family"
    elif coexact_low_family and flux_sensitive:
        observation = "periodic topology produces a mostly coexact low family, and flux moves it, but the family is not cleanly isolated across all cases"
        conclusion = "the edge sector is more than scalar leakage, but the current branch is still partially contaminated by harmonic or mixed modes"
    else:
        observation = "periodic and twisted L1 does not isolate a stable low coexact family in the tested range"
        conclusion = "the current gauge story still looks like scalar contamination, and the operator needs a deeper change"

    verdict = {
        "low_family_mostly_coexact": coexact_low_family,
        "responds_cleanly_to_flux": flux_sensitive,
        "persists_across_sizes": persists_across_sizes,
        "plausibly_vector_like": vector_like,
    }
    return observation, conclusion, verdict


def sanitize_case(case: dict[str, Any]) -> dict[str, Any]:
    return {
        "label": case["label"],
        "config": case["config"],
        "face_holonomy_phase": case["face_holonomy_phase"],
        "full_spectrum": case["full_spectrum"],
        "projected_spectrum": case["projected_spectrum"],
        "full_modes": case["full_modes"],
        "projected_modes": case["projected_modes"],
    }


def save_results(result: dict[str, Any], plot_paths: list[str]) -> tuple[Path, list[str], str]:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output = dict(result)
    output["plots"] = plot_paths
    payload = json.dumps(output, indent=2)
    stamped = RESULTS / f"{timestamp}_periodic_twisted_l1.json"
    latest = RESULTS / "periodic_twisted_l1_latest.json"
    stamped.write_text(payload, encoding="utf-8")
    latest.write_text(payload, encoding="utf-8")

    stamped_plots: list[str] = []
    for rel_path in plot_paths:
        src = REPO_ROOT / rel_path
        dst = PLOTS / f"{timestamp}_periodic_twisted_l1_{src.name}"
        shutil.copy2(src, dst)
        stamped_plots.append(str(dst.relative_to(REPO_ROOT)))
    return stamped, stamped_plots, timestamp


def run_periodic_twisted_l1_experiment(
    config: dict[str, Any] | None = None,
    config_path: Path | None = None,
    plot_name: str = "periodic_twisted_l1",
) -> tuple[dict[str, Any], Path, list[str], str]:
    cfg = load_config(config, config_path=config_path)
    epsilon = float(cfg.get("epsilon", 0.2))
    experiment_cfg = cfg["periodic_twisted_l1"]
    sizes = [int(value) for value in experiment_cfg.get("sizes", [4, 5])]
    flux_quanta = [int(value) for value in experiment_cfg.get("flux_quanta", [0, 1])]
    open_compare_side = int(experiment_cfg.get("open_compare_side", max(sizes)))

    cases: dict[str, Any] = {}

    open_complex = build_cubic_lattice_complex(open_compare_side, epsilon)
    open_complex["flux_quanta"] = 0
    open_complex["plaquette_angle"] = 0.0
    open_complex["face_holonomies"] = np.ones(len(open_complex["face_weights"]), dtype=complex)
    open_case_label = f"open_n{open_compare_side}"
    cases[open_case_label] = analyze_complex(open_complex, open_case_label)

    for n_side in sizes:
        for flux in flux_quanta:
            label = f"periodic_n{n_side}_flux{flux}"
            periodic_complex = build_periodic_twisted_complex(n_side=n_side, epsilon=epsilon, flux_quanta=flux)
            cases[label] = analyze_complex(periodic_complex, label)

    observation, conclusion, verdict = summarize_verdict(cases, sizes, flux_quanta, open_case_label)
    result = {
        "experiment": "periodic_twisted_l1",
        "config": {
            "epsilon": epsilon,
            "sizes": sizes,
            "flux_quanta": flux_quanta,
            "open_compare_side": open_compare_side,
        },
        "open_case": open_case_label,
        "cases": {label: sanitize_case(case) for label, case in cases.items()},
        "verdict": verdict,
        "observation": observation,
        "conclusion": conclusion,
    }
    plot_paths = make_plots(cases=cases, open_label=open_case_label, plot_name=plot_name)
    result_path, stamped_plots, timestamp = save_results(result, plot_paths)
    return result, result_path, stamped_plots, timestamp


def main() -> None:
    result, result_path, stamped_plots, _ = run_periodic_twisted_l1_experiment()
    print(f"Saved: {result_path.relative_to(REPO_ROOT)}")
    print(f"Plots: {', '.join(stamped_plots)}")
    print(f"Observation: {result['observation']}")
    print(f"Conclusion: {result['conclusion']}")
    print(f"Verdict: {result['verdict']}")


if __name__ == "__main__":
    main()
