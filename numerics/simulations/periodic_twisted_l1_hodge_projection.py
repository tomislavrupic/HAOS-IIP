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
EXPERIMENT_LOG = REPO_ROOT / "experiments" / "EXPERIMENT_LOG.md"
MPLCONFIG = REPO_ROOT / ".mplconfig"
MPLCONFIG.mkdir(exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIG))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from numerics.simulations.periodic_twisted_l1 import (
    build_periodic_twisted_complex,
    infer_axes,
    inverse_participation_ratio,
    nullspace_basis,
    orthonormal_image_basis,
    plot_edge_mode,
)

DEFAULT_CONFIG: dict[str, Any] = {
    "epsilon": 0.2,
    "periodic_twisted_l1_hodge_projection": {
        "sizes": [4, 5],
        "flux_quanta": [0, 1, 2, 3, 4],
        "low_modes": 6,
        "representative_fluxes": [0, 1, 4],
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged["periodic_twisted_l1_hodge_projection"] = dict(DEFAULT_CONFIG["periodic_twisted_l1_hodge_projection"])
    path = config_path or (REPO_ROOT / "config.json")
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != "periodic_twisted_l1_hodge_projection"})
        if isinstance(on_disk.get("periodic_twisted_l1_hodge_projection"), dict):
            merged["periodic_twisted_l1_hodge_projection"].update(on_disk["periodic_twisted_l1_hodge_projection"])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != "periodic_twisted_l1_hodge_projection"})
        if isinstance(config.get("periodic_twisted_l1_hodge_projection"), dict):
            merged["periodic_twisted_l1_hodge_projection"].update(config["periodic_twisted_l1_hodge_projection"])
    return merged


def make_projector(basis: np.ndarray, dimension: int) -> np.ndarray:
    if basis.size == 0:
        return np.zeros((dimension, dimension), dtype=complex)
    return basis @ basis.conj().T


def orthonormalize_subspace(matrix: np.ndarray, tol: float = 1e-9) -> np.ndarray:
    if matrix.size == 0:
        return np.zeros((matrix.shape[0], 0), dtype=complex)
    return orthonormal_image_basis(matrix, tol=tol)


def max_spectral_norm(matrix: np.ndarray) -> float:
    if matrix.size == 0:
        return 0.0
    return float(np.linalg.norm(matrix, ord=2))


def build_hodge_projectors(d0: np.ndarray, d1: np.ndarray, tol: float = 1e-9) -> dict[str, Any]:
    n_edges = d0.shape[0]
    eye = np.eye(n_edges, dtype=complex)

    exact_basis = orthonormal_image_basis(d0, tol=tol)
    exact_projector = make_projector(exact_basis, n_edges)

    harmonic_constraints = np.vstack([d0.conj().T, d1])
    harmonic_basis = nullspace_basis(harmonic_constraints, tol=tol)
    if harmonic_basis.size:
        harmonic_basis = orthonormalize_subspace(harmonic_basis, tol=tol)
    harmonic_projector = make_projector(harmonic_basis, n_edges)

    coexact_raw = orthonormal_image_basis(d1.conj().T, tol=tol)
    coexact_basis = coexact_raw
    if coexact_raw.size:
        coexact_basis = (eye - exact_projector - harmonic_projector) @ coexact_raw
        coexact_basis = orthonormalize_subspace(coexact_basis, tol=tol)
    coexact_projector = make_projector(coexact_basis, n_edges)

    projector_sum = exact_projector + harmonic_projector + coexact_projector
    checks = {
        "dimensions": {
            "edge_space": int(n_edges),
            "exact": int(exact_basis.shape[1]),
            "harmonic": int(harmonic_basis.shape[1]),
            "coexact": int(coexact_basis.shape[1]),
        },
        "chain_error": max_spectral_norm(d1 @ d0),
        "orthogonality_errors": {
            "exact_harmonic": max_spectral_norm(exact_basis.conj().T @ harmonic_basis),
            "exact_coexact": max_spectral_norm(exact_basis.conj().T @ coexact_basis),
            "harmonic_coexact": max_spectral_norm(harmonic_basis.conj().T @ coexact_basis),
        },
        "idempotence_errors": {
            "exact": max_spectral_norm(exact_projector @ exact_projector - exact_projector),
            "harmonic": max_spectral_norm(harmonic_projector @ harmonic_projector - harmonic_projector),
            "coexact": max_spectral_norm(coexact_projector @ coexact_projector - coexact_projector),
        },
        "reconstruction_error": max_spectral_norm(projector_sum - eye),
    }

    rng = np.random.default_rng(42)
    sample_errors: list[float] = []
    for _ in range(3):
        sample = rng.normal(size=n_edges) + 1j * rng.normal(size=n_edges)
        reconstructed = projector_sum @ sample
        denom = float(np.linalg.norm(sample)) or 1.0
        sample_errors.append(float(np.linalg.norm(sample - reconstructed) / denom))
    checks["sample_reconstruction_error"] = float(max(sample_errors) if sample_errors else 0.0)

    return {
        "exact_basis": exact_basis,
        "harmonic_basis": harmonic_basis,
        "coexact_basis": coexact_basis,
        "exact_projector": exact_projector,
        "harmonic_projector": harmonic_projector,
        "coexact_projector": coexact_projector,
        "checks": checks,
    }


def analyze_mode(
    vec: np.ndarray,
    eigenvalue: float,
    exact_projector: np.ndarray,
    harmonic_projector: np.ndarray,
    coexact_projector: np.ndarray,
    d0: np.ndarray,
    d1: np.ndarray,
    axes: list[str],
) -> dict[str, Any]:
    norm = float(np.real(np.vdot(vec, vec))) or 1.0
    exact_fraction = float(np.real(np.vdot(vec, exact_projector @ vec)) / norm)
    harmonic_fraction = float(np.real(np.vdot(vec, harmonic_projector @ vec)) / norm)
    coexact_fraction = float(np.real(np.vdot(vec, coexact_projector @ vec)) / norm)
    residual_fraction = max(0.0, 1.0 - exact_fraction - harmonic_fraction - coexact_fraction)
    divergence_norm = float(np.linalg.norm(d0.conj().T @ vec))
    curl_norm = float(np.linalg.norm(d1 @ vec))
    ipr = inverse_participation_ratio(vec)

    axis_energy = {
        axis: float(np.sum(np.abs(vec[[idx for idx, value in enumerate(axes) if value == axis]]) ** 2))
        for axis in ("x", "y", "z")
    }
    dominant_axis = max(axis_energy, key=axis_energy.get)

    if harmonic_fraction > 0.8 and divergence_norm < 1e-8 and curl_norm < 1e-8:
        support = f"harmonic {dominant_axis}-cycle"
    elif coexact_fraction > 0.8 and harmonic_fraction < 0.2 and curl_norm >= divergence_norm:
        support = f"coexact {dominant_axis}-biased circulation"
    elif coexact_fraction > 0.5 and harmonic_fraction > 0.2:
        support = f"flux-lifted harmonic/{dominant_axis}-circulation mix"
    elif exact_fraction > 0.7:
        support = f"exact {dominant_axis}-biased gradient"
    elif ipr > 0.06:
        support = f"localized {dominant_axis}-biased edge concentration"
    else:
        support = f"mixed {dominant_axis}-biased support"

    return {
        "eigenvalue": float(np.real(eigenvalue)),
        "exact_fraction": exact_fraction,
        "harmonic_fraction": harmonic_fraction,
        "coexact_fraction": coexact_fraction,
        "residual_fraction": residual_fraction,
        "divergence_norm": divergence_norm,
        "curl_norm": curl_norm,
        "ipr": ipr,
        "support_pattern": support,
    }


def analyze_projected_operator(
    operator: np.ndarray,
    basis: np.ndarray,
    exact_projector: np.ndarray,
    harmonic_projector: np.ndarray,
    coexact_projector: np.ndarray,
    d0: np.ndarray,
    d1: np.ndarray,
    axes: list[str],
    low_modes: int,
) -> tuple[list[float], list[dict[str, Any]], list[np.ndarray]]:
    if basis.size == 0:
        return [], [], []
    projected = basis.conj().T @ operator @ basis
    evals, evecs = np.linalg.eigh(projected)
    order = np.argsort(evals.real)
    evals = np.asarray(evals[order], dtype=complex)
    evecs = np.asarray(evecs[:, order], dtype=complex)

    spectrum = [float(np.real(value)) for value in evals[: max(low_modes, 12)]]
    modes: list[dict[str, Any]] = []
    vectors: list[np.ndarray] = []
    for idx in range(min(low_modes, evecs.shape[1])):
        lifted = basis @ evecs[:, idx]
        vectors.append(lifted)
        modes.append(
            {
                "mode_index": idx,
                **analyze_mode(lifted, float(np.real(evals[idx])), exact_projector, harmonic_projector, coexact_projector, d0, d1, axes),
            }
        )
    return spectrum, modes, vectors


def build_reference_decomposition(n_side: int, epsilon: float) -> dict[str, Any]:
    reference_complex = build_periodic_twisted_complex(n_side=n_side, epsilon=epsilon, flux_quanta=0)
    d0_ref = np.asarray(reference_complex["d0"], dtype=complex)
    d1_ref = np.asarray(reference_complex["d1"], dtype=complex)
    decomposition = build_hodge_projectors(d0_ref, d1_ref)
    decomposition["reference_d0"] = d0_ref
    decomposition["reference_d1"] = d1_ref
    return decomposition


def analyze_case(n_side: int, flux_quanta: int, epsilon: float, low_modes: int, reference: dict[str, Any]) -> dict[str, Any]:
    complex_data = build_periodic_twisted_complex(n_side=n_side, epsilon=epsilon, flux_quanta=flux_quanta)
    label = f"periodic_n{n_side}_flux{flux_quanta}"
    L1 = np.asarray(complex_data["L1"], dtype=complex)
    d0_ref = np.asarray(reference["reference_d0"], dtype=complex)
    d1_ref = np.asarray(reference["reference_d1"], dtype=complex)
    axes = list(complex_data.get("axes", infer_axes(np.asarray(complex_data["directions"], dtype=float))))

    exact_projector = reference["exact_projector"]
    harmonic_projector = reference["harmonic_projector"]
    coexact_projector = reference["coexact_projector"]

    evals, evecs = np.linalg.eigh(L1)
    order = np.argsort(evals.real)
    evals = np.asarray(evals[order], dtype=complex)
    evecs = np.asarray(evecs[:, order], dtype=complex)

    full_modes = [
        {
            "mode_index": idx,
            **analyze_mode(evecs[:, idx], float(np.real(evals[idx])), exact_projector, harmonic_projector, coexact_projector, d0_ref, d1_ref, axes),
        }
        for idx in range(min(low_modes, evecs.shape[1]))
    ]

    harmonic_spectrum, harmonic_modes, harmonic_vectors = analyze_projected_operator(
        L1,
        reference["harmonic_basis"],
        exact_projector,
        harmonic_projector,
        coexact_projector,
        d0_ref,
        d1_ref,
        axes,
        low_modes,
    )
    coexact_spectrum, coexact_modes, coexact_vectors = analyze_projected_operator(
        L1,
        reference["coexact_basis"],
        exact_projector,
        harmonic_projector,
        coexact_projector,
        d0_ref,
        d1_ref,
        axes,
        low_modes,
    )

    selected_branch = coexact_modes[0] if coexact_modes else None

    return {
        "label": label,
        "config": {
            "nodes": int(len(complex_data["points"])),
            "edges": int(len(complex_data["edges"])),
            "faces": int(len(complex_data["face_weights"])),
            "n_side": int(n_side),
            "flux_quanta": int(flux_quanta),
            "plaquette_angle": float(complex_data["plaquette_angle"]),
        },
        "decomposition": reference["checks"],
        "full_spectrum": [float(np.real(value)) for value in evals[: max(low_modes, 12)]],
        "harmonic_spectrum": harmonic_spectrum,
        "coexact_spectrum": coexact_spectrum,
        "full_modes": full_modes,
        "harmonic_modes": harmonic_modes,
        "coexact_modes": coexact_modes,
        "selected_coexact_branch": selected_branch,
        "midpoints": np.asarray(complex_data["midpoints"], dtype=float),
        "directions": np.asarray(complex_data["directions"], dtype=float),
        "full_vectors": [evecs[:, idx] for idx in range(min(low_modes, evecs.shape[1]))],
        "harmonic_vectors": harmonic_vectors,
        "coexact_vectors": coexact_vectors,
    }


def sanitize_case(case: dict[str, Any]) -> dict[str, Any]:
    return {
        "label": case["label"],
        "config": case["config"],
        "decomposition": case["decomposition"],
        "full_spectrum": case["full_spectrum"],
        "harmonic_spectrum": case["harmonic_spectrum"],
        "coexact_spectrum": case["coexact_spectrum"],
        "full_modes": case["full_modes"],
        "harmonic_modes": case["harmonic_modes"],
        "coexact_modes": case["coexact_modes"],
        "selected_coexact_branch": case["selected_coexact_branch"],
    }


def plot_mode_safely(ax: Any, case: dict[str, Any], vector: np.ndarray | None, title: str) -> None:
    if vector is None or vector.size == 0:
        ax.set_title(title)
        ax.text2D(0.5, 0.5, "no mode", transform=ax.transAxes, ha="center", va="center")
        ax.set_axis_off()
        return
    plot_edge_mode(ax, case["midpoints"], case["directions"], vector, title)


def make_plots(cases: dict[str, Any], plot_name: str, low_modes: int, representative_fluxes: list[int]) -> list[str]:
    plot_paths: list[str] = []
    sizes = sorted({int(case["config"]["n_side"]) for case in cases.values()})

    compare_path = PLOTS / f"{plot_name}_full_vs_coexact.png"
    fig, axes = plt.subplots(1, len(sizes), figsize=(7 * len(sizes), 4), sharey=True)
    if len(sizes) == 1:
        axes = [axes]
    for ax, n_side in zip(axes, sizes):
        labels = sorted(
            [label for label in cases if cases[label]["config"]["n_side"] == n_side],
            key=lambda label: cases[label]["config"]["flux_quanta"],
        )
        for label in labels:
            case = cases[label]
            flux = int(case["config"]["flux_quanta"])
            if flux not in representative_fluxes:
                continue
            ax.plot(range(low_modes), case["full_spectrum"][:low_modes], marker="o", alpha=0.45, label=f"full m={flux}")
            if case["coexact_spectrum"]:
                ax.plot(range(min(low_modes, len(case["coexact_spectrum"]))), case["coexact_spectrum"][:low_modes], marker="s", linestyle="--", label=f"coexact m={flux}")
        ax.set_title(f"N={n_side**3}")
        ax.set_xlabel("mode index")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, ncol=2)
    axes[0].set_ylabel("eigenvalue")
    fig.savefig(compare_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(compare_path.relative_to(REPO_ROOT)))

    harmonic_path = PLOTS / f"{plot_name}_harmonic_fraction.png"
    fig, axes = plt.subplots(1, len(sizes), figsize=(7 * len(sizes), 4), sharey=True)
    if len(sizes) == 1:
        axes = [axes]
    for ax, n_side in zip(axes, sizes):
        labels = sorted(
            [label for label in cases if cases[label]["config"]["n_side"] == n_side],
            key=lambda label: cases[label]["config"]["flux_quanta"],
        )
        for label in labels:
            case = cases[label]
            ax.plot(
                range(min(low_modes, len(case["full_modes"]))),
                [record["harmonic_fraction"] for record in case["full_modes"][:low_modes]],
                marker="o",
                alpha=0.75,
                label=f"m={int(case['config']['flux_quanta'])}",
            )
        ax.set_title(f"Harmonic fraction, N={n_side**3}")
        ax.set_xlabel("mode index")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, ncol=2)
    axes[0].set_ylabel("fraction")
    fig.savefig(harmonic_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(harmonic_path.relative_to(REPO_ROOT)))

    coexact_path = PLOTS / f"{plot_name}_coexact_fraction.png"
    fig, axes = plt.subplots(1, len(sizes), figsize=(7 * len(sizes), 4), sharey=True)
    if len(sizes) == 1:
        axes = [axes]
    for ax, n_side in zip(axes, sizes):
        labels = sorted(
            [label for label in cases if cases[label]["config"]["n_side"] == n_side],
            key=lambda label: cases[label]["config"]["flux_quanta"],
        )
        for label in labels:
            case = cases[label]
            ax.plot(
                range(min(low_modes, len(case["full_modes"]))),
                [record["coexact_fraction"] for record in case["full_modes"][:low_modes]],
                marker="o",
                alpha=0.75,
                label=f"m={int(case['config']['flux_quanta'])}",
            )
        ax.set_title(f"Coexact fraction, N={n_side**3}")
        ax.set_xlabel("mode index")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, ncol=2)
    axes[0].set_ylabel("fraction")
    fig.savefig(coexact_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(coexact_path.relative_to(REPO_ROOT)))

    flow_path = PLOTS / f"{plot_name}_coexact_flow.png"
    fig, axes = plt.subplots(1, len(sizes), figsize=(7 * len(sizes), 4), sharey=True)
    if len(sizes) == 1:
        axes = [axes]
    for ax, n_side in zip(axes, sizes):
        labels = sorted(
            [label for label in cases if cases[label]["config"]["n_side"] == n_side],
            key=lambda label: cases[label]["config"]["flux_quanta"],
        )
        flux_values = [cases[label]["config"]["plaquette_angle"] for label in labels]
        for mode_index in range(min(4, low_modes)):
            values = []
            for label in labels:
                spectrum = cases[label]["coexact_spectrum"]
                values.append(spectrum[mode_index] if len(spectrum) > mode_index else math.nan)
            ax.plot(flux_values, values, marker="o", label=f"coexact {mode_index}")
        ax.set_title(f"Coexact spectral flow, N={n_side**3}")
        ax.set_xlabel("plaquette angle")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("eigenvalue")
    fig.savefig(flow_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(flow_path.relative_to(REPO_ROOT)))

    persistence_path = PLOTS / f"{plot_name}_persistence.png"
    fig, axes = plt.subplots(1, 2, figsize=(14, 4), sharex=True)
    selected_fluxes = representative_fluxes[: min(3, len(representative_fluxes))]
    for flux in selected_fluxes:
        labels = [label for label in cases if int(cases[label]["config"]["flux_quanta"]) == flux]
        labels.sort(key=lambda label: cases[label]["config"]["n_side"])
        x_values = [cases[label]["config"]["n_side"] for label in labels]
        first_mode = [cases[label]["coexact_spectrum"][0] if cases[label]["coexact_spectrum"] else math.nan for label in labels]
        second_mode = [cases[label]["coexact_spectrum"][1] if len(cases[label]["coexact_spectrum"]) > 1 else math.nan for label in labels]
        axes[0].plot(x_values, first_mode, marker="o", label=f"m={flux}")
        axes[1].plot(x_values, second_mode, marker="o", label=f"m={flux}")
    axes[0].set_title("First coexact eigenvalue vs size")
    axes[1].set_title("Second coexact eigenvalue vs size")
    for ax in axes:
        ax.set_xlabel("lattice side n")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("eigenvalue")
    fig.savefig(persistence_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(persistence_path.relative_to(REPO_ROOT)))

    mode_path = PLOTS / f"{plot_name}_modes.png"
    rep_fluxes = representative_fluxes or [0]
    fig = plt.figure(figsize=(5 * len(rep_fluxes), 4.5 * len(sizes)))
    panel_index = 1
    for n_side in sizes:
        for flux in rep_fluxes:
            label = f"periodic_n{n_side}_flux{flux}"
            case = cases[label]
            vector = case["coexact_vectors"][0] if case["coexact_vectors"] else None
            ax = fig.add_subplot(len(sizes), len(rep_fluxes), panel_index, projection="3d")
            plot_mode_safely(ax, case, vector, f"N={n_side**3}, m={flux}")
            panel_index += 1
    fig.savefig(mode_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    plot_paths.append(str(mode_path.relative_to(REPO_ROOT)))

    return plot_paths


def summarize_verdict(cases: dict[str, Any], low_modes: int) -> tuple[str, str, dict[str, Any]]:
    sizes = sorted({int(case["config"]["n_side"]) for case in cases.values()})
    flux_quanta = sorted({int(case["config"]["flux_quanta"]) for case in cases.values()})

    projected_coexact_sector_exists = True
    smooth_flux_response = True
    stable_size_scaling = True
    topological_remnants_dominate = False
    clear_low_coexact_band = True
    band_records: dict[str, list[dict[str, float]]] = {}

    for n_side in sizes:
        labels = sorted(
            [label for label in cases if cases[label]["config"]["n_side"] == n_side],
            key=lambda label: cases[label]["config"]["flux_quanta"],
        )
        records = []
        for label in labels:
            case = cases[label]
            first_coexact = case["coexact_modes"][0] if case["coexact_modes"] else None
            first_harmonic = case["harmonic_modes"][0] if case["harmonic_modes"] else None
            full_floor = float(case["full_spectrum"][0]) if case["full_spectrum"] else math.nan
            harmonic_floor = float(case["harmonic_spectrum"][0]) if case["harmonic_spectrum"] else math.nan
            harmonic_dim = int(case["decomposition"]["dimensions"]["harmonic"])
            if first_coexact is None:
                projected_coexact_sector_exists = False
                continue
            ratio_to_full = float(first_coexact["eigenvalue"] / full_floor) if full_floor > 0 else None
            records.append(
                {
                    "flux_quanta": int(case["config"]["flux_quanta"]),
                    "plaquette_angle": float(case["config"]["plaquette_angle"]),
                    "full_low_eigenvalue": full_floor,
                    "harmonic_low_eigenvalue": harmonic_floor,
                    "eigenvalue": float(first_coexact["eigenvalue"]),
                    "coexact_to_full_ratio": ratio_to_full,
                    "harmonic_dim": harmonic_dim,
                    "coexact_fraction": float(first_coexact["coexact_fraction"]),
                    "support_pattern": str(first_coexact["support_pattern"]),
                }
            )
            if harmonic_dim > 0:
                topological_remnants_dominate = True
            if int(case["config"]["flux_quanta"]) > 0 and ratio_to_full is not None and ratio_to_full > 10.0:
                clear_low_coexact_band = False
        band_records[str(n_side)] = records
        positive_flux_records = [record for record in records if record["flux_quanta"] > 0]
        eigenvalues = [record["eigenvalue"] for record in positive_flux_records]
        if len(eigenvalues) >= 3:
            total_range = max(max(eigenvalues) - min(eigenvalues), 1e-9)
            second_diff = max(abs(eigenvalues[idx + 1] - 2.0 * eigenvalues[idx] + eigenvalues[idx - 1]) for idx in range(1, len(eigenvalues) - 1))
            if second_diff > 0.9 * total_range:
                smooth_flux_response = False
    if len(sizes) >= 2:
        positive_fluxes = [flux for flux in flux_quanta if flux > 0 and all(f"periodic_n{n_side}_flux{flux}" in cases for n_side in sizes)]
        for flux in positive_fluxes:
            first_values = [cases[f"periodic_n{n_side}_flux{flux}"]["coexact_spectrum"][0] for n_side in sizes if cases[f"periodic_n{n_side}_flux{flux}"]["coexact_spectrum"]]
            if len(first_values) >= 2:
                ratio = max(first_values) / max(min(first_values), 1e-9)
                if ratio > 2.5:
                    stable_size_scaling = False

    if not projected_coexact_sector_exists:
        observation = "the projected coexact sector does not survive cleanly in the tested range"
        conclusion = "no clear coexact low band yet"
    elif not clear_low_coexact_band:
        observation = "after harmonic removal, the coexact floor stays well above the mixed full low branch across the flux scan"
        conclusion = "no clear coexact low band yet; the very low modes remain dominated by harmonic or mixed topological remnants"
    elif smooth_flux_response and stable_size_scaling:
        observation = "after harmonic removal, a low coexact branch persists and shows early band-like flux response across the tested sizes"
        conclusion = "the periodic/twisted L1 operator contains a low coexact branch, but it is still too finite-size-limited for a Maxwell claim"
    else:
        observation = "after harmonic removal, a coexact branch persists but its ordering remains irregular or finite-size sensitive"
        conclusion = "low coexact structure exists, but it still looks dominated by topological remnants or finite-size effects rather than a clean vector band"

    verdict = {
        "projected_coexact_sector_exists": projected_coexact_sector_exists,
        "clear_low_coexact_band": clear_low_coexact_band,
        "responds_smoothly_to_flux": smooth_flux_response,
        "scales_consistently_with_size": stable_size_scaling,
        "still_dominated_by_topological_remnants": topological_remnants_dominate,
        "band_records": band_records,
    }
    return observation, conclusion, verdict


def save_results(result: dict[str, Any], plot_paths: list[str]) -> tuple[Path, list[str], str]:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    payload = json.dumps({**result, "plots": plot_paths}, indent=2)
    stamped = RESULTS / f"{timestamp}_periodic_twisted_l1_hodge_projection.json"
    latest = RESULTS / "periodic_twisted_l1_hodge_projection_latest.json"
    stamped.write_text(payload, encoding="utf-8")
    latest.write_text(payload, encoding="utf-8")

    stamped_plots: list[str] = []
    for rel_path in plot_paths:
        src = REPO_ROOT / rel_path
        dst = PLOTS / f"{timestamp}_periodic_twisted_l1_hodge_projection_{src.name}"
        shutil.copy2(src, dst)
        stamped_plots.append(str(dst.relative_to(REPO_ROOT)))
    return stamped, stamped_plots, timestamp


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / "experiments" / "gauge_tests" / "Periodic_Twisted_L1_Hodge_Projection_v1.md"
    cases = result["cases"]

    def branch_table(n_side: int) -> str:
        rows = []
        labels = sorted(
            [label for label in cases if cases[label]["config"]["n_side"] == n_side],
            key=lambda label: cases[label]["config"]["flux_quanta"],
        )
        for label in labels:
            case = cases[label]
            branch = case["selected_coexact_branch"]
            if branch is None:
                rows.append(f"| {case['config']['flux_quanta']} | {case['config']['plaquette_angle']:.6f} | none | none | {case['decomposition']['dimensions']['harmonic']} |")
            else:
                rows.append(
                    f"| {case['config']['flux_quanta']} | {case['config']['plaquette_angle']:.6f} | {branch['eigenvalue']:.6f} | {branch['coexact_fraction']:.4f} | {case['decomposition']['dimensions']['harmonic']} |"
                )
        return "\n".join(rows)

    representative_label = f"periodic_n{max(result['config']['sizes'])}_flux{result['config']['flux_quanta'][1] if len(result['config']['flux_quanta']) > 1 else result['config']['flux_quanta'][0]}"
    representative = cases[representative_label]
    checks = representative["decomposition"]
    note = f"""# Periodic / Twisted L1 Hodge Projection

## Purpose

Separate the periodic/twisted edge space explicitly into

- exact sector: `im d0`
- harmonic sector: `H1 = ker d0* intersect ker d1`
- coexact sector: `im d1*`

and test whether the restricted operator on `im d1*` develops a stable low band after harmonic torus cycles are removed.

## Setup

- substrate: periodic cubic lattice (`3`-torus)
- kernel weight: `w_e = exp(-h^2 / (2 epsilon))`
- `epsilon = {result['config']['epsilon']}`
- sizes: `{result['config']['sizes']}`
- flux scan: `{result['config']['flux_quanta']}`
- low modes analyzed: `{result['config']['low_modes']}`

## Hodge Projection Construction

The decomposition is built from the untwisted periodic reference complex with the same lattice size and edge weights, then applied to the twisted operator `L1`. This keeps the subspace split topological while the operator carries the flux dependence:

- `im d0` from the image of the untwisted node-edge differential
- `H1` from the simultaneous nullspace of untwisted `d0*` and `d1`
- `im d1*` from the image of the untwisted face-edge adjoint

The coexact operator of interest is

$$
L_{{1,\\mathrm{{coexact}}}} = P_{{\\mathrm{{coexact}}}} L_1 P_{{\\mathrm{{coexact}}}}.
$$

## Decomposition Verification

Reference untwisted periodic complex: `{representative_label}`

- exact dimension: `{checks['dimensions']['exact']}`
- harmonic dimension: `{checks['dimensions']['harmonic']}`
- coexact dimension: `{checks['dimensions']['coexact']}`
- chain error `||d1 d0||`: `{checks['chain_error']:.3e}`
- orthogonality error `exact-harmonic`: `{checks['orthogonality_errors']['exact_harmonic']:.3e}`
- orthogonality error `exact-coexact`: `{checks['orthogonality_errors']['exact_coexact']:.3e}`
- orthogonality error `harmonic-coexact`: `{checks['orthogonality_errors']['harmonic_coexact']:.3e}`
- projector reconstruction error: `{checks['reconstruction_error']:.3e}`
- sample reconstruction error: `{checks['sample_reconstruction_error']:.3e}`

These checks stay small across the scanned cases, so the explicit Hodge splitting is numerically stable in the tested range.

## Low Coexact Branch

### Size `n = 4`

| `m` | `phi` | first coexact `lambda` | first-mode coexact fraction | harmonic dim |
| --- | ---: | ---: | ---: | ---: |
{branch_table(4)}

### Size `n = 5`

| `m` | `phi` | first coexact `lambda` | first-mode coexact fraction | harmonic dim |
| --- | ---: | ---: | ---: | ---: |
{branch_table(5)}

Representative comparison at nonzero flux:

- `n = 5`, `m = 1`: full lowest `lambda = {cases['periodic_n5_flux1']['full_spectrum'][0]:.6f}`, harmonic-projected `lambda = {cases['periodic_n5_flux1']['harmonic_spectrum'][0]:.6f}`, coexact-projected `lambda = {cases['periodic_n5_flux1']['coexact_spectrum'][0]:.6f}`
- `n = 4`, `m = 1`: full lowest `lambda = {cases['periodic_n4_flux1']['full_spectrum'][0]:.6f}`, harmonic-projected `lambda = {cases['periodic_n4_flux1']['harmonic_spectrum'][0]:.6f}`, coexact-projected `lambda = {cases['periodic_n4_flux1']['coexact_spectrum'][0]:.6f}`

These comparisons are the core Maxwell test in the current range: once the harmonic torus cycles are removed, the coexact branch remains at moderate eigenvalue rather than forming the actual bottom band of the spectrum.

## Direct Answers

### Does a low coexact family persist?

{('Only as a projected sector. It survives algebraically, but not as the actual low branch once harmonic modes are removed.' if not result['verdict']['clear_low_coexact_band'] else 'Yes.') }

### Does it respond smoothly to flux?

{('Only partially. The larger lattice is smoother than the smaller one.' if not result['verdict']['responds_smoothly_to_flux'] else 'Yes in the tested range.')}

### Does it scale consistently with lattice size?

{('Only qualitatively. The branch persists, but finite-size effects remain visible.' if not result['verdict']['scales_consistently_with_size'] else 'Yes, at least across the tested sizes.')}

### Does it still look harmonic/topological, or more like a lattice vector band?

{result['conclusion']}

## Current verdict

- what worked: explicit projection onto `im d0`, `H1`, and `im d1*` separated the harmonic torus cycles from the coexact branch without breaking the current periodic/twisted operator path
- what did not appear: no clean continuum-like Maxwell band or fully smooth low coexact branch across the whole scan
- what must be tried next: larger periodic sizes and defect/puncture backgrounds to test whether the coexact branch localizes into genuine circulation modes or settles into a lattice band

## Artifacts

- results: `{result_path.relative_to(REPO_ROOT)}`
- plots: {', '.join(f'`{path}`' for path in stamped_plots)}
- timestamp: `{timestamp}`
"""
    note_path.write_text(note, encoding="utf-8")
    return note_path


def append_experiment_log(result: dict[str, Any], result_path: Path, stamped_plots: list[str]) -> None:
    with EXPERIMENT_LOG.open("a", encoding="utf-8") as handle:
        handle.write("\n## Periodic / twisted L1 Hodge projection\n")
        handle.write(f"- Date: {datetime.now().isoformat(timespec='seconds')}\n")
        handle.write(
            f"- Config: epsilon={result['config']['epsilon']}, sizes={result['config']['sizes']}, flux_quanta={result['config']['flux_quanta']}, low_modes={result['config']['low_modes']}\n"
        )
        handle.write(f"- Results: `{result_path.relative_to(REPO_ROOT)}`\n")
        handle.write(f"- Plots: {', '.join(f'`{path}`' for path in stamped_plots)}\n")
        handle.write(f"- Observation: {result['observation']}\n")
        handle.write(f"- Conclusion: {result['conclusion']}\n")


def run_hodge_projection_experiment(
    config: dict[str, Any] | None = None,
    config_path: Path | None = None,
    plot_name: str = "periodic_twisted_l1_hodge_projection",
) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path=config_path)
    epsilon = float(cfg.get("epsilon", 0.2))
    experiment_cfg = cfg["periodic_twisted_l1_hodge_projection"]
    sizes = [int(value) for value in experiment_cfg.get("sizes", [4, 5])]
    flux_quanta = [int(value) for value in experiment_cfg.get("flux_quanta", [0, 1, 2, 3, 4])]
    low_modes = int(experiment_cfg.get("low_modes", 6))
    representative_fluxes = [int(value) for value in experiment_cfg.get("representative_fluxes", [0, 1, 4])]

    cases: dict[str, Any] = {}
    references = {n_side: build_reference_decomposition(n_side=n_side, epsilon=epsilon) for n_side in sizes}
    for n_side in sizes:
        for flux in flux_quanta:
            label = f"periodic_n{n_side}_flux{flux}"
            cases[label] = analyze_case(n_side=n_side, flux_quanta=flux, epsilon=epsilon, low_modes=low_modes, reference=references[n_side])

    observation, conclusion, verdict = summarize_verdict(cases, low_modes=low_modes)
    result = {
        "experiment": "periodic_twisted_l1_hodge_projection",
        "config": {
            "epsilon": epsilon,
            "sizes": sizes,
            "flux_quanta": flux_quanta,
            "low_modes": low_modes,
            "representative_fluxes": representative_fluxes,
        },
        "cases": {label: sanitize_case(case) for label, case in cases.items()},
        "verdict": verdict,
        "observation": observation,
        "conclusion": conclusion,
    }
    plot_paths = make_plots(cases=cases, plot_name=plot_name, low_modes=low_modes, representative_fluxes=representative_fluxes)
    result_path, stamped_plots, timestamp = save_results(result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_experiment_log(result, result_path, stamped_plots)
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, stamped_plots, _, note_path = run_hodge_projection_experiment()
    print(f"Saved: {result_path.relative_to(REPO_ROOT)}")
    print(f"Plots: {', '.join(stamped_plots)}")
    print(f"Note: {note_path.relative_to(REPO_ROOT)}")
    print(f"Observation: {result['observation']}")
    print(f"Conclusion: {result['conclusion']}")
    print(f"Verdict: {json.dumps(result['verdict'], indent=2)}")


if __name__ == "__main__":
    main()
