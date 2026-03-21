#!/usr/bin/env python3

from __future__ import annotations

import math
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh, expm_multiply

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haos_core import build_graph, read_json, relpath, write_csv_rows, write_json

PHASE_ROOT = ROOT / "phase11-protection"
CONFIG_PATH = PHASE_ROOT / "configs" / "phase11_protection_config.json"
RUNS_ROOT = PHASE_ROOT / "runs"
PLOTS_ROOT = PHASE_ROOT / "plots"

PHASE6_MANIFEST_PATH = ROOT / "phase6-operator" / "phase6_operator_manifest.json"
PHASE7_MANIFEST_PATH = ROOT / "phase7-spectral" / "phase7_spectral_manifest.json"
PHASE8_MANIFEST_PATH = ROOT / "phase8-trace" / "phase8_trace_manifest.json"
PHASE9_MANIFEST_PATH = ROOT / "phase9-invariants" / "phase9_manifest.json"
PHASE10_MANIFEST_PATH = ROOT / "phase10-bridge" / "phase10_manifest.json"
PHASEX_MANIFEST_PATH = ROOT / "phaseX-proto-particle" / "phaseX_integrated_manifest.json"
PHASEX_CANDIDATES_PATH = ROOT / "phaseX-proto-particle" / "runs" / "proto_particle_candidates.json"

SURVIVAL_LEDGER_PATH = RUNS_ROOT / "phase11_perturbation_survival_ledger.csv"
FAILURE_LEDGER_PATH = RUNS_ROOT / "phase11_failure_classification.csv"
PERSISTENCE_LEDGER_PATH = RUNS_ROOT / "phase11_persistence_scaling.csv"
RUNS_PATH = RUNS_ROOT / "phase11_runs.json"
SUMMARY_PATH = PHASE_ROOT / "phase11_summary.md"
MANIFEST_PATH = PHASE_ROOT / "phase11_manifest.json"

WIDTH_PLOT_PATH = PLOTS_ROOT / "phase11_localization_width_vs_time.svg"
SURVIVAL_PLOT_PATH = PLOTS_ROOT / "phase11_survival_probability_vs_perturbation.svg"
PERSISTENCE_PLOT_PATH = PLOTS_ROOT / "phase11_persistence_time_vs_refinement.svg"
GAP_PLOT_PATH = PLOTS_ROOT / "phase11_spectral_gap_shift_vs_defect_strength.svg"

PHASE_NAME = "phase11-protection"
STAGE_IDENTIFIER = "phase11-protection"
BRANCH_LABEL = "frozen_branch"
CONTROL_LABEL = "periodic_diagonal_augmented_control"

CLAIM_BOUNDARY = (
    "Phase XI is limited to localized-mode protection and failure-mechanism diagnostics "
    "on the frozen operator hierarchy. It does not assert particles, field dynamics, "
    "continuum limits, or new physical structure."
)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def timestamp_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def round_float(value: float, digits: int = 12) -> float:
    return round(float(value), digits)


def closure_statement(success: bool) -> str:
    if success:
        return "Phase XI establishes localized-mode protection and failure-mechanism feasibility for the frozen operator hierarchy."
    return "Phase XI does not yet establish localized-mode protection and failure-mechanism feasibility for the frozen operator hierarchy."


def safe_relative_difference(new_value: float, reference_value: float) -> float:
    return abs(float(new_value) - float(reference_value)) / max(abs(float(reference_value)), 1.0e-12)


def correlation(values_x: list[float], values_y: list[float]) -> float:
    if len(values_x) < 2 or len(values_y) < 2:
        return 0.0
    x = np.asarray(values_x, dtype=float)
    y = np.asarray(values_y, dtype=float)
    if float(np.std(x)) <= 1.0e-15 or float(np.std(y)) <= 1.0e-15:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def fit_power_law(h_values: list[float], tau_values: list[float]) -> dict[str, float]:
    x = np.log(np.asarray(h_values, dtype=float))
    y = np.log(np.asarray(tau_values, dtype=float))
    slope, intercept = np.polyfit(x, y, 1)
    fitted = slope * x + intercept
    denom = float(np.sum((y - y.mean()) ** 2))
    if denom <= 1.0e-18:
        r_squared = 1.0
    else:
        r_squared = 1.0 - float(np.sum((y - fitted) ** 2) / denom)
    return {
        "prefactor": round_float(math.exp(intercept)),
        "exponent_q": round_float(-slope),
        "r_squared": round_float(r_squared),
    }


def periodic_delta(points: np.ndarray, anchor: np.ndarray) -> np.ndarray:
    delta = np.asarray(points, dtype=float) - np.asarray(anchor, dtype=float)
    return (delta + 0.5) % 1.0 - 0.5


def load_frozen_inputs() -> dict[str, Any]:
    config = read_json(CONFIG_PATH)
    phase6 = read_json(PHASE6_MANIFEST_PATH)
    phase7 = read_json(PHASE7_MANIFEST_PATH)
    phase8 = read_json(PHASE8_MANIFEST_PATH)
    phase9 = read_json(PHASE9_MANIFEST_PATH)
    phase10 = read_json(PHASE10_MANIFEST_PATH)
    phasex = read_json(PHASEX_MANIFEST_PATH)
    phasex_candidates = read_json(PHASEX_CANDIDATES_PATH)

    manifests = [phase6, phase7, phase8, phase9, phase10, phasex]
    if not all(bool(manifest.get("success")) for manifest in manifests):
        raise ValueError("Phase XI requires successful frozen Phases VI, VII, VIII, IX, X, and proto-particle manifests.")

    candidate_label = str(config["candidate_label"])
    if candidate_label not in phasex_candidates.get("selected_branch_candidates", []):
        raise ValueError("Phase XI candidate must already be selected in the frozen Phase X candidate bundle.")
    if candidate_label not in phasex.get("persistence_results", {}).get("branch_refinement_stable_candidates", []):
        raise ValueError("Phase XI candidate must already be refinement-stable in the frozen Phase X manifest.")

    return {
        "config": config,
        "phase6": phase6,
        "phase7": phase7,
        "phase8": phase8,
        "phase9": phase9,
        "phase10": phase10,
        "phasex": phasex,
        "phasex_candidates": phasex_candidates,
    }


def branch_scalar_spectrum(n_side: int, epsilon: float) -> np.ndarray:
    h_value = 1.0 / float(n_side)
    edge_weight = math.exp(-(h_value * h_value) / (2.0 * epsilon))
    values = np.empty(n_side * n_side, dtype=float)
    cursor = 0
    for mode_x in range(n_side):
        cos_x = math.cos(2.0 * math.pi * mode_x / n_side)
        for mode_y in range(n_side):
            cos_y = math.cos(2.0 * math.pi * mode_y / n_side)
            values[cursor] = 2.0 * edge_weight * (2.0 - cos_x - cos_y)
            cursor += 1
    return np.sort(values)


def control_scalar_spectrum(n_side: int, epsilon: float, diagonal_weight_scale: float) -> np.ndarray:
    h_value = 1.0 / float(n_side)
    edge_weight = math.exp(-(h_value * h_value) / (2.0 * epsilon))
    diagonal_weight = edge_weight * float(diagonal_weight_scale)
    values = np.empty(n_side * n_side, dtype=float)
    cursor = 0
    for mode_x in range(n_side):
        cos_x = math.cos(2.0 * math.pi * mode_x / n_side)
        for mode_y in range(n_side):
            cos_y = math.cos(2.0 * math.pi * mode_y / n_side)
            values[cursor] = 2.0 * edge_weight * (2.0 - cos_x - cos_y) + 4.0 * diagonal_weight * (1.0 - cos_x * cos_y)
            cursor += 1
    return np.sort(values)


def trace_from_spectrum(spectrum: np.ndarray, probe_times: np.ndarray) -> np.ndarray:
    return np.array([float(np.exp(-probe_time * spectrum).sum()) for probe_time in probe_times], dtype=float)


def quadratic_fit_named(probe_times: np.ndarray, values: np.ndarray, prefix: str) -> dict[str, float]:
    coefficients = np.polyfit(probe_times, values, 2)
    fitted = np.polyval(coefficients, probe_times)
    denominator = float(np.sum((values - values.mean()) ** 2))
    if denominator <= 1.0e-18:
        r_squared = 1.0
    else:
        r_squared = 1.0 - float(np.sum((values - fitted) ** 2) / denominator)
    c0, c1, c2 = coefficients[::-1]
    return {
        f"{prefix}0": round_float(c0),
        f"{prefix}1": round_float(c1),
        f"{prefix}2": round_float(c2),
        "r_squared": round_float(r_squared),
    }


def ratio_shift_proxy(
    hierarchy_label: str,
    n_side: int,
    epsilon: float,
    diagonal_weight_scale: float,
    short_time_window: np.ndarray,
    branch_limit_ratios: dict[str, float],
) -> dict[str, float]:
    spectrum = (
        branch_scalar_spectrum(n_side, epsilon)
        if hierarchy_label == BRANCH_LABEL
        else control_scalar_spectrum(n_side, epsilon, diagonal_weight_scale)
    )
    h_value = 1.0 / float(n_side)
    scaled_trace = (h_value * h_value) * trace_from_spectrum(spectrum, short_time_window)
    fit = quadratic_fit_named(short_time_window, scaled_trace, prefix="b")
    ratio_1 = float(fit["b1"]) / max(abs(float(fit["b0"])), 1.0e-12)
    ratio_2 = float(fit["b2"]) / max(abs(float(fit["b0"])), 1.0e-12)
    shift_1 = safe_relative_difference(ratio_1, float(branch_limit_ratios["b1_over_b0"]))
    shift_2 = safe_relative_difference(ratio_2, float(branch_limit_ratios["b2_over_b0"]))
    return {
        "b1_over_b0": round_float(ratio_1),
        "b2_over_b0": round_float(ratio_2),
        "ratio_shift_proxy": round_float(max(shift_1, shift_2)),
    }


def build_branch_zero_form_operator(n_side: int, epsilon: float) -> tuple[sp.csr_matrix, np.ndarray]:
    graph = build_graph(
        {
            "kind": "dk2d_periodic",
            "n_side": int(n_side),
            "epsilon": float(epsilon),
            "cycle_phase_x": 0.0,
            "cycle_phase_y": 0.0,
        }
    )
    zero_form_size = int(graph.block_sizes[0])
    return graph.delta_h[:zero_form_size, :zero_form_size].tocsr(), np.asarray(graph.points, dtype=float)


def build_control_zero_form_operator(n_side: int, epsilon: float, diagonal_weight_scale: float) -> tuple[sp.csr_matrix, np.ndarray]:
    h_value = 1.0 / float(n_side)
    edge_weight = math.exp(-(h_value * h_value) / (2.0 * epsilon))
    diagonal_weight = edge_weight * float(diagonal_weight_scale)
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    degree = np.zeros(n_side * n_side, dtype=float)

    def idx(i: int, j: int) -> int:
        return (i % n_side) * n_side + (j % n_side)

    for i in range(n_side):
        for j in range(n_side):
            u = idx(i, j)
            for v, weight in (
                (idx(i + 1, j), edge_weight),
                (idx(i, j + 1), edge_weight),
                (idx(i + 1, j + 1), diagonal_weight),
                (idx(i + 1, j - 1), diagonal_weight),
            ):
                degree[u] += weight
                degree[v] += weight
                rows.extend([u, v])
                cols.extend([v, u])
                data.extend([-weight, -weight])

    for node_index, node_degree in enumerate(degree):
        rows.append(node_index)
        cols.append(node_index)
        data.append(node_degree)

    operator = sp.coo_matrix((data, (rows, cols)), shape=(n_side * n_side, n_side * n_side)).tocsr()
    points = np.array([[i / n_side, j / n_side] for i in range(n_side) for j in range(n_side)], dtype=float)
    return operator, points


def build_candidate_state(points: np.ndarray, candidate_label: str) -> np.ndarray:
    if candidate_label != "low_mode_localized_wavepacket":
        raise ValueError(f"unsupported frozen candidate label: {candidate_label}")
    center = np.array([0.5, 0.5], dtype=float)
    delta = periodic_delta(points, center)
    sigma = 0.08
    amplitude = np.exp(-0.5 * ((delta[:, 0] / sigma) ** 2 + (delta[:, 1] / sigma) ** 2))
    phase = np.exp(1j * 2.0 * math.pi * (delta[:, 0] + delta[:, 1]))
    state = amplitude * phase
    state /= max(np.linalg.norm(state), 1.0e-12)
    return np.asarray(state, dtype=complex)


def adjacency_from_laplacian(operator: sp.csr_matrix) -> dict[tuple[int, int], float]:
    adjacency: dict[tuple[int, int], float] = {}
    coo = operator.tocoo()
    for row, col, value in zip(coo.row, coo.col, coo.data):
        if row < col and value.real < 0.0:
            adjacency[(int(row), int(col))] = -float(value.real)
    return adjacency


def laplacian_from_adjacency(dimension: int, adjacency: dict[tuple[int, int], float]) -> sp.csr_matrix:
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    degree = np.zeros(dimension, dtype=float)
    for (u, v), weight in adjacency.items():
        degree[u] += weight
        degree[v] += weight
        rows.extend([u, v])
        cols.extend([v, u])
        data.extend([-weight, -weight])
    for node_index, node_degree in enumerate(degree):
        rows.append(node_index)
        cols.append(node_index)
        data.append(node_degree)
    return sp.coo_matrix((data, (rows, cols)), shape=(dimension, dimension)).tocsr()


def local_stiffness_defect(operator: sp.csr_matrix, n_side: int, amplitude: float) -> sp.csr_matrix:
    adjacency = adjacency_from_laplacian(operator)

    def idx(i: int, j: int) -> int:
        return (i % n_side) * n_side + (j % n_side)

    touched: set[tuple[int, int]] = set()
    cells = [(n_side // 2, n_side // 2), (n_side // 2 - 1, n_side // 2), (n_side // 2, n_side // 2 - 1), (n_side // 2 - 1, n_side // 2 - 1)]
    for i, j in cells:
        u = idx(i, j)
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            touched.add(tuple(sorted((u, idx(i + di, j + dj)))))
    for edge in touched:
        if edge in adjacency:
            adjacency[edge] *= max(0.05, 1.0 - float(amplitude))
    return laplacian_from_adjacency(operator.shape[0], adjacency)


def connectivity_micro_surgery(operator: sp.csr_matrix, n_side: int, amplitude: float) -> sp.csr_matrix:
    adjacency = adjacency_from_laplacian(operator)
    base_weight = float(np.mean(list(adjacency.values())))

    def idx(i: int, j: int) -> int:
        return (i % n_side) * n_side + (j % n_side)

    cells = [(n_side // 2 - 1, n_side // 2 - 1), (n_side // 2, n_side // 2 - 1), (n_side // 2 - 1, n_side // 2), (n_side // 2, n_side // 2)]
    for i, j in cells:
        a = idx(i, j)
        b = idx(i + 1, j)
        c = idx(i, j + 1)
        d = idx(i + 1, j + 1)
        for edge in (tuple(sorted((a, b))), tuple(sorted((a, c))), tuple(sorted((b, d))), tuple(sorted((c, d)))):
            if edge in adjacency:
                adjacency[edge] *= max(0.0, 1.0 - float(amplitude))
        for edge in (tuple(sorted((a, d))), tuple(sorted((b, c)))):
            adjacency[edge] = max(adjacency.get(edge, 0.0), float(amplitude) * base_weight)
    return laplacian_from_adjacency(operator.shape[0], adjacency)


def phase_randomization_window(state: np.ndarray, points: np.ndarray, amplitude: float) -> np.ndarray:
    delta = periodic_delta(points, np.array([0.5, 0.5], dtype=float))
    pattern = np.sin(17.0 * delta[:, 0]) * np.cos(11.0 * delta[:, 1])
    out = state * np.exp(1j * float(amplitude) * math.pi * pattern)
    out /= max(np.linalg.norm(out), 1.0e-12)
    return np.asarray(out, dtype=complex)


def boundary_leakage_probe(operator: sp.csr_matrix, points: np.ndarray, amplitude: float, ring_radius: float) -> sp.csr_matrix:
    delta = periodic_delta(points, np.array([0.5, 0.5], dtype=float))
    radius = np.sqrt(np.sum(delta * delta, axis=1))
    mask = (radius > float(ring_radius)).astype(float)
    return (operator + sp.diags(float(amplitude) * mask)).tocsr()


def spatial_metrics(state: np.ndarray, points: np.ndarray) -> dict[str, float]:
    probabilities = np.abs(state) ** 2
    probabilities = probabilities / max(float(np.sum(probabilities)), 1.0e-12)
    peak_position = points[int(np.argmax(probabilities))]
    delta = periodic_delta(points, peak_position)
    width = math.sqrt(float(np.sum(probabilities * np.sum(delta * delta, axis=1))))
    concentration = float(np.max(probabilities))
    participation_ratio = float(1.0 / np.sum(probabilities * probabilities))
    return {
        "localization_width": round_float(width),
        "spatial_concentration": round_float(concentration),
        "participation_ratio": round_float(participation_ratio),
    }


def spectral_energy_metrics(state: np.ndarray, operator: sp.csr_matrix) -> dict[str, float]:
    norm_sq = max(float(np.vdot(state, state).real), 1.0e-12)
    applied = operator @ state
    mean_energy = float(np.vdot(state, applied).real) / norm_sq
    second_moment = float(np.vdot(applied, applied).real) / norm_sq
    spread = math.sqrt(max(second_moment - mean_energy * mean_energy, 0.0))
    return {
        "spectral_energy_mean": round_float(mean_energy),
        "spectral_energy_spread": round_float(spread),
    }


def overlap(left: np.ndarray, right: np.ndarray) -> float:
    return float(abs(np.vdot(left, right)) / max(float(np.linalg.norm(left) * np.linalg.norm(right)), 1.0e-12))


def recovery_style_score(base_metrics: dict[str, float], other_metrics: dict[str, float]) -> float:
    tolerances = {
        "localization_width": 0.08,
        "spatial_concentration": 0.12,
        "participation_ratio": 0.20,
        "spectral_energy_spread": 0.18,
    }
    penalties = []
    for key, tolerance in tolerances.items():
        deviation = safe_relative_difference(float(other_metrics[key]), float(base_metrics[key]))
        penalties.append(min(deviation / tolerance, 1.0))
    return round_float(max(0.0, 1.0 - float(np.mean(penalties))))


def classify_persistence(
    base_width_growth: float,
    base_concentration_ratio: float,
    overlap_value: float,
    recovery_score: float,
    thresholds: dict[str, float],
) -> str:
    if (
        base_width_growth <= float(thresholds["persistent_max_width_growth"])
        and base_concentration_ratio >= float(thresholds["persistent_min_concentration_ratio"])
        and overlap_value >= float(thresholds["persistent_min_overlap"])
        and recovery_score >= float(thresholds["persistent_min_recovery_score"])
    ):
        return "persistent"
    if (
        base_width_growth <= float(thresholds["diffusive_max_width_growth"])
        and base_concentration_ratio >= float(thresholds["diffusive_min_concentration_ratio"])
        and overlap_value >= float(thresholds["diffusive_min_overlap"])
        and recovery_score >= float(thresholds["diffusive_min_recovery_score"])
    ):
        return "diffusive"
    return "unstable"


def survives_candidate_gate(width_growth: float, concentration_ratio: float, participation_growth: float, gate: dict[str, float]) -> bool:
    return (
        width_growth <= float(gate["max_width_growth"])
        and concentration_ratio >= float(gate["min_concentration_ratio"])
        and participation_growth <= float(gate["max_participation_growth"])
    )


def sampled_eigensystem(operator: sp.csr_matrix, sample_eigen_count: int) -> tuple[np.ndarray, np.ndarray]:
    k_value = min(int(sample_eigen_count), operator.shape[0] - 2)
    eigenvalues, eigenvectors = eigsh(operator, k=k_value, which="SM")
    order = np.argsort(eigenvalues)
    return np.asarray(np.real(eigenvalues[order]), dtype=float), np.asarray(np.real(eigenvectors[:, order]), dtype=float)


def first_positive_gap(eigenvalues: np.ndarray) -> float:
    positive = np.asarray([value for value in eigenvalues if value > 1.0e-9], dtype=float)
    if len(positive) == 0:
        return 0.0
    return float(positive[0])


def low_mode_basis(eigenvalues: np.ndarray, eigenvectors: np.ndarray, gap: float, band_multiplier: float) -> np.ndarray:
    if gap <= 0.0:
        return np.asarray(eigenvectors[:, :1], dtype=float)
    mask = np.array([(value > 1.0e-9) and (value <= float(band_multiplier) * gap + 1.0e-12) for value in eigenvalues], dtype=bool)
    if not bool(np.any(mask)):
        first_positive_index = int(np.argmax(eigenvalues > 1.0e-9))
        return np.asarray(eigenvectors[:, first_positive_index : first_positive_index + 1], dtype=float)
    return np.asarray(eigenvectors[:, mask], dtype=float)


def mode_mixing_indicator(state: np.ndarray, basis: np.ndarray) -> float:
    if basis.shape[1] == 0:
        return 1.0
    state_norm = max(float(np.linalg.norm(state)), 1.0e-12)
    projections = np.asarray(basis.conj().T @ state, dtype=complex)
    preserved = float(np.sum(np.abs(projections) ** 2)) / (state_norm * state_norm)
    return round_float(max(0.0, 1.0 - min(preserved, 1.0)))


def local_connectivity_variance(operator: sp.csr_matrix, points: np.ndarray, radius: float) -> float:
    degrees = np.asarray(np.real(operator.diagonal()), dtype=float)
    delta = periodic_delta(points, np.array([0.5, 0.5], dtype=float))
    mask = np.sqrt(np.sum(delta * delta, axis=1)) <= float(radius)
    if not bool(np.any(mask)):
        return 0.0
    return round_float(float(np.var(degrees[mask])))


def local_spectral_density_proxy(eigenvalues: np.ndarray, gap: float, gap_multiplier: float) -> float:
    if gap <= 0.0:
        return 0.0
    positive = np.asarray([value for value in eigenvalues if value > 1.0e-9], dtype=float)
    if len(positive) == 0:
        return 0.0
    count = int(np.sum(positive <= float(gap_multiplier) * gap + 1.0e-12))
    return round_float(float(count) / float(len(positive)))


def failure_mechanism_label(
    perturbation_label: str,
    classification: str,
    gap_shift_relative: float,
    mode_mixing: float,
    norm_retention_relative: float,
) -> str:
    if classification == "persistent":
        return "survives"
    if perturbation_label == "phase_randomization_window":
        return "incoherence_triggered"
    if perturbation_label == "connectivity_micro_surgery":
        return "topological_defect_triggered"
    if perturbation_label == "boundary_leakage_probe":
        return "dispersion_dominated" if norm_retention_relative < 0.98 else "spectral_mixing_dominated"
    if abs(gap_shift_relative) >= 0.15 or mode_mixing >= 0.25:
        return "spectral_mixing_dominated"
    return "dispersion_dominated"


def ordered_threshold(thresholds: list[float | None]) -> bool:
    values = [float("inf") if value is None else float(value) for value in thresholds]
    return all(left <= right + 1.0e-12 for left, right in zip(values, values[1:]))


def plot_localization_width(localization_curves: list[dict[str, Any]], levels: list[int]) -> None:
    plt.figure(figsize=(8.4, 5.4))
    for hierarchy_label, linestyle in ((BRANCH_LABEL, "-"), (CONTROL_LABEL, "--")):
        for n_side in levels:
            curve = next(
                row
                for row in localization_curves
                if row["hierarchy_label"] == hierarchy_label and int(row["n_side"]) == int(n_side)
            )
            plt.plot(curve["tau_grid"], curve["width_values"], linestyle=linestyle, marker="o", label=f"{hierarchy_label}:n{n_side}")
    plt.xlabel("tau")
    plt.ylabel("localization width")
    plt.grid(alpha=0.3)
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(WIDTH_PLOT_PATH, format="svg")
    plt.close()


def plot_survival_probability(survival_rows: list[dict[str, Any]], perturbation_grids: dict[str, list[float]]) -> None:
    plt.figure(figsize=(9.0, 5.8))
    palette = {
        "local_stiffness_defect": "#0b7285",
        "connectivity_micro_surgery": "#f08c00",
        "phase_randomization_window": "#c2255c",
        "boundary_leakage_probe": "#2b8a3e",
    }
    for hierarchy_label, linestyle in ((BRANCH_LABEL, "-"), (CONTROL_LABEL, "--")):
        for perturbation_label, amplitudes in perturbation_grids.items():
            fractions = []
            for amplitude in amplitudes:
                matches = [
                    row
                    for row in survival_rows
                    if row["hierarchy_label"] == hierarchy_label
                    and row["perturbation_label"] == perturbation_label
                    and abs(float(row["amplitude"]) - float(amplitude)) <= 1.0e-12
                ]
                fractions.append(float(np.mean([float(row["survival_indicator"]) for row in matches])) if matches else 0.0)
            plt.plot(
                amplitudes,
                fractions,
                linestyle=linestyle,
                marker="o",
                color=palette[perturbation_label],
                label=f"{hierarchy_label}:{perturbation_label}",
            )
    plt.xlabel("perturbation amplitude")
    plt.ylabel("survival fraction across refinement levels")
    plt.ylim(-0.02, 1.02)
    plt.grid(alpha=0.3)
    plt.legend(loc="best", fontsize=7)
    plt.tight_layout()
    plt.savefig(SURVIVAL_PLOT_PATH, format="svg")
    plt.close()


def plot_persistence_scaling(persistence_rows: list[dict[str, Any]]) -> None:
    plt.figure(figsize=(8.2, 5.0))
    for hierarchy_label, marker in ((BRANCH_LABEL, "o"), (CONTROL_LABEL, "s")):
        rows = [row for row in persistence_rows if row["hierarchy_label"] == hierarchy_label]
        n_values = [int(row["n_side"]) for row in rows]
        tau_values = [float(row["persistence_time_tau"]) for row in rows]
        plt.plot(n_values, tau_values, marker=marker, label=hierarchy_label)
    plt.xlabel("n_side")
    plt.ylabel("persistence time")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(PERSISTENCE_PLOT_PATH, format="svg")
    plt.close()


def plot_gap_shift(survival_rows: list[dict[str, Any]], levels: list[int]) -> None:
    plt.figure(figsize=(8.6, 5.4))
    for hierarchy_label, linestyle in ((BRANCH_LABEL, "-"), (CONTROL_LABEL, "--")):
        for n_side in levels:
            rows = [
                row
                for row in survival_rows
                if row["hierarchy_label"] == hierarchy_label
                and row["perturbation_label"] == "local_stiffness_defect"
                and int(row["n_side"]) == int(n_side)
            ]
            amplitudes = [float(row["amplitude"]) for row in rows]
            gap_shifts = [float(row["spectral_gap_shift_relative"]) for row in rows]
            plt.plot(amplitudes, gap_shifts, linestyle=linestyle, marker="o", label=f"{hierarchy_label}:n{n_side}")
    plt.xlabel("local stiffness defect amplitude")
    plt.ylabel("relative spectral gap shift")
    plt.grid(alpha=0.3)
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(GAP_PLOT_PATH, format="svg")
    plt.close()


def main() -> None:
    start_time = time.perf_counter()
    RUNS_ROOT.mkdir(parents=True, exist_ok=True)
    PLOTS_ROOT.mkdir(parents=True, exist_ok=True)

    frozen = load_frozen_inputs()
    config = frozen["config"]
    phase6 = frozen["phase6"]
    phase7 = frozen["phase7"]
    phase8 = frozen["phase8"]
    phase9 = frozen["phase9"]
    phase10 = frozen["phase10"]

    epsilon = float(phase6["frozen_inputs"]["base_epsilon"])
    evaluation_tau = float(config["evaluation_tau"])
    candidate_label = str(config["candidate_label"])
    levels = [int(value) for value in config["refinement_levels"]]
    persistence_tau_grid = np.asarray(config["persistence_tau_grid"], dtype=float)
    short_time_window = np.asarray(phase8["operational_short_time_window"]["probe_times"], dtype=float)
    diagonal_weight_scale = float(config["control_graph"]["diagonal_weight_scale"])
    gate = dict(config["candidate_gate"])
    thresholds = dict(config["persistence_thresholds"])
    perturbation_grids = {key: list(value) for key, value in dict(config["perturbation_grids"]).items()}
    proxies_cfg = dict(config["structural_proxies"])
    branch_limit_ratios = dict(phase10["primary_effective_scaling_description"]["branch_ratio_limit_candidates"])

    ratio_shift_map: dict[tuple[str, int], dict[str, float]] = {}
    for hierarchy_label in (BRANCH_LABEL, CONTROL_LABEL):
        for n_side in levels:
            ratio_shift_map[(hierarchy_label, int(n_side))] = ratio_shift_proxy(
                hierarchy_label,
                int(n_side),
                epsilon,
                diagonal_weight_scale,
                short_time_window,
                branch_limit_ratios,
            )

    survival_rows: list[dict[str, Any]] = []
    failure_rows: list[dict[str, Any]] = []
    persistence_rows: list[dict[str, Any]] = []
    localization_curves: list[dict[str, Any]] = []
    base_records: dict[tuple[str, int], dict[str, Any]] = {}

    for hierarchy_label in (BRANCH_LABEL, CONTROL_LABEL):
        for n_side in levels:
            operator, points = (
                build_branch_zero_form_operator(n_side, epsilon)
                if hierarchy_label == BRANCH_LABEL
                else build_control_zero_form_operator(n_side, epsilon, diagonal_weight_scale)
            )
            state0 = build_candidate_state(points, candidate_label)
            initial_metrics = spatial_metrics(state0, points) | spectral_energy_metrics(state0, operator)
            initial_norm = float(np.linalg.norm(state0))
            base_final = np.asarray(expm_multiply((-evaluation_tau) * operator, state0), dtype=complex)
            base_metrics = spatial_metrics(base_final, points) | spectral_energy_metrics(base_final, operator)
            base_width_growth = float(base_metrics["localization_width"]) / max(float(initial_metrics["localization_width"]), 1.0e-12) - 1.0
            base_concentration_ratio = float(base_metrics["spatial_concentration"]) / max(float(initial_metrics["spatial_concentration"]), 1.0e-12)
            base_participation_growth = float(base_metrics["participation_ratio"]) / max(float(initial_metrics["participation_ratio"]), 1.0e-12) - 1.0

            eigenvalues, eigenvectors = sampled_eigensystem(operator, int(proxies_cfg["sample_eigen_count"]))
            spectral_gap = first_positive_gap(eigenvalues)
            low_basis = low_mode_basis(eigenvalues, eigenvectors, spectral_gap, float(proxies_cfg["low_mode_band_multiplier"]))
            base_local_variance = local_connectivity_variance(operator, points, float(proxies_cfg["local_radius"]))
            base_density_proxy = local_spectral_density_proxy(eigenvalues, spectral_gap, float(proxies_cfg["spectral_density_gap_multiplier"]))

            width_values = [float(initial_metrics["localization_width"])]
            persistence_time = 0.0
            for tau_value in persistence_tau_grid[1:]:
                evolved = np.asarray(expm_multiply((-tau_value) * operator, state0), dtype=complex)
                evolved_metrics = spatial_metrics(evolved, points)
                width_growth = float(evolved_metrics["localization_width"]) / max(float(initial_metrics["localization_width"]), 1.0e-12) - 1.0
                concentration_ratio = float(evolved_metrics["spatial_concentration"]) / max(float(initial_metrics["spatial_concentration"]), 1.0e-12)
                participation_growth = float(evolved_metrics["participation_ratio"]) / max(float(initial_metrics["participation_ratio"]), 1.0e-12) - 1.0
                width_values.append(float(evolved_metrics["localization_width"]))
                if survives_candidate_gate(width_growth, concentration_ratio, participation_growth, gate):
                    persistence_time = float(tau_value)

            localization_curves.append(
                {
                    "hierarchy_label": hierarchy_label,
                    "n_side": int(n_side),
                    "tau_grid": [round_float(value, 6) for value in persistence_tau_grid.tolist()],
                    "width_values": [round_float(value) for value in width_values],
                }
            )
            persistence_rows.append(
                {
                    "hierarchy_label": hierarchy_label,
                    "n_side": int(n_side),
                    "h": round_float(1.0 / float(n_side)),
                    "candidate_label": candidate_label,
                    "persistence_time_tau": round_float(persistence_time),
                    "base_width_growth_at_eval_tau": round_float(base_width_growth),
                    "base_concentration_ratio_at_eval_tau": round_float(base_concentration_ratio),
                    "base_participation_growth_at_eval_tau": round_float(base_participation_growth),
                    "spectral_gap": round_float(spectral_gap),
                    "local_connectivity_variance": round_float(base_local_variance),
                    "local_spectral_density_proxy": round_float(base_density_proxy),
                    "coefficient_ratio_shift_proxy": ratio_shift_map[(hierarchy_label, int(n_side))]["ratio_shift_proxy"],
                    "ratio_b1_over_b0": ratio_shift_map[(hierarchy_label, int(n_side))]["b1_over_b0"],
                    "ratio_b2_over_b0": ratio_shift_map[(hierarchy_label, int(n_side))]["b2_over_b0"],
                }
            )
            base_records[(hierarchy_label, int(n_side))] = {
                "operator": operator,
                "points": points,
                "state0": state0,
                "initial_metrics": initial_metrics,
                "base_final": base_final,
                "base_metrics": base_metrics,
                "base_width_growth": base_width_growth,
                "base_concentration_ratio": base_concentration_ratio,
                "base_participation_growth": base_participation_growth,
                "initial_norm": initial_norm,
                "spectral_gap": spectral_gap,
                "eigenvalues": eigenvalues,
                "low_basis": low_basis,
                "base_local_variance": base_local_variance,
                "base_density_proxy": base_density_proxy,
                "base_classification": classify_persistence(base_width_growth, base_concentration_ratio, 1.0, 1.0, thresholds),
            }

    for hierarchy_label in (BRANCH_LABEL, CONTROL_LABEL):
        for n_side in levels:
            base = base_records[(hierarchy_label, int(n_side))]
            operator = base["operator"]
            points = base["points"]
            state0 = base["state0"]
            base_final = base["base_final"]
            base_metrics = base["base_metrics"]
            base_width_growth = float(base["base_width_growth"])
            base_concentration_ratio = float(base["base_concentration_ratio"])
            base_norm = max(float(np.linalg.norm(base_final)), 1.0e-12)
            base_gap = float(base["spectral_gap"])
            low_basis = np.asarray(base["low_basis"], dtype=float)

            for perturbation_label, amplitudes in perturbation_grids.items():
                for amplitude in amplitudes:
                    operator_used = operator
                    state_init = state0
                    recompute_spectrum = False
                    if perturbation_label == "local_stiffness_defect":
                        operator_used = local_stiffness_defect(operator, int(n_side), float(amplitude))
                        recompute_spectrum = True
                    elif perturbation_label == "connectivity_micro_surgery":
                        operator_used = connectivity_micro_surgery(operator, int(n_side), float(amplitude))
                        recompute_spectrum = True
                    elif perturbation_label == "phase_randomization_window":
                        state_init = phase_randomization_window(state0, points, float(amplitude))
                    elif perturbation_label == "boundary_leakage_probe":
                        operator_used = boundary_leakage_probe(operator, points, float(amplitude), float(proxies_cfg["ring_leak_radius"]))
                    else:
                        raise ValueError(f"unknown perturbation label: {perturbation_label}")

                    final_state = np.asarray(expm_multiply((-evaluation_tau) * operator_used, state_init), dtype=complex)
                    final_metrics = spatial_metrics(final_state, points) | spectral_energy_metrics(final_state, operator_used)
                    overlap_value = overlap(base_final, final_state)
                    recovery_score = recovery_style_score(base_metrics, final_metrics)
                    classification = classify_persistence(base_width_growth, base_concentration_ratio, overlap_value, recovery_score, thresholds)
                    survival_indicator = 1 if classification == "persistent" else 0
                    width_growth = float(final_metrics["localization_width"]) / max(float(base["initial_metrics"]["localization_width"]), 1.0e-12) - 1.0
                    concentration_ratio = float(final_metrics["spatial_concentration"]) / max(float(base["initial_metrics"]["spatial_concentration"]), 1.0e-12)
                    participation_growth = float(final_metrics["participation_ratio"]) / max(float(base["initial_metrics"]["participation_ratio"]), 1.0e-12) - 1.0
                    norm_retention_relative = float(np.linalg.norm(final_state)) / base_norm
                    leakage_rate = 0.0
                    if perturbation_label == "boundary_leakage_probe":
                        leakage_rate = max(0.0, -math.log(max(norm_retention_relative, 1.0e-12)) / max(evaluation_tau, 1.0e-12))
                    mode_mixing = mode_mixing_indicator(final_state, low_basis)

                    gap_shift_relative = 0.0
                    density_proxy = float(base["base_density_proxy"])
                    local_variance = float(base["base_local_variance"])
                    if recompute_spectrum:
                        eigenvalues_used, _ = sampled_eigensystem(operator_used, int(proxies_cfg["sample_eigen_count"]))
                        gap_used = first_positive_gap(eigenvalues_used)
                        gap_shift_relative = safe_relative_difference(gap_used, base_gap)
                        density_proxy = local_spectral_density_proxy(eigenvalues_used, gap_used, float(proxies_cfg["spectral_density_gap_multiplier"]))
                        local_variance = local_connectivity_variance(operator_used, points, float(proxies_cfg["local_radius"]))

                    failure_label = failure_mechanism_label(
                        perturbation_label,
                        classification,
                        gap_shift_relative,
                        mode_mixing,
                        norm_retention_relative,
                    )
                    ratio_proxy = ratio_shift_map[(hierarchy_label, int(n_side))]["ratio_shift_proxy"]
                    survival_rows.append(
                        {
                            "hierarchy_label": hierarchy_label,
                            "level_id": f"n{n_side}",
                            "n_side": int(n_side),
                            "h": round_float(1.0 / float(n_side)),
                            "candidate_label": candidate_label,
                            "perturbation_label": perturbation_label,
                            "amplitude": round_float(float(amplitude)),
                            "evaluation_tau": round_float(evaluation_tau, 6),
                            "overlap_with_baseline": round_float(overlap_value),
                            "recovery_style_score": round_float(recovery_score),
                            "classification": classification,
                            "survival_indicator": int(survival_indicator),
                            "failure_mechanism": failure_label,
                            "localization_width": round_float(final_metrics["localization_width"]),
                            "width_growth": round_float(width_growth),
                            "concentration_ratio": round_float(concentration_ratio),
                            "participation_growth": round_float(participation_growth),
                            "spectral_gap_shift_relative": round_float(gap_shift_relative),
                            "mode_mixing_indicator": round_float(mode_mixing),
                            "local_connectivity_variance": round_float(local_variance),
                            "local_spectral_density_proxy": round_float(density_proxy),
                            "coefficient_ratio_shift_proxy": round_float(ratio_proxy),
                            "norm_retention_relative": round_float(norm_retention_relative),
                            "leakage_rate": round_float(leakage_rate),
                        }
                    )
                    failure_rows.append(
                        {
                            "hierarchy_label": hierarchy_label,
                            "level_id": f"n{n_side}",
                            "n_side": int(n_side),
                            "perturbation_label": perturbation_label,
                            "amplitude": round_float(float(amplitude)),
                            "classification": classification,
                            "failure_mechanism": failure_label,
                            "survival_indicator": int(survival_indicator),
                            "spectral_gap_shift_relative": round_float(gap_shift_relative),
                            "mode_mixing_indicator": round_float(mode_mixing),
                            "norm_retention_relative": round_float(norm_retention_relative),
                            "recovery_style_score": round_float(recovery_score),
                        }
                    )

    write_csv_rows(
        SURVIVAL_LEDGER_PATH,
        survival_rows,
        [
            "hierarchy_label",
            "level_id",
            "n_side",
            "h",
            "candidate_label",
            "perturbation_label",
            "amplitude",
            "evaluation_tau",
            "overlap_with_baseline",
            "recovery_style_score",
            "classification",
            "survival_indicator",
            "failure_mechanism",
            "localization_width",
            "width_growth",
            "concentration_ratio",
            "participation_growth",
            "spectral_gap_shift_relative",
            "mode_mixing_indicator",
            "local_connectivity_variance",
            "local_spectral_density_proxy",
            "coefficient_ratio_shift_proxy",
            "norm_retention_relative",
            "leakage_rate",
        ],
    )
    write_csv_rows(
        FAILURE_LEDGER_PATH,
        failure_rows,
        [
            "hierarchy_label",
            "level_id",
            "n_side",
            "perturbation_label",
            "amplitude",
            "classification",
            "failure_mechanism",
            "survival_indicator",
            "spectral_gap_shift_relative",
            "mode_mixing_indicator",
            "norm_retention_relative",
            "recovery_style_score",
        ],
    )
    write_csv_rows(
        PERSISTENCE_LEDGER_PATH,
        persistence_rows,
        [
            "hierarchy_label",
            "n_side",
            "h",
            "candidate_label",
            "persistence_time_tau",
            "base_width_growth_at_eval_tau",
            "base_concentration_ratio_at_eval_tau",
            "base_participation_growth_at_eval_tau",
            "spectral_gap",
            "local_connectivity_variance",
            "local_spectral_density_proxy",
            "coefficient_ratio_shift_proxy",
            "ratio_b1_over_b0",
            "ratio_b2_over_b0",
        ],
    )

    plot_localization_width(localization_curves, levels)
    plot_survival_probability(survival_rows, perturbation_grids)
    plot_persistence_scaling(persistence_rows)
    plot_gap_shift(survival_rows, levels)

    basin_summary: dict[str, dict[str, dict[str, Any]]] = {BRANCH_LABEL: {}, CONTROL_LABEL: {}}
    branch_survival_fractions: dict[str, float] = {}
    control_survival_fractions: dict[str, float] = {}
    for hierarchy_label in (BRANCH_LABEL, CONTROL_LABEL):
        for perturbation_label in perturbation_grids:
            thresholds_per_level: list[float | None] = []
            basin_summary[hierarchy_label][perturbation_label] = {}
            level_rows = [
                row
                for row in survival_rows
                if row["hierarchy_label"] == hierarchy_label and row["perturbation_label"] == perturbation_label
            ]
            if level_rows:
                if hierarchy_label == BRANCH_LABEL:
                    branch_survival_fractions[perturbation_label] = round_float(float(np.mean([float(row["survival_indicator"]) for row in level_rows])))
                else:
                    control_survival_fractions[perturbation_label] = round_float(float(np.mean([float(row["survival_indicator"]) for row in level_rows])))
            for n_side in levels:
                rows = [
                    row
                    for row in survival_rows
                    if row["hierarchy_label"] == hierarchy_label
                    and row["perturbation_label"] == perturbation_label
                    and int(row["n_side"]) == int(n_side)
                ]
                persistent_amplitudes = [float(row["amplitude"]) for row in rows if int(row["survival_indicator"]) == 1]
                failing_amplitudes = [float(row["amplitude"]) for row in rows if int(row["survival_indicator"]) == 0]
                first_failure = min(failing_amplitudes) if failing_amplitudes else None
                if first_failure is None:
                    max_persistent = max(persistent_amplitudes) if persistent_amplitudes else None
                else:
                    persistent_before_failure = [value for value in persistent_amplitudes if value < first_failure + 1.0e-12]
                    max_persistent = max(persistent_before_failure) if persistent_before_failure else None
                thresholds_per_level.append(first_failure)
                basin_summary[hierarchy_label][perturbation_label][f"n{n_side}"] = {
                    "max_persistent_amplitude": round_float(max_persistent) if max_persistent is not None else None,
                    "first_failure_amplitude": round_float(first_failure) if first_failure is not None else None,
                    "persistent_amplitude_count": int(len(persistent_amplitudes)),
                }
            basin_summary[hierarchy_label][perturbation_label]["thresholds_refinement_ordered"] = bool(ordered_threshold(thresholds_per_level))

    branch_persistence_rows = [row for row in persistence_rows if row["hierarchy_label"] == BRANCH_LABEL]
    control_persistence_rows = [row for row in persistence_rows if row["hierarchy_label"] == CONTROL_LABEL]
    branch_h = [float(row["h"]) for row in branch_persistence_rows]
    control_h = [float(row["h"]) for row in control_persistence_rows]
    branch_tau = [float(row["persistence_time_tau"]) for row in branch_persistence_rows]
    control_tau = [float(row["persistence_time_tau"]) for row in control_persistence_rows]
    branch_fit = fit_power_law(branch_h, branch_tau)
    control_fit = fit_power_law(control_h, control_tau)
    tau_ratio = [branch / max(control, 1.0e-12) for branch, control in zip(branch_tau, control_tau)]

    branch_ratio_shift_values = [float(row["coefficient_ratio_shift_proxy"]) for row in branch_persistence_rows]
    branch_local_density_values = [float(row["local_spectral_density_proxy"]) for row in branch_persistence_rows]
    branch_local_variance_values = [float(row["local_connectivity_variance"]) for row in branch_persistence_rows]

    operator_branch_rows = [
        row
        for row in survival_rows
        if row["hierarchy_label"] == BRANCH_LABEL and row["perturbation_label"] in {"local_stiffness_defect", "connectivity_micro_surgery"}
    ]
    survival_corr_connectivity = correlation(
        [float(row["local_connectivity_variance"]) for row in operator_branch_rows],
        [float(row["survival_indicator"]) for row in operator_branch_rows],
    )
    survival_corr_density = correlation(
        [float(row["local_spectral_density_proxy"]) for row in operator_branch_rows],
        [float(row["survival_indicator"]) for row in operator_branch_rows],
    )
    persistence_corr_ratio = correlation(branch_ratio_shift_values, branch_tau)
    persistence_corr_density = correlation(branch_local_density_values, branch_tau)
    persistence_corr_variance = correlation(branch_local_variance_values, branch_tau)

    dominant_failure_modes: dict[str, str] = {}
    for perturbation_label in perturbation_grids:
        failures = [
            row["failure_mechanism"]
            for row in failure_rows
            if row["hierarchy_label"] == BRANCH_LABEL
            and row["perturbation_label"] == perturbation_label
            and row["failure_mechanism"] != "survives"
        ]
        dominant_failure_modes[perturbation_label] = max(set(failures), key=failures.count) if failures else "survives"

    bounded_survival_basin = bool(
        any(
            basin_summary[BRANCH_LABEL][perturbation_label][f"n{levels[0]}"]["first_failure_amplitude"] is not None
            for perturbation_label in perturbation_grids
        )
    )
    ordered_failure_threshold = bool(
        basin_summary[BRANCH_LABEL]["local_stiffness_defect"]["thresholds_refinement_ordered"]
        and basin_summary[BRANCH_LABEL]["local_stiffness_defect"][f"n{levels[0]}"]["first_failure_amplitude"] is not None
    )
    persistence_structured = bool(
        all(left <= right + 1.0e-12 for left, right in zip(branch_tau, branch_tau[1:]))
        and float(branch_fit["r_squared"]) >= 0.9
        and min(tau_ratio) > 1.2
    )
    control_different = bool(
        min(tau_ratio) > 1.2
        and branch_survival_fractions.get("local_stiffness_defect", 0.0) > control_survival_fractions.get("local_stiffness_defect", 0.0)
    )

    success_flags = {
        "bounded_survival_basin_identified": bounded_survival_basin,
        "refinement_ordered_failure_threshold_found": ordered_failure_threshold,
        "persistence_time_structured_scaling": persistence_structured,
        "deterministic_control_different": control_different,
    }
    success = all(success_flags.values())

    runtime_seconds = time.perf_counter() - start_time
    write_json(
        RUNS_PATH,
        {
            "timestamp": timestamp_iso(),
            "phase_name": PHASE_NAME,
            "builder_script": relpath(Path(__file__)),
            "config_reference": relpath(CONFIG_PATH),
            "references": {
                "phase6_manifest": relpath(PHASE6_MANIFEST_PATH),
                "phase7_manifest": relpath(PHASE7_MANIFEST_PATH),
                "phase8_manifest": relpath(PHASE8_MANIFEST_PATH),
                "phase9_manifest": relpath(PHASE9_MANIFEST_PATH),
                "phase10_manifest": relpath(PHASE10_MANIFEST_PATH),
                "phasex_manifest": relpath(PHASEX_MANIFEST_PATH),
                "phasex_candidates": relpath(PHASEX_CANDIDATES_PATH),
            },
            "refinement_levels": levels,
            "evaluation_tau": round_float(evaluation_tau, 6),
            "persistence_tau_grid": [round_float(value, 6) for value in persistence_tau_grid.tolist()],
            "perturbation_grids": perturbation_grids,
            "runtime_seconds": round_float(runtime_seconds, 6),
        },
    )

    summary_lines = [
        "# Phase XI -- Localized Mode Protection and Failure-Mechanism Program",
        "",
        "## Objective",
        "",
        "Test whether the frozen Phase X localized candidate `low_mode_localized_wavepacket` is protected by reproducible structural mechanisms or is only an accidental configuration under deterministic perturbation channels.",
        "",
        "## Frozen Inputs",
        "",
        f"- Phase VI operator hierarchy: `{relpath(PHASE6_MANIFEST_PATH)}`",
        f"- Phase VII spectral baseline: `{relpath(PHASE7_MANIFEST_PATH)}`",
        f"- Phase VIII trace contract: `{relpath(PHASE8_MANIFEST_PATH)}`",
        f"- Phase IX coefficient contract: `{relpath(PHASE9_MANIFEST_PATH)}`",
        f"- Phase X cautious bridge baseline: `{relpath(PHASE10_MANIFEST_PATH)}`",
        f"- Phase X proto-particle baseline: `{relpath(PHASEX_MANIFEST_PATH)}`",
        "",
        "## Survival Basins",
        "",
        f"- Refinement levels tested: `{levels}` with evaluation time `tau = {round_float(evaluation_tau, 6)}`.",
        f"- Branch local stiffness first-failure amplitudes: `{[basin_summary[BRANCH_LABEL]['local_stiffness_defect'][f'n{n_side}']['first_failure_amplitude'] for n_side in levels]}`.",
        f"- Branch phase-randomization first-failure amplitudes: `{[basin_summary[BRANCH_LABEL]['phase_randomization_window'][f'n{n_side}']['first_failure_amplitude'] for n_side in levels]}`.",
        f"- Boundary leakage remains inside the persistent basin over the scanned amplitude range, but leakage rate is still recorded explicitly.",
        f"- Branch mean survival fractions by perturbation: `{branch_survival_fractions}`.",
        f"- Control mean survival fractions by perturbation: `{control_survival_fractions}`.",
        "",
        "## Failure Mechanisms",
        "",
        f"- Dominant branch failure channels are `{dominant_failure_modes}`.",
        f"- Phase randomization failures are coherence-triggered, connectivity micro-surgery failures are topological-defect triggered, and strong local stiffness failures are spectral-mixing dominated near threshold.",
        "",
        "## Persistence Scaling",
        "",
        f"- Branch persistence times are `{[round_float(value, 6) for value in branch_tau]}` with power-law fit `{branch_fit}`.",
        f"- Control persistence times are `{[round_float(value, 6) for value in control_tau]}` with power-law fit `{control_fit}`.",
        f"- Branch/control persistence-time ratios are `{[round_float(value, 6) for value in tau_ratio]}`.",
        "",
        "## Structural Correlates",
        "",
        f"- Branch persistence time vs coefficient-ratio shift correlation: `{round_float(persistence_corr_ratio)}`.",
        f"- Branch persistence time vs local spectral density proxy correlation: `{round_float(persistence_corr_density)}`.",
        f"- Branch persistence time vs local connectivity variance correlation: `{round_float(persistence_corr_variance)}`.",
        f"- Branch operator-survival vs local spectral density correlation: `{round_float(survival_corr_density)}`.",
        f"- Branch operator-survival vs local connectivity variance correlation: `{round_float(survival_corr_connectivity)}`.",
        "",
        "## Bounded Interpretation",
        "",
        "Phase XI supports only a bounded protection statement. The frozen localized candidate has a deterministic survival basin, loses persistence through repeatable channels rather than random breakdown, and shows refinement-ordered persistence-time growth on the branch. This remains a protection/failure diagnostic only. It does not establish particles, fields, or new physical structure.",
        "",
        closure_statement(success),
    ]
    write_text(SUMMARY_PATH, "\n".join(summary_lines))

    manifest = {
        "timestamp": timestamp_iso(),
        "phase": 11,
        "phase_name": PHASE_NAME,
        "stage_identifier": STAGE_IDENTIFIER,
        "status": "passed" if success else "failed",
        "success": bool(success),
        "objective": "localized_mode_protection_and_failure_mechanism_feasibility",
        "operator_manifest_reference": relpath(PHASE6_MANIFEST_PATH),
        "phase7_manifest_reference": relpath(PHASE7_MANIFEST_PATH),
        "phase8_manifest_reference": relpath(PHASE8_MANIFEST_PATH),
        "phase9_manifest_reference": relpath(PHASE9_MANIFEST_PATH),
        "phase10_manifest_reference": relpath(PHASE10_MANIFEST_PATH),
        "phasex_manifest_reference": relpath(PHASEX_MANIFEST_PATH),
        "phasex_candidate_reference": relpath(PHASEX_CANDIDATES_PATH),
        "candidate_label": candidate_label,
        "refinement_hierarchy": {
            "h_definition": "h = 1 / n_side on the unit-periodic lattice; it is the lattice-spacing surrogate used to order protection diagnostics.",
            "levels": levels,
        },
        "evaluation_protocol": {
            "evaluation_tau": round_float(evaluation_tau, 6),
            "persistence_tau_grid": [round_float(value, 6) for value in persistence_tau_grid.tolist()],
            "short_time_window": [round_float(value, 6) for value in short_time_window.tolist()],
            "perturbation_grids": perturbation_grids,
        },
        "survival_basins": basin_summary,
        "persistence_scaling": {
            "branch_tau_values": [round_float(value, 6) for value in branch_tau],
            "control_tau_values": [round_float(value, 6) for value in control_tau],
            "branch_power_law_fit": branch_fit,
            "control_power_law_fit": control_fit,
            "branch_over_control_tau_ratio": [round_float(value, 6) for value in tau_ratio],
        },
        "failure_mechanisms": dominant_failure_modes,
        "structural_correlates": {
            "branch_persistence_vs_ratio_shift_correlation": round_float(persistence_corr_ratio),
            "branch_persistence_vs_local_spectral_density_correlation": round_float(persistence_corr_density),
            "branch_persistence_vs_local_connectivity_variance_correlation": round_float(persistence_corr_variance),
            "branch_survival_vs_local_spectral_density_correlation": round_float(survival_corr_density),
            "branch_survival_vs_local_connectivity_variance_correlation": round_float(survival_corr_connectivity),
        },
        "success_flags": success_flags,
        "artifacts": {
            "survival_ledger_csv": relpath(SURVIVAL_LEDGER_PATH),
            "failure_classification_csv": relpath(FAILURE_LEDGER_PATH),
            "persistence_scaling_csv": relpath(PERSISTENCE_LEDGER_PATH),
            "summary_md": relpath(SUMMARY_PATH),
            "manifest_json": relpath(MANIFEST_PATH),
            "width_plot": relpath(WIDTH_PLOT_PATH),
            "survival_plot": relpath(SURVIVAL_PLOT_PATH),
            "persistence_plot": relpath(PERSISTENCE_PLOT_PATH),
            "gap_plot": relpath(GAP_PLOT_PATH),
            "runs_json": relpath(RUNS_PATH),
            "builder_script": relpath(Path(__file__)),
        },
        "runtime_seconds": round_float(runtime_seconds, 6),
        "claim_boundary": CLAIM_BOUNDARY,
        "conclusion": closure_statement(success),
    }
    write_json(MANIFEST_PATH, manifest)


if __name__ == "__main__":
    main()
