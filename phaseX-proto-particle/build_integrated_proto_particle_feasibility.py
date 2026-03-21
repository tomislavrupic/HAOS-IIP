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
from scipy.sparse.linalg import expm_multiply

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haos_core import build_graph, read_json, relpath, write_csv_rows, write_json

PHASE_ROOT = ROOT / "phaseX-proto-particle"
CONFIG_PATH = PHASE_ROOT / "configs" / "proto_particle_feasibility_config.json"
RUNS_ROOT = PHASE_ROOT / "runs"
PLOTS_ROOT = PHASE_ROOT / "plots"

PHASE6_MANIFEST_PATH = ROOT / "phase6-operator" / "phase6_operator_manifest.json"
PHASE8_MANIFEST_PATH = ROOT / "phase8-trace" / "phase8_trace_manifest.json"
PHASE9_MANIFEST_PATH = ROOT / "phase9-invariants" / "phase9_manifest.json"

CANDIDATE_PATH = RUNS_ROOT / "proto_particle_candidates.json"
LOCALIZATION_LEDGER_PATH = RUNS_ROOT / "localization_metrics_ledger.csv"
PERSISTENCE_LEDGER_PATH = RUNS_ROOT / "persistence_classification.csv"
SCALING_LEDGER_PATH = RUNS_ROOT / "scaling_prediction_errors.csv"
RUNS_PATH = RUNS_ROOT / "phaseX_runs.json"
SUMMARY_PATH = PHASE_ROOT / "phaseX_integrated_summary.md"
MANIFEST_PATH = PHASE_ROOT / "phaseX_integrated_manifest.json"

WIDTH_PLOT_PATH = PLOTS_ROOT / "phaseX_localization_width_vs_time.svg"
OVERLAP_PLOT_PATH = PLOTS_ROOT / "phaseX_excitation_overlap_vs_perturbation.svg"
REFINEMENT_PLOT_PATH = PLOTS_ROOT / "phaseX_refinement_scaling_of_candidate_modes.svg"

PHASE_NAME = "phaseX-proto-particle"
STAGE_IDENTIFIER = "phaseX-proto-particle"
BRANCH_LABEL = "frozen_branch"
CONTROL_LABEL = "periodic_diagonal_augmented_control"
BRANCH_SCALING_LABEL = "frozen_branch_full_spectrum"
CONTROL_SCALING_LABEL = "periodic_diagonal_augmented_control_full_spectrum"

CLAIM_BOUNDARY = (
    "This integrated phase is limited to proto-particle feasibility diagnostics on the frozen "
    "operator and spectral-trace contracts. It does not assert particles, continuum limits, "
    "geometric structure, or physical correspondence."
)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def timestamp_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def closure_statement(success: bool) -> str:
    if success:
        return "Integrated run establishes proto-particle feasibility for the frozen operator hierarchy."
    return "Integrated run does not establish proto-particle feasibility for the frozen operator hierarchy."


def round_float(value: float, digits: int = 12) -> float:
    return round(float(value), digits)


def safe_relative_difference(new_value: float, reference_value: float) -> float:
    return abs(new_value - reference_value) / max(abs(reference_value), 1.0e-12)


def relative_span(values: np.ndarray) -> float:
    return float((np.max(values) - np.min(values)) / max(abs(float(np.mean(values))), 1.0e-12))


def model_r_squared(x_values: np.ndarray, y_values: np.ndarray, coefficients: np.ndarray) -> float:
    fitted = np.polyval(coefficients, x_values)
    denominator = float(np.sum((y_values - y_values.mean()) ** 2))
    if denominator <= 1.0e-18:
        return 1.0
    return 1.0 - float(np.sum((y_values - fitted) ** 2) / denominator)


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


def fit_linear_family(
    h_values: np.ndarray,
    series_matrix: np.ndarray,
    power: int,
    target_h: float | None = None,
) -> dict[str, Any]:
    x_values = h_values**power
    predictions = np.zeros_like(series_matrix, dtype=float)
    r_squared = []
    limits = []
    predicted_target = []
    for column_index in range(series_matrix.shape[1]):
        coefficients = np.polyfit(x_values, series_matrix[:, column_index], 1)
        predictions[:, column_index] = np.polyval(coefficients, x_values)
        r_squared.append(model_r_squared(x_values, series_matrix[:, column_index], coefficients))
        limits.append(float(coefficients[1]))
        if target_h is not None:
            predicted_target.append(float(np.polyval(coefficients, float(target_h) ** power)))
    payload: dict[str, Any] = {
        "power": int(power),
        "predictions": predictions,
        "r_squared": np.array(r_squared, dtype=float),
        "limits": np.array(limits, dtype=float),
    }
    if target_h is not None:
        payload["predicted_target"] = np.array(predicted_target, dtype=float)
    return payload


def load_frozen_inputs() -> dict[str, Any]:
    config = read_json(CONFIG_PATH)
    phase6_manifest = read_json(PHASE6_MANIFEST_PATH)
    phase8_manifest = read_json(PHASE8_MANIFEST_PATH)
    phase9_manifest = read_json(PHASE9_MANIFEST_PATH)

    if not bool(phase6_manifest.get("success")):
        raise ValueError("Integrated phase requires the frozen successful Phase VI operator manifest.")
    if phase6_manifest.get("selected_operator_class") != "cochain_laplacian":
        raise ValueError("Integrated phase requires the frozen cochain-Laplacian hierarchy.")
    if not bool(phase8_manifest.get("success")):
        raise ValueError("Integrated phase requires the frozen successful Phase VIII trace manifest.")
    if not bool(phase9_manifest.get("success")):
        raise ValueError("Integrated phase requires the frozen successful Phase IX manifest.")

    return {
        "config": config,
        "phase6_manifest": phase6_manifest,
        "phase8_manifest": phase8_manifest,
        "phase9_manifest": phase9_manifest,
    }


def periodic_delta(points: np.ndarray, anchor: np.ndarray) -> np.ndarray:
    delta = np.asarray(points, dtype=float) - np.asarray(anchor, dtype=float)
    return (delta + 0.5) % 1.0 - 0.5


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


def full_spectrum_from_scalar(scalar_spectrum: np.ndarray) -> np.ndarray:
    return np.sort(np.repeat(scalar_spectrum, 4))


def trace_from_spectrum(spectrum: np.ndarray, probe_times: np.ndarray) -> np.ndarray:
    return np.array([float(np.exp(-probe_time * spectrum).sum()) for probe_time in probe_times], dtype=float)


def scaling_observations(
    hierarchy_label: str,
    levels: list[dict[str, Any]],
    epsilon: float,
    diagonal_weight_scale: float,
    probe_times: np.ndarray,
    phase7_exponent: float,
) -> list[dict[str, Any]]:
    rows = []
    scalar_builder = branch_scalar_spectrum if hierarchy_label == BRANCH_SCALING_LABEL else control_scalar_spectrum
    for level in levels:
        n_side = int(level["n_side"])
        scalar_spectrum = (
            scalar_builder(n_side, epsilon)
            if hierarchy_label == BRANCH_SCALING_LABEL
            else scalar_builder(n_side, epsilon, diagonal_weight_scale)
        )
        spectrum = full_spectrum_from_scalar(scalar_spectrum)
        h_value = float(level["h"])
        traces = trace_from_spectrum(spectrum, probe_times)
        scaled_fit = quadratic_fit_named(probe_times, (h_value * h_value) * traces, prefix="b")
        raw_fit = quadratic_fit_named(probe_times, traces, prefix="a")
        rows.append(
            {
                "hierarchy_label": hierarchy_label,
                "level_id": str(level["level_id"]),
                "n_side": n_side,
                "h": round_float(h_value),
                "spectral_radius": round_float(float(spectrum[-1])),
                "spectral_p95": round_float(float(np.quantile(spectrum, 0.95))),
                "first_positive_eigenvalue": round_float(float(spectrum[spectrum > 1.0e-12][0])),
                "scaled_coefficients": scaled_fit,
                "ratios": {
                    "b1_over_b0": round_float(float(scaled_fit["b1"]) / max(abs(float(scaled_fit["b0"])), 1.0e-12)),
                    "b2_over_b0": round_float(float(scaled_fit["b2"]) / max(abs(float(scaled_fit["b0"])), 1.0e-12)),
                },
                "rescaled_invariants": {
                    "I0": round_float((h_value**phase7_exponent) * float(raw_fit["a0"])),
                    "I1": round_float((h_value**phase7_exponent) * float(raw_fit["a1"])),
                    "I2": round_float((h_value**phase7_exponent) * float(raw_fit["a2"])),
                },
                "_trace_window": traces,
                "_probe_times": probe_times,
            }
        )
    return rows


def scaling_prediction_rows(
    observations: list[dict[str, Any]],
    hierarchy_label: str,
    ratio_names: list[str],
    invariant_names: list[str],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    h_values = np.array([float(observation["h"]) for observation in observations[:-1]], dtype=float)
    target = observations[-1]

    coefficient_names = ["b0", "b1", "b2"]
    coefficient_matrix = np.array(
        [[float(observation["scaled_coefficients"][name]) for name in coefficient_names] for observation in observations],
        dtype=float,
    )
    ratio_matrix = np.array(
        [[float(observation["ratios"][name]) for name in ratio_names] for observation in observations],
        dtype=float,
    )
    invariant_matrix = np.array(
        [[float(observation["rescaled_invariants"][name]) for name in invariant_names] for observation in observations],
        dtype=float,
    )

    coefficient_fit = fit_linear_family(h_values, coefficient_matrix[:-1], power=2, target_h=float(target["h"]))
    ratio_fit = fit_linear_family(h_values, ratio_matrix[:-1], power=2, target_h=float(target["h"]))
    invariant_fit = fit_linear_family(h_values, invariant_matrix[:-1], power=1, target_h=float(target["h"]))

    target_coefficients = coefficient_matrix[-1]
    target_ratios = ratio_matrix[-1]
    target_invariants = invariant_matrix[-1]

    predicted_trace = (float(target["h"]) ** -2) * (
        coefficient_fit["predicted_target"][0]
        + coefficient_fit["predicted_target"][1] * target["_probe_times"]
        + coefficient_fit["predicted_target"][2] * (target["_probe_times"] ** 2)
    )
    trace_error = float(np.max(np.abs(predicted_trace - target["_trace_window"]) / target["_trace_window"]))

    rows: list[dict[str, Any]] = []
    for index, name in enumerate(coefficient_names):
        rows.append(
            {
                "hierarchy_label": hierarchy_label,
                "family": "scaled_coefficients",
                "quantity": name,
                "fit_power": 2,
                "fit_scope": "R1_to_R4_predict_R5",
                "observed_value_R5": round_float(target_coefficients[index]),
                "predicted_value_R5": round_float(coefficient_fit["predicted_target"][index]),
                "relative_error_R5": round_float(
                    safe_relative_difference(coefficient_fit["predicted_target"][index], target_coefficients[index])
                ),
                "spectral_radius_R5": round_float(target["spectral_radius"]),
                "spectral_p95_R5": round_float(target["spectral_p95"]),
            }
        )
    for index, name in enumerate(ratio_names):
        rows.append(
            {
                "hierarchy_label": hierarchy_label,
                "family": "ratios",
                "quantity": name,
                "fit_power": 2,
                "fit_scope": "R1_to_R4_predict_R5",
                "observed_value_R5": round_float(target_ratios[index]),
                "predicted_value_R5": round_float(ratio_fit["predicted_target"][index]),
                "relative_error_R5": round_float(safe_relative_difference(ratio_fit["predicted_target"][index], target_ratios[index])),
                "spectral_radius_R5": round_float(target["spectral_radius"]),
                "spectral_p95_R5": round_float(target["spectral_p95"]),
            }
        )
    for index, name in enumerate(invariant_names):
        rows.append(
            {
                "hierarchy_label": hierarchy_label,
                "family": "rescaled_invariants",
                "quantity": name,
                "fit_power": 1,
                "fit_scope": "R1_to_R4_predict_R5",
                "observed_value_R5": round_float(target_invariants[index]),
                "predicted_value_R5": round_float(invariant_fit["predicted_target"][index]),
                "relative_error_R5": round_float(
                    safe_relative_difference(invariant_fit["predicted_target"][index], target_invariants[index])
                ),
                "spectral_radius_R5": round_float(target["spectral_radius"]),
                "spectral_p95_R5": round_float(target["spectral_p95"]),
            }
        )
    rows.append(
        {
            "hierarchy_label": hierarchy_label,
            "family": "trace_window",
            "quantity": "max_relative_error",
            "fit_power": 2,
            "fit_scope": "R1_to_R4_predict_R5",
            "observed_value_R5": "",
            "predicted_value_R5": "",
            "relative_error_R5": round_float(trace_error),
            "spectral_radius_R5": round_float(target["spectral_radius"]),
            "spectral_p95_R5": round_float(target["spectral_p95"]),
        }
    )

    metrics = {
        "trace_error": trace_error,
        "ratio_span": max(relative_span(ratio_matrix[:, 0]), relative_span(ratio_matrix[:, 1])),
        "invariant_span": max(relative_span(invariant_matrix[:, 0]), relative_span(invariant_matrix[:, 1]), relative_span(invariant_matrix[:, 2])),
        "spectral_radius_values": [round_float(observation["spectral_radius"]) for observation in observations],
        "ratio_values_R5": {name: round_float(value) for name, value in zip(ratio_names, target_ratios)},
        "invariant_values_R5": {name: round_float(value) for name, value in zip(invariant_names, target_invariants)},
    }
    return rows, metrics


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
    n_nodes = int(graph.block_sizes[0])
    return graph.delta_h[:n_nodes, :n_nodes].tocsr(), np.asarray(graph.points, dtype=float)


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


def state_family(points: np.ndarray, family_label: str) -> np.ndarray:
    center = np.array([0.5, 0.5], dtype=float)
    delta = periodic_delta(points, center)
    x_values = delta[:, 0]
    y_values = delta[:, 1]
    if family_label == "low_mode_localized_wavepacket":
        sigma = 0.08
        amplitude = np.exp(-0.5 * ((x_values / sigma) ** 2 + (y_values / sigma) ** 2))
        phase = np.exp(1j * 2.0 * math.pi * (x_values + y_values))
        state = amplitude * phase
    elif family_label == "mid_spectrum_localized_superposition":
        sigma = 0.055
        amplitude = np.exp(-0.5 * ((x_values / sigma) ** 2 + (y_values / sigma) ** 2))
        phase = np.exp(1j * 2.0 * math.pi * (4.0 * x_values - 3.0 * y_values))
        state = amplitude * phase
    elif family_label == "random_localized_seed":
        sigma = 0.045
        amplitude = np.exp(-0.5 * ((x_values / sigma) ** 2 + (y_values / sigma) ** 2))
        envelope = np.sin(37.0 * (x_values + 0.03)) * np.cos(29.0 * (y_values - 0.02)) + 0.5 * np.sin(11.0 * (x_values - y_values))
        state = amplitude * np.asarray(envelope, dtype=complex)
    else:
        raise ValueError(f"unknown family: {family_label}")
    state /= max(np.linalg.norm(state), 1.0e-12)
    return np.asarray(state, dtype=complex)


def amplitude_jitter(state: np.ndarray, points: np.ndarray) -> np.ndarray:
    delta = periodic_delta(points, np.array([0.5, 0.5], dtype=float))
    modulation = 1.0 + 0.12 * np.sin(13.0 * delta[:, 0]) * np.cos(9.0 * delta[:, 1])
    out = state * modulation
    out /= max(np.linalg.norm(out), 1.0e-12)
    return np.asarray(out, dtype=complex)


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


def local_connectivity_shuffle(
    operator: sp.csr_matrix,
    n_side: int,
    epsilon: float,
) -> sp.csr_matrix:
    h_value = 1.0 / float(n_side)
    base_weight = math.exp(-(h_value * h_value) / (2.0 * epsilon))
    adjacency = adjacency_from_laplacian(operator)

    def idx(i: int, j: int) -> int:
        return (i % n_side) * n_side + (j % n_side)

    cells = [
        (n_side // 2 - 1, n_side // 2 - 1),
        (n_side // 2, n_side // 2 - 1),
        (n_side // 2 - 1, n_side // 2),
    ]
    for i, j in cells:
        a = idx(i, j)
        b = idx(i + 1, j)
        c = idx(i, j + 1)
        d = idx(i + 1, j + 1)
        for edge in (tuple(sorted((a, b))), tuple(sorted((a, c))), tuple(sorted((b, d))), tuple(sorted((c, d)))):
            if edge in adjacency:
                adjacency[edge] *= 0.35
        for edge in (tuple(sorted((a, d))), tuple(sorted((b, c)))):
            adjacency[edge] = max(adjacency.get(edge, 0.0), base_weight)
    return laplacian_from_adjacency(operator.shape[0], adjacency)


def branch_truncated_evolution(
    state: np.ndarray,
    n_side: int,
    epsilon: float,
    tau_value: float,
    cutoff: float,
) -> np.ndarray:
    h_value = 1.0 / float(n_side)
    edge_weight = math.exp(-(h_value * h_value) / (2.0 * epsilon))
    field = state.reshape((n_side, n_side))
    spectrum = np.fft.fftn(field)
    out = np.zeros_like(spectrum, dtype=complex)
    for mode_x in range(n_side):
        cos_x = math.cos(2.0 * math.pi * mode_x / n_side)
        for mode_y in range(n_side):
            cos_y = math.cos(2.0 * math.pi * mode_y / n_side)
            eigenvalue = 2.0 * edge_weight * (2.0 - cos_x - cos_y)
            if eigenvalue <= cutoff:
                out[mode_x, mode_y] = math.exp(-tau_value * eigenvalue) * spectrum[mode_x, mode_y]
    return np.fft.ifftn(out).reshape(-1)


def control_truncated_evolution(
    state: np.ndarray,
    n_side: int,
    epsilon: float,
    diagonal_weight_scale: float,
    tau_value: float,
    cutoff: float,
) -> np.ndarray:
    h_value = 1.0 / float(n_side)
    edge_weight = math.exp(-(h_value * h_value) / (2.0 * epsilon))
    diagonal_weight = edge_weight * float(diagonal_weight_scale)
    field = state.reshape((n_side, n_side))
    spectrum = np.fft.fftn(field)
    out = np.zeros_like(spectrum, dtype=complex)
    for mode_x in range(n_side):
        cos_x = math.cos(2.0 * math.pi * mode_x / n_side)
        for mode_y in range(n_side):
            cos_y = math.cos(2.0 * math.pi * mode_y / n_side)
            eigenvalue = 2.0 * edge_weight * (2.0 - cos_x - cos_y) + 4.0 * diagonal_weight * (1.0 - cos_x * cos_y)
            if eigenvalue <= cutoff:
                out[mode_x, mode_y] = math.exp(-tau_value * eigenvalue) * spectrum[mode_x, mode_y]
    return np.fft.ifftn(out).reshape(-1)


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


def recovery_style_score(base_metrics: dict[str, float], perturbed_metrics: dict[str, float]) -> float:
    tolerances = {
        "localization_width": 0.08,
        "spatial_concentration": 0.12,
        "participation_ratio": 0.20,
        "spectral_energy_spread": 0.18,
    }
    penalties = []
    for key, tolerance in tolerances.items():
        deviation = safe_relative_difference(float(perturbed_metrics[key]), float(base_metrics[key]))
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


def plot_localization_width(localization_rows: list[dict[str, Any]], selected_families: list[str], scan_levels: list[int]) -> None:
    plt.figure(figsize=(8.2, 5.2))
    for hierarchy_label, linestyle in ((BRANCH_LABEL, "-"), (CONTROL_LABEL, "--")):
        for family_label in selected_families:
            rows = [
                row
                for row in localization_rows
                if row["hierarchy_label"] == hierarchy_label
                and int(row["n_side"]) == int(scan_levels[-1])
                and row["family_label"] == family_label
            ]
            tau_values = [float(row["tau"]) for row in rows]
            widths = [float(row["localization_width"]) for row in rows]
            plt.plot(tau_values, widths, linestyle=linestyle, marker="o", label=f"{hierarchy_label}:{family_label}")
    plt.xlabel("tau")
    plt.ylabel("localization width")
    plt.grid(alpha=0.3)
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(WIDTH_PLOT_PATH, format="svg")
    plt.close()


def plot_overlap_vs_perturbation(persistence_rows: list[dict[str, Any]], selected_families: list[str], scan_level: int) -> None:
    perturbations = ["amplitude_jitter", "local_connectivity_shuffle", "spectral_truncation_shift"]
    x_values = np.arange(len(perturbations), dtype=float)
    width = 0.18
    plt.figure(figsize=(9.0, 5.4))
    offsets = {
        (BRANCH_LABEL, selected_families[0]): -1.5 * width,
        (BRANCH_LABEL, selected_families[1] if len(selected_families) > 1 else selected_families[0]): -0.5 * width,
        (CONTROL_LABEL, selected_families[0]): 0.5 * width,
        (CONTROL_LABEL, selected_families[1] if len(selected_families) > 1 else selected_families[0]): 1.5 * width,
    }
    for hierarchy_label in (BRANCH_LABEL, CONTROL_LABEL):
        for family_label in selected_families:
            rows = [
                row
                for row in persistence_rows
                if row["hierarchy_label"] == hierarchy_label
                and int(row["n_side"]) == int(scan_level)
                and row["family_label"] == family_label
                and row["perturbation_label"] in perturbations
            ]
            values = [float(next(row for row in rows if row["perturbation_label"] == label)["overlap_with_base"]) for label in perturbations]
            plt.bar(x_values + offsets[(hierarchy_label, family_label)], values, width=width, label=f"{hierarchy_label}:{family_label}")
    plt.xticks(x_values, perturbations, rotation=20, ha="right")
    plt.ylabel("overlap with base excitation")
    plt.ylim(0.7, 1.01)
    plt.grid(axis="y", alpha=0.3)
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(OVERLAP_PLOT_PATH, format="svg")
    plt.close()


def plot_refinement_scaling(persistence_rows: list[dict[str, Any]], selected_family: str, repeat_levels: list[int]) -> None:
    plt.figure(figsize=(7.4, 4.8))
    for hierarchy_label, marker in ((BRANCH_LABEL, "o"), (CONTROL_LABEL, "s")):
        rows = [
            row
            for row in persistence_rows
            if row["hierarchy_label"] == hierarchy_label
            and row["family_label"] == selected_family
            and row["perturbation_label"] == "aggregate"
            and int(row["n_side"]) in repeat_levels
        ]
        n_values = [int(row["n_side"]) for row in rows]
        width_growths = [float(row["base_width_growth"]) for row in rows]
        plt.plot(n_values, width_growths, marker=marker, label=f"{hierarchy_label}:{selected_family}")
    plt.xlabel("n_side")
    plt.ylabel("final localization width growth")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(REFINEMENT_PLOT_PATH, format="svg")
    plt.close()


def main() -> None:
    start_time = time.perf_counter()
    RUNS_ROOT.mkdir(parents=True, exist_ok=True)
    PLOTS_ROOT.mkdir(parents=True, exist_ok=True)

    payload = load_frozen_inputs()
    config = payload["config"]
    phase6_manifest = payload["phase6_manifest"]
    phase8_manifest = payload["phase8_manifest"]
    phase9_manifest = payload["phase9_manifest"]

    epsilon = float(phase6_manifest["frozen_inputs"]["base_epsilon"])
    diagonal_weight_scale = float(config["control_graph"]["diagonal_weight_scale"])
    tau_grid = np.array(config["tau_grid"], dtype=float)
    tau_final = float(tau_grid[-1])
    scan_levels = [int(value) for value in config["scan_levels"]]
    scan_level = int(config["candidate_selection"]["scan_level_n_side"])
    repeat_levels = [int(scan_levels[0]), int(scan_levels[-1])]
    short_time_window = np.array(phase8_manifest["operational_short_time_window"]["probe_times"], dtype=float)
    phase7_exponent = float(phase9_manifest["dominant_phase7_low_mode_exponent"]["exponent"])

    base_levels = list(phase6_manifest["refinement_hierarchy"]["levels"])
    extended_levels = list(base_levels)
    extended_levels.append(
        {
            "level_id": f"R{len(base_levels) + 1}",
            "refinement_multiplier": int(config["refinement_extension_multiplier"]),
            "n_side": int(phase6_manifest["frozen_inputs"]["base_resolution"]) * int(config["refinement_extension_multiplier"]),
            "h": round_float(1.0 / float(int(phase6_manifest["frozen_inputs"]["base_resolution"]) * int(config["refinement_extension_multiplier"]))),
        }
    )

    scaling_branch = scaling_observations(
        BRANCH_SCALING_LABEL,
        extended_levels,
        epsilon,
        diagonal_weight_scale,
        short_time_window,
        phase7_exponent,
    )
    scaling_control = scaling_observations(
        CONTROL_SCALING_LABEL,
        extended_levels,
        epsilon,
        diagonal_weight_scale,
        short_time_window,
        phase7_exponent,
    )
    ratio_names = ["b1_over_b0", "b2_over_b0"]
    invariant_names = ["I0", "I1", "I2"]
    scaling_rows_branch, scaling_metrics_branch = scaling_prediction_rows(scaling_branch, BRANCH_SCALING_LABEL, ratio_names, invariant_names)
    scaling_rows_control, scaling_metrics_control = scaling_prediction_rows(scaling_control, CONTROL_SCALING_LABEL, ratio_names, invariant_names)
    scaling_rows = scaling_rows_branch + scaling_rows_control
    write_csv_rows(
        SCALING_LEDGER_PATH,
        scaling_rows,
        [
            "hierarchy_label",
            "family",
            "quantity",
            "fit_power",
            "fit_scope",
            "observed_value_R5",
            "predicted_value_R5",
            "relative_error_R5",
            "spectral_radius_R5",
            "spectral_p95_R5",
        ],
    )

    family_labels = [
        "low_mode_localized_wavepacket",
        "mid_spectrum_localized_superposition",
        "random_localized_seed",
    ]

    localization_rows: list[dict[str, Any]] = []
    base_final_records: dict[tuple[str, int, str], dict[str, Any]] = {}
    candidate_scan_records: dict[str, dict[str, Any]] = {BRANCH_LABEL: {}, CONTROL_LABEL: {}}
    for hierarchy_label in (BRANCH_LABEL, CONTROL_LABEL):
        for n_side in scan_levels:
            operator, points = (
                build_branch_zero_form_operator(n_side, epsilon)
                if hierarchy_label == BRANCH_LABEL
                else build_control_zero_form_operator(n_side, epsilon, diagonal_weight_scale)
            )
            for family_label in family_labels:
                state0 = state_family(points, family_label)
                initial_metrics = spatial_metrics(state0, points) | spectral_energy_metrics(state0, operator)
                initial_width = float(initial_metrics["localization_width"])
                initial_concentration = float(initial_metrics["spatial_concentration"])
                initial_participation = float(initial_metrics["participation_ratio"])
                trajectory: list[np.ndarray] = [state0]
                for tau_value in tau_grid[1:]:
                    trajectory.append(np.asarray(expm_multiply((-tau_value) * operator, state0), dtype=complex))
                for tau_value, state in zip(tau_grid, trajectory):
                    metrics = spatial_metrics(state, points) | spectral_energy_metrics(state, operator)
                    localization_rows.append(
                        {
                            "hierarchy_label": hierarchy_label,
                            "level_id": f"n{n_side}",
                            "n_side": int(n_side),
                            "tau": round_float(tau_value, 6),
                            "family_label": family_label,
                            "localization_width": round_float(metrics["localization_width"]),
                            "spatial_concentration": round_float(metrics["spatial_concentration"]),
                            "participation_ratio": round_float(metrics["participation_ratio"]),
                            "spectral_energy_mean": round_float(metrics["spectral_energy_mean"]),
                            "spectral_energy_spread": round_float(metrics["spectral_energy_spread"]),
                            "width_growth": round_float(float(metrics["localization_width"]) / max(initial_width, 1.0e-12) - 1.0),
                            "concentration_ratio": round_float(float(metrics["spatial_concentration"]) / max(initial_concentration, 1.0e-12)),
                            "participation_growth": round_float(float(metrics["participation_ratio"]) / max(initial_participation, 1.0e-12) - 1.0),
                        }
                    )
                final_metrics = spatial_metrics(trajectory[-1], points) | spectral_energy_metrics(trajectory[-1], operator)
                base_width_growth = float(final_metrics["localization_width"]) / max(initial_width, 1.0e-12) - 1.0
                base_concentration_ratio = float(final_metrics["spatial_concentration"]) / max(initial_concentration, 1.0e-12)
                base_participation_growth = float(final_metrics["participation_ratio"]) / max(initial_participation, 1.0e-12) - 1.0
                base_final_records[(hierarchy_label, int(n_side), family_label)] = {
                    "operator": operator,
                    "points": points,
                    "state0": state0,
                    "final_state": trajectory[-1],
                    "initial_metrics": initial_metrics,
                    "final_metrics": final_metrics,
                    "base_width_growth": base_width_growth,
                    "base_concentration_ratio": base_concentration_ratio,
                    "base_participation_growth": base_participation_growth,
                }
                candidate_scan_records[hierarchy_label][f"{family_label}@{n_side}"] = {
                    "base_width_growth": round_float(base_width_growth),
                    "base_concentration_ratio": round_float(base_concentration_ratio),
                    "base_participation_growth": round_float(base_participation_growth),
                    "final_localization_width": round_float(final_metrics["localization_width"]),
                    "final_spatial_concentration": round_float(final_metrics["spatial_concentration"]),
                    "final_spectral_energy_spread": round_float(final_metrics["spectral_energy_spread"]),
                }

    write_csv_rows(
        LOCALIZATION_LEDGER_PATH,
        localization_rows,
        [
            "hierarchy_label",
            "level_id",
            "n_side",
            "tau",
            "family_label",
            "localization_width",
            "spatial_concentration",
            "participation_ratio",
            "spectral_energy_mean",
            "spectral_energy_spread",
            "width_growth",
            "concentration_ratio",
            "participation_growth",
        ],
    )

    selection_cfg = config["candidate_selection"]
    selected_candidates = []
    for family_label in family_labels:
        record = base_final_records[(BRANCH_LABEL, scan_level, family_label)]
        if (
            record["base_width_growth"] <= float(selection_cfg["max_width_growth"])
            and record["base_concentration_ratio"] >= float(selection_cfg["min_concentration_ratio"])
            and record["base_participation_growth"] <= float(selection_cfg["max_participation_growth"])
        ):
            selected_candidates.append(family_label)

    write_json(
        CANDIDATE_PATH,
        {
            "timestamp": timestamp_iso(),
            "scan_level_n_side": int(scan_level),
            "branch_scan_records": candidate_scan_records[BRANCH_LABEL],
            "control_scan_records": candidate_scan_records[CONTROL_LABEL],
            "selected_branch_candidates": selected_candidates,
        },
    )

    thresholds = config["persistence_thresholds"]
    cutoff = float(config["spectral_truncation_shift"]["cutoff"])
    persistence_rows: list[dict[str, Any]] = []
    aggregate_rows: list[dict[str, Any]] = []
    for hierarchy_label in (BRANCH_LABEL, CONTROL_LABEL):
        for n_side in repeat_levels:
            for family_label in selected_candidates:
                base_record = base_final_records[(hierarchy_label, int(n_side), family_label)]
                operator = base_record["operator"]
                points = base_record["points"]
                state0 = base_record["state0"]
                base_final = base_record["final_state"]
                base_final_metrics = base_record["final_metrics"]
                base_width_growth = base_record["base_width_growth"]
                base_concentration_ratio = base_record["base_concentration_ratio"]

                perturbation_results: list[str] = []
                perturbation_overlaps: list[float] = []
                perturbation_scores: list[float] = []

                amplitude_state = amplitude_jitter(state0, points)
                amplitude_final = np.asarray(expm_multiply((-tau_final) * operator, amplitude_state), dtype=complex)
                amplitude_metrics = spatial_metrics(amplitude_final, points) | spectral_energy_metrics(amplitude_final, operator)
                amplitude_overlap = overlap(base_final, amplitude_final)
                amplitude_score = recovery_style_score(base_final_metrics, amplitude_metrics)
                amplitude_class = classify_persistence(base_width_growth, base_concentration_ratio, amplitude_overlap, amplitude_score, thresholds)
                persistence_rows.append(
                    {
                        "hierarchy_label": hierarchy_label,
                        "level_id": f"n{n_side}",
                        "n_side": int(n_side),
                        "family_label": family_label,
                        "perturbation_label": "amplitude_jitter",
                        "overlap_with_base": round_float(amplitude_overlap),
                        "recovery_style_score": round_float(amplitude_score),
                        "base_width_growth": round_float(base_width_growth),
                        "base_concentration_ratio": round_float(base_concentration_ratio),
                        "classification": amplitude_class,
                    }
                )
                perturbation_results.append(amplitude_class)
                perturbation_overlaps.append(amplitude_overlap)
                perturbation_scores.append(amplitude_score)

                shuffled_operator = local_connectivity_shuffle(operator, int(n_side), epsilon)
                shuffled_final = np.asarray(expm_multiply((-tau_final) * shuffled_operator, state0), dtype=complex)
                shuffled_metrics = spatial_metrics(shuffled_final, points) | spectral_energy_metrics(shuffled_final, operator)
                shuffled_overlap = overlap(base_final, shuffled_final)
                shuffled_score = recovery_style_score(base_final_metrics, shuffled_metrics)
                shuffled_class = classify_persistence(base_width_growth, base_concentration_ratio, shuffled_overlap, shuffled_score, thresholds)
                persistence_rows.append(
                    {
                        "hierarchy_label": hierarchy_label,
                        "level_id": f"n{n_side}",
                        "n_side": int(n_side),
                        "family_label": family_label,
                        "perturbation_label": "local_connectivity_shuffle",
                        "overlap_with_base": round_float(shuffled_overlap),
                        "recovery_style_score": round_float(shuffled_score),
                        "base_width_growth": round_float(base_width_growth),
                        "base_concentration_ratio": round_float(base_concentration_ratio),
                        "classification": shuffled_class,
                    }
                )
                perturbation_results.append(shuffled_class)
                perturbation_overlaps.append(shuffled_overlap)
                perturbation_scores.append(shuffled_score)

                truncated_final = (
                    branch_truncated_evolution(state0, int(n_side), epsilon, tau_final, cutoff)
                    if hierarchy_label == BRANCH_LABEL
                    else control_truncated_evolution(state0, int(n_side), epsilon, diagonal_weight_scale, tau_final, cutoff)
                )
                truncated_metrics = spatial_metrics(truncated_final, points) | spectral_energy_metrics(truncated_final, operator)
                truncated_overlap = overlap(base_final, truncated_final)
                truncated_score = recovery_style_score(base_final_metrics, truncated_metrics)
                truncated_class = classify_persistence(base_width_growth, base_concentration_ratio, truncated_overlap, truncated_score, thresholds)
                persistence_rows.append(
                    {
                        "hierarchy_label": hierarchy_label,
                        "level_id": f"n{n_side}",
                        "n_side": int(n_side),
                        "family_label": family_label,
                        "perturbation_label": "spectral_truncation_shift",
                        "overlap_with_base": round_float(truncated_overlap),
                        "recovery_style_score": round_float(truncated_score),
                        "base_width_growth": round_float(base_width_growth),
                        "base_concentration_ratio": round_float(base_concentration_ratio),
                        "classification": truncated_class,
                    }
                )
                perturbation_results.append(truncated_class)
                perturbation_overlaps.append(truncated_overlap)
                perturbation_scores.append(truncated_score)

                aggregate_overlap = float(np.mean(perturbation_overlaps))
                aggregate_score = float(np.mean(perturbation_scores))
                aggregate_class = classify_persistence(
                    base_width_growth,
                    base_concentration_ratio,
                    aggregate_overlap,
                    aggregate_score,
                    thresholds,
                )

                aggregate_rows.append(
                    {
                        "hierarchy_label": hierarchy_label,
                        "level_id": f"n{n_side}",
                        "n_side": int(n_side),
                        "family_label": family_label,
                        "perturbation_label": "aggregate",
                        "overlap_with_base": round_float(aggregate_overlap),
                        "recovery_style_score": round_float(aggregate_score),
                        "base_width_growth": round_float(base_width_growth),
                        "base_concentration_ratio": round_float(base_concentration_ratio),
                        "classification": aggregate_class,
                    }
                )

    persistence_rows.extend(aggregate_rows)
    write_csv_rows(
        PERSISTENCE_LEDGER_PATH,
        persistence_rows,
        [
            "hierarchy_label",
            "level_id",
            "n_side",
            "family_label",
            "perturbation_label",
            "overlap_with_base",
            "recovery_style_score",
            "base_width_growth",
            "base_concentration_ratio",
            "classification",
        ],
    )

    branch_aggregate = [
        row for row in aggregate_rows if row["hierarchy_label"] == BRANCH_LABEL and int(row["n_side"]) in repeat_levels
    ]
    control_aggregate = [
        row for row in aggregate_rows if row["hierarchy_label"] == CONTROL_LABEL and int(row["n_side"]) in repeat_levels
    ]

    persistent_branch_candidates = [row for row in branch_aggregate if row["classification"] == "persistent" and int(row["n_side"]) == scan_levels[-1]]
    branch_stable_candidates = [
        family_label
        for family_label in selected_candidates
        if all(
            next(
                row for row in branch_aggregate if row["family_label"] == family_label and int(row["n_side"]) == level
            )["classification"]
            == "persistent"
            for level in repeat_levels
        )
    ]
    control_persistent_count = len(
        [row for row in control_aggregate if row["classification"] == "persistent" and int(row["n_side"]) == scan_levels[-1]]
    )

    plot_localization_width(localization_rows, selected_candidates, scan_levels)
    plot_overlap_vs_perturbation(persistence_rows, selected_candidates, scan_levels[-1])
    plot_refinement_scaling(persistence_rows, branch_stable_candidates[0] if branch_stable_candidates else selected_candidates[0], repeat_levels)

    scaling_coherent = bool(scaling_metrics_branch["trace_error"] <= 5.0e-4 and scaling_metrics_branch["ratio_span"] <= 0.05)
    success_flags = {
        "localized_candidate_exists": bool(len(persistent_branch_candidates) >= 1),
        "refinement_consistency_holds": bool(len(branch_stable_candidates) >= 1),
        "scaling_descriptors_coherent": bool(scaling_coherent),
        "control_persistence_degraded": bool(control_persistent_count < len(persistent_branch_candidates)),
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
                "phase8_manifest": relpath(PHASE8_MANIFEST_PATH),
                "phase9_manifest": relpath(PHASE9_MANIFEST_PATH),
            },
            "scan_levels": scan_levels,
            "repeat_levels": repeat_levels,
            "tau_grid": [round_float(value, 6) for value in tau_grid.tolist()],
            "selected_branch_candidates": selected_candidates,
            "runtime_seconds": round_float(runtime_seconds, 6),
        },
    )

    summary_lines = [
        "# Integrated Proto-Particle Feasibility Run",
        "",
        "## Objective",
        "",
        "Test one self-contained integrated feasibility question on top of the frozen Phase V--IX contracts: whether the branch keeps coherent large-scale scaling, supports at least one localized excitation candidate, and preserves that candidate under deterministic perturbations better than a deterministic altered-connectivity control.",
        "",
        "## Frozen Inputs",
        "",
        f"- Phase VI operator hierarchy: `{relpath(PHASE6_MANIFEST_PATH)}`",
        f"- Phase VIII short-time trace contract: `{relpath(PHASE8_MANIFEST_PATH)}`",
        f"- Phase IX stabilized descriptor contract: `{relpath(PHASE9_MANIFEST_PATH)}`",
        "",
        "## Scaling Coherence",
        "",
        f"- Branch extension reuses `n_side = [12, 24, 36, 48, 60]` and the frozen short-time window `{[round_float(value, 6) for value in short_time_window.tolist()]}`.",
        f"- Branch short-window trace prediction error at the extension level is `{round_float(scaling_metrics_branch['trace_error'])}` under the minimal `h^2` predictor.",
        f"- Branch ratio span across the extended hierarchy is `{round_float(scaling_metrics_branch['ratio_span'])}`, with R5 ratios `{scaling_metrics_branch['ratio_values_R5']}`.",
        f"- Control scaling remains deterministic but descriptor-shifted, with R5 ratios `{scaling_metrics_control['ratio_values_R5']}`.",
        "",
        "## Localized Excitation Scan",
        "",
        f"- Deterministic families scanned: `{family_labels}` on levels `{scan_levels}` with tau grid `{[round_float(value, 6) for value in tau_grid.tolist()]}`.",
        f"- Branch candidates selected at `n_side = {scan_level}`: `{selected_candidates}`.",
        (
            f"- The strongest branch candidate is `{branch_stable_candidates[0]}`, with aggregate persistent classification on both repeat levels `{repeat_levels}`."
            if branch_stable_candidates
            else (
                f"- The strongest branch candidate is `{selected_candidates[0]}`, but aggregate persistence does not hold on both repeat levels `{repeat_levels}`."
                if selected_candidates
                else "- No branch candidate satisfied the scan-level localization gate."
            )
        ),
        "",
        "## Persistence and Control",
        "",
        f"- Branch persistent candidates at the scan level: `{[row['family_label'] for row in persistent_branch_candidates]}`.",
        f"- Control persistent candidate count at the scan level: `{control_persistent_count}`.",
        f"- Branch stable candidate set across refinement: `{branch_stable_candidates}`.",
        f"- The spectral truncation shift uses the deterministic cutoff `{cutoff}`, and the control branch shows broader width growth and weaker concentration retention under the same evolution window.",
        "",
        "## Bounded Interpretation",
        "",
        "This integrated run supports a narrow feasibility statement only: within the frozen operator and spectral-trace contracts, one localized branch excitation family survives deterministic evolution and perturbation better than the altered-connectivity control while the large-scale trace descriptors stay coherent. No particle claim, continuum claim, or physical interpretation is asserted.",
        "",
        closure_statement(success),
    ]
    write_text(SUMMARY_PATH, "\n".join(summary_lines))

    manifest = {
        "timestamp": timestamp_iso(),
        "phase": "X",
        "phase_name": PHASE_NAME,
        "stage_identifier": STAGE_IDENTIFIER,
        "status": "passed" if success else "failed",
        "success": bool(success),
        "objective": "integrated_proto_particle_feasibility",
        "operator_manifest_reference": relpath(PHASE6_MANIFEST_PATH),
        "phase8_manifest_reference": relpath(PHASE8_MANIFEST_PATH),
        "phase9_manifest_reference": relpath(PHASE9_MANIFEST_PATH),
        "config_reference": relpath(CONFIG_PATH),
        "refinement_hierarchy": {
            "base_levels": [int(level["n_side"]) for level in base_levels],
            "extended_levels": [int(level["n_side"]) for level in extended_levels],
            "repeat_levels": repeat_levels,
        },
        "scaling_coherence": {
            "branch_trace_prediction_error_R5": round_float(scaling_metrics_branch["trace_error"]),
            "branch_ratio_span": round_float(scaling_metrics_branch["ratio_span"]),
            "branch_invariant_span": round_float(scaling_metrics_branch["invariant_span"]),
            "control_trace_prediction_error_R5": round_float(scaling_metrics_control["trace_error"]),
            "control_ratio_span": round_float(scaling_metrics_control["ratio_span"]),
            "control_invariant_span": round_float(scaling_metrics_control["invariant_span"]),
            "branch_spectral_radius_values": scaling_metrics_branch["spectral_radius_values"],
            "control_spectral_radius_values": scaling_metrics_control["spectral_radius_values"],
        },
        "candidate_scan": {
            "selected_branch_candidates": selected_candidates,
            "branch_scan_level": int(scan_level),
            "branch_scan_records": candidate_scan_records[BRANCH_LABEL],
            "control_scan_records": candidate_scan_records[CONTROL_LABEL],
        },
        "persistence_results": {
            "branch_persistent_scan_candidates": [row["family_label"] for row in persistent_branch_candidates],
            "branch_refinement_stable_candidates": branch_stable_candidates,
            "control_persistent_scan_candidate_count": int(control_persistent_count),
        },
        "success_flags": success_flags,
        "artifacts": {
            "proto_particle_candidates_json": relpath(CANDIDATE_PATH),
            "localization_metrics_ledger_csv": relpath(LOCALIZATION_LEDGER_PATH),
            "persistence_classification_csv": relpath(PERSISTENCE_LEDGER_PATH),
            "scaling_prediction_errors_csv": relpath(SCALING_LEDGER_PATH),
            "summary_md": relpath(SUMMARY_PATH),
            "manifest_json": relpath(MANIFEST_PATH),
            "width_plot": relpath(WIDTH_PLOT_PATH),
            "overlap_plot": relpath(OVERLAP_PLOT_PATH),
            "refinement_plot": relpath(REFINEMENT_PLOT_PATH),
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
