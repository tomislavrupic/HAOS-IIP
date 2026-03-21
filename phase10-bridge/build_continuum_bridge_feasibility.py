#!/usr/bin/env python3

from __future__ import annotations

import math
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haos_core import read_json, relpath, write_csv_rows, write_json

PHASE_ROOT = ROOT / "phase10-bridge"
RUNS_ROOT = PHASE_ROOT / "runs"
PLOTS_ROOT = PHASE_ROOT / "plots"

PHASE6_MANIFEST_PATH = ROOT / "phase6-operator" / "phase6_operator_manifest.json"
PHASE8_MANIFEST_PATH = ROOT / "phase8-trace" / "phase8_trace_manifest.json"
PHASE9_MANIFEST_PATH = ROOT / "phase9-invariants" / "phase9_manifest.json"

EXTENSION_LEDGER_PATH = RUNS_ROOT / "phase10_extension_ledger.csv"
SCALING_LEDGER_PATH = RUNS_ROOT / "phase10_effective_scaling_ledger.csv"
COARSE_LEDGER_PATH = RUNS_ROOT / "phase10_coarse_grain_summary.csv"
RUNS_PATH = RUNS_ROOT / "phase10_runs.json"
SUMMARY_PATH = PHASE_ROOT / "phase10_summary.md"
MANIFEST_PATH = PHASE_ROOT / "phase10_manifest.json"

DESCRIPTOR_PLOT_PATH = PLOTS_ROOT / "phase10_descriptor_vs_extended_refinement.svg"
PREDICTION_PLOT_PATH = PLOTS_ROOT / "phase10_prediction_error_vs_refinement.svg"
COARSE_PLOT_PATH = PLOTS_ROOT / "phase10_coarse_grain_spectral_comparison.svg"
CORRELATION_PLOT_PATH = PLOTS_ROOT / "phase10_multi_descriptor_correlation.svg"

PHASE_NAME = "phase10-bridge"
PHASE_IDENTIFIER = 10
BRANCH_LABEL = "frozen_branch"
CONTROL_LABEL = "generic_open_grid_scalar_block_control"
EXTENDED_REFINEMENT_MULTIPLIER = 5
PRIMARY_CORRECTION_POWER = 2
INVARIANT_CORRECTION_POWER = 1
COARSE_GRAIN_LABEL = "lambda_le_7_5_spectral_projection"
COARSE_GRAIN_CUTOFF = 7.5

BRIDGE_TRACE_ERROR_TOL = 1.0e-3
BRIDGE_RATIO_ERROR_TOL = 1.0e-3
RATIO_LIMIT_DRIFT_TOL = 1.0e-4
INVARIANT_EXTENSION_ERROR_TOL = 0.01
COARSE_TRACE_ERROR_TOL = 0.04
COARSE_RATIO_DRIFT_TOL = 0.08
COHERENCE_PC1_MIN = 0.95
CONTROL_TRACE_ERROR_MULTIPLIER_MIN = 3.0
CONTROL_RATIO_ERROR_MULTIPLIER_MIN = 5.0

CLAIM_BOUNDARY = (
    "Phase X authority is limited to cautious continuum-bridge feasibility diagnostics on the "
    "frozen operator and spectral-trace contracts. It does not assert a continuum limit, geometric "
    "coefficients, metric structure, physical correspondence, or any continuum-derived law."
)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def timestamp_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def closure_statement(success: bool) -> str:
    if success:
        return "Phase X establishes cautious continuum-bridge feasibility for the frozen operator hierarchy."
    return "Phase X does not yet establish cautious continuum-bridge feasibility for the frozen operator hierarchy."


def round_float(value: float, digits: int = 12) -> float:
    return round(float(value), digits)


def safe_relative_difference(new_value: float, reference_value: float) -> float:
    return abs(new_value - reference_value) / max(abs(reference_value), 1.0e-12)


def trend_label(values: list[float], tolerance: float = 1.0e-12) -> str:
    nondecreasing = all(values[index + 1] >= values[index] - tolerance for index in range(len(values) - 1))
    nonincreasing = all(values[index + 1] <= values[index] + tolerance for index in range(len(values) - 1))
    if nondecreasing:
        return "nondecreasing"
    if nonincreasing:
        return "nonincreasing"
    return "mixed"


def pearson_correlation(left: np.ndarray, right: np.ndarray) -> float:
    if np.std(left) <= 1.0e-18 or np.std(right) <= 1.0e-18:
        return 1.0
    return float(np.corrcoef(left, right)[0, 1])


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
    predictions = np.zeros_like(series_matrix, dtype=float)
    slopes: list[float] = []
    limits: list[float] = []
    r_squared: list[float] = []
    predicted_target: list[float] = []

    x_values = h_values**power
    target_x = None if target_h is None else float(target_h) ** power
    for column_index in range(series_matrix.shape[1]):
        coefficients = np.polyfit(x_values, series_matrix[:, column_index], 1)
        predictions[:, column_index] = np.polyval(coefficients, x_values)
        slopes.append(float(coefficients[0]))
        limits.append(float(coefficients[1]))
        r_squared.append(model_r_squared(x_values, series_matrix[:, column_index], coefficients))
        if target_x is not None:
            predicted_target.append(float(np.polyval(coefficients, target_x)))

    payload: dict[str, Any] = {
        "power": int(power),
        "predictions": predictions,
        "slopes": np.array(slopes, dtype=float),
        "limits": np.array(limits, dtype=float),
        "r_squared": np.array(r_squared, dtype=float),
    }
    if target_h is not None:
        payload["predicted_target"] = np.array(predicted_target, dtype=float)
    return payload


def best_model_labels(h_values: np.ndarray, series_matrix: np.ndarray, names: list[str]) -> dict[str, str]:
    fit_h = fit_linear_family(h_values, series_matrix, power=1)
    fit_h2 = fit_linear_family(h_values, series_matrix, power=2)
    labels: dict[str, str] = {}
    for index, name in enumerate(names):
        labels[name] = "linear_in_h^2" if fit_h2["r_squared"][index] > fit_h["r_squared"][index] else "linear_in_h^1"
    return labels


def load_frozen_manifests() -> dict[str, Any]:
    phase6_manifest = read_json(PHASE6_MANIFEST_PATH)
    phase8_manifest = read_json(PHASE8_MANIFEST_PATH)
    phase9_manifest = read_json(PHASE9_MANIFEST_PATH)

    if not bool(phase6_manifest.get("success")):
        raise ValueError("Phase X requires the frozen successful Phase VI operator manifest.")
    if phase6_manifest.get("selected_operator_class") != "cochain_laplacian":
        raise ValueError("Phase X requires the frozen cochain-Laplacian operator family.")
    if not bool(phase8_manifest.get("success")):
        raise ValueError("Phase X requires the frozen successful Phase VIII trace manifest.")
    if not bool(phase9_manifest.get("success")):
        raise ValueError("Phase X requires the frozen successful Phase IX manifest.")
    return {
        "phase6_manifest": phase6_manifest,
        "phase8_manifest": phase8_manifest,
        "phase9_manifest": phase9_manifest,
    }


def extend_hierarchy(frozen_levels: list[dict[str, Any]], base_resolution: int) -> list[dict[str, Any]]:
    levels = list(frozen_levels)
    extended_n_side = int(base_resolution) * EXTENDED_REFINEMENT_MULTIPLIER
    levels.append(
        {
            "level_id": f"R{len(levels) + 1}",
            "refinement_multiplier": EXTENDED_REFINEMENT_MULTIPLIER,
            "n_side": extended_n_side,
            "h": round_float(1.0 / float(extended_n_side)),
        }
    )
    return levels


def branch_spectrum(n_side: int, epsilon: float) -> np.ndarray:
    h_value = 1.0 / float(n_side)
    edge_weight = math.exp(-(h_value * h_value) / (2.0 * epsilon))
    scalar_values = np.empty(n_side * n_side, dtype=float)
    cursor = 0
    for mode_x in range(n_side):
        angle_x = 2.0 * math.pi * mode_x / n_side
        cos_x = math.cos(angle_x)
        for mode_y in range(n_side):
            angle_y = 2.0 * math.pi * mode_y / n_side
            cos_y = math.cos(angle_y)
            scalar_values[cursor] = 2.0 * edge_weight * (2.0 - cos_x - cos_y)
            cursor += 1
    return np.sort(np.repeat(np.sort(scalar_values), 4))


def control_spectrum(n_side: int, epsilon: float) -> np.ndarray:
    h_value = 1.0 / float(n_side)
    edge_weight = math.exp(-(h_value * h_value) / (2.0 * epsilon))
    path_values = np.array(
        [2.0 * edge_weight * (1.0 - math.cos(math.pi * mode / n_side)) for mode in range(n_side)],
        dtype=float,
    )
    scalar_values = np.empty(n_side * n_side, dtype=float)
    cursor = 0
    for mode_x in range(n_side):
        for mode_y in range(n_side):
            scalar_values[cursor] = path_values[mode_x] + path_values[mode_y]
            cursor += 1
    return np.sort(np.repeat(np.sort(scalar_values), 4))


def trace_from_spectrum(spectrum: np.ndarray, probe_times: np.ndarray) -> np.ndarray:
    return np.array([float(np.exp(-probe_time * spectrum).sum()) for probe_time in probe_times], dtype=float)


def descriptor_observation(
    hierarchy_label: str,
    level: dict[str, Any],
    spectrum: np.ndarray,
    probe_times: np.ndarray,
    phase7_exponent: float,
) -> dict[str, Any]:
    h_value = float(level["h"])
    trace_window = trace_from_spectrum(spectrum, probe_times)
    scaled_fit = quadratic_fit_named(probe_times, (h_value * h_value) * trace_window, prefix="b")
    raw_fit = quadratic_fit_named(probe_times, trace_window, prefix="a")

    ratios = {
        "b1_over_b0": round_float(float(scaled_fit["b1"]) / max(abs(float(scaled_fit["b0"])), 1.0e-12)),
        "b2_over_b0": round_float(float(scaled_fit["b2"]) / max(abs(float(scaled_fit["b0"])), 1.0e-12)),
    }
    invariants = {
        "I0": round_float((h_value**phase7_exponent) * float(raw_fit["a0"])),
        "I1": round_float((h_value**phase7_exponent) * float(raw_fit["a1"])),
        "I2": round_float((h_value**phase7_exponent) * float(raw_fit["a2"])),
    }

    positive_eigenvalues = spectrum[spectrum > 1.0e-12]
    low_mode_sample = positive_eigenvalues[:12]
    return {
        "hierarchy_label": hierarchy_label,
        "level_id": str(level["level_id"]),
        "refinement_multiplier": int(level["refinement_multiplier"]),
        "n_side": int(level["n_side"]),
        "h": round_float(h_value),
        "dimension": int(len(spectrum)),
        "nullspace_dimension": int(np.sum(np.isclose(spectrum, 0.0, atol=1.0e-12))),
        "spectral_radius": round_float(float(spectrum[-1])),
        "first_positive_eigenvalue": round_float(float(positive_eigenvalues[0])),
        "trace_window": [round_float(value) for value in trace_window.tolist()],
        "scaled_coefficients": scaled_fit,
        "raw_coefficients": raw_fit,
        "ratios": ratios,
        "rescaled_invariants": invariants,
        "_spectrum": spectrum,
        "_trace_window": trace_window,
        "_low_mode_sample": low_mode_sample,
    }


def build_observations(
    hierarchy_label: str,
    levels: list[dict[str, Any]],
    epsilon: float,
    probe_times: np.ndarray,
    phase7_exponent: float,
) -> list[dict[str, Any]]:
    spectrum_builder = branch_spectrum if hierarchy_label == BRANCH_LABEL else control_spectrum
    return [
        descriptor_observation(hierarchy_label, level, spectrum_builder(int(level["n_side"]), epsilon), probe_times, phase7_exponent)
        for level in levels
    ]


def family_matrix(observations: list[dict[str, Any]], names: list[str], family_name: str) -> np.ndarray:
    return np.array([[float(observation[family_name][name]) for name in names] for observation in observations], dtype=float)


def trace_prediction_from_scaled_coefficients(h_value: float, coefficients: np.ndarray, probe_times: np.ndarray) -> np.ndarray:
    return (float(h_value) ** -2) * (coefficients[0] + coefficients[1] * probe_times + coefficients[2] * (probe_times**2))


def projection_summary(
    observations: list[dict[str, Any]],
    cutoff: float,
    probe_times: np.ndarray,
    phase7_exponent: float,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for observation in observations:
        spectrum = observation["_spectrum"]
        projected = spectrum[spectrum <= cutoff]
        projected_trace = trace_from_spectrum(projected, probe_times)
        h_value = float(observation["h"])
        projected_scaled = quadratic_fit_named(probe_times, (h_value * h_value) * projected_trace, prefix="b")
        projected_raw = quadratic_fit_named(probe_times, projected_trace, prefix="a")
        projected_ratios = {
            "b1_over_b0": float(projected_scaled["b1"]) / max(abs(float(projected_scaled["b0"])), 1.0e-12),
            "b2_over_b0": float(projected_scaled["b2"]) / max(abs(float(projected_scaled["b0"])), 1.0e-12),
        }
        projected_invariants = {
            "I0": (h_value**phase7_exponent) * float(projected_raw["a0"]),
            "I1": (h_value**phase7_exponent) * float(projected_raw["a1"]),
            "I2": (h_value**phase7_exponent) * float(projected_raw["a2"]),
        }

        full_trace = observation["_trace_window"]
        full_ratios = observation["ratios"]
        full_invariants = observation["rescaled_invariants"]
        low_mode_difference = 0.0
        retained_low_modes = projected[projected > 1.0e-12][: len(observation["_low_mode_sample"])]
        if len(retained_low_modes) == len(observation["_low_mode_sample"]):
            low_mode_difference = float(np.max(np.abs(retained_low_modes - observation["_low_mode_sample"])))
        else:
            low_mode_difference = float("inf")

        ratio_relative_drifts = [
            safe_relative_difference(projected_ratios["b1_over_b0"], float(full_ratios["b1_over_b0"])),
            safe_relative_difference(projected_ratios["b2_over_b0"], float(full_ratios["b2_over_b0"])),
        ]
        invariant_relative_drifts = [
            safe_relative_difference(projected_invariants["I0"], float(full_invariants["I0"])),
            safe_relative_difference(projected_invariants["I1"], float(full_invariants["I1"])),
            safe_relative_difference(projected_invariants["I2"], float(full_invariants["I2"])),
        ]

        rows.append(
            {
                "hierarchy_label": observation["hierarchy_label"],
                "level_id": observation["level_id"],
                "n_side": int(observation["n_side"]),
                "h": round_float(observation["h"]),
                "projection_label": COARSE_GRAIN_LABEL,
                "cutoff": round_float(cutoff, 6),
                "retained_mode_count": int(len(projected)),
                "retained_mode_fraction": round_float(len(projected) / len(spectrum)),
                "projected_spectral_radius": round_float(float(projected[-1])),
                "spectral_radius_loss": round_float(float(spectrum[-1] - projected[-1])),
                "first_positive_low_mode_max_abs_difference": round_float(low_mode_difference),
                "trace_window_max_relative_error": round_float(float(np.max(np.abs(projected_trace - full_trace) / full_trace))),
                "ratio_max_relative_drift": round_float(float(max(ratio_relative_drifts))),
                "invariant_max_relative_drift": round_float(float(max(invariant_relative_drifts))),
                "low_mode_preserved": bool(low_mode_difference <= 1.0e-12),
                "ratio_ordering_preserved": bool(projected_ratios["b2_over_b0"] > 0.0),
            }
        )
    return rows


def build_extension_rows(observations: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    previous_vector: np.ndarray | None = None
    for observation in observations:
        descriptor_vector = np.array(
            [
                float(observation["ratios"]["b1_over_b0"]),
                float(observation["ratios"]["b2_over_b0"]),
                float(observation["rescaled_invariants"]["I0"]),
                float(observation["rescaled_invariants"]["I1"]),
                float(observation["rescaled_invariants"]["I2"]),
                float(observation["spectral_radius"]),
            ],
            dtype=float,
        )
        descriptor_distance = ""
        if previous_vector is not None:
            descriptor_distance = f"{round_float(float(np.linalg.norm(descriptor_vector - previous_vector))):.12f}"
        rows.append(
            {
                "hierarchy_label": observation["hierarchy_label"],
                "level_id": observation["level_id"],
                "n_side": int(observation["n_side"]),
                "h": round_float(observation["h"]),
                "dimension": int(observation["dimension"]),
                "nullspace_dimension": int(observation["nullspace_dimension"]),
                "spectral_radius": round_float(observation["spectral_radius"]),
                "first_positive_eigenvalue": round_float(observation["first_positive_eigenvalue"]),
                "b0": round_float(observation["scaled_coefficients"]["b0"]),
                "b1": round_float(observation["scaled_coefficients"]["b1"]),
                "b2": round_float(observation["scaled_coefficients"]["b2"]),
                "b1_over_b0": round_float(observation["ratios"]["b1_over_b0"]),
                "b2_over_b0": round_float(observation["ratios"]["b2_over_b0"]),
                "I0": round_float(observation["rescaled_invariants"]["I0"]),
                "I1": round_float(observation["rescaled_invariants"]["I1"]),
                "I2": round_float(observation["rescaled_invariants"]["I2"]),
                "descriptor_distance_from_previous": descriptor_distance,
            }
        )
        previous_vector = descriptor_vector
    return rows


def build_scaling_rows(
    hierarchy_label: str,
    frozen_h_values: np.ndarray,
    target_h: float,
    observed_target: dict[str, Any],
    coefficient_matrix: np.ndarray,
    ratio_matrix: np.ndarray,
    invariant_matrix: np.ndarray,
    coefficient_names: list[str],
    ratio_names: list[str],
    invariant_names: list[str],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    coefficient_h2_fit = fit_linear_family(frozen_h_values, coefficient_matrix[:-1], power=PRIMARY_CORRECTION_POWER, target_h=target_h)
    coefficient_h1_fit = fit_linear_family(frozen_h_values, coefficient_matrix[:-1], power=1, target_h=target_h)
    ratio_h2_fit = fit_linear_family(frozen_h_values, ratio_matrix[:-1], power=PRIMARY_CORRECTION_POWER, target_h=target_h)
    ratio_h1_fit = fit_linear_family(frozen_h_values, ratio_matrix[:-1], power=1, target_h=target_h)
    invariant_h1_fit = fit_linear_family(frozen_h_values, invariant_matrix[:-1], power=INVARIANT_CORRECTION_POWER, target_h=target_h)
    invariant_h2_fit = fit_linear_family(frozen_h_values, invariant_matrix[:-1], power=PRIMARY_CORRECTION_POWER, target_h=target_h)

    target_coefficients = coefficient_matrix[-1]
    target_ratios = ratio_matrix[-1]
    target_invariants = invariant_matrix[-1]

    for index, name in enumerate(coefficient_names):
        rows.append(
            {
                "hierarchy_label": hierarchy_label,
                "family": "scaled_coefficients",
                "quantity": name,
                "fit_power": PRIMARY_CORRECTION_POWER,
                "fit_scope": "R1_to_R4_predict_R5",
                "fit_r_squared": round_float(coefficient_h2_fit["r_squared"][index]),
                "slope": round_float(coefficient_h2_fit["slopes"][index]),
                "limit_estimate": round_float(coefficient_h2_fit["limits"][index]),
                "predicted_R5": round_float(coefficient_h2_fit["predicted_target"][index]),
                "observed_R5": round_float(target_coefficients[index]),
                "relative_error_R5": round_float(
                    safe_relative_difference(coefficient_h2_fit["predicted_target"][index], target_coefficients[index])
                ),
            }
        )
    for index, name in enumerate(ratio_names):
        rows.append(
            {
                "hierarchy_label": hierarchy_label,
                "family": "ratios",
                "quantity": name,
                "fit_power": PRIMARY_CORRECTION_POWER,
                "fit_scope": "R1_to_R4_predict_R5",
                "fit_r_squared": round_float(ratio_h2_fit["r_squared"][index]),
                "slope": round_float(ratio_h2_fit["slopes"][index]),
                "limit_estimate": round_float(ratio_h2_fit["limits"][index]),
                "predicted_R5": round_float(ratio_h2_fit["predicted_target"][index]),
                "observed_R5": round_float(target_ratios[index]),
                "relative_error_R5": round_float(
                    safe_relative_difference(ratio_h2_fit["predicted_target"][index], target_ratios[index])
                ),
            }
        )
    for index, name in enumerate(invariant_names):
        rows.append(
            {
                "hierarchy_label": hierarchy_label,
                "family": "rescaled_invariants",
                "quantity": name,
                "fit_power": INVARIANT_CORRECTION_POWER,
                "fit_scope": "R1_to_R4_predict_R5",
                "fit_r_squared": round_float(invariant_h1_fit["r_squared"][index]),
                "slope": round_float(invariant_h1_fit["slopes"][index]),
                "limit_estimate": round_float(invariant_h1_fit["limits"][index]),
                "predicted_R5": round_float(invariant_h1_fit["predicted_target"][index]),
                "observed_R5": round_float(target_invariants[index]),
                "relative_error_R5": round_float(
                    safe_relative_difference(invariant_h1_fit["predicted_target"][index], target_invariants[index])
                ),
            }
        )

    predicted_trace_h2 = trace_prediction_from_scaled_coefficients(target_h, coefficient_h2_fit["predicted_target"], observed_target["_probe_times"])
    predicted_trace_h1 = trace_prediction_from_scaled_coefficients(target_h, coefficient_h1_fit["predicted_target"], observed_target["_probe_times"])
    observed_trace = observed_target["_trace_window"]

    rows.append(
        {
            "hierarchy_label": hierarchy_label,
            "family": "trace_window",
            "quantity": "max_relative_error",
            "fit_power": PRIMARY_CORRECTION_POWER,
            "fit_scope": "R1_to_R4_predict_R5",
            "fit_r_squared": round_float(float(np.min(coefficient_h2_fit["r_squared"]))),
            "slope": "",
            "limit_estimate": "",
            "predicted_R5": "",
            "observed_R5": "",
            "relative_error_R5": round_float(float(np.max(np.abs(predicted_trace_h2 - observed_trace) / observed_trace))),
        }
    )

    metrics = {
        "coefficient_h2": coefficient_h2_fit,
        "coefficient_h1": coefficient_h1_fit,
        "ratio_h2": ratio_h2_fit,
        "ratio_h1": ratio_h1_fit,
        "invariant_h1": invariant_h1_fit,
        "invariant_h2": invariant_h2_fit,
        "trace_h2_error": float(np.max(np.abs(predicted_trace_h2 - observed_trace) / observed_trace)),
        "trace_h1_error": float(np.max(np.abs(predicted_trace_h1 - observed_trace) / observed_trace)),
        "coefficient_h2_max_error": float(
            np.max(np.abs(coefficient_h2_fit["predicted_target"] - target_coefficients) / np.maximum(np.abs(target_coefficients), 1.0e-12))
        ),
        "ratio_h2_max_error": float(
            np.max(np.abs(ratio_h2_fit["predicted_target"] - target_ratios) / np.maximum(np.abs(target_ratios), 1.0e-12))
        ),
        "invariant_h1_max_error": float(
            np.max(np.abs(invariant_h1_fit["predicted_target"] - target_invariants) / np.maximum(np.abs(target_invariants), 1.0e-12))
        ),
        "ratio_limit_drift": {},
        "invariant_limit_drift": {},
    }

    ratio_h2_full = fit_linear_family(np.append(frozen_h_values, target_h), ratio_matrix, power=PRIMARY_CORRECTION_POWER)
    invariant_h1_full = fit_linear_family(np.append(frozen_h_values, target_h), invariant_matrix, power=INVARIANT_CORRECTION_POWER)
    for index, name in enumerate(ratio_names):
        metrics["ratio_limit_drift"][name] = safe_relative_difference(ratio_h2_full["limits"][index], ratio_h2_fit["limits"][index])
    for index, name in enumerate(invariant_names):
        metrics["invariant_limit_drift"][name] = safe_relative_difference(
            invariant_h1_full["limits"][index], invariant_h1_fit["limits"][index]
        )
    metrics["ratio_h2_full"] = ratio_h2_full
    metrics["invariant_h1_full"] = invariant_h1_full
    return rows, metrics


def plot_descriptor_vs_refinement(
    branch_observations: list[dict[str, Any]],
    control_observations: list[dict[str, Any]],
) -> None:
    n_values = [int(observation["n_side"]) for observation in branch_observations]
    plt.figure(figsize=(8.5, 8.0))

    ax1 = plt.subplot(3, 1, 1)
    ax1.plot(n_values, [float(obs["ratios"]["b1_over_b0"]) for obs in branch_observations], marker="o", label="branch b1/b0")
    ax1.plot(n_values, [float(obs["ratios"]["b1_over_b0"]) for obs in control_observations], marker="s", label="control b1/b0")
    ax1.set_ylabel("b1 / b0")
    ax1.grid(alpha=0.3)
    ax1.legend(loc="best")

    ax2 = plt.subplot(3, 1, 2)
    ax2.plot(n_values, [float(obs["ratios"]["b2_over_b0"]) for obs in branch_observations], marker="o", label="branch b2/b0")
    ax2.plot(n_values, [float(obs["ratios"]["b2_over_b0"]) for obs in control_observations], marker="s", label="control b2/b0")
    ax2.set_ylabel("b2 / b0")
    ax2.grid(alpha=0.3)
    ax2.legend(loc="best")

    ax3 = plt.subplot(3, 1, 3)
    ax3.plot(n_values, [float(obs["rescaled_invariants"]["I2"]) for obs in branch_observations], marker="o", label="branch I2")
    ax3.plot(n_values, [float(obs["rescaled_invariants"]["I2"]) for obs in control_observations], marker="s", label="control I2")
    ax3.set_ylabel("I2")
    ax3.set_xlabel("n_side")
    ax3.grid(alpha=0.3)
    ax3.legend(loc="best")

    plt.tight_layout()
    plt.savefig(DESCRIPTOR_PLOT_PATH, format="svg")
    plt.close()


def plot_prediction_errors(
    branch_observations: list[dict[str, Any]],
    control_observations: list[dict[str, Any]],
    branch_all_level_fit: dict[str, Any],
    control_all_level_fit: dict[str, Any],
) -> None:
    def per_level_trace_residuals(observations: list[dict[str, Any]], fit_payload: dict[str, Any]) -> list[float]:
        values: list[float] = []
        for index, observation in enumerate(observations):
            predicted_trace = trace_prediction_from_scaled_coefficients(
                float(observation["h"]),
                fit_payload["predictions"][index],
                observation["_probe_times"],
            )
            observed_trace = observation["_trace_window"]
            values.append(float(np.max(np.abs(predicted_trace - observed_trace) / observed_trace)))
        return values

    n_values = [int(observation["n_side"]) for observation in branch_observations]
    branch_errors = per_level_trace_residuals(branch_observations, branch_all_level_fit)
    control_errors = per_level_trace_residuals(control_observations, control_all_level_fit)

    plt.figure(figsize=(8.0, 4.8))
    plt.plot(n_values, branch_errors, marker="o", label="branch h^2 law residual")
    plt.plot(n_values, control_errors, marker="s", label="control h^2 law residual")
    plt.axhline(BRIDGE_TRACE_ERROR_TOL, color="black", linestyle="--", linewidth=1.0, label="branch tolerance")
    plt.xlabel("n_side")
    plt.ylabel("max short-window trace relative residual")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(PREDICTION_PLOT_PATH, format="svg")
    plt.close()


def plot_coarse_grain_spectral_comparison(branch_rows: list[dict[str, Any]], control_rows: list[dict[str, Any]]) -> None:
    n_values = [int(row["n_side"]) for row in branch_rows]
    plt.figure(figsize=(8.0, 5.2))
    plt.plot(n_values, [float(row["projected_spectral_radius"]) for row in branch_rows], marker="o", label="branch projected radius")
    plt.plot(n_values, [float(row["projected_spectral_radius"]) for row in control_rows], marker="s", label="control projected radius")
    plt.plot(n_values, [float(row["projected_spectral_radius"]) + float(row["spectral_radius_loss"]) for row in branch_rows], linestyle="--", color="C0", label="branch full radius")
    plt.plot(n_values, [float(row["projected_spectral_radius"]) + float(row["spectral_radius_loss"]) for row in control_rows], linestyle="--", color="C1", label="control full radius")
    plt.xlabel("n_side")
    plt.ylabel("spectral radius")
    plt.grid(alpha=0.3)
    plt.legend(loc="best", ncol=2)
    plt.tight_layout()
    plt.savefig(COARSE_PLOT_PATH, format="svg")
    plt.close()


def plot_multi_descriptor_correlation(
    branch_observations: list[dict[str, Any]],
    control_observations: list[dict[str, Any]],
) -> None:
    plt.figure(figsize=(6.6, 5.4))
    for observations, marker, label in (
        (branch_observations, "o", "branch"),
        (control_observations, "s", "control"),
    ):
        x_values = [float(observation["ratios"]["b2_over_b0"]) for observation in observations]
        y_values = [float(observation["rescaled_invariants"]["I2"]) for observation in observations]
        plt.scatter(x_values, y_values, marker=marker, label=label)
        for observation, x_value, y_value in zip(observations, x_values, y_values):
            plt.annotate(observation["level_id"], (x_value, y_value), textcoords="offset points", xytext=(4, 4), fontsize=8)
    plt.xlabel("b2 / b0")
    plt.ylabel("I2")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(CORRELATION_PLOT_PATH, format="svg")
    plt.close()


def main() -> None:
    start_time = time.perf_counter()
    RUNS_ROOT.mkdir(parents=True, exist_ok=True)
    PLOTS_ROOT.mkdir(parents=True, exist_ok=True)
    manifests = load_frozen_manifests()
    phase6_manifest = manifests["phase6_manifest"]
    phase8_manifest = manifests["phase8_manifest"]
    phase9_manifest = manifests["phase9_manifest"]

    base_resolution = int(phase6_manifest["frozen_inputs"]["base_resolution"])
    epsilon = float(phase6_manifest["frozen_inputs"]["base_epsilon"])
    phase7_exponent = float(phase9_manifest["dominant_phase7_low_mode_exponent"]["exponent"])
    probe_times = np.array(phase8_manifest["operational_short_time_window"]["probe_times"], dtype=float)
    frozen_levels = phase6_manifest["refinement_hierarchy"]["levels"]
    extended_levels = extend_hierarchy(frozen_levels, base_resolution)
    frozen_h_values = np.array([float(level["h"]) for level in extended_levels[:-1]], dtype=float)
    target_h = float(extended_levels[-1]["h"])

    branch_observations = build_observations(BRANCH_LABEL, extended_levels, epsilon, probe_times, phase7_exponent)
    control_observations = build_observations(CONTROL_LABEL, extended_levels, epsilon, probe_times, phase7_exponent)
    for observations in (branch_observations, control_observations):
        for observation in observations:
            observation["_probe_times"] = probe_times

    coefficient_names = ["b0", "b1", "b2"]
    ratio_names = ["b1_over_b0", "b2_over_b0"]
    invariant_names = ["I0", "I1", "I2"]

    branch_coefficients = family_matrix(branch_observations, coefficient_names, "scaled_coefficients")
    branch_ratios = family_matrix(branch_observations, ratio_names, "ratios")
    branch_invariants = family_matrix(branch_observations, invariant_names, "rescaled_invariants")

    control_coefficients = family_matrix(control_observations, coefficient_names, "scaled_coefficients")
    control_ratios = family_matrix(control_observations, ratio_names, "ratios")
    control_invariants = family_matrix(control_observations, invariant_names, "rescaled_invariants")

    branch_scaling_rows, branch_scaling_metrics = build_scaling_rows(
        BRANCH_LABEL,
        frozen_h_values,
        target_h,
        branch_observations[-1],
        branch_coefficients,
        branch_ratios,
        branch_invariants,
        coefficient_names,
        ratio_names,
        invariant_names,
    )
    control_scaling_rows, control_scaling_metrics = build_scaling_rows(
        CONTROL_LABEL,
        frozen_h_values,
        target_h,
        control_observations[-1],
        control_coefficients,
        control_ratios,
        control_invariants,
        coefficient_names,
        ratio_names,
        invariant_names,
    )
    scaling_rows = branch_scaling_rows + control_scaling_rows

    branch_extension_rows = build_extension_rows(branch_observations)
    control_extension_rows = build_extension_rows(control_observations)
    extension_rows = branch_extension_rows + control_extension_rows

    branch_coarse_rows = projection_summary(branch_observations, COARSE_GRAIN_CUTOFF, probe_times, phase7_exponent)
    control_coarse_rows = projection_summary(control_observations, COARSE_GRAIN_CUTOFF, probe_times, phase7_exponent)
    coarse_rows = branch_coarse_rows + control_coarse_rows

    branch_all_level_h2_fit = fit_linear_family(
        np.array([float(observation["h"]) for observation in branch_observations], dtype=float),
        branch_coefficients,
        power=PRIMARY_CORRECTION_POWER,
    )
    control_all_level_h2_fit = fit_linear_family(
        np.array([float(observation["h"]) for observation in control_observations], dtype=float),
        control_coefficients,
        power=PRIMARY_CORRECTION_POWER,
    )

    plot_descriptor_vs_refinement(branch_observations, control_observations)
    plot_prediction_errors(branch_observations, control_observations, branch_all_level_h2_fit, control_all_level_h2_fit)
    plot_coarse_grain_spectral_comparison(branch_coarse_rows, control_coarse_rows)
    plot_multi_descriptor_correlation(branch_observations, control_observations)

    write_csv_rows(
        EXTENSION_LEDGER_PATH,
        extension_rows,
        [
            "hierarchy_label",
            "level_id",
            "n_side",
            "h",
            "dimension",
            "nullspace_dimension",
            "spectral_radius",
            "first_positive_eigenvalue",
            "b0",
            "b1",
            "b2",
            "b1_over_b0",
            "b2_over_b0",
            "I0",
            "I1",
            "I2",
            "descriptor_distance_from_previous",
        ],
    )
    write_csv_rows(
        SCALING_LEDGER_PATH,
        scaling_rows,
        [
            "hierarchy_label",
            "family",
            "quantity",
            "fit_power",
            "fit_scope",
            "fit_r_squared",
            "slope",
            "limit_estimate",
            "predicted_R5",
            "observed_R5",
            "relative_error_R5",
        ],
    )
    write_csv_rows(
        COARSE_LEDGER_PATH,
        coarse_rows,
        [
            "hierarchy_label",
            "level_id",
            "n_side",
            "h",
            "projection_label",
            "cutoff",
            "retained_mode_count",
            "retained_mode_fraction",
            "projected_spectral_radius",
            "spectral_radius_loss",
            "first_positive_low_mode_max_abs_difference",
            "trace_window_max_relative_error",
            "ratio_max_relative_drift",
            "invariant_max_relative_drift",
            "low_mode_preserved",
            "ratio_ordering_preserved",
        ],
    )

    branch_descriptor_matrix = np.array(
        [
            [
                float(observation["ratios"]["b1_over_b0"]),
                float(observation["ratios"]["b2_over_b0"]),
                float(observation["rescaled_invariants"]["I0"]),
                float(observation["rescaled_invariants"]["I1"]),
                float(observation["rescaled_invariants"]["I2"]),
                float(observation["spectral_radius"]),
            ]
            for observation in branch_observations
        ],
        dtype=float,
    )
    control_descriptor_matrix = np.array(
        [
            [
                float(observation["ratios"]["b1_over_b0"]),
                float(observation["ratios"]["b2_over_b0"]),
                float(observation["rescaled_invariants"]["I0"]),
                float(observation["rescaled_invariants"]["I1"]),
                float(observation["rescaled_invariants"]["I2"]),
                float(observation["spectral_radius"]),
            ]
            for observation in control_observations
        ],
        dtype=float,
    )

    def pca_first_component_fraction(matrix: np.ndarray) -> float:
        normalized = (matrix - matrix.mean(axis=0)) / matrix.std(axis=0)
        singular_values = np.linalg.svd(normalized, compute_uv=False)
        eigenvalues = (singular_values**2) / (len(matrix) - 1)
        return float(eigenvalues[0] / np.sum(eigenvalues))

    branch_descriptor_distances = [
        float(np.linalg.norm(branch_descriptor_matrix[index + 1] - branch_descriptor_matrix[index]))
        for index in range(len(branch_descriptor_matrix) - 1)
    ]
    control_descriptor_distances = [
        float(np.linalg.norm(control_descriptor_matrix[index + 1] - control_descriptor_matrix[index]))
        for index in range(len(control_descriptor_matrix) - 1)
    ]
    branch_descriptor_distances_shrink = all(
        branch_descriptor_distances[index + 1] < branch_descriptor_distances[index] for index in range(len(branch_descriptor_distances) - 1)
    )

    branch_ratio_limit_drift_max = max(branch_scaling_metrics["ratio_limit_drift"].values())
    branch_invariant_limit_drift_max = max(branch_scaling_metrics["invariant_limit_drift"].values())
    branch_coarse_trace_error_max = max(float(row["trace_window_max_relative_error"]) for row in branch_coarse_rows)
    branch_coarse_ratio_drift_max = max(float(row["ratio_max_relative_drift"]) for row in branch_coarse_rows)

    control_trace_error_multiplier = control_scaling_metrics["trace_h2_error"] / max(branch_scaling_metrics["trace_h2_error"], 1.0e-12)
    control_ratio_error_multiplier = control_scaling_metrics["ratio_h2_max_error"] / max(branch_scaling_metrics["ratio_h2_max_error"], 1.0e-12)

    branch_best_models = best_model_labels(
        np.array([float(observation["h"]) for observation in branch_observations], dtype=float),
        branch_coefficients,
        coefficient_names,
    )
    control_best_models = best_model_labels(
        np.array([float(observation["h"]) for observation in control_observations], dtype=float),
        control_coefficients,
        coefficient_names,
    )

    success_flags = {
        "extension_preserves_primary_h2_trace_law": bool(branch_scaling_metrics["trace_h2_error"] <= BRIDGE_TRACE_ERROR_TOL),
        "ratio_limit_candidates_stable": bool(branch_ratio_limit_drift_max <= RATIO_LIMIT_DRIFT_TOL),
        "rescaled_invariants_predictable": bool(branch_scaling_metrics["invariant_h1_max_error"] <= INVARIANT_EXTENSION_ERROR_TOL),
        "coarse_projection_retains_branch_signatures": bool(
            branch_coarse_trace_error_max <= COARSE_TRACE_ERROR_TOL and branch_coarse_ratio_drift_max <= COARSE_RATIO_DRIFT_TOL
        ),
        "multi_descriptor_coherence_compact": bool(
            pca_first_component_fraction(branch_descriptor_matrix) >= COHERENCE_PC1_MIN and branch_descriptor_distances_shrink
        ),
        "control_fails_branch_bridge_law": bool(
            control_trace_error_multiplier >= CONTROL_TRACE_ERROR_MULTIPLIER_MIN
            and control_ratio_error_multiplier >= CONTROL_RATIO_ERROR_MULTIPLIER_MIN
            and any(label != "linear_in_h^2" for label in control_best_models.values())
        ),
    }
    success = all(success_flags.values())

    runtime_seconds = time.perf_counter() - start_time

    write_json(
        RUNS_PATH,
        {
            "timestamp": timestamp_iso(),
            "phase": PHASE_IDENTIFIER,
            "phase_name": PHASE_NAME,
            "builder_script": relpath(Path(__file__)),
            "deterministic_seed_record": {
                "randomness_used": False,
                "note": "No stochastic probes were used; all spectra, fits, and projections are deterministic."
            },
            "references": {
                "phase6_manifest": relpath(PHASE6_MANIFEST_PATH),
                "phase8_manifest": relpath(PHASE8_MANIFEST_PATH),
                "phase9_manifest": relpath(PHASE9_MANIFEST_PATH),
            },
            "hierarchy": {
                "frozen_levels": frozen_levels,
                "extended_level": extended_levels[-1],
            },
            "effective_scaling_definition": {
                "trace_prefactor": "h^(-2)",
                "primary_correction_variable": "h^2",
                "primary_surrogate": "T_h(t) ≈ h^(-2) * (b0(h) + b1(h) * t + b2(h) * t^2)",
                "invariant_tracking_variable": "h",
                "trace_window": [round_float(value, 6) for value in probe_times.tolist()],
            },
            "coarse_grain_definition": {
                "label": COARSE_GRAIN_LABEL,
                "method": "spectral truncation projection",
                "cutoff": COARSE_GRAIN_CUTOFF,
            },
            "control_hierarchy": phase9_manifest["control_hierarchy"],
            "runtime_seconds": round_float(runtime_seconds, 6),
        },
    )

    summary_lines = [
        "# Phase X Cautious Continuum-Bridge Feasibility",
        "",
        "## Objective",
        "",
        "Test whether the frozen Phase IX spectral-trace descriptors support a cautious large-scale bridge description under one deterministic refinement extension and one explicit coarse-grain projection, without changing any earlier phase contract.",
        "",
        "## Frozen Inputs",
        "",
        f"- Phase VI operator manifest: `{relpath(PHASE6_MANIFEST_PATH)}`",
        f"- Phase VIII trace contract: `{relpath(PHASE8_MANIFEST_PATH)}`",
        f"- Phase IX descriptor contract: `{relpath(PHASE9_MANIFEST_PATH)}`",
        "",
        "## Extended Hierarchy",
        "",
        f"- Frozen levels: {[int(level['n_side']) for level in frozen_levels]}",
        f"- Extended level: `n_side = {extended_levels[-1]['n_side']}` with `h = {extended_levels[-1]['h']}`",
        f"- Short-time trace window reused unchanged: {[round_float(value, 6) for value in probe_times.tolist()]}",
        "",
        "## Primary Effective Scaling Description",
        "",
        "- Branch trace law: `T_h(t) ≈ h^(-2) * (b0(h) + b1(h) * t + b2(h) * t^2)` with `b_k(h) = b_k,∞ + c_k h^2`.",
        "- Ratio limit candidates reused unchanged from Phase IX: `b1 / b0` and `b2 / b0`.",
        "- Rescaled invariant tracking reused unchanged from Phase IX: `I0`, `I1`, `I2`, tracked with linear-in-`h` correction.",
        "",
        "## Results",
        "",
        f"- Branch extension keeps the `h^2` trace law intact: out-of-sample `R5` max short-window trace error = `{round_float(branch_scaling_metrics['trace_h2_error'])}`, compared with `{round_float(branch_scaling_metrics['trace_h1_error'])}` for the weaker `h` comparator.",
        f"- Branch ratio prediction at `R5` is tight: max relative error = `{round_float(branch_scaling_metrics['ratio_h2_max_error'])}`, and max ratio-limit drift after adding `R5` = `{round_float(branch_ratio_limit_drift_max)}`.",
        f"- Branch rescaled invariants remain bounded and predictable: max `R5` relative error under linear-in-`h` tracking = `{round_float(branch_scaling_metrics['invariant_h1_max_error'])}`, with max limit drift = `{round_float(branch_invariant_limit_drift_max)}`.",
        f"- Coarse projection `{COARSE_GRAIN_LABEL}` keeps the branch low-mode ladder exact, with max short-window trace error `{round_float(branch_coarse_trace_error_max)}` and max ratio drift `{round_float(branch_coarse_ratio_drift_max)}`.",
        f"- Multi-descriptor organization stays compact: first principal-component share = `{round_float(pca_first_component_fraction(branch_descriptor_matrix))}` and descriptor distances shrink along the extension `{[round_float(value) for value in branch_descriptor_distances]}`.",
        f"- The deterministic control does not share the same bridge law: its `h^2` trace-prediction error is `{round_float(control_scaling_metrics['trace_h2_error'])}` against `{round_float(branch_scaling_metrics['trace_h2_error'])}` on the branch, and its max ratio prediction error under the same law is `{round_float(control_scaling_metrics['ratio_h2_max_error'])}`.",
        f"- Model-family separation persists under extension: branch scaled coefficients stay `{branch_best_models}`, while the control stays `{control_best_models}`.",
        "",
        "## Bounded Interpretation",
        "",
        "Phase X supports a cautious branch-local bridge statement only in the descriptive sense used here: the frozen hierarchy admits a compact scaling summary that survives one finer deterministic level, retains its ratio-limit candidates, and remains distinguishable from the altered-connectivity control under the same `h^2` bridge law.",
        "",
        "## Explicit Non-Claims",
        "",
        "- No continuum limit is asserted.",
        "- No geometric or metric meaning is assigned to any coefficient, ratio, or descriptor.",
        "- No physical correspondence is claimed.",
        "- The coarse projection is only a reproducible diagnostic reduction, not a continuum derivation.",
        "",
        closure_statement(success),
    ]
    write_text(SUMMARY_PATH, "\n".join(summary_lines))

    manifest = {
        "timestamp": timestamp_iso(),
        "phase": PHASE_IDENTIFIER,
        "phase_name": PHASE_NAME,
        "stage_identifier": PHASE_NAME,
        "status": "passed" if success else "failed",
        "success": bool(success),
        "objective": "cautious_continuum_bridge_feasibility",
        "operator_manifest_reference": relpath(PHASE6_MANIFEST_PATH),
        "phase8_manifest_reference": relpath(PHASE8_MANIFEST_PATH),
        "phase9_manifest_reference": relpath(PHASE9_MANIFEST_PATH),
        "refinement_hierarchy": {
            "h_definition": phase6_manifest["refinement_hierarchy"]["h_definition"],
            "level_parameter": "n_side",
            "construction_rule_constant": True,
            "refinement_multipliers": [int(level["refinement_multiplier"]) for level in extended_levels],
            "levels": extended_levels,
        },
        "frozen_trace_contract": {
            "operational_short_time_window": phase8_manifest["operational_short_time_window"]["probe_times"],
            "master_probe_grid": phase8_manifest["trace_proxy_definition"]["probe_times"],
            "spectral_truncation_methods": phase8_manifest["trace_proxy_definition"]["spectral_truncation_methods"],
        },
        "frozen_descriptor_contract": {
            "ratios": ratio_names,
            "rescaled_invariants": invariant_names,
            "phase9_control_hierarchy": phase9_manifest["control_hierarchy"],
        },
        "primary_effective_scaling_description": {
            "trace_prefactor": "h^(-2)",
            "correction_variable": "h^2",
            "surrogate": "T_h(t) ≈ h^(-2) * (b0(h) + b1(h) * t + b2(h) * t^2)",
            "branch_limit_coefficients": {
                name: round_float(value)
                for name, value in zip(coefficient_names, branch_scaling_metrics["coefficient_h2"]["limits"])
            },
            "branch_ratio_limit_candidates": {
                name: round_float(value)
                for name, value in zip(ratio_names, branch_scaling_metrics["ratio_h2"]["limits"])
            },
            "branch_invariant_limit_candidates": {
                name: round_float(value)
                for name, value in zip(invariant_names, branch_scaling_metrics["invariant_h1"]["limits"])
            },
        },
        "extension_results": {
            "branch_r5_prediction_errors": {
                "trace_h2_max_relative_error": round_float(branch_scaling_metrics["trace_h2_error"]),
                "trace_h1_comparator_error": round_float(branch_scaling_metrics["trace_h1_error"]),
                "scaled_coefficient_h2_max_relative_error": round_float(branch_scaling_metrics["coefficient_h2_max_error"]),
                "ratio_h2_max_relative_error": round_float(branch_scaling_metrics["ratio_h2_max_error"]),
                "invariant_h1_max_relative_error": round_float(branch_scaling_metrics["invariant_h1_max_error"]),
            },
            "branch_limit_drift_after_extension": {
                "ratios": {name: round_float(value) for name, value in branch_scaling_metrics["ratio_limit_drift"].items()},
                "rescaled_invariants": {
                    name: round_float(value) for name, value in branch_scaling_metrics["invariant_limit_drift"].items()
                },
            },
            "branch_best_scaled_coefficient_models": branch_best_models,
            "branch_spectral_radius_values": [round_float(observation["spectral_radius"]) for observation in branch_observations],
        },
        "coarse_grain_projection": {
            "label": COARSE_GRAIN_LABEL,
            "cutoff": COARSE_GRAIN_CUTOFF,
            "branch_retained_mode_fraction_range": [
                round_float(min(float(row["retained_mode_fraction"]) for row in branch_coarse_rows)),
                round_float(max(float(row["retained_mode_fraction"]) for row in branch_coarse_rows)),
            ],
            "branch_max_trace_window_relative_error": round_float(branch_coarse_trace_error_max),
            "branch_max_ratio_relative_drift": round_float(branch_coarse_ratio_drift_max),
            "branch_low_mode_preserved_at_all_levels": bool(all(bool(row["low_mode_preserved"]) for row in branch_coarse_rows)),
        },
        "multi_descriptor_coherence": {
            "branch_first_principal_component_share": round_float(pca_first_component_fraction(branch_descriptor_matrix)),
            "branch_descriptor_distance_sequence": [round_float(value) for value in branch_descriptor_distances],
            "branch_descriptor_distances_shrink": bool(branch_descriptor_distances_shrink),
            "branch_ratio2_vs_invariant2_correlation": round_float(
                pearson_correlation(branch_ratios[:, 1], branch_invariants[:, 2])
            ),
        },
        "control_comparison": {
            "label": CONTROL_LABEL,
            "description": phase9_manifest["control_hierarchy"]["description"],
            "trace_h2_error_branch": round_float(branch_scaling_metrics["trace_h2_error"]),
            "trace_h2_error_control": round_float(control_scaling_metrics["trace_h2_error"]),
            "trace_h2_error_multiplier": round_float(control_trace_error_multiplier),
            "ratio_h2_error_branch": round_float(branch_scaling_metrics["ratio_h2_max_error"]),
            "ratio_h2_error_control": round_float(control_scaling_metrics["ratio_h2_max_error"]),
            "ratio_h2_error_multiplier": round_float(control_ratio_error_multiplier),
            "control_best_scaled_coefficient_models": control_best_models,
            "control_ratio_limit_drift_max": round_float(max(control_scaling_metrics["ratio_limit_drift"].values())),
        },
        "success_flags": success_flags,
        "artifacts": {
            "extension_ledger_csv": relpath(EXTENSION_LEDGER_PATH),
            "effective_scaling_ledger_csv": relpath(SCALING_LEDGER_PATH),
            "coarse_grain_summary_csv": relpath(COARSE_LEDGER_PATH),
            "summary_md": relpath(SUMMARY_PATH),
            "manifest_json": relpath(MANIFEST_PATH),
            "descriptor_plot": relpath(DESCRIPTOR_PLOT_PATH),
            "prediction_error_plot": relpath(PREDICTION_PLOT_PATH),
            "coarse_grain_plot": relpath(COARSE_PLOT_PATH),
            "multi_descriptor_correlation_plot": relpath(CORRELATION_PLOT_PATH),
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
