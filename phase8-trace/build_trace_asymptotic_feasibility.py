#!/usr/bin/env python3

from __future__ import annotations

import itertools
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

from haos_core import build_graph, read_json, relpath, write_csv_rows, write_json

PHASE_ROOT = ROOT / "phase8-trace"
RUNS_ROOT = PHASE_ROOT / "runs"
PLOTS_ROOT = PHASE_ROOT / "plots"

PHASE5_MANIFEST_PATH = ROOT / "phase5-readout" / "phase5_authoritative_manifest.json"
PHASE6_MANIFEST_PATH = ROOT / "phase6-operator" / "phase6_operator_manifest.json"
PHASE7_MANIFEST_PATH = ROOT / "phase7-spectral" / "phase7_spectral_manifest.json"

LEDGER_PATH = RUNS_ROOT / "phase8_trace_scaling_ledger.csv"
RUNS_PATH = RUNS_ROOT / "phase8_runs.json"
SUMMARY_PATH = PHASE_ROOT / "phase8_trace_summary.md"
MANIFEST_PATH = PHASE_ROOT / "phase8_trace_manifest.json"

TRACE_PLOT_PATH = PLOTS_ROOT / "phase8_trace_vs_t.svg"
SLOPE_PLOT_PATH = PLOTS_ROOT / "phase8_log_trace_slope_vs_t.svg"
COEFF_PLOT_PATH = PLOTS_ROOT / "phase8_coefficient_like_terms_vs_refinement.svg"

PHASE_NAME = "phase8-trace"
PHASE_IDENTIFIER = 8
TRACE_METHOD_LABEL = "exact_torus_mode_spectrum"
TRACE_REDUCED_LABEL = "kernel_subtracted"
TRUNCATION_METHODS = (
    {"label": "full_exact", "cutoff": None},
    {"label": "lambda_le_7_5", "cutoff": 7.5},
    {"label": "lambda_le_7_8", "cutoff": 7.8},
)
MASTER_PROBE_TIMES = (
    0.01,
    0.015,
    0.02,
    0.03,
    0.05,
    0.075,
    0.1,
    0.15,
    0.2,
    0.3,
    0.5,
    0.75,
    1.0,
)
MIN_WINDOW_POINTS = 5
SLOPE_DISPERSION_TOL = 0.01
NORMALIZED_TRACE_REL_RANGE_TOL = 0.01
KERNEL_SHARE_TOL = 0.02
KERNEL_SLOPE_DIFF_TOL = 0.005
TRUNCATION_REL_ERROR_TOLERANCES = {
    "lambda_le_7_5": 0.04,
    "lambda_le_7_8": 0.015,
}
TRUNCATION_SLOPE_DIFF_TOL = 0.012
TRUNCATION_COEFFICIENT_REL_DIFF_TOL = 0.12
COEFFICIENT_R2_TOL = 0.99998
COEFFICIENT_SUBGRID_REL_DRIFT_TOL = 0.03
COEFFICIENT_REL_SPAN_TOL = 0.03
CLAIM_BOUNDARY = (
    "Phase VIII authority is limited to short-time trace asymptotic feasibility diagnostics on "
    "the frozen branch-local cochain-Laplacian hierarchy. It does not assert continuum limits, "
    "geometric coefficients, universal asymptotic coefficients, Seeley-DeWitt structure, or any "
    "physical correspondence."
)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def timestamp_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def closure_statement(success: bool) -> str:
    if success:
        return "Phase VIII establishes short-time trace asymptotic feasibility for the frozen operator hierarchy."
    return "Phase VIII does not yet establish short-time trace asymptotic feasibility for the frozen operator hierarchy."


def round_float(value: float, digits: int = 12) -> float:
    return round(float(value), digits)


def trend_label(values: list[float], tolerance: float = 1.0e-12) -> str:
    nondecreasing = all(values[index + 1] >= values[index] - tolerance for index in range(len(values) - 1))
    nonincreasing = all(values[index + 1] <= values[index] + tolerance for index in range(len(values) - 1))
    if nondecreasing:
        return "nondecreasing"
    if nonincreasing:
        return "nonincreasing"
    return "mixed"


def load_frozen_manifests() -> dict[str, Any]:
    phase5_manifest = read_json(PHASE5_MANIFEST_PATH)
    phase6_manifest = read_json(PHASE6_MANIFEST_PATH)
    phase7_manifest = read_json(PHASE7_MANIFEST_PATH)

    if not bool(phase5_manifest.get("success")):
        raise ValueError("Phase VIII requires the frozen successful Phase V authority manifest.")
    if not bool(phase6_manifest.get("success")):
        raise ValueError("Phase VIII requires the frozen successful Phase VI operator manifest.")
    if phase6_manifest.get("selected_operator_class") != "cochain_laplacian":
        raise ValueError("Phase VIII requires the frozen cochain-Laplacian operator hierarchy.")
    if not bool(phase7_manifest.get("success")):
        raise ValueError("Phase VIII requires the frozen successful Phase VII spectral manifest.")

    return {
        "phase5_manifest": phase5_manifest,
        "phase6_manifest": phase6_manifest,
        "phase7_manifest": phase7_manifest,
    }


def scalar_mode_spectrum(n_side: int, epsilon: float) -> tuple[np.ndarray, float]:
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
    return np.sort(scalar_values), edge_weight


def exact_cochain_spectrum(n_side: int, epsilon: float) -> tuple[np.ndarray, float]:
    scalar_values, edge_weight = scalar_mode_spectrum(n_side, epsilon)
    return np.sort(np.repeat(scalar_values, 4)), edge_weight


def trace_from_spectrum(spectrum: np.ndarray, probe_times: np.ndarray) -> np.ndarray:
    return np.array([float(np.exp(-probe_time * spectrum).sum()) for probe_time in probe_times], dtype=float)


def local_log_slope(probe_times: np.ndarray, traces: np.ndarray) -> np.ndarray:
    return np.gradient(np.log(traces), np.log(probe_times))


def quadratic_fit(probe_times: np.ndarray, values: np.ndarray) -> dict[str, Any]:
    coefficients = np.polyfit(probe_times, values, 2)
    fitted = np.polyval(coefficients, probe_times)
    denominator = float(np.sum((values - values.mean()) ** 2))
    if denominator <= 1.0e-18:
        r_squared = 1.0
    else:
        r_squared = 1.0 - float(np.sum((values - fitted) ** 2) / denominator)
    c0, c1, c2 = coefficients[::-1]
    return {
        "b0": round_float(c0),
        "b1": round_float(c1),
        "b2": round_float(c2),
        "r_squared": round_float(r_squared),
    }


def selected_window_mask(probe_times: np.ndarray, selected_times: np.ndarray) -> np.ndarray:
    return np.array(
        [
            any(math.isclose(float(probe_time), float(candidate), rel_tol=0.0, abs_tol=1.0e-12) for candidate in selected_times)
            for probe_time in probe_times
        ],
        dtype=bool,
    )


def window_local_log_slope(probe_times: np.ndarray, traces: np.ndarray, selected_times: np.ndarray) -> np.ndarray:
    mask = selected_window_mask(probe_times, selected_times)
    return local_log_slope(selected_times, traces[mask])


def direct_dense_validation(n_side: int, epsilon: float, probe_times: np.ndarray) -> dict[str, Any]:
    graph = build_graph(
        {
            "kind": "dk2d_periodic",
            "n_side": int(n_side),
            "epsilon": float(epsilon),
            "cycle_phase_x": 0.0,
            "cycle_phase_y": 0.0,
        }
    )
    dense_eigenvalues = np.linalg.eigvalsh(graph.delta_h.toarray().real)
    exact_eigenvalues, _ = exact_cochain_spectrum(n_side, epsilon)
    max_abs_eigenvalue_difference = float(np.max(np.abs(np.sort(dense_eigenvalues) - exact_eigenvalues)))
    dense_traces = trace_from_spectrum(dense_eigenvalues, probe_times)
    exact_traces = trace_from_spectrum(exact_eigenvalues, probe_times)
    max_abs_trace_difference = float(np.max(np.abs(dense_traces - exact_traces)))
    return {
        "validated_level_n_side": int(n_side),
        "max_abs_eigenvalue_difference": round_float(max_abs_eigenvalue_difference),
        "max_abs_trace_difference": round_float(max_abs_trace_difference),
    }


def compute_level_summary(level_index: int, level: dict[str, Any], epsilon: float, probe_times: np.ndarray) -> dict[str, Any]:
    n_side = int(level["n_side"])
    h_value = float(level["h"])
    spectrum, edge_weight = exact_cochain_spectrum(n_side, epsilon)
    full_traces = trace_from_spectrum(spectrum, probe_times)
    reduced_traces = full_traces - 4.0

    methods: dict[str, dict[str, Any]] = {}
    for method in TRUNCATION_METHODS:
        label = str(method["label"])
        cutoff = method["cutoff"]
        truncated_spectrum = spectrum if cutoff is None else spectrum[spectrum <= float(cutoff)]
        traces = trace_from_spectrum(truncated_spectrum, probe_times)
        methods[label] = {
            "cutoff": None if cutoff is None else round_float(float(cutoff), 6),
            "included_mode_count": int(len(truncated_spectrum)),
            "traces": [round_float(value) for value in traces.tolist()],
            "log_trace_slopes": [round_float(value) for value in local_log_slope(probe_times, traces).tolist()],
        }

    kernel_share = 4.0 / full_traces
    reduced_slope = local_log_slope(probe_times, reduced_traces)
    full_slope = np.array(methods["full_exact"]["log_trace_slopes"], dtype=float)

    return {
        "level_id": f"R{level_index}",
        "n_side": int(n_side),
        "h": round_float(h_value),
        "dimension": int(4 * n_side * n_side),
        "edge_weight": round_float(edge_weight),
        "nullspace_dimension": 4,
        "spectrum_max": round_float(float(spectrum[-1])),
        "trace_monotone_in_t": bool(all(full_traces[index + 1] < full_traces[index] for index in range(len(full_traces) - 1))),
        "probe_times": [round_float(value, 6) for value in probe_times.tolist()],
        "methods": methods,
        TRACE_REDUCED_LABEL: {
            "traces": [round_float(value) for value in reduced_traces.tolist()],
            "log_trace_slopes": [round_float(value) for value in reduced_slope.tolist()],
            "kernel_share": [round_float(value) for value in kernel_share.tolist()],
            "max_kernel_share": round_float(float(np.max(kernel_share))),
            "max_log_slope_difference_vs_full": round_float(float(np.max(np.abs(reduced_slope - full_slope)))),
        },
    }


def collect_time_metrics(levels: list[dict[str, Any]], probe_times: np.ndarray) -> dict[str, Any]:
    full_traces = {
        level["level_id"]: np.array(level["methods"]["full_exact"]["traces"], dtype=float)
        for level in levels
    }
    full_slopes = {
        level["level_id"]: np.array(level["methods"]["full_exact"]["log_trace_slopes"], dtype=float)
        for level in levels
    }
    reduced_metrics = {
        level["level_id"]: level[TRACE_REDUCED_LABEL]
        for level in levels
    }

    per_time = []
    for index, probe_time in enumerate(probe_times):
        scaled_traces = []
        slope_values = []
        kernel_shares = []
        record: dict[str, Any] = {"probe_time": round_float(float(probe_time), 6)}
        for level in levels:
            level_id = str(level["level_id"])
            trace_value = full_traces[level_id][index]
            scaled_traces.append((float(level["h"]) ** 2) * trace_value)
            slope_values.append(full_slopes[level_id][index])
            kernel_shares.append(float(reduced_metrics[level_id]["kernel_share"][index]))
        record["max_slope_dispersion"] = round_float(float(max(slope_values) - min(slope_values)))
        record["max_normalized_trace_relative_range"] = round_float(
            float((max(scaled_traces) - min(scaled_traces)) / max(float(np.mean(scaled_traces)), 1.0e-12))
        )
        record["max_kernel_share"] = round_float(float(max(kernel_shares)))
        for method in TRUNCATION_METHODS[1:]:
            label = str(method["label"])
            cutoff_traces = {
                level["level_id"]: np.array(level["methods"][label]["traces"], dtype=float)
                for level in levels
            }
            cutoff_slopes = {
                level["level_id"]: np.array(level["methods"][label]["log_trace_slopes"], dtype=float)
                for level in levels
            }
            rel_errors = []
            slope_drifts = []
            for level in levels:
                level_id = str(level["level_id"])
                full_value = full_traces[level_id][index]
                truncated_value = cutoff_traces[level_id][index]
                rel_errors.append(abs(truncated_value - full_value) / full_value)
                slope_drifts.append(abs(cutoff_slopes[level_id][index] - full_slopes[level_id][index]))
            record[f"{label}_max_relative_trace_error"] = round_float(float(max(rel_errors)))
            record[f"{label}_max_log_slope_difference"] = round_float(float(max(slope_drifts)))
        per_time.append(record)
    return {"per_time": per_time}


def identify_operational_window(levels: list[dict[str, Any]], probe_times: np.ndarray, time_metrics: dict[str, Any]) -> dict[str, Any]:
    best_candidate: dict[str, Any] | None = None
    metric_by_time = {float(record["probe_time"]): record for record in time_metrics["per_time"]}

    for start_index in range(len(probe_times)):
        for end_index in range(start_index + MIN_WINDOW_POINTS - 1, len(probe_times)):
            window_times = probe_times[start_index : end_index + 1]
            window_metrics = [metric_by_time[float(round_float(probe_time, 6))] for probe_time in window_times]

            max_slope_dispersion = max(float(record["max_slope_dispersion"]) for record in window_metrics)
            max_normalized_range = max(float(record["max_normalized_trace_relative_range"]) for record in window_metrics)
            max_kernel_share = max(float(record["max_kernel_share"]) for record in window_metrics)
            max_trunc_75 = max(float(record["lambda_le_7_5_max_relative_trace_error"]) for record in window_metrics)
            max_trunc_78 = max(float(record["lambda_le_7_8_max_relative_trace_error"]) for record in window_metrics)

            if max_slope_dispersion > SLOPE_DISPERSION_TOL:
                continue
            if max_normalized_range > NORMALIZED_TRACE_REL_RANGE_TOL:
                continue
            if max_kernel_share > KERNEL_SHARE_TOL:
                continue
            if max_trunc_75 > TRUNCATION_REL_ERROR_TOLERANCES["lambda_le_7_5"]:
                continue
            if max_trunc_78 > TRUNCATION_REL_ERROR_TOLERANCES["lambda_le_7_8"]:
                continue

            subgrid_times = window_times[::2]
            min_fit_r_squared = 1.0
            max_subgrid_rel_drift = 0.0
            per_level_fits = []
            for level in levels:
                full_trace = np.array(level["methods"]["full_exact"]["traces"], dtype=float)
                full_window_values = (float(level["h"]) ** 2) * full_trace[start_index : end_index + 1]
                fit = quadratic_fit(window_times, full_window_values)
                min_fit_r_squared = min(min_fit_r_squared, float(fit["r_squared"]))

                subgrid_values = full_window_values[::2]
                subgrid_fit = quadratic_fit(subgrid_times, subgrid_values)
                coefficient_drifts = []
                for key in ("b0", "b1", "b2"):
                    denominator = max(abs(float(fit[key])), 1.0e-12)
                    coefficient_drifts.append(abs(float(subgrid_fit[key]) - float(fit[key])) / denominator)
                max_subgrid_rel_drift = max(max_subgrid_rel_drift, max(coefficient_drifts))
                per_level_fits.append(
                    {
                        "level_id": level["level_id"],
                        "fit": fit,
                        "subgrid_fit": subgrid_fit,
                        "max_subgrid_rel_drift": round_float(max(coefficient_drifts)),
                    }
                )

            if min_fit_r_squared < COEFFICIENT_R2_TOL:
                continue
            if max_subgrid_rel_drift > COEFFICIENT_SUBGRID_REL_DRIFT_TOL:
                continue

            score = (
                len(window_times),
                -max_subgrid_rel_drift,
                -max_slope_dispersion,
                -max_trunc_75,
            )
            candidate = {
                "score": score,
                "start_index": start_index,
                "end_index": end_index,
                "probe_times": [round_float(value, 6) for value in window_times.tolist()],
                "max_slope_dispersion": round_float(max_slope_dispersion),
                "max_normalized_trace_relative_range": round_float(max_normalized_range),
                "max_kernel_share": round_float(max_kernel_share),
                "max_lambda_le_7_5_relative_trace_error": round_float(max_trunc_75),
                "max_lambda_le_7_8_relative_trace_error": round_float(max_trunc_78),
                "min_full_fit_r_squared": round_float(min_fit_r_squared),
                "max_full_fit_subgrid_rel_drift": round_float(max_subgrid_rel_drift),
                "full_fit_window": per_level_fits,
            }
            if best_candidate is None or candidate["score"] > best_candidate["score"]:
                best_candidate = candidate

    if best_candidate is None:
        raise ValueError("Phase VIII could not identify an operational short-time window from the frozen probe grid.")
    return best_candidate


def fit_all_methods(levels: list[dict[str, Any]], probe_times: np.ndarray, operational_window: dict[str, Any]) -> dict[str, Any]:
    selected_times = np.array(operational_window["probe_times"], dtype=float)
    selected_mask = selected_window_mask(probe_times, selected_times)

    method_results: dict[str, Any] = {}
    full_coefficients: dict[str, np.ndarray] = {}

    for method in TRUNCATION_METHODS:
        label = str(method["label"])
        per_level = []
        coefficient_rows = []
        for level in levels:
            traces = np.array(level["methods"][label]["traces"], dtype=float)
            scaled_values = (float(level["h"]) ** 2) * traces[selected_mask]
            fit = quadratic_fit(selected_times, scaled_values)
            coefficient_rows.append([float(fit["b0"]), float(fit["b1"]), float(fit["b2"])])
            per_level.append(
                {
                    "level_id": level["level_id"],
                    "n_side": level["n_side"],
                    "h": level["h"],
                    "fit": fit,
                }
            )

        coefficient_matrix = np.array(coefficient_rows, dtype=float)
        relative_spans = []
        trends = {}
        for index, key in enumerate(("b0", "b1", "b2")):
            values = coefficient_matrix[:, index].tolist()
            trends[key] = trend_label(values)
            relative_spans.append(
                float((np.max(coefficient_matrix[:, index]) - np.min(coefficient_matrix[:, index])) / max(abs(float(np.mean(coefficient_matrix[:, index]))), 1.0e-12))
            )
        summary = {
            "per_level": per_level,
            "relative_spans": {
                "b0": round_float(relative_spans[0]),
                "b1": round_float(relative_spans[1]),
                "b2": round_float(relative_spans[2]),
            },
            "trends": trends,
            "max_relative_span": round_float(max(relative_spans)),
            "min_r_squared": round_float(min(float(item["fit"]["r_squared"]) for item in per_level)),
        }
        method_results[label] = summary
        if label == "full_exact":
            full_coefficients = {
                level["level_id"]: np.array(
                    [float(item["fit"]["b0"]), float(item["fit"]["b1"]), float(item["fit"]["b2"])],
                    dtype=float,
                )
                for level, item in zip(levels, per_level)
            }

    robustness = {}
    for method in TRUNCATION_METHODS[1:]:
        label = str(method["label"])
        max_coefficient_rel_diff = 0.0
        max_slope_diff = 0.0
        for level, item in zip(levels, method_results[label]["per_level"]):
            level_id = str(level["level_id"])
            candidate_vector = np.array(
                [float(item["fit"]["b0"]), float(item["fit"]["b1"]), float(item["fit"]["b2"])],
                dtype=float,
            )
            full_vector = full_coefficients[level_id]
            rel_diff = np.max(np.abs((candidate_vector - full_vector) / np.maximum(np.abs(full_vector), 1.0e-12)))
            max_coefficient_rel_diff = max(max_coefficient_rel_diff, float(rel_diff))

            full_traces = np.array(level["methods"]["full_exact"]["traces"], dtype=float)
            candidate_traces = np.array(level["methods"][label]["traces"], dtype=float)
            full_slopes = window_local_log_slope(probe_times, full_traces, selected_times)
            candidate_slopes = window_local_log_slope(probe_times, candidate_traces, selected_times)
            slope_diff = np.max(np.abs(candidate_slopes - full_slopes))
            max_slope_diff = max(max_slope_diff, float(slope_diff))

        robustness[label] = {
            "cutoff": method["cutoff"],
            "max_relative_coefficient_difference_vs_full": round_float(max_coefficient_rel_diff),
            "max_log_slope_difference_vs_full": round_float(max_slope_diff),
        }
    return {"methods": method_results, "truncation_robustness": robustness}


def compute_window_kernel_sensitivity(levels: list[dict[str, Any]], probe_times: np.ndarray, operational_window: dict[str, Any]) -> dict[str, Any]:
    selected_times = np.array(operational_window["probe_times"], dtype=float)
    selected_mask = selected_window_mask(probe_times, selected_times)

    max_kernel_share = 0.0
    max_slope_diff = 0.0
    per_level = []
    for level in levels:
        full_traces = np.array(level["methods"]["full_exact"]["traces"], dtype=float)
        reduced_traces = np.array(level[TRACE_REDUCED_LABEL]["traces"], dtype=float)
        kernel_share = np.array(level[TRACE_REDUCED_LABEL]["kernel_share"], dtype=float)[selected_mask]
        full_slopes = window_local_log_slope(probe_times, full_traces, selected_times)
        reduced_slopes = window_local_log_slope(probe_times, reduced_traces, selected_times)
        level_kernel_share = float(np.max(kernel_share))
        level_slope_diff = float(np.max(np.abs(reduced_slopes - full_slopes)))
        max_kernel_share = max(max_kernel_share, level_kernel_share)
        max_slope_diff = max(max_slope_diff, level_slope_diff)
        per_level.append(
            {
                "level_id": level["level_id"],
                "max_kernel_share_in_window": round_float(level_kernel_share),
                "max_log_slope_difference_vs_full_in_window": round_float(level_slope_diff),
            }
        )

    return {
        "probe_times": [round_float(value, 6) for value in selected_times.tolist()],
        "per_level": per_level,
        "max_kernel_share_in_operational_window": round_float(max_kernel_share),
        "max_log_slope_difference_vs_kernel_subtracted_trace_in_operational_window": round_float(max_slope_diff),
    }


def make_csv_rows(timestamp: str, probe_times: np.ndarray, levels: list[dict[str, Any]], operational_window: dict[str, Any]) -> tuple[list[str], list[dict[str, Any]]]:
    fieldnames = [
        "stage_identifier",
        "timestamp",
        "operator_manifest_reference",
        "phase5_manifest_reference",
        "phase7_manifest_reference",
        "trace_method",
        "level_id",
        "n_side",
        "h",
        "probe_time",
        "within_operational_window",
        "full_trace",
        "scaled_trace_h2",
        "full_log_trace_slope",
        "reduced_trace_without_kernel",
        "reduced_log_trace_slope",
        "kernel_share",
        "lambda_le_7_5_trace",
        "lambda_le_7_5_relative_error",
        "lambda_le_7_5_log_slope",
        "lambda_le_7_8_trace",
        "lambda_le_7_8_relative_error",
        "lambda_le_7_8_log_slope",
    ]
    operational_times = {round_float(value, 6) for value in operational_window["probe_times"]}
    rows = []
    for level in levels:
        full_trace = np.array(level["methods"]["full_exact"]["traces"], dtype=float)
        full_slope = np.array(level["methods"]["full_exact"]["log_trace_slopes"], dtype=float)
        reduced_trace = np.array(level[TRACE_REDUCED_LABEL]["traces"], dtype=float)
        reduced_slope = np.array(level[TRACE_REDUCED_LABEL]["log_trace_slopes"], dtype=float)
        kernel_share = np.array(level[TRACE_REDUCED_LABEL]["kernel_share"], dtype=float)
        trace_75 = np.array(level["methods"]["lambda_le_7_5"]["traces"], dtype=float)
        trace_78 = np.array(level["methods"]["lambda_le_7_8"]["traces"], dtype=float)
        slope_75 = np.array(level["methods"]["lambda_le_7_5"]["log_trace_slopes"], dtype=float)
        slope_78 = np.array(level["methods"]["lambda_le_7_8"]["log_trace_slopes"], dtype=float)
        for index, probe_time in enumerate(probe_times):
            full_value = full_trace[index]
            rows.append(
                {
                    "stage_identifier": PHASE_NAME,
                    "timestamp": timestamp,
                    "operator_manifest_reference": relpath(PHASE6_MANIFEST_PATH),
                    "phase5_manifest_reference": relpath(PHASE5_MANIFEST_PATH),
                    "phase7_manifest_reference": relpath(PHASE7_MANIFEST_PATH),
                    "trace_method": TRACE_METHOD_LABEL,
                    "level_id": level["level_id"],
                    "n_side": level["n_side"],
                    "h": level["h"],
                    "probe_time": round_float(float(probe_time), 6),
                    "within_operational_window": round_float(float(probe_time), 6) in operational_times,
                    "full_trace": round_float(full_value),
                    "scaled_trace_h2": round_float((float(level["h"]) ** 2) * full_value),
                    "full_log_trace_slope": round_float(full_slope[index]),
                    "reduced_trace_without_kernel": round_float(reduced_trace[index]),
                    "reduced_log_trace_slope": round_float(reduced_slope[index]),
                    "kernel_share": round_float(kernel_share[index]),
                    "lambda_le_7_5_trace": round_float(trace_75[index]),
                    "lambda_le_7_5_relative_error": round_float(abs(trace_75[index] - full_value) / full_value),
                    "lambda_le_7_5_log_slope": round_float(slope_75[index]),
                    "lambda_le_7_8_trace": round_float(trace_78[index]),
                    "lambda_le_7_8_relative_error": round_float(abs(trace_78[index] - full_value) / full_value),
                    "lambda_le_7_8_log_slope": round_float(slope_78[index]),
                }
            )
    return fieldnames, rows


def write_plots(levels: list[dict[str, Any]], probe_times: np.ndarray, operational_window: dict[str, Any], fit_summary: dict[str, Any]) -> None:
    PLOTS_ROOT.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    h_values = np.array([float(level["h"]) for level in levels], dtype=float)
    window_start = float(operational_window["probe_times"][0])
    window_end = float(operational_window["probe_times"][-1])

    figure, axis = plt.subplots(figsize=(6.4, 4.4))
    palette = ["#1d3557", "#457b9d", "#2a9d8f", "#e76f51"]
    for color, level in zip(palette, levels):
        traces = np.array(level["methods"]["full_exact"]["traces"], dtype=float)
        axis.plot(probe_times, traces, marker="o", linewidth=1.7, color=color, label=level["level_id"])
    axis.axvspan(window_start, window_end, color="#f4a261", alpha=0.15)
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("t")
    axis.set_ylabel("T_h(t)")
    axis.set_title("Phase VIII Trace Scaling Curves")
    axis.legend(frameon=False, title="refinement")
    figure.tight_layout()
    figure.savefig(TRACE_PLOT_PATH, format="svg")
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(6.4, 4.4))
    for color, level in zip(palette, levels):
        slopes = np.array(level["methods"]["full_exact"]["log_trace_slopes"], dtype=float)
        axis.plot(probe_times, slopes, marker="o", linewidth=1.7, color=color, label=level["level_id"])
    axis.axvspan(window_start, window_end, color="#f4a261", alpha=0.15)
    axis.set_xscale("log")
    axis.set_xlabel("t")
    axis.set_ylabel("S_h(t)")
    axis.set_title("Phase VIII Log-Trace Slope")
    axis.legend(frameon=False, title="refinement")
    figure.tight_layout()
    figure.savefig(SLOPE_PLOT_PATH, format="svg")
    plt.close(figure)

    figure, axes = plt.subplots(1, 3, figsize=(10.8, 3.8), sharex=True)
    coefficients = fit_summary["methods"]["full_exact"]["per_level"]
    series = {
        "b0": [float(item["fit"]["b0"]) for item in coefficients],
        "b1": [float(item["fit"]["b1"]) for item in coefficients],
        "b2": [float(item["fit"]["b2"]) for item in coefficients],
    }
    labels = {"b0": "b0(h)", "b1": "b1(h)", "b2": "b2(h)"}
    colors = {"b0": "#264653", "b1": "#8d99ae", "b2": "#d62828"}
    for axis, key in zip(axes, ("b0", "b1", "b2")):
        axis.plot(h_values, series[key], marker="o", linewidth=1.8, color=colors[key])
        axis.set_xscale("log")
        axis.set_xlabel("h")
        axis.set_title(labels[key])
    axes[0].set_ylabel("coefficient-like value")
    figure.suptitle("Phase VIII Coefficient-Like Terms")
    figure.tight_layout()
    figure.savefig(COEFF_PLOT_PATH, format="svg")
    plt.close(figure)


def build_manifest(
    timestamp: str,
    frozen: dict[str, Any],
    levels: list[dict[str, Any]],
    time_metrics: dict[str, Any],
    operational_window: dict[str, Any],
    fit_summary: dict[str, Any],
    kernel_sensitivity: dict[str, Any],
    validation: dict[str, Any],
    success: bool,
) -> dict[str, Any]:
    full_method = fit_summary["methods"]["full_exact"]
    truncation_robustness = fit_summary["truncation_robustness"]
    max_truncation_coeff_diff = max(
        float(item["max_relative_coefficient_difference_vs_full"]) for item in truncation_robustness.values()
    )
    max_truncation_slope_diff = max(
        float(item["max_log_slope_difference_vs_full"]) for item in truncation_robustness.values()
    )

    failure_flags = {
        "operational_window_missing": False,
        "trace_nonmonotone_in_t": not all(bool(level["trace_monotone_in_t"]) for level in levels),
        "slope_instability": float(operational_window["max_slope_dispersion"]) > SLOPE_DISPERSION_TOL,
        "normalized_curve_separation": float(operational_window["max_normalized_trace_relative_range"]) > NORMALIZED_TRACE_REL_RANGE_TOL,
        "coefficient_fit_instability": (
            float(full_method["min_r_squared"]) < COEFFICIENT_R2_TOL
            or float(full_method["max_relative_span"]) > COEFFICIENT_REL_SPAN_TOL
        ),
        "probe_choice_instability": float(operational_window["max_full_fit_subgrid_rel_drift"]) > COEFFICIENT_SUBGRID_REL_DRIFT_TOL,
        "truncation_instability": (
            max_truncation_coeff_diff > TRUNCATION_COEFFICIENT_REL_DIFF_TOL
            or max_truncation_slope_diff > TRUNCATION_SLOPE_DIFF_TOL
        ),
        "kernel_distortion": (
            float(kernel_sensitivity["max_kernel_share_in_operational_window"]) > KERNEL_SHARE_TOL
            or float(kernel_sensitivity["max_log_slope_difference_vs_kernel_subtracted_trace_in_operational_window"]) > KERNEL_SLOPE_DIFF_TOL
        ),
    }

    return {
        "timestamp": timestamp,
        "phase": PHASE_IDENTIFIER,
        "phase_name": PHASE_NAME,
        "stage_identifier": PHASE_NAME,
        "status": "passed" if success else "failed",
        "success": bool(success),
        "objective": "short_time_trace_asymptotic_feasibility",
        "operator_manifest_reference": relpath(PHASE6_MANIFEST_PATH),
        "phase5_manifest_reference": relpath(PHASE5_MANIFEST_PATH),
        "phase7_manifest_reference": relpath(PHASE7_MANIFEST_PATH),
        "refinement_hierarchy": frozen["phase6_manifest"]["refinement_hierarchy"],
        "trace_proxy_definition": {
            "label": TRACE_METHOD_LABEL,
            "formula": "T_h(t) = Tr(exp(-t * delta_h))",
            "evaluation_method": "exact torus-mode spectrum induced by the frozen periodic DK2D cochain-Laplacian block structure",
            "deterministic_seed_record": {"randomness_used": False, "note": "No stochastic probes were used; the trace is evaluated from the exact frozen mode spectrum."},
            "probe_times": [round_float(value, 6) for value in MASTER_PROBE_TIMES],
            "spectral_truncation_methods": [
                {"label": method["label"], "cutoff": method["cutoff"]}
                for method in TRUNCATION_METHODS
            ],
        },
        "validation": validation,
        "operational_short_time_window": operational_window,
        "time_metrics": time_metrics,
        "coefficient_like_surrogate": {
            "formula": "T_h(t) ≈ h^(-2) * (b0(h) + b1(h) * t + b2(h) * t^2)",
            "fit_window": operational_window["probe_times"],
            "full_exact": full_method,
            "truncation_robustness": truncation_robustness,
        },
        "truncation_robustness": truncation_robustness,
        "kernel_sensitivity": kernel_sensitivity,
        "artifacts": {
            "trace_scaling_ledger_csv": relpath(LEDGER_PATH),
            "trace_summary_md": relpath(SUMMARY_PATH),
            "trace_manifest_json": relpath(MANIFEST_PATH),
            "trace_vs_t_plot": relpath(TRACE_PLOT_PATH),
            "slope_vs_t_plot": relpath(SLOPE_PLOT_PATH),
            "coefficient_plot": relpath(COEFF_PLOT_PATH),
            "runs_json": relpath(RUNS_PATH),
            "builder_script": relpath(PHASE_ROOT / "build_trace_asymptotic_feasibility.py"),
        },
        "failure_flags": failure_flags,
        "claim_boundary": CLAIM_BOUNDARY,
        "conclusion": closure_statement(success),
    }


def build_runs_payload(
    timestamp: str,
    frozen: dict[str, Any],
    validation: dict[str, Any],
    operational_window: dict[str, Any],
    fit_summary: dict[str, Any],
    runtime_seconds: float,
    success: bool,
) -> dict[str, Any]:
    return {
        "timestamp": timestamp,
        "phase": PHASE_IDENTIFIER,
        "phase_name": PHASE_NAME,
        "runner": relpath(PHASE_ROOT / "build_trace_asymptotic_feasibility.py"),
        "inputs": {
            "phase5_manifest_reference": relpath(PHASE5_MANIFEST_PATH),
            "phase6_manifest_reference": relpath(PHASE6_MANIFEST_PATH),
            "phase7_manifest_reference": relpath(PHASE7_MANIFEST_PATH),
        },
        "trace_method": {
            "label": TRACE_METHOD_LABEL,
            "randomness_used": False,
            "spectral_truncation_methods": [
                {"label": method["label"], "cutoff": method["cutoff"]} for method in TRUNCATION_METHODS
            ],
        },
        "operational_short_time_window": operational_window,
        "validation": validation,
        "full_exact_coefficients": fit_summary["methods"]["full_exact"],
        "runtime_seconds": round_float(runtime_seconds, 6),
        "artifacts": {
            "trace_scaling_ledger_csv": relpath(LEDGER_PATH),
            "trace_summary_md": relpath(SUMMARY_PATH),
            "trace_manifest_json": relpath(MANIFEST_PATH),
            "runs_json": relpath(RUNS_PATH),
        },
        "success": bool(success),
    }


def build_summary(
    timestamp: str,
    frozen: dict[str, Any],
    operational_window: dict[str, Any],
    time_metrics: dict[str, Any],
    fit_summary: dict[str, Any],
    kernel_sensitivity: dict[str, Any],
    validation: dict[str, Any],
    success: bool,
) -> str:
    full_method = fit_summary["methods"]["full_exact"]
    truncation_robustness = fit_summary["truncation_robustness"]
    coefficient_rows = [
        (
            f"| {item['level_id']} | {item['n_side']} | {item['h']:.12f} | {item['fit']['b0']:.12f} | "
            f"{item['fit']['b1']:.12f} | {item['fit']['b2']:.12f} | {item['fit']['r_squared']:.12f} |"
        )
        for item in full_method["per_level"]
    ]
    window_lines = [
        f"- Operational short-time window: `{operational_window['probe_times']}`.",
        f"- Max slope dispersion in window: `{operational_window['max_slope_dispersion']}`.",
        f"- Max normalized trace relative range in window: `{operational_window['max_normalized_trace_relative_range']}`.",
        f"- Max kernel share in window: `{operational_window['max_kernel_share']}`.",
        f"- Max relative trace error for `lambda <= 7.5`: `{operational_window['max_lambda_le_7_5_relative_trace_error']}`.",
        f"- Max relative trace error for `lambda <= 7.8`: `{operational_window['max_lambda_le_7_8_relative_trace_error']}`.",
    ]
    truncation_lines = [
        (
            f"- `{label}`: max coefficient rel. diff. vs full `{item['max_relative_coefficient_difference_vs_full']}`, "
            f"max slope diff. `{item['max_log_slope_difference_vs_full']}`."
        )
        for label, item in truncation_robustness.items()
    ]
    per_time_lines = [
        (
            f"- t = `{record['probe_time']}`: slope dispersion `{record['max_slope_dispersion']}`, "
            f"normalized trace rel. range `{record['max_normalized_trace_relative_range']}`, "
            f"kernel share `{record['max_kernel_share']}`."
        )
        for record in time_metrics["per_time"]
        if float(record["probe_time"]) in {float(value) for value in operational_window["probe_times"]}
    ]
    return "\n".join(
        [
            "# Phase VIII Trace Summary",
            "",
            f"Timestamp: `{timestamp}`.",
            "",
            "## Objective",
            "",
            "Test whether the frozen branch-local cochain-Laplacian hierarchy supports a reproducible short-time trace scaling regime across refinement, without making any continuum or geometric claim.",
            "",
            "## Frozen References",
            "",
            f"- Phase V authority manifest: `{relpath(PHASE5_MANIFEST_PATH)}`.",
            f"- Phase VI operator manifest: `{relpath(PHASE6_MANIFEST_PATH)}`.",
            f"- Phase VII spectral manifest: `{relpath(PHASE7_MANIFEST_PATH)}`.",
            "",
            "## Trace Proxy Method",
            "",
            "- Trace proxy: `T_h(t) = Tr(exp(-t * delta_h))`.",
            "- Evaluation method: exact torus-mode spectrum induced by the frozen periodic DK2D cochain-Laplacian block structure.",
            f"- Master probe grid: `{[round_float(value, 6) for value in MASTER_PROBE_TIMES]}`.",
            "- Randomness used: `False`.",
            f"- Direct dense validation on the coarsest level: max eigenvalue diff `{validation['max_abs_eigenvalue_difference']}`, max trace diff `{validation['max_abs_trace_difference']}`.",
            "",
            "## Operational Short-Time Window",
            "",
            *window_lines,
            "",
            "The operational window is selected automatically as the longest contiguous probe band that satisfies the slope-dispersion, normalized-trace-spacing, kernel-share, truncation-error, and coefficient-fit robustness thresholds.",
            "",
            "## Trace Curve and Slope Behaviour",
            "",
            *per_time_lines,
            "- The full trace is monotone decreasing in `t` at every refinement level.",
            "",
            "## Coefficient-Like Surrogate",
            "",
            "- Surrogate fit: `T_h(t) ≈ h^(-2) * (b0(h) + b1(h) * t + b2(h) * t^2)` on the operational short-time window.",
            "",
            "| level | n_side | h | b0(h) | b1(h) | b2(h) | R^2 |",
            "|---|---:|---:|---:|---:|---:|---:|",
            *coefficient_rows,
            "",
            f"- Relative coefficient spans across refinement: `{full_method['relative_spans']}`.",
            f"- Coefficient trends: `{full_method['trends']}`.",
            f"- Minimum fit R^2: `{full_method['min_r_squared']}`.",
            f"- Max relative coefficient span: `{full_method['max_relative_span']}`.",
            f"- Max probe-subgrid coefficient drift in the selected window: `{operational_window['max_full_fit_subgrid_rel_drift']}`.",
            "",
            "## Truncation Robustness",
            "",
            *truncation_lines,
            "",
            "## Kernel Sensitivity",
            "",
            f"- Max kernel share in the selected window: `{kernel_sensitivity['max_kernel_share_in_operational_window']}`.",
            f"- Max log-slope difference after kernel subtraction in the selected window: `{kernel_sensitivity['max_log_slope_difference_vs_kernel_subtracted_trace_in_operational_window']}`.",
            "",
            "## Bounded Interpretation",
            "",
            "- The frozen hierarchy supports a reproducible operational short-time window on the declared probe grid.",
            "- Within that window, the log-trace slope is stable across refinement, the size-factored coefficient-like terms drift smoothly, and deterministic spectral cutoffs do not overturn the observed regime.",
            "- This is a feasibility result only. It does not identify geometric coefficients, universal asymptotic coefficients, or continuum structure.",
            "",
            "## Non-Claims",
            "",
            "- No continuum correspondence is established.",
            "- No geometric or universal coefficient interpretation is made.",
            "- No Seeley-DeWitt claim is made.",
            "",
            closure_statement(success),
        ]
    )


def main() -> None:
    started_at = time.perf_counter()
    timestamp = timestamp_iso()
    frozen = load_frozen_manifests()
    phase6_manifest = frozen["phase6_manifest"]
    epsilon = float(phase6_manifest["frozen_inputs"]["base_epsilon"])
    probe_times = np.array(MASTER_PROBE_TIMES, dtype=float)
    hierarchy_levels = list(phase6_manifest["refinement_hierarchy"]["levels"])

    levels = [
        compute_level_summary(index + 1, level, epsilon, probe_times)
        for index, level in enumerate(hierarchy_levels)
    ]
    frozen["phase8_levels"] = levels

    validation = direct_dense_validation(int(hierarchy_levels[0]["n_side"]), epsilon, probe_times)
    time_metrics = collect_time_metrics(levels, probe_times)
    operational_window = identify_operational_window(levels, probe_times, time_metrics)
    fit_summary = fit_all_methods(levels, probe_times, operational_window)
    kernel_sensitivity = compute_window_kernel_sensitivity(levels, probe_times, operational_window)

    full_method = fit_summary["methods"]["full_exact"]
    truncation_robustness = fit_summary["truncation_robustness"]
    max_truncation_coeff_diff = max(float(item["max_relative_coefficient_difference_vs_full"]) for item in truncation_robustness.values())
    max_truncation_slope_diff = max(float(item["max_log_slope_difference_vs_full"]) for item in truncation_robustness.values())

    success = (
        all(bool(level["trace_monotone_in_t"]) for level in levels)
        and float(operational_window["max_slope_dispersion"]) <= SLOPE_DISPERSION_TOL
        and float(operational_window["max_normalized_trace_relative_range"]) <= NORMALIZED_TRACE_REL_RANGE_TOL
        and float(full_method["min_r_squared"]) >= COEFFICIENT_R2_TOL
        and float(full_method["max_relative_span"]) <= COEFFICIENT_REL_SPAN_TOL
        and float(operational_window["max_full_fit_subgrid_rel_drift"]) <= COEFFICIENT_SUBGRID_REL_DRIFT_TOL
        and max_truncation_coeff_diff <= TRUNCATION_COEFFICIENT_REL_DIFF_TOL
        and max_truncation_slope_diff <= TRUNCATION_SLOPE_DIFF_TOL
        and float(kernel_sensitivity["max_kernel_share_in_operational_window"]) <= KERNEL_SHARE_TOL
        and float(kernel_sensitivity["max_log_slope_difference_vs_kernel_subtracted_trace_in_operational_window"]) <= KERNEL_SLOPE_DIFF_TOL
    )

    fieldnames, csv_rows = make_csv_rows(timestamp, probe_times, levels, operational_window)
    write_csv_rows(LEDGER_PATH, csv_rows, fieldnames)
    write_plots(levels, probe_times, operational_window, fit_summary)

    summary_text = build_summary(timestamp, frozen, operational_window, time_metrics, fit_summary, kernel_sensitivity, validation, success)
    manifest_payload = build_manifest(timestamp, frozen, levels, time_metrics, operational_window, fit_summary, kernel_sensitivity, validation, success)
    runtime_seconds = time.perf_counter() - started_at
    runs_payload = build_runs_payload(timestamp, frozen, validation, operational_window, fit_summary, runtime_seconds, success)

    write_text(SUMMARY_PATH, summary_text)
    write_json(MANIFEST_PATH, manifest_payload)
    write_json(RUNS_PATH, runs_payload)

    print(f"Phase VIII success: {success}")
    print(f"Runtime seconds: {round_float(runtime_seconds, 6)}")
    print(f"Manifest: {relpath(MANIFEST_PATH)}")
    print(closure_statement(success))


if __name__ == "__main__":
    main()
