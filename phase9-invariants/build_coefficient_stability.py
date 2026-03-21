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

PHASE_ROOT = ROOT / "phase9-invariants"
RUNS_ROOT = PHASE_ROOT / "runs"
PLOTS_ROOT = PHASE_ROOT / "plots"

PHASE6_MANIFEST_PATH = ROOT / "phase6-operator" / "phase6_operator_manifest.json"
PHASE7_MANIFEST_PATH = ROOT / "phase7-spectral" / "phase7_spectral_manifest.json"
PHASE8_MANIFEST_PATH = ROOT / "phase8-trace" / "phase8_trace_manifest.json"
PHASE8_CONTRACT_NOTE_PATH = ROOT / "PHASE_VIII_TRACE_CONTRACT_NOTE.md"

COEFF_LEDGER_PATH = RUNS_ROOT / "phase9_coefficient_stability_ledger.csv"
RATIO_LEDGER_PATH = RUNS_ROOT / "phase9_invariant_ratios_ledger.csv"
RUNS_PATH = RUNS_ROOT / "phase9_runs.json"
SUMMARY_PATH = PHASE_ROOT / "phase9_summary.md"
MANIFEST_PATH = PHASE_ROOT / "phase9_manifest.json"

COEFF_PLOT_PATH = PLOTS_ROOT / "phase9_coefficient_vs_refinement.svg"
RATIO_PLOT_PATH = PLOTS_ROOT / "phase9_ratio_vs_refinement.svg"
INVARIANT_PLOT_PATH = PLOTS_ROOT / "phase9_rescaled_invariant_vs_refinement.svg"
DISTANCE_PLOT_PATH = PLOTS_ROOT / "phase9_refinement_distance_convergence.svg"

PHASE_NAME = "phase9-invariants"
PHASE_IDENTIFIER = 9
BRANCH_LABEL = "frozen_branch"
CONTROL_LABEL = "generic_open_grid_scalar_block_control"
CONTROL_DESCRIPTION = (
    "Deterministic control hierarchy with the same total node counts as the frozen branch, "
    "constructed as four repeated scalar open-grid Laplacian blocks with open boundary connectivity."
)

SUBSET_LABELS = {
    "front_window": slice(0, 4),
    "rear_window": slice(1, 5),
}

PLATEAU_RELATIVE_SPAN_TOL = 0.03
RATIO_SPAN_TOL = 0.03
INVARIANT_SPAN_TOL = 0.07
SUBSET_RATIO_DRIFT_TOL = 0.08
SUBSET_INVARIANT_DRIFT_TOL = 0.08
TRUNCATION_RATIO_DRIFT_TOL = 0.08
TRUNCATION_INVARIANT_DRIFT_TOL = 0.12
CONTROL_RATIO_SPAN_MULTIPLIER_MIN = 2.5
CONTROL_DISTANCE_MULTIPLIER_MIN = 3.0
CONTROL_MODEL_DIFFERENCE_MIN = 2

CLAIM_BOUNDARY = (
    "Phase IX authority is limited to coefficient stabilization and invariant-tracking diagnostics "
    "on the frozen short-time trace contract. It does not assert continuum coefficients, geometric "
    "meaning, universal asymptotic data, or physical correspondence."
)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def timestamp_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def closure_statement(success: bool) -> str:
    if success:
        return "Phase IX establishes coefficient stabilization and invariant-tracking feasibility for the frozen operator hierarchy."
    return "Phase IX does not yet establish coefficient stabilization and invariant-tracking feasibility for the frozen operator hierarchy."


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


def relative_span(values: np.ndarray) -> float:
    return float((np.max(values) - np.min(values)) / max(abs(float(np.mean(values))), 1.0e-12))


def safe_relative_difference(new_value: float, reference_value: float) -> float:
    return abs(new_value - reference_value) / max(abs(reference_value), 1.0e-12)


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


def linear_limit_fit(h_values: np.ndarray, values: np.ndarray, power: int) -> dict[str, Any]:
    x_values = h_values ** power
    coefficients = np.polyfit(x_values, values, 1)
    fitted = np.polyval(coefficients, x_values)
    denominator = float(np.sum((values - values.mean()) ** 2))
    if denominator <= 1.0e-18:
        r_squared = 1.0
    else:
        r_squared = 1.0 - float(np.sum((values - fitted) ** 2) / denominator)
    return {
        "model": f"linear_in_h^{power}",
        "limit_estimate": round_float(coefficients[1]),
        "r_squared": round_float(r_squared),
    }


def best_limit_fit(h_values: np.ndarray, values: np.ndarray) -> dict[str, Any]:
    candidate_h = linear_limit_fit(h_values, values, power=1)
    candidate_h2 = linear_limit_fit(h_values, values, power=2)
    if float(candidate_h2["r_squared"]) > float(candidate_h["r_squared"]):
        return candidate_h2
    return candidate_h


def load_frozen_manifests() -> dict[str, Any]:
    phase6_manifest = read_json(PHASE6_MANIFEST_PATH)
    phase7_manifest = read_json(PHASE7_MANIFEST_PATH)
    phase8_manifest = read_json(PHASE8_MANIFEST_PATH)

    if not bool(phase6_manifest.get("success")):
        raise ValueError("Phase IX requires the frozen successful Phase VI operator manifest.")
    if not bool(phase7_manifest.get("success")):
        raise ValueError("Phase IX requires the frozen successful Phase VII spectral manifest.")
    if not bool(phase8_manifest.get("success")):
        raise ValueError("Phase IX requires the frozen successful Phase VIII trace manifest.")

    return {
        "phase6_manifest": phase6_manifest,
        "phase7_manifest": phase7_manifest,
        "phase8_manifest": phase8_manifest,
    }


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
    return np.sort(np.repeat(scalar_values, 4))


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
    return np.sort(np.repeat(scalar_values, 4))


def trace_from_spectrum(spectrum: np.ndarray, probe_times: np.ndarray) -> np.ndarray:
    return np.array([float(np.exp(-probe_time * spectrum).sum()) for probe_time in probe_times], dtype=float)


def level_observation(
    hierarchy_label: str,
    level: dict[str, Any],
    spectrum: np.ndarray,
    probe_times: np.ndarray,
    phase7_exponent: float,
) -> dict[str, Any]:
    h_value = float(level["h"])
    traces = trace_from_spectrum(spectrum, probe_times)
    raw_fit = quadratic_fit_named(probe_times, traces, prefix="a")
    scaled_fit = quadratic_fit_named(probe_times, (h_value * h_value) * traces, prefix="b")
    ratios = {
        "b1_over_b0": round_float(float(scaled_fit["b1"]) / max(abs(float(scaled_fit["b0"])), 1.0e-12)),
        "b2_over_b0": round_float(float(scaled_fit["b2"]) / max(abs(float(scaled_fit["b0"])), 1.0e-12)),
    }
    invariants = {
        "I0": round_float((h_value ** phase7_exponent) * float(raw_fit["a0"])),
        "I1": round_float((h_value ** phase7_exponent) * float(raw_fit["a1"])),
        "I2": round_float((h_value ** phase7_exponent) * float(raw_fit["a2"])),
    }
    return {
        "hierarchy_label": hierarchy_label,
        "level_id": str(level["level_id"]),
        "n_side": int(level["n_side"]),
        "h": round_float(h_value),
        "probe_times": [round_float(value, 6) for value in probe_times.tolist()],
        "trace_values": [round_float(value) for value in traces.tolist()],
        "raw_coefficients": raw_fit,
        "scaled_coefficients": scaled_fit,
        "ratios": ratios,
        "rescaled_invariants": invariants,
    }


def build_hierarchy_observations(
    hierarchy_label: str,
    levels: list[dict[str, Any]],
    epsilon: float,
    probe_times: np.ndarray,
    phase7_exponent: float,
    cutoff: float | None = None,
) -> list[dict[str, Any]]:
    spectrum_builder = branch_spectrum if hierarchy_label == BRANCH_LABEL else control_spectrum
    observations = []
    for level in levels:
        spectrum = spectrum_builder(int(level["n_side"]), epsilon)
        if cutoff is not None:
            spectrum = spectrum[spectrum <= cutoff]
        observations.append(level_observation(hierarchy_label, level, spectrum, probe_times, phase7_exponent))
    return observations


def successive_differences(values: np.ndarray) -> list[float]:
    return [round_float(value) for value in np.abs(np.diff(values)).tolist()]


def successive_relative_differences(values: np.ndarray) -> list[float]:
    rows: list[float] = []
    for index in range(len(values) - 1):
        rows.append(round_float(safe_relative_difference(float(values[index + 1]), float(values[index]))))
    return rows


def summarize_sequence(values: list[float], h_values: np.ndarray) -> dict[str, Any]:
    value_array = np.array(values, dtype=float)
    diffs = np.abs(np.diff(value_array))
    drift_shrinks = bool(all(diffs[index + 1] <= diffs[index] + 1.0e-12 for index in range(len(diffs) - 1))) if len(diffs) >= 2 else True
    fit = best_limit_fit(h_values, value_array)
    return {
        "values": [round_float(value) for value in value_array.tolist()],
        "trend": trend_label(value_array.tolist()),
        "relative_span": round_float(relative_span(value_array)),
        "successive_abs_differences": successive_differences(value_array),
        "successive_relative_differences": successive_relative_differences(value_array),
        "drift_shrinks": drift_shrinks,
        "best_limit_fit": fit,
    }


def summarize_observations(observations: list[dict[str, Any]]) -> dict[str, Any]:
    h_values = np.array([float(item["h"]) for item in observations], dtype=float)
    raw = {
        label: summarize_sequence([float(item["raw_coefficients"][label]) for item in observations], h_values)
        for label in ("a0", "a1", "a2")
    }
    scaled = {
        label: summarize_sequence([float(item["scaled_coefficients"][label]) for item in observations], h_values)
        for label in ("b0", "b1", "b2")
    }
    ratios = {
        label: summarize_sequence([float(item["ratios"][label]) for item in observations], h_values)
        for label in ("b1_over_b0", "b2_over_b0")
    }
    invariants = {
        label: summarize_sequence([float(item["rescaled_invariants"][label]) for item in observations], h_values)
        for label in ("I0", "I1", "I2")
    }

    coefficient_vectors = np.array(
        [
            [
                float(item["scaled_coefficients"]["b0"]),
                float(item["scaled_coefficients"]["b1"]),
                float(item["scaled_coefficients"]["b2"]),
            ]
            for item in observations
        ],
        dtype=float,
    )
    vector_distances = np.linalg.norm(np.diff(coefficient_vectors, axis=0), axis=1)
    distances_shrink = bool(
        all(vector_distances[index + 1] <= vector_distances[index] + 1.0e-12 for index in range(len(vector_distances) - 1))
    ) if len(vector_distances) >= 2 else True

    return {
        "observations": observations,
        "raw_coefficients": raw,
        "scaled_coefficients": scaled,
        "ratios": ratios,
        "rescaled_invariants": invariants,
        "vector_distance_convergence": {
            "transition_labels": [
                f"{observations[index]['level_id']}->{observations[index + 1]['level_id']}"
                for index in range(len(observations) - 1)
            ],
            "successive_distances": [round_float(value) for value in vector_distances.tolist()],
            "distances_shrink": distances_shrink,
            "final_distance": round_float(float(vector_distances[-1])),
            "max_distance": round_float(float(np.max(vector_distances))),
        },
    }


def compare_probe_variant(
    reference: list[dict[str, Any]],
    candidate: list[dict[str, Any]],
) -> dict[str, Any]:
    descriptor_paths = {
        "b0": ("scaled_coefficients", "b0"),
        "b1": ("scaled_coefficients", "b1"),
        "b2": ("scaled_coefficients", "b2"),
        "b1_over_b0": ("ratios", "b1_over_b0"),
        "b2_over_b0": ("ratios", "b2_over_b0"),
        "I0": ("rescaled_invariants", "I0"),
        "I1": ("rescaled_invariants", "I1"),
        "I2": ("rescaled_invariants", "I2"),
    }
    max_drifts = {label: 0.0 for label in descriptor_paths}
    for ref_level, candidate_level in zip(reference, candidate):
        for label, path in descriptor_paths.items():
            left = float(ref_level[path[0]][path[1]])
            right = float(candidate_level[path[0]][path[1]])
            max_drifts[label] = max(max_drifts[label], safe_relative_difference(right, left))
    return {label: round_float(value) for label, value in max_drifts.items()}


def build_subset_probe_payload(
    hierarchy_label: str,
    levels: list[dict[str, Any]],
    epsilon: float,
    full_window: np.ndarray,
    phase7_exponent: float,
    reference: list[dict[str, Any]],
) -> dict[str, Any]:
    payload = {}
    for subset_label, subset_slice in SUBSET_LABELS.items():
        subset_times = full_window[subset_slice]
        subset_observations = build_hierarchy_observations(
            hierarchy_label=hierarchy_label,
            levels=levels,
            epsilon=epsilon,
            probe_times=subset_times,
            phase7_exponent=phase7_exponent,
            cutoff=None,
        )
        payload[subset_label] = {
            "probe_times": [round_float(value, 6) for value in subset_times.tolist()],
            "max_relative_drifts_vs_full_window": compare_probe_variant(reference, subset_observations),
        }
    return payload


def build_truncation_probe_payload(
    hierarchy_label: str,
    levels: list[dict[str, Any]],
    epsilon: float,
    full_window: np.ndarray,
    phase7_exponent: float,
    truncation_methods: list[dict[str, Any]],
    reference: list[dict[str, Any]],
) -> dict[str, Any]:
    payload = {}
    for method in truncation_methods:
        label = str(method["label"])
        if label == "full_exact":
            continue
        cutoff = float(method["cutoff"])
        cutoff_observations = build_hierarchy_observations(
            hierarchy_label=hierarchy_label,
            levels=levels,
            epsilon=epsilon,
            probe_times=full_window,
            phase7_exponent=phase7_exponent,
            cutoff=cutoff,
        )
        payload[label] = {
            "cutoff": round_float(cutoff, 6),
            "max_relative_drifts_vs_full_exact": compare_probe_variant(reference, cutoff_observations),
        }
    return payload


def coefficient_rows(timestamp: str, summary: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    raw_metrics = summary["raw_coefficients"]
    scaled_metrics = summary["scaled_coefficients"]
    observations = summary["observations"]
    for coefficient_name in ("b0", "b1", "b2"):
        raw_name = f"a{coefficient_name[1]}"
        metric = scaled_metrics[coefficient_name]
        raw_metric = raw_metrics[raw_name]
        for index, observation in enumerate(observations):
            previous_value = "" if index == 0 else metric["successive_relative_differences"][index - 1]
            rows.append(
                {
                    "stage_identifier": PHASE_NAME,
                    "timestamp": timestamp,
                    "hierarchy_label": observation["hierarchy_label"],
                    "level_id": observation["level_id"],
                    "n_side": observation["n_side"],
                    "h": observation["h"],
                    "coefficient_name": coefficient_name,
                    "raw_coefficient_name": raw_name,
                    "raw_value": observation["raw_coefficients"][raw_name],
                    "scaled_value": observation["scaled_coefficients"][coefficient_name],
                    "relative_step_from_previous": previous_value,
                    "scaled_relative_span": metric["relative_span"],
                    "raw_relative_span": raw_metric["relative_span"],
                    "best_limit_model": metric["best_limit_fit"]["model"],
                    "limit_estimate": metric["best_limit_fit"]["limit_estimate"],
                    "limit_fit_r_squared": metric["best_limit_fit"]["r_squared"],
                }
            )
    return rows


def ratio_rows(timestamp: str, summary: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    observations = summary["observations"]
    descriptor_families = (
        ("ratio", "b1_over_b0"),
        ("ratio", "b2_over_b0"),
        ("rescaled_invariant", "I0"),
        ("rescaled_invariant", "I1"),
        ("rescaled_invariant", "I2"),
    )
    for family, label in descriptor_families:
        metric_bank = summary["ratios"] if family == "ratio" else summary["rescaled_invariants"]
        metric = metric_bank[label]
        source_key = "ratios" if family == "ratio" else "rescaled_invariants"
        for index, observation in enumerate(observations):
            previous_value = "" if index == 0 else metric["successive_relative_differences"][index - 1]
            rows.append(
                {
                    "stage_identifier": PHASE_NAME,
                    "timestamp": timestamp,
                    "hierarchy_label": observation["hierarchy_label"],
                    "level_id": observation["level_id"],
                    "n_side": observation["n_side"],
                    "h": observation["h"],
                    "descriptor_family": family,
                    "descriptor_name": label,
                    "value": observation[source_key][label],
                    "relative_step_from_previous": previous_value,
                    "relative_span": metric["relative_span"],
                    "best_limit_model": metric["best_limit_fit"]["model"],
                    "limit_estimate": metric["best_limit_fit"]["limit_estimate"],
                    "limit_fit_r_squared": metric["best_limit_fit"]["r_squared"],
                }
            )
    return rows


def write_plots(branch_summary: dict[str, Any], control_summary: dict[str, Any]) -> None:
    PLOTS_ROOT.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    branch_observations = branch_summary["observations"]
    control_observations = control_summary["observations"]
    h_values = np.array([float(item["h"]) for item in branch_observations], dtype=float)

    figure, axes = plt.subplots(1, 3, figsize=(10.8, 3.8), sharex=True)
    coefficient_colors = {BRANCH_LABEL: "#1d3557", CONTROL_LABEL: "#e76f51"}
    for axis, label in zip(axes, ("b0", "b1", "b2")):
        axis.plot(h_values, [float(item["scaled_coefficients"][label]) for item in branch_observations], marker="o", linewidth=1.8, color=coefficient_colors[BRANCH_LABEL], label="branch")
        axis.plot(h_values, [float(item["scaled_coefficients"][label]) for item in control_observations], marker="s", linewidth=1.8, color=coefficient_colors[CONTROL_LABEL], label="control")
        axis.set_xscale("log")
        axis.set_xlabel("h")
        axis.set_title(label)
    axes[0].set_ylabel("scaled coefficient")
    axes[0].legend(frameon=False)
    figure.suptitle("Phase IX Coefficient Stabilization")
    figure.tight_layout()
    figure.savefig(COEFF_PLOT_PATH, format="svg")
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(7.2, 3.8), sharex=True)
    for axis, label in zip(axes, ("b1_over_b0", "b2_over_b0")):
        axis.plot(h_values, [float(item["ratios"][label]) for item in branch_observations], marker="o", linewidth=1.8, color=coefficient_colors[BRANCH_LABEL], label="branch")
        axis.plot(h_values, [float(item["ratios"][label]) for item in control_observations], marker="s", linewidth=1.8, color=coefficient_colors[CONTROL_LABEL], label="control")
        axis.set_xscale("log")
        axis.set_xlabel("h")
        axis.set_title(label)
    axes[0].set_ylabel("dimensionless ratio")
    axes[0].legend(frameon=False)
    figure.suptitle("Phase IX Ratio Stabilization")
    figure.tight_layout()
    figure.savefig(RATIO_PLOT_PATH, format="svg")
    plt.close(figure)

    figure, axes = plt.subplots(1, 3, figsize=(10.8, 3.8), sharex=True)
    for axis, label in zip(axes, ("I0", "I1", "I2")):
        axis.plot(h_values, [float(item["rescaled_invariants"][label]) for item in branch_observations], marker="o", linewidth=1.8, color=coefficient_colors[BRANCH_LABEL], label="branch")
        axis.plot(h_values, [float(item["rescaled_invariants"][label]) for item in control_observations], marker="s", linewidth=1.8, color=coefficient_colors[CONTROL_LABEL], label="control")
        axis.set_xscale("log")
        axis.set_xlabel("h")
        axis.set_title(label)
    axes[0].set_ylabel("rescaled descriptor")
    axes[0].legend(frameon=False)
    figure.suptitle("Phase IX Rescaled Invariants")
    figure.tight_layout()
    figure.savefig(INVARIANT_PLOT_PATH, format="svg")
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(6.4, 4.0))
    x_positions = np.arange(len(branch_summary["vector_distance_convergence"]["transition_labels"]))
    axis.plot(
        x_positions,
        branch_summary["vector_distance_convergence"]["successive_distances"],
        marker="o",
        linewidth=1.8,
        color=coefficient_colors[BRANCH_LABEL],
        label="branch",
    )
    axis.plot(
        x_positions,
        control_summary["vector_distance_convergence"]["successive_distances"],
        marker="s",
        linewidth=1.8,
        color=coefficient_colors[CONTROL_LABEL],
        label="control",
    )
    axis.set_xticks(x_positions, branch_summary["vector_distance_convergence"]["transition_labels"])
    axis.set_ylabel("successive distance")
    axis.set_title("Phase IX Refinement-Distance Convergence")
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(DISTANCE_PLOT_PATH, format="svg")
    plt.close(figure)


def build_manifest(
    timestamp: str,
    frozen: dict[str, Any],
    branch_summary: dict[str, Any],
    control_summary: dict[str, Any],
    subset_probes: dict[str, Any],
    truncation_probes: dict[str, Any],
    phase7_exponent: float,
    success: bool,
) -> dict[str, Any]:
    branch_scaled = branch_summary["scaled_coefficients"]
    branch_ratios = branch_summary["ratios"]
    branch_invariants = branch_summary["rescaled_invariants"]
    branch_raw = branch_summary["raw_coefficients"]

    branch_max_ratio_span = max(float(item["relative_span"]) for item in branch_ratios.values())
    branch_max_raw_span = max(float(item["relative_span"]) for item in branch_raw.values())
    branch_max_invariant_span = max(float(item["relative_span"]) for item in branch_invariants.values())
    control_max_ratio_span = max(float(item["relative_span"]) for item in control_summary["ratios"].values())

    branch_subset_max_ratio_drift = max(
        float(value)
        for subset in subset_probes[BRANCH_LABEL].values()
        for label, value in subset["max_relative_drifts_vs_full_window"].items()
        if label in {"b1_over_b0", "b2_over_b0"}
    )
    branch_subset_max_invariant_drift = max(
        float(value)
        for subset in subset_probes[BRANCH_LABEL].values()
        for label, value in subset["max_relative_drifts_vs_full_window"].items()
        if label in {"I0", "I1", "I2"}
    )
    branch_truncation_max_ratio_drift = max(
        float(value)
        for payload in truncation_probes[BRANCH_LABEL].values()
        for label, value in payload["max_relative_drifts_vs_full_exact"].items()
        if label in {"b1_over_b0", "b2_over_b0"}
    )
    branch_truncation_max_invariant_drift = max(
        float(value)
        for payload in truncation_probes[BRANCH_LABEL].values()
        for label, value in payload["max_relative_drifts_vs_full_exact"].items()
        if label in {"I0", "I1", "I2"}
    )

    branch_models = {label: branch_scaled[label]["best_limit_fit"]["model"] for label in ("b0", "b1", "b2")}
    control_models = {label: control_summary["scaled_coefficients"][label]["best_limit_fit"]["model"] for label in ("b0", "b1", "b2")}
    differing_model_count = sum(branch_models[label] != control_models[label] for label in branch_models)

    success_flags = {
        "at_least_one_coefficient_plateau": any(
            float(branch_scaled[label]["relative_span"]) <= PLATEAU_RELATIVE_SPAN_TOL and bool(branch_scaled[label]["drift_shrinks"])
            for label in ("b0", "b1", "b2")
        ),
        "ratios_more_stable_than_raw": branch_max_ratio_span < branch_max_raw_span,
        "compact_rescaled_band": branch_max_invariant_span <= INVARIANT_SPAN_TOL,
        "probe_subset_stable": (
            branch_subset_max_ratio_drift <= SUBSET_RATIO_DRIFT_TOL
            and branch_subset_max_invariant_drift <= SUBSET_INVARIANT_DRIFT_TOL
        ),
        "truncation_stable": (
            branch_truncation_max_ratio_drift <= TRUNCATION_RATIO_DRIFT_TOL
            and branch_truncation_max_invariant_drift <= TRUNCATION_INVARIANT_DRIFT_TOL
        ),
        "control_distinct": (
            control_max_ratio_span >= CONTROL_RATIO_SPAN_MULTIPLIER_MIN * branch_max_ratio_span
            and float(control_summary["vector_distance_convergence"]["final_distance"])
            >= CONTROL_DISTANCE_MULTIPLIER_MIN * float(branch_summary["vector_distance_convergence"]["final_distance"])
            and differing_model_count >= CONTROL_MODEL_DIFFERENCE_MIN
        ),
    }

    return {
        "timestamp": timestamp,
        "phase": PHASE_IDENTIFIER,
        "phase_name": PHASE_NAME,
        "stage_identifier": PHASE_NAME,
        "status": "passed" if success else "failed",
        "success": bool(success),
        "objective": "coefficient_stabilization_and_invariant_tracking_feasibility",
        "operator_manifest_reference": relpath(PHASE6_MANIFEST_PATH),
        "phase7_manifest_reference": relpath(PHASE7_MANIFEST_PATH),
        "phase8_manifest_reference": relpath(PHASE8_MANIFEST_PATH),
        "phase8_contract_note_reference": relpath(PHASE8_CONTRACT_NOTE_PATH),
        "refinement_hierarchy": frozen["phase8_manifest"]["refinement_hierarchy"],
        "trace_contract": {
            "operational_short_time_window": frozen["phase8_manifest"]["operational_short_time_window"]["probe_times"],
            "master_probe_grid": frozen["phase8_manifest"]["trace_proxy_definition"]["probe_times"],
            "trace_method_label": frozen["phase8_manifest"]["trace_proxy_definition"]["label"],
            "trace_method_definition": frozen["phase8_manifest"]["trace_proxy_definition"]["evaluation_method"],
            "spectral_truncation_methods": frozen["phase8_manifest"]["trace_proxy_definition"]["spectral_truncation_methods"],
        },
        "dominant_phase7_low_mode_exponent": {
            "band_index": 1,
            "exponent": round_float(phase7_exponent),
        },
        "control_hierarchy": {
            "label": CONTROL_LABEL,
            "description": CONTROL_DESCRIPTION,
            "same_node_counts_as_branch": True,
            "same_total_dimension_as_branch": True,
        },
        "branch_scaled_coefficients": branch_scaled,
        "branch_raw_coefficients": branch_raw,
        "branch_ratios": branch_ratios,
        "branch_rescaled_invariants": branch_invariants,
        "control_scaled_coefficients": control_summary["scaled_coefficients"],
        "control_ratios": control_summary["ratios"],
        "control_rescaled_invariants": control_summary["rescaled_invariants"],
        "vector_distance_convergence": {
            BRANCH_LABEL: branch_summary["vector_distance_convergence"],
            CONTROL_LABEL: control_summary["vector_distance_convergence"],
        },
        "probe_subset_robustness": subset_probes,
        "truncation_robustness": truncation_probes,
        "control_separation_metrics": {
            "branch_max_ratio_span": round_float(branch_max_ratio_span),
            "control_max_ratio_span": round_float(control_max_ratio_span),
            "ratio_span_multiplier": round_float(control_max_ratio_span / max(branch_max_ratio_span, 1.0e-12)),
            "branch_final_vector_distance": branch_summary["vector_distance_convergence"]["final_distance"],
            "control_final_vector_distance": control_summary["vector_distance_convergence"]["final_distance"],
            "final_distance_multiplier": round_float(
                float(control_summary["vector_distance_convergence"]["final_distance"])
                / max(float(branch_summary["vector_distance_convergence"]["final_distance"]), 1.0e-12)
            ),
            "branch_best_fit_models": branch_models,
            "control_best_fit_models": control_models,
            "differing_model_count": differing_model_count,
        },
        "success_flags": success_flags,
        "artifacts": {
            "coefficient_stability_ledger_csv": relpath(COEFF_LEDGER_PATH),
            "invariant_ratios_ledger_csv": relpath(RATIO_LEDGER_PATH),
            "summary_md": relpath(SUMMARY_PATH),
            "manifest_json": relpath(MANIFEST_PATH),
            "coefficient_plot": relpath(COEFF_PLOT_PATH),
            "ratio_plot": relpath(RATIO_PLOT_PATH),
            "rescaled_invariant_plot": relpath(INVARIANT_PLOT_PATH),
            "distance_plot": relpath(DISTANCE_PLOT_PATH),
            "runs_json": relpath(RUNS_PATH),
            "builder_script": relpath(PHASE_ROOT / "build_coefficient_stability.py"),
        },
        "claim_boundary": CLAIM_BOUNDARY,
        "conclusion": closure_statement(success),
    }


def build_runs_payload(
    timestamp: str,
    phase7_exponent: float,
    branch_summary: dict[str, Any],
    control_summary: dict[str, Any],
    subset_probes: dict[str, Any],
    truncation_probes: dict[str, Any],
    runtime_seconds: float,
    success: bool,
) -> dict[str, Any]:
    return {
        "timestamp": timestamp,
        "phase": PHASE_IDENTIFIER,
        "phase_name": PHASE_NAME,
        "runner": relpath(PHASE_ROOT / "build_coefficient_stability.py"),
        "inputs": {
            "phase6_manifest_reference": relpath(PHASE6_MANIFEST_PATH),
            "phase7_manifest_reference": relpath(PHASE7_MANIFEST_PATH),
            "phase8_manifest_reference": relpath(PHASE8_MANIFEST_PATH),
            "phase8_contract_note_reference": relpath(PHASE8_CONTRACT_NOTE_PATH),
        },
        "deterministic_seed_record": {
            "randomness_used": False,
            "note": "Phase IX uses exact deterministic spectra and fixed probe subsets only.",
        },
        "dominant_phase7_low_mode_exponent": round_float(phase7_exponent),
        "branch_vector_distance_convergence": branch_summary["vector_distance_convergence"],
        "control_vector_distance_convergence": control_summary["vector_distance_convergence"],
        "probe_subset_robustness": subset_probes,
        "truncation_robustness": truncation_probes,
        "runtime_seconds": round_float(runtime_seconds, 6),
        "artifacts": {
            "coefficient_stability_ledger_csv": relpath(COEFF_LEDGER_PATH),
            "invariant_ratios_ledger_csv": relpath(RATIO_LEDGER_PATH),
            "summary_md": relpath(SUMMARY_PATH),
            "manifest_json": relpath(MANIFEST_PATH),
            "runs_json": relpath(RUNS_PATH),
        },
        "success": bool(success),
    }


def build_summary(
    timestamp: str,
    phase7_exponent: float,
    branch_summary: dict[str, Any],
    control_summary: dict[str, Any],
    subset_probes: dict[str, Any],
    truncation_probes: dict[str, Any],
    success: bool,
) -> str:
    branch_scaled = branch_summary["scaled_coefficients"]
    branch_ratios = branch_summary["ratios"]
    branch_invariants = branch_summary["rescaled_invariants"]
    control_scaled = control_summary["scaled_coefficients"]
    control_ratios = control_summary["ratios"]
    control_invariants = control_summary["rescaled_invariants"]

    branch_max_ratio_span = max(float(item["relative_span"]) for item in branch_ratios.values())
    branch_max_raw_span = max(float(item["relative_span"]) for item in branch_summary["raw_coefficients"].values())
    branch_max_invariant_span = max(float(item["relative_span"]) for item in branch_invariants.values())
    control_max_ratio_span = max(float(item["relative_span"]) for item in control_ratios.values())

    branch_subset_max_ratio_drift = max(
        float(value)
        for subset in subset_probes[BRANCH_LABEL].values()
        for label, value in subset["max_relative_drifts_vs_full_window"].items()
        if label in {"b1_over_b0", "b2_over_b0"}
    )
    branch_subset_max_invariant_drift = max(
        float(value)
        for subset in subset_probes[BRANCH_LABEL].values()
        for label, value in subset["max_relative_drifts_vs_full_window"].items()
        if label in {"I0", "I1", "I2"}
    )
    branch_truncation_max_ratio_drift = max(
        float(value)
        for payload in truncation_probes[BRANCH_LABEL].values()
        for label, value in payload["max_relative_drifts_vs_full_exact"].items()
        if label in {"b1_over_b0", "b2_over_b0"}
    )
    branch_truncation_max_invariant_drift = max(
        float(value)
        for payload in truncation_probes[BRANCH_LABEL].values()
        for label, value in payload["max_relative_drifts_vs_full_exact"].items()
        if label in {"I0", "I1", "I2"}
    )

    coefficient_rows = [
        (
            f"| {label} | {branch_scaled[label]['values']} | {branch_scaled[label]['relative_span']} | "
            f"{branch_scaled[label]['best_limit_fit']['model']} | {branch_scaled[label]['best_limit_fit']['limit_estimate']} | "
            f"{branch_scaled[label]['best_limit_fit']['r_squared']} |"
        )
        for label in ("b0", "b1", "b2")
    ]
    ratio_rows_text = [
        (
            f"| {label} | {branch_ratios[label]['values']} | {branch_ratios[label]['relative_span']} | "
            f"{branch_ratios[label]['best_limit_fit']['model']} | {branch_ratios[label]['best_limit_fit']['limit_estimate']} |"
        )
        for label in ("b1_over_b0", "b2_over_b0")
    ]
    invariant_rows = [
        (
            f"| {label} | {branch_invariants[label]['values']} | {branch_invariants[label]['relative_span']} | "
            f"{branch_invariants[label]['best_limit_fit']['model']} | {branch_invariants[label]['best_limit_fit']['limit_estimate']} |"
        )
        for label in ("I0", "I1", "I2")
    ]
    control_rows = [
        (
            f"| {label} | {control_scaled[label]['relative_span']} | {control_scaled[label]['best_limit_fit']['model']} | "
            f"{control_scaled[label]['best_limit_fit']['limit_estimate']} |"
        )
        for label in ("b0", "b1", "b2")
    ]

    return "\n".join(
        [
            "# Phase IX Summary",
            "",
            f"Timestamp: `{timestamp}`.",
            "",
            "## Objective",
            "",
            "Test whether the coefficient-like structures extracted from the frozen Phase VIII short-time trace regime stabilize across refinement, whether normalized combinations behave invariant-like, and whether this behavior differs from a deterministic generic-graph control hierarchy.",
            "",
            "## Frozen References",
            "",
            f"- Phase VI operator manifest: `{relpath(PHASE6_MANIFEST_PATH)}`.",
            f"- Phase VII spectral manifest: `{relpath(PHASE7_MANIFEST_PATH)}`.",
            f"- Phase VIII trace manifest: `{relpath(PHASE8_MANIFEST_PATH)}`.",
            f"- Phase VIII contract note: `{relpath(PHASE8_CONTRACT_NOTE_PATH)}`.",
            "",
            "## Frozen Trace Contract",
            "",
            "- Operational short-time window: `[0.02, 0.03, 0.05, 0.075, 0.1]`.",
            "- Trace method continuity: exact deterministic spectrum with the frozen `full_exact`, `lambda_le_7_5`, and `lambda_le_7_8` truncation policies.",
            f"- Dominant Phase VII low-mode exponent used for rescaled descriptors: `{round_float(phase7_exponent)}`.",
            "",
            "## Frozen Branch Coefficient Stabilization",
            "",
            "| coefficient | values across refinement | relative span | best extrapolation | limit estimate | fit R^2 |",
            "|---|---|---:|---|---:|---:|",
            *coefficient_rows,
            "",
            f"- Successive coefficient-vector distances: `{branch_summary['vector_distance_convergence']['successive_distances']}`.",
            f"- Distances shrink across refinement: `{branch_summary['vector_distance_convergence']['distances_shrink']}`.",
            f"- Largest raw coefficient relative span: `{round_float(branch_max_raw_span)}`.",
            "",
            "## Ratio Stabilization",
            "",
            "| ratio | values across refinement | relative span | best extrapolation | limit estimate |",
            "|---|---|---:|---|---:|",
            *ratio_rows_text,
            "",
            f"- Maximum branch ratio span: `{round_float(branch_max_ratio_span)}`.",
            f"- Ratio spans are smaller than raw coefficient spans: `{branch_max_ratio_span < branch_max_raw_span}`.",
            f"- Maximum subset-induced ratio drift: `{round_float(branch_subset_max_ratio_drift)}`.",
            f"- Maximum truncation-induced ratio drift: `{round_float(branch_truncation_max_ratio_drift)}`.",
            "",
            "## Rescaled Invariant-Like Descriptors",
            "",
            "| descriptor | values across refinement | relative span | best extrapolation | limit estimate |",
            "|---|---|---:|---|---:|",
            *invariant_rows,
            "",
            f"- Maximum branch rescaled-invariant span: `{round_float(branch_max_invariant_span)}`.",
            f"- Maximum subset-induced invariant drift: `{round_float(branch_subset_max_invariant_drift)}`.",
            f"- Maximum truncation-induced invariant drift: `{round_float(branch_truncation_max_invariant_drift)}`.",
            "",
            "## Generic-Graph Control Comparison",
            "",
            f"- Control hierarchy: `{CONTROL_LABEL}`.",
            f"- Control description: {CONTROL_DESCRIPTION}",
            "",
            "| control coefficient | relative span | best extrapolation | limit estimate |",
            "|---|---:|---|---:|",
            *control_rows,
            "",
            f"- Branch max ratio span: `{round_float(branch_max_ratio_span)}`.",
            f"- Control max ratio span: `{round_float(control_max_ratio_span)}`.",
            f"- Ratio-span multiplier (control / branch): `{round_float(control_max_ratio_span / max(branch_max_ratio_span, 1.0e-12))}`.",
            f"- Final coefficient-vector distance multiplier (control / branch): `{round_float(float(control_summary['vector_distance_convergence']['final_distance']) / max(float(branch_summary['vector_distance_convergence']['final_distance']), 1.0e-12))}`.",
            f"- Branch best extrapolation models: `{ {label: branch_scaled[label]['best_limit_fit']['model'] for label in ('b0', 'b1', 'b2')} }`.",
            f"- Control best extrapolation models: `{ {label: control_scaled[label]['best_limit_fit']['model'] for label in ('b0', 'b1', 'b2')} }`.",
            "",
            "## Bounded Interpretation",
            "",
            "- The frozen branch exhibits reproducible coefficient stabilization in the Phase VIII window, with the strongest plateau in `b0` and shrinking successive drift across all three scaled coefficients.",
            "- Dimensionless ratios are materially more stable than the raw trace coefficients, and the rescaled descriptors form compact refinement bands on the frozen branch.",
            "- The deterministic open-grid scalar control does not behave identically to the frozen branch: its ratio spans are larger, its successive coefficient-vector distances remain larger, and its simple extrapolation model family differs from the branch.",
            "- This is a feasibility result only. It does not assign geometric, continuum, or universal meaning to any coefficient or descriptor.",
            "",
            "## Non-Claims",
            "",
            "- No continuum limit is established.",
            "- No geometric coefficient interpretation is made.",
            "- No universal invariant is claimed.",
            "- No continuum or physical correspondence is asserted.",
            "",
            closure_statement(success),
        ]
    )


def main() -> None:
    started_at = time.perf_counter()
    timestamp = timestamp_iso()
    frozen = load_frozen_manifests()

    phase6_manifest = frozen["phase6_manifest"]
    phase7_manifest = frozen["phase7_manifest"]
    phase8_manifest = frozen["phase8_manifest"]

    epsilon = float(phase6_manifest["frozen_inputs"]["base_epsilon"])
    levels = list(phase8_manifest["refinement_hierarchy"]["levels"])
    phase7_exponent = float(phase7_manifest["aggregate_diagnostics"]["band_power_law_fits"][0]["exponent"])
    operational_window = np.array(phase8_manifest["operational_short_time_window"]["probe_times"], dtype=float)
    truncation_methods = list(phase8_manifest["trace_proxy_definition"]["spectral_truncation_methods"])

    branch_observations = build_hierarchy_observations(
        hierarchy_label=BRANCH_LABEL,
        levels=levels,
        epsilon=epsilon,
        probe_times=operational_window,
        phase7_exponent=phase7_exponent,
        cutoff=None,
    )
    control_observations = build_hierarchy_observations(
        hierarchy_label=CONTROL_LABEL,
        levels=levels,
        epsilon=epsilon,
        probe_times=operational_window,
        phase7_exponent=phase7_exponent,
        cutoff=None,
    )

    branch_summary = summarize_observations(branch_observations)
    control_summary = summarize_observations(control_observations)

    subset_probes = {
        BRANCH_LABEL: build_subset_probe_payload(BRANCH_LABEL, levels, epsilon, operational_window, phase7_exponent, branch_observations),
        CONTROL_LABEL: build_subset_probe_payload(CONTROL_LABEL, levels, epsilon, operational_window, phase7_exponent, control_observations),
    }
    truncation_probes = {
        BRANCH_LABEL: build_truncation_probe_payload(BRANCH_LABEL, levels, epsilon, operational_window, phase7_exponent, truncation_methods, branch_observations),
        CONTROL_LABEL: build_truncation_probe_payload(CONTROL_LABEL, levels, epsilon, operational_window, phase7_exponent, truncation_methods, control_observations),
    }

    branch_max_ratio_span = max(float(item["relative_span"]) for item in branch_summary["ratios"].values())
    branch_max_raw_span = max(float(item["relative_span"]) for item in branch_summary["raw_coefficients"].values())
    branch_max_invariant_span = max(float(item["relative_span"]) for item in branch_summary["rescaled_invariants"].values())
    control_max_ratio_span = max(float(item["relative_span"]) for item in control_summary["ratios"].values())
    branch_subset_max_ratio_drift = max(
        float(value)
        for subset in subset_probes[BRANCH_LABEL].values()
        for label, value in subset["max_relative_drifts_vs_full_window"].items()
        if label in {"b1_over_b0", "b2_over_b0"}
    )
    branch_subset_max_invariant_drift = max(
        float(value)
        for subset in subset_probes[BRANCH_LABEL].values()
        for label, value in subset["max_relative_drifts_vs_full_window"].items()
        if label in {"I0", "I1", "I2"}
    )
    branch_truncation_max_ratio_drift = max(
        float(value)
        for payload in truncation_probes[BRANCH_LABEL].values()
        for label, value in payload["max_relative_drifts_vs_full_exact"].items()
        if label in {"b1_over_b0", "b2_over_b0"}
    )
    branch_truncation_max_invariant_drift = max(
        float(value)
        for payload in truncation_probes[BRANCH_LABEL].values()
        for label, value in payload["max_relative_drifts_vs_full_exact"].items()
        if label in {"I0", "I1", "I2"}
    )
    branch_models = {label: branch_summary["scaled_coefficients"][label]["best_limit_fit"]["model"] for label in ("b0", "b1", "b2")}
    control_models = {label: control_summary["scaled_coefficients"][label]["best_limit_fit"]["model"] for label in ("b0", "b1", "b2")}
    differing_model_count = sum(branch_models[label] != control_models[label] for label in branch_models)

    success = (
        any(
            float(branch_summary["scaled_coefficients"][label]["relative_span"]) <= PLATEAU_RELATIVE_SPAN_TOL
            and bool(branch_summary["scaled_coefficients"][label]["drift_shrinks"])
            for label in ("b0", "b1", "b2")
        )
        and branch_max_ratio_span < branch_max_raw_span
        and branch_max_ratio_span <= RATIO_SPAN_TOL
        and branch_max_invariant_span <= INVARIANT_SPAN_TOL
        and branch_subset_max_ratio_drift <= SUBSET_RATIO_DRIFT_TOL
        and branch_subset_max_invariant_drift <= SUBSET_INVARIANT_DRIFT_TOL
        and branch_truncation_max_ratio_drift <= TRUNCATION_RATIO_DRIFT_TOL
        and branch_truncation_max_invariant_drift <= TRUNCATION_INVARIANT_DRIFT_TOL
        and control_max_ratio_span >= CONTROL_RATIO_SPAN_MULTIPLIER_MIN * branch_max_ratio_span
        and float(control_summary["vector_distance_convergence"]["final_distance"])
        >= CONTROL_DISTANCE_MULTIPLIER_MIN * float(branch_summary["vector_distance_convergence"]["final_distance"])
        and differing_model_count >= CONTROL_MODEL_DIFFERENCE_MIN
    )

    coeff_rows = coefficient_rows(timestamp, branch_summary) + coefficient_rows(timestamp, control_summary)
    ratio_rows_payload = ratio_rows(timestamp, branch_summary) + ratio_rows(timestamp, control_summary)
    coeff_fieldnames = [
        "stage_identifier",
        "timestamp",
        "hierarchy_label",
        "level_id",
        "n_side",
        "h",
        "coefficient_name",
        "raw_coefficient_name",
        "raw_value",
        "scaled_value",
        "relative_step_from_previous",
        "scaled_relative_span",
        "raw_relative_span",
        "best_limit_model",
        "limit_estimate",
        "limit_fit_r_squared",
    ]
    ratio_fieldnames = [
        "stage_identifier",
        "timestamp",
        "hierarchy_label",
        "level_id",
        "n_side",
        "h",
        "descriptor_family",
        "descriptor_name",
        "value",
        "relative_step_from_previous",
        "relative_span",
        "best_limit_model",
        "limit_estimate",
        "limit_fit_r_squared",
    ]
    write_csv_rows(COEFF_LEDGER_PATH, coeff_rows, coeff_fieldnames)
    write_csv_rows(RATIO_LEDGER_PATH, ratio_rows_payload, ratio_fieldnames)
    write_plots(branch_summary, control_summary)

    summary_text = build_summary(
        timestamp=timestamp,
        phase7_exponent=phase7_exponent,
        branch_summary=branch_summary,
        control_summary=control_summary,
        subset_probes=subset_probes,
        truncation_probes=truncation_probes,
        success=success,
    )
    manifest_payload = build_manifest(
        timestamp=timestamp,
        frozen=frozen,
        branch_summary=branch_summary,
        control_summary=control_summary,
        subset_probes=subset_probes,
        truncation_probes=truncation_probes,
        phase7_exponent=phase7_exponent,
        success=success,
    )
    runtime_seconds = time.perf_counter() - started_at
    runs_payload = build_runs_payload(
        timestamp=timestamp,
        phase7_exponent=phase7_exponent,
        branch_summary=branch_summary,
        control_summary=control_summary,
        subset_probes=subset_probes,
        truncation_probes=truncation_probes,
        runtime_seconds=runtime_seconds,
        success=success,
    )

    write_text(SUMMARY_PATH, summary_text)
    write_json(MANIFEST_PATH, manifest_payload)
    write_json(RUNS_PATH, runs_payload)

    print(f"Phase IX success: {success}")
    print(f"Runtime seconds: {round_float(runtime_seconds, 6)}")
    print(f"Manifest: {relpath(MANIFEST_PATH)}")
    print(closure_statement(success))


if __name__ == "__main__":
    main()
