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

from haos_core import read_json, write_csv_rows, write_json

PHASE_ROOT = ROOT / "phase16-temporal-ordering"
CONFIG_PATH = PHASE_ROOT / "configs" / "phase16_temporal_config.json"
RUNS_ROOT = PHASE_ROOT / "runs"
PLOTS_ROOT = PHASE_ROOT / "plots"

PHASE6_MANIFEST_PATH = ROOT / "phase6-operator" / "phase6_operator_manifest.json"
PHASE7_MANIFEST_PATH = ROOT / "phase7-spectral" / "phase7_spectral_manifest.json"
PHASE8_MANIFEST_PATH = ROOT / "phase8-trace" / "phase8_trace_manifest.json"
PHASE9_MANIFEST_PATH = ROOT / "phase9-invariants" / "phase9_manifest.json"
PHASE10_MANIFEST_PATH = ROOT / "phase10-bridge" / "phase10_manifest.json"
PHASE11_MANIFEST_PATH = ROOT / "phase11-protection" / "phase11_manifest.json"
PHASE12_MANIFEST_PATH = ROOT / "phase12-interactions" / "phase12_manifest.json"
PHASE13_MANIFEST_PATH = ROOT / "phase13-sector-formation" / "phase13_manifest.json"
PHASE14_MANIFEST_PATH = ROOT / "phase14-collective-dynamics" / "phase14_manifest.json"
PHASE15_MANIFEST_PATH = ROOT / "phase15-propagation" / "phase15_manifest.json"
PHASE15_CONFIG_PATH = ROOT / "phase15-propagation" / "configs" / "phase15_propagation_config.json"

EVENT_LEDGER_PATH = RUNS_ROOT / "phase16_event_ordering_ledger.csv"
MONOTONIC_LEDGER_PATH = RUNS_ROOT / "phase16_monotonic_parameter_ledger.csv"
FRONT_LEDGER_PATH = RUNS_ROOT / "phase16_front_arrival_ordering.csv"
ROBUSTNESS_LEDGER_PATH = RUNS_ROOT / "phase16_ordering_robustness.csv"
RUNS_PATH = RUNS_ROOT / "phase16_runs.json"
SUMMARY_PATH = PHASE_ROOT / "phase16_summary.md"
MANIFEST_PATH = PHASE_ROOT / "phase16_manifest.json"

CONSISTENCY_PLOT_PATH = PLOTS_ROOT / "phase16_ordering_consistency_score_vs_refinement.svg"
MONOTONIC_PLOT_PATH = PLOTS_ROOT / "phase16_monotonic_evolution_variable_vs_progression.svg"
FRONT_VAR_PLOT_PATH = PLOTS_ROOT / "phase16_front_arrival_ordering_variance_vs_refinement.svg"
ROBUSTNESS_PLOT_PATH = PLOTS_ROOT / "phase16_ordering_graph_distance_under_perturbation.svg"

PHASE_NAME = "phase16-temporal-ordering"
STAGE_IDENTIFIER = "phase16-temporal-ordering"
BRANCH_LABEL = "frozen_branch"
CONTROL_LABEL = "periodic_diagonal_augmented_control"
CLAIM_BOUNDARY = (
    "Phase XVI is limited to emergent temporal-ordering feasibility diagnostics on the frozen "
    "collective sector and frozen operator hierarchy. It does not assert physical time, "
    "relativistic causality, metric structure, spacetime emergence, or continuum temporal models."
)
EVENT_LABELS = [
    "radius_half",
    "dispersion_half",
    "spectral_half",
    "low_k_half",
    "width_half",
    "front_near",
    "front_far",
]
PRIMARY_PARAMETER = "integrated_radius_progress"


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def timestamp_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def round_float(value: float, digits: int = 12) -> float:
    return round(float(value), digits)


def closure_statement(success: bool) -> str:
    if success:
        return "Phase XVI establishes emergent temporal-ordering feasibility for the frozen operator hierarchy."
    return "Phase XVI does not yet establish emergent temporal-ordering feasibility for the frozen operator hierarchy."


def repo_rel(path: Path) -> str:
    return str(Path(path).resolve().relative_to(ROOT))


def periodic_delta(points: np.ndarray, anchor: np.ndarray) -> np.ndarray:
    delta = np.asarray(points, dtype=float) - np.asarray(anchor, dtype=float)
    return (delta + 0.5) % 1.0 - 0.5


def torus_distance(left: np.ndarray, right: np.ndarray) -> float:
    delta = periodic_delta(np.asarray([left], dtype=float), np.asarray(right, dtype=float))[0]
    return float(np.sqrt(np.sum(delta * delta)))


def build_points(n_side: int) -> np.ndarray:
    return np.array([[i / n_side, j / n_side] for i in range(n_side) for j in range(n_side)], dtype=float)


def scalar_eigen_grid(hierarchy_label: str, n_side: int, epsilon: float, diagonal_weight_scale: float) -> np.ndarray:
    h_value = 1.0 / float(n_side)
    edge_weight = math.exp(-(h_value * h_value) / (2.0 * epsilon))
    mode_x = np.arange(n_side, dtype=float)
    mode_y = np.arange(n_side, dtype=float)
    cos_x = np.cos(2.0 * math.pi * mode_x / n_side)[:, None]
    cos_y = np.cos(2.0 * math.pi * mode_y / n_side)[None, :]
    base = 2.0 * edge_weight * (2.0 - cos_x - cos_y)
    if hierarchy_label == BRANCH_LABEL:
        return np.asarray(base, dtype=float)
    diagonal_weight = edge_weight * float(diagonal_weight_scale)
    return np.asarray(base + 4.0 * diagonal_weight * (1.0 - cos_x * cos_y), dtype=float)


def first_positive_gap(eigen_grid: np.ndarray) -> float:
    positive = np.sort(np.asarray(eigen_grid[eigen_grid > 1.0e-12], dtype=float))
    return float(positive[0]) if len(positive) else 0.0


def evolve_state(state: np.ndarray, eigen_grid: np.ndarray, tau_value: float) -> np.ndarray:
    n_side = int(round(math.sqrt(state.size)))
    field = state.reshape((n_side, n_side))
    spectrum = np.fft.fftn(field)
    evolved = np.fft.ifftn(np.exp((-tau_value) * eigen_grid) * spectrum)
    return np.asarray(evolved.reshape(-1), dtype=complex)


def phase_kick(points: np.ndarray, axis: str, bias_strength: float) -> np.ndarray:
    if axis == "x":
        coordinate = points[:, 0] - 0.5
    else:
        coordinate = points[:, 1] - 0.5
    return np.exp(1j * 2.0 * math.pi * float(bias_strength) * coordinate)


def evolve_with_bias(
    state0: np.ndarray,
    eigen_grid: np.ndarray,
    points: np.ndarray,
    tau_grid: list[float],
    bias_strength: float,
    bias_remove_tau: float,
    phase_axis: str,
) -> list[np.ndarray]:
    states = [np.asarray(state0, dtype=complex)]
    current = np.asarray(state0, dtype=complex)
    last_tau = float(tau_grid[0])
    kick = phase_kick(points, phase_axis, bias_strength)
    for tau_value in tau_grid[1:]:
        dt = float(tau_value) - float(last_tau)
        current = evolve_state(current, eigen_grid, dt)
        if float(bias_strength) > 0.0 and float(tau_value) <= float(bias_remove_tau) + 1.0e-12:
            current = np.asarray(current * kick, dtype=complex)
        states.append(np.asarray(current, dtype=complex))
        last_tau = float(tau_value)
    return states


def low_mode_occupancy(state: np.ndarray, eigen_grid: np.ndarray, band_multiplier: float) -> float:
    gap = first_positive_gap(eigen_grid)
    if gap <= 0.0:
        return 0.0
    n_side = int(round(math.sqrt(state.size)))
    spectrum = np.fft.fftn(state.reshape((n_side, n_side)))
    energy = np.abs(spectrum) ** 2
    mask = (eigen_grid > 1.0e-12) & (eigen_grid <= float(band_multiplier) * gap + 1.0e-12)
    total = max(float(np.sum(energy)), 1.0e-12)
    return float(np.sum(energy[mask]) / total)


def build_candidate_state(points: np.ndarray, center: np.ndarray) -> np.ndarray:
    delta = periodic_delta(points, np.asarray(center, dtype=float))
    sigma = 0.08
    amplitude = np.exp(-0.5 * ((delta[:, 0] / sigma) ** 2 + (delta[:, 1] / sigma) ** 2))
    phase = np.exp(1j * 2.0 * math.pi * (delta[:, 0] + delta[:, 1]))
    state = amplitude * phase
    state /= max(np.linalg.norm(state), 1.0e-12)
    return np.asarray(state, dtype=complex)


def gaussian_window(points: np.ndarray, center: np.ndarray, sigma: float) -> np.ndarray:
    delta = periodic_delta(points, np.asarray(center, dtype=float))
    squared_radius = np.sum(delta * delta, axis=1)
    return np.exp(-0.5 * squared_radius / max(float(sigma) * float(sigma), 1.0e-12))


def local_mass(state: np.ndarray, window: np.ndarray) -> float:
    probabilities = np.abs(state) ** 2
    total = max(float(np.sum(probabilities)), 1.0e-12)
    return float(np.sum(probabilities * window) / total)


def spatial_metrics(state: np.ndarray, points: np.ndarray) -> dict[str, float]:
    probabilities = np.abs(state) ** 2
    probabilities = probabilities / max(float(np.sum(probabilities)), 1.0e-12)
    peak_position = points[int(np.argmax(probabilities))]
    delta = periodic_delta(points, peak_position)
    width = math.sqrt(float(np.sum(probabilities * np.sum(delta * delta, axis=1))))
    concentration = float(np.max(probabilities))
    participation_ratio = float(1.0 / np.sum(probabilities * probabilities))
    return {
        "localization_width": float(width),
        "spatial_concentration": concentration,
        "participation_ratio": participation_ratio,
    }


def generate_layout(seed: int, ensemble_size: int, min_separation_physical: float, max_attempts: int) -> list[np.ndarray]:
    rng = np.random.default_rng(int(seed))
    centers: list[np.ndarray] = []
    attempts = 0
    while len(centers) < int(ensemble_size) and attempts < int(max_attempts):
        candidate = np.asarray(rng.random(2), dtype=float)
        if all(torus_distance(candidate, center) >= float(min_separation_physical) - 1.0e-12 for center in centers):
            centers.append(candidate)
        attempts += 1
    if len(centers) != int(ensemble_size):
        raise ValueError(f"could not place {ensemble_size} candidates with seed {seed}")
    return centers


def density_profile_x(centers: np.ndarray, weights: np.ndarray, bin_count: int) -> np.ndarray:
    profile = np.zeros(int(bin_count), dtype=float)
    for center, weight in zip(centers, weights):
        index = min(int(math.floor(float(center[0]) * bin_count)), int(bin_count) - 1)
        profile[index] += float(weight)
    total = max(float(np.sum(profile)), 1.0e-12)
    return np.asarray(profile / total, dtype=float)


def fluctuation_summary(profile: np.ndarray) -> dict[str, Any]:
    centered = np.asarray(profile, dtype=float) - float(np.mean(profile))
    spectrum = np.abs(np.fft.rfft(centered)) ** 2
    low_k_power = float(np.sum(spectrum[1:3])) if spectrum.size > 1 else 0.0
    total_power = float(np.sum(spectrum[1:])) if spectrum.size > 1 else 0.0
    return {
        "low_k_power": low_k_power,
        "low_k_fraction": low_k_power / max(total_power, 1.0e-12),
    }


def response_profile_x(centers: np.ndarray, baseline_weights: np.ndarray, disturbed_weights: np.ndarray, bin_count: int) -> np.ndarray:
    baseline = density_profile_x(centers, baseline_weights, bin_count)
    disturbed = density_profile_x(centers, disturbed_weights, bin_count)
    response = np.abs(disturbed - baseline)
    total = max(float(np.sum(response)), 1.0e-12)
    return np.asarray(response / total, dtype=float)


def initial_collective_state(
    ref_states: list[np.ndarray],
    centers: np.ndarray,
    density_gradient_strength: float,
) -> tuple[np.ndarray, np.ndarray]:
    gradient_weights = 1.0 + density_gradient_strength * (0.5 - centers[:, 0])
    ensemble_state0 = np.sum(np.vstack([weight * state for weight, state in zip(gradient_weights, ref_states)]), axis=0)
    ensemble_state0 = np.asarray(ensemble_state0 / max(np.linalg.norm(ensemble_state0), 1.0e-12), dtype=complex)
    return ensemble_state0, np.asarray(gradient_weights, dtype=float)


def build_disturbed_state(
    probe_name: str,
    base_weights: np.ndarray,
    ref_states: list[np.ndarray],
    anchor_index: int,
    pulse_amplitude: float,
    removal_amplitude: float,
) -> np.ndarray:
    weights = np.asarray(base_weights, dtype=float).copy()
    if probe_name == "density_pulse":
        weights[anchor_index] *= float(pulse_amplitude)
    elif probe_name == "candidate_removal":
        weights[anchor_index] *= float(removal_amplitude)
    else:
        raise ValueError(f"unsupported static disturbance probe: {probe_name}")
    disturbed = np.sum(np.vstack([weight * state for weight, state in zip(weights, ref_states)]), axis=0)
    disturbed = np.asarray(disturbed / max(np.linalg.norm(disturbed), 1.0e-12), dtype=complex)
    return disturbed


def cumulative_trapezoid(tau_grid: list[float], values: list[float]) -> list[float]:
    out = [0.0]
    total = 0.0
    for left, right, value_left, value_right in zip(tau_grid, tau_grid[1:], values, values[1:]):
        dt = float(right) - float(left)
        total += 0.5 * dt * (max(float(value_left), 0.0) + max(float(value_right), 0.0))
        out.append(float(total))
    return out


def threshold_crossing_time(tau_grid: list[float], values: list[float], fraction: float, absolute: bool = False) -> float:
    transformed = [abs(float(value)) if absolute else float(value) for value in values]
    target = float(fraction) * max(max(transformed), 1.0e-12)
    for tau_value, value in zip(tau_grid, transformed):
        if value >= target - 1.0e-12:
            return float(tau_value)
    return float(tau_grid[-1])


def chain_signature(event_times: dict[str, float]) -> str:
    ordered = [label for label, _ in sorted(event_times.items(), key=lambda item: (float(item[1]), item[0]))]
    return ">".join(ordered)


def event_graph(event_times: dict[str, float]) -> dict[tuple[str, str], float]:
    edges: dict[tuple[str, str], float] = {}
    for left in EVENT_LABELS:
        for right in EVENT_LABELS:
            if left == right:
                continue
            if float(event_times[right]) > float(event_times[left]) + 1.0e-12:
                edges[(left, right)] = float(event_times[right] - event_times[left])
    return edges


def chain_distance(left_signature: str, right_signature: str) -> float:
    left = left_signature.split(">")
    right = right_signature.split(">")
    left_pos = {label: index for index, label in enumerate(left)}
    right_pos = {label: index for index, label in enumerate(right)}
    if len(left) <= 1:
        return 0.0
    return float(np.mean([abs(left_pos[label] - right_pos[label]) for label in EVENT_LABELS]) / (len(EVENT_LABELS) - 1))


def event_time_distance(left_times: dict[str, float], right_times: dict[str, float]) -> float:
    return float(np.mean([abs(float(left_times[label]) - float(right_times[label])) for label in EVENT_LABELS]))


def graph_distance(left_graph: dict[tuple[str, str], float], right_graph: dict[tuple[str, str], float]) -> float:
    keys = sorted(set(left_graph) | set(right_graph))
    if not keys:
        return 0.0
    return float(np.mean([abs(float(left_graph.get(key, 0.0)) - float(right_graph.get(key, 0.0))) for key in keys]))


def monotonicity_score(values: list[float]) -> tuple[int, float]:
    reversals = sum(1 for left, right in zip(values, values[1:]) if float(right) < float(left) - 1.0e-12)
    if len(values) <= 1:
        return 0, 1.0
    score = float(sum(1 for left, right in zip(values, values[1:]) if float(right) >= float(left) - 1.0e-12) / (len(values) - 1))
    return reversals, score


def deterministic_edge_noise_factor(left_label: str, right_label: str, scale: float) -> float:
    checksum = sum(ord(char) for char in (left_label + ":" + right_label))
    direction = 1.0 if checksum % 2 else -1.0
    return float(1.0 + direction * float(scale))


def perturb_event_graph(base_graph: dict[tuple[str, str], float], scale: float) -> dict[tuple[str, str], float]:
    return {
        key: float(value) * deterministic_edge_noise_factor(key[0], key[1], scale)
        for key, value in base_graph.items()
    }


def compute_run(
    hierarchy_label: str,
    n_side: int,
    ensemble_size: int,
    seed: int,
    probe_name: str,
    *,
    epsilon: float,
    diagonal_weight_scale: float,
    tau_grid: list[float],
    density_gradient_strength: float,
    response_threshold_fraction: float,
    pulse_amplitude: float,
    removal_amplitude: float,
    bias_strength: float,
    bias_remove_tau: float,
    phase_axis: str,
    band_multiplier: float,
    x_bin_count: int,
    min_sep: float,
    max_attempts: int,
    mass_sigma: float,
    event_thresholds: dict[str, float],
) -> dict[str, Any]:
    points = build_points(n_side)
    eigen_grid = scalar_eigen_grid(hierarchy_label, n_side, epsilon, diagonal_weight_scale)
    centers = np.asarray(generate_layout(seed, ensemble_size, min_sep, max_attempts), dtype=float)
    ref_states = [build_candidate_state(points, center) for center in centers]
    mass_windows = [gaussian_window(points, center, mass_sigma) for center in centers]
    anchor_index = int(np.argmin(centers[:, 0]))
    anchor_center = np.asarray(centers[anchor_index], dtype=float)
    distances = np.asarray([torus_distance(center, anchor_center) for center in centers], dtype=float)

    baseline_state0, base_weights = initial_collective_state(ref_states, centers, density_gradient_strength)
    baseline_states = [evolve_state(baseline_state0, eigen_grid, tau_value) for tau_value in tau_grid]

    if probe_name == "bias_onset":
        disturbed_states = evolve_with_bias(
            baseline_state0,
            eigen_grid,
            points,
            tau_grid,
            bias_strength,
            bias_remove_tau,
            phase_axis,
        )
    else:
        disturbed_state0 = build_disturbed_state(
            probe_name,
            base_weights,
            ref_states,
            anchor_index,
            pulse_amplitude,
            removal_amplitude,
        )
        disturbed_states = [evolve_state(disturbed_state0, eigen_grid, tau_value) for tau_value in tau_grid]

    snapshots: list[dict[str, Any]] = []
    max_response = 0.0
    for tau_value, base_state, disturbed_state in zip(tau_grid, baseline_states, disturbed_states):
        base_local = np.asarray([local_mass(base_state, window) for window in mass_windows], dtype=float)
        disturbed_local = np.asarray([local_mass(disturbed_state, window) for window in mass_windows], dtype=float)
        response_vector = np.asarray(np.abs(disturbed_local - base_local), dtype=float)
        max_response = max(max_response, float(np.max(response_vector)))
        response_profile = response_profile_x(centers, base_local, disturbed_local, x_bin_count)
        fluctuation = fluctuation_summary(response_profile)
        response_dispersion = float(np.sum(response_vector * distances) / max(np.sum(response_vector), 1.0e-12))
        disturbed_width = spatial_metrics(disturbed_state, points)["localization_width"]
        baseline_width = spatial_metrics(base_state, points)["localization_width"]
        spectral_shift = float(
            low_mode_occupancy(disturbed_state, eigen_grid, band_multiplier)
            - low_mode_occupancy(base_state, eigen_grid, band_multiplier)
        )
        snapshots.append(
            {
                "tau": float(tau_value),
                "response_vector": response_vector,
                "response_dispersion": response_dispersion,
                "spectral_shift": spectral_shift,
                "width_shift": float(disturbed_width - baseline_width),
                "low_k_response_fraction": float(fluctuation["low_k_fraction"]),
            }
        )

    threshold_value = float(response_threshold_fraction * max(max_response, 1.0e-12))
    latencies = np.full(len(centers), float(tau_grid[-1]), dtype=float)
    for index in range(len(centers)):
        for snapshot in snapshots:
            if float(snapshot["response_vector"][index]) >= threshold_value - 1.0e-12:
                latencies[index] = float(snapshot["tau"])
                break

    disturbance_radius_values = []
    response_dispersion_values = []
    spectral_shift_values = []
    low_k_values = []
    width_shift_values = []
    for snapshot in snapshots:
        response_vector = np.asarray(snapshot["response_vector"], dtype=float)
        active = response_vector >= threshold_value - 1.0e-12
        active_distances = distances[active]
        disturbance_radius_values.append(float(np.max(active_distances)) if active_distances.size else 0.0)
        response_dispersion_values.append(float(snapshot["response_dispersion"]))
        spectral_shift_values.append(float(snapshot["spectral_shift"]))
        low_k_values.append(float(snapshot["low_k_response_fraction"]))
        width_shift_values.append(float(snapshot["width_shift"]))

    non_anchor_mask = distances > 1.0e-12
    non_anchor_distances = np.asarray(distances[non_anchor_mask], dtype=float)
    non_anchor_latencies = np.asarray(latencies[non_anchor_mask], dtype=float)
    if non_anchor_distances.size:
        median_distance = float(np.median(non_anchor_distances))
        near_mask = non_anchor_distances <= median_distance + 1.0e-12
        far_mask = non_anchor_distances >= median_distance - 1.0e-12
        near_latency = float(np.mean(non_anchor_latencies[near_mask]))
        far_latency = float(np.mean(non_anchor_latencies[far_mask]))
        near_std = float(np.std(non_anchor_latencies[near_mask]))
        far_std = float(np.std(non_anchor_latencies[far_mask]))
    else:
        near_latency = float(tau_grid[-1])
        far_latency = float(tau_grid[-1])
        near_std = 0.0
        far_std = 0.0

    event_times = {
        "radius_half": threshold_crossing_time(tau_grid, disturbance_radius_values, event_thresholds["radius_fraction"]),
        "dispersion_half": threshold_crossing_time(tau_grid, response_dispersion_values, event_thresholds["dispersion_fraction"]),
        "spectral_half": threshold_crossing_time(tau_grid, spectral_shift_values, event_thresholds["spectral_fraction"], absolute=True),
        "low_k_half": threshold_crossing_time(tau_grid, low_k_values, event_thresholds["low_k_fraction"]),
        "width_half": threshold_crossing_time(tau_grid, width_shift_values, event_thresholds["width_fraction"], absolute=True),
        "front_near": near_latency,
        "front_far": far_latency,
    }
    event_chain = chain_signature(event_times)
    graph = event_graph(event_times)

    monotonic_series = {
        "integrated_radius_progress": cumulative_trapezoid(tau_grid, disturbance_radius_values),
        "integrated_dispersion_progress": cumulative_trapezoid(tau_grid, response_dispersion_values),
        "integrated_low_k_progress": cumulative_trapezoid(tau_grid, low_k_values),
    }
    monotonic_summary = {}
    for parameter_name, values in monotonic_series.items():
        reversals, score = monotonicity_score(values)
        monotonic_summary[parameter_name] = {
            "reversal_count": int(reversals),
            "monotonicity_score": float(score),
            "final_value": float(values[-1]),
        }

    return {
        "hierarchy_label": hierarchy_label,
        "n_side": int(n_side),
        "ensemble_size": int(ensemble_size),
        "seed": int(seed),
        "probe_name": probe_name,
        "event_times": event_times,
        "event_chain": event_chain,
        "event_graph": graph,
        "monotonic_series": monotonic_series,
        "monotonic_summary": monotonic_summary,
        "front_summary": {
            "near_latency": near_latency,
            "far_latency": far_latency,
            "near_std": near_std,
            "far_std": far_std,
            "front_distance": float(abs(far_latency - near_latency)),
            "front_ordering_valid": bool(near_latency <= far_latency + 1.0e-12),
        },
        "support_observables": {
            "disturbance_radius_values": disturbance_radius_values,
            "response_dispersion_values": response_dispersion_values,
            "spectral_shift_values": spectral_shift_values,
            "low_k_values": low_k_values,
            "width_shift_values": width_shift_values,
        },
    }


def successive_event_distance(
    results: dict[tuple[str, int, int, int, str], dict[str, Any]],
    hierarchy_label: str,
    probe_name: str,
    levels: list[int],
    ensemble_sizes: list[int],
    seeds: list[int],
) -> tuple[float, dict[int, float]]:
    distances = []
    level_means: dict[int, list[float]] = {int(level): [] for level in levels}
    for ensemble_size in ensemble_sizes:
        for seed in seeds:
            previous = None
            for n_side in levels:
                current = results[(hierarchy_label, int(n_side), int(ensemble_size), int(seed), probe_name)]
                if previous is not None:
                    distance_value = event_time_distance(previous["event_times"], current["event_times"])
                    distances.append(distance_value)
                    level_means[int(n_side)].append(distance_value)
                previous = current
    if not distances:
        return 0.0, {int(level): 0.0 for level in levels}
    return float(np.mean(distances)), {int(level): float(np.mean(values)) if values else 0.0 for level, values in level_means.items()}


def successive_front_distance(
    results: dict[tuple[str, int, int, int, str], dict[str, Any]],
    hierarchy_label: str,
    probe_name: str,
    levels: list[int],
    ensemble_sizes: list[int],
    seeds: list[int],
) -> tuple[float, dict[int, float]]:
    distances = []
    level_means: dict[int, list[float]] = {int(level): [] for level in levels}
    for ensemble_size in ensemble_sizes:
        for seed in seeds:
            previous = None
            for n_side in levels:
                current = results[(hierarchy_label, int(n_side), int(ensemble_size), int(seed), probe_name)]
                pair = np.asarray(
                    [current["front_summary"]["near_latency"], current["front_summary"]["far_latency"]],
                    dtype=float,
                )
                if previous is not None:
                    previous_pair = np.asarray(
                        [previous["front_summary"]["near_latency"], previous["front_summary"]["far_latency"]],
                        dtype=float,
                    )
                    distance_value = float(np.mean(np.abs(pair - previous_pair)))
                    distances.append(distance_value)
                    level_means[int(n_side)].append(distance_value)
                previous = current
    if not distances:
        return 0.0, {int(level): 0.0 for level in levels}
    return float(np.mean(distances)), {int(level): float(np.mean(values)) if values else 0.0 for level, values in level_means.items()}


def front_variance_by_level(
    results: dict[tuple[str, int, int, int, str], dict[str, Any]],
    hierarchy_label: str,
    probe_name: str,
    levels: list[int],
    ensemble_sizes: list[int],
    seeds: list[int],
) -> dict[int, float]:
    level_values: dict[int, float] = {}
    for n_side in levels:
        near_values = []
        far_values = []
        for ensemble_size in ensemble_sizes:
            for seed in seeds:
                current = results[(hierarchy_label, int(n_side), int(ensemble_size), int(seed), probe_name)]
                near_values.append(float(current["front_summary"]["near_latency"]))
                far_values.append(float(current["front_summary"]["far_latency"]))
        level_values[int(n_side)] = float(0.5 * (np.std(near_values) + np.std(far_values)))
    return level_values


def robustness_summary(
    baseline_results: dict[tuple[str, int, int, int, str], dict[str, Any]],
    perturbed_results: dict[tuple[str, int, int, int, str], dict[str, Any]],
    hierarchy_label: str,
    probe_name: str,
    levels: list[int],
    ensemble_sizes: list[int],
    seeds: list[int],
    parameter_name: str,
) -> tuple[float, float, float, float]:
    event_distances = []
    chain_distances = []
    graph_distances = []
    parameter_distances = []
    for n_side in levels:
        for ensemble_size in ensemble_sizes:
            for seed in seeds:
                base = baseline_results[(hierarchy_label, int(n_side), int(ensemble_size), int(seed), probe_name)]
                perturbed = perturbed_results[(hierarchy_label, int(n_side), int(ensemble_size), int(seed), probe_name)]
                event_distances.append(event_time_distance(base["event_times"], perturbed["event_times"]))
                chain_distances.append(chain_distance(base["event_chain"], perturbed["event_chain"]))
                graph_distances.append(graph_distance(base["event_graph"], perturbed["event_graph"]))
                parameter_distances.append(
                    float(
                        np.mean(
                            np.abs(
                                np.asarray(base["monotonic_series"][parameter_name], dtype=float)
                                - np.asarray(perturbed["monotonic_series"][parameter_name], dtype=float)
                            )
                        )
                    )
                )
    return (
        float(np.mean(event_distances)) if event_distances else 0.0,
        float(np.mean(chain_distances)) if chain_distances else 0.0,
        float(np.mean(graph_distances)) if graph_distances else 0.0,
        float(np.mean(parameter_distances)) if parameter_distances else 0.0,
    )


def plot_ordering_consistency(consistency_scores: dict[str, dict[int, float]], levels: list[int]) -> None:
    plt.figure(figsize=(8.8, 5.4))
    for hierarchy_label, marker in ((BRANCH_LABEL, "o"), (CONTROL_LABEL, "s")):
        values = [float(consistency_scores[hierarchy_label][int(level)]) for level in levels]
        plt.plot(levels, values, marker=marker, label=hierarchy_label)
    plt.xlabel("n_side")
    plt.ylabel("ordering consistency score")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(CONSISTENCY_PLOT_PATH, format="svg")
    plt.close()


def plot_monotonic_series(
    monotonic_rows: list[dict[str, Any]],
    tau_grid: list[float],
    levels: list[int],
    primary_probe: str,
    parameter_name: str,
    ensemble_size: int,
) -> None:
    plt.figure(figsize=(8.8, 5.4))
    level = int(levels[-1])
    for hierarchy_label, linestyle in ((BRANCH_LABEL, "-"), (CONTROL_LABEL, "--")):
        series = []
        for tau_value in tau_grid:
            values = [
                float(row["value"])
                for row in monotonic_rows
                if row["hierarchy_label"] == hierarchy_label
                and row["probe_name"] == primary_probe
                and row["parameter_name"] == parameter_name
                and int(row["n_side"]) == level
                and int(row["ensemble_size"]) == int(ensemble_size)
                and abs(float(row["tau"]) - float(tau_value)) <= 1.0e-12
            ]
            series.append(float(np.mean(values)) if values else 0.0)
        plt.plot(tau_grid, series, linestyle=linestyle, marker="o", label=hierarchy_label)
    plt.xlabel("tau")
    plt.ylabel(parameter_name.replace("_", " "))
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(MONOTONIC_PLOT_PATH, format="svg")
    plt.close()


def plot_front_variance(front_variance: dict[str, dict[int, float]], levels: list[int]) -> None:
    plt.figure(figsize=(8.8, 5.4))
    for hierarchy_label, marker in ((BRANCH_LABEL, "o"), (CONTROL_LABEL, "s")):
        values = [float(front_variance[hierarchy_label][int(level)]) for level in levels]
        plt.plot(levels, values, marker=marker, label=hierarchy_label)
    plt.xlabel("n_side")
    plt.ylabel("front-arrival variance")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(FRONT_VAR_PLOT_PATH, format="svg")
    plt.close()


def plot_robustness(robustness_rows: list[dict[str, Any]], perturbation_kinds: list[str]) -> None:
    plt.figure(figsize=(9.2, 5.6))
    x_positions = np.arange(len(perturbation_kinds), dtype=float)
    width = 0.35
    for offset, hierarchy_label in ((-0.18, BRANCH_LABEL), (0.18, CONTROL_LABEL)):
        values = []
        for perturbation_kind in perturbation_kinds:
            matches = [
                float(row["event_time_distance"])
                for row in robustness_rows
                if row["hierarchy_label"] == hierarchy_label and row["perturbation_kind"] == perturbation_kind
            ]
            values.append(float(np.mean(matches)) if matches else 0.0)
        plt.bar(x_positions + offset, values, width=width, label=hierarchy_label)
    plt.xticks(x_positions, perturbation_kinds, rotation=15)
    plt.ylabel("ordering-graph distance proxy")
    plt.grid(axis="y", alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(ROBUSTNESS_PLOT_PATH, format="svg")
    plt.close()


def load_frozen_inputs() -> dict[str, Any]:
    config = read_json(CONFIG_PATH)
    phase6 = read_json(PHASE6_MANIFEST_PATH)
    phase7 = read_json(PHASE7_MANIFEST_PATH)
    phase8 = read_json(PHASE8_MANIFEST_PATH)
    phase9 = read_json(PHASE9_MANIFEST_PATH)
    phase10 = read_json(PHASE10_MANIFEST_PATH)
    phase11 = read_json(PHASE11_MANIFEST_PATH)
    phase12 = read_json(PHASE12_MANIFEST_PATH)
    phase13 = read_json(PHASE13_MANIFEST_PATH)
    phase14 = read_json(PHASE14_MANIFEST_PATH)
    phase15 = read_json(PHASE15_MANIFEST_PATH)
    phase15_config = read_json(PHASE15_CONFIG_PATH)

    manifests = [phase6, phase7, phase8, phase9, phase10, phase11, phase12, phase13, phase14, phase15]
    if not all(bool(manifest.get("success")) for manifest in manifests):
        raise ValueError("Phase XVI requires successful frozen Phases VI through XV.")

    candidate_label = str(config["candidate_label"])
    if candidate_label != str(phase15["candidate_label"]):
        raise ValueError("Phase XVI candidate label must match the frozen Phase XV candidate label.")
    if [int(value) for value in config["refinement_levels"]] != [int(value) for value in phase15["refinement_hierarchy"]["levels"]]:
        raise ValueError("Phase XVI refinement levels must match the frozen Phase XV refinement hierarchy.")
    if [int(value) for value in config["ensemble_sizes"]] != [int(value) for value in phase15["evaluation_protocol"]["ensemble_sizes"]]:
        raise ValueError("Phase XVI ensemble sizes must match the frozen Phase XV ensemble sizes.")
    if [int(value) for value in config["layout_seeds"]] != [int(value) for value in phase15["evaluation_protocol"]["layout_seeds"]]:
        raise ValueError("Phase XVI layout seeds must match the frozen Phase XV layout seeds.")
    if [float(value) for value in config["persistence_tau_grid"]] != [float(value) for value in phase15["evaluation_protocol"]["persistence_tau_grid"]]:
        raise ValueError("Phase XVI tau grid must match the frozen Phase XV tau grid.")
    if float(config["evaluation_tau"]) != float(phase15["evaluation_protocol"]["evaluation_tau"]):
        raise ValueError("Phase XVI evaluation tau must match the frozen Phase XV evaluation tau.")
    if float(config["layout_rule"]["min_separation_physical"]) != float(phase15["evaluation_protocol"]["minimum_physical_separation"]):
        raise ValueError("Phase XVI minimum separation must match the frozen Phase XV layout rule.")
    if float(config["transport_probe"]["response_threshold_fraction"]) != float(phase15["evaluation_protocol"]["response_threshold_fraction"]):
        raise ValueError("Phase XVI response threshold fraction must match the frozen Phase XV threshold.")
    if float(config["transport_probe"]["density_gradient_strength"]) != float(phase15["evaluation_protocol"]["transport_probe"]["density_gradient_strength"]):
        raise ValueError("Phase XVI density gradient must match the frozen Phase XV transport probe.")
    if float(config["disturbance_probes"]["density_pulse_amplitude"]) != float(phase15["evaluation_protocol"]["disturbance_probes"]["density_pulse_amplitude"]):
        raise ValueError("Phase XVI density pulse amplitude must match the frozen Phase XV disturbance probe.")
    if float(config["disturbance_probes"]["candidate_removal_amplitude"]) != float(phase15["evaluation_protocol"]["disturbance_probes"]["candidate_removal_amplitude"]):
        raise ValueError("Phase XVI removal amplitude must match the frozen Phase XV disturbance probe.")
    if float(config["disturbance_probes"]["bias_strength"]) != float(phase15["evaluation_protocol"]["disturbance_probes"]["bias_strength"]):
        raise ValueError("Phase XVI bias strength must match the frozen Phase XV disturbance probe.")
    if float(config["disturbance_probes"]["bias_remove_tau"]) != float(phase15["evaluation_protocol"]["disturbance_probes"]["bias_remove_tau"]):
        raise ValueError("Phase XVI bias removal tau must match the frozen Phase XV disturbance probe.")
    if str(config["disturbance_probes"]["phase_kick_axis"]) != str(phase15["evaluation_protocol"]["disturbance_probes"]["phase_kick_axis"]):
        raise ValueError("Phase XVI phase axis must match the frozen Phase XV disturbance probe.")
    if [str(value) for value in config["ordering_probes"]] != ["density_pulse", "candidate_removal", "bias_onset"]:
        raise ValueError("Phase XVI ordering probes must reuse the frozen Phase XV disturbance families.")
    if str(config["ordering_primary_probe"]) != "bias_onset":
        raise ValueError("Phase XVI primary ordering probe must be the frozen nontrivial bias-onset probe.")
    if float(config["control_graph"]["diagonal_weight_scale"]) != float(phase15_config["control_graph"]["diagonal_weight_scale"]):
        raise ValueError("Phase XVI control graph must match the frozen Phase XV control graph.")

    return {
        "config": config,
        "phase6": phase6,
        "phase7": phase7,
        "phase8": phase8,
        "phase9": phase9,
        "phase10": phase10,
        "phase11": phase11,
        "phase12": phase12,
        "phase13": phase13,
        "phase14": phase14,
        "phase15": phase15,
        "phase15_config": phase15_config,
    }


def main() -> None:
    start_time = time.perf_counter()
    RUNS_ROOT.mkdir(parents=True, exist_ok=True)
    PLOTS_ROOT.mkdir(parents=True, exist_ok=True)

    payload = load_frozen_inputs()
    config = payload["config"]
    phase6 = payload["phase6"]

    epsilon = float(phase6["frozen_inputs"]["base_epsilon"])
    diagonal_weight_scale = float(config["control_graph"]["diagonal_weight_scale"])
    levels = [int(value) for value in config["refinement_levels"]]
    ensemble_sizes = [int(value) for value in config["ensemble_sizes"]]
    seeds = [int(value) for value in config["layout_seeds"]]
    tau_grid = [float(value) for value in config["persistence_tau_grid"]]
    min_sep = float(config["layout_rule"]["min_separation_physical"])
    max_attempts = int(config["layout_rule"]["max_attempts"])
    mass_sigma = float(config["candidate_windows"]["local_mass_sigma"])
    band_multiplier = float(config["spectral_identity_proxy"]["low_mode_band_multiplier"])
    density_gradient_strength = float(config["transport_probe"]["density_gradient_strength"])
    response_threshold_fraction = float(config["transport_probe"]["response_threshold_fraction"])
    x_bin_count = int(config["transport_probe"]["x_bin_count"])
    pulse_amplitude = float(config["disturbance_probes"]["density_pulse_amplitude"])
    removal_amplitude = float(config["disturbance_probes"]["candidate_removal_amplitude"])
    bias_strength = float(config["disturbance_probes"]["bias_strength"])
    bias_remove_tau = float(config["disturbance_probes"]["bias_remove_tau"])
    phase_axis = str(config["disturbance_probes"]["phase_kick_axis"])
    probe_names = [str(value) for value in config["ordering_probes"]]
    primary_probe = str(config["ordering_primary_probe"])
    event_thresholds = {str(key): float(value) for key, value in config["event_thresholds"].items()}
    monotonic_parameters = [str(value) for value in config["monotonic_parameters"]]
    thresholds = {str(key): float(value) for key, value in config["success_thresholds"].items()}

    baseline_results: dict[tuple[str, int, int, int, str], dict[str, Any]] = {}
    event_rows: list[dict[str, Any]] = []
    monotonic_rows: list[dict[str, Any]] = []
    front_rows: list[dict[str, Any]] = []
    robustness_rows: list[dict[str, Any]] = []

    for hierarchy_label in (BRANCH_LABEL, CONTROL_LABEL):
        for n_side in levels:
            for ensemble_size in ensemble_sizes:
                for seed in seeds:
                    for probe_name in probe_names:
                        result = compute_run(
                            hierarchy_label,
                            n_side,
                            ensemble_size,
                            seed,
                            probe_name,
                            epsilon=epsilon,
                            diagonal_weight_scale=diagonal_weight_scale,
                            tau_grid=tau_grid,
                            density_gradient_strength=density_gradient_strength,
                            response_threshold_fraction=response_threshold_fraction,
                            pulse_amplitude=pulse_amplitude,
                            removal_amplitude=removal_amplitude,
                            bias_strength=bias_strength,
                            bias_remove_tau=bias_remove_tau,
                            phase_axis=phase_axis,
                            band_multiplier=band_multiplier,
                            x_bin_count=x_bin_count,
                            min_sep=min_sep,
                            max_attempts=max_attempts,
                            mass_sigma=mass_sigma,
                            event_thresholds=event_thresholds,
                        )
                        baseline_results[(hierarchy_label, int(n_side), int(ensemble_size), int(seed), probe_name)] = result

                        ordered_labels = [label for label, _ in sorted(result["event_times"].items(), key=lambda item: (float(item[1]), item[0]))]
                        rank_map = {label: index + 1 for index, label in enumerate(ordered_labels)}
                        for event_label in EVENT_LABELS:
                            event_rows.append(
                                {
                                    "hierarchy_label": hierarchy_label,
                                    "level_id": f"n{n_side}",
                                    "n_side": int(n_side),
                                    "h": round_float(1.0 / float(n_side)),
                                    "ensemble_size": int(ensemble_size),
                                    "seed": int(seed),
                                    "probe_name": probe_name,
                                    "event_label": event_label,
                                    "event_time": round_float(result["event_times"][event_label], 6),
                                    "event_rank": int(rank_map[event_label]),
                                    "chain_signature": result["event_chain"],
                                    "is_primary_probe": bool(probe_name == primary_probe),
                                }
                            )

                        for parameter_name in monotonic_parameters:
                            values = result["monotonic_series"][parameter_name]
                            reversal_count = int(result["monotonic_summary"][parameter_name]["reversal_count"])
                            monotonicity_value = float(result["monotonic_summary"][parameter_name]["monotonicity_score"])
                            for tau_value, value in zip(tau_grid, values):
                                monotonic_rows.append(
                                    {
                                        "hierarchy_label": hierarchy_label,
                                        "level_id": f"n{n_side}",
                                        "n_side": int(n_side),
                                        "h": round_float(1.0 / float(n_side)),
                                        "ensemble_size": int(ensemble_size),
                                        "seed": int(seed),
                                        "probe_name": probe_name,
                                        "parameter_name": parameter_name,
                                        "tau": round_float(tau_value, 6),
                                        "value": round_float(value),
                                        "reversal_count": int(reversal_count),
                                        "monotonicity_score": round_float(monotonicity_value),
                                        "is_primary_parameter": bool(parameter_name == PRIMARY_PARAMETER),
                                    }
                                )

                        for band_name, latency_value, band_std in (
                            ("near", result["front_summary"]["near_latency"], result["front_summary"]["near_std"]),
                            ("far", result["front_summary"]["far_latency"], result["front_summary"]["far_std"]),
                        ):
                            front_rows.append(
                                {
                                    "hierarchy_label": hierarchy_label,
                                    "level_id": f"n{n_side}",
                                    "n_side": int(n_side),
                                    "h": round_float(1.0 / float(n_side)),
                                    "ensemble_size": int(ensemble_size),
                                    "seed": int(seed),
                                    "probe_name": probe_name,
                                    "distance_band": band_name,
                                    "mean_latency": round_float(latency_value, 6),
                                    "latency_std": round_float(band_std),
                                    "front_distance": round_float(result["front_summary"]["front_distance"]),
                                    "front_ordering_valid": bool(result["front_summary"]["front_ordering_valid"]),
                                    "chain_signature": result["event_chain"],
                                }
                            )

    perturbation_specs = {
        "density_gradient_plus": {
            "density_gradient_strength": float(config["robustness_perturbations"]["density_gradient_plus"]["density_gradient_strength"]),
            "pulse_amplitude": pulse_amplitude,
            "bias_strength": bias_strength,
        },
        "amplitude_shift": {
            "density_gradient_strength": density_gradient_strength,
            "pulse_amplitude": float(config["robustness_perturbations"]["amplitude_shift"]["density_pulse_amplitude"]),
            "bias_strength": float(config["robustness_perturbations"]["amplitude_shift"]["bias_strength"]),
        },
    }

    perturbed_results: dict[str, dict[tuple[str, int, int, int, str], dict[str, Any]]] = {}
    for perturbation_kind, spec in perturbation_specs.items():
        perturbed_results[perturbation_kind] = {}
        for hierarchy_label in (BRANCH_LABEL, CONTROL_LABEL):
            for n_side in levels:
                for ensemble_size in ensemble_sizes:
                    for seed in seeds:
                        result = compute_run(
                            hierarchy_label,
                            n_side,
                            ensemble_size,
                            seed,
                            primary_probe,
                            epsilon=epsilon,
                            diagonal_weight_scale=diagonal_weight_scale,
                            tau_grid=tau_grid,
                            density_gradient_strength=float(spec["density_gradient_strength"]),
                            response_threshold_fraction=response_threshold_fraction,
                            pulse_amplitude=float(spec["pulse_amplitude"]),
                            removal_amplitude=removal_amplitude,
                            bias_strength=float(spec["bias_strength"]),
                            bias_remove_tau=bias_remove_tau,
                            phase_axis=phase_axis,
                            band_multiplier=band_multiplier,
                            x_bin_count=x_bin_count,
                            min_sep=min_sep,
                            max_attempts=max_attempts,
                            mass_sigma=mass_sigma,
                            event_thresholds=event_thresholds,
                        )
                        perturbed_results[perturbation_kind][(hierarchy_label, int(n_side), int(ensemble_size), int(seed), primary_probe)] = result

    noise_scale = float(config["robustness_perturbations"]["event_graph_connectivity_noise"]["edge_weight_noise_scale"])
    for hierarchy_label in (BRANCH_LABEL, CONTROL_LABEL):
        for n_side in levels:
            for ensemble_size in ensemble_sizes:
                for seed in seeds:
                    baseline = baseline_results[(hierarchy_label, int(n_side), int(ensemble_size), int(seed), primary_probe)]
                    for perturbation_kind, perturb_results in perturbed_results.items():
                        perturbed = perturb_results[(hierarchy_label, int(n_side), int(ensemble_size), int(seed), primary_probe)]
                        robustness_rows.append(
                            {
                                "hierarchy_label": hierarchy_label,
                                "level_id": f"n{n_side}",
                                "n_side": int(n_side),
                                "h": round_float(1.0 / float(n_side)),
                                "ensemble_size": int(ensemble_size),
                                "seed": int(seed),
                                "probe_name": primary_probe,
                                "perturbation_kind": perturbation_kind,
                                "event_time_distance": round_float(event_time_distance(baseline["event_times"], perturbed["event_times"])),
                                "chain_distance": round_float(chain_distance(baseline["event_chain"], perturbed["event_chain"])),
                                "graph_distance": round_float(graph_distance(baseline["event_graph"], perturbed["event_graph"])),
                                "front_latency_distance": round_float(
                                    0.5
                                    * (
                                        abs(float(baseline["front_summary"]["near_latency"]) - float(perturbed["front_summary"]["near_latency"]))
                                        + abs(float(baseline["front_summary"]["far_latency"]) - float(perturbed["front_summary"]["far_latency"]))
                                    )
                                ),
                                "primary_parameter_distance": round_float(
                                    float(
                                        np.mean(
                                            np.abs(
                                                np.asarray(baseline["monotonic_series"][PRIMARY_PARAMETER], dtype=float)
                                                - np.asarray(perturbed["monotonic_series"][PRIMARY_PARAMETER], dtype=float)
                                            )
                                        )
                                    )
                                ),
                            }
                        )
                    noisy_graph = perturb_event_graph(baseline["event_graph"], noise_scale)
                    robustness_rows.append(
                        {
                            "hierarchy_label": hierarchy_label,
                            "level_id": f"n{n_side}",
                            "n_side": int(n_side),
                            "h": round_float(1.0 / float(n_side)),
                            "ensemble_size": int(ensemble_size),
                            "seed": int(seed),
                            "probe_name": primary_probe,
                            "perturbation_kind": "event_graph_connectivity_noise",
                            "event_time_distance": 0.0,
                            "chain_distance": 0.0,
                            "graph_distance": round_float(graph_distance(baseline["event_graph"], noisy_graph)),
                            "front_latency_distance": 0.0,
                            "primary_parameter_distance": 0.0,
                        }
                    )

    branch_event_distance, branch_event_level_distance = successive_event_distance(
        baseline_results,
        BRANCH_LABEL,
        primary_probe,
        levels,
        ensemble_sizes,
        seeds,
    )
    control_event_distance, control_event_level_distance = successive_event_distance(
        baseline_results,
        CONTROL_LABEL,
        primary_probe,
        levels,
        ensemble_sizes,
        seeds,
    )
    branch_front_distance, branch_front_level_distance = successive_front_distance(
        baseline_results,
        BRANCH_LABEL,
        primary_probe,
        levels,
        ensemble_sizes,
        seeds,
    )
    control_front_distance, control_front_level_distance = successive_front_distance(
        baseline_results,
        CONTROL_LABEL,
        primary_probe,
        levels,
        ensemble_sizes,
        seeds,
    )

    branch_monotonic_scores = [
        float(result["monotonic_summary"][PRIMARY_PARAMETER]["monotonicity_score"])
        for key, result in baseline_results.items()
        if key[0] == BRANCH_LABEL and key[4] == primary_probe
    ]
    control_monotonic_scores = [
        float(result["monotonic_summary"][PRIMARY_PARAMETER]["monotonicity_score"])
        for key, result in baseline_results.items()
        if key[0] == CONTROL_LABEL and key[4] == primary_probe
    ]
    branch_max_reversals = max(
        int(result["monotonic_summary"][PRIMARY_PARAMETER]["reversal_count"])
        for key, result in baseline_results.items()
        if key[0] == BRANCH_LABEL and key[4] == primary_probe
    )
    control_max_reversals = max(
        int(result["monotonic_summary"][PRIMARY_PARAMETER]["reversal_count"])
        for key, result in baseline_results.items()
        if key[0] == CONTROL_LABEL and key[4] == primary_probe
    )
    branch_monotonic_mean = float(np.mean(branch_monotonic_scores))
    control_monotonic_mean = float(np.mean(control_monotonic_scores))

    branch_front_variance = front_variance_by_level(baseline_results, BRANCH_LABEL, primary_probe, levels, ensemble_sizes, seeds)
    control_front_variance = front_variance_by_level(baseline_results, CONTROL_LABEL, primary_probe, levels, ensemble_sizes, seeds)
    branch_front_variance_mean = float(np.mean(list(branch_front_variance.values())))
    control_front_variance_mean = float(np.mean(list(control_front_variance.values())))

    branch_robustness_distance = float(
        np.mean(
            [
                float(row["event_time_distance"])
                for row in robustness_rows
                if row["hierarchy_label"] == BRANCH_LABEL and row["perturbation_kind"] in ("density_gradient_plus", "amplitude_shift")
            ]
        )
    )
    control_robustness_distance = float(
        np.mean(
            [
                float(row["event_time_distance"])
                for row in robustness_rows
                if row["hierarchy_label"] == CONTROL_LABEL and row["perturbation_kind"] in ("density_gradient_plus", "amplitude_shift")
            ]
        )
    )
    branch_noise_graph_distance = float(
        np.mean(
            [
                float(row["graph_distance"])
                for row in robustness_rows
                if row["hierarchy_label"] == BRANCH_LABEL and row["perturbation_kind"] == "event_graph_connectivity_noise"
            ]
        )
    )
    control_noise_graph_distance = float(
        np.mean(
            [
                float(row["graph_distance"])
                for row in robustness_rows
                if row["hierarchy_label"] == CONTROL_LABEL and row["perturbation_kind"] == "event_graph_connectivity_noise"
            ]
        )
    )

    branch_event_score = float(1.0 - branch_event_distance / max(float(tau_grid[-1]), 1.0e-12))
    control_event_score = float(1.0 - control_event_distance / max(float(tau_grid[-1]), 1.0e-12))

    success_flags = {
        "event_ordering_chains_stable": bool(branch_event_distance <= thresholds["max_primary_event_distance"]),
        "monotonic_evolution_parameter_identified": bool(
            branch_monotonic_mean >= thresholds["min_primary_monotonicity_score"] and branch_max_reversals == 0
        ),
        "front_ordering_dispersion_bounded": bool(branch_front_distance <= thresholds["max_primary_front_distance"]),
        "ordering_robustness_holds": bool(branch_robustness_distance <= thresholds["max_primary_robustness_distance"]),
    }
    control_flags = {
        "event_ordering_chains_stable": bool(control_event_distance <= thresholds["max_primary_event_distance"]),
        "monotonic_evolution_parameter_identified": bool(
            control_monotonic_mean >= thresholds["min_primary_monotonicity_score"] and control_max_reversals == 0
        ),
        "front_ordering_dispersion_bounded": bool(control_front_distance <= thresholds["max_primary_front_distance"]),
        "ordering_robustness_holds": bool(control_robustness_distance <= thresholds["max_primary_robustness_distance"]),
    }
    success_flags["control_hierarchy_different"] = bool(not all(control_flags.values()))
    success = all(success_flags.values())

    consistency_scores = {
        BRANCH_LABEL: {int(levels[0]): 1.0} | {level: float(1.0 - value / max(float(tau_grid[-1]), 1.0e-12)) for level, value in branch_event_level_distance.items() if int(level) != int(levels[0])},
        CONTROL_LABEL: {int(levels[0]): 1.0} | {level: float(1.0 - value / max(float(tau_grid[-1]), 1.0e-12)) for level, value in control_event_level_distance.items() if int(level) != int(levels[0])},
    }

    plot_ordering_consistency(consistency_scores, levels)
    plot_monotonic_series(monotonic_rows, tau_grid, levels, primary_probe, PRIMARY_PARAMETER, ensemble_sizes[1])
    plot_front_variance({BRANCH_LABEL: branch_front_variance, CONTROL_LABEL: control_front_variance}, levels)
    plot_robustness(robustness_rows, ["density_gradient_plus", "amplitude_shift", "event_graph_connectivity_noise"])

    write_csv_rows(
        EVENT_LEDGER_PATH,
        event_rows,
        [
            "hierarchy_label",
            "level_id",
            "n_side",
            "h",
            "ensemble_size",
            "seed",
            "probe_name",
            "event_label",
            "event_time",
            "event_rank",
            "chain_signature",
            "is_primary_probe",
        ],
    )
    write_csv_rows(
        MONOTONIC_LEDGER_PATH,
        monotonic_rows,
        [
            "hierarchy_label",
            "level_id",
            "n_side",
            "h",
            "ensemble_size",
            "seed",
            "probe_name",
            "parameter_name",
            "tau",
            "value",
            "reversal_count",
            "monotonicity_score",
            "is_primary_parameter",
        ],
    )
    write_csv_rows(
        FRONT_LEDGER_PATH,
        front_rows,
        [
            "hierarchy_label",
            "level_id",
            "n_side",
            "h",
            "ensemble_size",
            "seed",
            "probe_name",
            "distance_band",
            "mean_latency",
            "latency_std",
            "front_distance",
            "front_ordering_valid",
            "chain_signature",
        ],
    )
    write_csv_rows(
        ROBUSTNESS_LEDGER_PATH,
        robustness_rows,
        [
            "hierarchy_label",
            "level_id",
            "n_side",
            "h",
            "ensemble_size",
            "seed",
            "probe_name",
            "perturbation_kind",
            "event_time_distance",
            "chain_distance",
            "graph_distance",
            "front_latency_distance",
            "primary_parameter_distance",
        ],
    )

    runtime_seconds = float(time.perf_counter() - start_time)
    runs_payload = {
        "timestamp": timestamp_iso(),
        "phase": 16,
        "stage_identifier": STAGE_IDENTIFIER,
        "deterministic_seed_record": config["deterministic_seed_record"],
        "refinement_levels": levels,
        "ensemble_sizes": ensemble_sizes,
        "layout_seeds": seeds,
        "tau_grid": tau_grid,
        "ordering_probes": probe_names,
        "primary_ordering_probe": primary_probe,
        "monotonic_parameters": monotonic_parameters,
        "robustness_perturbations": config["robustness_perturbations"],
        "runtime_seconds": round_float(runtime_seconds, 6),
    }
    write_json(RUNS_PATH, runs_payload)

    summary_lines = [
        "# Phase XVI Summary",
        "",
        "Phase XVI tests whether the frozen collective sector supports a stable internal ordering structure derived from the already frozen Phase XV disturbance regime.",
        "",
        "## Primary Ordering Probe",
        "",
        f"- Primary probe: `{primary_probe}`",
        f"- Branch primary event distance: `{round_float(branch_event_distance)}`",
        f"- Control primary event distance: `{round_float(control_event_distance)}`",
        f"- Branch primary ordering score: `{round_float(branch_event_score)}`",
        f"- Control primary ordering score: `{round_float(control_event_score)}`",
        "",
        "## Aggregate Metrics",
        "",
        f"- Branch front-order distance: `{round_float(branch_front_distance)}`",
        f"- Control front-order distance: `{round_float(control_front_distance)}`",
        f"- Branch monotonicity mean: `{round_float(branch_monotonic_mean)}`",
        f"- Control monotonicity mean: `{round_float(control_monotonic_mean)}`",
        f"- Branch robustness distance: `{round_float(branch_robustness_distance)}`",
        f"- Control robustness distance: `{round_float(control_robustness_distance)}`",
        f"- Branch event-graph noise distance: `{round_float(branch_noise_graph_distance)}`",
        f"- Control event-graph noise distance: `{round_float(control_noise_graph_distance)}`",
        f"- Branch front variance mean: `{round_float(branch_front_variance_mean)}`",
        f"- Control front variance mean: `{round_float(control_front_variance_mean)}`",
        "",
        "## Success Flags",
        "",
        f"- Event ordering chains stable: `{success_flags['event_ordering_chains_stable']}`",
        f"- Monotonic evolution parameter identified: `{success_flags['monotonic_evolution_parameter_identified']}`",
        f"- Front ordering dispersion bounded: `{success_flags['front_ordering_dispersion_bounded']}`",
        f"- Ordering robustness holds: `{success_flags['ordering_robustness_holds']}`",
        f"- Control hierarchy different: `{success_flags['control_hierarchy_different']}`",
        "",
        "## Control Flags",
        "",
        f"- Control event ordering chains stable: `{control_flags['event_ordering_chains_stable']}`",
        f"- Control monotonic evolution parameter identified: `{control_flags['monotonic_evolution_parameter_identified']}`",
        f"- Control front ordering dispersion bounded: `{control_flags['front_ordering_dispersion_bounded']}`",
        f"- Control ordering robustness holds: `{control_flags['ordering_robustness_holds']}`",
        "",
        closure_statement(success),
    ]
    write_text(SUMMARY_PATH, "\n".join(summary_lines))

    manifest = {
        "timestamp": timestamp_iso(),
        "phase": 16,
        "phase_name": PHASE_NAME,
        "stage_identifier": STAGE_IDENTIFIER,
        "status": "passed" if success else "failed",
        "success": bool(success),
        "objective": "emergent_temporal_ordering_feasibility",
        "operator_manifest_reference": repo_rel(PHASE6_MANIFEST_PATH),
        "phase7_manifest_reference": repo_rel(PHASE7_MANIFEST_PATH),
        "phase8_manifest_reference": repo_rel(PHASE8_MANIFEST_PATH),
        "phase9_manifest_reference": repo_rel(PHASE9_MANIFEST_PATH),
        "phase10_manifest_reference": repo_rel(PHASE10_MANIFEST_PATH),
        "phase11_manifest_reference": repo_rel(PHASE11_MANIFEST_PATH),
        "phase12_manifest_reference": repo_rel(PHASE12_MANIFEST_PATH),
        "phase13_manifest_reference": repo_rel(PHASE13_MANIFEST_PATH),
        "phase14_manifest_reference": repo_rel(PHASE14_MANIFEST_PATH),
        "phase15_manifest_reference": repo_rel(PHASE15_MANIFEST_PATH),
        "candidate_label": config["candidate_label"],
        "primary_ordering_probe": primary_probe,
        "refinement_hierarchy": {
            "h_definition": "h = 1 / n_side on the unit-periodic lattice; it is the refinement-order parameter inherited from the frozen operator hierarchy.",
            "levels": levels,
        },
        "evaluation_protocol": {
            "persistence_tau_grid": [round_float(value, 6) for value in tau_grid],
            "ensemble_sizes": ensemble_sizes,
            "layout_seeds": seeds,
            "minimum_physical_separation": round_float(min_sep),
            "ordering_probes": probe_names,
            "event_thresholds": config["event_thresholds"],
            "robustness_perturbations": config["robustness_perturbations"],
        },
        "aggregate_metrics": {
            "branch_primary_event_distance": round_float(branch_event_distance),
            "control_primary_event_distance": round_float(control_event_distance),
            "branch_primary_ordering_score": round_float(branch_event_score),
            "control_primary_ordering_score": round_float(control_event_score),
            "branch_primary_front_distance": round_float(branch_front_distance),
            "control_primary_front_distance": round_float(control_front_distance),
            "branch_primary_monotonicity_mean": round_float(branch_monotonic_mean),
            "control_primary_monotonicity_mean": round_float(control_monotonic_mean),
            "branch_primary_max_reversals": int(branch_max_reversals),
            "control_primary_max_reversals": int(control_max_reversals),
            "branch_primary_robustness_distance": round_float(branch_robustness_distance),
            "control_primary_robustness_distance": round_float(control_robustness_distance),
            "branch_front_variance_mean": round_float(branch_front_variance_mean),
            "control_front_variance_mean": round_float(control_front_variance_mean),
            "branch_event_graph_noise_distance": round_float(branch_noise_graph_distance),
            "control_event_graph_noise_distance": round_float(control_noise_graph_distance),
        },
        "success_flags": success_flags,
        "control_flags": control_flags,
        "artifacts": {
            "event_ordering_ledger_csv": repo_rel(EVENT_LEDGER_PATH),
            "monotonic_parameter_ledger_csv": repo_rel(MONOTONIC_LEDGER_PATH),
            "front_arrival_ordering_csv": repo_rel(FRONT_LEDGER_PATH),
            "ordering_robustness_csv": repo_rel(ROBUSTNESS_LEDGER_PATH),
            "runs_json": repo_rel(RUNS_PATH),
            "summary_md": repo_rel(SUMMARY_PATH),
            "manifest_json": repo_rel(MANIFEST_PATH),
            "ordering_consistency_plot": repo_rel(CONSISTENCY_PLOT_PATH),
            "monotonic_plot": repo_rel(MONOTONIC_PLOT_PATH),
            "front_variance_plot": repo_rel(FRONT_VAR_PLOT_PATH),
            "robustness_plot": repo_rel(ROBUSTNESS_PLOT_PATH),
            "builder_script": repo_rel(Path(__file__)),
        },
        "runtime_seconds": round_float(runtime_seconds, 6),
        "claim_boundary": CLAIM_BOUNDARY,
        "conclusion": closure_statement(success),
    }
    write_json(MANIFEST_PATH, manifest)


if __name__ == "__main__":
    main()
