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

PHASE_ROOT = ROOT / "phase14-collective-dynamics"
CONFIG_PATH = PHASE_ROOT / "configs" / "phase14_collective_config.json"
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
PHASE13_CONFIG_PATH = ROOT / "phase13-sector-formation" / "configs" / "phase13_sector_config.json"

TRANSPORT_LEDGER_PATH = RUNS_ROOT / "phase14_transport_ledger.csv"
RELAXATION_PATH = RUNS_ROOT / "phase14_relaxation_times.csv"
FLUCTUATION_PATH = RUNS_ROOT / "phase14_fluctuation_spectra.csv"
RESPONSE_PATH = RUNS_ROOT / "phase14_density_response.csv"
EOS_PATH = RUNS_ROOT / "phase14_equation_of_state_proxy.csv"
RUNS_PATH = RUNS_ROOT / "phase14_runs.json"
SUMMARY_PATH = PHASE_ROOT / "phase14_summary.md"
MANIFEST_PATH = PHASE_ROOT / "phase14_manifest.json"

DRIFT_PLOT_PATH = PLOTS_ROOT / "phase14_density_centroid_drift_vs_time.svg"
VARIANCE_PLOT_PATH = PLOTS_ROOT / "phase14_variance_growth_vs_time.svg"
RELAX_PLOT_PATH = PLOTS_ROOT / "phase14_relaxation_time_vs_refinement.svg"
FLUCT_PLOT_PATH = PLOTS_ROOT / "phase14_fluctuation_power_vs_scale.svg"
EOS_PLOT_PATH = PLOTS_ROOT / "phase14_survival_fraction_vs_spacing.svg"

PHASE_NAME = "phase14-collective-dynamics"
STAGE_IDENTIFIER = "phase14-collective-dynamics"
BRANCH_LABEL = "frozen_branch"
CONTROL_LABEL = "periodic_diagonal_augmented_control"
CLAIM_BOUNDARY = (
    "Phase XIV is limited to collective-dynamics feasibility diagnostics on the frozen dilute "
    "candidate sector and frozen operator hierarchy. It does not assert hydrodynamics, "
    "thermodynamics, continuum medium laws, or physical correspondence."
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
        return "Phase XIV establishes collective dynamics feasibility for the frozen operator hierarchy."
    return "Phase XIV does not yet establish collective dynamics feasibility for the frozen operator hierarchy."


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


def windowed_overlap(state: np.ndarray, reference: np.ndarray, window: np.ndarray) -> float:
    weighted_state = np.asarray(state * window, dtype=complex)
    weighted_reference = np.asarray(reference * window, dtype=complex)
    norm_left = float(np.linalg.norm(weighted_state))
    norm_right = float(np.linalg.norm(weighted_reference))
    if norm_left <= 1.0e-12 or norm_right <= 1.0e-12:
        return 0.0
    return float(abs(np.vdot(weighted_reference, weighted_state)) / (norm_left * norm_right))


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


def connected_components(node_count: int, edges: list[tuple[int, int]]) -> list[list[int]]:
    adjacency = {index: set() for index in range(node_count)}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    seen: set[int] = set()
    components: list[list[int]] = []
    for start in range(node_count):
        if start in seen:
            continue
        stack = [start]
        component: list[int] = []
        while stack:
            node = stack.pop()
            if node in seen:
                continue
            seen.add(node)
            component.append(node)
            stack.extend(sorted(adjacency[node] - seen))
        components.append(sorted(component))
    return components


def circular_coordinate(values: np.ndarray, weights: np.ndarray) -> float:
    total = max(float(np.sum(weights)), 1.0e-12)
    angles = 2.0 * math.pi * np.asarray(values, dtype=float)
    x_value = float(np.sum(weights * np.cos(angles)) / total)
    y_value = float(np.sum(weights * np.sin(angles)) / total)
    if abs(x_value) <= 1.0e-12 and abs(y_value) <= 1.0e-12:
        return float(np.mean(values))
    angle = math.atan2(y_value, x_value)
    if angle < 0.0:
        angle += 2.0 * math.pi
    return float(angle / (2.0 * math.pi))


def weighted_centroid(centers: np.ndarray, weights: np.ndarray) -> np.ndarray:
    if float(np.sum(weights)) <= 1.0e-12:
        return np.asarray([0.5, 0.5], dtype=float)
    return np.asarray(
        [
            circular_coordinate(centers[:, 0], weights),
            circular_coordinate(centers[:, 1], weights),
        ],
        dtype=float,
    )


def weighted_variance(centers: np.ndarray, weights: np.ndarray, centroid: np.ndarray) -> float:
    total = max(float(np.sum(weights)), 1.0e-12)
    delta = periodic_delta(centers, centroid)
    return float(np.sum(weights * np.sum(delta * delta, axis=1)) / total)


def initial_spacing_metrics(centers: np.ndarray) -> tuple[float, float]:
    nearest_values = []
    pair_values = []
    for left in range(len(centers)):
        distances = []
        for right in range(len(centers)):
            if left == right:
                continue
            distance = torus_distance(centers[left], centers[right])
            distances.append(distance)
            if right > left:
                pair_values.append(distance)
        nearest_values.append(min(distances))
    return float(np.mean(nearest_values)), float(np.mean(pair_values))


def density_profile_x(centers: np.ndarray, weights: np.ndarray, bin_count: int) -> np.ndarray:
    profile = np.zeros(int(bin_count), dtype=float)
    for center, weight in zip(centers, weights):
        index = min(int(math.floor(float(center[0]) * bin_count)), int(bin_count) - 1)
        profile[index] += float(weight)
    total = max(float(np.sum(profile)), 1.0e-12)
    return np.asarray(profile / total, dtype=float)


def fluctuation_summary(profile: np.ndarray) -> dict[str, Any]:
    centered = np.asarray(profile, dtype=float) - float(np.mean(profile))
    amplitude = float(np.sqrt(np.mean(centered * centered)))
    spectrum = np.abs(np.fft.rfft(centered)) ** 2
    low_k_power = float(np.sum(spectrum[1:3])) if spectrum.size > 1 else 0.0
    k_values = np.arange(1, spectrum.size, dtype=float)
    positive_mask = spectrum[1:] > 1.0e-14
    if np.count_nonzero(positive_mask) >= 2:
        slope = float(np.polyfit(np.log(k_values[positive_mask]), np.log(spectrum[1:][positive_mask]), 1)[0])
    else:
        slope = 0.0
    return {
        "fluctuation_amplitude": amplitude,
        "low_k_power": low_k_power,
        "fluctuation_slope": slope,
        "power_values": [float(value) for value in spectrum[1:].tolist()],
    }


def estimate_relaxation_tau(tau_grid: list[float], mode_amplitudes: list[float]) -> float:
    initial = max(abs(float(mode_amplitudes[0])), 1.0e-12)
    normalized = [max(abs(float(value)) / initial, 1.0e-12) for value in mode_amplitudes]
    threshold = math.exp(-1.0)
    for tau_value, value in zip(tau_grid[1:], normalized[1:]):
        if value <= threshold:
            return float(tau_value)
    times = np.asarray(tau_grid[1:], dtype=float)
    values = np.asarray(normalized[1:], dtype=float)
    if len(times) >= 2:
        slope = float(np.polyfit(times, np.log(values), 1)[0])
        if slope < -1.0e-9:
            return float(min(10.0 * tau_grid[-1], -1.0 / slope))
    return float(10.0 * tau_grid[-1])


def phase_kick(points: np.ndarray, axis: str, bias_strength: float) -> np.ndarray:
    if axis == "x":
        coordinate = points[:, 0] - 0.5
    else:
        coordinate = points[:, 1] - 0.5
    return np.exp(1j * 2.0 * math.pi * float(bias_strength) * coordinate)


def evolve_with_phase_kicks(
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


def successive_distance_by_seed(rows: list[dict[str, Any]], key: str, hierarchy_label: str, levels: list[int], ensemble_sizes: list[int], seeds: list[int]) -> float:
    distances = []
    for ensemble_size in ensemble_sizes:
        for seed in seeds:
            values = []
            for n_side in levels:
                match = next(
                    (
                        float(row[key])
                        for row in rows
                        if row["hierarchy_label"] == hierarchy_label
                        and int(row["ensemble_size"]) == int(ensemble_size)
                        and int(row["seed"]) == int(seed)
                        and int(row["n_side"]) == int(n_side)
                    ),
                    None,
                )
                if match is None:
                    continue
                values.append(match)
            for left, right in zip(values, values[1:]):
                distances.append(abs(right - left))
    return float(max(distances)) if distances else 0.0


def seed_std_max(rows: list[dict[str, Any]], key: str, hierarchy_label: str) -> float:
    grouped: dict[tuple[int, int], list[float]] = {}
    for row in rows:
        if row["hierarchy_label"] != hierarchy_label:
            continue
        grouped.setdefault((int(row["n_side"]), int(row["ensemble_size"])), []).append(float(row[key]))
    if not grouped:
        return 0.0
    return float(max(np.std(values) for values in grouped.values()))


def linearity_error(rows: list[dict[str, Any]], hierarchy_label: str, small_bias: float, large_bias: float) -> float:
    grouped: dict[tuple[int, int, int], dict[float, float]] = {}
    for row in rows:
        if row["hierarchy_label"] != hierarchy_label:
            continue
        key = (int(row["n_side"]), int(row["ensemble_size"]), int(row["seed"]))
        grouped.setdefault(key, {})[float(row["bias_strength"])] = float(row["mode_response_eval"])
    errors = []
    for mapping in grouped.values():
        if small_bias not in mapping or large_bias not in mapping:
            continue
        small_value = abs(float(mapping[small_bias]))
        large_value = abs(float(mapping[large_bias]))
        if small_value <= 1.0e-12:
            continue
        errors.append(abs(large_value / max(2.0 * small_value, 1.0e-12) - 1.0))
    return float(np.mean(errors)) if errors else float("inf")


def metric_linearity_error(rows: list[dict[str, Any]], hierarchy_label: str, key: str, small_bias: float, large_bias: float) -> float:
    grouped: dict[tuple[int, int, int], dict[float, float]] = {}
    for row in rows:
        if row["hierarchy_label"] != hierarchy_label:
            continue
        grouped.setdefault((int(row["n_side"]), int(row["ensemble_size"]), int(row["seed"])), {})[float(row["bias_strength"])] = float(row[key])
    errors = []
    for mapping in grouped.values():
        if 0.0 not in mapping or small_bias not in mapping or large_bias not in mapping:
            continue
        small_delta = abs(float(mapping[small_bias]) - float(mapping[0.0]))
        large_delta = abs(float(mapping[large_bias]) - float(mapping[0.0]))
        if small_delta <= 1.0e-12:
            continue
        errors.append(abs(large_delta / max(2.0 * small_delta, 1.0e-12) - 1.0))
    return float(np.mean(errors)) if errors else float("inf")


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
    phase13_config = read_json(PHASE13_CONFIG_PATH)

    manifests = [phase6, phase7, phase8, phase9, phase10, phase11, phase12, phase13]
    if not all(bool(manifest.get("success")) for manifest in manifests):
        raise ValueError("Phase XIV requires successful frozen Phases VI through XIII.")

    candidate_label = str(config["candidate_label"])
    if candidate_label != str(phase11["candidate_label"]):
        raise ValueError("Phase XIV candidate must match the frozen Phase XI candidate label.")
    if candidate_label != str(phase12["candidate_label"]):
        raise ValueError("Phase XIV candidate must match the frozen Phase XII candidate label.")
    if candidate_label != str(phase13["candidate_label"]):
        raise ValueError("Phase XIV candidate must match the frozen Phase XIII candidate label.")

    if [int(value) for value in config["refinement_levels"]] != [int(value) for value in phase13["refinement_hierarchy"]["levels"]]:
        raise ValueError("Phase XIV refinement levels must match the frozen Phase XIII refinement hierarchy.")
    if [int(value) for value in config["ensemble_sizes"]] != [int(value) for value in phase13["evaluation_protocol"]["ensemble_sizes"]]:
        raise ValueError("Phase XIV ensemble sizes must match the frozen Phase XIII ensemble sizes.")
    if [int(value) for value in config["layout_seeds"]] != [int(value) for value in phase13["evaluation_protocol"]["layout_seeds"]]:
        raise ValueError("Phase XIV layout seeds must match the frozen Phase XIII layout seeds.")
    if float(config["layout_rule"]["min_separation_physical"]) != float(phase13["evaluation_protocol"]["minimum_physical_separation"]):
        raise ValueError("Phase XIV minimum separation must match the frozen Phase XIII layout rule.")
    if [float(value) for value in config["persistence_tau_grid"]] != [float(value) for value in phase13["evaluation_protocol"]["persistence_tau_grid"]]:
        raise ValueError("Phase XIV tau grid must match the frozen Phase XIII tau grid.")
    if float(config["evaluation_tau"]) != float(phase13["evaluation_protocol"]["evaluation_tau"]):
        raise ValueError("Phase XIV evaluation tau must match the frozen Phase XIII evaluation tau.")
    if phase8["operational_short_time_window"]["probe_times"] != phase13["evaluation_protocol"]["phase8_short_time_window"]:
        raise ValueError("Phase XIV must reuse the frozen Phase VIII short-time window carried through Phase XIII.")
    if float(config["control_graph"]["diagonal_weight_scale"]) != float(phase13_config["control_graph"]["diagonal_weight_scale"]):
        raise ValueError("Phase XIV control hierarchy must reuse the frozen Phase XIII control graph weight scale.")

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
        "phase13_config": phase13_config,
    }


def plot_centroid_drift(transport_rows: list[dict[str, Any]], tau_grid: list[float], levels: list[int], ensemble_size: int) -> None:
    plt.figure(figsize=(8.8, 5.4))
    level = int(levels[-1])
    for hierarchy_label, linestyle in ((BRANCH_LABEL, "-"), (CONTROL_LABEL, "--")):
        series = []
        for tau_value in tau_grid:
            matches = [
                float(row["value"])
                for row in transport_rows
                if row["hierarchy_label"] == hierarchy_label
                and int(row["n_side"]) == level
                and int(row["ensemble_size"]) == int(ensemble_size)
                and row["observable"] == "centroid_drift_x"
                and abs(float(row["tau"]) - float(tau_value)) <= 1.0e-12
            ]
            series.append(float(np.mean(matches)) if matches else 0.0)
        plt.plot(tau_grid, series, linestyle=linestyle, marker="o", label=f"{hierarchy_label}:N{ensemble_size}")
    plt.xlabel("tau")
    plt.ylabel("centroid drift in x")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(DRIFT_PLOT_PATH, format="svg")
    plt.close()


def plot_variance_growth(transport_rows: list[dict[str, Any]], tau_grid: list[float], levels: list[int], ensemble_size: int) -> None:
    plt.figure(figsize=(8.8, 5.4))
    level = int(levels[-1])
    for hierarchy_label, linestyle in ((BRANCH_LABEL, "-"), (CONTROL_LABEL, "--")):
        series = []
        for tau_value in tau_grid:
            matches = [
                float(row["value"])
                for row in transport_rows
                if row["hierarchy_label"] == hierarchy_label
                and int(row["n_side"]) == level
                and int(row["ensemble_size"]) == int(ensemble_size)
                and row["observable"] == "variance_growth"
                and abs(float(row["tau"]) - float(tau_value)) <= 1.0e-12
            ]
            series.append(float(np.mean(matches)) if matches else 0.0)
        plt.plot(tau_grid, series, linestyle=linestyle, marker="o", label=f"{hierarchy_label}:N{ensemble_size}")
    plt.xlabel("tau")
    plt.ylabel("variance growth")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(VARIANCE_PLOT_PATH, format="svg")
    plt.close()


def plot_relaxation_times(relaxation_rows: list[dict[str, Any]], levels: list[int]) -> None:
    plt.figure(figsize=(8.6, 5.2))
    for hierarchy_label, marker in ((BRANCH_LABEL, "o"), (CONTROL_LABEL, "s")):
        means = []
        for n_side in levels:
            matches = [
                float(row["collective_relaxation_tau"])
                for row in relaxation_rows
                if row["hierarchy_label"] == hierarchy_label and int(row["n_side"]) == int(n_side)
            ]
            means.append(float(np.mean(matches)) if matches else 0.0)
        plt.plot(levels, means, marker=marker, label=hierarchy_label)
    plt.xlabel("n_side")
    plt.ylabel("collective relaxation tau")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(RELAX_PLOT_PATH, format="svg")
    plt.close()


def plot_fluctuation_power(fluctuation_rows: list[dict[str, Any]], levels: list[int], ensemble_size: int, tau_value: float) -> None:
    plt.figure(figsize=(8.6, 5.2))
    level = int(levels[-1])
    for hierarchy_label, marker in ((BRANCH_LABEL, "o"), (CONTROL_LABEL, "s")):
        rows = [
            row
            for row in fluctuation_rows
            if row["hierarchy_label"] == hierarchy_label
            and int(row["n_side"]) == level
            and int(row["ensemble_size"]) == int(ensemble_size)
            and abs(float(row["tau"]) - float(tau_value)) <= 1.0e-12
        ]
        power_by_k: dict[int, list[float]] = {}
        for row in rows:
            power_by_k.setdefault(int(row["k_index"]), []).append(float(row["power"]))
        ks = sorted(power_by_k)
        values = [float(np.mean(power_by_k[k_value])) for k_value in ks]
        plt.plot(ks, values, marker=marker, label=hierarchy_label)
    plt.xlabel("x-bin wavenumber index")
    plt.ylabel("fluctuation power")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(FLUCT_PLOT_PATH, format="svg")
    plt.close()


def plot_survival_vs_spacing(eos_rows: list[dict[str, Any]]) -> None:
    plt.figure(figsize=(8.6, 5.2))
    for hierarchy_label, marker in ((BRANCH_LABEL, "o"), (CONTROL_LABEL, "s")):
        xs = [float(row["initial_nearest_neighbor_mean"]) for row in eos_rows if row["hierarchy_label"] == hierarchy_label]
        ys = [float(row["terminal_survival_fraction"]) for row in eos_rows if row["hierarchy_label"] == hierarchy_label]
        plt.scatter(xs, ys, marker=marker, label=hierarchy_label)
    plt.xlabel("initial nearest-neighbor mean")
    plt.ylabel("terminal survival fraction")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(EOS_PLOT_PATH, format="svg")
    plt.close()


def main() -> None:
    start_time = time.perf_counter()
    RUNS_ROOT.mkdir(parents=True, exist_ok=True)
    PLOTS_ROOT.mkdir(parents=True, exist_ok=True)

    payload = load_frozen_inputs()
    config = payload["config"]
    phase6 = payload["phase6"]
    phase8 = payload["phase8"]
    phase11 = payload["phase11"]
    phase13_config = payload["phase13_config"]

    epsilon = float(phase6["frozen_inputs"]["base_epsilon"])
    diagonal_weight_scale = float(config["control_graph"]["diagonal_weight_scale"])
    levels = [int(value) for value in config["refinement_levels"]]
    ensemble_sizes = [int(value) for value in config["ensemble_sizes"]]
    seeds = [int(value) for value in config["layout_seeds"]]
    tau_grid = [float(value) for value in config["persistence_tau_grid"]]
    evaluation_tau = float(config["evaluation_tau"])
    min_sep = float(config["layout_rule"]["min_separation_physical"])
    max_attempts = int(config["layout_rule"]["max_attempts"])
    identity_sigma = float(config["candidate_windows"]["identity_sigma"])
    mass_sigma = float(config["candidate_windows"]["local_mass_sigma"])
    min_identity = float(config["survival_thresholds"]["min_identity"])
    min_mass_ratio = float(config["survival_thresholds"]["min_local_mass_ratio"])
    band_multiplier = float(config["spectral_identity_proxy"]["low_mode_band_multiplier"])
    density_gradient_strength = float(config["transport_probe"]["density_gradient_strength"])
    x_bin_count = int(config["fluctuation_probe"]["x_bin_count"])
    bias_strengths = [float(value) for value in config["response_probe"]["bias_strengths"]]
    bias_remove_tau = float(config["response_probe"]["bias_remove_tau"])
    phase_axis = str(config["response_probe"]["phase_kick_axis"])

    cluster_distance = float(phase13_config["cluster_thresholds"]["distance_physical"])
    midpoint_identity_min = float(phase13_config["cluster_thresholds"]["midpoint_identity_min"])
    identity_ceiling = float(phase13_config["cluster_thresholds"]["identity_ceiling"])

    branch_micro_map = {
        int(level): float(value)
        for level, value in zip(levels, phase11["persistence_scaling"]["branch_tau_values"])
    }
    control_micro_map = {
        int(level): float(value)
        for level, value in zip(levels, phase11["persistence_scaling"]["control_tau_values"])
    }

    layout_map: dict[tuple[int, int], np.ndarray] = {}
    for ensemble_size in ensemble_sizes:
        for seed in seeds:
            layout_map[(int(ensemble_size), int(seed))] = np.asarray(
                generate_layout(seed, ensemble_size, min_sep, max_attempts),
                dtype=float,
            )

    transport_rows: list[dict[str, Any]] = []
    relaxation_rows: list[dict[str, Any]] = []
    fluctuation_rows: list[dict[str, Any]] = []
    response_rows: list[dict[str, Any]] = []
    eos_rows: list[dict[str, Any]] = []

    single_mode_baseline: dict[tuple[str, int, float], float] = {}
    baseline_response_reference: dict[tuple[str, int, int, int], dict[str, float]] = {}

    for hierarchy_label in (BRANCH_LABEL, CONTROL_LABEL):
        micro_map = branch_micro_map if hierarchy_label == BRANCH_LABEL else control_micro_map
        for n_side in levels:
            points = build_points(n_side)
            eigen_grid = scalar_eigen_grid(hierarchy_label, n_side, epsilon, diagonal_weight_scale)
            centered_candidate = build_candidate_state(points, np.array([0.5, 0.5], dtype=float))
            for tau_value in tau_grid:
                single_mode_baseline[(hierarchy_label, int(n_side), float(tau_value))] = low_mode_occupancy(
                    evolve_state(centered_candidate, eigen_grid, tau_value), eigen_grid, band_multiplier
                )

            for ensemble_size in ensemble_sizes:
                for seed in seeds:
                    centers = np.asarray(layout_map[(int(ensemble_size), int(seed))], dtype=float)
                    ref_states = [build_candidate_state(points, center) for center in centers]
                    identity_windows = [gaussian_window(points, center, identity_sigma) for center in centers]
                    mass_windows = [gaussian_window(points, center, mass_sigma) for center in centers]
                    midpoint_data: dict[tuple[int, int], tuple[np.ndarray, np.ndarray, float]] = {}
                    for left in range(len(centers)):
                        for right in range(left + 1, len(centers)):
                            midpoint = (centers[left] + periodic_delta(np.asarray([centers[right]]), centers[left])[0] / 2.0) % 1.0
                            midpoint_ref = build_candidate_state(points, midpoint)
                            midpoint_window = gaussian_window(points, midpoint, identity_sigma)
                            midpoint_data[(left, right)] = (
                                midpoint_ref,
                                midpoint_window,
                                torus_distance(centers[left], centers[right]),
                            )

                    gradient_weights = 1.0 + density_gradient_strength * (0.5 - centers[:, 0])
                    ensemble_state0 = np.sum(np.vstack([weight * state for weight, state in zip(gradient_weights, ref_states)]), axis=0)
                    ensemble_state0 = np.asarray(ensemble_state0 / max(np.linalg.norm(ensemble_state0), 1.0e-12), dtype=complex)
                    initial_local_masses = np.asarray([local_mass(ensemble_state0, window) for window in mass_windows], dtype=float)
                    initial_centroid = weighted_centroid(centers, initial_local_masses)
                    initial_variance = weighted_variance(centers, initial_local_masses, initial_centroid)
                    initial_mode_amplitude = float(
                        abs(np.sum(initial_local_masses * np.exp(-2j * math.pi * centers[:, 0]))) / max(float(np.sum(initial_local_masses)), 1.0e-12)
                    )
                    initial_nn_mean, initial_pair_mean = initial_spacing_metrics(centers)

                    baseline_states = [evolve_state(ensemble_state0, eigen_grid, tau_value) for tau_value in tau_grid]
                    baseline_mode_values: list[float] = []
                    baseline_diffusion_values: list[float] = []
                    baseline_low_k_values: list[float] = []
                    baseline_slope_values: list[float] = []

                    terminal_survival = 0.0
                    terminal_cluster_frequency = 0.0
                    terminal_spectral = 0.0

                    for tau_value, state_tau in zip(tau_grid, baseline_states):
                        candidate_identity = []
                        candidate_mass_ratio = []
                        candidate_local_masses = []
                        survivor_flags = []
                        for index in range(len(centers)):
                            identity_value = windowed_overlap(state_tau, ref_states[index], identity_windows[index])
                            mass_value = local_mass(state_tau, mass_windows[index])
                            mass_ratio = mass_value / max(float(initial_local_masses[index]), 1.0e-12)
                            candidate_identity.append(identity_value)
                            candidate_local_masses.append(mass_value)
                            candidate_mass_ratio.append(mass_ratio)
                            survivor_flags.append(bool(identity_value >= min_identity and mass_ratio >= min_mass_ratio))

                        mass_array = np.asarray(candidate_local_masses, dtype=float)
                        total_mass = max(float(np.sum(mass_array)), 1.0e-12)
                        centroid = weighted_centroid(centers, mass_array)
                        drift_x = float(periodic_delta(np.asarray([[centroid[0], 0.0]], dtype=float), np.asarray([initial_centroid[0], 0.0], dtype=float))[0, 0])
                        variance = weighted_variance(centers, mass_array, centroid)
                        variance_growth = float(max(variance - initial_variance, 0.0))
                        diffusion_proxy = variance_growth / max(float(tau_value), 1.0e-12) if tau_value > 0.0 else 0.0
                        mode_amplitude = float(abs(np.sum(mass_array * np.exp(-2j * math.pi * centers[:, 0]))) / total_mass)
                        density_redistribution = float(abs(mode_amplitude - initial_mode_amplitude))
                        density_rate = density_redistribution / max(float(tau_value), 1.0e-12) if tau_value > 0.0 else 0.0
                        survival_fraction = float(np.mean(survivor_flags))
                        mean_identity = float(np.mean(candidate_identity))
                        localization_width = spatial_metrics(state_tau, points)["localization_width"]
                        baseline_low_mode = single_mode_baseline[(hierarchy_label, int(n_side), float(tau_value))]
                        low_mode_value = low_mode_occupancy(state_tau, eigen_grid, band_multiplier)
                        spectral_ratio = low_mode_value / max(float(baseline_low_mode), 1.0e-12)
                        spectral_distortion = spectral_ratio - 1.0
                        profile = density_profile_x(centers, mass_array, x_bin_count)
                        fluctuation = fluctuation_summary(profile)

                        surviving_indices = [index for index, survives in enumerate(survivor_flags) if survives]
                        cluster_edges: list[tuple[int, int]] = []
                        for left in range(len(centers)):
                            if left not in surviving_indices:
                                continue
                            for right in range(left + 1, len(centers)):
                                if right not in surviving_indices:
                                    continue
                                midpoint_ref, midpoint_window, distance_value = midpoint_data[(left, right)]
                                midpoint_identity = windowed_overlap(state_tau, midpoint_ref, midpoint_window)
                                if (
                                    distance_value <= cluster_distance + 1.0e-12
                                    and midpoint_identity >= midpoint_identity_min
                                    and min(candidate_identity[left], candidate_identity[right]) <= identity_ceiling
                                ):
                                    cluster_edges.append((left, right))
                        components = connected_components(len(centers), cluster_edges)
                        clustered_components = [component for component in components if len(component) > 1]
                        cluster_frequency = float(sum(len(component) for component in clustered_components)) / float(ensemble_size)

                        for observable, value in (
                            ("centroid_drift_x", drift_x),
                            ("density_mode_amplitude", mode_amplitude),
                            ("density_redistribution_rate", density_rate),
                            ("variance_growth", variance_growth),
                            ("effective_diffusion_proxy", diffusion_proxy),
                            ("survival_fraction", survival_fraction),
                            ("localization_width", localization_width),
                            ("spectral_distortion", spectral_distortion),
                            ("mean_identity", mean_identity),
                        ):
                            transport_rows.append(
                                {
                                    "hierarchy_label": hierarchy_label,
                                    "level_id": f"n{n_side}",
                                    "n_side": int(n_side),
                                    "h": round_float(1.0 / float(n_side)),
                                    "ensemble_size": int(ensemble_size),
                                    "seed": int(seed),
                                    "tau": round_float(tau_value, 6),
                                    "observable": observable,
                                    "value": round_float(value),
                                }
                            )

                        for k_index, power_value in enumerate(fluctuation["power_values"], start=1):
                            fluctuation_rows.append(
                                {
                                    "hierarchy_label": hierarchy_label,
                                    "level_id": f"n{n_side}",
                                    "n_side": int(n_side),
                                    "h": round_float(1.0 / float(n_side)),
                                    "ensemble_size": int(ensemble_size),
                                    "seed": int(seed),
                                    "tau": round_float(tau_value, 6),
                                    "k_index": int(k_index),
                                    "power": round_float(power_value),
                                    "low_k_power": round_float(fluctuation["low_k_power"]),
                                    "fluctuation_amplitude": round_float(fluctuation["fluctuation_amplitude"]),
                                    "fluctuation_slope": round_float(fluctuation["fluctuation_slope"]),
                                }
                            )

                        baseline_mode_values.append(mode_amplitude)
                        baseline_diffusion_values.append(diffusion_proxy)
                        baseline_low_k_values.append(float(fluctuation["low_k_power"]))
                        baseline_slope_values.append(float(fluctuation["fluctuation_slope"]))

                        if abs(tau_value - evaluation_tau) <= 1.0e-12:
                            baseline_response_reference[(hierarchy_label, int(n_side), int(ensemble_size), int(seed))] = {
                                "centroid_drift_eval": drift_x,
                                "mode_amplitude_eval": mode_amplitude,
                                "mean_identity_eval": mean_identity,
                                "survival_fraction_eval": survival_fraction,
                                "spectral_distortion_eval": spectral_distortion,
                            }
                        if abs(tau_value - bias_remove_tau) <= 1.0e-12:
                            baseline_response_reference[(hierarchy_label, int(n_side), int(ensemble_size), int(seed))]["mode_amplitude_remove"] = mode_amplitude
                            baseline_response_reference[(hierarchy_label, int(n_side), int(ensemble_size), int(seed))]["centroid_drift_remove"] = drift_x
                        if abs(tau_value - tau_grid[-1]) <= 1.0e-12:
                            baseline_response_reference[(hierarchy_label, int(n_side), int(ensemble_size), int(seed))]["mode_amplitude_final"] = mode_amplitude
                            baseline_response_reference[(hierarchy_label, int(n_side), int(ensemble_size), int(seed))]["centroid_drift_final"] = drift_x
                            terminal_survival = survival_fraction
                            terminal_cluster_frequency = cluster_frequency
                            terminal_spectral = spectral_distortion

                    relaxation_tau = estimate_relaxation_tau(tau_grid, baseline_mode_values)
                    microscopic_tau = float(micro_map[int(n_side)])
                    relaxation_ratio = relaxation_tau / max(float(microscopic_tau), 1.0e-12)
                    relaxation_rows.append(
                        {
                            "hierarchy_label": hierarchy_label,
                            "level_id": f"n{n_side}",
                            "n_side": int(n_side),
                            "h": round_float(1.0 / float(n_side)),
                            "ensemble_size": int(ensemble_size),
                            "seed": int(seed),
                            "collective_relaxation_tau": round_float(relaxation_tau),
                            "microscopic_persistence_tau": round_float(microscopic_tau),
                            "relaxation_ratio": round_float(relaxation_ratio),
                            "initial_density_mode_amplitude": round_float(initial_mode_amplitude),
                            "terminal_density_mode_amplitude": round_float(baseline_mode_values[-1]),
                            "terminal_diffusion_proxy": round_float(baseline_diffusion_values[-1]),
                            "terminal_low_k_power": round_float(baseline_low_k_values[-1]),
                            "terminal_fluctuation_slope": round_float(baseline_slope_values[-1]),
                        }
                    )

                    eos_rows.append(
                        {
                            "hierarchy_label": hierarchy_label,
                            "level_id": f"n{n_side}",
                            "n_side": int(n_side),
                            "h": round_float(1.0 / float(n_side)),
                            "ensemble_size": int(ensemble_size),
                            "seed": int(seed),
                            "initial_nearest_neighbor_mean": round_float(initial_nn_mean),
                            "initial_pair_distance_mean": round_float(initial_pair_mean),
                            "density_proxy": round_float(1.0 / max(initial_nn_mean, 1.0e-12)),
                            "terminal_survival_fraction": round_float(terminal_survival),
                            "terminal_cluster_frequency": round_float(terminal_cluster_frequency),
                            "terminal_spectral_distortion": round_float(terminal_spectral),
                        }
                    )

                    for bias_strength in bias_strengths:
                        driven_states = evolve_with_phase_kicks(
                            ensemble_state0,
                            eigen_grid,
                            points,
                            tau_grid,
                            bias_strength,
                            bias_remove_tau,
                            phase_axis,
                        )
                        eval_snapshot: dict[str, float] | None = None
                        remove_snapshot: dict[str, float] | None = None
                        final_snapshot: dict[str, float] | None = None
                        for tau_value, state_tau in zip(tau_grid, driven_states):
                            candidate_identity = []
                            candidate_local_masses = []
                            survivor_flags = []
                            for index in range(len(centers)):
                                identity_value = windowed_overlap(state_tau, ref_states[index], identity_windows[index])
                                mass_value = local_mass(state_tau, mass_windows[index])
                                mass_ratio = mass_value / max(float(initial_local_masses[index]), 1.0e-12)
                                candidate_identity.append(identity_value)
                                candidate_local_masses.append(mass_value)
                                survivor_flags.append(bool(identity_value >= min_identity and mass_ratio >= min_mass_ratio))
                            mass_array = np.asarray(candidate_local_masses, dtype=float)
                            total_mass = max(float(np.sum(mass_array)), 1.0e-12)
                            centroid = weighted_centroid(centers, mass_array)
                            drift_x = float(periodic_delta(np.asarray([[centroid[0], 0.0]], dtype=float), np.asarray([initial_centroid[0], 0.0], dtype=float))[0, 0])
                            mode_amplitude = float(abs(np.sum(mass_array * np.exp(-2j * math.pi * centers[:, 0]))) / total_mass)
                            mean_identity = float(np.mean(candidate_identity))
                            survival_fraction = float(np.mean(survivor_flags))
                            baseline_low_mode = single_mode_baseline[(hierarchy_label, int(n_side), float(tau_value))]
                            spectral_ratio = low_mode_occupancy(state_tau, eigen_grid, band_multiplier) / max(float(baseline_low_mode), 1.0e-12)
                            snapshot = {
                                "centroid_drift": drift_x,
                                "mode_amplitude": mode_amplitude,
                                "mean_identity": mean_identity,
                                "survival_fraction": survival_fraction,
                                "spectral_distortion": spectral_ratio - 1.0,
                            }
                            if abs(tau_value - evaluation_tau) <= 1.0e-12:
                                eval_snapshot = dict(snapshot)
                            if abs(tau_value - bias_remove_tau) <= 1.0e-12:
                                remove_snapshot = dict(snapshot)
                            if abs(tau_value - tau_grid[-1]) <= 1.0e-12:
                                final_snapshot = dict(snapshot)

                        if eval_snapshot is None or remove_snapshot is None or final_snapshot is None:
                            raise ValueError("Phase XIV density-response snapshots were not captured on the frozen tau grid.")

                        baseline_key = (hierarchy_label, int(n_side), int(ensemble_size), int(seed))
                        baseline_snapshot = baseline_response_reference[baseline_key]
                        centroid_response_eval = float(eval_snapshot["centroid_drift"] - baseline_snapshot["centroid_drift_eval"])
                        mode_response_eval = float(eval_snapshot["mode_amplitude"] - baseline_snapshot["mode_amplitude_eval"])
                        mode_response_remove = float(remove_snapshot["mode_amplitude"] - baseline_snapshot["mode_amplitude_remove"])
                        mode_response_final = float(final_snapshot["mode_amplitude"] - baseline_snapshot["mode_amplitude_final"])
                        memory_ratio = (
                            abs(mode_response_final) / max(abs(mode_response_remove), 1.0e-12)
                            if float(bias_strength) > 0.0
                            else 0.0
                        )
                        response_rows.append(
                            {
                                "hierarchy_label": hierarchy_label,
                                "level_id": f"n{n_side}",
                                "n_side": int(n_side),
                                "h": round_float(1.0 / float(n_side)),
                                "ensemble_size": int(ensemble_size),
                                "seed": int(seed),
                                "bias_strength": round_float(bias_strength, 6),
                                "bias_remove_tau": round_float(bias_remove_tau, 6),
                                "centroid_response_eval": round_float(centroid_response_eval),
                                "mode_response_eval": round_float(mode_response_eval),
                                "mode_response_final": round_float(mode_response_final),
                                "mean_identity_eval": round_float(eval_snapshot["mean_identity"]),
                                "survival_fraction_eval": round_float(eval_snapshot["survival_fraction"]),
                                "spectral_distortion_eval": round_float(eval_snapshot["spectral_distortion"]),
                                "memory_ratio": round_float(memory_ratio),
                                "linearity_error": "",
                            }
                        )

    small_bias = float(sorted([value for value in bias_strengths if value > 0.0])[0])
    large_bias = float(sorted([value for value in bias_strengths if value > 0.0])[-1])
    branch_linearity_error = linearity_error(response_rows, BRANCH_LABEL, small_bias, large_bias)
    control_linearity_error = linearity_error(response_rows, CONTROL_LABEL, small_bias, large_bias)
    for row in response_rows:
        if row["hierarchy_label"] == BRANCH_LABEL:
            row["linearity_error"] = round_float(branch_linearity_error)
        else:
            row["linearity_error"] = round_float(control_linearity_error)

    write_csv_rows(
        TRANSPORT_LEDGER_PATH,
        transport_rows,
        ["hierarchy_label", "level_id", "n_side", "h", "ensemble_size", "seed", "tau", "observable", "value"],
    )
    write_csv_rows(
        RELAXATION_PATH,
        relaxation_rows,
        [
            "hierarchy_label",
            "level_id",
            "n_side",
            "h",
            "ensemble_size",
            "seed",
            "collective_relaxation_tau",
            "microscopic_persistence_tau",
            "relaxation_ratio",
            "initial_density_mode_amplitude",
            "terminal_density_mode_amplitude",
            "terminal_diffusion_proxy",
            "terminal_low_k_power",
            "terminal_fluctuation_slope",
        ],
    )
    write_csv_rows(
        FLUCTUATION_PATH,
        fluctuation_rows,
        [
            "hierarchy_label",
            "level_id",
            "n_side",
            "h",
            "ensemble_size",
            "seed",
            "tau",
            "k_index",
            "power",
            "low_k_power",
            "fluctuation_amplitude",
            "fluctuation_slope",
        ],
    )
    write_csv_rows(
        RESPONSE_PATH,
        response_rows,
        [
            "hierarchy_label",
            "level_id",
            "n_side",
            "h",
            "ensemble_size",
            "seed",
            "bias_strength",
            "bias_remove_tau",
            "centroid_response_eval",
            "mode_response_eval",
            "mode_response_final",
            "mean_identity_eval",
            "survival_fraction_eval",
            "spectral_distortion_eval",
            "memory_ratio",
            "linearity_error",
        ],
    )
    write_csv_rows(
        EOS_PATH,
        eos_rows,
        [
            "hierarchy_label",
            "level_id",
            "n_side",
            "h",
            "ensemble_size",
            "seed",
            "initial_nearest_neighbor_mean",
            "initial_pair_distance_mean",
            "density_proxy",
            "terminal_survival_fraction",
            "terminal_cluster_frequency",
            "terminal_spectral_distortion",
        ],
    )

    plot_centroid_drift(transport_rows, tau_grid, levels, int(ensemble_sizes[-1]))
    plot_variance_growth(transport_rows, tau_grid, levels, int(ensemble_sizes[-1]))
    plot_relaxation_times(relaxation_rows, levels)
    plot_fluctuation_power(fluctuation_rows, levels, int(ensemble_sizes[-1]), float(tau_grid[-1]))
    plot_survival_vs_spacing(eos_rows)

    terminal_transport_rows = [
        row
        for row in transport_rows
        if abs(float(row["tau"]) - float(tau_grid[-1])) <= 1.0e-12
        and row["observable"] in {"effective_diffusion_proxy", "density_redistribution_rate"}
    ]
    terminal_fluct_rows = [
        row
        for row in relaxation_rows
        if True
    ]
    max_bias_rows = [row for row in response_rows if abs(float(row["bias_strength"]) - float(large_bias)) <= 1.0e-12]

    branch_transport_drift = successive_distance_by_seed(
        [row for row in terminal_transport_rows if row["observable"] == "effective_diffusion_proxy"],
        "value",
        BRANCH_LABEL,
        levels,
        ensemble_sizes,
        seeds,
    )
    branch_transport_seed_std = seed_std_max(
        [row for row in terminal_transport_rows if row["observable"] == "effective_diffusion_proxy"],
        "value",
        BRANCH_LABEL,
    )
    control_transport_drift = successive_distance_by_seed(
        [row for row in terminal_transport_rows if row["observable"] == "effective_diffusion_proxy"],
        "value",
        CONTROL_LABEL,
        levels,
        ensemble_sizes,
        seeds,
    )
    control_transport_seed_std = seed_std_max(
        [row for row in terminal_transport_rows if row["observable"] == "effective_diffusion_proxy"],
        "value",
        CONTROL_LABEL,
    )

    branch_relax_mean = float(
        np.mean([float(row["relaxation_ratio"]) for row in relaxation_rows if row["hierarchy_label"] == BRANCH_LABEL])
    )
    control_relax_mean = float(
        np.mean([float(row["relaxation_ratio"]) for row in relaxation_rows if row["hierarchy_label"] == CONTROL_LABEL])
    )
    branch_relax_span = float(
        (
            np.max([float(row["relaxation_ratio"]) for row in relaxation_rows if row["hierarchy_label"] == BRANCH_LABEL])
            - np.min([float(row["relaxation_ratio"]) for row in relaxation_rows if row["hierarchy_label"] == BRANCH_LABEL])
        )
        / max(branch_relax_mean, 1.0e-12)
    )
    control_relax_span = float(
        (
            np.max([float(row["relaxation_ratio"]) for row in relaxation_rows if row["hierarchy_label"] == CONTROL_LABEL])
            - np.min([float(row["relaxation_ratio"]) for row in relaxation_rows if row["hierarchy_label"] == CONTROL_LABEL])
        )
        / max(control_relax_mean, 1.0e-12)
    )
    branch_relax_tau_drift = successive_distance_by_seed(
        relaxation_rows,
        "collective_relaxation_tau",
        BRANCH_LABEL,
        levels,
        ensemble_sizes,
        seeds,
    )
    control_relax_tau_drift = successive_distance_by_seed(
        relaxation_rows,
        "collective_relaxation_tau",
        CONTROL_LABEL,
        levels,
        ensemble_sizes,
        seeds,
    )
    branch_relax_ratio_min = float(
        np.min([float(row["relaxation_ratio"]) for row in relaxation_rows if row["hierarchy_label"] == BRANCH_LABEL])
    )
    control_relax_ratio_min = float(
        np.min([float(row["relaxation_ratio"]) for row in relaxation_rows if row["hierarchy_label"] == CONTROL_LABEL])
    )

    branch_fluct_drift = successive_distance_by_seed(
        relaxation_rows,
        "terminal_low_k_power",
        BRANCH_LABEL,
        levels,
        ensemble_sizes,
        seeds,
    )
    control_fluct_drift = successive_distance_by_seed(
        relaxation_rows,
        "terminal_low_k_power",
        CONTROL_LABEL,
        levels,
        ensemble_sizes,
        seeds,
    )
    branch_fluct_seed_std = seed_std_max(relaxation_rows, "terminal_low_k_power", BRANCH_LABEL)
    control_fluct_seed_std = seed_std_max(relaxation_rows, "terminal_low_k_power", CONTROL_LABEL)
    branch_slope_seed_std = seed_std_max(relaxation_rows, "terminal_fluctuation_slope", BRANCH_LABEL)
    control_slope_seed_std = seed_std_max(relaxation_rows, "terminal_fluctuation_slope", CONTROL_LABEL)
    branch_slope_drift = successive_distance_by_seed(
        relaxation_rows,
        "terminal_fluctuation_slope",
        BRANCH_LABEL,
        levels,
        ensemble_sizes,
        seeds,
    )
    control_slope_drift = successive_distance_by_seed(
        relaxation_rows,
        "terminal_fluctuation_slope",
        CONTROL_LABEL,
        levels,
        ensemble_sizes,
        seeds,
    )

    branch_identity_min_max_bias = float(
        np.min([float(row["mean_identity_eval"]) for row in max_bias_rows if row["hierarchy_label"] == BRANCH_LABEL])
    )
    branch_survival_min_max_bias = float(
        np.min([float(row["survival_fraction_eval"]) for row in max_bias_rows if row["hierarchy_label"] == BRANCH_LABEL])
    )
    control_identity_min_max_bias = float(
        np.min([float(row["mean_identity_eval"]) for row in max_bias_rows if row["hierarchy_label"] == CONTROL_LABEL])
    )
    control_survival_min_max_bias = float(
        np.min([float(row["survival_fraction_eval"]) for row in max_bias_rows if row["hierarchy_label"] == CONTROL_LABEL])
    )

    branch_density_response_mean = float(
        np.mean([float(row["mode_response_eval"]) for row in max_bias_rows if row["hierarchy_label"] == BRANCH_LABEL])
    )
    control_density_response_mean = float(
        np.mean([float(row["mode_response_eval"]) for row in max_bias_rows if row["hierarchy_label"] == CONTROL_LABEL])
    )
    branch_identity_mean_max_bias = float(
        np.mean([float(row["mean_identity_eval"]) for row in max_bias_rows if row["hierarchy_label"] == BRANCH_LABEL])
    )
    control_identity_mean_max_bias = float(
        np.mean([float(row["mean_identity_eval"]) for row in max_bias_rows if row["hierarchy_label"] == CONTROL_LABEL])
    )
    branch_survival_mean_max_bias = float(
        np.mean([float(row["survival_fraction_eval"]) for row in max_bias_rows if row["hierarchy_label"] == BRANCH_LABEL])
    )
    control_survival_mean_max_bias = float(
        np.mean([float(row["survival_fraction_eval"]) for row in max_bias_rows if row["hierarchy_label"] == CONTROL_LABEL])
    )

    def safe_corr(xs: list[float], ys: list[float]) -> float:
        if len(xs) < 2:
            return 0.0
        x_array = np.asarray(xs, dtype=float)
        y_array = np.asarray(ys, dtype=float)
        if np.std(x_array) <= 1.0e-12 or np.std(y_array) <= 1.0e-12:
            return 0.0
        return float(np.corrcoef(x_array, y_array)[0, 1])

    branch_eos_rows = [row for row in eos_rows if row["hierarchy_label"] == BRANCH_LABEL]
    control_eos_rows = [row for row in eos_rows if row["hierarchy_label"] == CONTROL_LABEL]
    branch_spacing_survival_corr = safe_corr(
        [float(row["initial_nearest_neighbor_mean"]) for row in branch_eos_rows],
        [float(row["terminal_survival_fraction"]) for row in branch_eos_rows],
    )
    control_spacing_survival_corr = safe_corr(
        [float(row["initial_nearest_neighbor_mean"]) for row in control_eos_rows],
        [float(row["terminal_survival_fraction"]) for row in control_eos_rows],
    )

    branch_identity_linearity_error = metric_linearity_error(response_rows, BRANCH_LABEL, "mean_identity_eval", small_bias, large_bias)
    control_identity_linearity_error = metric_linearity_error(response_rows, CONTROL_LABEL, "mean_identity_eval", small_bias, large_bias)
    branch_flags = {
        "transport_descriptor_stable": bool(branch_transport_drift <= 0.02 and branch_transport_seed_std <= 0.02),
        "collective_relaxation_slower_than_microscopic": bool(branch_relax_tau_drift <= 1.0 and branch_relax_ratio_min >= 1.1),
        "fluctuation_scaling_consistent": bool(branch_fluct_drift <= 0.02 and branch_slope_drift <= 0.08),
        "weak_bias_response_linear_identity_preserving": bool(
            branch_identity_linearity_error <= 0.9
            and branch_identity_mean_max_bias >= 0.96
            and branch_survival_mean_max_bias >= 0.96
        ),
    }
    control_flags = {
        "transport_descriptor_stable": bool(control_transport_drift <= 0.02 and control_transport_seed_std <= 0.02),
        "collective_relaxation_slower_than_microscopic": bool(control_relax_tau_drift <= 1.0 and control_relax_ratio_min >= 1.1),
        "fluctuation_scaling_consistent": bool(control_fluct_drift <= 0.02 and control_slope_drift <= 0.08),
        "weak_bias_response_linear_identity_preserving": bool(
            control_identity_linearity_error <= 0.9
            and control_identity_mean_max_bias >= 0.96
            and control_survival_mean_max_bias >= 0.96
        ),
    }
    success_flags = dict(branch_flags)
    success_flags["control_hierarchy_different"] = bool(not all(control_flags.values()))
    success = bool(all(success_flags.values()))

    runtime_seconds = round_float(time.perf_counter() - start_time, 6)

    write_json(
        RUNS_PATH,
        {
            "timestamp": timestamp_iso(),
            "phase_name": PHASE_NAME,
            "stage_identifier": STAGE_IDENTIFIER,
            "config_reference": repo_rel(CONFIG_PATH),
            "phase13_manifest_reference": repo_rel(PHASE13_MANIFEST_PATH),
            "phase11_manifest_reference": repo_rel(PHASE11_MANIFEST_PATH),
            "phase8_manifest_reference": repo_rel(PHASE8_MANIFEST_PATH),
            "candidate_label": str(config["candidate_label"]),
            "refinement_levels": levels,
            "ensemble_sizes": ensemble_sizes,
            "layout_seeds": seeds,
            "layout_rule": dict(config["layout_rule"]),
            "evaluation_tau": round_float(evaluation_tau, 6),
            "persistence_tau_grid": [round_float(value, 6) for value in tau_grid],
            "phase8_short_time_window": phase8["operational_short_time_window"]["probe_times"],
            "transport_probe": dict(config["transport_probe"]),
            "response_probe": dict(config["response_probe"]),
            "fluctuation_probe": dict(config["fluctuation_probe"]),
            "deterministic_seed_record": config["deterministic_seed_record"],
            "runtime_seconds": runtime_seconds,
        },
    )

    summary_lines = [
        "# Phase XIV - Collective Dynamics and Mesoscopic Transport",
        "",
        "## Objective",
        "",
        "Test whether the frozen dilute candidate sector supports bounded collective transport, collective relaxation, fluctuation scaling, and weak-bias response without modifying any earlier frozen contract.",
        "",
        "## Frozen Inputs",
        "",
        f"- Candidate: `{config['candidate_label']}`",
        f"- Refinement levels: `{levels}` with `h = 1 / n_side`",
        f"- Ensemble sizes: `{ensemble_sizes}`",
        f"- Layout seeds: `{seeds}` with minimum physical separation `{min_sep}`",
        f"- Tau grid: `{tau_grid}`",
        f"- Phase VIII short-time window retained by contract: `{phase8['operational_short_time_window']['probe_times']}`",
        "",
        "## Key Results",
        "",
        f"- Branch terminal diffusion-proxy fixed-layout drift: `{round_float(branch_transport_drift)}`.",
        f"- Control terminal diffusion-proxy fixed-layout drift: `{round_float(control_transport_drift)}`.",
        f"- Branch mean collective/microscopic relaxation ratio: `{round_float(branch_relax_mean)}` with min ratio `{round_float(branch_relax_ratio_min)}` and fixed-layout tau drift `{round_float(branch_relax_tau_drift)}`.",
        f"- Control mean collective/microscopic relaxation ratio: `{round_float(control_relax_mean)}` with min ratio `{round_float(control_relax_ratio_min)}` and fixed-layout tau drift `{round_float(control_relax_tau_drift)}`.",
        f"- Branch terminal low-k fluctuation drift: `{round_float(branch_fluct_drift)}` and slope drift `{round_float(branch_slope_drift)}`.",
        f"- Control terminal low-k fluctuation drift: `{round_float(control_fluct_drift)}` and slope drift `{round_float(control_slope_drift)}`.",
        f"- Branch weak-bias identity-linearity error: `{round_float(branch_identity_linearity_error)}` with mean identity `{round_float(branch_identity_mean_max_bias)}` and mean survival `{round_float(branch_survival_mean_max_bias)}` under max bias.",
        f"- Control weak-bias identity-linearity error: `{round_float(control_identity_linearity_error)}` with mean identity `{round_float(control_identity_mean_max_bias)}` and mean survival `{round_float(control_survival_mean_max_bias)}` under max bias.",
        f"- Mean maximum-bias mode response (branch, control): `{round_float(branch_density_response_mean)}`, `{round_float(control_density_response_mean)}`.",
        f"- Spacing-survival correlation (branch, control): `{round_float(branch_spacing_survival_corr)}`, `{round_float(control_spacing_survival_corr)}`.",
        "",
        "## Bounded Interpretation",
        "",
        "Phase XIV treats collective feasibility as established only when at least one transport descriptor stays bounded across refinement within each fixed seeded layout family, collective relaxation remains slower than microscopic persistence with stable fixed-layout timescales, fluctuation scaling stays coherent in low-k power and slope, weak bias produces approximately linear aggregate identity loss while preserving mean candidate identity and survival, and the control hierarchy fails at least one of those same tests. These diagnostics do not assert hydrodynamics, thermodynamics, field equations, or continuum medium laws.",
        "",
        closure_statement(success),
    ]
    write_text(SUMMARY_PATH, "\n".join(summary_lines))

    manifest = {
        "timestamp": timestamp_iso(),
        "phase": 14,
        "phase_name": PHASE_NAME,
        "stage_identifier": STAGE_IDENTIFIER,
        "status": "passed" if success else "failed",
        "success": bool(success),
        "objective": "collective_dynamics_feasibility",
        "operator_manifest_reference": repo_rel(PHASE6_MANIFEST_PATH),
        "phase7_manifest_reference": repo_rel(PHASE7_MANIFEST_PATH),
        "phase8_manifest_reference": repo_rel(PHASE8_MANIFEST_PATH),
        "phase9_manifest_reference": repo_rel(PHASE9_MANIFEST_PATH),
        "phase10_manifest_reference": repo_rel(PHASE10_MANIFEST_PATH),
        "phase11_manifest_reference": repo_rel(PHASE11_MANIFEST_PATH),
        "phase12_manifest_reference": repo_rel(PHASE12_MANIFEST_PATH),
        "phase13_manifest_reference": repo_rel(PHASE13_MANIFEST_PATH),
        "candidate_label": str(config["candidate_label"]),
        "refinement_hierarchy": {
            "h_definition": "h = 1 / n_side on the unit-periodic lattice; it is the refinement-order parameter inherited from the frozen operator hierarchy.",
            "levels": levels,
        },
        "evaluation_protocol": {
            "evaluation_tau": round_float(evaluation_tau, 6),
            "persistence_tau_grid": [round_float(value, 6) for value in tau_grid],
            "ensemble_sizes": ensemble_sizes,
            "layout_seeds": seeds,
            "minimum_physical_separation": round_float(min_sep),
            "phase8_short_time_window": phase8["operational_short_time_window"]["probe_times"],
            "transport_probe": dict(config["transport_probe"]),
            "response_probe": dict(config["response_probe"]),
            "fluctuation_probe": dict(config["fluctuation_probe"]),
        },
        "aggregate_metrics": {
            "branch_transport_fixed_layout_drift": round_float(branch_transport_drift),
            "branch_transport_seed_std": round_float(branch_transport_seed_std),
            "control_transport_fixed_layout_drift": round_float(control_transport_drift),
            "control_transport_seed_std": round_float(control_transport_seed_std),
            "branch_relaxation_ratio_mean": round_float(branch_relax_mean),
            "branch_relaxation_ratio_relative_span": round_float(branch_relax_span),
            "branch_relaxation_ratio_min": round_float(branch_relax_ratio_min),
            "branch_relaxation_tau_fixed_layout_drift": round_float(branch_relax_tau_drift),
            "control_relaxation_ratio_mean": round_float(control_relax_mean),
            "control_relaxation_ratio_relative_span": round_float(control_relax_span),
            "control_relaxation_ratio_min": round_float(control_relax_ratio_min),
            "control_relaxation_tau_fixed_layout_drift": round_float(control_relax_tau_drift),
            "branch_low_k_power_drift": round_float(branch_fluct_drift),
            "branch_low_k_power_seed_std": round_float(branch_fluct_seed_std),
            "branch_fluctuation_slope_seed_std": round_float(branch_slope_seed_std),
            "branch_fluctuation_slope_drift": round_float(branch_slope_drift),
            "control_low_k_power_drift": round_float(control_fluct_drift),
            "control_low_k_power_seed_std": round_float(control_fluct_seed_std),
            "control_fluctuation_slope_seed_std": round_float(control_slope_seed_std),
            "control_fluctuation_slope_drift": round_float(control_slope_drift),
            "branch_linearity_error": round_float(branch_linearity_error),
            "control_linearity_error": round_float(control_linearity_error),
            "branch_identity_linearity_error": round_float(branch_identity_linearity_error),
            "control_identity_linearity_error": round_float(control_identity_linearity_error),
            "branch_identity_min_max_bias": round_float(branch_identity_min_max_bias),
            "branch_survival_min_max_bias": round_float(branch_survival_min_max_bias),
            "control_identity_min_max_bias": round_float(control_identity_min_max_bias),
            "control_survival_min_max_bias": round_float(control_survival_min_max_bias),
            "branch_identity_mean_max_bias": round_float(branch_identity_mean_max_bias),
            "branch_survival_mean_max_bias": round_float(branch_survival_mean_max_bias),
            "control_identity_mean_max_bias": round_float(control_identity_mean_max_bias),
            "control_survival_mean_max_bias": round_float(control_survival_mean_max_bias),
            "branch_max_bias_mode_response_mean": round_float(branch_density_response_mean),
            "control_max_bias_mode_response_mean": round_float(control_density_response_mean),
            "branch_spacing_survival_correlation": round_float(branch_spacing_survival_corr),
            "control_spacing_survival_correlation": round_float(control_spacing_survival_corr),
        },
        "success_flags": success_flags,
        "control_flags": control_flags,
        "artifacts": {
            "transport_ledger_csv": repo_rel(TRANSPORT_LEDGER_PATH),
            "relaxation_times_csv": repo_rel(RELAXATION_PATH),
            "fluctuation_spectra_csv": repo_rel(FLUCTUATION_PATH),
            "density_response_csv": repo_rel(RESPONSE_PATH),
            "equation_of_state_proxy_csv": repo_rel(EOS_PATH),
            "runs_json": repo_rel(RUNS_PATH),
            "summary_md": repo_rel(SUMMARY_PATH),
            "manifest_json": repo_rel(MANIFEST_PATH),
            "drift_plot": repo_rel(DRIFT_PLOT_PATH),
            "variance_plot": repo_rel(VARIANCE_PLOT_PATH),
            "relaxation_plot": repo_rel(RELAX_PLOT_PATH),
            "fluctuation_plot": repo_rel(FLUCT_PLOT_PATH),
            "equation_of_state_plot": repo_rel(EOS_PLOT_PATH),
            "builder_script": repo_rel(Path(__file__)),
        },
        "runtime_seconds": runtime_seconds,
        "claim_boundary": CLAIM_BOUNDARY,
        "conclusion": closure_statement(success),
    }
    write_json(MANIFEST_PATH, manifest)


if __name__ == "__main__":
    main()
