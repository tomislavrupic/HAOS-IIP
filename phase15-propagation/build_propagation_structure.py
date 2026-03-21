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

PHASE_ROOT = ROOT / "phase15-propagation"
CONFIG_PATH = PHASE_ROOT / "configs" / "phase15_propagation_config.json"
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
PHASE13_CONFIG_PATH = ROOT / "phase13-sector-formation" / "configs" / "phase13_sector_config.json"
PHASE14_CONFIG_PATH = ROOT / "phase14-collective-dynamics" / "configs" / "phase14_collective_config.json"

PROPAGATION_LEDGER_PATH = RUNS_ROOT / "phase15_propagation_ledger.csv"
SPEED_LEDGER_PATH = RUNS_ROOT / "phase15_effective_speed_ledger.csv"
INFLUENCE_LEDGER_PATH = RUNS_ROOT / "phase15_influence_range_ledger.csv"
TRANSPORT_LEDGER_PATH = RUNS_ROOT / "phase15_transport_descriptor_ledger.csv"
RUNS_PATH = RUNS_ROOT / "phase15_runs.json"
SUMMARY_PATH = PHASE_ROOT / "phase15_summary.md"
MANIFEST_PATH = PHASE_ROOT / "phase15_manifest.json"

RADIUS_PLOT_PATH = PLOTS_ROOT / "phase15_disturbance_radius_vs_time.svg"
SPEED_PLOT_PATH = PLOTS_ROOT / "phase15_effective_speed_vs_refinement.svg"
INFLUENCE_PLOT_PATH = PLOTS_ROOT / "phase15_influence_range_envelope_vs_time.svg"
LOWK_PLOT_PATH = PLOTS_ROOT / "phase15_low_k_transport_descriptor_vs_refinement.svg"

PHASE_NAME = "phase15-propagation"
STAGE_IDENTIFIER = "phase15-propagation"
BRANCH_LABEL = "frozen_branch"
CONTROL_LABEL = "periodic_diagonal_augmented_control"
CLAIM_BOUNDARY = (
    "Phase XV is limited to effective propagation-structure feasibility diagnostics on the frozen "
    "collective sector and frozen operator hierarchy. It does not assert metric structure, "
    "relativistic light-cones, spacetime emergence, continuum field equations, or physical correspondence."
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
        return "Phase XV establishes effective propagation-structure feasibility for the frozen operator hierarchy."
    return "Phase XV does not yet establish effective propagation-structure feasibility for the frozen operator hierarchy."


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
    total_power = float(np.sum(spectrum[1:])) if spectrum.size > 1 else 0.0
    return {
        "fluctuation_amplitude": amplitude,
        "low_k_power": low_k_power,
        "fluctuation_slope": slope,
        "low_k_fraction": low_k_power / max(total_power, 1.0e-12),
        "power_values": [float(value) for value in spectrum[1:].tolist()],
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


def successive_distance_by_seed(
    rows: list[dict[str, Any]],
    key: str,
    hierarchy_label: str,
    probe_name: str,
    levels: list[int],
    ensemble_sizes: list[int],
    seeds: list[int],
) -> float:
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
                        and row["probe_name"] == probe_name
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


def seed_std_max(
    rows: list[dict[str, Any]],
    key: str,
    hierarchy_label: str,
    probe_name: str,
) -> float:
    grouped: dict[tuple[int, int], list[float]] = {}
    for row in rows:
        if row["hierarchy_label"] != hierarchy_label or row["probe_name"] != probe_name:
            continue
        grouped.setdefault((int(row["n_side"]), int(row["ensemble_size"])), []).append(float(row[key]))
    if not grouped:
        return 0.0
    return float(max(np.std(values) for values in grouped.values()))


def influence_order_score(distances: np.ndarray, latencies: np.ndarray) -> float:
    pairs = 0
    violations = 0
    for left in range(len(distances)):
        for right in range(len(distances)):
            if distances[left] + 1.0e-12 < distances[right]:
                pairs += 1
                if latencies[left] > latencies[right] + 1.0e-12:
                    violations += 1
    if pairs == 0:
        return 1.0
    return float(1.0 - violations / pairs)


def latency_shell_envelope(distances: np.ndarray, latencies: np.ndarray) -> tuple[list[float], list[float]]:
    unique_distances = sorted({round(float(value), 6) for value in distances if float(value) > 1.0e-12})
    shell_distances = []
    shell_latencies = []
    for distance_key in unique_distances:
        shell_values = [float(latency) for distance, latency in zip(distances, latencies) if abs(float(distance) - float(distance_key)) <= 1.0e-6]
        if not shell_values:
            continue
        shell_distances.append(float(distance_key))
        shell_latencies.append(float(np.mean(shell_values)))
    return shell_distances, shell_latencies


def fit_effective_speed(distances: np.ndarray, latencies: np.ndarray) -> float:
    mask = (latencies > 1.0e-12) & (distances > 1.0e-12)
    if np.count_nonzero(mask) == 0:
        return 0.0
    speed_values = np.asarray(distances[mask] / latencies[mask], dtype=float)
    return float(np.mean(speed_values))


def plot_radius(propagation_rows: list[dict[str, Any]], tau_grid: list[float], levels: list[int], ensemble_size: int) -> None:
    plt.figure(figsize=(8.8, 5.4))
    level = int(levels[-1])
    for hierarchy_label, linestyle in ((BRANCH_LABEL, "-"), (CONTROL_LABEL, "--")):
        matches = [
            row
            for row in propagation_rows
            if row["probe_name"] == "density_pulse"
            and row["hierarchy_label"] == hierarchy_label
            and int(row["n_side"]) == level
            and int(row["ensemble_size"]) == int(ensemble_size)
            and row["observable"] == "disturbance_radius"
        ]
        series = []
        for tau_value in tau_grid:
            values = [float(row["value"]) for row in matches if abs(float(row["tau"]) - float(tau_value)) <= 1.0e-12]
            series.append(float(np.mean(values)) if values else 0.0)
        plt.plot(tau_grid, series, linestyle=linestyle, marker="o", label=hierarchy_label)
    plt.xlabel("tau")
    plt.ylabel("disturbance radius")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(RADIUS_PLOT_PATH, format="svg")
    plt.close()


def plot_speed(speed_rows: list[dict[str, Any]], levels: list[int]) -> None:
    plt.figure(figsize=(8.6, 5.2))
    for hierarchy_label, marker in ((BRANCH_LABEL, "o"), (CONTROL_LABEL, "s")):
        means = []
        for n_side in levels:
            values = [
                float(row["effective_speed"])
                for row in speed_rows
                if row["probe_name"] == "density_pulse"
                and row["hierarchy_label"] == hierarchy_label
                and int(row["n_side"]) == int(n_side)
            ]
            means.append(float(np.mean(values)) if values else 0.0)
        plt.plot(levels, means, marker=marker, label=hierarchy_label)
    plt.xlabel("n_side")
    plt.ylabel("effective propagation speed")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(SPEED_PLOT_PATH, format="svg")
    plt.close()


def plot_influence(influence_rows: list[dict[str, Any]], levels: list[int], ensemble_size: int) -> None:
    plt.figure(figsize=(8.8, 5.4))
    level = int(levels[-1])
    for hierarchy_label, linestyle in ((BRANCH_LABEL, "-"), (CONTROL_LABEL, "--")):
        rows = [
            row
            for row in influence_rows
            if row["probe_name"] == "density_pulse"
            and row["hierarchy_label"] == hierarchy_label
            and int(row["n_side"]) == level
            and int(row["ensemble_size"]) == int(ensemble_size)
        ]
        grouped: dict[float, list[float]] = {}
        for row in rows:
            grouped.setdefault(float(row["distance"]), []).append(float(row["latency"]))
        xs = sorted(grouped)
        ys = [float(np.mean(grouped[value])) for value in xs]
        plt.plot(xs, ys, linestyle=linestyle, marker="o", label=hierarchy_label)
    plt.xlabel("distance from anchor")
    plt.ylabel("first-response latency")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(INFLUENCE_PLOT_PATH, format="svg")
    plt.close()


def plot_low_k(transport_rows: list[dict[str, Any]], levels: list[int]) -> None:
    plt.figure(figsize=(8.6, 5.2))
    for hierarchy_label, marker in ((BRANCH_LABEL, "o"), (CONTROL_LABEL, "s")):
        means = []
        for n_side in levels:
            values = [
                float(row["terminal_low_k_fraction"])
                for row in transport_rows
                if row["probe_name"] == "density_pulse"
                and row["hierarchy_label"] == hierarchy_label
                and int(row["n_side"]) == int(n_side)
            ]
            means.append(float(np.mean(values)) if values else 0.0)
        plt.plot(levels, means, marker=marker, label=hierarchy_label)
    plt.xlabel("n_side")
    plt.ylabel("terminal low-k response fraction")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(LOWK_PLOT_PATH, format="svg")
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
    phase13_config = read_json(PHASE13_CONFIG_PATH)
    phase14_config = read_json(PHASE14_CONFIG_PATH)

    manifests = [phase6, phase7, phase8, phase9, phase10, phase11, phase12, phase13, phase14]
    if not all(bool(manifest.get("success")) for manifest in manifests):
        raise ValueError("Phase XV requires successful frozen Phases VI through XIV.")

    candidate_label = str(config["candidate_label"])
    for manifest in (phase11, phase12, phase13, phase14):
        if candidate_label != str(manifest["candidate_label"]):
            raise ValueError("Phase XV candidate label must match the frozen earlier candidate label.")

    if [int(value) for value in config["refinement_levels"]] != [int(value) for value in phase13["refinement_hierarchy"]["levels"]]:
        raise ValueError("Phase XV refinement levels must match the frozen Phase XIII refinement hierarchy.")
    if [int(value) for value in config["ensemble_sizes"]] != [int(value) for value in phase13["evaluation_protocol"]["ensemble_sizes"]]:
        raise ValueError("Phase XV ensemble sizes must match the frozen Phase XIII ensemble sizes.")
    if [int(value) for value in config["layout_seeds"]] != [int(value) for value in phase13["evaluation_protocol"]["layout_seeds"]]:
        raise ValueError("Phase XV layout seeds must match the frozen Phase XIII layout seeds.")
    if [float(value) for value in config["persistence_tau_grid"]] != [float(value) for value in phase14["evaluation_protocol"]["persistence_tau_grid"]]:
        raise ValueError("Phase XV tau grid must match the frozen Phase XIV tau grid.")
    if float(config["evaluation_tau"]) != float(phase14["evaluation_protocol"]["evaluation_tau"]):
        raise ValueError("Phase XV evaluation tau must match the frozen Phase XIV evaluation tau.")
    if float(config["layout_rule"]["min_separation_physical"]) != float(phase13["evaluation_protocol"]["minimum_physical_separation"]):
        raise ValueError("Phase XV minimum separation must match the frozen Phase XIII layout rule.")
    if phase8["operational_short_time_window"]["probe_times"] != phase14["evaluation_protocol"]["phase8_short_time_window"]:
        raise ValueError("Phase XV must reuse the frozen Phase VIII short-time window inherited by Phase XIV.")
    if float(config["control_graph"]["diagonal_weight_scale"]) != float(phase14_config["control_graph"]["diagonal_weight_scale"]):
        raise ValueError("Phase XV control hierarchy must reuse the frozen Phase XIV control graph weight scale.")
    if float(config["transport_probe"]["density_gradient_strength"]) != float(phase14_config["transport_probe"]["density_gradient_strength"]):
        raise ValueError("Phase XV density gradient must match the frozen Phase XIV transport probe.")
    if float(config["disturbance_probes"]["bias_remove_tau"]) != float(phase14_config["response_probe"]["bias_remove_tau"]):
        raise ValueError("Phase XV bias removal tau must match the frozen Phase XIV response probe.")
    if str(config["disturbance_probes"]["phase_kick_axis"]) != str(phase14_config["response_probe"]["phase_kick_axis"]):
        raise ValueError("Phase XV phase-kick axis must match the frozen Phase XIV response probe.")

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
        "phase13_config": phase13_config,
        "phase14_config": phase14_config,
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
    evaluation_tau = float(config["evaluation_tau"])
    min_sep = float(config["layout_rule"]["min_separation_physical"])
    max_attempts = int(config["layout_rule"]["max_attempts"])
    identity_sigma = float(config["candidate_windows"]["identity_sigma"])
    mass_sigma = float(config["candidate_windows"]["local_mass_sigma"])
    min_identity = float(config["survival_thresholds"]["min_identity"])
    min_mass_ratio = float(config["survival_thresholds"]["min_local_mass_ratio"])
    band_multiplier = float(config["spectral_identity_proxy"]["low_mode_band_multiplier"])
    density_gradient_strength = float(config["transport_probe"]["density_gradient_strength"])
    response_threshold_fraction = float(config["transport_probe"]["response_threshold_fraction"])
    x_bin_count = int(config["transport_probe"]["x_bin_count"])
    pulse_amplitude = float(config["disturbance_probes"]["density_pulse_amplitude"])
    removal_amplitude = float(config["disturbance_probes"]["candidate_removal_amplitude"])
    bias_strength = float(config["disturbance_probes"]["bias_strength"])
    bias_remove_tau = float(config["disturbance_probes"]["bias_remove_tau"])
    phase_axis = str(config["disturbance_probes"]["phase_kick_axis"])

    layout_map: dict[tuple[int, int], np.ndarray] = {}
    for ensemble_size in ensemble_sizes:
        for seed in seeds:
            layout_map[(int(ensemble_size), int(seed))] = np.asarray(
                generate_layout(seed, ensemble_size, min_sep, max_attempts),
                dtype=float,
            )

    propagation_rows: list[dict[str, Any]] = []
    speed_rows: list[dict[str, Any]] = []
    influence_rows: list[dict[str, Any]] = []
    transport_rows: list[dict[str, Any]] = []

    probe_names = ("density_pulse", "candidate_removal", "bias_onset")
    branch_gate_rows: dict[str, list[dict[str, Any]]] = {name: [] for name in probe_names}
    control_gate_rows: dict[str, list[dict[str, Any]]] = {name: [] for name in probe_names}

    for hierarchy_label in (BRANCH_LABEL, CONTROL_LABEL):
        for n_side in levels:
            points = build_points(n_side)
            eigen_grid = scalar_eigen_grid(hierarchy_label, n_side, epsilon, diagonal_weight_scale)

            for ensemble_size in ensemble_sizes:
                for seed in seeds:
                    centers = np.asarray(layout_map[(int(ensemble_size), int(seed))], dtype=float)
                    ref_states = [build_candidate_state(points, center) for center in centers]
                    identity_windows = [gaussian_window(points, center, identity_sigma) for center in centers]
                    mass_windows = [gaussian_window(points, center, mass_sigma) for center in centers]

                    anchor_index = int(np.argmin(centers[:, 0]))
                    anchor_center = np.asarray(centers[anchor_index], dtype=float)
                    distances = np.asarray([torus_distance(center, anchor_center) for center in centers], dtype=float)

                    baseline_state0, base_weights = initial_collective_state(ref_states, centers, density_gradient_strength)
                    baseline_states = [evolve_state(baseline_state0, eigen_grid, tau_value) for tau_value in tau_grid]
                    baseline_local_masses = []
                    baseline_survival = []
                    baseline_identity = []
                    baseline_profiles = []
                    baseline_spectral = []
                    baseline_widths = []
                    for state_tau in baseline_states:
                        local_values = np.asarray([local_mass(state_tau, window) for window in mass_windows], dtype=float)
                        identities = np.asarray(
                            [windowed_overlap(state_tau, ref_state, window) for ref_state, window in zip(ref_states, identity_windows)],
                            dtype=float,
                        )
                        mass_ratios = np.asarray(
                            [
                                local_values[index] / max(local_mass(baseline_state0, mass_windows[index]), 1.0e-12)
                                for index in range(len(centers))
                            ],
                            dtype=float,
                        )
                        survivors = np.asarray(
                            [bool(identities[index] >= min_identity and mass_ratios[index] >= min_mass_ratio) for index in range(len(centers))],
                            dtype=bool,
                        )
                        baseline_local_masses.append(local_values)
                        baseline_survival.append(float(np.mean(survivors)))
                        baseline_identity.append(float(np.mean(identities)))
                        baseline_profiles.append(density_profile_x(centers, local_values, x_bin_count))
                        baseline_spectral.append(low_mode_occupancy(state_tau, eigen_grid, band_multiplier))
                        baseline_widths.append(spatial_metrics(state_tau, points)["localization_width"])

                    probe_state0_map = {
                        "density_pulse": build_disturbed_state(
                            "density_pulse",
                            base_weights,
                            ref_states,
                            anchor_index,
                            pulse_amplitude,
                            removal_amplitude,
                        ),
                        "candidate_removal": build_disturbed_state(
                            "candidate_removal",
                            base_weights,
                            ref_states,
                            anchor_index,
                            pulse_amplitude,
                            removal_amplitude,
                        ),
                    }
                    driven_states = evolve_with_bias(
                        baseline_state0,
                        eigen_grid,
                        points,
                        tau_grid,
                        bias_strength,
                        bias_remove_tau,
                        phase_axis,
                    )

                    for probe_name in probe_names:
                        if probe_name == "bias_onset":
                            disturbed_states = driven_states
                        else:
                            disturbed_states = [evolve_state(probe_state0_map[probe_name], eigen_grid, tau_value) for tau_value in tau_grid]

                        snapshots: list[dict[str, Any]] = []
                        max_response_value = 0.0
                        for tau_value, state_tau, base_local, base_profile, base_spectral, base_width, base_survival_value, base_identity_value in zip(
                            tau_grid,
                            disturbed_states,
                            baseline_local_masses,
                            baseline_profiles,
                            baseline_spectral,
                            baseline_widths,
                            baseline_survival,
                            baseline_identity,
                        ):
                            disturbed_local = np.asarray([local_mass(state_tau, window) for window in mass_windows], dtype=float)
                            disturbed_identity = np.asarray(
                                [windowed_overlap(state_tau, ref_state, window) for ref_state, window in zip(ref_states, identity_windows)],
                                dtype=float,
                            )
                            disturbed_mass_ratios = np.asarray(
                                [
                                    disturbed_local[index] / max(local_mass(baseline_state0, mass_windows[index]), 1.0e-12)
                                    for index in range(len(centers))
                                ],
                                dtype=float,
                            )
                            disturbed_survivors = np.asarray(
                                [
                                    bool(disturbed_identity[index] >= min_identity and disturbed_mass_ratios[index] >= min_mass_ratio)
                                    for index in range(len(centers))
                                ],
                                dtype=bool,
                            )

                            response_vector = np.asarray(np.abs(disturbed_local - base_local), dtype=float)
                            response_profile = response_profile_x(centers, base_local, disturbed_local, x_bin_count)
                            response_fluctuation = fluctuation_summary(response_profile)
                            response_dispersion = float(np.sum(response_vector * distances) / max(np.sum(response_vector), 1.0e-12))
                            disturbed_spectral = low_mode_occupancy(state_tau, eigen_grid, band_multiplier)
                            spectral_shift = float(disturbed_spectral - base_spectral)
                            disturbed_width = spatial_metrics(state_tau, points)["localization_width"]
                            width_shift = float(disturbed_width - base_width)
                            survival_shift = float(np.mean(disturbed_survivors) - base_survival_value)
                            identity_shift = float(np.mean(disturbed_identity) - base_identity_value)

                            max_response_value = max(max_response_value, float(np.max(response_vector)))
                            snapshots.append(
                                {
                                    "tau": float(tau_value),
                                    "response_vector": response_vector,
                                    "response_profile": response_profile,
                                    "response_dispersion": response_dispersion,
                                    "spectral_shift": spectral_shift,
                                    "width_shift": width_shift,
                                    "survival_shift": survival_shift,
                                    "identity_shift": identity_shift,
                                    "low_k_response_fraction": float(response_fluctuation["low_k_fraction"]),
                                }
                            )

                        threshold_value = float(response_threshold_fraction * max(max_response_value, 1.0e-12))
                        response_vectors = [np.asarray(snapshot["response_vector"], dtype=float) for snapshot in snapshots]
                        response_profiles = [np.asarray(snapshot["response_profile"], dtype=float) for snapshot in snapshots]

                        for snapshot in snapshots:
                            response_vector = np.asarray(snapshot["response_vector"], dtype=float)
                            active = response_vector >= threshold_value - 1.0e-12
                            active_distances = distances[active]
                            disturbance_radius = float(np.max(active_distances)) if active_distances.size else 0.0
                            active_count = int(np.count_nonzero(active))
                            front_sharpness = float(np.max(response_vector) / max(np.mean(response_vector), 1.0e-12))

                            for observable, value in (
                                ("disturbance_radius", disturbance_radius),
                                ("response_dispersion", snapshot["response_dispersion"]),
                                ("front_sharpness", front_sharpness),
                                ("active_candidate_fraction", active_count / float(len(centers))),
                                ("survival_shift", snapshot["survival_shift"]),
                                ("identity_shift", snapshot["identity_shift"]),
                                ("spectral_shift", snapshot["spectral_shift"]),
                                ("width_shift", snapshot["width_shift"]),
                                ("low_k_response_fraction", snapshot["low_k_response_fraction"]),
                            ):
                                propagation_rows.append(
                                    {
                                        "probe_name": probe_name,
                                        "hierarchy_label": hierarchy_label,
                                        "level_id": f"n{n_side}",
                                        "n_side": int(n_side),
                                        "h": round_float(1.0 / float(n_side)),
                                        "ensemble_size": int(ensemble_size),
                                        "seed": int(seed),
                                        "tau": round_float(snapshot["tau"], 6),
                                        "observable": observable,
                                        "value": round_float(value),
                                    }
                                )

                        latencies = np.full(len(centers), float(tau_grid[-1]), dtype=float)
                        for index in range(len(centers)):
                            found = False
                            for tau_value, response_vector in zip(tau_grid, response_vectors):
                                if float(response_vector[index]) >= threshold_value - 1.0e-12:
                                    latencies[index] = float(tau_value)
                                    found = True
                                    break
                            if not found:
                                latencies[index] = float(tau_grid[-1])

                        shell_distances, shell_latencies = latency_shell_envelope(distances, latencies)
                        effective_speed = fit_effective_speed(distances, latencies)
                        influence_order = influence_order_score(distances, latencies)
                        terminal_profile = response_profiles[-1]
                        terminal_fluctuation = fluctuation_summary(terminal_profile)
                        transport_character = float(
                            np.mean(
                                [
                                    float(row["value"])
                                    for row in propagation_rows[-(len(tau_grid) * 9):]
                                    if row["observable"] == "disturbance_radius"
                                ]
                            )
                            / max(np.mean(latencies[distances > 1.0e-12]) if np.count_nonzero(distances > 1.0e-12) else 1.0, 1.0e-12)
                        )
                        terminal_radius = float(
                            next(
                                float(row["value"])
                                for row in reversed(propagation_rows)
                                if row["probe_name"] == probe_name
                                and row["hierarchy_label"] == hierarchy_label
                                and int(row["n_side"]) == int(n_side)
                                and int(row["ensemble_size"]) == int(ensemble_size)
                                and int(row["seed"]) == int(seed)
                                and row["observable"] == "disturbance_radius"
                                and abs(float(row["tau"]) - float(tau_grid[-1])) <= 1.0e-12
                            )
                        )

                        speed_row = {
                            "probe_name": probe_name,
                            "hierarchy_label": hierarchy_label,
                            "level_id": f"n{n_side}",
                            "n_side": int(n_side),
                            "h": round_float(1.0 / float(n_side)),
                            "ensemble_size": int(ensemble_size),
                            "seed": int(seed),
                            "anchor_index": int(anchor_index),
                            "effective_speed": round_float(effective_speed),
                            "terminal_radius": round_float(terminal_radius),
                            "mean_latency": round_float(float(np.mean(latencies[distances > 1.0e-12])) if np.count_nonzero(distances > 1.0e-12) else 0.0),
                            "influence_order_score": round_float(influence_order),
                            "response_threshold": round_float(threshold_value),
                        }
                        speed_rows.append(speed_row)

                        for distance_value, latency_value in zip(distances, latencies):
                            influence_rows.append(
                                {
                                    "probe_name": probe_name,
                                    "hierarchy_label": hierarchy_label,
                                    "level_id": f"n{n_side}",
                                    "n_side": int(n_side),
                                    "h": round_float(1.0 / float(n_side)),
                                    "ensemble_size": int(ensemble_size),
                                    "seed": int(seed),
                                    "distance": round_float(distance_value),
                                    "latency": round_float(latency_value),
                                    "normalized_latency": round_float(latency_value / max(tau_grid[-1], 1.0e-12)),
                                    "anchor_index": int(anchor_index),
                                }
                            )

                        transport_row = {
                            "probe_name": probe_name,
                            "hierarchy_label": hierarchy_label,
                            "level_id": f"n{n_side}",
                            "n_side": int(n_side),
                            "h": round_float(1.0 / float(n_side)),
                            "ensemble_size": int(ensemble_size),
                            "seed": int(seed),
                            "terminal_low_k_fraction": round_float(terminal_fluctuation["low_k_fraction"]),
                            "terminal_low_k_power": round_float(terminal_fluctuation["low_k_power"]),
                            "terminal_fluctuation_slope": round_float(terminal_fluctuation["fluctuation_slope"]),
                            "transport_character_index": round_float(transport_character),
                            "shell_distance_count": int(len(shell_distances)),
                            "shell_latency_span": round_float(max(shell_latencies) - min(shell_latencies) if shell_latencies else 0.0),
                        }
                        transport_rows.append(transport_row)

                        target_rows = branch_gate_rows if hierarchy_label == BRANCH_LABEL else control_gate_rows
                        target_rows[probe_name].append(speed_row | transport_row)

    branch_speed_drift = successive_distance_by_seed(speed_rows, "effective_speed", BRANCH_LABEL, "density_pulse", levels, ensemble_sizes, seeds)
    control_speed_drift = successive_distance_by_seed(speed_rows, "effective_speed", CONTROL_LABEL, "density_pulse", levels, ensemble_sizes, seeds)
    branch_speed_seed_std = seed_std_max(speed_rows, "effective_speed", BRANCH_LABEL, "density_pulse")
    control_speed_seed_std = seed_std_max(speed_rows, "effective_speed", CONTROL_LABEL, "density_pulse")
    branch_order_drift = successive_distance_by_seed(speed_rows, "influence_order_score", BRANCH_LABEL, "density_pulse", levels, ensemble_sizes, seeds)
    control_order_drift = successive_distance_by_seed(speed_rows, "influence_order_score", CONTROL_LABEL, "density_pulse", levels, ensemble_sizes, seeds)
    branch_low_k_drift = successive_distance_by_seed(transport_rows, "terminal_low_k_fraction", BRANCH_LABEL, "density_pulse", levels, ensemble_sizes, seeds)
    control_low_k_drift = successive_distance_by_seed(transport_rows, "terminal_low_k_fraction", CONTROL_LABEL, "density_pulse", levels, ensemble_sizes, seeds)

    branch_speed_mean = float(
        np.mean([float(row["effective_speed"]) for row in speed_rows if row["probe_name"] == "density_pulse" and row["hierarchy_label"] == BRANCH_LABEL])
    )
    control_speed_mean = float(
        np.mean([float(row["effective_speed"]) for row in speed_rows if row["probe_name"] == "density_pulse" and row["hierarchy_label"] == CONTROL_LABEL])
    )
    branch_order_mean = float(
        np.mean([float(row["influence_order_score"]) for row in speed_rows if row["probe_name"] == "density_pulse" and row["hierarchy_label"] == BRANCH_LABEL])
    )
    control_order_mean = float(
        np.mean([float(row["influence_order_score"]) for row in speed_rows if row["probe_name"] == "density_pulse" and row["hierarchy_label"] == CONTROL_LABEL])
    )
    branch_low_k_mean = float(
        np.mean([float(row["terminal_low_k_fraction"]) for row in transport_rows if row["probe_name"] == "density_pulse" and row["hierarchy_label"] == BRANCH_LABEL])
    )
    control_low_k_mean = float(
        np.mean([float(row["terminal_low_k_fraction"]) for row in transport_rows if row["probe_name"] == "density_pulse" and row["hierarchy_label"] == CONTROL_LABEL])
    )
    branch_shell_span_mean = float(
        np.mean([float(row["shell_latency_span"]) for row in transport_rows if row["probe_name"] == "density_pulse" and row["hierarchy_label"] == BRANCH_LABEL])
    )
    control_shell_span_mean = float(
        np.mean([float(row["shell_latency_span"]) for row in transport_rows if row["probe_name"] == "density_pulse" and row["hierarchy_label"] == CONTROL_LABEL])
    )

    probe_metrics: dict[str, dict[str, float]] = {}
    for probe_name in probe_names:
        branch_probe_speed = [
            float(row["effective_speed"])
            for row in speed_rows
            if row["probe_name"] == probe_name and row["hierarchy_label"] == BRANCH_LABEL
        ]
        control_probe_speed = [
            float(row["effective_speed"])
            for row in speed_rows
            if row["probe_name"] == probe_name and row["hierarchy_label"] == CONTROL_LABEL
        ]
        branch_probe_order = [
            float(row["influence_order_score"])
            for row in speed_rows
            if row["probe_name"] == probe_name and row["hierarchy_label"] == BRANCH_LABEL
        ]
        control_probe_order = [
            float(row["influence_order_score"])
            for row in speed_rows
            if row["probe_name"] == probe_name and row["hierarchy_label"] == CONTROL_LABEL
        ]
        branch_probe_low_k = [
            float(row["terminal_low_k_fraction"])
            for row in transport_rows
            if row["probe_name"] == probe_name and row["hierarchy_label"] == BRANCH_LABEL
        ]
        control_probe_low_k = [
            float(row["terminal_low_k_fraction"])
            for row in transport_rows
            if row["probe_name"] == probe_name and row["hierarchy_label"] == CONTROL_LABEL
        ]
        probe_metrics[probe_name] = {
            "branch_speed_mean": round_float(np.mean(branch_probe_speed) if branch_probe_speed else 0.0),
            "control_speed_mean": round_float(np.mean(control_probe_speed) if control_probe_speed else 0.0),
            "branch_order_mean": round_float(np.mean(branch_probe_order) if branch_probe_order else 0.0),
            "control_order_mean": round_float(np.mean(control_probe_order) if control_probe_order else 0.0),
            "branch_low_k_mean": round_float(np.mean(branch_probe_low_k) if branch_probe_low_k else 0.0),
            "control_low_k_mean": round_float(np.mean(control_probe_low_k) if control_probe_low_k else 0.0),
        }

    success_flags = {
        "propagation_drift_bounded": bool(branch_speed_drift <= 0.03),
        "effective_speed_band_compact": bool(branch_speed_seed_std <= 0.05),
        "influence_ordering_stable": bool(branch_order_mean >= 0.8 and branch_order_drift <= 0.2),
        "low_k_transport_coherent": bool(branch_low_k_drift <= 0.02),
    }
    control_flags = {
        "propagation_drift_bounded": bool(control_speed_drift <= 0.03),
        "effective_speed_band_compact": bool(control_speed_seed_std <= 0.05),
        "influence_ordering_stable": bool(control_order_mean >= 0.8 and control_order_drift <= 0.2),
        "low_k_transport_coherent": bool(control_low_k_drift <= 0.02),
    }
    success_flags["control_hierarchy_different"] = bool(not all(control_flags.values()))
    success = all(success_flags.values())

    write_csv_rows(
        PROPAGATION_LEDGER_PATH,
        propagation_rows,
        ["probe_name", "hierarchy_label", "level_id", "n_side", "h", "ensemble_size", "seed", "tau", "observable", "value"],
    )
    write_csv_rows(
        SPEED_LEDGER_PATH,
        speed_rows,
        [
            "probe_name",
            "hierarchy_label",
            "level_id",
            "n_side",
            "h",
            "ensemble_size",
            "seed",
            "anchor_index",
            "effective_speed",
            "terminal_radius",
            "mean_latency",
            "influence_order_score",
            "response_threshold",
        ],
    )
    write_csv_rows(
        INFLUENCE_LEDGER_PATH,
        influence_rows,
        [
            "probe_name",
            "hierarchy_label",
            "level_id",
            "n_side",
            "h",
            "ensemble_size",
            "seed",
            "distance",
            "latency",
            "normalized_latency",
            "anchor_index",
        ],
    )
    write_csv_rows(
        TRANSPORT_LEDGER_PATH,
        transport_rows,
        [
            "probe_name",
            "hierarchy_label",
            "level_id",
            "n_side",
            "h",
            "ensemble_size",
            "seed",
            "terminal_low_k_fraction",
            "terminal_low_k_power",
            "terminal_fluctuation_slope",
            "transport_character_index",
            "shell_distance_count",
            "shell_latency_span",
        ],
    )

    plot_radius(propagation_rows, tau_grid, levels, ensemble_sizes[-1])
    plot_speed(speed_rows, levels)
    plot_influence(influence_rows, levels, ensemble_sizes[-1])
    plot_low_k(transport_rows, levels)

    runtime_seconds = float(time.perf_counter() - start_time)

    runs_payload = {
        "timestamp": timestamp_iso(),
        "phase": 15,
        "stage_identifier": STAGE_IDENTIFIER,
        "deterministic_seed_record": config["deterministic_seed_record"],
        "refinement_levels": levels,
        "ensemble_sizes": ensemble_sizes,
        "layout_seeds": seeds,
        "tau_grid": tau_grid,
        "probe_names": list(probe_names),
        "control_graph": config["control_graph"],
        "disturbance_probes": config["disturbance_probes"],
        "runtime_seconds": round_float(runtime_seconds, 6),
    }
    write_json(RUNS_PATH, runs_payload)

    summary_lines = [
        "# Phase XV Summary",
        "",
        "Phase XV tests whether the frozen dilute collective sector supports bounded propagation-level structure under deterministic disturbances.",
        "",
        "## Frozen Inputs",
        "",
        f"- Candidate label: `{config['candidate_label']}`",
        f"- Refinement levels: `{levels}`",
        f"- Ensemble sizes: `{ensemble_sizes}`",
        f"- Layout seeds: `{seeds}`",
        f"- Tau grid: `{tau_grid}`",
        "",
        "## Aggregate Metrics",
        "",
        f"- Branch effective-speed drift: `{round_float(branch_speed_drift)}`",
        f"- Control effective-speed drift: `{round_float(control_speed_drift)}`",
        f"- Branch mean effective speed: `{round_float(branch_speed_mean)}`",
        f"- Control mean effective speed: `{round_float(control_speed_mean)}`",
        f"- Branch influence-order mean: `{round_float(branch_order_mean)}`",
        f"- Control influence-order mean: `{round_float(control_order_mean)}`",
        f"- Branch low-k drift: `{round_float(branch_low_k_drift)}`",
        f"- Control low-k drift: `{round_float(control_low_k_drift)}`",
        f"- Branch shell-latency mean span: `{round_float(branch_shell_span_mean)}`",
        f"- Control shell-latency mean span: `{round_float(control_shell_span_mean)}`",
        "",
        "## Probe Means",
        "",
        f"- Density pulse branch/control speed means: `{probe_metrics['density_pulse']['branch_speed_mean']}` / `{probe_metrics['density_pulse']['control_speed_mean']}`",
        f"- Candidate removal branch/control speed means: `{probe_metrics['candidate_removal']['branch_speed_mean']}` / `{probe_metrics['candidate_removal']['control_speed_mean']}`",
        f"- Bias onset branch/control speed means: `{probe_metrics['bias_onset']['branch_speed_mean']}` / `{probe_metrics['bias_onset']['control_speed_mean']}`",
        "",
        "## Success Flags",
        "",
        f"- Propagation drift bounded: `{success_flags['propagation_drift_bounded']}`",
        f"- Effective speed band compact: `{success_flags['effective_speed_band_compact']}`",
        f"- Influence ordering stable: `{success_flags['influence_ordering_stable']}`",
        f"- Low-k transport coherent: `{success_flags['low_k_transport_coherent']}`",
        f"- Control hierarchy different: `{success_flags['control_hierarchy_different']}`",
        "",
        "## Control Flags",
        "",
        f"- Control propagation drift bounded: `{control_flags['propagation_drift_bounded']}`",
        f"- Control effective speed band compact: `{control_flags['effective_speed_band_compact']}`",
        f"- Control influence ordering stable: `{control_flags['influence_ordering_stable']}`",
        f"- Control low-k transport coherent: `{control_flags['low_k_transport_coherent']}`",
        "",
        closure_statement(success),
    ]
    write_text(SUMMARY_PATH, "\n".join(summary_lines))

    manifest = {
        "timestamp": timestamp_iso(),
        "phase": 15,
        "phase_name": PHASE_NAME,
        "stage_identifier": STAGE_IDENTIFIER,
        "status": "passed" if success else "failed",
        "success": bool(success),
        "objective": "effective_propagation_structure_feasibility",
        "operator_manifest_reference": repo_rel(PHASE6_MANIFEST_PATH),
        "phase7_manifest_reference": repo_rel(PHASE7_MANIFEST_PATH),
        "phase8_manifest_reference": repo_rel(PHASE8_MANIFEST_PATH),
        "phase9_manifest_reference": repo_rel(PHASE9_MANIFEST_PATH),
        "phase10_manifest_reference": repo_rel(PHASE10_MANIFEST_PATH),
        "phase11_manifest_reference": repo_rel(PHASE11_MANIFEST_PATH),
        "phase12_manifest_reference": repo_rel(PHASE12_MANIFEST_PATH),
        "phase13_manifest_reference": repo_rel(PHASE13_MANIFEST_PATH),
        "phase14_manifest_reference": repo_rel(PHASE14_MANIFEST_PATH),
        "candidate_label": config["candidate_label"],
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
            "response_threshold_fraction": round_float(response_threshold_fraction),
            "transport_probe": config["transport_probe"],
            "disturbance_probes": config["disturbance_probes"],
        },
        "aggregate_metrics": {
            "branch_effective_speed_mean": round_float(branch_speed_mean),
            "control_effective_speed_mean": round_float(control_speed_mean),
            "branch_effective_speed_drift": round_float(branch_speed_drift),
            "control_effective_speed_drift": round_float(control_speed_drift),
            "branch_effective_speed_seed_std": round_float(branch_speed_seed_std),
            "control_effective_speed_seed_std": round_float(control_speed_seed_std),
            "branch_influence_order_mean": round_float(branch_order_mean),
            "control_influence_order_mean": round_float(control_order_mean),
            "branch_influence_order_drift": round_float(branch_order_drift),
            "control_influence_order_drift": round_float(control_order_drift),
            "branch_low_k_transport_mean": round_float(branch_low_k_mean),
            "control_low_k_transport_mean": round_float(control_low_k_mean),
            "branch_low_k_transport_drift": round_float(branch_low_k_drift),
            "control_low_k_transport_drift": round_float(control_low_k_drift),
            "branch_shell_latency_span_mean": round_float(branch_shell_span_mean),
            "control_shell_latency_span_mean": round_float(control_shell_span_mean),
        },
        "probe_metrics": probe_metrics,
        "success_flags": success_flags,
        "control_flags": control_flags,
        "artifacts": {
            "propagation_ledger_csv": repo_rel(PROPAGATION_LEDGER_PATH),
            "effective_speed_ledger_csv": repo_rel(SPEED_LEDGER_PATH),
            "influence_range_ledger_csv": repo_rel(INFLUENCE_LEDGER_PATH),
            "transport_descriptor_ledger_csv": repo_rel(TRANSPORT_LEDGER_PATH),
            "runs_json": repo_rel(RUNS_PATH),
            "summary_md": repo_rel(SUMMARY_PATH),
            "manifest_json": repo_rel(MANIFEST_PATH),
            "radius_plot": repo_rel(RADIUS_PLOT_PATH),
            "speed_plot": repo_rel(SPEED_PLOT_PATH),
            "influence_plot": repo_rel(INFLUENCE_PLOT_PATH),
            "low_k_plot": repo_rel(LOWK_PLOT_PATH),
            "builder_script": repo_rel(Path(__file__)),
        },
        "runtime_seconds": round_float(runtime_seconds, 6),
        "claim_boundary": CLAIM_BOUNDARY,
        "conclusion": closure_statement(success),
    }
    write_json(MANIFEST_PATH, manifest)


if __name__ == "__main__":
    main()
