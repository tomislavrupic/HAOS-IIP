#!/usr/bin/env python3

from __future__ import annotations

import json
import math
import os
import sys
import time
from collections import Counter
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

PHASE_ROOT = ROOT / "phase13-sector-formation"
CONFIG_PATH = PHASE_ROOT / "configs" / "phase13_sector_config.json"
RUNS_ROOT = PHASE_ROOT / "runs"
PLOTS_ROOT = PHASE_ROOT / "plots"

PHASE6_MANIFEST_PATH = ROOT / "phase6-operator" / "phase6_operator_manifest.json"
PHASE7_MANIFEST_PATH = ROOT / "phase7-spectral" / "phase7_spectral_manifest.json"
PHASE8_MANIFEST_PATH = ROOT / "phase8-trace" / "phase8_trace_manifest.json"
PHASE9_MANIFEST_PATH = ROOT / "phase9-invariants" / "phase9_manifest.json"
PHASE10_MANIFEST_PATH = ROOT / "phase10-bridge" / "phase10_manifest.json"
PHASEX_MANIFEST_PATH = ROOT / "phaseX-proto-particle" / "phaseX_integrated_manifest.json"
PHASE11_MANIFEST_PATH = ROOT / "phase11-protection" / "phase11_manifest.json"
PHASE12_MANIFEST_PATH = ROOT / "phase12-interactions" / "phase12_manifest.json"

SURVIVAL_LEDGER_PATH = RUNS_ROOT / "phase13_population_survival_ledger.csv"
SPACING_LEDGER_PATH = RUNS_ROOT / "phase13_spacing_statistics_ledger.csv"
CLUSTER_LEDGER_PATH = RUNS_ROOT / "phase13_cluster_metrics.csv"
SPECTRAL_LEDGER_PATH = RUNS_ROOT / "phase13_spectral_ensemble_proxy.csv"
RUNS_PATH = RUNS_ROOT / "phase13_runs.json"
SUMMARY_PATH = PHASE_ROOT / "phase13_summary.md"
MANIFEST_PATH = PHASE_ROOT / "phase13_manifest.json"

SURVIVAL_PLOT_PATH = PLOTS_ROOT / "phase13_survival_fraction_vs_time.svg"
PAIR_HIST_PLOT_PATH = PLOTS_ROOT / "phase13_pair_distance_histogram_vs_refinement.svg"
NN_PLOT_PATH = PLOTS_ROOT / "phase13_nearest_neighbor_scale_vs_refinement.svg"
CLUSTER_PLOT_PATH = PLOTS_ROOT / "phase13_cluster_count_vs_time.svg"
SPECTRAL_PLOT_PATH = PLOTS_ROOT / "phase13_ensemble_spectral_proxy_vs_population.svg"

PHASE_NAME = "phase13-sector-formation"
STAGE_IDENTIFIER = "phase13-sector-formation"
BRANCH_LABEL = "frozen_branch"
CONTROL_LABEL = "periodic_diagonal_augmented_control"
CLAIM_BOUNDARY = (
    "Phase XIII is limited to multi-mode statistical sector formation diagnostics on the frozen "
    "localized candidate and frozen operator hierarchy. It does not assert thermodynamics, particles, "
    "gases, kinetic laws, or physical correspondence."
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
        return "Phase XIII establishes multi-mode statistical sector formation feasibility for the frozen operator hierarchy."
    return "Phase XIII does not yet establish multi-mode statistical sector formation feasibility for the frozen operator hierarchy."


def repo_rel(path: Path) -> str:
    return str(Path(path).resolve().relative_to(ROOT))


def periodic_delta(points: np.ndarray, anchor: np.ndarray) -> np.ndarray:
    delta = np.asarray(points, dtype=float) - np.asarray(anchor, dtype=float)
    return (delta + 0.5) % 1.0 - 0.5


def torus_distance(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.sqrt(np.sum(periodic_delta(np.asarray([left], dtype=float), np.asarray(right, dtype=float))[0] ** 2)))


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
    control = base + 4.0 * diagonal_weight * (1.0 - cos_x * cos_y)
    return np.asarray(control, dtype=float)


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
    delta = periodic_delta(points, center)
    squared_radius = np.sum(delta * delta, axis=1)
    return np.exp(-0.5 * squared_radius / (float(sigma) * float(sigma)))


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


def expected_uniform_shell_density(bin_edges: np.ndarray) -> np.ndarray:
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    widths = np.diff(bin_edges)
    shell_areas = 2.0 * math.pi * centers * widths
    total = max(float(np.sum(shell_areas)), 1.0e-12)
    return np.asarray(shell_areas / total, dtype=float)


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


def load_frozen_inputs() -> dict[str, Any]:
    config = read_json(CONFIG_PATH)
    phase6 = read_json(PHASE6_MANIFEST_PATH)
    phase7 = read_json(PHASE7_MANIFEST_PATH)
    phase8 = read_json(PHASE8_MANIFEST_PATH)
    phase9 = read_json(PHASE9_MANIFEST_PATH)
    phase10 = read_json(PHASE10_MANIFEST_PATH)
    phasex = read_json(PHASEX_MANIFEST_PATH)
    phase11 = read_json(PHASE11_MANIFEST_PATH)
    phase12 = read_json(PHASE12_MANIFEST_PATH)

    manifests = [phase6, phase7, phase8, phase9, phase10, phasex, phase11, phase12]
    if not all(bool(manifest.get("success")) for manifest in manifests):
        raise ValueError("Phase XIII requires successful frozen Phases VI, VII, VIII, IX, X, proto-particle, XI, and XII manifests.")

    candidate_label = str(config["candidate_label"])
    if candidate_label != str(phase11["candidate_label"]):
        raise ValueError("Phase XIII candidate must match the frozen Phase XI candidate label.")
    if candidate_label != str(phase12["candidate_label"]):
        raise ValueError("Phase XIII candidate must match the frozen Phase XII candidate label.")
    if candidate_label not in phasex.get("persistence_results", {}).get("branch_refinement_stable_candidates", []):
        raise ValueError("Phase XIII candidate must already be refinement-stable in the frozen proto-particle manifest.")

    tau_grid = [float(value) for value in config["persistence_tau_grid"]]
    if tau_grid != [float(value) for value in phase11["evaluation_protocol"]["persistence_tau_grid"]]:
        raise ValueError("Phase XIII tau grid must match the frozen Phase XI persistence tau grid.")
    if float(config["evaluation_tau"]) != float(phase11["evaluation_protocol"]["evaluation_tau"]):
        raise ValueError("Phase XIII evaluation tau must match the frozen Phase XI evaluation tau.")
    if phase8["operational_short_time_window"]["probe_times"] != phase11["evaluation_protocol"]["short_time_window"]:
        raise ValueError("Phase VIII and XI frozen short-time windows must agree.")

    return {
        "config": config,
        "phase6": phase6,
        "phase7": phase7,
        "phase8": phase8,
        "phase9": phase9,
        "phase10": phase10,
        "phasex": phasex,
        "phase11": phase11,
        "phase12": phase12,
    }


def plot_survival_fraction(survival_rows: list[dict[str, Any]], tau_grid: list[float], ensemble_sizes: list[int], levels: list[int]) -> None:
    plt.figure(figsize=(9.0, 5.6))
    level = int(levels[-1])
    for hierarchy_label, linestyle in ((BRANCH_LABEL, "-"), (CONTROL_LABEL, "--")):
        for ensemble_size in ensemble_sizes:
            rows = [
                row
                for row in survival_rows
                if row["hierarchy_label"] == hierarchy_label
                and int(row["n_side"]) == level
                and int(row["ensemble_size"]) == int(ensemble_size)
            ]
            grouped = []
            for tau_value in tau_grid:
                matches = [float(row["survival_fraction"]) for row in rows if abs(float(row["tau"]) - float(tau_value)) <= 1.0e-12]
                grouped.append(float(np.mean(matches)) if matches else 0.0)
            plt.plot(tau_grid, grouped, linestyle=linestyle, marker="o", label=f"{hierarchy_label}:N{ensemble_size}")
    plt.xlabel("tau")
    plt.ylabel("survival fraction")
    plt.ylim(-0.02, 1.02)
    plt.grid(alpha=0.3)
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(SURVIVAL_PLOT_PATH, format="svg")
    plt.close()


def plot_pair_histogram(spacing_rows: list[dict[str, Any]], ensemble_size: int, seed: int, levels: list[int], bin_centers: np.ndarray) -> None:
    plt.figure(figsize=(8.8, 5.4))
    for hierarchy_label, linestyle in ((BRANCH_LABEL, "-"), (CONTROL_LABEL, "--")):
        for n_side in levels:
            row = next(
                row
                for row in spacing_rows
                if row["hierarchy_label"] == hierarchy_label
                and int(row["n_side"]) == int(n_side)
                and int(row["ensemble_size"]) == int(ensemble_size)
                and int(row["seed"]) == int(seed)
                and abs(float(row["tau"]) - 0.8) <= 1.0e-12
            )
            hist = np.asarray(json.loads(str(row["pair_distance_histogram_density_json"])), dtype=float)
            plt.plot(bin_centers, hist, linestyle=linestyle, marker="o", label=f"{hierarchy_label}:n{n_side}")
    plt.xlabel("pair distance")
    plt.ylabel("normalized pair-distance density")
    plt.grid(alpha=0.3)
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(PAIR_HIST_PLOT_PATH, format="svg")
    plt.close()


def plot_nearest_neighbor_scale(spacing_rows: list[dict[str, Any]], ensemble_sizes: list[int], levels: list[int]) -> None:
    plt.figure(figsize=(8.6, 5.2))
    for hierarchy_label, marker in ((BRANCH_LABEL, "o"), (CONTROL_LABEL, "s")):
        for ensemble_size in ensemble_sizes:
            values = []
            for n_side in levels:
                matches = [
                    float(row["nearest_neighbor_mean"])
                    for row in spacing_rows
                    if row["hierarchy_label"] == hierarchy_label
                    and int(row["n_side"]) == int(n_side)
                    and int(row["ensemble_size"]) == int(ensemble_size)
                    and abs(float(row["tau"]) - 0.8) <= 1.0e-12
                ]
                values.append(float(np.mean(matches)) if matches else 0.0)
            plt.plot(levels, values, marker=marker, label=f"{hierarchy_label}:N{ensemble_size}")
    plt.xlabel("n_side")
    plt.ylabel("nearest-neighbor mean")
    plt.grid(alpha=0.3)
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(NN_PLOT_PATH, format="svg")
    plt.close()


def plot_cluster_count(cluster_rows: list[dict[str, Any]], tau_grid: list[float], levels: list[int], ensemble_size: int) -> None:
    plt.figure(figsize=(8.8, 5.4))
    level = int(levels[-1])
    for hierarchy_label, linestyle in ((BRANCH_LABEL, "-"), (CONTROL_LABEL, "--")):
        rows = [
            row
            for row in cluster_rows
            if row["hierarchy_label"] == hierarchy_label
            and int(row["n_side"]) == level
            and int(row["ensemble_size"]) == int(ensemble_size)
        ]
        means = []
        for tau_value in tau_grid:
            values = [float(row["cluster_count"]) for row in rows if abs(float(row["tau"]) - float(tau_value)) <= 1.0e-12]
            means.append(float(np.mean(values)) if values else 0.0)
        plt.plot(tau_grid, means, linestyle=linestyle, marker="o", label=f"{hierarchy_label}:N{ensemble_size}")
    plt.xlabel("tau")
    plt.ylabel("cluster count")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(CLUSTER_PLOT_PATH, format="svg")
    plt.close()


def plot_spectral_proxy(spectral_rows: list[dict[str, Any]], ensemble_sizes: list[int], levels: list[int]) -> None:
    plt.figure(figsize=(8.6, 5.2))
    level = int(levels[-1])
    for hierarchy_label, marker in ((BRANCH_LABEL, "o"), (CONTROL_LABEL, "s")):
        values = []
        for ensemble_size in ensemble_sizes:
            matches = [
                float(row["spectral_distortion_vs_single"])
                for row in spectral_rows
                if row["hierarchy_label"] == hierarchy_label
                and int(row["n_side"]) == level
                and int(row["ensemble_size"]) == int(ensemble_size)
                and abs(float(row["tau"]) - 0.8) <= 1.0e-12
            ]
            values.append(float(np.mean(matches)) if matches else 0.0)
        plt.plot(ensemble_sizes, values, marker=marker, label=hierarchy_label)
    plt.xlabel("ensemble size")
    plt.ylabel("spectral distortion vs single-mode baseline")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(SPECTRAL_PLOT_PATH, format="svg")
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
    phase12 = payload["phase12"]

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
    cluster_distance = float(config["cluster_thresholds"]["distance_physical"])
    midpoint_identity_min = float(config["cluster_thresholds"]["midpoint_identity_min"])
    identity_ceiling = float(config["cluster_thresholds"]["identity_ceiling"])
    hist_bins = np.asarray(config["spacing_statistics"]["histogram_bins"], dtype=float)
    hist_bin_centers = 0.5 * (hist_bins[:-1] + hist_bins[1:])
    expected_shell = expected_uniform_shell_density(hist_bins)
    band_multiplier = float(config["spectral_identity_proxy"]["low_mode_band_multiplier"])

    if abs(evaluation_tau - float(phase11["evaluation_protocol"]["evaluation_tau"])) > 1.0e-12:
        raise ValueError("Phase XIII evaluation tau must remain identical to the frozen Phase XI value.")
    if phase8["operational_short_time_window"]["probe_times"] != phase11["evaluation_protocol"]["short_time_window"]:
        raise ValueError("Phase VIII and XI short-time windows must stay aligned.")

    layout_map: dict[tuple[int, int], list[np.ndarray]] = {}
    for ensemble_size in ensemble_sizes:
        for seed in seeds:
            layout_map[(int(ensemble_size), int(seed))] = generate_layout(seed, ensemble_size, min_sep, max_attempts)

    survival_rows: list[dict[str, Any]] = []
    spacing_rows: list[dict[str, Any]] = []
    cluster_rows: list[dict[str, Any]] = []
    spectral_rows: list[dict[str, Any]] = []

    single_mode_baseline: dict[tuple[str, int, float], float] = {}

    for hierarchy_label in (BRANCH_LABEL, CONTROL_LABEL):
        for n_side in levels:
            points = build_points(n_side)
            eigen_grid = scalar_eigen_grid(hierarchy_label, n_side, epsilon, diagonal_weight_scale)
            baseline_state = build_candidate_state(points, np.array([0.5, 0.5], dtype=float))
            for tau_value in tau_grid:
                single_mode_baseline[(hierarchy_label, int(n_side), float(tau_value))] = low_mode_occupancy(
                    evolve_state(baseline_state, eigen_grid, float(tau_value)), eigen_grid, band_multiplier
                )

            for ensemble_size in ensemble_sizes:
                for seed in seeds:
                    centers = layout_map[(int(ensemble_size), int(seed))]
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

                    ensemble_state0 = np.sum(np.vstack(ref_states), axis=0)
                    ensemble_state0 = np.asarray(ensemble_state0 / max(np.linalg.norm(ensemble_state0), 1.0e-12), dtype=complex)
                    initial_local_mass = [local_mass(ensemble_state0, window) for window in mass_windows]
                    prior_survival = [True for _ in centers]
                    first_failure_tau: list[float | None] = [None for _ in centers]

                    for tau_value in tau_grid:
                        state_tau = evolve_state(ensemble_state0, eigen_grid, float(tau_value))
                        candidate_identity = []
                        candidate_mass_ratio = []
                        survivor_flags = []
                        for index, center in enumerate(centers):
                            identity_value = windowed_overlap(state_tau, ref_states[index], identity_windows[index])
                            mass_value = local_mass(state_tau, mass_windows[index])
                            mass_ratio = mass_value / max(float(initial_local_mass[index]), 1.0e-12)
                            candidate_identity.append(identity_value)
                            candidate_mass_ratio.append(mass_ratio)
                            survives = bool(identity_value >= min_identity and mass_ratio >= min_mass_ratio)
                            survivor_flags.append(survives)
                            if prior_survival[index] and not survives and first_failure_tau[index] is None:
                                first_failure_tau[index] = float(tau_value)
                        new_failures = sum(
                            1 for old_value, new_value in zip(prior_survival, survivor_flags) if old_value and not new_value
                        )
                        prior_survival = list(survivor_flags)

                        surviving_indices = [index for index, survives in enumerate(survivor_flags) if survives]
                        surviving_centers = [centers[index] for index in surviving_indices]
                        surviving_count = int(len(surviving_indices))
                        survival_fraction = float(surviving_count) / float(ensemble_size)
                        decay_times = [value for value in first_failure_tau if value is not None]
                        mean_decay_time = float(np.mean(decay_times)) if decay_times else ""

                        pair_distances: list[float] = []
                        nearest_neighbor_values: list[float] = []
                        if surviving_count >= 2:
                            for local_index, left_index in enumerate(surviving_indices):
                                distances = []
                                for right_index in surviving_indices:
                                    if left_index == right_index:
                                        continue
                                    distance_value = torus_distance(centers[left_index], centers[right_index])
                                    distances.append(distance_value)
                                    if right_index > left_index:
                                        pair_distances.append(distance_value)
                                nearest_neighbor_values.append(min(distances))
                        pair_hist_counts, _ = np.histogram(pair_distances, bins=hist_bins)
                        pair_hist_density = (
                            pair_hist_counts.astype(float) / max(float(np.sum(pair_hist_counts)), 1.0e-12)
                        )
                        g_r_proxy = pair_hist_density / np.maximum(expected_shell, 1.0e-12)
                        nearest_neighbor_mean = float(np.mean(nearest_neighbor_values)) if nearest_neighbor_values else 0.0
                        nearest_neighbor_std = float(np.std(nearest_neighbor_values)) if nearest_neighbor_values else 0.0
                        effective_exclusion_radius = float(min(nearest_neighbor_values)) if nearest_neighbor_values else 0.0

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
                        cluster_count = int(len(clustered_components))
                        largest_cluster_size = int(max((len(component) for component in clustered_components), default=1 if surviving_count else 0))
                        merged_candidate_count = int(sum(len(component) for component in clustered_components))
                        cluster_frequency = float(merged_candidate_count) / float(ensemble_size) if ensemble_size else 0.0

                        low_mode_value = low_mode_occupancy(state_tau, eigen_grid, band_multiplier)
                        baseline_low_mode = single_mode_baseline[(hierarchy_label, int(n_side), float(tau_value))]
                        spectral_ratio = low_mode_value / max(float(baseline_low_mode), 1.0e-12)
                        spectral_distortion = spectral_ratio - 1.0

                        survival_rows.append(
                            {
                                "hierarchy_label": hierarchy_label,
                                "level_id": f"n{n_side}",
                                "n_side": int(n_side),
                                "h": round_float(1.0 / float(n_side)),
                                "ensemble_size": int(ensemble_size),
                                "seed": int(seed),
                                "tau": round_float(tau_value, 6),
                                "surviving_count": int(surviving_count),
                                "survival_fraction": round_float(survival_fraction),
                                "new_failures": int(new_failures),
                                "cumulative_failures": int(ensemble_size - surviving_count),
                                "mean_decay_event_tau": round_float(mean_decay_time) if mean_decay_time != "" else "",
                            }
                        )
                        spacing_rows.append(
                            {
                                "hierarchy_label": hierarchy_label,
                                "level_id": f"n{n_side}",
                                "n_side": int(n_side),
                                "h": round_float(1.0 / float(n_side)),
                                "ensemble_size": int(ensemble_size),
                                "seed": int(seed),
                                "tau": round_float(tau_value, 6),
                                "surviving_count": int(surviving_count),
                                "pair_distance_histogram_counts_json": json.dumps(pair_hist_counts.tolist()),
                                "pair_distance_histogram_density_json": json.dumps([round_float(value) for value in pair_hist_density.tolist()]),
                                "nearest_neighbor_mean": round_float(nearest_neighbor_mean),
                                "nearest_neighbor_std": round_float(nearest_neighbor_std),
                                "effective_exclusion_radius": round_float(effective_exclusion_radius),
                                "g_r_proxy_json": json.dumps([round_float(value) for value in g_r_proxy.tolist()]),
                            }
                        )
                        cluster_rows.append(
                            {
                                "hierarchy_label": hierarchy_label,
                                "level_id": f"n{n_side}",
                                "n_side": int(n_side),
                                "h": round_float(1.0 / float(n_side)),
                                "ensemble_size": int(ensemble_size),
                                "seed": int(seed),
                                "tau": round_float(tau_value, 6),
                                "cluster_count": int(cluster_count),
                                "largest_cluster_size": int(largest_cluster_size),
                                "merged_candidate_count": int(merged_candidate_count),
                                "cluster_frequency": round_float(cluster_frequency),
                                "effective_exclusion_radius": round_float(effective_exclusion_radius),
                            }
                        )
                        spectral_rows.append(
                            {
                                "hierarchy_label": hierarchy_label,
                                "level_id": f"n{n_side}",
                                "n_side": int(n_side),
                                "h": round_float(1.0 / float(n_side)),
                                "ensemble_size": int(ensemble_size),
                                "seed": int(seed),
                                "tau": round_float(tau_value, 6),
                                "low_mode_occupancy": round_float(low_mode_value),
                                "single_mode_baseline_occupancy": round_float(baseline_low_mode),
                                "low_mode_occupancy_ratio_vs_single": round_float(spectral_ratio),
                                "spectral_distortion_vs_single": round_float(spectral_distortion),
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
            "ensemble_size",
            "seed",
            "tau",
            "surviving_count",
            "survival_fraction",
            "new_failures",
            "cumulative_failures",
            "mean_decay_event_tau",
        ],
    )
    write_csv_rows(
        SPACING_LEDGER_PATH,
        spacing_rows,
        [
            "hierarchy_label",
            "level_id",
            "n_side",
            "h",
            "ensemble_size",
            "seed",
            "tau",
            "surviving_count",
            "pair_distance_histogram_counts_json",
            "pair_distance_histogram_density_json",
            "nearest_neighbor_mean",
            "nearest_neighbor_std",
            "effective_exclusion_radius",
            "g_r_proxy_json",
        ],
    )
    write_csv_rows(
        CLUSTER_LEDGER_PATH,
        cluster_rows,
        [
            "hierarchy_label",
            "level_id",
            "n_side",
            "h",
            "ensemble_size",
            "seed",
            "tau",
            "cluster_count",
            "largest_cluster_size",
            "merged_candidate_count",
            "cluster_frequency",
            "effective_exclusion_radius",
        ],
    )
    write_csv_rows(
        SPECTRAL_LEDGER_PATH,
        spectral_rows,
        [
            "hierarchy_label",
            "level_id",
            "n_side",
            "h",
            "ensemble_size",
            "seed",
            "tau",
            "low_mode_occupancy",
            "single_mode_baseline_occupancy",
            "low_mode_occupancy_ratio_vs_single",
            "spectral_distortion_vs_single",
        ],
    )

    plot_survival_fraction(survival_rows, tau_grid, ensemble_sizes, levels)
    plot_pair_histogram(spacing_rows, int(ensemble_sizes[-1]), int(seeds[0]), levels, hist_bin_centers)
    plot_nearest_neighbor_scale(spacing_rows, ensemble_sizes, levels)
    plot_cluster_count(cluster_rows, tau_grid, levels, int(ensemble_sizes[-1]))
    plot_spectral_proxy(spectral_rows, ensemble_sizes, levels)

    evaluation_survival_rows = [row for row in survival_rows if abs(float(row["tau"]) - evaluation_tau) <= 1.0e-12]
    evaluation_spacing_rows = [row for row in spacing_rows if abs(float(row["tau"]) - evaluation_tau) <= 1.0e-12]
    evaluation_cluster_rows = [row for row in cluster_rows if abs(float(row["tau"]) - evaluation_tau) <= 1.0e-12]
    evaluation_spectral_rows = [row for row in spectral_rows if abs(float(row["tau"]) - evaluation_tau) <= 1.0e-12]

    def mean_metric(rows: list[dict[str, Any]], key: str, hierarchy_label: str, ensemble_size: int) -> float:
        values = [float(row[key]) for row in rows if row["hierarchy_label"] == hierarchy_label and int(row["ensemble_size"]) == int(ensemble_size)]
        return float(np.mean(values)) if values else 0.0

    def seed_std_max(rows: list[dict[str, Any]], key: str, hierarchy_label: str) -> float:
        grouped: dict[tuple[int, int], list[float]] = {}
        for row in rows:
            if row["hierarchy_label"] != hierarchy_label:
                continue
            grouped.setdefault((int(row["n_side"]), int(row["ensemble_size"])), []).append(float(row[key]))
        if not grouped:
            return 0.0
        return float(max(np.std(values) for values in grouped.values()))

    def successive_distance(rows: list[dict[str, Any]], key: str, hierarchy_label: str) -> float:
        distances = []
        for ensemble_size in ensemble_sizes:
            means = []
            for n_side in levels:
                values = [
                    float(row[key])
                    for row in rows
                    if row["hierarchy_label"] == hierarchy_label
                    and int(row["ensemble_size"]) == int(ensemble_size)
                    and int(row["n_side"]) == int(n_side)
                ]
                means.append(float(np.mean(values)) if values else 0.0)
            for left, right in zip(means, means[1:]):
                distances.append(abs(right - left))
        return float(max(distances)) if distances else 0.0

    def successive_distance_by_seed(rows: list[dict[str, Any]], key: str, hierarchy_label: str) -> float:
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

    def max_histogram_l1_distance(rows: list[dict[str, Any]], hierarchy_label: str) -> float:
        distances = []
        for ensemble_size in ensemble_sizes:
            for seed in seeds:
                histograms = []
                for n_side in levels:
                    row = next(
                        (
                            row
                            for row in rows
                            if row["hierarchy_label"] == hierarchy_label
                            and int(row["ensemble_size"]) == int(ensemble_size)
                            and int(row["seed"]) == int(seed)
                            and int(row["n_side"]) == int(n_side)
                        ),
                        None,
                    )
                    if row is None:
                        continue
                    histograms.append(np.asarray(json.loads(str(row["pair_distance_histogram_density_json"])), dtype=float))
                for left, right in zip(histograms, histograms[1:]):
                    distances.append(float(np.sum(np.abs(right - left))))
        return float(max(distances)) if distances else 0.0

    def max_tau_gap(rows: list[dict[str, Any]], key: str, left_label: str, right_label: str) -> float:
        distances = []
        for tau_value in tau_grid:
            left_values = [
                float(row[key])
                for row in rows
                if row["hierarchy_label"] == left_label and abs(float(row["tau"]) - float(tau_value)) <= 1.0e-12
            ]
            right_values = [
                float(row[key])
                for row in rows
                if row["hierarchy_label"] == right_label and abs(float(row["tau"]) - float(tau_value)) <= 1.0e-12
            ]
            if not left_values or not right_values:
                continue
            distances.append(float(np.mean(left_values) - np.mean(right_values)))
        return float(max(distances)) if distances else 0.0

    def relative_span(values: list[float]) -> float:
        if not values:
            return float("inf")
        array = np.asarray(values, dtype=float)
        return float((np.max(array) - np.min(array)) / max(abs(float(np.mean(array))), 1.0e-12))

    branch_seed_std_survival = seed_std_max(evaluation_survival_rows, "survival_fraction", BRANCH_LABEL)
    branch_seed_std_spacing = seed_std_max(evaluation_spacing_rows, "nearest_neighbor_mean", BRANCH_LABEL)
    branch_seed_std_cluster = seed_std_max(evaluation_cluster_rows, "cluster_frequency", BRANCH_LABEL)
    branch_successive_survival = successive_distance(evaluation_survival_rows, "survival_fraction", BRANCH_LABEL)
    branch_successive_spacing = successive_distance(evaluation_spacing_rows, "nearest_neighbor_mean", BRANCH_LABEL)
    branch_successive_spacing_by_seed = successive_distance_by_seed(evaluation_spacing_rows, "nearest_neighbor_mean", BRANCH_LABEL)
    branch_successive_cluster_by_seed = successive_distance_by_seed(evaluation_cluster_rows, "cluster_frequency", BRANCH_LABEL)
    branch_successive_exclusion_by_seed = successive_distance_by_seed(evaluation_spacing_rows, "effective_exclusion_radius", BRANCH_LABEL)
    branch_histogram_l1 = max_histogram_l1_distance(evaluation_spacing_rows, BRANCH_LABEL)

    branch_nn_values = [float(row["nearest_neighbor_mean"]) for row in evaluation_spacing_rows if row["hierarchy_label"] == BRANCH_LABEL]
    control_nn_values = [float(row["nearest_neighbor_mean"]) for row in evaluation_spacing_rows if row["hierarchy_label"] == CONTROL_LABEL]
    branch_exclusion_values = [float(row["effective_exclusion_radius"]) for row in evaluation_spacing_rows if row["hierarchy_label"] == BRANCH_LABEL]
    control_exclusion_values = [float(row["effective_exclusion_radius"]) for row in evaluation_spacing_rows if row["hierarchy_label"] == CONTROL_LABEL]

    branch_survival_mean = float(np.mean([float(row["survival_fraction"]) for row in evaluation_survival_rows if row["hierarchy_label"] == BRANCH_LABEL]))
    control_survival_mean = float(np.mean([float(row["survival_fraction"]) for row in evaluation_survival_rows if row["hierarchy_label"] == CONTROL_LABEL]))
    branch_cluster_mean = float(np.mean([float(row["cluster_frequency"]) for row in evaluation_cluster_rows if row["hierarchy_label"] == BRANCH_LABEL]))
    control_cluster_mean = float(np.mean([float(row["cluster_frequency"]) for row in evaluation_cluster_rows if row["hierarchy_label"] == CONTROL_LABEL]))
    branch_spectral_mean = float(
        np.mean([float(row["spectral_distortion_vs_single"]) for row in evaluation_spectral_rows if row["hierarchy_label"] == BRANCH_LABEL])
    )
    control_spectral_mean = float(
        np.mean([float(row["spectral_distortion_vs_single"]) for row in evaluation_spectral_rows if row["hierarchy_label"] == CONTROL_LABEL])
    )
    branch_nn_mean = float(np.mean(branch_nn_values)) if branch_nn_values else 0.0
    control_nn_mean = float(np.mean(control_nn_values)) if control_nn_values else 0.0
    max_tau_survival_gap = max_tau_gap(survival_rows, "survival_fraction", BRANCH_LABEL, CONTROL_LABEL)
    max_tau_spectral_gap = max_tau_gap(spectral_rows, "spectral_distortion_vs_single", BRANCH_LABEL, CONTROL_LABEL)
    terminal_tau = float(tau_grid[-1])
    terminal_branch_survival = float(
        np.mean(
            [
                float(row["survival_fraction"])
                for row in survival_rows
                if row["hierarchy_label"] == BRANCH_LABEL and abs(float(row["tau"]) - terminal_tau) <= 1.0e-12
            ]
        )
    )
    terminal_control_survival = float(
        np.mean(
            [
                float(row["survival_fraction"])
                for row in survival_rows
                if row["hierarchy_label"] == CONTROL_LABEL and abs(float(row["tau"]) - terminal_tau) <= 1.0e-12
            ]
        )
    )

    success_flags = {
        "ensemble_survival_reproducible": bool(branch_seed_std_survival <= 0.08 and branch_successive_survival <= 0.15),
        "spacing_bands_compact": bool(branch_successive_spacing_by_seed <= 0.01 and branch_histogram_l1 <= 0.05),
        "cluster_behavior_stable": bool(branch_successive_cluster_by_seed <= 0.05 and branch_successive_exclusion_by_seed <= 0.05),
        "control_hierarchy_different": bool(
            max_tau_survival_gap >= 0.05
            or max_tau_spectral_gap >= 0.04
            or (terminal_branch_survival - terminal_control_survival) >= 0.05
        ),
    }
    success = bool(all(success_flags.values()))

    runtime_seconds = round_float(time.perf_counter() - start_time, 6)

    write_json(
        RUNS_PATH,
        {
            "timestamp": timestamp_iso(),
            "phase_name": PHASE_NAME,
            "stage_identifier": STAGE_IDENTIFIER,
            "config_reference": repo_rel(CONFIG_PATH),
            "operator_manifest_reference": repo_rel(PHASE6_MANIFEST_PATH),
            "phase8_manifest_reference": repo_rel(PHASE8_MANIFEST_PATH),
            "phase11_manifest_reference": repo_rel(PHASE11_MANIFEST_PATH),
            "phase12_manifest_reference": repo_rel(PHASE12_MANIFEST_PATH),
            "candidate_label": str(config["candidate_label"]),
            "refinement_levels": levels,
            "ensemble_sizes": ensemble_sizes,
            "layout_seeds": seeds,
            "layout_rule": dict(config["layout_rule"]),
            "evaluation_tau": round_float(evaluation_tau, 6),
            "persistence_tau_grid": [round_float(value, 6) for value in tau_grid],
            "phase8_short_time_window": phase8["operational_short_time_window"]["probe_times"],
            "deterministic_seed_record": config["deterministic_seed_record"],
            "runtime_seconds": runtime_seconds,
        },
    )

    summary_lines = [
        "# Phase XIII - Multi-Mode Statistical Sector Formation",
        "",
        "## Objective",
        "",
        "Test whether dilute populations of the frozen persistent localized candidate preserve approximate identity, develop reproducible spacing structure, and remain distinguishable from the deterministic control hierarchy.",
        "",
        "## Frozen Inputs",
        "",
        f"- Candidate: `{config['candidate_label']}` from Phase XI and Phase XII",
        f"- Refinement levels: `{levels}` with `h = 1 / n_side`",
        f"- Ensemble sizes: `{ensemble_sizes}`",
        f"- Layout seeds: `{seeds}` with minimum physical separation `{min_sep}`",
        f"- Evaluation tau: `{evaluation_tau}`",
        f"- Persistence tau grid: `{tau_grid}`",
        f"- Phase VIII short-time window retained by contract: `{phase8['operational_short_time_window']['probe_times']}`",
        "",
        "## Key Results",
        "",
        f"- Mean branch evaluation survival fraction: `{round_float(branch_survival_mean)}`.",
        f"- Mean control evaluation survival fraction: `{round_float(control_survival_mean)}`.",
        f"- Mean branch evaluation cluster frequency: `{round_float(branch_cluster_mean)}`.",
        f"- Mean control evaluation cluster frequency: `{round_float(control_cluster_mean)}`.",
        f"- Mean branch evaluation spectral distortion: `{round_float(branch_spectral_mean)}`.",
        f"- Mean control evaluation spectral distortion: `{round_float(control_spectral_mean)}`.",
        f"- Mean branch nearest-neighbor scale: `{round_float(branch_nn_mean)}`.",
        f"- Mean control nearest-neighbor scale: `{round_float(control_nn_mean)}`.",
        f"- Branch seed std (survival, spacing, cluster): `{round_float(branch_seed_std_survival)}`, `{round_float(branch_seed_std_spacing)}`, `{round_float(branch_seed_std_cluster)}`.",
        f"- Branch successive-refinement distances (survival, spacing): `{round_float(branch_successive_survival)}`, `{round_float(branch_successive_spacing)}`.",
        f"- Branch refinement drift within fixed seeded layouts (spacing, cluster, exclusion, pair-hist L1): `{round_float(branch_successive_spacing_by_seed)}`, `{round_float(branch_successive_cluster_by_seed)}`, `{round_float(branch_successive_exclusion_by_seed)}`, `{round_float(branch_histogram_l1)}`.",
        f"- Maximum branch-control gaps across the tau grid (survival, spectral distortion): `{round_float(max_tau_survival_gap)}`, `{round_float(max_tau_spectral_gap)}`.",
        f"- Terminal tau `{round_float(terminal_tau, 6)}` survival fractions (branch, control): `{round_float(terminal_branch_survival)}`, `{round_float(terminal_control_survival)}`.",
        "",
        "## Bounded Interpretation",
        "",
        "The branch is treated as feasible only if each deterministic layout family remains compact across refinement, and if the same seeded layouts on the control hierarchy separate under time-resolved survival or spectral occupancy diagnostics. Seed-to-seed spacing variation is treated as deterministic geometry variation, not as refinement instability. These results do not assert thermodynamics, particle ensembles, or kinetic laws.",
        "",
        closure_statement(success),
    ]
    write_text(SUMMARY_PATH, "\n".join(summary_lines))

    manifest = {
        "timestamp": timestamp_iso(),
        "phase": 13,
        "phase_name": PHASE_NAME,
        "stage_identifier": STAGE_IDENTIFIER,
        "status": "passed" if success else "failed",
        "success": bool(success),
        "objective": "multi_mode_statistical_sector_formation_feasibility",
        "operator_manifest_reference": repo_rel(PHASE6_MANIFEST_PATH),
        "phase7_manifest_reference": repo_rel(PHASE7_MANIFEST_PATH),
        "phase8_manifest_reference": repo_rel(PHASE8_MANIFEST_PATH),
        "phase9_manifest_reference": repo_rel(PHASE9_MANIFEST_PATH),
        "phase10_manifest_reference": repo_rel(PHASE10_MANIFEST_PATH),
        "phasex_manifest_reference": repo_rel(PHASEX_MANIFEST_PATH),
        "phase11_manifest_reference": repo_rel(PHASE11_MANIFEST_PATH),
        "phase12_manifest_reference": repo_rel(PHASE12_MANIFEST_PATH),
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
        },
        "evaluation_means": {
            "branch_survival_fraction": round_float(branch_survival_mean),
            "control_survival_fraction": round_float(control_survival_mean),
            "branch_cluster_frequency": round_float(branch_cluster_mean),
            "control_cluster_frequency": round_float(control_cluster_mean),
            "branch_spectral_distortion": round_float(branch_spectral_mean),
            "control_spectral_distortion": round_float(control_spectral_mean),
            "branch_nearest_neighbor_mean": round_float(branch_nn_mean),
            "control_nearest_neighbor_mean": round_float(control_nn_mean),
        },
        "stability_metrics": {
            "branch_seed_std_survival": round_float(branch_seed_std_survival),
            "branch_seed_std_spacing": round_float(branch_seed_std_spacing),
            "branch_seed_std_cluster": round_float(branch_seed_std_cluster),
            "branch_successive_refinement_survival_distance": round_float(branch_successive_survival),
            "branch_successive_refinement_spacing_distance": round_float(branch_successive_spacing),
            "branch_successive_refinement_spacing_distance_fixed_layout": round_float(branch_successive_spacing_by_seed),
            "branch_successive_refinement_cluster_distance_fixed_layout": round_float(branch_successive_cluster_by_seed),
            "branch_successive_refinement_exclusion_distance_fixed_layout": round_float(branch_successive_exclusion_by_seed),
            "branch_pair_histogram_l1_distance_fixed_layout": round_float(branch_histogram_l1),
            "branch_nearest_neighbor_relative_span": round_float(relative_span(branch_nn_values)),
            "branch_exclusion_radius_relative_span": round_float(relative_span(branch_exclusion_values)),
            "control_nearest_neighbor_relative_span": round_float(relative_span(control_nn_values)),
            "control_exclusion_radius_relative_span": round_float(relative_span(control_exclusion_values)),
            "max_tau_survival_gap_branch_minus_control": round_float(max_tau_survival_gap),
            "max_tau_spectral_gap_branch_minus_control": round_float(max_tau_spectral_gap),
            "terminal_tau": round_float(terminal_tau, 6),
            "terminal_branch_survival_fraction": round_float(terminal_branch_survival),
            "terminal_control_survival_fraction": round_float(terminal_control_survival),
        },
        "success_flags": success_flags,
        "artifacts": {
            "population_survival_ledger_csv": repo_rel(SURVIVAL_LEDGER_PATH),
            "spacing_statistics_ledger_csv": repo_rel(SPACING_LEDGER_PATH),
            "cluster_metrics_csv": repo_rel(CLUSTER_LEDGER_PATH),
            "spectral_ensemble_proxy_csv": repo_rel(SPECTRAL_LEDGER_PATH),
            "summary_md": repo_rel(SUMMARY_PATH),
            "manifest_json": repo_rel(MANIFEST_PATH),
            "survival_plot": repo_rel(SURVIVAL_PLOT_PATH),
            "pair_hist_plot": repo_rel(PAIR_HIST_PLOT_PATH),
            "nearest_neighbor_plot": repo_rel(NN_PLOT_PATH),
            "cluster_plot": repo_rel(CLUSTER_PLOT_PATH),
            "spectral_plot": repo_rel(SPECTRAL_PLOT_PATH),
            "runs_json": repo_rel(RUNS_PATH),
            "builder_script": repo_rel(Path(__file__)),
        },
        "runtime_seconds": runtime_seconds,
        "claim_boundary": CLAIM_BOUNDARY,
        "conclusion": closure_statement(success),
    }
    write_json(MANIFEST_PATH, manifest)


if __name__ == "__main__":
    main()
