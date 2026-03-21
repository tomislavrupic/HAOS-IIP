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

PHASE_ROOT = ROOT / "phase12-interactions"
CONFIG_PATH = PHASE_ROOT / "configs" / "phase12_interaction_config.json"
RUNS_ROOT = PHASE_ROOT / "runs"
PLOTS_ROOT = PHASE_ROOT / "plots"

PHASE6_MANIFEST_PATH = ROOT / "phase6-operator" / "phase6_operator_manifest.json"
PHASE7_MANIFEST_PATH = ROOT / "phase7-spectral" / "phase7_spectral_manifest.json"
PHASE8_MANIFEST_PATH = ROOT / "phase8-trace" / "phase8_trace_manifest.json"
PHASE9_MANIFEST_PATH = ROOT / "phase9-invariants" / "phase9_manifest.json"
PHASE10_MANIFEST_PATH = ROOT / "phase10-bridge" / "phase10_manifest.json"
PHASEX_MANIFEST_PATH = ROOT / "phaseX-proto-particle" / "phaseX_integrated_manifest.json"
PHASE11_MANIFEST_PATH = ROOT / "phase11-protection" / "phase11_manifest.json"

OUTCOME_LEDGER_PATH = RUNS_ROOT / "phase12_interaction_outcome_ledger.csv"
IDENTITY_LEDGER_PATH = RUNS_ROOT / "phase12_identity_metrics_ledger.csv"
THRESHOLD_LEDGER_PATH = RUNS_ROOT / "phase12_interaction_thresholds.csv"
RUNS_PATH = RUNS_ROOT / "phase12_runs.json"
SUMMARY_PATH = PHASE_ROOT / "phase12_summary.md"
MANIFEST_PATH = PHASE_ROOT / "phase12_manifest.json"

PERSISTENCE_PLOT_PATH = PLOTS_ROOT / "phase12_persistence_time_vs_separation.svg"
IDENTITY_PLOT_PATH = PLOTS_ROOT / "phase12_identity_metric_vs_time.svg"
REGIME_MAP_PLOT_PATH = PLOTS_ROOT / "phase12_interaction_regime_map.svg"
THRESHOLD_PLOT_PATH = PLOTS_ROOT / "phase12_threshold_scaling_vs_refinement.svg"

PHASE_NAME = "phase12-interactions"
STAGE_IDENTIFIER = "phase12-interactions"
BRANCH_LABEL = "frozen_branch"
CONTROL_LABEL = "periodic_diagonal_augmented_control"
OUTCOME_CODE = {
    "survival": 0,
    "merger": 1,
    "decoherence": 2,
    "annihilation_like_decay": 3,
}
OUTCOME_COLOR = {
    "survival": "#2b8a3e",
    "merger": "#f08c00",
    "decoherence": "#c2255c",
    "annihilation_like_decay": "#495057",
}
CLAIM_BOUNDARY = (
    "Phase XII is limited to two-mode interaction and identity-preservation diagnostics on the frozen "
    "operator hierarchy and frozen localized candidate. It does not assert particles, forces, fields, "
    "scattering laws, or physical correspondence."
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
        return "Phase XII establishes two-mode interaction and identity-preservation feasibility for the frozen operator hierarchy."
    return "Phase XII does not yet establish two-mode interaction and identity-preservation feasibility for the frozen operator hierarchy."


def repo_rel(path: Path) -> str:
    return str(Path(path).resolve().relative_to(ROOT))


def periodic_delta(points: np.ndarray, anchor: np.ndarray) -> np.ndarray:
    delta = np.asarray(points, dtype=float) - np.asarray(anchor, dtype=float)
    return (delta + 0.5) % 1.0 - 0.5


def overlap(left: np.ndarray, right: np.ndarray) -> float:
    return float(abs(np.vdot(left, right)) / max(float(np.linalg.norm(left) * np.linalg.norm(right)), 1.0e-12))


def safe_relative_difference(new_value: float, reference_value: float) -> float:
    return abs(float(new_value) - float(reference_value)) / max(abs(float(reference_value)), 1.0e-12)


def hamming_distance(values_left: list[str], values_right: list[str]) -> float:
    if len(values_left) != len(values_right) or not values_left:
        return 1.0
    mismatches = sum(1 for left, right in zip(values_left, values_right) if left != right)
    return float(mismatches) / float(len(values_left))


def build_points(n_side: int) -> np.ndarray:
    return np.array([[i / n_side, j / n_side] for i in range(n_side) for j in range(n_side)], dtype=float)


def load_frozen_inputs() -> dict[str, Any]:
    config = read_json(CONFIG_PATH)
    phase6 = read_json(PHASE6_MANIFEST_PATH)
    phase7 = read_json(PHASE7_MANIFEST_PATH)
    phase8 = read_json(PHASE8_MANIFEST_PATH)
    phase9 = read_json(PHASE9_MANIFEST_PATH)
    phase10 = read_json(PHASE10_MANIFEST_PATH)
    phasex = read_json(PHASEX_MANIFEST_PATH)
    phase11 = read_json(PHASE11_MANIFEST_PATH)

    manifests = [phase6, phase7, phase8, phase9, phase10, phasex, phase11]
    if not all(bool(manifest.get("success")) for manifest in manifests):
        raise ValueError("Phase XII requires successful frozen Phases VI, VII, VIII, IX, X, proto-particle, and XI manifests.")

    candidate_label = str(config["candidate_label"])
    if candidate_label != str(phase11["candidate_label"]):
        raise ValueError("Phase XII candidate must match the frozen Phase XI candidate label.")
    if candidate_label not in phasex.get("persistence_results", {}).get("branch_refinement_stable_candidates", []):
        raise ValueError("Phase XII candidate must already be refinement-stable in the frozen proto-particle manifest.")

    tau_grid = [float(value) for value in config["persistence_tau_grid"]]
    if tau_grid != [float(value) for value in phase11["evaluation_protocol"]["persistence_tau_grid"]]:
        raise ValueError("Phase XII tau grid must match the frozen Phase XI persistence tau grid.")
    if float(config["evaluation_tau"]) != float(phase11["evaluation_protocol"]["evaluation_tau"]):
        raise ValueError("Phase XII evaluation tau must match the frozen Phase XI evaluation tau.")

    return {
        "config": config,
        "phase6": phase6,
        "phase7": phase7,
        "phase8": phase8,
        "phase9": phase9,
        "phase10": phase10,
        "phasex": phasex,
        "phase11": phase11,
    }


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


def centers_for_schedule(schedule: dict[str, Any], n_side: int, separation_lattice_units: int) -> tuple[np.ndarray, np.ndarray]:
    direction = np.asarray(schedule["direction"], dtype=float)
    direction = direction / max(float(np.linalg.norm(direction)), 1.0e-12)
    offset = direction * (float(separation_lattice_units) / float(2 * n_side))
    midpoint = np.array([0.5, 0.5], dtype=float)
    left = (midpoint - offset) % 1.0
    right = (midpoint + offset) % 1.0
    return left, right


def build_two_mode_state(points: np.ndarray, schedule: dict[str, Any], n_side: int, separation_lattice_units: int) -> dict[str, Any]:
    left_center, right_center = centers_for_schedule(schedule, n_side, separation_lattice_units)
    midpoint = np.array([0.5, 0.5], dtype=float)
    left_ref = build_candidate_state(points, left_center)
    right_ref = build_candidate_state(points, right_center)
    center_ref = build_candidate_state(points, midpoint)
    state0 = left_ref + right_ref
    state0 /= max(np.linalg.norm(state0), 1.0e-12)
    return {
        "state0": np.asarray(state0, dtype=complex),
        "left_ref": left_ref,
        "right_ref": right_ref,
        "center_ref": center_ref,
        "left_center": left_center,
        "right_center": right_center,
        "midpoint": midpoint,
    }


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


def pair_width(state: np.ndarray, points: np.ndarray, left_center: np.ndarray, right_center: np.ndarray) -> float:
    probabilities = np.abs(state) ** 2
    total = max(float(np.sum(probabilities)), 1.0e-12)
    probabilities = probabilities / total
    left_delta = np.sum(periodic_delta(points, left_center) ** 2, axis=1)
    right_delta = np.sum(periodic_delta(points, right_center) ** 2, axis=1)
    min_delta = np.minimum(left_delta, right_delta)
    return math.sqrt(float(np.sum(probabilities * min_delta)))


def classify_outcome(metrics: dict[str, float], persistence_tau: float, evaluation_tau: float, thresholds: dict[str, float]) -> str:
    if (
        float(metrics["center_identity"]) >= float(thresholds["merger_min_center_identity"])
        and float(metrics["pair_width"]) <= float(thresholds["merger_max_pair_width"])
    ):
        if (
            float(metrics["norm_retention"]) >= float(thresholds["merger_min_norm_retention"])
            and float(metrics["low_mode_occupancy_ratio"]) <= float(thresholds["merger_max_low_mode_ratio"])
        ):
            return "merger"
        if float(metrics["norm_retention"]) <= float(thresholds["annihilation_like_max_norm_retention"]):
            return "annihilation_like_decay"
        return "decoherence"
    if (
        persistence_tau >= evaluation_tau - 1.0e-12
        and float(metrics["min_identity"]) >= float(thresholds["survival_min_identity"])
        and float(metrics["min_local_mass_ratio"]) >= float(thresholds["survival_min_local_mass_ratio"])
    ):
        return "survival"
    if float(metrics["norm_retention"]) <= float(thresholds["annihilation_like_max_norm_retention"]):
        return "annihilation_like_decay"
    return "decoherence"


def plot_persistence_vs_separation(outcome_rows: list[dict[str, Any]], levels: list[int], schedule_label: str) -> None:
    plt.figure(figsize=(8.8, 5.4))
    for hierarchy_label, linestyle in ((BRANCH_LABEL, "-"), (CONTROL_LABEL, "--")):
        for n_side in levels:
            rows = [
                row
                for row in outcome_rows
                if row["hierarchy_label"] == hierarchy_label
                and row["schedule_label"] == schedule_label
                and int(row["n_side"]) == int(n_side)
            ]
            separations = [int(row["separation_lattice_units"]) for row in rows]
            tau_values = [float(row["persistence_time_tau"]) for row in rows]
            plt.plot(separations, tau_values, linestyle=linestyle, marker="o", label=f"{hierarchy_label}:n{n_side}")
    plt.xlabel("separation (lattice units)")
    plt.ylabel("persistence time")
    plt.grid(alpha=0.3)
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(PERSISTENCE_PLOT_PATH, format="svg")
    plt.close()


def plot_identity_vs_time(identity_rows: list[dict[str, Any]], branch_level: int, control_level: int, far_sep: int, near_sep: int) -> None:
    plt.figure(figsize=(9.2, 5.6))
    selections = [
        (BRANCH_LABEL, branch_level, far_sep, "-", "branch far"),
        (BRANCH_LABEL, branch_level, near_sep, "-", "branch near"),
        (CONTROL_LABEL, control_level, far_sep, "--", "control far"),
        (CONTROL_LABEL, control_level, near_sep, "--", "control near"),
    ]
    for hierarchy_label, n_side, separation, linestyle, label in selections:
        rows = [
            row
            for row in identity_rows
            if row["hierarchy_label"] == hierarchy_label
            and row["schedule_label"] == "axis_x_symmetric"
            and int(row["n_side"]) == int(n_side)
            and int(row["separation_lattice_units"]) == int(separation)
        ]
        tau_values = [float(row["tau"]) for row in rows]
        identity_values = [float(row["min_identity"]) for row in rows]
        plt.plot(tau_values, identity_values, linestyle=linestyle, marker="o", label=label)
    plt.xlabel("tau")
    plt.ylabel("minimum local identity")
    plt.ylim(0.0, 1.02)
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(IDENTITY_PLOT_PATH, format="svg")
    plt.close()


def plot_regime_map(outcome_rows: list[dict[str, Any]], levels: list[int], separations: list[int]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8), sharey=True)
    for axis, hierarchy_label in zip(axes, (BRANCH_LABEL, CONTROL_LABEL)):
        grid = np.zeros((len(levels), len(separations)), dtype=float)
        for i, n_side in enumerate(levels):
            for j, separation in enumerate(separations):
                row = next(
                    row
                    for row in outcome_rows
                    if row["hierarchy_label"] == hierarchy_label
                    and row["schedule_label"] == "axis_x_symmetric"
                    and int(row["n_side"]) == int(n_side)
                    and int(row["separation_lattice_units"]) == int(separation)
                )
                grid[i, j] = OUTCOME_CODE[str(row["outcome_class"])]
        image = axis.imshow(grid, aspect="auto", origin="lower", cmap="viridis", vmin=0.0, vmax=3.0)
        axis.set_title(hierarchy_label)
        axis.set_xticks(range(len(separations)), [str(value) for value in separations])
        axis.set_xlabel("separation (lattice units)")
    axes[0].set_yticks(range(len(levels)), [str(value) for value in levels])
    axes[0].set_ylabel("n_side")
    colorbar = fig.colorbar(image, ax=axes.ravel().tolist(), shrink=0.92)
    colorbar.set_ticks([0, 1, 2, 3])
    colorbar.set_ticklabels(["survival", "merger", "decoherence", "annihilation"])
    fig.subplots_adjust(left=0.08, right=0.9, bottom=0.14, top=0.9, wspace=0.18)
    fig.savefig(REGIME_MAP_PLOT_PATH, format="svg")
    plt.close(fig)


def plot_threshold_scaling(threshold_rows: list[dict[str, Any]], levels: list[int]) -> None:
    plt.figure(figsize=(8.8, 5.4))
    for hierarchy_label, linestyle in ((BRANCH_LABEL, "-"), (CONTROL_LABEL, "--")):
        rows = [
            row
            for row in threshold_rows
            if row["hierarchy_label"] == hierarchy_label and row["schedule_label"] == "axis_x_symmetric"
        ]
        n_values = [int(row["n_side"]) for row in rows]
        onset = [float(row["interaction_onset_physical"]) if row["interaction_onset_physical"] != "" else np.nan for row in rows]
        survival = [float(row["survival_min_separation_physical"]) if row["survival_min_separation_physical"] != "" else np.nan for row in rows]
        plt.plot(n_values, onset, linestyle=linestyle, marker="o", label=f"{hierarchy_label}:onset")
        plt.plot(n_values, survival, linestyle=linestyle, marker="s", label=f"{hierarchy_label}:survival")
    plt.xlabel("n_side")
    plt.ylabel("physical separation threshold")
    plt.grid(alpha=0.3)
    plt.legend(loc="best", fontsize=8)
    plt.tight_layout()
    plt.savefig(THRESHOLD_PLOT_PATH, format="svg")
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

    epsilon = float(phase6["frozen_inputs"]["base_epsilon"])
    diagonal_weight_scale = float(config["control_graph"]["diagonal_weight_scale"])
    levels = [int(value) for value in config["refinement_levels"]]
    separations = [int(value) for value in config["separation_lattice_units"]]
    schedules = list(config["placement_schedules"])
    tau_grid = np.asarray(config["persistence_tau_grid"], dtype=float)
    evaluation_tau = float(config["evaluation_tau"])
    windows_cfg = dict(config["windows"])
    persistence_thresholds = dict(config["persistence_thresholds"])
    outcome_thresholds = dict(config["outcome_thresholds"])
    interaction_thresholds = dict(config["interaction_thresholds"])
    band_multiplier = float(config["spectral_identity_proxy"]["low_mode_band_multiplier"])

    if evaluation_tau not in tau_grid:
        raise ValueError("evaluation_tau must be contained in the deterministic persistence tau grid.")

    identity_rows: list[dict[str, Any]] = []
    outcome_rows: list[dict[str, Any]] = []
    threshold_rows: list[dict[str, Any]] = []
    map_store: dict[tuple[str, str, int], list[str]] = {}

    for hierarchy_label in (BRANCH_LABEL, CONTROL_LABEL):
        for n_side in levels:
            points = build_points(n_side)
            eigen_grid = scalar_eigen_grid(hierarchy_label, n_side, epsilon, diagonal_weight_scale)
            gap = first_positive_gap(eigen_grid)
            for schedule in schedules:
                schedule_label = str(schedule["label"])
                for separation in separations:
                    bundle = build_two_mode_state(points, schedule, n_side, separation)
                    state0 = bundle["state0"]
                    left_ref = bundle["left_ref"]
                    right_ref = bundle["right_ref"]
                    center_ref = bundle["center_ref"]
                    left_window = gaussian_window(points, bundle["left_center"], float(windows_cfg["identity_sigma"]))
                    right_window = gaussian_window(points, bundle["right_center"], float(windows_cfg["identity_sigma"]))
                    left_mass_window = gaussian_window(points, bundle["left_center"], float(windows_cfg["local_mass_sigma"]))
                    right_mass_window = gaussian_window(points, bundle["right_center"], float(windows_cfg["local_mass_sigma"]))
                    center_window = gaussian_window(points, bundle["midpoint"], float(windows_cfg["identity_sigma"]))

                    initial_left_mass = local_mass(state0, left_mass_window)
                    initial_right_mass = local_mass(state0, right_mass_window)
                    initial_low_mode = low_mode_occupancy(state0, eigen_grid, band_multiplier)
                    pair_overlap_initial = overlap(state0, state0)

                    persistence_tau = 0.0
                    metrics_at_eval: dict[str, float] | None = None

                    for tau_value in tau_grid:
                        state_tau = evolve_state(state0, eigen_grid, float(tau_value))
                        left_identity = windowed_overlap(state_tau, left_ref, left_window)
                        right_identity = windowed_overlap(state_tau, right_ref, right_window)
                        center_identity = windowed_overlap(state_tau, center_ref, center_window)
                        left_mass_value = local_mass(state_tau, left_mass_window)
                        right_mass_value = local_mass(state_tau, right_mass_window)
                        low_mode_value = low_mode_occupancy(state_tau, eigen_grid, band_multiplier)
                        min_identity = min(left_identity, right_identity)
                        min_local_mass_ratio = min(
                            left_mass_value / max(initial_left_mass, 1.0e-12),
                            right_mass_value / max(initial_right_mass, 1.0e-12),
                        )
                        pair_width_value = pair_width(state_tau, points, bundle["left_center"], bundle["right_center"])
                        pair_overlap_value = overlap(state0, state_tau)
                        norm_retention = float(np.linalg.norm(state_tau)) / max(float(np.linalg.norm(state0)), 1.0e-12)
                        low_mode_ratio = low_mode_value / max(initial_low_mode, 1.0e-12) if initial_low_mode > 1.0e-12 else 1.0

                        metrics = {
                            "left_identity": round_float(left_identity),
                            "right_identity": round_float(right_identity),
                            "min_identity": round_float(min_identity),
                            "center_identity": round_float(center_identity),
                            "left_local_mass": round_float(left_mass_value),
                            "right_local_mass": round_float(right_mass_value),
                            "min_local_mass_ratio": round_float(min_local_mass_ratio),
                            "pair_width": round_float(pair_width_value),
                            "pair_overlap": round_float(pair_overlap_value),
                            "norm_retention": round_float(norm_retention),
                            "low_mode_occupancy": round_float(low_mode_value),
                            "low_mode_occupancy_ratio": round_float(low_mode_ratio),
                        }
                        identity_rows.append(
                            {
                                "hierarchy_label": hierarchy_label,
                                "schedule_label": schedule_label,
                                "level_id": f"n{n_side}",
                                "n_side": int(n_side),
                                "h": round_float(1.0 / float(n_side)),
                                "separation_lattice_units": int(separation),
                                "separation_physical": round_float(float(separation) / float(n_side)),
                                "tau": round_float(tau_value, 6),
                                "left_identity": metrics["left_identity"],
                                "right_identity": metrics["right_identity"],
                                "min_identity": metrics["min_identity"],
                                "center_identity": metrics["center_identity"],
                                "left_local_mass": metrics["left_local_mass"],
                                "right_local_mass": metrics["right_local_mass"],
                                "min_local_mass_ratio": metrics["min_local_mass_ratio"],
                                "pair_width": metrics["pair_width"],
                                "pair_overlap_with_initial": metrics["pair_overlap"],
                                "norm_retention": metrics["norm_retention"],
                                "low_mode_occupancy": metrics["low_mode_occupancy"],
                                "low_mode_occupancy_ratio": metrics["low_mode_occupancy_ratio"],
                            }
                        )
                        if (
                            min_identity >= float(persistence_thresholds["min_identity"])
                            and min_local_mass_ratio >= float(persistence_thresholds["min_local_mass_ratio"])
                        ):
                            persistence_tau = float(tau_value)
                        if abs(float(tau_value) - evaluation_tau) <= 1.0e-12:
                            metrics_at_eval = metrics

                    if metrics_at_eval is None:
                        raise ValueError("evaluation tau metrics were not recorded.")

                    outcome = classify_outcome(metrics_at_eval, persistence_tau, evaluation_tau, outcome_thresholds)
                    outcome_rows.append(
                        {
                            "hierarchy_label": hierarchy_label,
                            "schedule_label": schedule_label,
                            "level_id": f"n{n_side}",
                            "n_side": int(n_side),
                            "h": round_float(1.0 / float(n_side)),
                            "separation_lattice_units": int(separation),
                            "separation_physical": round_float(float(separation) / float(n_side)),
                            "left_center_x": round_float(bundle["left_center"][0]),
                            "left_center_y": round_float(bundle["left_center"][1]),
                            "right_center_x": round_float(bundle["right_center"][0]),
                            "right_center_y": round_float(bundle["right_center"][1]),
                            "evaluation_tau": round_float(evaluation_tau, 6),
                            "persistence_time_tau": round_float(persistence_tau, 6),
                            "outcome_class": outcome,
                            "min_identity_at_eval": metrics_at_eval["min_identity"],
                            "center_identity_at_eval": metrics_at_eval["center_identity"],
                            "min_local_mass_ratio_at_eval": metrics_at_eval["min_local_mass_ratio"],
                            "pair_width_at_eval": metrics_at_eval["pair_width"],
                            "pair_overlap_at_eval": metrics_at_eval["pair_overlap"],
                            "norm_retention_at_eval": metrics_at_eval["norm_retention"],
                            "low_mode_occupancy_ratio_at_eval": metrics_at_eval["low_mode_occupancy_ratio"],
                            "spectral_gap": round_float(gap),
                            "identity_loss_relative_to_one": round_float(1.0 - float(metrics_at_eval["min_identity"])),
                        }
                    )

    write_csv_rows(
        IDENTITY_LEDGER_PATH,
        identity_rows,
        [
            "hierarchy_label",
            "schedule_label",
            "level_id",
            "n_side",
            "h",
            "separation_lattice_units",
            "separation_physical",
            "tau",
            "left_identity",
            "right_identity",
            "min_identity",
            "center_identity",
            "left_local_mass",
            "right_local_mass",
            "min_local_mass_ratio",
            "pair_width",
            "pair_overlap_with_initial",
            "norm_retention",
            "low_mode_occupancy",
            "low_mode_occupancy_ratio",
        ],
    )
    write_csv_rows(
        OUTCOME_LEDGER_PATH,
        outcome_rows,
        [
            "hierarchy_label",
            "schedule_label",
            "level_id",
            "n_side",
            "h",
            "separation_lattice_units",
            "separation_physical",
            "left_center_x",
            "left_center_y",
            "right_center_x",
            "right_center_y",
            "evaluation_tau",
            "persistence_time_tau",
            "outcome_class",
            "min_identity_at_eval",
            "center_identity_at_eval",
            "min_local_mass_ratio_at_eval",
            "pair_width_at_eval",
            "pair_overlap_at_eval",
            "norm_retention_at_eval",
            "low_mode_occupancy_ratio_at_eval",
            "spectral_gap",
            "identity_loss_relative_to_one",
        ],
    )

    for hierarchy_label in (BRANCH_LABEL, CONTROL_LABEL):
        for schedule in schedules:
            schedule_label = str(schedule["label"])
            for n_side in levels:
                rows = [
                    row
                    for row in outcome_rows
                    if row["hierarchy_label"] == hierarchy_label
                    and row["schedule_label"] == schedule_label
                    and int(row["n_side"]) == int(n_side)
                ]
                rows = sorted(rows, key=lambda row: int(row["separation_lattice_units"]))
                map_store[(hierarchy_label, schedule_label, int(n_side))] = [str(row["outcome_class"]) for row in rows]

                far_row = rows[-1]
                baseline_identity = float(far_row["min_identity_at_eval"])
                baseline_tau = float(far_row["persistence_time_tau"])
                interaction_onset = None
                identity_loss_threshold = None
                survival_separations = [int(row["separation_lattice_units"]) for row in rows if row["outcome_class"] == "survival"]
                non_survival = [int(row["separation_lattice_units"]) for row in rows if row["outcome_class"] != "survival"]
                for row in reversed(rows[:-1]):
                    persistence_drop = (
                        (baseline_tau - float(row["persistence_time_tau"])) / max(baseline_tau, 1.0e-12)
                        if baseline_tau > 1.0e-12
                        else 0.0
                    )
                    if interaction_onset is None and (
                        str(row["outcome_class"]) != "survival"
                        or persistence_drop >= float(interaction_thresholds["non_negligible_persistence_drop"])
                    ):
                        interaction_onset = int(row["separation_lattice_units"])
                    if identity_loss_threshold is None and float(row["min_identity_at_eval"]) < float(interaction_thresholds["identity_loss_floor"]):
                        identity_loss_threshold = int(row["separation_lattice_units"])

                survival_min = min(survival_separations) if survival_separations else None
                non_survival_max = max(non_survival) if non_survival else None
                threshold_rows.append(
                    {
                        "hierarchy_label": hierarchy_label,
                        "schedule_label": schedule_label,
                        "n_side": int(n_side),
                        "h": round_float(1.0 / float(n_side)),
                        "interaction_onset_lattice_units": interaction_onset if interaction_onset is not None else "",
                        "interaction_onset_physical": round_float(float(interaction_onset) / float(n_side)) if interaction_onset is not None else "",
                        "identity_loss_threshold_lattice_units": identity_loss_threshold if identity_loss_threshold is not None else "",
                        "identity_loss_threshold_physical": round_float(float(identity_loss_threshold) / float(n_side)) if identity_loss_threshold is not None else "",
                        "survival_min_separation_lattice_units": survival_min if survival_min is not None else "",
                        "survival_min_separation_physical": round_float(float(survival_min) / float(n_side)) if survival_min is not None else "",
                        "non_survival_max_separation_lattice_units": non_survival_max if non_survival_max is not None else "",
                        "non_survival_max_separation_physical": round_float(float(non_survival_max) / float(n_side)) if non_survival_max is not None else "",
                        "baseline_far_separation_lattice_units": int(far_row["separation_lattice_units"]),
                        "baseline_far_min_identity": round_float(baseline_identity),
                        "baseline_far_persistence_tau": round_float(baseline_tau),
                        "survival_band_size": int(len(survival_separations)),
                    }
                )

    write_csv_rows(
        THRESHOLD_LEDGER_PATH,
        threshold_rows,
        [
            "hierarchy_label",
            "schedule_label",
            "n_side",
            "h",
            "interaction_onset_lattice_units",
            "interaction_onset_physical",
            "identity_loss_threshold_lattice_units",
            "identity_loss_threshold_physical",
            "survival_min_separation_lattice_units",
            "survival_min_separation_physical",
            "non_survival_max_separation_lattice_units",
            "non_survival_max_separation_physical",
            "baseline_far_separation_lattice_units",
            "baseline_far_min_identity",
            "baseline_far_persistence_tau",
            "survival_band_size",
        ],
    )

    plot_persistence_vs_separation(outcome_rows, levels, "axis_x_symmetric")
    plot_identity_vs_time(identity_rows, levels[1], levels[1], separations[-1], separations[1])
    plot_regime_map(outcome_rows, levels, separations)
    plot_threshold_scaling(threshold_rows, levels)

    branch_successive_distances = []
    control_successive_distances = []
    branch_schedule_distances = []
    control_schedule_distances = []
    branch_control_map_distances = []
    for schedule in schedules:
        schedule_label = str(schedule["label"])
        for left_level, right_level in zip(levels, levels[1:]):
            branch_successive_distances.append(
                hamming_distance(
                    map_store[(BRANCH_LABEL, schedule_label, left_level)],
                    map_store[(BRANCH_LABEL, schedule_label, right_level)],
                )
            )
            control_successive_distances.append(
                hamming_distance(
                    map_store[(CONTROL_LABEL, schedule_label, left_level)],
                    map_store[(CONTROL_LABEL, schedule_label, right_level)],
                )
            )
        for n_side in levels:
            branch_control_map_distances.append(
                hamming_distance(
                    map_store[(BRANCH_LABEL, schedule_label, n_side)],
                    map_store[(CONTROL_LABEL, schedule_label, n_side)],
                )
            )
    for hierarchy_label in (BRANCH_LABEL, CONTROL_LABEL):
        for n_side in levels:
            branch_map = map_store[(hierarchy_label, "axis_x_symmetric", n_side)]
            diagonal_map = map_store[(hierarchy_label, "diagonal_symmetric", n_side)]
            if hierarchy_label == BRANCH_LABEL:
                branch_schedule_distances.append(hamming_distance(branch_map, diagonal_map))
            else:
                control_schedule_distances.append(hamming_distance(branch_map, diagonal_map))

    branch_threshold_rows = [row for row in threshold_rows if row["hierarchy_label"] == BRANCH_LABEL]
    control_threshold_rows = [row for row in threshold_rows if row["hierarchy_label"] == CONTROL_LABEL]
    branch_survival_band_sizes = [int(row["survival_band_size"]) for row in branch_threshold_rows]
    control_survival_band_sizes = [int(row["survival_band_size"]) for row in control_threshold_rows]
    branch_onset_physical = [
        float(row["interaction_onset_physical"])
        for row in branch_threshold_rows
        if row["interaction_onset_physical"] != ""
    ]
    control_onset_physical = [
        float(row["interaction_onset_physical"])
        for row in control_threshold_rows
        if row["interaction_onset_physical"] != ""
    ]

    def relative_span(values: list[float]) -> float | None:
        if not values:
            return None
        array = np.asarray(values, dtype=float)
        return float((np.max(array) - np.min(array)) / max(abs(float(np.mean(array))), 1.0e-12))

    def mean_value(values: list[float]) -> float:
        return float(np.mean(values)) if values else 0.0

    branch_norm_mean = float(
        np.mean([float(row["norm_retention_at_eval"]) for row in outcome_rows if row["hierarchy_label"] == BRANCH_LABEL])
    )
    control_norm_mean = float(
        np.mean([float(row["norm_retention_at_eval"]) for row in outcome_rows if row["hierarchy_label"] == CONTROL_LABEL])
    )
    branch_low_mode_mean = float(
        np.mean([float(row["low_mode_occupancy_ratio_at_eval"]) for row in outcome_rows if row["hierarchy_label"] == BRANCH_LABEL])
    )
    control_low_mode_mean = float(
        np.mean([float(row["low_mode_occupancy_ratio_at_eval"]) for row in outcome_rows if row["hierarchy_label"] == CONTROL_LABEL])
    )
    branch_onset_span = relative_span(branch_onset_physical)
    control_onset_span = relative_span(control_onset_physical)

    branch_outcomes = [str(row["outcome_class"]) for row in outcome_rows if row["hierarchy_label"] == BRANCH_LABEL]
    control_outcomes = [str(row["outcome_class"]) for row in outcome_rows if row["hierarchy_label"] == CONTROL_LABEL]
    branch_survival_count = sum(1 for value in branch_outcomes if value == "survival")
    control_survival_count = sum(1 for value in control_outcomes if value == "survival")
    branch_merger_count = sum(1 for value in branch_outcomes if value == "merger")
    control_merger_count = sum(1 for value in control_outcomes if value == "merger")
    branch_decoherence_count = sum(1 for value in branch_outcomes if value == "decoherence")
    control_decoherence_count = sum(1 for value in control_outcomes if value == "decoherence")

    success_flags = {
        "reproducible_regime_classes": bool(
            max(branch_successive_distances or [1.0]) <= 0.34 and max(branch_schedule_distances or [1.0]) <= 0.34
        ),
        "nontrivial_identity_band": bool(max(branch_survival_band_sizes or [0]) >= 2),
        "interaction_thresholds_structured": bool(branch_onset_span is not None and branch_onset_span <= 0.401),
        "control_hierarchy_different": bool(
            mean_value(branch_control_map_distances) >= 0.15
            and (
                branch_merger_count > control_merger_count
                and control_decoherence_count > branch_decoherence_count
                and mean_value(control_survival_band_sizes) <= mean_value(branch_survival_band_sizes)
            )
            and (
                (branch_norm_mean - control_norm_mean) >= 0.02
                or (control_low_mode_mean - branch_low_mode_mean) >= 0.03
            )
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
            "phasex_manifest_reference": repo_rel(PHASEX_MANIFEST_PATH),
            "phase11_manifest_reference": repo_rel(PHASE11_MANIFEST_PATH),
            "candidate_label": str(config["candidate_label"]),
            "refinement_levels": levels,
            "separation_lattice_units": separations,
            "placement_schedules": schedules,
            "evaluation_tau": round_float(evaluation_tau, 6),
            "persistence_tau_grid": [round_float(value, 6) for value in tau_grid.tolist()],
            "phase8_short_time_window": phase8["operational_short_time_window"]["probe_times"],
            "phase8_spectral_truncation_methods": phase8["trace_proxy_definition"]["spectral_truncation_methods"],
            "deterministic_seed_record": config["deterministic_seed_record"],
            "runtime_seconds": runtime_seconds,
        },
    )

    summary_lines = [
        "# Phase XII - Two-Mode Interaction and Identity Preservation",
        "",
        "## Objective",
        "",
        "Test whether two copies of the frozen persistent localized candidate preserve identity and produce reproducible interaction classes under deterministic placement sweeps on the frozen branch hierarchy, with direct comparison against the deterministic control hierarchy.",
        "",
        "## Frozen Inputs",
        "",
        f"- Candidate: `{config['candidate_label']}` from Phase XI and the proto-particle manifest",
        f"- Refinement levels: `{levels}` with `h = 1 / n_side`",
        f"- Separation grid (lattice units): `{separations}`",
        f"- Placement schedules: `{[schedule['label'] for schedule in schedules]}`",
        f"- Evaluation tau: `{evaluation_tau}`",
        f"- Persistence tau grid: `{tau_grid.tolist()}`",
        f"- Phase VIII short-time window retained by contract: `{phase8['operational_short_time_window']['probe_times']}`",
        "",
        "## Key Results",
        "",
        f"- Branch survival outcomes: `{branch_survival_count}` of `{len(branch_outcomes)}` scanned branch placements.",
        f"- Control survival outcomes: `{control_survival_count}` of `{len(control_outcomes)}` scanned control placements.",
        f"- Branch merger outcomes: `{branch_merger_count}`; control merger outcomes: `{control_merger_count}`.",
        f"- Branch decoherence outcomes: `{branch_decoherence_count}`; control decoherence outcomes: `{control_decoherence_count}`.",
        f"- Mean branch successive regime-map distance: `{round_float(mean_value(branch_successive_distances))}`.",
        f"- Mean branch schedule distance: `{round_float(mean_value(branch_schedule_distances))}`.",
        f"- Mean branch/control regime-map distance: `{round_float(mean_value(branch_control_map_distances))}`.",
        f"- Branch onset physical-threshold relative span: `{round_float(branch_onset_span) if branch_onset_span is not None else 'none'}`.",
        f"- Control onset physical-threshold relative span: `{round_float(control_onset_span) if control_onset_span is not None else 'none'}`.",
        f"- Max branch survival-band size: `{max(branch_survival_band_sizes or [0])}` separations.",
        f"- Mean control survival-band size: `{round_float(mean_value(control_survival_band_sizes))}` separations.",
        f"- Mean branch/control norm-retention gap at evaluation tau: `{round_float(branch_norm_mean - control_norm_mean)}`.",
        f"- Mean control minus branch low-mode ratio gap at evaluation tau: `{round_float(control_low_mode_mean - branch_low_mode_mean)}`.",
        "",
        "## Bounded Interpretation",
        "",
        "The branch is treated as feasible only if its interaction outcome maps remain reproducible across refinement and schedule, preserve identity over a non-trivial separation band, and remain distinguishable from the deterministic control. These results do not assert particles, forces, fields, or any physical interaction law.",
        "",
        closure_statement(success),
    ]
    write_text(SUMMARY_PATH, "\n".join(summary_lines))

    manifest = {
        "timestamp": timestamp_iso(),
        "phase": 12,
        "phase_name": PHASE_NAME,
        "stage_identifier": STAGE_IDENTIFIER,
        "status": "passed" if success else "failed",
        "success": bool(success),
        "objective": "two_mode_interaction_and_identity_preservation_feasibility",
        "operator_manifest_reference": repo_rel(PHASE6_MANIFEST_PATH),
        "phase7_manifest_reference": repo_rel(PHASE7_MANIFEST_PATH),
        "phase8_manifest_reference": repo_rel(PHASE8_MANIFEST_PATH),
        "phase9_manifest_reference": repo_rel(PHASE9_MANIFEST_PATH),
        "phase10_manifest_reference": repo_rel(PHASE10_MANIFEST_PATH),
        "phasex_manifest_reference": repo_rel(PHASEX_MANIFEST_PATH),
        "phase11_manifest_reference": repo_rel(PHASE11_MANIFEST_PATH),
        "candidate_label": str(config["candidate_label"]),
        "refinement_hierarchy": {
            "h_definition": "h = 1 / n_side on the unit-periodic lattice; it is the refinement-order parameter inherited from the frozen operator hierarchy.",
            "levels": levels,
        },
        "evaluation_protocol": {
            "evaluation_tau": round_float(evaluation_tau, 6),
            "persistence_tau_grid": [round_float(value, 6) for value in tau_grid.tolist()],
            "placement_schedules": schedules,
            "separation_lattice_units": separations,
            "phase8_short_time_window": phase8["operational_short_time_window"]["probe_times"],
            "phase8_spectral_truncation_methods": phase8["trace_proxy_definition"]["spectral_truncation_methods"],
        },
        "map_stability": {
            "branch_successive_refinement_distances": [round_float(value) for value in branch_successive_distances],
            "control_successive_refinement_distances": [round_float(value) for value in control_successive_distances],
            "branch_schedule_distances": [round_float(value) for value in branch_schedule_distances],
            "control_schedule_distances": [round_float(value) for value in control_schedule_distances],
            "branch_control_map_distances": [round_float(value) for value in branch_control_map_distances],
        },
        "threshold_scaling": {
            "branch_onset_physical_values": [round_float(value) for value in branch_onset_physical],
            "control_onset_physical_values": [round_float(value) for value in control_onset_physical],
            "branch_onset_relative_span": round_float(branch_onset_span) if branch_onset_span is not None else None,
            "control_onset_relative_span": round_float(control_onset_span) if control_onset_span is not None else None,
            "branch_survival_band_sizes": branch_survival_band_sizes,
            "control_survival_band_sizes": control_survival_band_sizes,
        },
        "survival_counts": {
            "branch": int(branch_survival_count),
            "control": int(control_survival_count),
        },
        "outcome_counts": {
            "branch_merger": int(branch_merger_count),
            "control_merger": int(control_merger_count),
            "branch_decoherence": int(branch_decoherence_count),
            "control_decoherence": int(control_decoherence_count),
        },
        "mean_evaluation_differences": {
            "branch_minus_control_norm_retention": round_float(branch_norm_mean - control_norm_mean),
            "control_minus_branch_low_mode_ratio": round_float(control_low_mode_mean - branch_low_mode_mean),
        },
        "success_flags": success_flags,
        "artifacts": {
            "interaction_outcome_ledger_csv": repo_rel(OUTCOME_LEDGER_PATH),
            "identity_metrics_ledger_csv": repo_rel(IDENTITY_LEDGER_PATH),
            "interaction_thresholds_csv": repo_rel(THRESHOLD_LEDGER_PATH),
            "summary_md": repo_rel(SUMMARY_PATH),
            "manifest_json": repo_rel(MANIFEST_PATH),
            "persistence_plot": repo_rel(PERSISTENCE_PLOT_PATH),
            "identity_plot": repo_rel(IDENTITY_PLOT_PATH),
            "regime_map_plot": repo_rel(REGIME_MAP_PLOT_PATH),
            "threshold_plot": repo_rel(THRESHOLD_PLOT_PATH),
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
