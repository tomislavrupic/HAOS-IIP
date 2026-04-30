#!/usr/bin/env python3
"""Phase 62 multi-pole / supertranslation toy probe.

This sidecar extends the Phase 61 toy memory benchmark from one permanent
memory-like deformation to a structured sequence of low-l mode shifts on S2.
The target is a synthetic "supertranslation-like" memory function, represented
as ordered l=2,3,4 mode packets. The probe asks whether HAOS-style spectral
telemetry can recover the cumulative field, mode address, event ordering, and
permanence better than controls that preserve smoothness or partial statistics.

It does not test real BMS supertranslations, gravitational memory, celestial
holography, Ward identities, or S-matrix data.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
ROOT = Path(__file__).resolve().parent

CSV_PATH = ROOT / "multipole_diagnostics.csv"
STATUS_PATH = ROOT / "bridge_status.json"
REPORT_PATH = ROOT / "multipole_supertranslation_report.md"

FIELDNAMES = [
    "sample",
    "sample_role",
    "control_type",
    "supertranslation_specificity_score",
    "final_recoverability",
    "multipole_address_score",
    "event_order_score",
    "permanence_score",
    "leakage_control_score",
    "generic_smoothness_score",
    "recovered_target_correlation",
    "recovered_amplitude",
    "address_correlation",
    "band_power_score",
    "event_time_score",
    "plateau_stability",
    "off_band_leakage",
]


@dataclass(frozen=True)
class ProbeConfig:
    resolution: int = 192
    time_count: int = 120
    l_max: int = 4
    k_neighbors: int = 10
    seed: int = 6201
    event_times: tuple[float, float, float] = (0.26, 0.48, 0.70)
    event_width: float = 0.028
    noise_scale: float = 0.010
    pass_min_specificity_score: float = 0.86
    pass_min_recoverability: float = 0.90
    pass_min_address: float = 0.88
    pass_min_event_order: float = 0.82
    pass_min_permanence: float = 0.82
    pass_min_margin: float = 0.18
    pass_min_separation_z: float = 1.50
    pass_max_control_false_positive_rate: float = 0.0


@dataclass(frozen=True)
class TargetSpec:
    points: np.ndarray
    component_fields: list[np.ndarray]
    component_units: list[np.ndarray]
    total_field: np.ndarray
    memory_basis: np.ndarray
    band_bases: dict[int, np.ndarray]
    target_coefficients: np.ndarray
    target_band_powers: dict[int, float]


@dataclass(frozen=True)
class SampleSpec:
    sample: str
    sample_role: str
    control_type: str
    series: np.ndarray


@dataclass(frozen=True)
class SampleSummary:
    sample: str
    sample_role: str
    control_type: str
    supertranslation_specificity_score: float
    final_recoverability: float
    multipole_address_score: float
    event_order_score: float
    permanence_score: float
    leakage_control_score: float
    generic_smoothness_score: float
    recovered_target_correlation: float
    recovered_amplitude: float
    address_correlation: float
    band_power_score: float
    event_time_score: float
    plateau_stability: float
    off_band_leakage: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 62 multi-pole supertranslation toy probe.")
    parser.add_argument("--resolution", type=int, default=192)
    parser.add_argument("--time-count", type=int, default=120)
    parser.add_argument("--seed", type=int, default=6201)
    parser.add_argument("--noise-scale", type=float, default=0.010)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT,
        help="Directory for generated CSV, JSON, markdown report, and figures.",
    )
    return parser.parse_args()


def make_config(args: argparse.Namespace) -> ProbeConfig:
    return ProbeConfig(
        resolution=int(args.resolution),
        time_count=int(args.time_count),
        seed=int(args.seed),
        noise_scale=float(args.noise_scale),
    )


def fibonacci_sphere(n: int) -> np.ndarray:
    indices = np.arange(n, dtype=float)
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    z = 1.0 - 2.0 * (indices + 0.5) / float(n)
    radius = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    phi = golden_angle * indices
    return np.column_stack((radius * np.cos(phi), radius * np.sin(phi), z)).astype(float)


def pairwise_distances(points: np.ndarray) -> np.ndarray:
    diff = points[:, None, :] - points[None, :, :]
    return np.sqrt(np.sum(diff * diff, axis=2))


def knn_graph(points: np.ndarray, k_neighbors: int) -> np.ndarray:
    n = points.shape[0]
    distances = pairwise_distances(points)
    kth = np.partition(distances, min(k_neighbors, n - 1), axis=1)[:, min(k_neighbors, n - 1)]
    scale = max(float(np.median(kth)), 1.0e-12)
    adjacency = np.zeros((n, n), dtype=float)
    for i in range(n):
        neighbors = np.argsort(distances[i])[1 : k_neighbors + 1]
        adjacency[i, neighbors] = np.exp(-((distances[i, neighbors] / scale) ** 2))
    adjacency = np.maximum(adjacency, adjacency.T)
    np.fill_diagonal(adjacency, 0.0)
    return adjacency


def raw_harmonic_columns(points: np.ndarray, l_value: int) -> np.ndarray:
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    if l_value == 0:
        return np.column_stack([np.ones_like(x)])
    if l_value == 1:
        return np.column_stack([x, y, z])
    if l_value == 2:
        return np.column_stack([x * y, y * z, z * x, x * x - y * y, 3.0 * z * z - 1.0])
    if l_value == 3:
        return np.column_stack(
            [
                z * (5.0 * z * z - 3.0),
                x * (5.0 * z * z - 1.0),
                y * (5.0 * z * z - 1.0),
                z * (x * x - y * y),
                2.0 * x * y * z,
                x * (x * x - 3.0 * y * y),
                y * (3.0 * x * x - y * y),
            ]
        )
    if l_value == 4:
        return np.column_stack(
            [
                35.0 * z**4 - 30.0 * z * z + 3.0,
                x * z * (7.0 * z * z - 3.0),
                y * z * (7.0 * z * z - 3.0),
                (x * x - y * y) * (7.0 * z * z - 1.0),
                2.0 * x * y * (7.0 * z * z - 1.0),
                x * z * (x * x - 3.0 * y * y),
                y * z * (3.0 * x * x - y * y),
                x**4 - 6.0 * x * x * y * y + y**4,
                4.0 * x * y * (x * x - y * y),
            ]
        )
    raise ValueError(f"closed-form real harmonic basis implemented only for l <= 4, got {l_value}")


def orthonormalize(columns: np.ndarray, against: np.ndarray | None = None) -> np.ndarray:
    residual = columns.astype(float).copy()
    if against is not None and against.size:
        residual = residual - against @ (against.T @ residual)
    q, r = np.linalg.qr(residual)
    keep = np.abs(np.diag(r)) > 1.0e-10
    return q[:, keep]


def harmonic_subspaces(points: np.ndarray, l_max: int) -> dict[int, np.ndarray]:
    subspaces: dict[int, np.ndarray] = {}
    previous = np.zeros((points.shape[0], 0), dtype=float)
    for l_value in range(l_max + 1):
        basis = orthonormalize(raw_harmonic_columns(points, l_value), previous)
        subspaces[l_value] = basis
        previous = np.column_stack([previous, basis])
    return subspaces


def normalize_field(field: np.ndarray) -> np.ndarray:
    centered = field.astype(float) - float(np.mean(field))
    rms = float(np.sqrt(np.mean(centered * centered)))
    if rms <= 1.0e-12:
        return np.zeros_like(centered)
    return centered / rms


def cosine_similarity(left: np.ndarray, right: np.ndarray) -> float:
    denom = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denom <= 1.0e-12:
        return 0.0
    return float(np.dot(left, right) / denom)


def pearson_correlation(left: np.ndarray, right: np.ndarray) -> float:
    return cosine_similarity(left - float(np.mean(left)), right - float(np.mean(right)))


def sigmoid_gate(time: np.ndarray, event_time: float, event_width: float) -> np.ndarray:
    return 0.5 * (1.0 + np.tanh((time - event_time) / max(event_width, 1.0e-9)))


def pre_post_indices(config: ProbeConfig) -> tuple[np.ndarray, np.ndarray]:
    time = np.linspace(0.0, 1.0, config.time_count)
    pre = np.flatnonzero(time <= max(0.10, config.event_times[0] - 0.14))
    post = np.flatnonzero(time >= min(0.88, config.event_times[-1] + 0.18))
    return pre, post


def recovered_memory(series: np.ndarray, config: ProbeConfig) -> np.ndarray:
    pre, post = pre_post_indices(config)
    return np.mean(series[post], axis=0) - np.mean(series[pre], axis=0)


def projection_trace(series: np.ndarray, field: np.ndarray) -> np.ndarray:
    return series @ field / max(float(np.dot(field, field)), 1.0e-12)


def target_spec(points: np.ndarray, l_max: int) -> TargetSpec:
    subspaces = harmonic_subspaces(points, l_max)
    band_coeffs = {
        2: np.array([0.76, -0.31, 0.46, 0.27, -0.38], dtype=float),
        3: np.array([0.35, -0.44, 0.22, 0.31, -0.27, 0.18, 0.29], dtype=float),
        4: np.array([0.24, -0.18, 0.35, -0.22, 0.16, 0.29, -0.14, 0.21, -0.19], dtype=float),
    }
    weights = {2: 1.00, 3: 0.74, 4: 0.52}
    component_fields: list[np.ndarray] = []
    component_units: list[np.ndarray] = []
    for band in (2, 3, 4):
        unit = normalize_field(subspaces[band] @ band_coeffs[band])
        component_units.append(unit)
        component_fields.append(weights[band] * unit)
    total = np.sum(np.vstack(component_fields), axis=0)
    memory_basis = np.column_stack([subspaces[2], subspaces[3], subspaces[4]])
    total_coefficients = memory_basis.T @ total
    target_band_powers = band_power_fractions(total, {band: subspaces[band] for band in (2, 3, 4)})
    return TargetSpec(
        points=points,
        component_fields=component_fields,
        component_units=component_units,
        total_field=total,
        memory_basis=memory_basis,
        band_bases={band: subspaces[band] for band in (2, 3, 4)},
        target_coefficients=total_coefficients,
        target_band_powers=target_band_powers,
    )


def random_low_mode_pattern(
    rng: np.random.Generator,
    basis: np.ndarray,
    target: np.ndarray | None = None,
) -> np.ndarray:
    coeffs = rng.normal(size=basis.shape[1])
    field = basis @ coeffs
    if target is not None:
        projection = float(np.dot(field, target) / max(float(np.dot(target, target)), 1.0e-12))
        field = field - projection * target
    return normalize_field(field)


def build_series_from_components(
    config: ProbeConfig,
    component_fields: list[np.ndarray],
    rng: np.random.Generator,
    noise_scale: float | None = None,
    event_times: tuple[float, float, float] | None = None,
) -> np.ndarray:
    time = np.linspace(0.0, 1.0, config.time_count)
    active_times = event_times if event_times is not None else config.event_times
    series = np.zeros((config.time_count, component_fields[0].shape[0]), dtype=float)
    for event_time, field in zip(active_times, component_fields):
        series += sigmoid_gate(time, event_time, config.event_width)[:, None] * field[None, :]
    scale = config.noise_scale if noise_scale is None else noise_scale
    series += rng.normal(scale=scale, size=series.shape)
    return series


def build_samples(config: ProbeConfig) -> tuple[np.ndarray, TargetSpec, list[SampleSpec]]:
    rng = np.random.default_rng(config.seed)
    points = fibonacci_sphere(config.resolution)
    target = target_spec(points, config.l_max)
    time = np.linspace(0.0, 1.0, config.time_count)
    target_series = build_series_from_components(config, target.component_fields, rng)

    single_basis = target.band_bases[2]
    single_field = random_low_mode_pattern(rng, single_basis, target=target.total_field)
    single_field *= float(np.linalg.norm(target.total_field) / max(np.linalg.norm(single_field), 1.0e-12))
    single_series = build_series_from_components(config, [single_field, 0.0 * single_field, 0.0 * single_field], rng)

    shuffled_order = [target.component_fields[1], target.component_fields[2], target.component_fields[0]]
    event_order_shuffle = build_series_from_components(config, shuffled_order, rng)

    scrambled_components: list[np.ndarray] = []
    for band, reference in zip((2, 3, 4), target.component_fields):
        scrambled = random_low_mode_pattern(rng, target.band_bases[band], target=reference)
        scrambled *= float(np.linalg.norm(reference) / max(np.linalg.norm(scrambled), 1.0e-12))
        scrambled_components.append(scrambled)
    coefficient_scramble = build_series_from_components(config, scrambled_components, rng)

    permutation = rng.permutation(config.resolution)
    coordinate_permutation_components = [field[permutation] for field in target.component_fields]
    coordinate_permutation = build_series_from_components(config, coordinate_permutation_components, rng)

    transient = np.zeros_like(target_series)
    for event_time, field in zip(config.event_times, target.component_fields):
        burst = np.exp(-((time - event_time) / 0.050) ** 2)
        transient += burst[:, None] * field[None, :]
    transient += rng.normal(scale=config.noise_scale, size=transient.shape)

    decoy_basis = np.column_stack([harmonic_subspaces(points, config.l_max)[1], target.band_bases[2]])
    decoy_components = []
    for reference in target.component_fields:
        decoy = random_low_mode_pattern(rng, decoy_basis, target=target.total_field)
        decoy *= float(np.linalg.norm(reference) / max(np.linalg.norm(decoy), 1.0e-12))
        decoy_components.append(decoy)
    smooth_decoy = build_series_from_components(config, decoy_components, rng)

    drift_coeffs = rng.normal(scale=0.017, size=(config.time_count, target.memory_basis.shape[1]))
    drift_coeffs = np.cumsum(drift_coeffs, axis=0)
    drift = drift_coeffs @ target.memory_basis.T
    drift = drift / max(float(np.std(drift)), 1.0e-12) * 0.58
    drift += rng.normal(scale=config.noise_scale, size=drift.shape)

    return time, target, [
        SampleSpec("multipole_supertranslation_target", "target", "none", target_series),
        SampleSpec("single_pole_memory_control", "negative_control", "single_l2_step", single_series),
        SampleSpec("event_order_shuffle_control", "negative_control", "right_modes_wrong_order", event_order_shuffle),
        SampleSpec("coefficient_scramble_control", "negative_control", "right_bands_wrong_coefficients", coefficient_scramble),
        SampleSpec("coordinate_permutation_control", "negative_control", "node_identity_permutation", coordinate_permutation),
        SampleSpec("transient_multipole_control", "negative_control", "multipole_bursts_without_memory", transient),
        SampleSpec("smooth_supertranslation_decoy", "negative_control", "smooth_wrong_band_step_sequence", smooth_decoy),
        SampleSpec("stochastic_multipole_drift_control", "negative_control", "smooth_random_walk_no_event_lock", drift),
    ]


def graph_smoothness_score(adjacency: np.ndarray, field: np.ndarray) -> float:
    upper = np.argwhere(np.triu(adjacency, k=1) > 1.0e-12)
    if len(upper) == 0:
        return 0.0
    weights = np.array([adjacency[i, j] for i, j in upper], dtype=float)
    diffs = np.array([field[i] - field[j] for i, j in upper], dtype=float)
    energy = float(np.sum(weights * diffs * diffs) / max(float(np.sum(weights)), 1.0e-12))
    variance = float(np.var(field))
    normalized = energy / max(variance, 1.0e-12)
    return float(np.clip(1.0 - normalized / 0.66, 0.0, 1.0))


def final_recoverability_score(recovered: np.ndarray, truth: np.ndarray) -> tuple[float, float, float]:
    corr = cosine_similarity(recovered, truth)
    corr_score = float(np.clip(max(corr, 0.0), 0.0, 1.0))
    amplitude = float(np.linalg.norm(recovered) / max(float(np.linalg.norm(truth)), 1.0e-12))
    amplitude_score = float(np.clip(1.0 - abs(amplitude - 1.0) / 0.35, 0.0, 1.0))
    return float(0.74 * corr_score + 0.26 * amplitude_score), corr, amplitude


def band_power_fractions(field: np.ndarray, band_bases: dict[int, np.ndarray]) -> dict[int, float]:
    energies: dict[int, float] = {}
    for band, basis in band_bases.items():
        projected = basis @ (basis.T @ field)
        energies[band] = float(np.dot(projected, projected))
    total = max(sum(energies.values()), 1.0e-12)
    return {band: value / total for band, value in energies.items()}


def multipole_address_score(recovered: np.ndarray, target: TargetSpec) -> tuple[float, float, float, float, float]:
    coefficients = target.memory_basis.T @ recovered
    projected = target.memory_basis @ coefficients
    recovered_energy = max(float(np.dot(recovered, recovered)), 1.0e-12)
    in_band_energy = float(np.dot(projected, projected))
    off_band_leakage = float(np.clip(1.0 - in_band_energy / recovered_energy, 0.0, 1.0))
    leakage_score = float(np.clip(1.0 - off_band_leakage / 0.16, 0.0, 1.0))
    address_corr = cosine_similarity(coefficients, target.target_coefficients)
    address_score = float(np.clip(max(address_corr, 0.0), 0.0, 1.0))
    recovered_powers = band_power_fractions(recovered, target.band_bases)
    l1_error = sum(abs(recovered_powers[band] - target.target_band_powers[band]) for band in target.band_bases)
    band_power_score = float(np.clip(1.0 - l1_error / 0.58, 0.0, 1.0))
    score = float(0.68 * address_score + 0.22 * band_power_score + 0.10 * leakage_score)
    return score, address_corr, band_power_score, leakage_score, off_band_leakage


def permanence_score(config: ProbeConfig, series: np.ndarray, truth: np.ndarray, recovered: np.ndarray) -> tuple[float, float]:
    pre, post = pre_post_indices(config)
    trace = projection_trace(series, truth)
    shift = abs(float(np.mean(trace[post]) - np.mean(trace[pre])))
    post_stability = float(np.std(trace[post]) / max(shift, 1.0e-12))
    pre_stability = float(np.std(trace[pre]) / max(shift, 1.0e-12))
    stability = float(np.clip(1.0 - (post_stability + pre_stability) / 0.13, 0.0, 1.0))
    residual_norm = float(np.linalg.norm(recovered) / max(float(np.linalg.norm(truth)), 1.0e-12))
    residual_score = float(np.clip(residual_norm / 0.82, 0.0, 1.0))
    contrast_score = float(np.clip(shift / 0.78, 0.0, 1.0))
    return float(0.48 * stability + 0.30 * contrast_score + 0.22 * residual_score), stability


def event_order_score(config: ProbeConfig, series: np.ndarray, target: TargetSpec) -> tuple[float, float]:
    time = np.linspace(0.0, 1.0, config.time_count)
    component_scores: list[float] = []
    time_scores: list[float] = []
    for event_time, unit in zip(config.event_times, target.component_units):
        gate = sigmoid_gate(time, event_time, config.event_width)
        trace = projection_trace(series, unit)
        corr_score = float(np.clip(max(pearson_correlation(trace, gate), 0.0), 0.0, 1.0))
        before = trace[time <= max(0.05, event_time - 0.10)]
        after = trace[time >= min(0.95, event_time + 0.12)]
        contrast = float(np.mean(after) - np.mean(before))
        contrast_score = float(np.clip(contrast / 0.56, 0.0, 1.0))
        derivative = np.gradient(trace, time)
        detected_time = float(time[int(np.argmax(derivative))])
        time_score = float(np.clip(1.0 - abs(detected_time - event_time) / 0.10, 0.0, 1.0))
        # Phase 62 is specifically about ordered multi-pole binding. A control
        # with the right final field but wrong event order should not pass just
        # because cumulative step traces remain smooth and permanent.
        component_scores.append(0.25 * corr_score + 0.20 * contrast_score + 0.55 * time_score)
        time_scores.append(time_score)
    return float(np.mean(component_scores)), float(np.mean(time_scores))


def evaluate_sample(
    config: ProbeConfig,
    sample: SampleSpec,
    target: TargetSpec,
    adjacency: np.ndarray,
) -> SampleSummary:
    recovered = recovered_memory(sample.series, config)
    recoverability, corr, amplitude = final_recoverability_score(recovered, target.total_field)
    address, address_corr, band_power, leakage_score, off_band_leakage = multipole_address_score(
        recovered,
        target,
    )
    event_order, event_time = event_order_score(config, sample.series, target)
    permanence, plateau_stability = permanence_score(config, sample.series, target.total_field, recovered)
    smoothness = graph_smoothness_score(adjacency, recovered)
    specificity = float(
        0.24 * recoverability
        + 0.24 * address
        + 0.34 * event_order
        + 0.13 * permanence
        + 0.03 * leakage_score
        + 0.02 * smoothness
    )
    return SampleSummary(
        sample=sample.sample,
        sample_role=sample.sample_role,
        control_type=sample.control_type,
        supertranslation_specificity_score=specificity,
        final_recoverability=recoverability,
        multipole_address_score=address,
        event_order_score=event_order,
        permanence_score=permanence,
        leakage_control_score=leakage_score,
        generic_smoothness_score=smoothness,
        recovered_target_correlation=corr,
        recovered_amplitude=amplitude,
        address_correlation=address_corr,
        band_power_score=band_power,
        event_time_score=event_time,
        plateau_stability=plateau_stability,
        off_band_leakage=off_band_leakage,
    )


def run_probe(config: ProbeConfig) -> tuple[np.ndarray, TargetSpec, list[SampleSpec], list[SampleSummary], dict[str, Any]]:
    time, target, samples = build_samples(config)
    adjacency = knn_graph(target.points, config.k_neighbors)
    summaries = [evaluate_sample(config, sample, target, adjacency) for sample in samples]
    status = build_status(config, summaries)
    return time, target, samples, summaries, status


def build_status(config: ProbeConfig, summaries: list[SampleSummary]) -> dict[str, Any]:
    target = next(row for row in summaries if row.sample_role == "target")
    controls = [row for row in summaries if row.sample_role == "negative_control"]
    control_scores = np.array([row.supertranslation_specificity_score for row in controls], dtype=float)
    control_threshold = target.supertranslation_specificity_score - config.pass_min_margin
    control_false_positive_rate = float(np.mean(control_scores >= control_threshold))
    control_max = max(controls, key=lambda row: row.supertranslation_specificity_score)
    score_margin = float(target.supertranslation_specificity_score - control_max.supertranslation_specificity_score)
    separation_z = float(
        (target.supertranslation_specificity_score - float(np.mean(control_scores)))
        / max(float(np.std(control_scores)), 1.0e-9)
    )
    target_passes = bool(
        target.supertranslation_specificity_score >= config.pass_min_specificity_score
        and target.final_recoverability >= config.pass_min_recoverability
        and target.multipole_address_score >= config.pass_min_address
        and target.event_order_score >= config.pass_min_event_order
        and target.permanence_score >= config.pass_min_permanence
    )
    controls_pass = bool(
        control_false_positive_rate <= config.pass_max_control_false_positive_rate
        and score_margin >= config.pass_min_margin
        and separation_z >= config.pass_min_separation_z
    )
    if target_passes and controls_pass:
        bridge_status = "PASS"
        failure_mode = "NONE"
    elif target.supertranslation_specificity_score >= 0.72 and separation_z >= 1.0:
        bridge_status = "MARGINAL"
        failure_mode = "MULTIPOLE_CONTROL_SEPARATION_INCOMPLETE"
    else:
        bridge_status = "OPEN"
        failure_mode = "SUPERTRANSLATION_TOY_SIGNATURE_NOT_DISTINCT"
    return {
        "phase": "62",
        "bridge_name": "multipole_supertranslation_probe",
        "bridge_status": bridge_status,
        "failure_mode": failure_mode,
        "success_definition": (
            "PASS requires a synthetic ordered l=2,3,4 permanent memory target to recover "
            "the cumulative field, multipole address, event order, and permanence while "
            "smoothness-preserving and partial-statistic controls fail."
        ),
        "config": asdict(config),
        "target_supertranslation_specificity_score": float(target.supertranslation_specificity_score),
        "target_final_recoverability": float(target.final_recoverability),
        "target_multipole_address_score": float(target.multipole_address_score),
        "target_event_order_score": float(target.event_order_score),
        "target_permanence_score": float(target.permanence_score),
        "target_leakage_control_score": float(target.leakage_control_score),
        "target_generic_smoothness_score": float(target.generic_smoothness_score),
        "control_max_score": float(control_max.supertranslation_specificity_score),
        "control_max_sample": control_max.sample,
        "control_false_positive_threshold": float(control_threshold),
        "control_false_positive_rate": control_false_positive_rate,
        "score_margin_over_best_control": score_margin,
        "separation_z": separation_z,
        "core_modified": False,
        "safe_to_continue": True,
        "next_phase": "63_claim_gated_celestial_toy_summary_or_bms_audit",
        "non_claims": [
            "no real BMS supertranslation recovery claim",
            "no real gravitational-memory observable claim",
            "no BMS charge or soft-hair recovery claim",
            "no celestial holography recovery claim",
            "no soft theorem or S-matrix reconstruction claim",
            "synthetic ordered multi-pole S2 toy data only",
        ],
    }


def write_csv(path: Path, rows: list[SampleSummary]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)


def write_status(path: Path, status: dict[str, Any]) -> None:
    path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")


def write_report(path: Path, summaries: list[SampleSummary], status: dict[str, Any]) -> None:
    summary_rows = "\n".join(
        "| {sample} | {sample_role} | {control_type} | {supertranslation_specificity_score:.6f} | {final_recoverability:.6f} | {multipole_address_score:.6f} | {event_order_score:.6f} | {permanence_score:.6f} | {generic_smoothness_score:.6f} |".format(
            **asdict(row)
        )
        for row in summaries
    )
    controls = [row for row in summaries if row.sample_role == "negative_control"]
    leakage_rows = "\n".join(
        f"- `{row.sample}`: supertranslation={row.supertranslation_specificity_score:.6f}, "
        f"final={row.final_recoverability:.6f}, address={row.multipole_address_score:.6f}, "
        f"event={row.event_order_score:.6f}, permanence={row.permanence_score:.6f}"
        for row in sorted(controls, key=lambda item: item.supertranslation_specificity_score, reverse=True)
    )
    report = f"""# Phase 62 Multi-Pole Supertranslation Toy Probe

## Scope

This file is generated by `experiments/physics_bridge/multipole_supertranslation_probe/run_multipole_supertranslation_probe.py`.

This is a synthetic ordered multi-pole memory benchmark on a sampled S2 boundary. The target is an l=2,3,4 sequence of permanent mode shifts. It does not use real BMS charges, gravitational-wave memory data, celestial amplitudes, Ward identities, or S-matrix inputs.

The narrow question is:

> Can HAOS-style spectral telemetry distinguish a structured supertranslation-like toy deformation from smooth decoys, wrong multipole addresses, wrong event ordering, transient bursts, and stochastic drift?

## Bridge Status

- `bridge_status`: `{status['bridge_status']}`
- `failure_mode`: `{status['failure_mode']}`
- target supertranslation score: `{status['target_supertranslation_specificity_score']:.6f}`
- target final recoverability: `{status['target_final_recoverability']:.6f}`
- target multipole address: `{status['target_multipole_address_score']:.6f}`
- target event order: `{status['target_event_order_score']:.6f}`
- target permanence: `{status['target_permanence_score']:.6f}`
- target leakage control: `{status['target_leakage_control_score']:.6f}`
- target generic smoothness: `{status['target_generic_smoothness_score']:.6f}`
- control max sample: `{status['control_max_sample']}`
- control max score: `{status['control_max_score']:.6f}`
- control false-positive threshold: `{status['control_false_positive_threshold']:.6f}`
- control false-positive rate: `{status['control_false_positive_rate']:.6f}`
- separation z: `{status['separation_z']:.6f}`
- margin over best control: `{status['score_margin_over_best_control']:.6f}`

## Summary Table

| sample | role | control | supertranslation_score | final | address | event_order | permanence | smoothness |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
{summary_rows}

## Control Leakage Read

{leakage_rows}

## 62 Hardening Note

The strongest control is expected to be `event_order_shuffle_control`: it keeps
the same cumulative final field and nearly the same multipole address while
binding the mode packets to the wrong event order. Phase 62 therefore weights
ordered address-binding more strongly than Phase 61 weighted a single memory
step. This is intentional. A multi-pole toy supertranslation benchmark should
not pass on final-field recovery alone.

## Interpretation

`PASS` means the toy target is separable from controls on cumulative permanent field recovery, multi-pole harmonic address, event order, and permanence.

`MARGINAL` means the toy target is detected, but at least one control remains too competitive.

`OPEN` means the telemetry is mostly measuring smoothness, low-mode content, or generic step structure rather than the ordered multi-pole target.

## Claim Gate

This sidecar remains inside the Phase 60 claim gate. Even under `PASS`, the established physics rows for BMS supertranslations, BMS charge recovery, real gravitational memory, celestial holography, and real soft theorems remain `OPEN`.

Allowed language:

- multi-pole supertranslation toy proxy
- synthetic ordered S2 memory deformation
- toy multi-mode memory-address recovery

Disallowed language:

- BMS supertranslations recovered
- gravitational memory recovered
- soft hair detected
- celestial holography validated

## Authority

- CSV: `experiments/physics_bridge/multipole_supertranslation_probe/multipole_diagnostics.csv`
- JSON: `experiments/physics_bridge/multipole_supertranslation_probe/bridge_status.json`
- Figures: `experiments/physics_bridge/multipole_supertranslation_probe/figures/`
- Generator: `experiments/physics_bridge/multipole_supertranslation_probe/run_multipole_supertranslation_probe.py`
"""
    path.write_text(report, encoding="utf-8")


def write_plots(
    output_dir: Path,
    config: ProbeConfig,
    time: np.ndarray,
    target: TargetSpec,
    samples: list[SampleSpec],
    summaries: list[SampleSummary],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    labels = [row.sample for row in summaries]
    x = np.arange(len(summaries))

    fig, ax = plt.subplots(figsize=(12, 5))
    colors = ["tab:green" if row.sample_role == "target" else "tab:orange" for row in summaries]
    ax.bar(labels, [row.supertranslation_specificity_score for row in summaries], color=colors)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("supertranslation toy specificity")
    ax.set_title("Phase 62 multi-pole target vs controls")
    ax.tick_params(axis="x", rotation=25, labelsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "control_comparison.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, 5))
    width = 0.15
    ax.bar(x - 2.0 * width, [row.final_recoverability for row in summaries], width, label="final")
    ax.bar(x - width, [row.multipole_address_score for row in summaries], width, label="address")
    ax.bar(x, [row.event_order_score for row in summaries], width, label="event")
    ax.bar(x + width, [row.permanence_score for row in summaries], width, label="permanence")
    ax.bar(x + 2.0 * width, [row.generic_smoothness_score for row in summaries], width, label="smoothness")
    ax.set_xticks(x, labels)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("score")
    ax.set_title("Toy supertranslation metrics vs generic smoothness")
    ax.tick_params(axis="x", rotation=25, labelsize=8)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "metric_breakdown.png", dpi=180)
    plt.close(fig)

    target_sample = next(sample for sample in samples if sample.sample_role == "target")
    best_control = max(
        [row for row in summaries if row.sample_role == "negative_control"],
        key=lambda row: row.supertranslation_specificity_score,
    )
    control_sample = next(sample for sample in samples if sample.sample == best_control.sample)
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    for axis, unit, event_time in zip(axes, target.component_units, config.event_times):
        axis.plot(time, projection_trace(target_sample.series, unit), color="tab:green", label="target")
        axis.plot(time, projection_trace(control_sample.series, unit), color="tab:orange", label=best_control.sample)
        axis.axvline(event_time, color="tab:blue", linestyle="--", linewidth=1)
        axis.set_ylabel("projection")
        axis.legend(fontsize=7)
    axes[-1].set_xlabel("normalized time")
    fig.suptitle("Ordered event projections")
    fig.tight_layout()
    fig.savefig(figure_dir / "event_projection_traces.png", dpi=180)
    plt.close(fig)

    recovered_target = recovered_memory(target_sample.series, config)
    recovered_control = recovered_memory(control_sample.series, config)
    target_powers = target.target_band_powers
    recovered_powers = band_power_fractions(recovered_target, target.band_bases)
    control_powers = band_power_fractions(recovered_control, target.band_bases)
    bands = sorted(target.band_bases)
    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.24
    idx = np.arange(len(bands))
    ax.bar(idx - width, [target_powers[band] for band in bands], width, label="true")
    ax.bar(idx, [recovered_powers[band] for band in bands], width, label="target recovered")
    ax.bar(idx + width, [control_powers[band] for band in bands], width, label=best_control.sample)
    ax.set_xticks(idx, [f"l={band}" for band in bands])
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("band power fraction")
    ax.set_title("Multi-pole address power signature")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "band_power_signature.png", dpi=180)
    plt.close(fig)

    lon = np.arctan2(target.points[:, 1], target.points[:, 0])
    lat = np.arcsin(np.clip(target.points[:, 2], -1.0, 1.0))
    plots = [
        ("true total", target.total_field),
        ("target recovered", recovered_target),
        (f"{best_control.sample} recovered", recovered_control),
    ]
    vmax = max(float(np.max(np.abs(values))) for _, values in plots)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharex=True, sharey=True)
    for axis, (title, values) in zip(axes, plots):
        sc = axis.scatter(lon, lat, c=values, cmap="coolwarm", vmin=-vmax, vmax=vmax, s=18)
        axis.set_title(title, fontsize=9)
        axis.set_xlabel("longitude")
    axes[0].set_ylabel("latitude")
    fig.colorbar(sc, ax=axes.ravel().tolist(), shrink=0.85, label="displacement")
    fig.suptitle("Toy supertranslation field identity")
    fig.savefig(figure_dir / "field_identity.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT))


def main() -> None:
    args = parse_args()
    config = make_config(args)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    time, target, samples, summaries, status = run_probe(config)
    write_csv(output_dir / CSV_PATH.name, summaries)
    write_status(output_dir / STATUS_PATH.name, status)
    write_report(output_dir / REPORT_PATH.name, summaries, status)
    write_plots(output_dir, config, time, target, samples, summaries)

    print("Phase 62 multi-pole supertranslation toy probe generated.")
    print(f"CSV: {relative(output_dir / CSV_PATH.name)}")
    print(f"JSON: {relative(output_dir / STATUS_PATH.name)}")
    print(f"Report: {relative(output_dir / REPORT_PATH.name)}")
    print(f"Figures: {relative(output_dir / 'figures')}")
    print(
        f"bridge_status: {status['bridge_status']} ({status['failure_mode']}); "
        f"target={status['target_supertranslation_specificity_score']:.6f}, "
        f"control_max={status['control_max_score']:.6f}"
    )


if __name__ == "__main__":
    main()
