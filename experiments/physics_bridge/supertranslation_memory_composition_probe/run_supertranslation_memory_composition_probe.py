#!/usr/bin/env python3
"""Phase 63 supertranslation + memory composition toy probe.

This sidecar composes two previous toy ideas under a stricter claim gate:

1. a synthetic multi-pole supertranslation-like shift on S2;
2. a later synthetic memory response induced by that shift through a fixed toy
   response map in low-l harmonic coefficient space.

The probe asks whether HAOS-style spectral telemetry can recover the two fields,
their temporal order, and the induced response relation better than controls
that preserve only pieces of the construction.

It does not test real BMS supertranslations, gravitational memory, BMS charges,
celestial holography, Ward identities, or S-matrix data.
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

CSV_PATH = ROOT / "composition_diagnostics.csv"
STATUS_PATH = ROOT / "bridge_status.json"
REPORT_PATH = ROOT / "supertranslation_memory_composition_report.md"

FIELDNAMES = [
    "sample",
    "sample_role",
    "control_type",
    "composition_specificity_score",
    "shift_recoverability",
    "memory_recoverability",
    "composition_relation_score",
    "causal_order_score",
    "permanence_score",
    "address_binding_score",
    "generic_smoothness_score",
    "shift_correlation",
    "memory_correlation",
    "response_correlation",
    "response_amplitude_score",
    "shift_event_time_score",
    "memory_event_time_score",
    "plateau_stability",
]


@dataclass(frozen=True)
class ProbeConfig:
    resolution: int = 192
    time_count: int = 128
    l_max: int = 4
    k_neighbors: int = 10
    seed: int = 6301
    shift_time: float = 0.30
    memory_time: float = 0.66
    event_width: float = 0.030
    memory_gain: float = 0.78
    noise_scale: float = 0.010
    pass_min_specificity_score: float = 0.86
    pass_min_shift_recoverability: float = 0.90
    pass_min_memory_recoverability: float = 0.90
    pass_min_composition: float = 0.88
    pass_min_order: float = 0.82
    pass_min_margin: float = 0.18
    pass_min_separation_z: float = 1.50
    pass_max_control_false_positive_rate: float = 0.0


@dataclass(frozen=True)
class TargetSpec:
    points: np.ndarray
    shift_field: np.ndarray
    memory_field: np.ndarray
    memory_basis: np.ndarray
    band_bases: dict[int, np.ndarray]
    shift_coefficients: np.ndarray
    memory_coefficients: np.ndarray


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
    composition_specificity_score: float
    shift_recoverability: float
    memory_recoverability: float
    composition_relation_score: float
    causal_order_score: float
    permanence_score: float
    address_binding_score: float
    generic_smoothness_score: float
    shift_correlation: float
    memory_correlation: float
    response_correlation: float
    response_amplitude_score: float
    shift_event_time_score: float
    memory_event_time_score: float
    plateau_stability: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 63 supertranslation-memory composition toy probe.")
    parser.add_argument("--resolution", type=int, default=192)
    parser.add_argument("--time-count", type=int, default=128)
    parser.add_argument("--seed", type=int, default=6301)
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


def rms(field: np.ndarray) -> float:
    centered = field - float(np.mean(field))
    return float(np.sqrt(np.mean(centered * centered)))


def cosine_similarity(left: np.ndarray, right: np.ndarray) -> float:
    denom = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denom <= 1.0e-12:
        return 0.0
    return float(np.dot(left, right) / denom)


def pearson_correlation(left: np.ndarray, right: np.ndarray) -> float:
    return cosine_similarity(left - float(np.mean(left)), right - float(np.mean(right)))


def sigmoid_gate(time: np.ndarray, event_time: float, event_width: float) -> np.ndarray:
    return 0.5 * (1.0 + np.tanh((time - event_time) / max(event_width, 1.0e-9)))


def window_indices(config: ProbeConfig) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    time = np.linspace(0.0, 1.0, config.time_count)
    pre = np.flatnonzero(time <= max(0.10, config.shift_time - 0.13))
    mid = np.flatnonzero((time >= config.shift_time + 0.13) & (time <= config.memory_time - 0.13))
    post = np.flatnonzero(time >= min(0.90, config.memory_time + 0.16))
    return pre, mid, post


def recover_shift_and_memory(series: np.ndarray, config: ProbeConfig) -> tuple[np.ndarray, np.ndarray]:
    pre, mid, post = window_indices(config)
    shift = np.mean(series[mid], axis=0) - np.mean(series[pre], axis=0)
    memory = np.mean(series[post], axis=0) - np.mean(series[mid], axis=0)
    return shift, memory


def projection_trace(series: np.ndarray, field: np.ndarray) -> np.ndarray:
    return series @ field / max(float(np.dot(field, field)), 1.0e-12)


def response_coefficients(coefficients: np.ndarray) -> np.ndarray:
    c2 = coefficients[:5]
    c3 = coefficients[5:12]
    c4 = coefficients[12:21]
    r2 = 0.22 * c2 + 0.31 * c3[:5] - 0.24 * c4[:5]
    r3 = 0.19 * c3.copy()
    r3[:5] += 0.34 * c2
    r3 += 0.17 * c4[:7]
    r4 = 0.16 * c4.copy()
    r4[:7] += 0.29 * c3
    r4[:5] -= 0.21 * c2
    return np.concatenate([r2, r3, r4])


def response_from_shift(shift_field: np.ndarray, target: TargetSpec, config: ProbeConfig) -> np.ndarray:
    coefficients = target.memory_basis.T @ shift_field
    raw = target.memory_basis @ response_coefficients(coefficients)
    raw = raw - float(np.dot(raw, shift_field) / max(float(np.dot(shift_field, shift_field)), 1.0e-12)) * shift_field
    response = normalize_field(raw)
    return response * (config.memory_gain * max(rms(shift_field), 1.0e-12))


def target_spec(points: np.ndarray, config: ProbeConfig) -> TargetSpec:
    subspaces = harmonic_subspaces(points, config.l_max)
    band_bases = {band: subspaces[band] for band in (2, 3, 4)}
    memory_basis = np.column_stack([band_bases[2], band_bases[3], band_bases[4]])
    c2 = np.array([0.68, -0.35, 0.52, 0.29, -0.44], dtype=float)
    c3 = np.array([0.31, -0.41, 0.24, 0.35, -0.29, 0.17, 0.27], dtype=float)
    c4 = np.array([0.22, -0.16, 0.34, -0.24, 0.19, 0.31, -0.13, 0.23, -0.18], dtype=float)
    shift_raw = band_bases[2] @ c2 + 0.78 * (band_bases[3] @ c3) + 0.55 * (band_bases[4] @ c4)
    shift = normalize_field(shift_raw)
    placeholder = TargetSpec(
        points=points,
        shift_field=shift,
        memory_field=np.zeros_like(shift),
        memory_basis=memory_basis,
        band_bases=band_bases,
        shift_coefficients=memory_basis.T @ shift,
        memory_coefficients=np.zeros(memory_basis.shape[1], dtype=float),
    )
    memory = response_from_shift(shift, placeholder, config)
    return TargetSpec(
        points=points,
        shift_field=shift,
        memory_field=memory,
        memory_basis=memory_basis,
        band_bases=band_bases,
        shift_coefficients=memory_basis.T @ shift,
        memory_coefficients=memory_basis.T @ memory,
    )


def random_low_mode_pattern(
    rng: np.random.Generator,
    basis: np.ndarray,
    target: np.ndarray | None = None,
) -> np.ndarray:
    coefficients = rng.normal(size=basis.shape[1])
    field = basis @ coefficients
    if target is not None:
        field = field - float(np.dot(field, target) / max(float(np.dot(target, target)), 1.0e-12)) * target
    return normalize_field(field)


def build_series(
    config: ProbeConfig,
    shift_field: np.ndarray,
    memory_field: np.ndarray,
    rng: np.random.Generator,
    shift_time: float | None = None,
    memory_time: float | None = None,
    noise_scale: float | None = None,
) -> np.ndarray:
    time = np.linspace(0.0, 1.0, config.time_count)
    shift_gate = sigmoid_gate(time, config.shift_time if shift_time is None else shift_time, config.event_width)
    memory_gate = sigmoid_gate(time, config.memory_time if memory_time is None else memory_time, config.event_width)
    series = shift_gate[:, None] * shift_field[None, :] + memory_gate[:, None] * memory_field[None, :]
    scale = config.noise_scale if noise_scale is None else noise_scale
    series += rng.normal(scale=scale, size=series.shape)
    return series


def build_samples(config: ProbeConfig) -> tuple[np.ndarray, TargetSpec, list[SampleSpec]]:
    rng = np.random.default_rng(config.seed)
    points = fibonacci_sphere(config.resolution)
    target = target_spec(points, config)
    time = np.linspace(0.0, 1.0, config.time_count)
    target_series = build_series(config, target.shift_field, target.memory_field, rng)

    wrong_memory = random_low_mode_pattern(rng, target.memory_basis, target=target.memory_field)
    wrong_memory *= rms(target.memory_field) / max(rms(wrong_memory), 1.0e-12)
    independent_response = build_series(config, target.shift_field, wrong_memory, rng)

    wrong_shift = random_low_mode_pattern(rng, target.memory_basis, target=target.shift_field)
    wrong_shift *= rms(target.shift_field) / max(rms(wrong_shift), 1.0e-12)
    right_memory_wrong_shift = build_series(config, wrong_shift, target.memory_field, rng)

    simultaneous_final = build_series(
        config,
        target.shift_field + target.memory_field,
        np.zeros_like(target.memory_field),
        rng,
        shift_time=config.shift_time,
        memory_time=config.memory_time,
    )

    reversed_order = build_series(
        config,
        target.memory_field,
        target.shift_field,
        rng,
        shift_time=config.shift_time,
        memory_time=config.memory_time,
    )

    shift_only = build_series(config, target.shift_field, np.zeros_like(target.memory_field), rng)
    memory_only = build_series(config, np.zeros_like(target.shift_field), target.memory_field, rng)

    permutation = rng.permutation(config.resolution)
    coordinate_permutation = build_series(
        config,
        target.shift_field[permutation],
        target.memory_field[permutation],
        rng,
    )

    transient = np.zeros_like(target_series)
    for event_time, field in ((config.shift_time, target.shift_field), (config.memory_time, target.memory_field)):
        burst = np.exp(-((time - event_time) / 0.055) ** 2)
        transient += burst[:, None] * field[None, :]
    transient += rng.normal(scale=config.noise_scale, size=transient.shape)

    smooth_basis = np.column_stack([harmonic_subspaces(points, config.l_max)[1], target.band_bases[2]])
    smooth_shift = random_low_mode_pattern(rng, smooth_basis, target=target.shift_field)
    smooth_memory = random_low_mode_pattern(rng, smooth_basis, target=target.memory_field)
    smooth_shift *= rms(target.shift_field) / max(rms(smooth_shift), 1.0e-12)
    smooth_memory *= rms(target.memory_field) / max(rms(smooth_memory), 1.0e-12)
    smooth_decoy = build_series(config, smooth_shift, smooth_memory, rng)

    return time, target, [
        SampleSpec("composition_target", "target", "none", target_series),
        SampleSpec("independent_memory_control", "negative_control", "right_shift_wrong_response", independent_response),
        SampleSpec("wrong_shift_right_memory_control", "negative_control", "wrong_shift_right_memory", right_memory_wrong_shift),
        SampleSpec("simultaneous_final_control", "negative_control", "final_field_without_two_stage_composition", simultaneous_final),
        SampleSpec("reversed_order_control", "negative_control", "memory_before_shift", reversed_order),
        SampleSpec("shift_only_control", "negative_control", "supertranslation_without_memory", shift_only),
        SampleSpec("memory_only_control", "negative_control", "memory_without_shift", memory_only),
        SampleSpec("coordinate_permutation_control", "negative_control", "node_identity_permutation", coordinate_permutation),
        SampleSpec("transient_composition_control", "negative_control", "bursts_without_permanent_composition", transient),
        SampleSpec("smooth_two_step_decoy", "negative_control", "smooth_wrong_address_two_step", smooth_decoy),
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


def field_recoverability(recovered: np.ndarray, truth: np.ndarray) -> tuple[float, float, float]:
    corr = cosine_similarity(recovered, truth)
    corr_score = float(np.clip(max(corr, 0.0), 0.0, 1.0))
    amplitude = float(np.linalg.norm(recovered) / max(float(np.linalg.norm(truth)), 1.0e-12))
    amplitude_score = float(np.clip(1.0 - abs(amplitude - 1.0) / 0.35, 0.0, 1.0))
    return float(0.74 * corr_score + 0.26 * amplitude_score), corr, amplitude_score


def composition_relation_score(
    recovered_shift: np.ndarray,
    recovered_memory: np.ndarray,
    target: TargetSpec,
    config: ProbeConfig,
) -> tuple[float, float, float]:
    predicted = response_from_shift(recovered_shift, target, config)
    corr = cosine_similarity(recovered_memory, predicted)
    corr_score = float(np.clip(max(corr, 0.0), 0.0, 1.0))
    amplitude = float(np.linalg.norm(recovered_memory) / max(float(np.linalg.norm(predicted)), 1.0e-12))
    amplitude_score = float(np.clip(1.0 - abs(amplitude - 1.0) / 0.35, 0.0, 1.0))
    return float(0.78 * corr_score + 0.22 * amplitude_score), corr, amplitude_score


def address_binding_score(
    recovered_shift: np.ndarray,
    recovered_memory: np.ndarray,
    target: TargetSpec,
) -> float:
    shift_coeffs = target.memory_basis.T @ recovered_shift
    memory_coeffs = target.memory_basis.T @ recovered_memory
    shift_score = float(np.clip(max(cosine_similarity(shift_coeffs, target.shift_coefficients), 0.0), 0.0, 1.0))
    memory_score = float(np.clip(max(cosine_similarity(memory_coeffs, target.memory_coefficients), 0.0), 0.0, 1.0))
    return float(0.52 * shift_score + 0.48 * memory_score)


def causal_order_score(config: ProbeConfig, series: np.ndarray, target: TargetSpec) -> tuple[float, float, float]:
    time = np.linspace(0.0, 1.0, config.time_count)
    shift_gate = sigmoid_gate(time, config.shift_time, config.event_width)
    memory_gate = sigmoid_gate(time, config.memory_time, config.event_width)
    shift_trace = projection_trace(series, target.shift_field)
    memory_trace = projection_trace(series, target.memory_field)
    shift_derivative = np.gradient(shift_trace, time)
    memory_derivative = np.gradient(memory_trace, time)
    shift_detected = float(time[int(np.argmax(shift_derivative))])
    memory_detected = float(time[int(np.argmax(memory_derivative))])
    shift_time_score = float(np.clip(1.0 - abs(shift_detected - config.shift_time) / 0.10, 0.0, 1.0))
    memory_time_score = float(np.clip(1.0 - abs(memory_detected - config.memory_time) / 0.10, 0.0, 1.0))
    order_gap = memory_detected - shift_detected
    gap_score = float(np.clip((order_gap - 0.18) / 0.14, 0.0, 1.0))
    shift_corr = float(np.clip(max(pearson_correlation(shift_trace, shift_gate), 0.0), 0.0, 1.0))
    memory_corr = float(np.clip(max(pearson_correlation(memory_trace, memory_gate), 0.0), 0.0, 1.0))
    score = float(
        0.26 * shift_time_score
        + 0.26 * memory_time_score
        + 0.22 * gap_score
        + 0.13 * shift_corr
        + 0.13 * memory_corr
    )
    return score, shift_time_score, memory_time_score


def permanence_score(config: ProbeConfig, series: np.ndarray, target: TargetSpec) -> tuple[float, float]:
    pre, mid, post = window_indices(config)
    shift_trace = projection_trace(series, target.shift_field)
    memory_trace = projection_trace(series, target.memory_field)
    shift_step = abs(float(np.mean(shift_trace[mid]) - np.mean(shift_trace[pre])))
    memory_step = abs(float(np.mean(memory_trace[post]) - np.mean(memory_trace[mid])))
    shift_stability = float(np.std(shift_trace[mid]) / max(shift_step, 1.0e-12))
    post_shift_stability = float(np.std(shift_trace[post]) / max(shift_step, 1.0e-12))
    memory_stability = float(np.std(memory_trace[post]) / max(memory_step, 1.0e-12))
    stability = float(np.clip(1.0 - (shift_stability + post_shift_stability + memory_stability) / 0.20, 0.0, 1.0))
    contrast = float(np.clip(0.5 * shift_step + 0.5 * memory_step / max(config.memory_gain, 1.0e-12), 0.0, 1.0))
    return float(0.72 * stability + 0.28 * contrast), stability


def evaluate_sample(
    config: ProbeConfig,
    sample: SampleSpec,
    target: TargetSpec,
    adjacency: np.ndarray,
) -> SampleSummary:
    recovered_shift, recovered_memory = recover_shift_and_memory(sample.series, config)
    shift_score, shift_corr, _ = field_recoverability(recovered_shift, target.shift_field)
    memory_score, memory_corr, _ = field_recoverability(recovered_memory, target.memory_field)
    composition, response_corr, response_amp = composition_relation_score(
        recovered_shift,
        recovered_memory,
        target,
        config,
    )
    order, shift_time, memory_time = causal_order_score(config, sample.series, target)
    permanence, plateau = permanence_score(config, sample.series, target)
    address = address_binding_score(recovered_shift, recovered_memory, target)
    smoothness = graph_smoothness_score(adjacency, recovered_shift + recovered_memory)
    specificity = float(
        0.19 * shift_score
        + 0.19 * memory_score
        + 0.27 * composition
        + 0.18 * order
        + 0.10 * address
        + 0.05 * permanence
        + 0.02 * smoothness
    )
    return SampleSummary(
        sample=sample.sample,
        sample_role=sample.sample_role,
        control_type=sample.control_type,
        composition_specificity_score=specificity,
        shift_recoverability=shift_score,
        memory_recoverability=memory_score,
        composition_relation_score=composition,
        causal_order_score=order,
        permanence_score=permanence,
        address_binding_score=address,
        generic_smoothness_score=smoothness,
        shift_correlation=shift_corr,
        memory_correlation=memory_corr,
        response_correlation=response_corr,
        response_amplitude_score=response_amp,
        shift_event_time_score=shift_time,
        memory_event_time_score=memory_time,
        plateau_stability=plateau,
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
    control_scores = np.array([row.composition_specificity_score for row in controls], dtype=float)
    control_threshold = target.composition_specificity_score - config.pass_min_margin
    control_false_positive_rate = float(np.mean(control_scores >= control_threshold))
    control_max = max(controls, key=lambda row: row.composition_specificity_score)
    score_margin = float(target.composition_specificity_score - control_max.composition_specificity_score)
    separation_z = float(
        (target.composition_specificity_score - float(np.mean(control_scores)))
        / max(float(np.std(control_scores)), 1.0e-9)
    )
    target_passes = bool(
        target.composition_specificity_score >= config.pass_min_specificity_score
        and target.shift_recoverability >= config.pass_min_shift_recoverability
        and target.memory_recoverability >= config.pass_min_memory_recoverability
        and target.composition_relation_score >= config.pass_min_composition
        and target.causal_order_score >= config.pass_min_order
    )
    controls_pass = bool(
        control_false_positive_rate <= config.pass_max_control_false_positive_rate
        and score_margin >= config.pass_min_margin
        and separation_z >= config.pass_min_separation_z
    )
    if target_passes and controls_pass:
        bridge_status = "PASS"
        failure_mode = "NONE"
    elif target.composition_specificity_score >= 0.72 and separation_z >= 1.0:
        bridge_status = "MARGINAL"
        failure_mode = "COMPOSITION_CONTROL_SEPARATION_INCOMPLETE"
    else:
        bridge_status = "OPEN"
        failure_mode = "COMPOSED_TOY_SIGNATURE_NOT_DISTINCT"
    return {
        "phase": "63",
        "bridge_name": "supertranslation_memory_composition_probe",
        "bridge_status": bridge_status,
        "failure_mode": failure_mode,
        "success_definition": (
            "PASS requires a synthetic supertranslation-like shift and later induced memory "
            "response to recover both fields, their temporal order, and the fixed toy response "
            "relation while controls that preserve only pieces fail."
        ),
        "config": asdict(config),
        "target_composition_specificity_score": float(target.composition_specificity_score),
        "target_shift_recoverability": float(target.shift_recoverability),
        "target_memory_recoverability": float(target.memory_recoverability),
        "target_composition_relation_score": float(target.composition_relation_score),
        "target_causal_order_score": float(target.causal_order_score),
        "target_permanence_score": float(target.permanence_score),
        "target_address_binding_score": float(target.address_binding_score),
        "target_generic_smoothness_score": float(target.generic_smoothness_score),
        "control_max_score": float(control_max.composition_specificity_score),
        "control_max_sample": control_max.sample,
        "control_false_positive_threshold": float(control_threshold),
        "control_false_positive_rate": control_false_positive_rate,
        "score_margin_over_best_control": score_margin,
        "separation_z": separation_z,
        "core_modified": False,
        "safe_to_continue": True,
        "next_phase": "64_claim_gated_celestial_toy_summary_or_asymptotic_scaling",
        "non_claims": [
            "no real BMS supertranslation recovery claim",
            "no real gravitational-memory observable claim",
            "no BMS charge or soft-hair recovery claim",
            "no celestial holography recovery claim",
            "no soft theorem or S-matrix reconstruction claim",
            "synthetic two-stage composition data only",
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
        "| {sample} | {sample_role} | {control_type} | {composition_specificity_score:.6f} | {shift_recoverability:.6f} | {memory_recoverability:.6f} | {composition_relation_score:.6f} | {causal_order_score:.6f} | {address_binding_score:.6f} |".format(
            **asdict(row)
        )
        for row in summaries
    )
    controls = [row for row in summaries if row.sample_role == "negative_control"]
    leakage_rows = "\n".join(
        f"- `{row.sample}`: composition={row.composition_specificity_score:.6f}, "
        f"shift={row.shift_recoverability:.6f}, memory={row.memory_recoverability:.6f}, "
        f"relation={row.composition_relation_score:.6f}, order={row.causal_order_score:.6f}"
        for row in sorted(controls, key=lambda item: item.composition_specificity_score, reverse=True)
    )
    report = f"""# Phase 63 Supertranslation + Memory Composition Toy Probe

## Scope

This file is generated by `experiments/physics_bridge/supertranslation_memory_composition_probe/run_supertranslation_memory_composition_probe.py`.

This is a synthetic two-stage benchmark on a sampled S2 boundary. A multi-pole supertranslation-like shift occurs first; a later memory-like response is induced by that shift through a fixed toy response map in harmonic coefficient space.

It does not use real BMS charges, gravitational-wave memory data, celestial amplitudes, Ward identities, or S-matrix inputs.

The narrow question is:

> Can HAOS-style spectral telemetry distinguish a composed shift -> memory relation from controls that preserve only the shift, only the memory, the final field, generic smoothness, or the wrong response binding?

## Bridge Status

- `bridge_status`: `{status['bridge_status']}`
- `failure_mode`: `{status['failure_mode']}`
- target composition score: `{status['target_composition_specificity_score']:.6f}`
- target shift recoverability: `{status['target_shift_recoverability']:.6f}`
- target memory recoverability: `{status['target_memory_recoverability']:.6f}`
- target composition relation: `{status['target_composition_relation_score']:.6f}`
- target causal order: `{status['target_causal_order_score']:.6f}`
- target permanence: `{status['target_permanence_score']:.6f}`
- target address binding: `{status['target_address_binding_score']:.6f}`
- target generic smoothness: `{status['target_generic_smoothness_score']:.6f}`
- control max sample: `{status['control_max_sample']}`
- control max score: `{status['control_max_score']:.6f}`
- control false-positive threshold: `{status['control_false_positive_threshold']:.6f}`
- control false-positive rate: `{status['control_false_positive_rate']:.6f}`
- separation z: `{status['separation_z']:.6f}`
- margin over best control: `{status['score_margin_over_best_control']:.6f}`

## Summary Table

| sample | role | control | composition_score | shift | memory | relation | order | address |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
{summary_rows}

## Control Leakage Read

{leakage_rows}

## Composition Gate

The decisive row is `composition_relation_score`: the recovered memory must be
predicted by the recovered shift under the fixed toy response map. A control can
therefore preserve smoothness, final-field energy, or one of the two fields and
still fail if the shift -> memory binding is absent.

## Interpretation

`PASS` means the toy target is separable from controls on shift recovery, memory recovery, temporal order, and induced response relation.

`MARGINAL` means the composed target is detected, but at least one control remains too competitive.

`OPEN` means the telemetry is mostly measuring smoothness, low-mode content, or final-field recovery rather than composed dynamics.

## Claim Gate

This sidecar remains inside the Phase 60 claim gate. Even under `PASS`, the established physics rows for BMS supertranslations, BMS charge recovery, real gravitational memory, celestial holography, and real soft theorems remain `OPEN`.

Allowed language:

- supertranslation-memory composition toy proxy
- synthetic two-stage S2 deformation benchmark
- toy shift-to-memory response recovery

Disallowed language:

- BMS supertranslations recovered
- gravitational memory recovered
- soft hair detected
- celestial holography validated

## Authority

- CSV: `experiments/physics_bridge/supertranslation_memory_composition_probe/composition_diagnostics.csv`
- JSON: `experiments/physics_bridge/supertranslation_memory_composition_probe/bridge_status.json`
- Figures: `experiments/physics_bridge/supertranslation_memory_composition_probe/figures/`
- Generator: `experiments/physics_bridge/supertranslation_memory_composition_probe/run_supertranslation_memory_composition_probe.py`
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

    fig, ax = plt.subplots(figsize=(13, 5))
    colors = ["tab:green" if row.sample_role == "target" else "tab:orange" for row in summaries]
    ax.bar(labels, [row.composition_specificity_score for row in summaries], color=colors)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("composition specificity")
    ax.set_title("Phase 63 composed target vs controls")
    ax.tick_params(axis="x", rotation=25, labelsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "control_comparison.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(13, 5))
    width = 0.14
    ax.bar(x - 2.5 * width, [row.shift_recoverability for row in summaries], width, label="shift")
    ax.bar(x - 1.5 * width, [row.memory_recoverability for row in summaries], width, label="memory")
    ax.bar(x - 0.5 * width, [row.composition_relation_score for row in summaries], width, label="relation")
    ax.bar(x + 0.5 * width, [row.causal_order_score for row in summaries], width, label="order")
    ax.bar(x + 1.5 * width, [row.address_binding_score for row in summaries], width, label="address")
    ax.bar(x + 2.5 * width, [row.generic_smoothness_score for row in summaries], width, label="smoothness")
    ax.set_xticks(x, labels)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("score")
    ax.set_title("Composition metrics vs generic smoothness")
    ax.tick_params(axis="x", rotation=25, labelsize=8)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "metric_breakdown.png", dpi=180)
    plt.close(fig)

    target_sample = next(sample for sample in samples if sample.sample_role == "target")
    best_control = max(
        [row for row in summaries if row.sample_role == "negative_control"],
        key=lambda row: row.composition_specificity_score,
    )
    control_sample = next(sample for sample in samples if sample.sample == best_control.sample)
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    axes[0].plot(time, projection_trace(target_sample.series, target.shift_field), color="tab:green", label="target")
    axes[0].plot(time, projection_trace(control_sample.series, target.shift_field), color="tab:orange", label=best_control.sample)
    axes[0].axvline(config.shift_time, color="tab:blue", linestyle="--", linewidth=1)
    axes[0].set_ylabel("shift projection")
    axes[0].legend(fontsize=8)
    axes[1].plot(time, projection_trace(target_sample.series, target.memory_field), color="tab:green", label="target")
    axes[1].plot(time, projection_trace(control_sample.series, target.memory_field), color="tab:orange", label=best_control.sample)
    axes[1].axvline(config.memory_time, color="tab:blue", linestyle="--", linewidth=1)
    axes[1].set_ylabel("memory projection")
    axes[1].set_xlabel("normalized time")
    axes[1].legend(fontsize=8)
    fig.suptitle("Two-stage composition traces")
    fig.tight_layout()
    fig.savefig(figure_dir / "composition_projection_traces.png", dpi=180)
    plt.close(fig)

    recovered_shift, recovered_memory = recover_shift_and_memory(target_sample.series, config)
    predicted_memory = response_from_shift(recovered_shift, target, config)
    control_shift, control_memory = recover_shift_and_memory(control_sample.series, config)
    control_predicted = response_from_shift(control_shift, target, config)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(predicted_memory, recovered_memory, s=12, alpha=0.75, label="target")
    ax.scatter(control_predicted, control_memory, s=12, alpha=0.55, label=best_control.sample)
    lim = max(
        float(np.max(np.abs(predicted_memory))),
        float(np.max(np.abs(recovered_memory))),
        float(np.max(np.abs(control_predicted))),
        float(np.max(np.abs(control_memory))),
    )
    ax.plot([-lim, lim], [-lim, lim], color="black", linewidth=1, linestyle="--")
    ax.set_xlabel("predicted memory from recovered shift")
    ax.set_ylabel("recovered memory")
    ax.set_title("Composition response relation")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "composition_relation.png", dpi=180)
    plt.close(fig)

    lon = np.arctan2(target.points[:, 1], target.points[:, 0])
    lat = np.arcsin(np.clip(target.points[:, 2], -1.0, 1.0))
    plots = [
        ("true shift", target.shift_field),
        ("true memory", target.memory_field),
        ("target recovered memory", recovered_memory),
        (f"{best_control.sample} memory", control_memory),
    ]
    vmax = max(float(np.max(np.abs(values))) for _, values in plots)
    fig, axes = plt.subplots(1, 4, figsize=(15, 4), sharex=True, sharey=True)
    for axis, (title, values) in zip(axes, plots):
        sc = axis.scatter(lon, lat, c=values, cmap="coolwarm", vmin=-vmax, vmax=vmax, s=16)
        axis.set_title(title, fontsize=9)
        axis.set_xlabel("longitude")
    axes[0].set_ylabel("latitude")
    fig.colorbar(sc, ax=axes.ravel().tolist(), shrink=0.82, label="field")
    fig.suptitle("Shift/memory field identity")
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

    print("Phase 63 supertranslation-memory composition toy probe generated.")
    print(f"CSV: {relative(output_dir / CSV_PATH.name)}")
    print(f"JSON: {relative(output_dir / STATUS_PATH.name)}")
    print(f"Report: {relative(output_dir / REPORT_PATH.name)}")
    print(f"Figures: {relative(output_dir / 'figures')}")
    print(
        f"bridge_status: {status['bridge_status']} ({status['failure_mode']}); "
        f"target={status['target_composition_specificity_score']:.6f}, "
        f"control_max={status['control_max_score']:.6f}"
    )


if __name__ == "__main__":
    main()
