#!/usr/bin/env python3
"""Phase 61 toy gravitational-memory / soft-hair probe.

This sidecar is deliberately claim-gated. It creates a synthetic permanent
displacement-memory pattern on a sampled S2 boundary and asks whether
HAOS-style spectral telemetry can distinguish that known target from controls
that preserve smoothness, time-series scale, or low-mode content while breaking
the specific memory signature.

It does not test a real gravitational-memory observable, BMS charge, celestial
amplitude, or S-matrix statement.
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

CSV_PATH = ROOT / "memory_diagnostics.csv"
STATUS_PATH = ROOT / "bridge_status.json"
REPORT_PATH = ROOT / "gravitational_memory_toy_report.md"
FIGURE_DIR = ROOT / "figures"

FIELDNAMES = [
    "sample",
    "sample_role",
    "control_type",
    "memory_specificity_score",
    "memory_recoverability",
    "permanence_score",
    "temporal_step_score",
    "mode_signature_score",
    "generic_smoothness_score",
    "recovered_target_correlation",
    "recovered_amplitude",
    "plateau_stability",
    "low_mode_energy_fraction",
    "target_mode_correlation",
]


@dataclass(frozen=True)
class MemoryConfig:
    resolution: int = 144
    time_count: int = 96
    l_max: int = 4
    k_neighbors: int = 10
    seed: int = 6101
    event_time: float = 0.42
    event_width: float = 0.035
    noise_scale: float = 0.012
    pass_min_memory_score: float = 0.86
    pass_min_recoverability: float = 0.90
    pass_min_permanence: float = 0.82
    pass_min_temporal_step: float = 0.86
    pass_min_mode_signature: float = 0.84
    pass_min_margin: float = 0.18
    pass_min_separation_z: float = 1.50
    pass_max_control_false_positive_rate: float = 0.0


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
    memory_specificity_score: float
    memory_recoverability: float
    permanence_score: float
    temporal_step_score: float
    mode_signature_score: float
    generic_smoothness_score: float
    recovered_target_correlation: float
    recovered_amplitude: float
    plateau_stability: float
    low_mode_energy_fraction: float
    target_mode_correlation: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 61 toy gravitational-memory probe.")
    parser.add_argument("--resolution", type=int, default=144)
    parser.add_argument("--time-count", type=int, default=96)
    parser.add_argument("--seed", type=int, default=6101)
    parser.add_argument("--noise-scale", type=float, default=0.012)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT,
        help="Directory for generated CSV, JSON, markdown report, and figures.",
    )
    return parser.parse_args()


def make_config(args: argparse.Namespace) -> MemoryConfig:
    return MemoryConfig(
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
    points = np.column_stack((radius * np.cos(phi), radius * np.sin(phi), z))
    return points.astype(float)


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
        weights = np.exp(-((distances[i, neighbors] / scale) ** 2))
        adjacency[i, neighbors] = weights
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


def target_memory_pattern(points: np.ndarray, l_max: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    subspaces = harmonic_subspaces(points, l_max)
    basis = np.column_stack([subspaces[2], subspaces[3]])
    coeffs = np.array(
        [
            0.82,
            -0.36,
            0.28,
            0.51,
            -0.43,
            0.31,
            -0.24,
            0.18,
            -0.29,
            0.22,
            0.16,
            -0.20,
        ],
        dtype=float,
    )
    field = normalize_field(basis @ coeffs)
    target_coeffs = basis.T @ field
    return field, basis, target_coeffs


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


def projection_trace(series: np.ndarray, truth: np.ndarray) -> np.ndarray:
    return series @ truth / max(float(np.dot(truth, truth)), 1.0e-12)


def pre_post_indices(config: MemoryConfig) -> tuple[np.ndarray, np.ndarray]:
    time = np.linspace(0.0, 1.0, config.time_count)
    pre = np.flatnonzero(time <= max(0.16, config.event_time - 0.20))
    post = np.flatnonzero(time >= min(0.84, config.event_time + 0.32))
    return pre, post


def recovered_memory(series: np.ndarray, config: MemoryConfig) -> np.ndarray:
    pre, post = pre_post_indices(config)
    return np.mean(series[post], axis=0) - np.mean(series[pre], axis=0)


def build_target_series(
    config: MemoryConfig,
    rng: np.random.Generator,
    truth: np.ndarray,
    high_mode_basis: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    time = np.linspace(0.0, 1.0, config.time_count)
    gate = sigmoid_gate(time, config.event_time, config.event_width)
    burst = np.exp(-((time - config.event_time) / 0.055) ** 2)
    burst_field = random_low_mode_pattern(rng, high_mode_basis, target=truth)
    noise = rng.normal(scale=config.noise_scale, size=(config.time_count, truth.shape[0]))
    series = gate[:, None] * truth[None, :] + 0.10 * burst[:, None] * burst_field[None, :] + noise
    return time, gate, series


def build_samples(
    config: MemoryConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[SampleSpec]]:
    rng = np.random.default_rng(config.seed)
    points = fibonacci_sphere(config.resolution)
    truth, memory_basis, target_coeffs = target_memory_pattern(points, config.l_max)
    subspaces = harmonic_subspaces(points, config.l_max)
    high_mode_basis = subspaces[4]
    time, gate, target_series = build_target_series(config, rng, truth, high_mode_basis)

    target_projection = projection_trace(target_series, truth)

    temporal_shuffle = target_series.copy()
    for node in range(temporal_shuffle.shape[1]):
        temporal_shuffle[:, node] = temporal_shuffle[rng.permutation(config.time_count), node]

    transient_field = normalize_field(truth + 0.30 * random_low_mode_pattern(rng, high_mode_basis, target=truth))
    transient = np.exp(-((time - config.event_time) / 0.075) ** 2)
    transient_decay = transient[:, None] * transient_field[None, :]
    transient_decay += rng.normal(scale=config.noise_scale, size=target_series.shape)

    permuted_truth = truth[rng.permutation(len(truth))]
    coordinate_permutation = gate[:, None] * permuted_truth[None, :]
    coordinate_permutation += rng.normal(scale=config.noise_scale, size=target_series.shape)

    mode_scramble = random_low_mode_pattern(rng, memory_basis, target=truth)
    mode_scramble_series = gate[:, None] * mode_scramble[None, :]
    mode_scramble_series += rng.normal(scale=config.noise_scale, size=target_series.shape)

    low_basis = np.column_stack([subspaces[1], subspaces[2]])
    low_mode_decoy = random_low_mode_pattern(rng, low_basis, target=truth)
    low_mode_decoy_series = gate[:, None] * low_mode_decoy[None, :]
    low_mode_decoy_series += rng.normal(scale=config.noise_scale, size=target_series.shape)

    drift_coeffs = rng.normal(scale=0.020, size=(config.time_count, memory_basis.shape[1]))
    drift_coeffs = np.cumsum(drift_coeffs, axis=0)
    drift = drift_coeffs @ memory_basis.T
    drift = drift / max(float(np.std(drift)), 1.0e-12) * 0.55
    drift += 0.30 * target_projection[:, None] * random_low_mode_pattern(rng, memory_basis, target=truth)[None, :]

    samples = [
        SampleSpec("toy_memory_target", "target", "none", target_series),
        SampleSpec(
            "temporal_shuffle_control",
            "negative_control",
            "per_node_time_shuffle",
            temporal_shuffle,
        ),
        SampleSpec(
            "transient_decay_control",
            "negative_control",
            "burst_without_permanent_shift",
            transient_decay,
        ),
        SampleSpec(
            "coordinate_permutation_control",
            "negative_control",
            "node_identity_permutation",
            coordinate_permutation,
        ),
        SampleSpec(
            "mode_scramble_control",
            "negative_control",
            "same_low_band_wrong_memory_coefficients",
            mode_scramble_series,
        ),
        SampleSpec(
            "low_mode_decoy_control",
            "negative_control",
            "smooth_step_wrong_l_band",
            low_mode_decoy_series,
        ),
        SampleSpec(
            "stochastic_drift_control",
            "negative_control",
            "smooth_random_walk_no_event_lock",
            drift,
        ),
    ]
    return points, time, gate, truth, target_coeffs, samples


def graph_smoothness_score(adjacency: np.ndarray, field: np.ndarray) -> float:
    upper = np.argwhere(np.triu(adjacency, k=1) > 1.0e-12)
    if len(upper) == 0:
        return 0.0
    weights = np.array([adjacency[i, j] for i, j in upper], dtype=float)
    diffs = np.array([field[i] - field[j] for i, j in upper], dtype=float)
    energy = float(np.sum(weights * diffs * diffs) / max(float(np.sum(weights)), 1.0e-12))
    variance = float(np.var(field))
    normalized = energy / max(variance, 1.0e-12)
    return float(np.clip(1.0 - normalized / 0.60, 0.0, 1.0))


def memory_recoverability_score(recovered: np.ndarray, truth: np.ndarray) -> tuple[float, float, float]:
    corr = cosine_similarity(recovered, truth)
    corr_score = float(np.clip(max(corr, 0.0), 0.0, 1.0))
    amplitude = float(np.linalg.norm(recovered) / max(float(np.linalg.norm(truth)), 1.0e-12))
    amplitude_score = float(np.clip(1.0 - abs(amplitude - 1.0) / 0.35, 0.0, 1.0))
    score = float(0.75 * corr_score + 0.25 * amplitude_score)
    return score, corr, amplitude


def permanence_score(
    config: MemoryConfig,
    series: np.ndarray,
    recovered: np.ndarray,
    truth: np.ndarray,
) -> tuple[float, float]:
    pre, post = pre_post_indices(config)
    trace = projection_trace(series, truth)
    shift = abs(float(np.mean(trace[post]) - np.mean(trace[pre])))
    post_stability = float(np.std(trace[post]) / max(shift, 1.0e-12))
    pre_stability = float(np.std(trace[pre]) / max(shift, 1.0e-12))
    stability = float(np.clip(1.0 - (post_stability + pre_stability) / 0.12, 0.0, 1.0))
    contrast = float(np.clip(shift / 0.82, 0.0, 1.0))
    residual_norm = float(np.linalg.norm(recovered) / max(float(np.linalg.norm(truth)), 1.0e-12))
    residual_score = float(np.clip(residual_norm / 0.82, 0.0, 1.0))
    score = float(0.50 * stability + 0.30 * contrast + 0.20 * residual_score)
    return score, stability


def temporal_step_score(
    series: np.ndarray,
    truth: np.ndarray,
    gate: np.ndarray,
    config: MemoryConfig,
) -> float:
    pre, post = pre_post_indices(config)
    trace = projection_trace(series, truth)
    corr = pearson_correlation(trace, gate)
    contrast = float(np.mean(trace[post]) - np.mean(trace[pre]))
    contrast_score = float(np.clip(contrast / 0.82, 0.0, 1.0))
    return float(0.72 * np.clip(max(corr, 0.0), 0.0, 1.0) + 0.28 * contrast_score)


def mode_signature_score(
    recovered: np.ndarray,
    memory_basis: np.ndarray,
    target_coeffs: np.ndarray,
) -> tuple[float, float, float]:
    coeffs = memory_basis.T @ recovered
    projected = memory_basis @ coeffs
    low_mode_fraction = float(
        np.linalg.norm(projected) / max(float(np.linalg.norm(recovered)), 1.0e-12)
    )
    target_corr = cosine_similarity(coeffs, target_coeffs)
    corr_score = float(np.clip(max(target_corr, 0.0), 0.0, 1.0))
    fraction_score = float(np.clip(low_mode_fraction, 0.0, 1.0))
    score = float(0.72 * corr_score + 0.28 * fraction_score)
    return score, low_mode_fraction, target_corr


def evaluate_sample(
    config: MemoryConfig,
    sample: SampleSpec,
    truth: np.ndarray,
    gate: np.ndarray,
    memory_basis: np.ndarray,
    target_coeffs: np.ndarray,
    adjacency: np.ndarray,
) -> SampleSummary:
    recovered = recovered_memory(sample.series, config)
    recoverability, corr, amplitude = memory_recoverability_score(recovered, truth)
    permanence, plateau_stability = permanence_score(config, sample.series, recovered, truth)
    temporal = temporal_step_score(sample.series, truth, gate, config)
    mode_score, low_mode_fraction, target_mode_corr = mode_signature_score(
        recovered,
        memory_basis,
        target_coeffs,
    )
    smoothness = graph_smoothness_score(adjacency, recovered)
    specificity = float(
        0.34 * recoverability
        + 0.24 * permanence
        + 0.22 * temporal
        + 0.17 * mode_score
        + 0.03 * smoothness
    )
    return SampleSummary(
        sample=sample.sample,
        sample_role=sample.sample_role,
        control_type=sample.control_type,
        memory_specificity_score=specificity,
        memory_recoverability=recoverability,
        permanence_score=permanence,
        temporal_step_score=temporal,
        mode_signature_score=mode_score,
        generic_smoothness_score=smoothness,
        recovered_target_correlation=corr,
        recovered_amplitude=amplitude,
        plateau_stability=plateau_stability,
        low_mode_energy_fraction=low_mode_fraction,
        target_mode_correlation=target_mode_corr,
    )


def run_probe(
    config: MemoryConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[SampleSpec], list[SampleSummary], dict[str, Any]]:
    points, time, gate, truth, target_coeffs, samples = build_samples(config)
    memory_basis = np.column_stack([harmonic_subspaces(points, config.l_max)[2], harmonic_subspaces(points, config.l_max)[3]])
    adjacency = knn_graph(points, config.k_neighbors)
    summaries = [
        evaluate_sample(config, sample, truth, gate, memory_basis, target_coeffs, adjacency)
        for sample in samples
    ]
    status = build_status(config, summaries)
    return points, time, gate, truth, samples, summaries, status


def build_status(config: MemoryConfig, summaries: list[SampleSummary]) -> dict[str, Any]:
    target = next(row for row in summaries if row.sample_role == "target")
    controls = [row for row in summaries if row.sample_role == "negative_control"]
    control_scores = np.array([row.memory_specificity_score for row in controls], dtype=float)
    control_threshold = target.memory_specificity_score - config.pass_min_margin
    control_false_positive_rate = float(np.mean(control_scores >= control_threshold))
    control_max = max(controls, key=lambda row: row.memory_specificity_score)
    score_margin = float(target.memory_specificity_score - control_max.memory_specificity_score)
    separation_z = float(
        (target.memory_specificity_score - float(np.mean(control_scores)))
        / max(float(np.std(control_scores)), 1.0e-9)
    )
    target_passes = bool(
        target.memory_specificity_score >= config.pass_min_memory_score
        and target.memory_recoverability >= config.pass_min_recoverability
        and target.permanence_score >= config.pass_min_permanence
        and target.temporal_step_score >= config.pass_min_temporal_step
        and target.mode_signature_score >= config.pass_min_mode_signature
    )
    controls_pass = bool(
        control_false_positive_rate <= config.pass_max_control_false_positive_rate
        and score_margin >= config.pass_min_margin
        and separation_z >= config.pass_min_separation_z
    )
    if target_passes and controls_pass:
        bridge_status = "PASS"
        failure_mode = "NONE"
    elif target.memory_specificity_score >= 0.72 and separation_z >= 1.0:
        bridge_status = "MARGINAL"
        failure_mode = "MEMORY_CONTROL_SEPARATION_INCOMPLETE"
    else:
        bridge_status = "OPEN"
        failure_mode = "TOY_MEMORY_SIGNATURE_NOT_DISTINCT"

    return {
        "phase": "61",
        "bridge_name": "gravitational_memory_toy_probe",
        "bridge_status": bridge_status,
        "failure_mode": failure_mode,
        "success_definition": (
            "PASS requires a synthetic permanent displacement-memory target on S2 to recover "
            "the known spatial mode, permanent step, temporal profile, and target harmonic "
            "coefficients while controls that preserve smoothness or time-series scale fail."
        ),
        "config": asdict(config),
        "target_memory_specificity_score": float(target.memory_specificity_score),
        "target_memory_recoverability": float(target.memory_recoverability),
        "target_permanence_score": float(target.permanence_score),
        "target_temporal_step_score": float(target.temporal_step_score),
        "target_mode_signature_score": float(target.mode_signature_score),
        "target_generic_smoothness_score": float(target.generic_smoothness_score),
        "control_max_score": float(control_max.memory_specificity_score),
        "control_max_sample": control_max.sample,
        "control_false_positive_threshold": float(control_threshold),
        "control_false_positive_rate": control_false_positive_rate,
        "score_margin_over_best_control": score_margin,
        "separation_z": separation_z,
        "core_modified": False,
        "safe_to_continue": True,
        "next_phase": "62_claim_gated_memory_update_or_real_memory_audit",
        "non_claims": [
            "no real gravitational-memory observable claim",
            "no BMS charge or supertranslation charge recovery claim",
            "no celestial holography recovery claim",
            "no soft theorem or S-matrix reconstruction claim",
            "synthetic displacement-memory toy data only",
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
        "| {sample} | {sample_role} | {control_type} | {memory_specificity_score:.6f} | {memory_recoverability:.6f} | {permanence_score:.6f} | {temporal_step_score:.6f} | {mode_signature_score:.6f} | {generic_smoothness_score:.6f} |".format(
            **asdict(row)
        )
        for row in summaries
    )
    controls = [row for row in summaries if row.sample_role == "negative_control"]
    leakage_rows = "\n".join(
        f"- `{row.sample}`: memory={row.memory_specificity_score:.6f}, "
        f"recoverability={row.memory_recoverability:.6f}, permanence={row.permanence_score:.6f}, "
        f"temporal={row.temporal_step_score:.6f}, mode={row.mode_signature_score:.6f}"
        for row in sorted(controls, key=lambda item: item.memory_specificity_score, reverse=True)
    )
    report = f"""# Phase 61 Gravitational Memory Toy Probe

## Scope

This file is generated by `experiments/physics_bridge/gravitational_memory_toy_probe/run_gravitational_memory_toy_probe.py`.

This is a synthetic displacement-memory benchmark on a sampled S2 boundary. It does not use real gravitational-wave memory data, BMS charges, celestial amplitudes, Ward identities, or S-matrix inputs.

The narrow question is:

> Can HAOS-style spectral telemetry distinguish a known permanent memory-like deformation from transient bursts, node permutations, smooth low-mode decoys, and stochastic drift controls?

## Bridge Status

- `bridge_status`: `{status['bridge_status']}`
- `failure_mode`: `{status['failure_mode']}`
- target memory score: `{status['target_memory_specificity_score']:.6f}`
- target recoverability: `{status['target_memory_recoverability']:.6f}`
- target permanence: `{status['target_permanence_score']:.6f}`
- target temporal step: `{status['target_temporal_step_score']:.6f}`
- target mode signature: `{status['target_mode_signature_score']:.6f}`
- target generic smoothness: `{status['target_generic_smoothness_score']:.6f}`
- control max sample: `{status['control_max_sample']}`
- control max score: `{status['control_max_score']:.6f}`
- control false-positive threshold: `{status['control_false_positive_threshold']:.6f}`
- control false-positive rate: `{status['control_false_positive_rate']:.6f}`
- separation z: `{status['separation_z']:.6f}`
- margin over best control: `{status['score_margin_over_best_control']:.6f}`

## Summary Table

| sample | role | control | memory_score | recoverability | permanence | temporal_step | mode_signature | smoothness |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
{summary_rows}

## Control Leakage Read

{leakage_rows}

## Interpretation

`PASS` means the toy target is separable from controls on permanent displacement recovery, event-locked temporal profile, and target harmonic-mode identity.

`MARGINAL` means the toy memory pattern is detected, but at least one control remains too competitive.

`OPEN` means the telemetry is mostly measuring smoothness, low-mode content, or generic step structure rather than the known memory signature.

## Claim Gate

This sidecar is inside the Phase 60 claim gate. Even under `PASS`, the established physics row `gravitational_memory_observable` remains `OPEN`.

Allowed language:

- toy displacement-memory proxy
- synthetic permanent S2 deformation benchmark
- memory-like step recovery under controls

Disallowed language:

- gravitational memory recovered
- BMS soft hair recovered
- real supertranslation charge detected
- celestial holography validated

## Authority

- CSV: `experiments/physics_bridge/gravitational_memory_toy_probe/memory_diagnostics.csv`
- JSON: `experiments/physics_bridge/gravitational_memory_toy_probe/bridge_status.json`
- Figures: `experiments/physics_bridge/gravitational_memory_toy_probe/figures/`
- Generator: `experiments/physics_bridge/gravitational_memory_toy_probe/run_gravitational_memory_toy_probe.py`
"""
    path.write_text(report, encoding="utf-8")


def write_plots(
    output_dir: Path,
    config: MemoryConfig,
    points: np.ndarray,
    time: np.ndarray,
    truth: np.ndarray,
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

    fig, ax = plt.subplots(figsize=(11, 5))
    colors = ["tab:green" if row.sample_role == "target" else "tab:orange" for row in summaries]
    ax.bar(labels, [row.memory_specificity_score for row in summaries], color=colors)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("memory specificity score")
    ax.set_title("Phase 61 toy memory target vs controls")
    ax.tick_params(axis="x", rotation=25, labelsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "control_comparison.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 5))
    width = 0.16
    ax.bar(x - 2.0 * width, [row.memory_recoverability for row in summaries], width, label="recoverability")
    ax.bar(x - 1.0 * width, [row.permanence_score for row in summaries], width, label="permanence")
    ax.bar(x, [row.temporal_step_score for row in summaries], width, label="temporal")
    ax.bar(x + 1.0 * width, [row.mode_signature_score for row in summaries], width, label="mode")
    ax.bar(x + 2.0 * width, [row.generic_smoothness_score for row in summaries], width, label="smoothness")
    ax.set_xticks(x, labels)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("score")
    ax.set_title("Toy memory metrics vs generic smoothness")
    ax.tick_params(axis="x", rotation=25, labelsize=8)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "metric_breakdown.png", dpi=180)
    plt.close(fig)

    best_control = max(
        [row for row in summaries if row.sample_role == "negative_control"],
        key=lambda row: row.memory_specificity_score,
    )
    target_sample = next(sample for sample in samples if sample.sample_role == "target")
    control_sample = next(sample for sample in samples if sample.sample == best_control.sample)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(time, projection_trace(target_sample.series, truth), label="target projection", color="tab:green")
    ax.plot(time, projection_trace(control_sample.series, truth), label=best_control.sample, color="tab:orange")
    ax.set_xlabel("normalized time")
    ax.set_ylabel("projection onto target memory mode")
    ax.set_title("Event-locked permanent memory trace")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "memory_projection_trace.png", dpi=180)
    plt.close(fig)

    lon = np.arctan2(points[:, 1], points[:, 0])
    lat = np.arcsin(np.clip(points[:, 2], -1.0, 1.0))
    target_recovered = recovered_memory(target_sample.series, config)
    control_recovered = recovered_memory(control_sample.series, config)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharex=True, sharey=True)
    plots = [
        ("true memory", truth),
        ("target recovered", target_recovered),
        (f"{best_control.sample} recovered", control_recovered),
    ]
    vmax = max(float(np.max(np.abs(values))) for _, values in plots)
    for ax, (title, values) in zip(axes, plots):
        sc = ax.scatter(lon, lat, c=values, cmap="coolwarm", vmin=-vmax, vmax=vmax, s=20)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("longitude")
    axes[0].set_ylabel("latitude")
    fig.colorbar(sc, ax=axes.ravel().tolist(), shrink=0.85, label="displacement")
    fig.suptitle("Toy memory field identity")
    fig.savefig(figure_dir / "memory_field_identity.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT))


def main() -> None:
    args = parse_args()
    config = make_config(args)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    points, time, _gate, truth, samples, summaries, status = run_probe(config)
    write_csv(output_dir / CSV_PATH.name, summaries)
    write_status(output_dir / STATUS_PATH.name, status)
    write_report(output_dir / REPORT_PATH.name, summaries, status)
    write_plots(output_dir, config, points, time, truth, samples, summaries)

    print("Phase 61 gravitational-memory toy probe generated.")
    print(f"CSV: {relative(output_dir / CSV_PATH.name)}")
    print(f"JSON: {relative(output_dir / STATUS_PATH.name)}")
    print(f"Report: {relative(output_dir / REPORT_PATH.name)}")
    print(f"Figures: {relative(output_dir / 'figures')}")
    print(
        f"bridge_status: {status['bridge_status']} ({status['failure_mode']}); "
        f"target={status['target_memory_specificity_score']:.6f}, "
        f"control_max={status['control_max_score']:.6f}"
    )


if __name__ == "__main__":
    main()
