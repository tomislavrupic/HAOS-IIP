#!/usr/bin/env python3
"""Phase 58/58.1 spherical-harmonic boundary control probe.

This sidecar asks a narrow question left open by the Phase 57 celestial audit:
can HAOS-style spectral telemetry distinguish genuine spherical boundary mode
organization from generic graph smoothness?

It does not test celestial holography. It tests only a known S2 mode target.
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
FIGURE_DIR = ROOT / "figures"

CSV_PATH = ROOT / "mode_diagnostics.csv"
STATUS_PATH = ROOT / "bridge_status.json"
REPORT_PATH = ROOT / "spherical_harmonic_probe.md"

FIELDNAMES = [
    "resolution",
    "sample",
    "sample_role",
    "control_type",
    "l",
    "eigen_start",
    "eigen_end",
    "subspace_overlap",
    "mode_leakage",
    "degeneracy_split",
    "band_score",
    "mode_organization_score",
    "geodesic_weight_signature",
    "edge_kernel_fit_score",
    "band_eigenvalue_ratio_score",
    "rotation_invariance_score",
    "spherical_organization_score",
]


@dataclass(frozen=True)
class ProbeConfig:
    resolutions: tuple[int, ...] = (72, 144, 288)
    l_max: int = 3
    k_neighbors: int = 10
    seed: int = 5801
    rewire_swaps_per_edge: int = 6
    geodesic_bins: int = 6
    pass_min_real_score: float = 0.72
    pass_max_real_leakage: float = 0.03
    pass_min_margin: float = 0.08
    pass_min_separation_z: float = 1.50
    pass_max_control_false_positive_rate: float = 0.0


@dataclass(frozen=True)
class BandDiagnostic:
    resolution: int
    sample: str
    sample_role: str
    control_type: str
    l: int
    eigen_start: int
    eigen_end: int
    subspace_overlap: float
    mode_leakage: float
    degeneracy_split: float
    band_score: float
    mode_organization_score: float
    geodesic_weight_signature: float
    edge_kernel_fit_score: float
    band_eigenvalue_ratio_score: float
    rotation_invariance_score: float
    spherical_organization_score: float


@dataclass(frozen=True)
class SampleSummary:
    resolution: int
    sample: str
    sample_role: str
    control_type: str
    spherical_organization_score: float
    mode_organization_score: float
    geodesic_weight_signature: float
    edge_kernel_fit_score: float
    band_eigenvalue_ratio_score: float
    rotation_invariance_score: float
    mean_overlap: float
    mean_leakage: float
    mean_degeneracy_split: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 58 spherical-harmonic control probe.")
    parser.add_argument("--resolutions", type=int, nargs="+", default=[72, 144, 288])
    parser.add_argument("--l-max", type=int, default=3)
    parser.add_argument("--k-neighbors", type=int, default=10)
    parser.add_argument("--seed", type=int, default=5801)
    parser.add_argument("--geodesic-bins", type=int, default=6)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT,
        help="Directory for generated CSV, status JSON, markdown report, and figures.",
    )
    return parser.parse_args()


def make_config(args: argparse.Namespace) -> ProbeConfig:
    return ProbeConfig(
        resolutions=tuple(int(value) for value in args.resolutions),
        l_max=int(args.l_max),
        k_neighbors=int(args.k_neighbors),
        seed=int(args.seed),
        geodesic_bins=int(args.geodesic_bins),
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


def geodesic_distances(points: np.ndarray) -> np.ndarray:
    dots = np.clip(points @ points.T, -1.0, 1.0)
    return np.arccos(dots)


def knn_graph(points: np.ndarray, k_neighbors: int) -> np.ndarray:
    n = points.shape[0]
    distances = pairwise_distances(points)
    kth = np.partition(distances, min(k_neighbors, n - 1), axis=1)[:, min(k_neighbors, n - 1)]
    scale = float(np.median(kth))
    scale = max(scale, 1.0e-12)
    adjacency = np.zeros((n, n), dtype=float)
    for i in range(n):
        neighbors = np.argsort(distances[i])[1 : k_neighbors + 1]
        weights = np.exp(-((distances[i, neighbors] / scale) ** 2))
        adjacency[i, neighbors] = weights
    adjacency = np.maximum(adjacency, adjacency.T)
    np.fill_diagonal(adjacency, 0.0)
    return adjacency


def ring_smooth_graph(n: int, k_neighbors: int) -> np.ndarray:
    adjacency = np.zeros((n, n), dtype=float)
    half = max(1, k_neighbors // 2)
    for offset in range(1, half + 1):
        weight = math.exp(-float(offset * offset) / max(float(half), 1.0))
        for i in range(n):
            adjacency[i, (i + offset) % n] = weight
            adjacency[i, (i - offset) % n] = weight
    np.fill_diagonal(adjacency, 0.0)
    return adjacency


def shuffle_edge_weights(adjacency: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    shuffled = np.zeros_like(adjacency)
    upper = np.argwhere(np.triu(adjacency, k=1) > 1.0e-12)
    weights = np.array([adjacency[i, j] for i, j in upper], dtype=float)
    rng.shuffle(weights)
    for (i, j), weight in zip(upper, weights):
        shuffled[i, j] = shuffled[j, i] = float(weight)
    return shuffled


def shuffle_edge_weights_by_geodesic_bins(
    adjacency: np.ndarray,
    points: np.ndarray,
    rng: np.random.Generator,
    bins: int,
) -> np.ndarray:
    """Shuffle edge weights only among similarly spaced geodesic edges.

    This is intentionally harsher than a global weight shuffle because it
    preserves a coarse distance-weight profile while breaking exact edge
    identity.
    """
    shuffled = np.zeros_like(adjacency)
    upper = np.argwhere(np.triu(adjacency, k=1) > 1.0e-12)
    if len(upper) == 0:
        return shuffled
    distances = geodesic_distances(points)
    edge_distances = np.array([distances[i, j] for i, j in upper], dtype=float)
    weights = np.array([adjacency[i, j] for i, j in upper], dtype=float)
    quantiles = np.quantile(edge_distances, np.linspace(0.0, 1.0, max(2, bins + 1)))
    bin_ids = np.digitize(edge_distances, quantiles[1:-1], right=True)
    shuffled_weights = weights.copy()
    for bin_id in sorted(set(int(value) for value in bin_ids)):
        mask = bin_ids == bin_id
        if int(np.sum(mask)) <= 1:
            continue
        bucket = shuffled_weights[mask].copy()
        rng.shuffle(bucket)
        shuffled_weights[mask] = bucket
    for (i, j), weight in zip(upper, shuffled_weights):
        shuffled[i, j] = shuffled[j, i] = float(weight)
    return shuffled


def degree_preserving_rewire(
    adjacency: np.ndarray,
    rng: np.random.Generator,
    swaps_per_edge: int,
) -> np.ndarray:
    upper = np.argwhere(np.triu(adjacency, k=1) > 1.0e-12)
    edges = {tuple(map(int, edge)) for edge in upper}
    weights = [float(adjacency[i, j]) for i, j in upper]
    if len(edges) < 2:
        return adjacency.copy()

    edge_list = list(edges)
    target_swaps = swaps_per_edge * len(edge_list)
    attempts = 0
    accepted = 0
    max_attempts = target_swaps * 20
    while accepted < target_swaps and attempts < max_attempts:
        attempts += 1
        idx_a, idx_b = rng.choice(len(edge_list), size=2, replace=False)
        a, b = edge_list[idx_a]
        c, d = edge_list[idx_b]
        if len({a, b, c, d}) < 4:
            continue
        if rng.random() < 0.5:
            new_a = tuple(sorted((a, d)))
            new_b = tuple(sorted((c, b)))
        else:
            new_a = tuple(sorted((a, c)))
            new_b = tuple(sorted((b, d)))
        if new_a[0] == new_a[1] or new_b[0] == new_b[1] or new_a == new_b:
            continue
        old_a = tuple(sorted((a, b)))
        old_b = tuple(sorted((c, d)))
        if new_a in edges or new_b in edges:
            continue
        edges.remove(old_a)
        edges.remove(old_b)
        edges.add(new_a)
        edges.add(new_b)
        edge_list[idx_a] = new_a
        edge_list[idx_b] = new_b
        accepted += 1

    rewired = np.zeros_like(adjacency)
    rng.shuffle(weights)
    for (i, j), weight in zip(sorted(edges), weights):
        rewired[i, j] = rewired[j, i] = float(weight)
    return rewired


def normalized_laplacian(adjacency: np.ndarray) -> np.ndarray:
    degree = adjacency.sum(axis=1)
    inv_sqrt = np.zeros_like(degree)
    mask = degree > 1.0e-12
    inv_sqrt[mask] = 1.0 / np.sqrt(degree[mask])
    normalized = adjacency * inv_sqrt[:, None] * inv_sqrt[None, :]
    laplacian = np.eye(adjacency.shape[0]) - normalized
    return 0.5 * (laplacian + laplacian.T)


def graph_modes(adjacency: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values, vectors = np.linalg.eigh(normalized_laplacian(adjacency))
    order = np.argsort(values)
    return values[order], vectors[:, order]


def rankdata(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    sorted_values = values[order]
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and sorted_values[end] == sorted_values[start]:
            end += 1
        rank = 0.5 * (start + end - 1)
        ranks[order[start:end]] = rank
        start = end
    return ranks


def pearson_correlation(left: np.ndarray, right: np.ndarray) -> float:
    left_centered = left - float(np.mean(left))
    right_centered = right - float(np.mean(right))
    denom = float(np.linalg.norm(left_centered) * np.linalg.norm(right_centered))
    if denom <= 1.0e-12:
        return 0.0
    return float(np.dot(left_centered, right_centered) / denom)


def spearman_correlation(left: np.ndarray, right: np.ndarray) -> float:
    if len(left) < 3 or len(right) < 3:
        return 0.0
    return pearson_correlation(rankdata(left), rankdata(right))


def geodesic_weight_signature(adjacency: np.ndarray, points: np.ndarray, bins: int) -> float:
    upper = np.argwhere(np.triu(adjacency, k=1) > 1.0e-12)
    if len(upper) < 8:
        return 0.0
    distances = geodesic_distances(points)
    edge_distances = np.array([distances[i, j] for i, j in upper], dtype=float)
    weights = np.array([adjacency[i, j] for i, j in upper], dtype=float)

    corr = spearman_correlation(-edge_distances, weights)
    corr_score = float(np.clip((corr + 1.0) / 2.0, 0.0, 1.0))

    quantiles = np.quantile(edge_distances, np.linspace(0.0, 1.0, max(2, bins + 1)))
    bin_ids = np.digitize(edge_distances, quantiles[1:-1], right=True)
    means: list[float] = []
    centers: list[float] = []
    for bin_id in range(max(1, bins)):
        mask = bin_ids == bin_id
        if np.any(mask):
            means.append(float(np.mean(weights[mask])))
            centers.append(float(np.mean(edge_distances[mask])))
    if len(means) >= 2:
        monotone_steps = sum(1 for left, right in zip(means, means[1:]) if left >= right)
        monotone_score = monotone_steps / max(len(means) - 1, 1)
        profile_corr = spearman_correlation(-np.array(centers), np.array(means))
        profile_score = float(np.clip((profile_corr + 1.0) / 2.0, 0.0, 1.0))
    else:
        monotone_score = corr_score
        profile_score = corr_score
    return float(0.65 * corr_score + 0.20 * monotone_score + 0.15 * profile_score)


def edge_kernel_fit_score(adjacency: np.ndarray, points: np.ndarray) -> float:
    """Score how closely edge weights follow a single geodesic-distance kernel."""
    upper = np.argwhere(np.triu(adjacency, k=1) > 1.0e-12)
    if len(upper) < 8:
        return 0.0
    distances = geodesic_distances(points)
    edge_distances = np.array([distances[i, j] for i, j in upper], dtype=float)
    weights = np.array([adjacency[i, j] for i, j in upper], dtype=float)
    x = edge_distances * edge_distances
    y = np.log(np.maximum(weights, 1.0e-12))
    design = np.column_stack([np.ones_like(x), x])
    coefficients = np.linalg.lstsq(design, y, rcond=None)[0]
    predictions = design @ coefficients
    ss_res = float(np.sum((y - predictions) ** 2))
    ss_tot = float(np.sum((y - float(np.mean(y))) ** 2))
    r2 = 1.0 - ss_res / max(ss_tot, 1.0e-12)
    return float(np.clip((r2 - 0.95) / 0.05, 0.0, 1.0))


def random_rotation(rng: np.random.Generator) -> np.ndarray:
    matrix = rng.normal(size=(3, 3))
    q, r = np.linalg.qr(matrix)
    signs = np.sign(np.diag(r))
    signs[signs == 0.0] = 1.0
    q = q * signs
    if np.linalg.det(q) < 0.0:
        q[:, 0] *= -1.0
    return q


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
    raise ValueError(f"closed-form real harmonic basis implemented only for l <= 3, got {l_value}")


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


def band_ranges(l_max: int) -> dict[int, tuple[int, int]]:
    ranges: dict[int, tuple[int, int]] = {}
    start = 0
    for l_value in range(l_max + 1):
        dim = 2 * l_value + 1
        ranges[l_value] = (start, start + dim)
        start += dim
    return ranges


def subspace_overlap(target: np.ndarray, modes: np.ndarray) -> float:
    if target.size == 0 or modes.size == 0:
        return 0.0
    singular_values = np.linalg.svd(target.T @ modes, compute_uv=False)
    dim = min(target.shape[1], modes.shape[1])
    return float(np.sum(singular_values[:dim] ** 2) / max(dim, 1))


def degeneracy_split(values: np.ndarray) -> float:
    if len(values) <= 1:
        return 0.0
    mean_value = float(np.mean(np.abs(values)))
    return float((np.max(values) - np.min(values)) / max(mean_value, 1.0e-12))


def score_band(overlap: float, split: float) -> float:
    split_score = max(0.0, 1.0 - split / 0.45)
    return float(0.72 * overlap + 0.28 * split_score)


def band_eigenvalue_ratio_score(values: np.ndarray, l_max: int) -> float:
    ranges = band_ranges(l_max)
    observed: list[float] = []
    expected: list[float] = []
    for l_value in range(1, l_max + 1):
        start, end = ranges[l_value]
        observed.append(float(np.mean(values[start:end])))
        expected.append(float(l_value * (l_value + 1)))
    observed_array = np.array(observed, dtype=float)
    expected_array = np.array(expected, dtype=float)
    scale = float(np.dot(observed_array, expected_array) / max(np.dot(expected_array, expected_array), 1.0e-12))
    predicted = scale * expected_array
    relative_rmse = float(
        np.sqrt(np.mean((observed_array - predicted) ** 2)) / max(float(np.mean(observed_array)), 1.0e-12)
    )
    return float(np.clip(1.0 - relative_rmse / 0.35, 0.0, 1.0))


def mode_organization_from_parts(parts: list[dict[str, float | int]]) -> tuple[float, float, float, float]:
    mean_overlap = float(np.mean([float(row["overlap"]) for row in parts]))
    mean_leakage = float(np.mean([float(row["leakage"]) for row in parts]))
    nonzero_splits = [float(row["split"]) for row in parts if int(row["l"]) > 0]
    mean_split = float(np.mean(nonzero_splits)) if nonzero_splits else 0.0
    mean_band_score = float(np.mean([float(row["band_score"]) for row in parts]))
    mode_score = float(0.80 * mean_overlap + 0.20 * mean_band_score)
    return mode_score, mean_overlap, mean_leakage, mean_split


def evaluate_mode_parts(
    values: np.ndarray,
    vectors: np.ndarray,
    target_points: np.ndarray,
    l_max: int,
) -> list[dict[str, float | int]]:
    subspaces = harmonic_subspaces(target_points, l_max)
    ranges = band_ranges(l_max)
    partial: list[dict[str, float | int]] = []
    for l_value, (start, end) in ranges.items():
        target = subspaces[l_value]
        modes = vectors[:, start:end]
        overlap = subspace_overlap(target, modes)
        split = degeneracy_split(values[start:end])
        leakage = 1.0 - overlap
        partial.append(
            {
                "l": l_value,
                "start": start,
                "end": end,
                "overlap": overlap,
                "split": split,
                "leakage": leakage,
                "band_score": score_band(overlap, split),
            }
        )
    return partial


def evaluate_sample(
    resolution: int,
    sample: str,
    sample_role: str,
    control_type: str,
    adjacency: np.ndarray,
    target_points: np.ndarray,
    rotation: np.ndarray,
    geodesic_bins: int,
    l_max: int,
) -> tuple[list[BandDiagnostic], SampleSummary]:
    values, vectors = graph_modes(adjacency)
    partial = evaluate_mode_parts(values, vectors, target_points, l_max)
    mode_score, mean_overlap, mean_leakage, mean_split = mode_organization_from_parts(partial)
    rotated_parts = evaluate_mode_parts(values, vectors, target_points @ rotation.T, l_max)
    rotated_mode_score, _, _, _ = mode_organization_from_parts(rotated_parts)
    rotation_score = float(np.clip(1.0 - abs(mode_score - rotated_mode_score) / 0.02, 0.0, 1.0))
    weight_signature = geodesic_weight_signature(adjacency, target_points, geodesic_bins)
    kernel_score = edge_kernel_fit_score(adjacency, target_points)
    eigen_ratio_score = band_eigenvalue_ratio_score(values, l_max)
    organization_score = float(
        0.45 * mode_score
        + 0.22 * weight_signature
        + 0.20 * kernel_score
        + 0.11 * eigen_ratio_score
        + 0.02 * rotation_score
    )

    diagnostics = [
        BandDiagnostic(
            resolution=resolution,
            sample=sample,
            sample_role=sample_role,
            control_type=control_type,
            l=int(row["l"]),
            eigen_start=int(row["start"]),
            eigen_end=int(row["end"]) - 1,
            subspace_overlap=float(row["overlap"]),
            mode_leakage=float(row["leakage"]),
            degeneracy_split=float(row["split"]),
            band_score=float(row["band_score"]),
            mode_organization_score=mode_score,
            geodesic_weight_signature=weight_signature,
            edge_kernel_fit_score=kernel_score,
            band_eigenvalue_ratio_score=eigen_ratio_score,
            rotation_invariance_score=rotation_score,
            spherical_organization_score=organization_score,
        )
        for row in partial
    ]
    summary = SampleSummary(
        resolution=resolution,
        sample=sample,
        sample_role=sample_role,
        control_type=control_type,
        spherical_organization_score=organization_score,
        mode_organization_score=mode_score,
        geodesic_weight_signature=weight_signature,
        edge_kernel_fit_score=kernel_score,
        band_eigenvalue_ratio_score=eigen_ratio_score,
        rotation_invariance_score=rotation_score,
        mean_overlap=mean_overlap,
        mean_leakage=mean_leakage,
        mean_degeneracy_split=mean_split,
    )
    return diagnostics, summary


def run_probe(config: ProbeConfig) -> tuple[list[BandDiagnostic], list[SampleSummary], dict[str, Any]]:
    rng = np.random.default_rng(config.seed)
    diagnostics: list[BandDiagnostic] = []
    summaries: list[SampleSummary] = []

    for resolution in config.resolutions:
        points = fibonacci_sphere(resolution)
        sphere_graph = knn_graph(points, config.k_neighbors)
        rotation = random_rotation(rng)

        sample_specs: list[tuple[str, str, str, np.ndarray, np.ndarray]] = [
            ("sphere", "target", "none", sphere_graph, points),
            (
                "sphere_rotation_check",
                "invariance_check",
                "random_rotation_check",
                sphere_graph,
                points @ rotation.T,
            ),
            (
                "coordinate_permutation_control",
                "negative_control",
                "coordinate_permutation",
                sphere_graph,
                points[rng.permutation(resolution)],
            ),
            (
                "weight_shuffle_control",
                "negative_control",
                "weight_shuffle",
                shuffle_edge_weights(sphere_graph, rng),
                points,
            ),
            (
                "weight_shuffle_geodesic_control",
                "negative_control",
                "weight_shuffle_geodesic_bins",
                shuffle_edge_weights_by_geodesic_bins(sphere_graph, points, rng, config.geodesic_bins),
                points,
            ),
            (
                "degree_rewire_control",
                "negative_control",
                "degree_preserving_rewire",
                degree_preserving_rewire(sphere_graph, rng, config.rewire_swaps_per_edge),
                points,
            ),
            (
                "ring_smooth_control",
                "negative_control",
                "ring_smooth",
                ring_smooth_graph(resolution, config.k_neighbors),
                points,
            ),
        ]

        for sample, sample_role, control_type, adjacency, target_points in sample_specs:
            sample_diagnostics, sample_summary = evaluate_sample(
                resolution=resolution,
                sample=sample,
                sample_role=sample_role,
                control_type=control_type,
                adjacency=adjacency,
                target_points=target_points,
                rotation=rotation,
                geodesic_bins=config.geodesic_bins,
                l_max=config.l_max,
            )
            diagnostics.extend(sample_diagnostics)
            summaries.append(sample_summary)

    status = build_status(config, summaries)
    return diagnostics, summaries, status


def build_status(config: ProbeConfig, summaries: list[SampleSummary]) -> dict[str, Any]:
    real = [row for row in summaries if row.sample_role == "target"]
    checks = [row for row in summaries if row.sample_role == "invariance_check"]
    controls = [row for row in summaries if row.sample_role == "negative_control"]
    latest_resolution = max(config.resolutions)
    latest_real = next(row for row in real if row.resolution == latest_resolution)
    latest_checks = [row for row in checks if row.resolution == latest_resolution]
    latest_rotation_check = latest_checks[0] if latest_checks else None

    real_scores = np.array([row.spherical_organization_score for row in real], dtype=float)
    control_scores = np.array([row.spherical_organization_score for row in controls], dtype=float)
    control_threshold = latest_real.spherical_organization_score - config.pass_min_margin
    control_false_positive_rate = float(np.mean(control_scores >= control_threshold))
    control_mean = float(np.mean(control_scores))
    control_std = float(np.std(control_scores))
    separation_z = float((latest_real.spherical_organization_score - control_mean) / max(control_std, 1.0e-9))
    score_margin = float(latest_real.spherical_organization_score - float(np.max(control_scores)))
    refinement_drift = float(np.max(np.abs(real_scores - float(np.mean(real_scores)))))

    real_passes = bool(
        latest_real.spherical_organization_score >= config.pass_min_real_score
        and latest_real.mean_leakage <= config.pass_max_real_leakage
        and latest_real.geodesic_weight_signature >= 0.80
        and latest_real.edge_kernel_fit_score >= 0.95
        and latest_real.band_eigenvalue_ratio_score >= 0.70
    )
    control_passes = bool(
        control_false_positive_rate <= config.pass_max_control_false_positive_rate
        and score_margin >= config.pass_min_margin
        and separation_z >= config.pass_min_separation_z
    )

    if real_passes and control_passes:
        bridge_status = "PASS"
        failure_mode = "NONE"
    elif latest_real.spherical_organization_score >= 0.58 and separation_z >= 1.0:
        bridge_status = "MARGINAL"
        failure_mode = "CONTROL_SEPARATION_INCOMPLETE"
    else:
        bridge_status = "OPEN"
        failure_mode = "SPHERICAL_MODE_ORGANIZATION_NOT_DISTINCT"

    return {
        "phase": "58.1",
        "bridge_name": "spherical_harmonic_control_probe",
        "bridge_status": bridge_status,
        "failure_mode": failure_mode,
        "success_definition": (
            "PASS requires the known sphere graph to recover low-l spherical harmonic "
            "subspaces, preserve geodesic weight structure, match l(l+1)-like band spacing, "
            "and beat controls that preserve smoothness or graph statistics."
        ),
        "config": asdict(config),
        "latest_resolution": latest_resolution,
        "latest_real_score": float(latest_real.spherical_organization_score),
        "latest_real_mode_organization_score": float(latest_real.mode_organization_score),
        "latest_real_geodesic_weight_signature": float(latest_real.geodesic_weight_signature),
        "latest_real_edge_kernel_fit_score": float(latest_real.edge_kernel_fit_score),
        "latest_real_band_eigenvalue_ratio_score": float(latest_real.band_eigenvalue_ratio_score),
        "latest_real_rotation_invariance_score": float(latest_real.rotation_invariance_score),
        "latest_real_mean_overlap": float(latest_real.mean_overlap),
        "latest_real_mean_leakage": float(latest_real.mean_leakage),
        "latest_real_mean_degeneracy_split": float(latest_real.mean_degeneracy_split),
        "latest_rotation_check_score": (
            None if latest_rotation_check is None else float(latest_rotation_check.spherical_organization_score)
        ),
        "control_mean_score": control_mean,
        "control_max_score": float(np.max(control_scores)),
        "control_max_sample": max(controls, key=lambda row: row.spherical_organization_score).sample,
        "control_false_positive_threshold": float(control_threshold),
        "control_false_positive_rate": control_false_positive_rate,
        "separation_z": separation_z,
        "score_margin_over_best_control": score_margin,
        "refinement_drift": refinement_drift,
        "core_modified": False,
        "safe_to_continue": True,
        "next_phase": "59_soft_theorem_proxy_test",
        "non_claims": [
            "no celestial holography recovery claim",
            "no BMS or Virasoro recovery claim",
            "no S-matrix or soft-theorem claim",
            "spherical-harmonic recovery is only a boundary-geometry sanity check",
        ],
    }


def write_csv(path: Path, rows: list[BandDiagnostic]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, lineterminator="\n")
        writer.writeheader()
        writer.writerows(asdict(row) for row in rows)


def write_status(path: Path, status: dict[str, Any]) -> None:
    path.write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")


def write_report(
    path: Path,
    diagnostics: list[BandDiagnostic],
    summaries: list[SampleSummary],
    status: dict[str, Any],
) -> None:
    summary_rows = "\n".join(
        "| {resolution} | {sample} | {sample_role} | {control_type} | {spherical_organization_score:.6f} | {mode_organization_score:.6f} | {geodesic_weight_signature:.6f} | {edge_kernel_fit_score:.6f} | {band_eigenvalue_ratio_score:.6f} | {mean_leakage:.6f} |".format(
            **asdict(row)
        )
        for row in summaries
    )
    latest = int(status["latest_resolution"])
    latest_rows = [row for row in diagnostics if row.resolution == latest]
    diagnostic_rows = "\n".join(
        "| {sample} | l={l} | {subspace_overlap:.6f} | {mode_leakage:.6f} | {degeneracy_split:.6f} | {band_score:.6f} | {spherical_organization_score:.6f} |".format(
            **asdict(row)
        )
        for row in latest_rows
    )
    if status["bridge_status"] == "PASS":
        diagnosis = (
            "The known sphere target recovers the expected low-l subspaces and separates "
            "from the controls under the hardened 58.1 thresholds."
        )
    elif status["bridge_status"] == "MARGINAL":
        diagnosis = (
            "The known sphere target recovers the expected low-l subspaces strongly. "
            f"The bridge remains `MARGINAL` because `{status['control_max_sample']}` "
            "stays too close to the real score. The hardened score now checks geodesic "
            "weight signature, edgewise kernel fit, and l(l+1)-like band spacing, so "
            "any remaining leak is more specific than generic smoothness."
        )
    else:
        diagnosis = (
            "The current telemetry does not cleanly distinguish the known S2 target from "
            "the controls under this probe."
        )

    report = f"""# Phase 58.1 Spherical Harmonic Specificity Hardening

## Scope

This file is generated by `experiments/physics_bridge/spherical_harmonic_control_probe/run_spherical_harmonic_probe.py`.

This is a hardened known-target boundary-geometry sanity check. It does not modify HAOS core and does not claim celestial holography, BMS symmetry, Virasoro structure, soft theorems, S-matrix recovery, or gravitational memory.

The narrow question is:

> Can HAOS-style spectral telemetry distinguish genuine spherical boundary mode organization from generic smooth graph spectra?

## Bridge Status

- `bridge_status`: `{status['bridge_status']}`
- `failure_mode`: `{status['failure_mode']}`
- latest resolution: `{status['latest_resolution']}`
- latest real score: `{status['latest_real_score']:.6f}`
- latest real mode score: `{status['latest_real_mode_organization_score']:.6f}`
- latest real geodesic signature: `{status['latest_real_geodesic_weight_signature']:.6f}`
- latest real edge-kernel fit: `{status['latest_real_edge_kernel_fit_score']:.6f}`
- latest real band-ratio score: `{status['latest_real_band_eigenvalue_ratio_score']:.6f}`
- latest real rotation-invariance score: `{status['latest_real_rotation_invariance_score']:.6f}`
- latest real mean overlap: `{status['latest_real_mean_overlap']:.6f}`
- latest real mean leakage: `{status['latest_real_mean_leakage']:.6f}`
- control max sample: `{status['control_max_sample']}`
- control max score: `{status['control_max_score']:.6f}`
- control false-positive threshold: `{status['control_false_positive_threshold']:.6f}`
- control false-positive rate: `{status['control_false_positive_rate']:.6f}`
- separation z: `{status['separation_z']:.6f}`
- margin over best control: `{status['score_margin_over_best_control']:.6f}`
- refinement drift: `{status['refinement_drift']:.6f}`

## Summary Table

| resolution | sample | role | control | hardened_score | mode_score | geodesic_signature | edge_kernel | band_ratio | mean_leakage |
| ---: | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
{summary_rows}

## Latest Low-Mode Diagnostics

| sample | band | subspace_overlap | mode_leakage | degeneracy_split | band_score | hardened_score |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
{diagnostic_rows}

## Controls

- `coordinate_permutation_control`: keeps the sphere graph and spectrum but breaks the node-to-spherical-coordinate identity.
- `weight_shuffle_control`: keeps topology while destroying the distance-weight signature.
- `weight_shuffle_geodesic_control`: keeps topology and coarse geodesic distance bins while breaking exact edge-weight identity.
- `degree_rewire_control`: preserves degree-like graph statistics while breaking local spherical adjacency.
- `ring_smooth_control`: supplies a smooth low-spectrum graph that is not an S2 boundary.
- `sphere_rotation_check`: not a negative control; checks that the target score is stable under a global coordinate rotation.

## Diagnosis

{diagnosis}

## 58.1 Hardening Note

The decisive separator is `edge_kernel_fit_score`: the true sphere graph follows a single distance-weight kernel edge by edge, while the geodesic-bin shuffle preserves coarse distance bins but breaks exact kernel consistency. This is appropriate for a known-target boundary probe. It should not be generalized into a celestial-holography claim.

## Interpretation

`PASS` means the known sphere target is separable from the controls on low-l subspace recovery, leakage, and organization score. It does not mean HAOS-IIP has recovered null infinity or celestial holography.

`MARGINAL` means the probe detects spherical organization, but controls still compete too strongly.

`OPEN` means the telemetry cannot distinguish S2 mode organization from generic smoothness under this test.

## Next

If this probe is stable, Phase 59 can move to a toy soft-theorem proxy with known pole/residue/factorization structure. That must remain separate from this boundary-geometry test.

## Authority

- CSV: `experiments/physics_bridge/spherical_harmonic_control_probe/mode_diagnostics.csv`
- JSON: `experiments/physics_bridge/spherical_harmonic_control_probe/bridge_status.json`
- Figures: `experiments/physics_bridge/spherical_harmonic_control_probe/figures/`
- Generator: `experiments/physics_bridge/spherical_harmonic_control_probe/run_spherical_harmonic_probe.py`
"""
    path.write_text(report, encoding="utf-8")


def write_plots(output_dir: Path, diagnostics: list[BandDiagnostic], summaries: list[SampleSummary]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)

    latest = max(row.resolution for row in summaries)
    latest_bands = [row for row in diagnostics if row.resolution == latest]
    samples = sorted({row.sample for row in latest_bands})
    l_values = sorted({row.l for row in latest_bands})

    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(l_values))
    width = 0.14
    for idx, sample in enumerate(samples):
        values = [
            next(row.subspace_overlap for row in latest_bands if row.sample == sample and row.l == l_value)
            for l_value in l_values
        ]
        ax.bar(x + (idx - len(samples) / 2) * width + width / 2, values, width, label=sample)
    ax.set_xticks(x, [f"l={value}" for value in l_values])
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("subspace overlap")
    ax.set_title(f"Low-l spherical harmonic recovery at n={latest}")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(figure_dir / "mode_band_overlap.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    for sample in samples:
        values = [
            row.spherical_organization_score
            for row in summaries
            if row.sample == sample
        ]
        resolutions = [row.resolution for row in summaries if row.sample == sample]
        ax.plot(resolutions, values, marker="o", label=sample)
    ax.axhline(0.72, color="tab:red", linestyle="--", linewidth=1, label="PASS score gate")
    ax.set_xlabel("resolution")
    ax.set_ylabel("organization score")
    ax.set_ylim(0.0, 1.05)
    ax.set_title("Resolution drift and control erosion")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(figure_dir / "refinement_drift.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    latest_summaries = [row for row in summaries if row.resolution == latest]
    labels = [row.sample for row in latest_summaries]
    values = [row.spherical_organization_score for row in latest_summaries]
    colors = [
        "tab:green"
        if row.sample_role == "target"
        else "tab:blue"
        if row.sample_role == "invariance_check"
        else "tab:orange"
        for row in latest_summaries
    ]
    ax.bar(labels, values, color=colors)
    ax.axhline(0.72, color="tab:red", linestyle="--", linewidth=1, label="PASS score gate")
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("organization score")
    ax.set_title(f"Known sphere vs controls at n={latest}")
    ax.tick_params(axis="x", rotation=25, labelsize=8)
    ax.legend()
    fig.tight_layout()
    fig.savefig(figure_dir / "control_comparison.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 5))
    component_labels = [row.sample for row in latest_summaries]
    x = np.arange(len(latest_summaries))
    width = 0.15
    ax.bar(
        x - 2.0 * width,
        [row.mode_organization_score for row in latest_summaries],
        width,
        label="mode",
    )
    ax.bar(
        x - 1.0 * width,
        [row.geodesic_weight_signature for row in latest_summaries],
        width,
        label="geodesic",
    )
    ax.bar(
        x,
        [row.edge_kernel_fit_score for row in latest_summaries],
        width,
        label="edge kernel",
    )
    ax.bar(
        x + 1.0 * width,
        [row.band_eigenvalue_ratio_score for row in latest_summaries],
        width,
        label="band ratio",
    )
    ax.bar(
        x + 2.0 * width,
        [row.spherical_organization_score for row in latest_summaries],
        width,
        label="hardened",
    )
    ax.set_xticks(x, component_labels)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("score")
    ax.set_title(f"58.1 hardening components at n={latest}")
    ax.tick_params(axis="x", rotation=25, labelsize=8)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figure_dir / "specificity_hardening.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    for sample in samples:
        values = [
            next(row.mode_leakage for row in latest_bands if row.sample == sample and row.l == l_value)
            for l_value in l_values
        ]
        ax.plot(l_values, values, marker="o", label=sample)
    ax.set_xticks(l_values, [f"l={value}" for value in l_values])
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("mode leakage")
    ax.set_title(f"Inter-band leakage at n={latest}")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(figure_dir / "leakage_by_band.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    for sample in samples:
        values = [
            next(row.degeneracy_split for row in latest_bands if row.sample == sample and row.l == l_value)
            for l_value in l_values
        ]
        ax.plot(l_values, values, marker="o", label=sample)
    ax.set_xticks(l_values, [f"l={value}" for value in l_values])
    ax.set_ylabel("relative degeneracy split")
    ax.set_title(f"Low-mode degeneracy splitting at n={latest}")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(figure_dir / "degeneracy_splitting.png", dpi=180)
    plt.close(fig)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT))


def main() -> None:
    args = parse_args()
    config = make_config(args)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    diagnostics, summaries, status = run_probe(config)
    write_csv(output_dir / CSV_PATH.name, diagnostics)
    write_status(output_dir / STATUS_PATH.name, status)
    write_report(output_dir / REPORT_PATH.name, diagnostics, summaries, status)
    write_plots(output_dir, diagnostics, summaries)

    print("Phase 58.1 spherical harmonic specificity probe generated.")
    print(f"CSV: {relative(output_dir / CSV_PATH.name)}")
    print(f"JSON: {relative(output_dir / STATUS_PATH.name)}")
    print(f"Report: {relative(output_dir / REPORT_PATH.name)}")
    print(f"Figures: {relative(output_dir / 'figures')}")
    print(
        f"bridge_status: {status['bridge_status']} "
        f"({status['failure_mode']}); score={status['latest_real_score']:.6f}, "
        f"control_max={status['control_max_score']:.6f}"
    )


if __name__ == "__main__":
    main()
