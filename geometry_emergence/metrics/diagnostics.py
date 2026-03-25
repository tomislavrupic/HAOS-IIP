from __future__ import annotations

from typing import Optional

import numpy as np

from geometry_emergence.operators import InteractionGraph, TransportOperator


EPSILON = 1.0e-12


def _pairwise_feature_distances(features: np.ndarray) -> np.ndarray:
    delta = features[:, None, :] - features[None, :, :]
    return np.linalg.norm(delta, axis=-1)


def _mean(values: list[float]) -> float:
    if not values:
        return 0.0
    return float(sum(values) / len(values))


def compute_transport_arrival_times(
    transport_operator: TransportOperator,
    max_steps: int = 12,
    threshold: float = 0.01,
) -> np.ndarray:
    transition = np.asarray(transport_operator.transition, dtype=float)
    n_nodes = int(transition.shape[0])
    current = np.eye(n_nodes, dtype=float)
    arrival = np.full((n_nodes, n_nodes), np.inf, dtype=float)
    np.fill_diagonal(arrival, 0.0)

    for step in range(1, int(max_steps) + 1):
        current = current @ transition
        reached = (current >= float(threshold)) & np.isinf(arrival)
        arrival[reached] = float(step)
    return arrival


def compute_path_distortion(
    interaction_graph: InteractionGraph,
    transport_operator: TransportOperator,
    max_steps: int = 12,
    threshold: float = 0.01,
) -> dict[str, object]:
    graph_distances = interaction_graph.shortest_path_distances()
    transport_distances = compute_transport_arrival_times(
        transport_operator=transport_operator,
        max_steps=max_steps,
        threshold=threshold,
    )

    pair_mask = ~np.eye(interaction_graph.n_nodes, dtype=bool)
    graph_reachable = np.isfinite(graph_distances) & pair_mask
    transport_reachable = np.isfinite(transport_distances) & pair_mask
    jointly_reachable = graph_reachable & transport_reachable

    graph_pair_count = int(np.sum(graph_reachable))
    joint_pair_count = int(np.sum(jointly_reachable))
    reachable_fraction = (
        float(joint_pair_count / graph_pair_count)
        if graph_pair_count > 0
        else 0.0
    )

    if joint_pair_count == 0:
        return {
            "distortion_score": 1.0,
            "reachable_fraction": 0.0,
            "graph_reachable_pairs": graph_pair_count,
            "joint_reachable_pairs": 0,
        }

    graph_values = graph_distances[jointly_reachable]
    transport_values = transport_distances[jointly_reachable]
    graph_scale = max(float(np.max(graph_values)), EPSILON)
    transport_scale = max(float(np.max(transport_values)), EPSILON)
    graph_values = graph_values / graph_scale
    transport_values = transport_values / transport_scale

    relative_gap = np.abs(graph_values - transport_values) / (graph_values + transport_values + EPSILON)
    distortion_score = float(np.clip(relative_gap.mean() + 0.5 * (1.0 - reachable_fraction), 0.0, 1.0))

    return {
        "distortion_score": distortion_score,
        "reachable_fraction": reachable_fraction,
        "graph_reachable_pairs": graph_pair_count,
        "joint_reachable_pairs": joint_pair_count,
    }


def _top_k_indices(values: np.ndarray, k: int, exclude_index: int) -> np.ndarray:
    ranking = np.lexsort((np.arange(values.size), -values))
    filtered = ranking[ranking != int(exclude_index)]
    return filtered[: int(k)]


def compute_neighborhood_persistence(
    transport_operator: TransportOperator,
    k_nearest: int = 6,
    steps: int = 12,
    support_threshold: float = 0.01,
) -> dict[str, object]:
    transition = np.asarray(transport_operator.transition, dtype=float)
    n_nodes = int(transition.shape[0])
    current = np.eye(n_nodes, dtype=float)
    prior_neighbors: Optional[list[np.ndarray]] = None
    prior_coverages: Optional[list[float]] = None
    step_retention: list[float] = []
    step_support_coverage: list[float] = []

    for _ in range(int(steps)):
        current = current @ transition
        current_coverages: list[float] = []
        current_neighbors = [
            _top_k_indices(current[row], min(int(k_nearest), n_nodes - 1), row)
            for row in range(n_nodes)
        ]
        for row in range(n_nodes):
            support_mask = current[row] >= float(support_threshold)
            support_mask[row] = False
            support_count = int(np.sum(support_mask))
            denominator = max(n_nodes - 1, 1)
            current_coverages.append(float(min(support_count, denominator) / denominator))
        step_support_coverage.append(_mean(current_coverages))
        if prior_neighbors is not None:
            retention_values: list[float] = []
            for previous, updated, previous_coverage, current_coverage in zip(
                prior_neighbors,
                current_neighbors,
                prior_coverages or [],
                current_coverages,
            ):
                overlap = np.intersect1d(previous, updated, assume_unique=False).size
                retention = float(overlap / max(previous.size, 1))
                retention_values.append(retention * min(previous_coverage, current_coverage))
            step_retention.append(_mean(retention_values))
        prior_neighbors = current_neighbors
        prior_coverages = current_coverages

    return {
        "neighborhood_retention": _mean(step_retention),
        "step_retention": step_retention,
        "support_coverage": _mean(step_support_coverage),
    }


def _transport_signature_matrix(
    transport_operator: TransportOperator,
    horizon: int,
) -> np.ndarray:
    transition = np.asarray(transport_operator.transition, dtype=float)
    n_nodes = int(transition.shape[0])
    current = np.eye(n_nodes, dtype=float)
    blocks = []
    for _ in range(int(horizon)):
        current = current @ transition
        blocks.append(current.copy())
    signatures = np.hstack(blocks)
    row_norms = np.linalg.norm(signatures, axis=1, keepdims=True)
    return signatures / np.maximum(row_norms, EPSILON)


def _local_pca_dimension(
    signatures: np.ndarray,
    local_k: int,
    variance_capture: float = 0.9,
) -> tuple[float, list[float], np.ndarray]:
    pairwise_distances = _pairwise_feature_distances(signatures)
    local_dimensions: list[float] = []

    for index in range(signatures.shape[0]):
        ordering = np.argsort(pairwise_distances[index], kind="mergesort")
        neighborhood = ordering[: min(local_k + 1, signatures.shape[0])]
        centered = signatures[neighborhood] - signatures[neighborhood].mean(axis=0, keepdims=True)
        singular_values = np.linalg.svd(centered, full_matrices=False, compute_uv=False)
        energy = singular_values ** 2
        total_energy = float(np.sum(energy))
        if total_energy <= EPSILON:
            local_dimensions.append(0.0)
            continue
        cumulative_energy = np.cumsum(energy) / total_energy
        local_dimensions.append(float(np.searchsorted(cumulative_energy, variance_capture) + 1))

    return _mean(local_dimensions), local_dimensions, pairwise_distances


def _correlation_dimension(signatures: np.ndarray) -> float:
    pairwise_distances = _pairwise_feature_distances(signatures)
    upper = pairwise_distances[np.triu_indices(signatures.shape[0], k=1)]
    upper = upper[upper > EPSILON]
    if upper.size < 8:
        return 0.0

    lower = float(np.quantile(upper, 0.2))
    upper_q = float(np.quantile(upper, 0.8))
    if upper_q <= lower + EPSILON:
        return 0.0

    radii = np.linspace(lower, upper_q, 8)
    counts = np.array([(upper <= radius).mean() for radius in radii], dtype=float)
    valid = (counts > 0.0) & (counts < 1.0)
    if int(np.sum(valid)) < 3:
        return 0.0

    slope, _ = np.polyfit(np.log(radii[valid]), np.log(counts[valid]), 1)
    return float(np.clip(slope, 0.0, signatures.shape[0] - 1))


def compute_effective_dimension(
    transport_operator: TransportOperator,
    signature_horizon: int = 6,
    local_k: int = 8,
) -> dict[str, object]:
    signatures = _transport_signature_matrix(transport_operator, horizon=signature_horizon)
    local_pca_dimension, local_dimensions, _ = _local_pca_dimension(
        signatures=signatures,
        local_k=local_k,
    )
    correlation_dimension = _correlation_dimension(signatures)
    effective_dimension = float((local_pca_dimension + correlation_dimension) / 2.0)

    return {
        "local_pca_dimension": local_pca_dimension,
        "correlation_dimension": correlation_dimension,
        "effective_dimension": effective_dimension,
        "local_dimension_trace": local_dimensions,
    }


def compute_flow_bending_index(
    transport_operator: TransportOperator,
    steps: int = 12,
) -> dict[str, object]:
    transition = np.asarray(transport_operator.transition, dtype=float)
    n_nodes = int(transition.shape[0])
    current = np.eye(n_nodes, dtype=float)
    step_scores: list[float] = []

    for _ in range(int(steps)):
        updated = current @ transition
        divergence = updated - current
        smoothed_divergence = divergence @ transition
        numerator = np.mean(np.abs(divergence - smoothed_divergence), axis=1)
        denominator = np.mean(np.abs(divergence), axis=1) + EPSILON
        source_scores = numerator / denominator
        raw_score = float(np.mean(source_scores))
        step_scores.append(float(raw_score / (1.0 + raw_score)))
        current = updated

    return {
        "flow_bending_index": _mean(step_scores),
        "step_flow_bending": step_scores,
    }


def compute_recoverability(
    transport_operator: TransportOperator,
    reachable_fraction: float,
    seed: int = 17,
    burn_in_steps: int = 6,
    recovery_steps: int = 10,
    perturbation_scale: float = 0.15,
) -> dict[str, object]:
    transition = np.asarray(transport_operator.transition, dtype=float)
    n_nodes = int(transition.shape[0])
    rng = np.random.default_rng(seed)

    stationary = np.full(n_nodes, 1.0 / max(n_nodes, 1), dtype=float)
    for _ in range(max(int(burn_in_steps) + int(recovery_steps), 32)):
        stationary = stationary @ transition

    selected_nodes = rng.choice(n_nodes, size=min(4, n_nodes), replace=False)
    perturbed = np.zeros(n_nodes, dtype=float)
    perturbed[selected_nodes] = float(perturbation_scale)
    perturbed = perturbed / np.maximum(np.sum(perturbed), EPSILON)

    initial_gap = float(np.mean(np.abs(perturbed - stationary)))
    if initial_gap <= EPSILON:
        return {
            "recoverability_score": 0.0,
            "gap_decay_score": 0.0,
            "recoverability_trace": [],
        }

    gap_trace: list[float] = []
    perturbed_state = perturbed.copy()

    for _ in range(int(recovery_steps)):
        perturbed_state = perturbed_state @ transition
        gap_ratio = float(np.mean(np.abs(perturbed_state - stationary)) / initial_gap)
        gap_trace.append(gap_ratio)

    final_gap = gap_trace[-1] if gap_trace else 1.0
    gap_decay_score = float(np.clip(1.0 - final_gap, 0.0, 1.0))
    recoverability_score = float(np.sqrt(max(reachable_fraction, 0.0) * gap_decay_score))

    return {
        "recoverability_score": recoverability_score,
        "gap_decay_score": gap_decay_score,
        "recoverability_trace": gap_trace,
    }
