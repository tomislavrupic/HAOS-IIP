from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class GraphMetrics:
    algebraic_connectivity: float
    effective_conductance: float


def ring_with_backups(
    node_count: int,
    backup_fraction: float,
    backup_modifiers: np.ndarray | None = None,
) -> np.ndarray:
    if node_count < 8:
        raise ValueError("node_count must be at least 8")
    if not 0.0 <= backup_fraction < 1.0:
        raise ValueError("backup_fraction must be in [0, 1)")
    modifiers = np.ones(node_count, dtype=float) if backup_modifiers is None else np.asarray(backup_modifiers, dtype=float)
    if modifiers.shape != (node_count,) or np.any(modifiers < 0.0) or float(modifiers.sum()) <= 0.0:
        raise ValueError("backup_modifiers must be nonnegative and match node_count")
    modifiers = modifiers * (node_count / float(modifiers.sum()))

    weights = np.zeros((node_count, node_count), dtype=float)
    primary_weight = 1.0 - backup_fraction
    for index in range(node_count):
        neighbor = (index + 1) % node_count
        weights[index, neighbor] = primary_weight
        weights[neighbor, index] = primary_weight
    for index, modifier in enumerate(modifiers):
        neighbor = (index + 2) % node_count
        backup_weight = backup_fraction * float(modifier)
        weights[index, neighbor] = backup_weight
        weights[neighbor, index] = backup_weight
    return weights


def interaction_budget(weights: np.ndarray) -> float:
    return float(np.triu(weights, k=1).sum())


def graph_metrics(weights: np.ndarray) -> GraphMetrics:
    laplacian = np.diag(weights.sum(axis=1)) - weights
    eigenvalues = np.linalg.eigvalsh(laplacian)
    algebraic_connectivity = float(max(eigenvalues[1], 0.0))
    unseen = set(range(weights.shape[0]))
    components: list[list[int]] = []
    while unseen:
        start = unseen.pop()
        component = [start]
        frontier = [start]
        while frontier:
            node = frontier.pop()
            neighbors = [int(index) for index in np.flatnonzero(weights[node] > 1.0e-15)]
            for neighbor in neighbors:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    component.append(neighbor)
                    frontier.append(neighbor)
        components.append(sorted(component))

    component_by_node: dict[int, int] = {}
    component_pseudoinverses: dict[int, tuple[list[int], np.ndarray]] = {}
    for component_index, component in enumerate(components):
        for node in component:
            component_by_node[node] = component_index
        if len(component) > 1:
            subgraph = weights[np.ix_(component, component)]
            sublaplacian = np.diag(subgraph.sum(axis=1)) - subgraph
            component_pseudoinverses[component_index] = (component, np.linalg.pinv(sublaplacian))

    inverse_resistances: list[float] = []
    for left in range(weights.shape[0]):
        for right in range(left + 1, weights.shape[0]):
            component_index = component_by_node[left]
            if component_index != component_by_node[right]:
                inverse_resistances.append(0.0)
                continue
            component, pseudoinverse = component_pseudoinverses[component_index]
            local_left = component.index(left)
            local_right = component.index(right)
            resistance = float(
                pseudoinverse[local_left, local_left]
                + pseudoinverse[local_right, local_right]
                - 2.0 * pseudoinverse[local_left, local_right]
            )
            if resistance <= 0.0:
                raise ValueError("within-component effective resistance must be positive")
            inverse_resistances.append(1.0 / resistance)
    return GraphMetrics(
        algebraic_connectivity=algebraic_connectivity,
        effective_conductance=float(np.mean(inverse_resistances)),
    )


def recovery_capacity(metrics: GraphMetrics, reference: GraphMetrics) -> float:
    connectivity_ratio = metrics.algebraic_connectivity / reference.algebraic_connectivity
    conductance_ratio = metrics.effective_conductance / reference.effective_conductance
    if connectivity_ratio + conductance_ratio <= 0.0:
        return 0.0
    return float(2.0 * connectivity_ratio * conductance_ratio / (connectivity_ratio + conductance_ratio))


def operator_identity(left: np.ndarray, right: np.ndarray) -> float:
    numerator = float(np.sum(left * right))
    denominator = float(np.sqrt(np.sum(left * left) * np.sum(right * right)))
    if denominator <= 0.0:
        raise ValueError("operator identity is undefined for a zero graph")
    return numerator / denominator


def primary_capacity_retention(backup_fraction: float, baseline_backup_fraction: float) -> float:
    return (1.0 - backup_fraction) / (1.0 - baseline_backup_fraction)


def adaptive_backup_fraction(
    baseline_backup_fraction: float,
    observed_exposure_severity: float,
    adaptation_gain: float,
    cap: float,
) -> float:
    if not 0.0 <= observed_exposure_severity <= 1.0:
        raise ValueError("observed_exposure_severity must be in [0, 1]")
    return min(cap, baseline_backup_fraction + adaptation_gain * observed_exposure_severity)


def random_backup_modifiers(node_count: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.lognormal(mean=0.0, sigma=1.0, size=node_count)


def stress_edge_indices(node_count: int, family: str, start: int, block_size: int) -> tuple[int, list[int]]:
    if family in {"contiguous_primary_damage", "displaced_primary_block"}:
        return 1, [int((start + offset) % node_count) for offset in range(block_size)]
    if family == "distributed_primary_damage":
        stride = max(node_count // block_size, 1)
        return 1, [int((start + offset * stride) % node_count) for offset in range(block_size)]
    if family == "backup_channel_damage":
        return 2, [int((start + offset) % node_count) for offset in range(block_size)]
    raise ValueError(f"unknown stress family: {family}")


def apply_stress(
    weights: np.ndarray,
    family: str,
    start: int,
    block_size: int,
    severity: float,
) -> np.ndarray:
    if not 0.0 <= severity <= 1.0:
        raise ValueError("severity must be in [0, 1]")
    stressed = weights.copy()
    distance, indices = stress_edge_indices(weights.shape[0], family, start, block_size)
    for index in indices:
        neighbor = (index + distance) % weights.shape[0]
        stressed[index, neighbor] *= 1.0 - severity
        stressed[neighbor, index] *= 1.0 - severity
    return stressed


def disjoint_primary_locations(
    node_count: int,
    exposure_start: int,
    test_start: int,
    block_size: int,
    test_family: str,
) -> bool:
    if test_family == "backup_channel_damage":
        return True
    _, exposure = stress_edge_indices(node_count, "contiguous_primary_damage", exposure_start, block_size)
    _, held_out = stress_edge_indices(node_count, test_family, test_start, block_size)
    return set(exposure).isdisjoint(held_out)
