from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np


EPSILON = 1.0e-12


def _pairwise_distances(positions: np.ndarray) -> np.ndarray:
    delta = positions[:, None, :] - positions[None, :, :]
    return np.linalg.norm(delta, axis=-1)


@dataclass(frozen=True)
class KernelField:
    kernel_width: float
    locality_radius: float

    def build_affinity(self, positions: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if self.kernel_width <= 0.0:
            raise ValueError("kernel_width must be positive")
        if self.locality_radius <= 0.0:
            raise ValueError("locality_radius must be positive")

        distances = _pairwise_distances(positions)
        affinity = np.exp(-(distances ** 2) / (2.0 * self.kernel_width ** 2))
        affinity[distances > self.locality_radius] = 0.0
        np.fill_diagonal(affinity, 0.0)
        return affinity, distances


@dataclass(frozen=True)
class InteractionGraph:
    positions: np.ndarray
    affinity: np.ndarray
    distances: np.ndarray
    kernel_width: float
    locality_radius: float

    @property
    def n_nodes(self) -> int:
        return int(self.positions.shape[0])

    def degrees(self) -> np.ndarray:
        return np.sum(self.affinity > 0.0, axis=1).astype(int)

    def mean_degree(self) -> float:
        return float(np.mean(self.degrees()))

    def edge_density(self) -> float:
        n_nodes = self.n_nodes
        if n_nodes <= 1:
            return 0.0
        upper = np.triu(self.affinity > 0.0, k=1)
        total_pairs = n_nodes * (n_nodes - 1) / 2.0
        return float(np.sum(upper) / total_pairs)

    def connected_components(self) -> List[List[int]]:
        adjacency = self.affinity > 0.0
        seen = np.zeros(self.n_nodes, dtype=bool)
        components: List[List[int]] = []

        for start in range(self.n_nodes):
            if seen[start]:
                continue
            stack = [start]
            seen[start] = True
            component: List[int] = []
            while stack:
                node = stack.pop()
                component.append(node)
                neighbors = np.flatnonzero(adjacency[node])
                for neighbor in neighbors.tolist():
                    if not seen[neighbor]:
                        seen[neighbor] = True
                        stack.append(int(neighbor))
            components.append(sorted(component))
        return components

    def connected_component_count(self) -> int:
        return len(self.connected_components())

    def edge_costs(self) -> np.ndarray:
        costs = np.full_like(self.affinity, np.inf, dtype=float)
        mask = self.affinity > 0.0
        costs[mask] = -np.log(np.clip(self.affinity[mask], EPSILON, 1.0))
        np.fill_diagonal(costs, 0.0)
        return costs

    def shortest_path_distances(self) -> np.ndarray:
        distances = self.edge_costs().copy()
        for pivot in range(self.n_nodes):
            distances = np.minimum(
                distances,
                distances[:, [pivot]] + distances[[pivot], :],
            )
        return distances


@dataclass(frozen=True)
class TransportOperator:
    transition: np.ndarray

    def apply_transport_step(self, state: np.ndarray) -> np.ndarray:
        array = np.asarray(state, dtype=float)
        if array.ndim == 1:
            return self.transition.T @ array
        if array.ndim == 2:
            return array @ self.transition
        raise ValueError("state must be one-dimensional or two-dimensional")


def _resolve_locality_radius(
    kernel_width: float,
    locality_radius: Optional[float],
    embedding_dim: int,
) -> float:
    if locality_radius is not None:
        return float(locality_radius)
    return float(min(3.0 * kernel_width, np.sqrt(float(embedding_dim))))


def _build_graph_from_positions(
    positions: np.ndarray,
    kernel_width: float,
    locality_radius: Optional[float],
) -> InteractionGraph:
    embedding_dim = int(positions.shape[1])
    resolved_radius = _resolve_locality_radius(kernel_width, locality_radius, embedding_dim)
    kernel_field = KernelField(
        kernel_width=float(kernel_width),
        locality_radius=resolved_radius,
    )
    affinity, distances = kernel_field.build_affinity(positions)
    return InteractionGraph(
        positions=np.asarray(positions, dtype=float),
        affinity=np.asarray(affinity, dtype=float),
        distances=np.asarray(distances, dtype=float),
        kernel_width=float(kernel_width),
        locality_radius=resolved_radius,
    )


def build_random_interaction_graph(
    n_nodes: int,
    kernel_width: float,
    locality_radius: Optional[float] = None,
    seed: int = 17,
    embedding_dim: int = 2,
) -> InteractionGraph:
    rng = np.random.default_rng(seed)
    positions = rng.uniform(0.0, 1.0, size=(int(n_nodes), int(embedding_dim)))
    return _build_graph_from_positions(positions, kernel_width, locality_radius)


def build_clustered_graph(
    n_nodes: int,
    kernel_width: float,
    locality_radius: Optional[float] = None,
    seed: int = 17,
    embedding_dim: int = 2,
    n_clusters: int = 3,
    cluster_spread: float = 0.08,
) -> InteractionGraph:
    rng = np.random.default_rng(seed)
    centers = rng.uniform(0.15, 0.85, size=(int(n_clusters), int(embedding_dim)))
    counts = np.full(int(n_clusters), int(n_nodes) // int(n_clusters), dtype=int)
    counts[: int(n_nodes) % int(n_clusters)] += 1
    assignments = np.repeat(np.arange(int(n_clusters)), counts)
    rng.shuffle(assignments)
    positions = centers[assignments] + rng.normal(
        loc=0.0,
        scale=float(cluster_spread),
        size=(int(n_nodes), int(embedding_dim)),
    )
    positions = np.clip(positions, 0.0, 1.0)
    return _build_graph_from_positions(positions, kernel_width, locality_radius)


def build_transport_operator(
    interaction_graph: InteractionGraph,
    self_weight: float = 0.05,
) -> TransportOperator:
    if not 0.0 <= self_weight < 1.0:
        raise ValueError("self_weight must lie in [0, 1)")

    affinity = np.asarray(interaction_graph.affinity, dtype=float)
    n_nodes = int(affinity.shape[0])
    transition = np.zeros_like(affinity, dtype=float)
    row_sums = affinity.sum(axis=1)

    nonzero_rows = row_sums > EPSILON
    transition[nonzero_rows] = affinity[nonzero_rows] / row_sums[nonzero_rows, None]
    for index in np.flatnonzero(~nonzero_rows).tolist():
        transition[index, index] = 1.0

    if self_weight > 0.0:
        transition = (1.0 - self_weight) * transition + self_weight * np.eye(n_nodes, dtype=float)

    transition = transition / np.maximum(transition.sum(axis=1, keepdims=True), EPSILON)
    return TransportOperator(transition=transition)


def apply_transport_step(
    state: np.ndarray,
    transport_operator: TransportOperator,
) -> np.ndarray:
    return transport_operator.apply_transport_step(state)

