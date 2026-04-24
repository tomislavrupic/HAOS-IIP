from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Iterable

import numpy as np


EPSILON = 1.0e-12


@dataclass(frozen=True)
class Edge:
    source: int
    target: int
    edge_class: str
    baseline_strength: float
    strength: float
    perturbation_factor: float


@dataclass(frozen=True)
class TissueLattice:
    grid_size: tuple[int, int, int]
    edges: tuple[Edge, ...]
    retention: float
    lesion_center: tuple[float, float, float]
    lesion_radius: float
    active_edge_threshold: float

    @property
    def node_count(self) -> int:
        x_size, y_size, z_size = self.grid_size
        return x_size * y_size * z_size

    def node_index(self, x: int, y: int, z: int) -> int:
        x_size, y_size, z_size = self.grid_size
        x = int(np.clip(x, 0, x_size - 1))
        y = int(np.clip(y, 0, y_size - 1))
        z = int(np.clip(z, 0, z_size - 1))
        return (x * y_size * z_size) + (y * z_size) + z

    def coordinates(self, node: int) -> tuple[int, int, int]:
        _, y_size, z_size = self.grid_size
        x = node // (y_size * z_size)
        rem = node % (y_size * z_size)
        y = rem // z_size
        z = rem % z_size
        return x, y, z


def build_edges(
    grid_size: tuple[int, int, int],
    *,
    local_neighbor_strength: float,
    diagonal_neighbor_strength: float,
    long_range_signal_strength: float,
) -> tuple[Edge, ...]:
    x_size, y_size, z_size = grid_size
    edges: list[Edge] = []

    def idx(x: int, y: int, z: int) -> int:
        return x * y_size * z_size + y * z_size + z

    def add_edge(source: int, target: int, edge_class: str, strength: float) -> None:
        edges.append(
            Edge(
                source=source,
                target=target,
                edge_class=edge_class,
                baseline_strength=strength,
                strength=strength,
                perturbation_factor=1.0,
            )
        )

    local_offsets = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    diagonal_offsets = (
        (1, 1, 0),
        (1, -1, 0),
        (1, 0, 1),
        (1, 0, -1),
        (0, 1, 1),
        (0, 1, -1),
    )

    for x in range(x_size):
        for y in range(y_size):
            for z in range(z_size):
                for dx, dy, dz in local_offsets:
                    nx, ny, nz = x + dx, y + dy, z + dz
                    if nx < x_size and ny < y_size and nz < z_size:
                        add_edge(idx(x, y, z), idx(nx, ny, nz), "local_neighbor", local_neighbor_strength)
                for dx, dy, dz in diagonal_offsets:
                    nx, ny, nz = x + dx, y + dy, z + dz
                    if 0 <= nx < x_size and 0 <= ny < y_size and 0 <= nz < z_size:
                        add_edge(idx(x, y, z), idx(nx, ny, nz), "diagonal_neighbor", diagonal_neighbor_strength)

    for x in range(0, x_size - 3, 2):
        for y in range(1, y_size - 2, 3):
            for z in range(0, z_size - 3, 3):
                add_edge(idx(x, y, z), idx(x + 3, y + 2, z + 3), "long_range_signal", long_range_signal_strength)

    return tuple(edges)


def build_tissue_lattice(
    *,
    grid_size: tuple[int, int, int] = (8, 8, 8),
    local_neighbor_strength: float = 1.0,
    diagonal_neighbor_strength: float = 0.28,
    long_range_signal_strength: float = 0.08,
    retention: float = 0.18,
    lesion_center: tuple[float, float, float] = (3.5, 3.5, 3.5),
    lesion_radius: float = 4.10,
    active_edge_threshold: float = 0.125,
) -> TissueLattice:
    return TissueLattice(
        grid_size=grid_size,
        edges=build_edges(
            grid_size,
            local_neighbor_strength=local_neighbor_strength,
            diagonal_neighbor_strength=diagonal_neighbor_strength,
            long_range_signal_strength=long_range_signal_strength,
        ),
        retention=retention,
        lesion_center=lesion_center,
        lesion_radius=lesion_radius,
        active_edge_threshold=active_edge_threshold,
    )


def lesion_weight(lattice: TissueLattice, edge: Edge) -> float:
    source = np.array(lattice.coordinates(edge.source), dtype=float)
    target = np.array(lattice.coordinates(edge.target), dtype=float)
    midpoint = 0.5 * (source + target)
    center = np.array(lattice.lesion_center, dtype=float)
    distance = float(np.linalg.norm(midpoint - center))
    if distance > lattice.lesion_radius:
        return 0.0
    shell = max(lattice.lesion_radius, EPSILON)
    return float(1.0 - 0.15 * (distance / shell))


def apply_local_lesion(lattice: TissueLattice, perturbation_level: float) -> TissueLattice:
    p = float(np.clip(perturbation_level, 0.0, 1.0))
    perturbed_edges: list[Edge] = []
    for edge in lattice.edges:
        exposure = lesion_weight(lattice, edge)
        if exposure <= EPSILON:
            perturbed_edges.append(edge)
            continue
        factor = max(0.0, 1.0 - p * exposure)
        perturbed_edges.append(
            replace(edge, strength=edge.baseline_strength * factor, perturbation_factor=factor)
        )
    return replace(lattice, edges=tuple(perturbed_edges))


def apply_signaling_weakening(lattice: TissueLattice, perturbation_level: float) -> TissueLattice:
    p = float(np.clip(perturbation_level, 0.0, 1.0))
    perturbed_edges: list[Edge] = []
    for edge in lattice.edges:
        if edge.edge_class == "long_range_signal":
            factor = 1.0 - p
            perturbed_edges.append(
                replace(edge, strength=edge.baseline_strength * factor, perturbation_factor=factor)
            )
        else:
            perturbed_edges.append(edge)
    return replace(lattice, edges=tuple(perturbed_edges))


def build_weighted_adjacency(lattice: TissueLattice) -> np.ndarray:
    adjacency = np.zeros((lattice.node_count, lattice.node_count), dtype=float)
    for edge in lattice.edges:
        if edge.strength <= EPSILON:
            continue
        adjacency[edge.source, edge.target] += edge.strength
        adjacency[edge.target, edge.source] += edge.strength
    return adjacency


def build_transition_matrix(lattice: TissueLattice) -> np.ndarray:
    adjacency = build_weighted_adjacency(lattice)
    degree = adjacency.sum(axis=1)
    transition = np.zeros_like(adjacency)
    active = degree > EPSILON
    transition[active] = adjacency[active] / degree[active, None]
    transition[~active, ~active] = 1.0
    return transition


def initial_signal(lattice: TissueLattice) -> np.ndarray:
    signal = np.zeros(lattice.node_count, dtype=float)
    x_size, y_size, z_size = lattice.grid_size
    del x_size
    for y in range(y_size):
        for z in range(z_size):
            signal[lattice.node_index(0, y, z)] = 1.0
    signal /= max(float(signal.sum()), EPSILON)
    return signal


def simulate_propagation(
    lattice: TissueLattice,
    *,
    steps: int = 40,
    initial_state: np.ndarray | None = None,
) -> np.ndarray:
    state = initial_signal(lattice) if initial_state is None else np.asarray(initial_state, dtype=float).copy()
    state = np.clip(state, 0.0, 1.0)
    total = float(state.sum())
    if total > EPSILON:
        state /= total

    transition = build_transition_matrix(lattice)
    trajectory = np.zeros((steps + 1, lattice.node_count), dtype=float)
    trajectory[0] = state
    for step in range(steps):
        neighbor_input = transition.T @ trajectory[step]
        next_state = lattice.retention * trajectory[step] + (1.0 - lattice.retention) * neighbor_input
        next_state = np.clip(next_state, 0.0, 1.0)
        total = float(next_state.sum())
        if total > EPSILON:
            next_state /= total
        trajectory[step + 1] = next_state
    return trajectory


def largest_connected_component_fraction(lattice: TissueLattice) -> float:
    adjacency = build_weighted_adjacency(lattice)
    active = adjacency >= lattice.active_edge_threshold
    seen = np.zeros(lattice.node_count, dtype=bool)
    largest = 0
    for node in range(lattice.node_count):
        if seen[node]:
            continue
        stack = [node]
        seen[node] = True
        size = 0
        while stack:
            current = stack.pop()
            size += 1
            neighbors = np.flatnonzero(active[current])
            for neighbor in neighbors:
                if not seen[neighbor]:
                    seen[neighbor] = True
                    stack.append(int(neighbor))
        largest = max(largest, size)
    return float(largest / lattice.node_count)


def weighted_degree_retention(lattice: TissueLattice, baseline_weighted_degree: np.ndarray) -> float:
    current = build_weighted_adjacency(lattice).sum(axis=1)
    baseline_mean = float(baseline_weighted_degree.mean())
    if baseline_mean <= EPSILON:
        return 0.0
    return float(np.clip(current.mean() / baseline_mean, 0.0, 1.0))


def compute_spatial_coherence(
    lattice: TissueLattice,
    baseline_weighted_degree: np.ndarray,
) -> float:
    current = build_weighted_adjacency(lattice).sum(axis=1)
    baseline = np.maximum(baseline_weighted_degree, EPSILON)
    local_weights = np.zeros(lattice.node_count, dtype=float)
    center = np.array(lattice.lesion_center, dtype=float)
    for node in range(lattice.node_count):
        coords = np.array(lattice.coordinates(node), dtype=float)
        distance = float(np.linalg.norm(coords - center))
        if distance <= lattice.lesion_radius + 1.0:
            local_weights[node] = 1.0 - 0.30 * min(distance / (lattice.lesion_radius + 1.0), 1.0)
    if local_weights.sum() <= EPSILON:
        return 1.0
    local_support = float(np.average(np.clip(current / baseline, 0.0, 1.0), weights=local_weights))
    # Morphological coherence is less brittle than raw local support loss.
    return float(np.clip(1.0 - 0.55 * (1.0 - local_support), 0.0, 1.0))


def propagation_retention(baseline_final: np.ndarray, current_final: np.ndarray) -> float:
    baseline = baseline_final / max(float(np.linalg.norm(baseline_final)), EPSILON)
    current = current_final / max(float(np.linalg.norm(current_final)), EPSILON)
    return float(np.clip(np.dot(baseline, current), 0.0, 1.0))


def find_k_star(
    recoverability: Iterable[float],
    *,
    collapse_threshold: float,
    sustain_steps: int,
) -> int | None:
    values = list(recoverability)
    for idx, value in enumerate(values):
        stop = idx + sustain_steps + 1
        if stop > len(values):
            break
        if value < collapse_threshold and all(v < collapse_threshold for v in values[idx + 1 : stop]):
            return idx
    return None


def compute_safety_margin(perturbation_level: float, p_at_k_star: float | None) -> float:
    if p_at_k_star is None:
        return 1.0 - perturbation_level
    return p_at_k_star - perturbation_level


def compute_effective_graph(lattice: TissueLattice) -> tuple[dict[str, float | int | str], ...]:
    edges: list[dict[str, float | int | str]] = []
    for edge in lattice.edges:
        if edge.strength <= EPSILON:
            continue
        edges.append(
            {
                "source": edge.source,
                "target": edge.target,
                "edge_class": edge.edge_class,
                "effective_weight": edge.strength,
                "perturbation_factor": edge.perturbation_factor,
            }
        )
    return tuple(edges)


def compute_recoverability_metrics(
    baseline_lattice: TissueLattice,
    perturbed_lattices: tuple[TissueLattice, ...],
    perturbation_levels: np.ndarray,
    *,
    collapse_threshold: float = 0.70,
    sustain_steps: int = 2,
    propagation_steps: int = 40,
) -> tuple[list[dict[str, float | int | bool | None]], dict[str, float | int | bool | None]]:
    baseline_adjacency = build_weighted_adjacency(baseline_lattice)
    baseline_weighted_degree = baseline_adjacency.sum(axis=1)
    baseline_trajectory = simulate_propagation(baseline_lattice, steps=propagation_steps)
    baseline_final = baseline_trajectory[-1]

    raw_rows: list[dict[str, float | bool]] = []
    recoverability_values: list[float] = []
    for lattice in perturbed_lattices:
        component_fraction = largest_connected_component_fraction(lattice)
        degree_retention = weighted_degree_retention(lattice, baseline_weighted_degree)
        spatial_retention = compute_spatial_coherence(lattice, baseline_weighted_degree)
        trajectory = simulate_propagation(lattice, steps=propagation_steps)
        signal_retention = propagation_retention(baseline_final, trajectory[-1])
        recoverability = float(
            np.clip(
                0.15 * component_fraction
                + 0.45 * degree_retention
                + 0.25 * spatial_retention
                + 0.15 * signal_retention,
                0.0,
                1.0,
            )
        )
        raw_rows.append(
            {
                "recoverability": recoverability,
                "largest_component_fraction": component_fraction,
                "weighted_degree_retention": degree_retention,
                "spatial_coherence_retention": spatial_retention,
                "propagation_retention": signal_retention,
                "visible_failure": bool(component_fraction < 0.85 or spatial_retention < 0.50),
            }
        )
        recoverability_values.append(recoverability)

    delta_persistence = [0.0]
    for idx in range(1, len(recoverability_values)):
        delta_persistence.append(recoverability_values[idx] - recoverability_values[idx - 1])

    k_star = find_k_star(
        recoverability_values,
        collapse_threshold=collapse_threshold,
        sustain_steps=sustain_steps,
    )
    p_at_k_star = None if k_star is None else float(perturbation_levels[k_star])

    rows: list[dict[str, float | int | bool | None]] = []
    for idx, raw in enumerate(raw_rows):
        p = float(perturbation_levels[idx])
        rows.append(
            {
                "perturbation_index": idx,
                "perturbation_level": p,
                "recoverability": raw["recoverability"],
                "delta_persistence": delta_persistence[idx],
                "k_star": k_star,
                "safety_margin": compute_safety_margin(p, p_at_k_star),
                "confidence": abs(min(delta_persistence[: idx + 1])),
                "largest_component_fraction": raw["largest_component_fraction"],
                "weighted_degree_retention": raw["weighted_degree_retention"],
                "spatial_coherence_retention": raw["spatial_coherence_retention"],
                "propagation_retention": raw["propagation_retention"],
                "visible_failure": raw["visible_failure"],
            }
        )

    summary = {
        "k_star": k_star,
        "p_at_k_star": p_at_k_star,
        "baseline_recoverability": recoverability_values[0],
        "baseline_largest_component_fraction": raw_rows[0]["largest_component_fraction"],
        "baseline_weighted_degree_retention": raw_rows[0]["weighted_degree_retention"],
        "baseline_spatial_coherence_retention": raw_rows[0]["spatial_coherence_retention"],
        "baseline_propagation_retention": raw_rows[0]["propagation_retention"],
        "min_recoverability": min(recoverability_values),
    }
    return rows, summary


def run_perturbation_sweep(
    lattice: TissueLattice,
    perturbation_levels: np.ndarray,
    *,
    mode: str = "local_lesion",
    collapse_threshold: float = 0.70,
    sustain_steps: int = 2,
) -> tuple[list[dict[str, float | int | bool | None]], dict[str, float | int | bool | None]]:
    perturbed: list[TissueLattice] = []
    for level in perturbation_levels:
        if mode == "local_lesion":
            perturbed.append(apply_local_lesion(lattice, float(level)))
        elif mode == "signaling_weakening":
            perturbed.append(apply_signaling_weakening(lattice, float(level)))
        elif mode == "lesion_plus_signal":
            perturbed.append(apply_signaling_weakening(apply_local_lesion(lattice, float(level)), float(level)))
        else:
            raise ValueError(f"unknown perturbation mode: {mode}")
    return compute_recoverability_metrics(
        lattice,
        tuple(perturbed),
        perturbation_levels,
        collapse_threshold=collapse_threshold,
        sustain_steps=sustain_steps,
    )
