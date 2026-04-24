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
class MicrotubuleLattice:
    protofilaments: int
    dimers_per_protofilament: int
    edges: tuple[Edge, ...]
    lambda_locality: float
    axial_spacing: float
    active_edge_threshold: float
    self_retention: float

    @property
    def node_count(self) -> int:
        return self.protofilaments * self.dimers_per_protofilament

    def node_index(self, protofilament: int, z_index: int) -> int:
        p = protofilament % self.protofilaments
        z = int(np.clip(z_index, 0, self.dimers_per_protofilament - 1))
        return p * self.dimers_per_protofilament + z

    def coordinates(self, node: int) -> tuple[int, int, float]:
        p = node // self.dimers_per_protofilament
        z = node % self.dimers_per_protofilament
        theta = 2.0 * np.pi * p / self.protofilaments
        return p, z, theta


def build_microtubule_lattice(
    *,
    protofilaments: int = 13,
    dimers_per_protofilament: int = 24,
    longitudinal_strength: float = 1.0,
    lateral_strength: float = 0.82,
    diagonal_strength: float = 0.06,
    weak_support_strength: float = 0.035,
    lambda_locality: float = 0.78,
    axial_spacing: float = 0.32,
    active_edge_threshold: float = 0.060,
    self_retention: float = 0.12,
) -> MicrotubuleLattice:
    edges: list[Edge] = []

    def idx(p: int, z: int) -> int:
        return (p % protofilaments) * dimers_per_protofilament + z

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

    for p in range(protofilaments):
        for z in range(dimers_per_protofilament - 1):
            add_edge(idx(p, z), idx(p, z + 1), "longitudinal", longitudinal_strength)

    for p in range(protofilaments):
        for z in range(dimers_per_protofilament):
            add_edge(idx(p, z), idx(p + 1, z), "lateral", lateral_strength)

    for p in range(protofilaments):
        for z in range(dimers_per_protofilament - 1):
            add_edge(idx(p, z), idx(p + 1, z + 1), "seam_or_diagonal", diagonal_strength)

    for p in range(0, protofilaments, 2):
        for z in range(0, dimers_per_protofilament - 4, 4):
            add_edge(idx(p, z), idx(p + 3, z + 4), "weak_support", weak_support_strength)

    return MicrotubuleLattice(
        protofilaments=protofilaments,
        dimers_per_protofilament=dimers_per_protofilament,
        edges=tuple(edges),
        lambda_locality=lambda_locality,
        axial_spacing=axial_spacing,
        active_edge_threshold=active_edge_threshold,
        self_retention=self_retention,
    )


def compute_cylindrical_distance(lattice: MicrotubuleLattice, left: int, right: int) -> float:
    left_p, left_z, left_theta = lattice.coordinates(left)
    right_p, right_z, right_theta = lattice.coordinates(right)
    del left_p, right_p
    angular = abs(left_theta - right_theta)
    angular = min(angular, 2.0 * np.pi - angular)
    axial = lattice.axial_spacing * abs(left_z - right_z)
    return float(np.sqrt(angular * angular + axial * axial))


def locality_stability(lattice: MicrotubuleLattice, edge: Edge) -> float:
    distance = compute_cylindrical_distance(lattice, edge.source, edge.target)
    return float(np.exp(-distance / lattice.lambda_locality))


def edge_score(lattice: MicrotubuleLattice, edge: Edge) -> float:
    return edge.baseline_strength * locality_stability(lattice, edge) * edge.perturbation_factor


def compute_filtration_order(lattice: MicrotubuleLattice) -> tuple[Edge, ...]:
    """Order strongest and most local interaction supports first."""
    return tuple(sorted(lattice.edges, key=lambda edge: edge_score(lattice, edge), reverse=True))


def apply_lateral_weakening(lattice: MicrotubuleLattice, perturbation_level: float) -> MicrotubuleLattice:
    p = float(np.clip(perturbation_level, 0.0, 1.0))
    weakened_edges: list[Edge] = []
    for edge in lattice.edges:
        if edge.edge_class == "lateral":
            factor = 1.0 - p
            weakened_edges.append(
                replace(edge, strength=edge.baseline_strength * factor, perturbation_factor=factor)
            )
        else:
            weakened_edges.append(edge)
    return replace(lattice, edges=tuple(weakened_edges))


def apply_defect_patch(
    lattice: MicrotubuleLattice,
    perturbation_level: float,
    *,
    protofilament_range: range = range(4, 7),
    z_range: range = range(10, 15),
) -> MicrotubuleLattice:
    p = float(np.clip(perturbation_level, 0.0, 1.0))
    patch_nodes = {lattice.node_index(protofilament, z) for protofilament in protofilament_range for z in z_range}
    weakened_edges: list[Edge] = []
    for edge in lattice.edges:
        if edge.source in patch_nodes or edge.target in patch_nodes:
            factor = 1.0 - p
            weakened_edges.append(
                replace(edge, strength=edge.strength * factor, perturbation_factor=edge.perturbation_factor * factor)
            )
        else:
            weakened_edges.append(edge)
    return replace(lattice, edges=tuple(weakened_edges))


def build_weighted_adjacency(lattice: MicrotubuleLattice) -> np.ndarray:
    adjacency = np.zeros((lattice.node_count, lattice.node_count), dtype=float)
    for edge in lattice.edges:
        score = edge_score(lattice, edge)
        if score <= EPSILON:
            continue
        adjacency[edge.source, edge.target] += score
        adjacency[edge.target, edge.source] += score
    return adjacency


def build_transition_matrix(lattice: MicrotubuleLattice) -> np.ndarray:
    adjacency = build_weighted_adjacency(lattice)
    degree = adjacency.sum(axis=1)
    transition = np.zeros_like(adjacency)
    active = degree > EPSILON
    transition[active] = adjacency[active] / degree[active, None]
    transition[~active, ~active] = 1.0
    return (1.0 - lattice.self_retention) * transition + lattice.self_retention * np.eye(lattice.node_count)


def initial_signal(lattice: MicrotubuleLattice) -> np.ndarray:
    signal = np.zeros(lattice.node_count, dtype=float)
    signal[lattice.node_index(0, 0)] = 1.0
    return signal


def simulate_propagation(
    lattice: MicrotubuleLattice,
    *,
    steps: int = 54,
    initial_state: np.ndarray | None = None,
) -> np.ndarray:
    state = initial_signal(lattice) if initial_state is None else np.asarray(initial_state, dtype=float).copy()
    state_sum = state.sum()
    if state_sum > EPSILON:
        state = state / state_sum
    transition = build_transition_matrix(lattice)
    trajectory = np.zeros((steps + 1, lattice.node_count), dtype=float)
    trajectory[0] = state
    for step in range(steps):
        trajectory[step + 1] = transition.T @ trajectory[step]
        total = trajectory[step + 1].sum()
        if total > EPSILON:
            trajectory[step + 1] /= total
    return trajectory


def weighted_degree_retention(lattice: MicrotubuleLattice, baseline_weighted_degree: np.ndarray) -> float:
    current = build_weighted_adjacency(lattice).sum(axis=1)
    baseline_mean = float(baseline_weighted_degree.mean())
    if baseline_mean <= EPSILON:
        return 0.0
    return float(np.clip(current.mean() / baseline_mean, 0.0, 1.0))


def largest_connected_component_fraction(lattice: MicrotubuleLattice) -> float:
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


def compute_safety_margin(
    perturbation_level: float,
    p_at_k_star: float | None,
) -> float:
    if p_at_k_star is None:
        return 1.0 - perturbation_level
    return p_at_k_star - perturbation_level


def compute_recoverability_metrics(
    baseline_lattice: MicrotubuleLattice,
    perturbed_lattices: tuple[MicrotubuleLattice, ...],
    perturbation_levels: np.ndarray,
    *,
    collapse_threshold: float = 0.70,
    sustain_steps: int = 2,
    propagation_steps: int = 54,
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
        trajectory = simulate_propagation(lattice, steps=propagation_steps)
        signal_retention = propagation_retention(baseline_final, trajectory[-1])
        recoverability = float(
            np.clip(
                0.05 * component_fraction + 0.90 * degree_retention + 0.05 * signal_retention,
                0.0,
                1.0,
            )
        )
        raw_rows.append(
            {
                "recoverability": recoverability,
                "largest_component_fraction": component_fraction,
                "weighted_degree_retention": degree_retention,
                "propagation_retention": signal_retention,
                "visible_failure": bool(component_fraction < 0.85 or signal_retention < 0.50),
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
        confidence = abs(min(delta_persistence[: idx + 1]))
        rows.append(
            {
                "perturbation_index": idx,
                "perturbation_level": p,
                "recoverability": raw["recoverability"],
                "delta_persistence": delta_persistence[idx],
                "k_star": k_star,
                "safety_margin": compute_safety_margin(p, p_at_k_star),
                "confidence": confidence,
                "largest_component_fraction": raw["largest_component_fraction"],
                "weighted_degree_retention": raw["weighted_degree_retention"],
                "propagation_retention": raw["propagation_retention"],
                "visible_failure": raw["visible_failure"],
            }
        )

    summary = {
        "k_star": k_star,
        "p_at_k_star": p_at_k_star,
        "min_recoverability": min(recoverability_values),
        "baseline_recoverability": recoverability_values[0],
        "baseline_largest_component_fraction": raw_rows[0]["largest_component_fraction"],
        "baseline_weighted_degree_retention": raw_rows[0]["weighted_degree_retention"],
        "baseline_propagation_retention": raw_rows[0]["propagation_retention"],
    }
    return rows, summary


def run_perturbation_sweep(
    lattice: MicrotubuleLattice,
    perturbation_levels: np.ndarray,
    *,
    mode: str = "lateral_weakening",
    collapse_threshold: float = 0.70,
    sustain_steps: int = 2,
) -> tuple[list[dict[str, float | int | bool | None]], dict[str, float | int | bool | None]]:
    perturbed: list[MicrotubuleLattice] = []
    for level in perturbation_levels:
        if mode == "lateral_weakening":
            perturbed.append(apply_lateral_weakening(lattice, float(level)))
        elif mode == "defect_patch":
            perturbed.append(apply_defect_patch(lattice, float(level)))
        elif mode == "lateral_plus_defect":
            perturbed.append(apply_defect_patch(apply_lateral_weakening(lattice, float(level)), float(level)))
        else:
            raise ValueError(f"unknown perturbation mode: {mode}")
    return compute_recoverability_metrics(
        lattice,
        tuple(perturbed),
        perturbation_levels,
        collapse_threshold=collapse_threshold,
        sustain_steps=sustain_steps,
    )
