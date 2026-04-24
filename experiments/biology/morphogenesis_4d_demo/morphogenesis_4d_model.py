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
class MorphogenesisModel:
    grid_size: tuple[int, int, int]
    developmental_steps: int
    edges: tuple[Edge, ...]
    retention: float
    neighbor_weight: float
    target_weight: float
    drift_level: float
    lesion_level: float
    drift_center: tuple[float, float, float]
    drift_radius: float
    delay_steps: int
    lesion_center: tuple[float, float, float]
    lesion_radius: float

    @property
    def node_count(self) -> int:
        x_size, y_size, z_size = self.grid_size
        return x_size * y_size * z_size

    def node_index(self, x: int, y: int, z: int) -> int:
        x_size, y_size, z_size = self.grid_size
        x = int(np.clip(x, 0, x_size - 1))
        y = int(np.clip(y, 0, y_size - 1))
        z = int(np.clip(z, 0, z_size - 1))
        return x * y_size * z_size + y * z_size + z

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
    developmental_signal_strength: float,
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
        for y in range(0, y_size - 2, 2):
            for z in range(0, z_size - 2, 2):
                add_edge(idx(x, y, z), idx(x + 3, y + 2, z + 2), "developmental_signal", developmental_signal_strength)
    return tuple(edges)


def build_tissue_lattice(
    *,
    grid_size: tuple[int, int, int] = (8, 8, 8),
    developmental_steps: int = 24,
    local_neighbor_strength: float = 1.0,
    diagonal_neighbor_strength: float = 0.22,
    developmental_signal_strength: float = 0.32,
    retention: float = 0.18,
    neighbor_weight: float = 0.17,
    target_weight: float = 0.65,
    drift_center: tuple[float, float, float] = (5.1, 3.5, 3.5),
    drift_radius: float = 4.4,
    delay_steps: int = 10,
    lesion_center: tuple[float, float, float] = (3.5, 3.5, 3.5),
    lesion_radius: float = 2.6,
) -> MorphogenesisModel:
    return MorphogenesisModel(
        grid_size=grid_size,
        developmental_steps=developmental_steps,
        edges=build_edges(
            grid_size,
            local_neighbor_strength=local_neighbor_strength,
            diagonal_neighbor_strength=diagonal_neighbor_strength,
            developmental_signal_strength=developmental_signal_strength,
        ),
        retention=retention,
        neighbor_weight=neighbor_weight,
        target_weight=target_weight,
        drift_level=0.0,
        lesion_level=0.0,
        drift_center=drift_center,
        drift_radius=drift_radius,
        delay_steps=delay_steps,
        lesion_center=lesion_center,
        lesion_radius=lesion_radius,
    )


def smoothstep(value: float) -> float:
    x = float(np.clip(value, 0.0, 1.0))
    return x * x * (3.0 - 2.0 * x)


def normalized_positions(model: MorphogenesisModel) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    coords = np.array([model.coordinates(node) for node in range(model.node_count)], dtype=float)
    scales = np.maximum(np.array(model.grid_size, dtype=float) - 1.0, 1.0)
    normalized = coords / scales
    centered = 2.0 * normalized - 1.0
    radius = np.linalg.norm(centered, axis=1) / np.sqrt(3.0)
    return normalized[:, 0], centered[:, 0], radius, coords


def compute_target_morphology(model: MorphogenesisModel, time_index: int) -> np.ndarray:
    if model.developmental_steps <= 1:
        tau = 1.0
    else:
        tau = float(np.clip(time_index / (model.developmental_steps - 1), 0.0, 1.0))
    x_norm, x_centered, radius, _ = normalized_positions(model)
    core = np.exp(-3.0 * radius * radius)
    shell = np.exp(-((radius - 0.62) ** 2) / (2.0 * 0.16 * 0.16))
    early = 0.16 + 0.78 * core
    mid = 0.14 + 0.52 * core + 0.30 * x_norm
    late = 0.12 + 0.42 * shell + 0.34 * x_norm + 0.08 * (1.0 - np.abs(x_centered))
    first_blend = smoothstep(tau / 0.52)
    second_blend = smoothstep((tau - 0.42) / 0.58)
    target = (1.0 - first_blend) * early + first_blend * mid
    target = (1.0 - second_blend) * target + second_blend * late
    return np.clip(target, 0.0, 1.0)


def region_exposure(
    model: MorphogenesisModel,
    node: int,
    *,
    center: tuple[float, float, float],
    radius: float,
) -> float:
    coords = np.array(model.coordinates(node), dtype=float)
    distance = float(np.linalg.norm(coords - np.array(center, dtype=float)))
    if distance > radius:
        return 0.0
    return float(1.0 - 0.20 * (distance / max(radius, EPSILON)))


def drift_mask(model: MorphogenesisModel) -> np.ndarray:
    return np.array(
        [
            region_exposure(model, node, center=model.drift_center, radius=model.drift_radius)
            for node in range(model.node_count)
        ],
        dtype=float,
    )


def lesion_edge_exposure(model: MorphogenesisModel, edge: Edge) -> float:
    source = np.array(model.coordinates(edge.source), dtype=float)
    target = np.array(model.coordinates(edge.target), dtype=float)
    midpoint = 0.5 * (source + target)
    distance = float(np.linalg.norm(midpoint - np.array(model.lesion_center, dtype=float)))
    if distance > model.lesion_radius:
        return 0.0
    return float(1.0 - 0.20 * (distance / max(model.lesion_radius, EPSILON)))


def drift_edge_exposure(model: MorphogenesisModel, edge: Edge) -> float:
    if edge.edge_class != "developmental_signal":
        return 0.0
    return max(
        region_exposure(model, edge.source, center=model.drift_center, radius=model.drift_radius),
        region_exposure(model, edge.target, center=model.drift_center, radius=model.drift_radius),
    )


def apply_developmental_drift(model: MorphogenesisModel, perturbation_level: float) -> MorphogenesisModel:
    p = float(np.clip(perturbation_level, 0.0, 1.0))
    perturbed_edges: list[Edge] = []
    for edge in model.edges:
        exposure = drift_edge_exposure(model, edge)
        if exposure <= EPSILON:
            perturbed_edges.append(edge)
            continue
        factor = max(0.0, 1.0 - 0.85 * p * exposure)
        perturbed_edges.append(
            replace(edge, strength=edge.baseline_strength * factor, perturbation_factor=factor)
        )
    return replace(model, edges=tuple(perturbed_edges), drift_level=p)


def apply_local_developmental_lesion(model: MorphogenesisModel, perturbation_level: float) -> MorphogenesisModel:
    p = float(np.clip(perturbation_level, 0.0, 1.0))
    perturbed_edges: list[Edge] = []
    for edge in model.edges:
        exposure = lesion_edge_exposure(model, edge)
        if exposure <= EPSILON:
            perturbed_edges.append(edge)
            continue
        factor = max(0.0, 1.0 - p * exposure)
        perturbed_edges.append(
            replace(edge, strength=edge.baseline_strength * factor, perturbation_factor=edge.perturbation_factor * factor)
        )
    return replace(model, edges=tuple(perturbed_edges), lesion_level=p)


def build_weighted_adjacency(model: MorphogenesisModel) -> np.ndarray:
    adjacency = np.zeros((model.node_count, model.node_count), dtype=float)
    for edge in model.edges:
        if edge.strength <= EPSILON:
            continue
        adjacency[edge.source, edge.target] += edge.strength
        adjacency[edge.target, edge.source] += edge.strength
    return adjacency


def build_transition_matrix(model: MorphogenesisModel) -> np.ndarray:
    adjacency = build_weighted_adjacency(model)
    degree = adjacency.sum(axis=1)
    transition = np.zeros_like(adjacency)
    active = degree > EPSILON
    transition[active] = adjacency[active] / degree[active, None]
    transition[~active, ~active] = 1.0
    return transition


def effective_target(model: MorphogenesisModel, time_index: int) -> np.ndarray:
    target = compute_target_morphology(model, time_index)
    if model.drift_level <= EPSILON:
        return target
    delayed_index = max(0, time_index - model.delay_steps)
    delayed_target = compute_target_morphology(model, delayed_index)
    neutral_target = np.full(model.node_count, 0.34, dtype=float)
    tau = time_index / max(model.developmental_steps - 1, 1)
    time_gate = smoothstep((tau - 0.18) / 0.62)
    drift_strength = np.clip(model.drift_level * time_gate * drift_mask(model), 0.0, 1.0)
    drifted_target = 0.72 * delayed_target + 0.28 * neutral_target
    return np.clip((1.0 - drift_strength) * target + drift_strength * drifted_target, 0.0, 1.0)


def simulate_development(model: MorphogenesisModel) -> np.ndarray:
    transition = build_transition_matrix(model)
    trajectory = np.zeros((model.developmental_steps, model.node_count), dtype=float)
    trajectory[0] = compute_target_morphology(model, 0)
    weight_sum = max(model.retention + model.neighbor_weight + model.target_weight, EPSILON)
    for time_index in range(model.developmental_steps - 1):
        neighbor_input = transition @ trajectory[time_index]
        target_pull = effective_target(model, time_index + 1)
        next_state = (
            model.retention * trajectory[time_index]
            + model.neighbor_weight * neighbor_input
            + model.target_weight * target_pull
        ) / weight_sum
        trajectory[time_index + 1] = np.clip(next_state, 0.0, 1.0)
    return trajectory


def match_score(left: np.ndarray, right: np.ndarray, *, scale: float = 0.30) -> float:
    rmse = float(np.sqrt(np.mean((left - right) ** 2)))
    return float(np.clip(1.0 - rmse / scale, 0.0, 1.0))


def compute_final_morphology_match(model: MorphogenesisModel, trajectory: np.ndarray) -> float:
    return match_score(trajectory[-1], compute_target_morphology(model, model.developmental_steps - 1))


def compute_trajectory_coherence(model: MorphogenesisModel, trajectory: np.ndarray) -> float:
    scores = [
        match_score(trajectory[time_index], compute_target_morphology(model, time_index))
        for time_index in range(model.developmental_steps)
    ]
    return float(np.mean(scores))


def compute_spatial_continuity(model: MorphogenesisModel, trajectory: np.ndarray, baseline_final_state: np.ndarray) -> float:
    final_state = trajectory[-1]
    baseline_gradient = local_gradient_mean(model, baseline_final_state)
    current_gradient = local_gradient_mean(model, final_state)
    if baseline_gradient <= EPSILON:
        return 1.0
    excess = max(0.0, current_gradient - baseline_gradient)
    return float(np.clip(1.0 - excess / 0.22, 0.0, 1.0))


def local_gradient_mean(model: MorphogenesisModel, state: np.ndarray) -> float:
    diffs: list[float] = []
    for edge in model.edges:
        if edge.edge_class == "local_neighbor":
            diffs.append(abs(float(state[edge.source] - state[edge.target])))
    return float(np.mean(diffs)) if diffs else 0.0


def compute_interaction_support_retention(model: MorphogenesisModel, baseline_weighted_degree: np.ndarray) -> float:
    current_degree = build_weighted_adjacency(model).sum(axis=1)
    degree_score = float(np.clip(current_degree.mean() / max(float(baseline_weighted_degree.mean()), EPSILON), 0.0, 1.0))
    target_score = float(np.clip(1.0 - 0.82 * model.drift_level, 0.0, 1.0))
    return float(np.clip(0.55 * degree_score + 0.45 * target_score, 0.0, 1.0))


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


def compute_recoverability_metrics(
    baseline_model: MorphogenesisModel,
    perturbed_models: tuple[MorphogenesisModel, ...],
    perturbation_levels: np.ndarray,
    *,
    collapse_threshold: float = 0.70,
    sustain_steps: int = 2,
) -> tuple[list[dict[str, float | int | bool | None]], dict[str, float | int | bool | None]]:
    baseline_weighted_degree = build_weighted_adjacency(baseline_model).sum(axis=1)
    baseline_trajectory = simulate_development(baseline_model)
    baseline_final_state = baseline_trajectory[-1]

    raw_rows: list[dict[str, float | bool]] = []
    recoverability_values: list[float] = []
    for model in perturbed_models:
        trajectory = simulate_development(model)
        final_match = compute_final_morphology_match(baseline_model, trajectory)
        trajectory_coherence = compute_trajectory_coherence(baseline_model, trajectory)
        support_retention = compute_interaction_support_retention(model, baseline_weighted_degree)
        spatial_continuity = compute_spatial_continuity(baseline_model, trajectory, baseline_final_state)
        recoverability = float(
            np.clip(
                0.35 * final_match
                + 0.32 * trajectory_coherence
                + 0.25 * support_retention
                + 0.08 * spatial_continuity,
                0.0,
                1.0,
            )
        )
        raw_rows.append(
            {
                "recoverability": recoverability,
                "final_morphology_match": final_match,
                "trajectory_coherence": trajectory_coherence,
                "interaction_support_retention": support_retention,
                "spatial_continuity_retention": spatial_continuity,
                "visible_failure": bool(final_match < 0.50 or spatial_continuity < 0.50),
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
                "final_morphology_match": raw["final_morphology_match"],
                "trajectory_coherence": raw["trajectory_coherence"],
                "interaction_support_retention": raw["interaction_support_retention"],
                "spatial_continuity_retention": raw["spatial_continuity_retention"],
                "visible_failure": raw["visible_failure"],
            }
        )

    summary = {
        "k_star": k_star,
        "p_at_k_star": p_at_k_star,
        "baseline_recoverability": recoverability_values[0],
        "baseline_final_morphology_match": raw_rows[0]["final_morphology_match"],
        "baseline_trajectory_coherence": raw_rows[0]["trajectory_coherence"],
        "baseline_interaction_support_retention": raw_rows[0]["interaction_support_retention"],
        "baseline_spatial_continuity_retention": raw_rows[0]["spatial_continuity_retention"],
        "min_recoverability": min(recoverability_values),
    }
    return rows, summary


def run_perturbation_sweep(
    model: MorphogenesisModel,
    perturbation_levels: np.ndarray,
    *,
    mode: str = "developmental_drift",
    collapse_threshold: float = 0.70,
    sustain_steps: int = 2,
) -> tuple[list[dict[str, float | int | bool | None]], dict[str, float | int | bool | None]]:
    perturbed: list[MorphogenesisModel] = []
    for level in perturbation_levels:
        if mode == "developmental_drift":
            perturbed.append(apply_developmental_drift(model, float(level)))
        elif mode == "local_developmental_lesion":
            perturbed.append(apply_local_developmental_lesion(model, float(level)))
        elif mode == "drift_plus_lesion":
            perturbed.append(apply_local_developmental_lesion(apply_developmental_drift(model, float(level)), float(level)))
        else:
            raise ValueError(f"unknown perturbation mode: {mode}")
    return compute_recoverability_metrics(
        model,
        tuple(perturbed),
        perturbation_levels,
        collapse_threshold=collapse_threshold,
        sustain_steps=sustain_steps,
    )
