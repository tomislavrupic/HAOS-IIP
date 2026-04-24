from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np


EPSILON = 1.0e-12


@dataclass(frozen=True)
class GeneNetwork:
    genes: tuple[str, ...]
    W: np.ndarray
    bias: np.ndarray
    decay: float
    hub_gene: str
    fragile_branch: tuple[str, ...]
    edge_weakening_targets: tuple[tuple[str, str], ...]

    @property
    def index(self) -> dict[str, int]:
        return {gene: idx for idx, gene in enumerate(self.genes)}


def sigmoid(z: np.ndarray | float) -> np.ndarray | float:
    return 1.0 / (1.0 + np.exp(-z))


def build_default_network() -> GeneNetwork:
    """Build a deterministic 12-gene toy regulatory network.

    W[source, target] stores directed regulatory influence:
    positive values activate, negative values inhibit.
    """
    genes = tuple(f"G{i}" for i in range(12))
    idx = {gene: i for i, gene in enumerate(genes)}
    W = np.zeros((len(genes), len(genes)), dtype=float)

    def edge(source: str, target: str, weight: float) -> None:
        W[idx[source], idx[target]] = weight

    # G0 is the hub regulator feeding the core and the fragile branch.
    edge("G0", "G1", 1.45)
    edge("G0", "G2", 1.05)
    edge("G0", "G3", 0.95)
    edge("G0", "G5", 1.05)
    edge("G0", "G6", 0.70)
    edge("G0", "G8", 6.00)

    # Negative feedback loop: G1 activates G2, G2 inhibits G1 and G0.
    edge("G1", "G2", 1.15)
    edge("G2", "G1", -1.35)
    edge("G2", "G0", -0.70)

    # Positive feedback motif: G3 and G4 mutually activate.
    edge("G3", "G4", 1.55)
    edge("G4", "G3", 1.35)
    edge("G2", "G4", -0.35)

    # Stabilizing downstream branch outside the fragile module.
    edge("G5", "G6", 1.10)
    edge("G6", "G7", 1.00)
    edge("G7", "G5", -0.85)

    # Fragile branch: weak external support plus internal positive maintenance.
    edge("G4", "G8", 1.50)
    edge("G5", "G9", 1.40)
    edge("G8", "G9", 3.40)
    edge("G9", "G10", 4.00)
    edge("G10", "G11", 4.00)
    edge("G8", "G10", 1.00)
    edge("G9", "G11", 1.00)
    edge("G11", "G9", 0.20)
    edge("G9", "G8", 0.20)

    bias = np.array(
        [
            0.52,   # G0
            -0.55,  # G1
            -0.38,  # G2
            -0.68,  # G3
            -0.70,  # G4
            -0.48,  # G5
            -0.56,  # G6
            -0.60,  # G7
            -2.60,  # G8
            -2.30,  # G9
            -2.40,  # G10
            -2.40,  # G11
        ],
        dtype=float,
    )

    return GeneNetwork(
        genes=genes,
        W=W,
        bias=bias,
        decay=0.36,
        hub_gene="G0",
        fragile_branch=("G8", "G9", "G10", "G11"),
        edge_weakening_targets=(("G0", "G8"), ("G4", "G8"), ("G5", "G9")),
    )


def fixed_initial_state(network: GeneNetwork) -> np.ndarray:
    return np.array(
        [0.72, 0.34, 0.31, 0.40, 0.38, 0.36, 0.32, 0.30, 0.55, 0.52, 0.50, 0.48],
        dtype=float,
    )


def simulate_network(
    network: GeneNetwork,
    *,
    steps: int = 260,
    initial_state: np.ndarray | None = None,
) -> np.ndarray:
    if initial_state is None:
        state = fixed_initial_state(network)
    else:
        state = np.asarray(initial_state, dtype=float).copy()

    trajectory = np.zeros((steps + 1, len(network.genes)), dtype=float)
    trajectory[0] = np.clip(state, 0.0, 1.0)
    for t in range(steps):
        inputs = network.bias + trajectory[t] @ network.W
        next_state = (1.0 - network.decay) * trajectory[t] + network.decay * sigmoid(inputs)
        trajectory[t + 1] = np.clip(next_state, 0.0, 1.0)
    return trajectory


def apply_edge_weakening(network: GeneNetwork, perturbation_level: float) -> GeneNetwork:
    p = float(np.clip(perturbation_level, 0.0, 1.0))
    W = network.W.copy()
    idx = network.index
    for source, target in network.edge_weakening_targets:
        W[idx[source], idx[target]] *= 1.0 - p
    return GeneNetwork(
        genes=network.genes,
        W=W,
        bias=network.bias.copy(),
        decay=network.decay,
        hub_gene=network.hub_gene,
        fragile_branch=network.fragile_branch,
        edge_weakening_targets=network.edge_weakening_targets,
    )


def apply_hub_knockdown(network: GeneNetwork, perturbation_level: float) -> GeneNetwork:
    p = float(np.clip(perturbation_level, 0.0, 1.0))
    W = network.W.copy()
    hub_idx = network.index[network.hub_gene]
    W[hub_idx, :] *= 1.0 - p
    return GeneNetwork(
        genes=network.genes,
        W=W,
        bias=network.bias.copy(),
        decay=network.decay,
        hub_gene=network.hub_gene,
        fragile_branch=network.fragile_branch,
        edge_weakening_targets=network.edge_weakening_targets,
    )


def final_window_mean(trajectory: np.ndarray, window: int = 30) -> np.ndarray:
    return trajectory[-window:].mean(axis=0)


def normalized_distance(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    return float(np.linalg.norm(left - right) / np.sqrt(left.size))


def compute_effective_graph(
    network: GeneNetwork,
    trajectory: np.ndarray,
    *,
    final_window: int = 30,
) -> tuple[dict[str, float | str], ...]:
    """Extract an activity-weighted interaction graph from the final trajectory.

    effective_weight(source, target) = abs(W[source, target]) * activity_factor
    where activity_factor is the source gene's mean activity over the final window.
    """
    source_activity = final_window_mean(trajectory, final_window)
    edges: list[dict[str, float | str]] = []
    for source_idx, source in enumerate(network.genes):
        for target_idx, target in enumerate(network.genes):
            weight = float(network.W[source_idx, target_idx])
            if abs(weight) <= EPSILON:
                continue
            activity_factor = float(source_activity[source_idx])
            edges.append(
                {
                    "source": source,
                    "target": target,
                    "raw_weight": weight,
                    "effective_weight": abs(weight) * activity_factor,
                    "activity_factor": activity_factor,
                }
            )
    return tuple(edges)


def first_sustained_collapse(
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


def compute_recoverability_metrics(
    baseline_trajectory: np.ndarray,
    perturbed_trajectories: tuple[np.ndarray, ...],
    perturbation_levels: np.ndarray,
    *,
    collapse_threshold: float = 0.65,
    sustain_steps: int = 2,
    visible_failure_threshold: float = 0.48,
    final_window: int = 30,
) -> tuple[list[dict[str, float | int | bool | None]], dict[str, float | int | None]]:
    """Compute lightweight external proxies aligned with HAOS-IIP stability logic.

    These are not frozen HAOS core metrics. Recoverability compares final-window
    expression states to the unperturbed baseline, delta_persistence is the step
    change in recoverability, k_star is the first sustained threshold crossing,
    safety_margin is distance in perturbation level to k_star, and confidence is
    the magnitude of the strongest negative recoverability step.
    """
    baseline_state = final_window_mean(baseline_trajectory, final_window)
    distances = [
        normalized_distance(final_window_mean(trajectory, final_window), baseline_state)
        for trajectory in perturbed_trajectories
    ]
    recoverability = [float(np.clip(1.0 - distance, 0.0, 1.0)) for distance in distances]
    delta_persistence = [0.0]
    for idx in range(1, len(recoverability)):
        delta_persistence.append(recoverability[idx] - recoverability[idx - 1])

    k_star = first_sustained_collapse(
        recoverability,
        collapse_threshold=collapse_threshold,
        sustain_steps=sustain_steps,
    )
    p_at_k_star = None if k_star is None else float(perturbation_levels[k_star])
    global_confidence = abs(min(delta_persistence)) if delta_persistence else 0.0

    rows: list[dict[str, float | int | bool | None]] = []
    for idx, p in enumerate(perturbation_levels):
        if p_at_k_star is None:
            safety_margin = 1.0 - float(p)
            confidence = abs(min(delta_persistence[: idx + 1]))
        else:
            safety_margin = p_at_k_star - float(p)
            confidence = global_confidence
        rows.append(
            {
                "perturbation_index": idx,
                "perturbation_level": float(p),
                "recoverability": recoverability[idx],
                "delta_persistence": delta_persistence[idx],
                "k_star": k_star,
                "safety_margin": safety_margin,
                "confidence": confidence,
                "visible_failure": bool(distances[idx] > visible_failure_threshold),
            }
        )

    summary = {
        "k_star": k_star,
        "p_at_k_star": p_at_k_star,
        "min_recoverability": min(recoverability),
        "max_distance": max(distances),
        "confidence": global_confidence,
    }
    return rows, summary
