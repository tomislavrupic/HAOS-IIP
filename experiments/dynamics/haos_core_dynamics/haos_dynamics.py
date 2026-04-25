"""Minimal native dynamics probe for HAOS-IIP sidecar experiments.

The rig evolves a scalar field on a HAOS-derived interaction graph using a
bounded local Laplacian update plus a reference recovery term. It measures
recoverability, invariant drift, causal spread, and deterministic controls.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parent.parent.parent
EPS = 1.0e-12
RECOVERABILITY_FLOOR = 0.72
INVARIANT_DRIFT_CEILING = 0.18
CONTROL_MARGIN = 0.06


@dataclass(frozen=True)
class DynamicsConfig:
    """Runtime configuration for the HAOS dynamics probe."""

    output_dir: Path = ROOT / "outputs"
    seed: int = 20260425
    steps: int = 160
    dt: float = 0.08
    diffusion: float = 0.55
    recovery_gain: float = 0.18
    perturbation_step: int = 24
    perturbation_scale: float = 0.75
    damage_fraction: float = 0.0
    max_nodes: int = 96


@dataclass(frozen=True)
class GraphBundle:
    adjacency: np.ndarray
    positions: np.ndarray
    source_label: str
    source_kind: str

    @property
    def node_count(self) -> int:
        return int(self.adjacency.shape[0])

    @property
    def edge_count(self) -> int:
        return int(np.count_nonzero(np.triu(self.adjacency, k=1)))


@dataclass(frozen=True)
class DynamicsRun:
    label: str
    states: np.ndarray
    metrics: list[dict[str, float | int | str]]
    summary: dict[str, float | int | bool | str | None]


@dataclass(frozen=True)
class DynamicsOutcome:
    status: str
    graph_source: str
    node_count: int
    edge_count: int
    summary: dict[str, Any]
    output_dir: Path


def run_dynamics_probe(config: DynamicsConfig) -> DynamicsOutcome:
    """Run observed dynamics plus deterministic controls and write outputs."""

    config.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(config.seed)
    graph = load_geometry_graph(config, rng)
    graph = damage_graph(graph, config.damage_fraction, rng)
    reference = build_reference_field(graph.positions)
    perturbation_nodes = choose_perturbation_nodes(graph.positions, rng)

    observed = run_single_dynamics("observed", graph.adjacency, reference, perturbation_nodes, config, rng)
    controls = [
        run_single_dynamics(
            "state_shuffle_control",
            graph.adjacency,
            reference,
            perturbation_nodes,
            config,
            rng,
            initial_override=rng.permutation(reference),
        ),
        run_single_dynamics("edge_shuffle_control", edge_shuffle(graph.adjacency, rng), reference, perturbation_nodes, config, rng),
        run_single_dynamics("topology_randomization_control", topology_randomize(graph.adjacency, rng), reference, perturbation_nodes, config, rng),
    ]
    status = classify_status(graph, observed, controls)
    write_outputs(config.output_dir, graph, observed, controls, status)
    outcome_summary = dict(observed.summary)
    outcome_summary["control_contrast"] = status["control_contrast"]
    outcome_summary["failure_mode"] = status["failure_mode"]
    outcome_summary["best_control_recoverability_score"] = status["best_control_recoverability_score"]
    return DynamicsOutcome(status["bridge_status"], graph.source_label, graph.node_count, graph.edge_count, outcome_summary, config.output_dir)


def load_geometry_graph(config: DynamicsConfig, rng: np.random.Generator) -> GraphBundle:
    """Reconstruct the committed geometry-emergence graph without cross-imports."""

    path = REPO_ROOT / "geometry_emergence" / "experiments" / "default_geometry_probe_config.json"
    if not path.exists():
        return synthetic_ring_graph(min(48, config.max_nodes))

    payload = json.loads(path.read_text(encoding="utf-8"))
    n_nodes = min(int(payload.get("n_nodes", 72)), config.max_nodes)
    dim = int(payload.get("embedding_dim", 2))
    clusters = int(payload.get("cluster_count", 3))
    spread = float(payload.get("cluster_spread", 0.075))
    widths = payload.get("kernel_widths", [0.13])
    kernel_width = float(widths[len(widths) // 2])
    locality_radius = float(payload.get("locality_radius_factor", 3.0)) * kernel_width
    local_rng = np.random.default_rng(int(payload.get("seed", 17)))

    centers = local_rng.uniform(0.15, 0.85, size=(clusters, dim))
    counts = np.full(clusters, n_nodes // clusters, dtype=int)
    counts[: n_nodes % clusters] += 1
    assignments = np.repeat(np.arange(clusters), counts)
    local_rng.shuffle(assignments)
    positions = centers[assignments] + local_rng.normal(0.0, spread, size=(n_nodes, dim))
    positions = np.clip(positions, 0.0, 1.0)

    distances = pairwise_distances(positions)
    adjacency = np.exp(-(distances ** 2) / max(2.0 * kernel_width * kernel_width, EPS))
    adjacency[distances > locality_radius] = 0.0
    np.fill_diagonal(adjacency, 0.0)
    return GraphBundle(normalize_adjacency(adjacency), positions, "geometry_emergence/default_geometry_probe_config.json", "haos_config_reconstruction")


def synthetic_ring_graph(node_count: int) -> GraphBundle:
    angles = np.linspace(0.0, 2.0 * math.pi, node_count, endpoint=False)
    positions = np.column_stack([np.cos(angles), np.sin(angles)])
    adjacency = np.zeros((node_count, node_count), dtype=float)
    for idx in range(node_count):
        adjacency[idx, (idx - 1) % node_count] = 1.0
        adjacency[idx, (idx + 1) % node_count] = 1.0
    return GraphBundle(adjacency, positions, "synthetic_ring_fallback", "OPEN_NO_DATA_SYNTHETIC")


def run_single_dynamics(
    label: str,
    adjacency: np.ndarray,
    reference: np.ndarray,
    perturbation_nodes: np.ndarray,
    config: DynamicsConfig,
    rng: np.random.Generator,
    initial_override: np.ndarray | None = None,
) -> DynamicsRun:
    """Evolve one deterministic field trajectory."""

    transition = build_transition(adjacency)
    laplacian = np.eye(adjacency.shape[0], dtype=float) - transition
    initial_seed = reference + 0.05 * rng.normal(size=reference.shape) if initial_override is None else initial_override
    initial = normalize_field(initial_seed, reference)
    states = np.zeros((config.steps + 1, reference.size), dtype=float)
    states[0] = initial
    perturbed = False

    for step in range(1, config.steps + 1):
        prior = states[step - 1].copy()
        if step == config.perturbation_step and not perturbed:
            prior[perturbation_nodes] += config.perturbation_scale * rng.normal(size=perturbation_nodes.size)
            prior = normalize_field(prior, reference)
            perturbed = True
        update = prior - config.dt * config.diffusion * (laplacian @ prior)
        update += config.dt * config.recovery_gain * (reference - update)
        states[step] = normalize_field(update, reference)

    metrics = compute_metrics(label, states, reference, adjacency, perturbation_nodes, config)
    summary = summarize_metrics(label, metrics, config.perturbation_step)
    return DynamicsRun(label, states, metrics, summary)


def build_reference_field(positions: np.ndarray) -> np.ndarray:
    x = positions[:, 0]
    y = positions[:, 1] if positions.shape[1] > 1 else np.zeros_like(x)
    field = np.sin(2.0 * math.pi * x) + 0.55 * np.cos(2.0 * math.pi * y)
    return normalize_field(field, field)


def choose_perturbation_nodes(positions: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    center = positions[int(rng.integers(0, positions.shape[0]))]
    distances = np.linalg.norm(positions - center, axis=1)
    count = max(4, int(round(0.12 * positions.shape[0])))
    return np.argsort(distances)[:count]


def compute_metrics(
    label: str,
    states: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    perturbation_nodes: np.ndarray,
    config: DynamicsConfig,
) -> list[dict[str, float | int | str]]:
    degree = np.maximum(adjacency.sum(axis=1), EPS)
    baseline_energy = dirichlet_energy(reference, adjacency)
    baseline_mean = float(np.mean(reference))
    baseline_norm = float(np.linalg.norm(reference))
    distances = graph_distances_from_sources(adjacency, perturbation_nodes)
    rows: list[dict[str, float | int | str]] = []
    previous_recoverability: float | None = None

    for step, state in enumerate(states):
        distance = float(np.linalg.norm(state - reference) / max(math.sqrt(reference.size), EPS))
        recoverability = float(np.clip(1.0 - distance, 0.0, 1.0))
        energy = dirichlet_energy(state, adjacency)
        invariant_drift = float(
            0.45 * abs(np.mean(state) - baseline_mean)
            + 0.35 * abs(np.linalg.norm(state) - baseline_norm) / max(baseline_norm, EPS)
            + 0.20 * abs(energy - baseline_energy) / max(abs(baseline_energy), EPS)
        )
        causal_spread = causal_spread_radius(state - reference, distances, degree)
        delta = 0.0 if previous_recoverability is None else recoverability - previous_recoverability
        previous_recoverability = recoverability
        rows.append(
            {
                "run": label,
                "step": step,
                "time": float(step * config.dt),
                "recoverability": recoverability,
                "delta_persistence": float(delta),
                "invariant_drift": invariant_drift,
                "dirichlet_energy": float(energy),
                "causal_spread": causal_spread,
            }
        )
    return rows


def summarize_metrics(
    label: str,
    rows: list[dict[str, float | int | str]],
    perturbation_step: int,
) -> dict[str, float | int | bool | str | None]:
    recoverability = np.array([float(row["recoverability"]) for row in rows], dtype=float)
    drift = np.array([float(row["invariant_drift"]) for row in rows], dtype=float)
    deltas = np.array([float(row["delta_persistence"]) for row in rows], dtype=float)
    start = min(max(int(perturbation_step), 0), len(recoverability) - 1)
    recovery_window = recoverability[max(start, len(recoverability) // 2) :]
    k_star = first_sustained_crossing(recoverability, RECOVERABILITY_FLOOR, sustain=3, start=start)
    return {
        "run": label,
        "recoverability_score": float(np.mean(recovery_window)),
        "final_recoverability": float(recoverability[-1]),
        "min_recoverability": float(np.min(recoverability)),
        "max_invariant_drift": float(np.max(drift)),
        "min_delta_persistence": float(np.min(deltas)),
        "k_star_time": None if k_star is None else float(rows[k_star]["time"]),
        "recovered": bool(float(recoverability[-1]) >= RECOVERABILITY_FLOOR),
        "invariants_bounded": bool(float(np.max(drift)) <= INVARIANT_DRIFT_CEILING),
    }


def classify_status(graph: GraphBundle, observed: DynamicsRun, controls: list[DynamicsRun]) -> dict[str, Any]:
    observed_score = float(observed.summary["recoverability_score"])
    best_control = max(float(run.summary["recoverability_score"]) for run in controls)
    control_contrast = observed_score >= best_control + CONTROL_MARGIN
    observed_recovered = bool(observed.summary["recovered"])
    observed_bounded = bool(observed.summary["invariants_bounded"])

    if graph.source_kind == "OPEN_NO_DATA_SYNTHETIC":
        bridge_status = "OPEN_NO_DATA_SYNTHETIC"
        failure_mode = "NO_HAOS_GRAPH_ARTIFACT_FOUND"
    elif observed_recovered and observed_bounded and control_contrast:
        bridge_status = "PASS"
        failure_mode = "RECOVERY_WITH_INVARIANT_CONTROL_CONTRAST"
    elif observed_recovered and observed_bounded:
        bridge_status = "MARGINAL"
        failure_mode = "RECOVERY_WITHOUT_CONTROL_CONTRAST"
    elif not observed_bounded:
        bridge_status = "FAIL"
        failure_mode = "INVARIANT_DRIFT_EXCEEDED"
    else:
        bridge_status = "FAIL"
        failure_mode = "RECOVERABILITY_COLLAPSE"

    return {
        "bridge_status": bridge_status,
        "failure_mode": failure_mode,
        "graph_source": graph.source_label,
        "graph_source_kind": graph.source_kind,
        "nodes": graph.node_count,
        "edges": graph.edge_count,
        "observed_recoverability_score": observed_score,
        "best_control_recoverability_score": best_control,
        "control_contrast": control_contrast,
        "recoverability_floor": RECOVERABILITY_FLOOR,
        "invariant_drift_ceiling": INVARIANT_DRIFT_CEILING,
        "observed_summary": observed.summary,
        "control_summaries": {run.label: run.summary for run in controls},
    }


def write_outputs(output_dir: Path, graph: GraphBundle, observed: DynamicsRun, controls: list[DynamicsRun], status: dict[str, Any]) -> None:
    write_json(output_dir / "bridge_status.json", status)
    write_timeseries(output_dir / "dynamics_timeseries.csv", [observed, *controls])
    write_report(output_dir / "dynamics_report.md", status, observed, controls)
    write_plots(output_dir, observed, controls)


def write_report(path: Path, status: dict[str, Any], observed: DynamicsRun, controls: list[DynamicsRun]) -> None:
    lines = [
        "# HAOS Core Dynamics Report",
        "",
        f"- bridge_status: {status['bridge_status']}",
        f"- failure_mode: {status['failure_mode']}",
        f"- graph_source: {status['graph_source']}",
        f"- graph_source_kind: {status['graph_source_kind']}",
        f"- nodes: {status['nodes']}",
        f"- edges: {status['edges']}",
        f"- control_contrast: {status['control_contrast']}",
        "",
        "| run | recoverability_score | final_recoverability | min_recoverability | max_invariant_drift | k_star_time | recovered | invariants_bounded |",
        "| --- | ---: | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    for run in [observed, *controls]:
        summary = run.summary
        lines.append(
            "| {run} | {score:.6f} | {final:.6f} | {minimum:.6f} | {drift:.6f} | {k_star} | {recovered} | {bounded} |".format(
                run=run.label,
                score=float(summary["recoverability_score"]),
                final=float(summary["final_recoverability"]),
                minimum=float(summary["min_recoverability"]),
                drift=float(summary["max_invariant_drift"]),
                k_star=format_optional(summary["k_star_time"]),
                recovered=summary["recovered"],
                bounded=summary["invariants_bounded"],
            )
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "This is scalar-field dynamics on a HAOS-derived graph, not a full cochain-Laplacian hierarchy claim.",
            "A PASS requires observed recovery, bounded invariant drift, and control contrast.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_timeseries(path: Path, runs: list[DynamicsRun]) -> None:
    fields = ["run", "step", "time", "recoverability", "delta_persistence", "invariant_drift", "dirichlet_energy", "causal_spread"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for run in runs:
            writer.writerows(run.metrics)


def write_plots(output_dir: Path, observed: DynamicsRun, controls: list[DynamicsRun]) -> None:
    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 4))
    for run in [observed, *controls]:
        plt.plot([float(row["time"]) for row in run.metrics], [float(row["recoverability"]) for row in run.metrics], label=run.label)
    plt.axhline(RECOVERABILITY_FLOOR, color="firebrick", linestyle="--", linewidth=1)
    plt.xlabel("time")
    plt.ylabel("recoverability")
    plt.title("HAOS dynamics recoverability")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "recoverability_timeseries.png", dpi=160)
    plt.close()

    plt.figure(figsize=(8, 4))
    for run in [observed, *controls]:
        plt.plot([float(row["time"]) for row in run.metrics], [float(row["invariant_drift"]) for row in run.metrics], label=run.label)
    plt.axhline(INVARIANT_DRIFT_CEILING, color="firebrick", linestyle="--", linewidth=1)
    plt.xlabel("time")
    plt.ylabel("invariant drift")
    plt.title("Invariant drift")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "invariant_drift.png", dpi=160)
    plt.close()

    labels = [run.label for run in [observed, *controls]]
    scores = [float(run.summary["recoverability_score"]) for run in [observed, *controls]]
    plt.figure(figsize=(8, 4))
    plt.bar(labels, scores)
    plt.xticks(rotation=20, ha="right")
    plt.ylabel("recoverability score")
    plt.title("Observed dynamics vs controls")
    plt.tight_layout()
    plt.savefig(output_dir / "control_comparison.png", dpi=160)
    plt.close()


def normalize_adjacency(adjacency: np.ndarray) -> np.ndarray:
    arr = np.asarray(adjacency, dtype=float)
    arr = 0.5 * (arr + arr.T)
    np.fill_diagonal(arr, 0.0)
    arr[arr < 0.0] = 0.0
    positive = arr[arr > 0.0]
    if positive.size:
        arr = arr / max(float(np.median(positive)), EPS)
    return arr


def build_transition(adjacency: np.ndarray) -> np.ndarray:
    row_sums = adjacency.sum(axis=1, keepdims=True)
    transition = np.divide(adjacency, np.maximum(row_sums, EPS), out=np.zeros_like(adjacency), where=row_sums > EPS)
    isolated = np.flatnonzero(row_sums[:, 0] <= EPS)
    transition[isolated, isolated] = 1.0
    return transition


def normalize_field(field: np.ndarray, reference: np.ndarray) -> np.ndarray:
    centered = np.asarray(field, dtype=float) - float(np.mean(field))
    target_norm = max(float(np.linalg.norm(reference - float(np.mean(reference)))), EPS)
    norm = max(float(np.linalg.norm(centered)), EPS)
    return centered * (target_norm / norm)


def dirichlet_energy(field: np.ndarray, adjacency: np.ndarray) -> float:
    diff = field[:, None] - field[None, :]
    return float(0.5 * np.sum(adjacency * diff * diff))


def graph_distances_from_sources(adjacency: np.ndarray, sources: np.ndarray) -> np.ndarray:
    distances = np.full(adjacency.shape[0], np.inf, dtype=float)
    distances[sources] = 0.0
    frontier = [int(node) for node in sources.tolist()]
    binary = adjacency > 0.0
    while frontier:
        node = frontier.pop(0)
        for neighbor in np.flatnonzero(binary[node]).tolist():
            if not np.isfinite(distances[neighbor]):
                distances[neighbor] = distances[node] + 1.0
                frontier.append(int(neighbor))
    finite = distances[np.isfinite(distances)]
    fill = float(np.max(finite) + 1.0) if finite.size else 0.0
    distances[~np.isfinite(distances)] = fill
    return distances


def causal_spread_radius(residual: np.ndarray, distances: np.ndarray, degree: np.ndarray) -> float:
    weights = np.abs(residual) * np.maximum(degree, 1.0)
    total = float(np.sum(weights))
    if total <= EPS:
        return 0.0
    return float(np.sum(weights * distances) / total)


def damage_graph(graph: GraphBundle, fraction: float, rng: np.random.Generator) -> GraphBundle:
    if fraction <= 0.0:
        return graph
    adjacency = graph.adjacency.copy()
    edge_i, edge_j = np.triu_indices(adjacency.shape[0], k=1)
    mask = adjacency[edge_i, edge_j] > 0.0
    edges = np.column_stack([edge_i[mask], edge_j[mask]])
    remove_count = int(round(float(fraction) * len(edges)))
    if remove_count > 0:
        selected = rng.choice(len(edges), size=min(remove_count, len(edges)), replace=False)
        for idx in selected:
            i, j = edges[int(idx)]
            adjacency[i, j] = 0.0
            adjacency[j, i] = 0.0
    return GraphBundle(adjacency, graph.positions, f"{graph.source_label}:damage={fraction:.3f}", graph.source_kind)


def edge_shuffle(adjacency: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    upper = adjacency[np.triu_indices(adjacency.shape[0], k=1)]
    shuffled = rng.permutation(upper)
    out = np.zeros_like(adjacency)
    out[np.triu_indices(adjacency.shape[0], k=1)] = shuffled
    return normalize_adjacency(out + out.T)


def topology_randomize(adjacency: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    n_nodes = adjacency.shape[0]
    edge_count = int(np.count_nonzero(np.triu(adjacency, k=1)))
    positive = adjacency[np.triu_indices(n_nodes, k=1)]
    weights = positive[positive > 0.0]
    pairs = [(i, j) for i in range(n_nodes) for j in range(i + 1, n_nodes)]
    selected = rng.choice(len(pairs), size=min(edge_count, len(pairs)), replace=False)
    sampled = rng.choice(weights, size=len(selected), replace=True) if weights.size else np.ones(len(selected))
    out = np.zeros_like(adjacency)
    for pair_idx, weight in zip(selected, sampled):
        i, j = pairs[int(pair_idx)]
        out[i, j] = float(weight)
        out[j, i] = float(weight)
    return normalize_adjacency(out)


def first_sustained_crossing(values: np.ndarray, threshold: float, sustain: int, start: int = 0) -> int | None:
    for idx in range(max(int(start), 0), len(values) - sustain):
        if bool(np.all(values[idx : idx + sustain + 1] >= threshold)):
            return idx
    return None


def pairwise_distances(points: np.ndarray) -> np.ndarray:
    delta = points[:, None, :] - points[None, :, :]
    return np.sqrt(np.sum(delta * delta, axis=2))


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(to_jsonable(data), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def to_jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(item) for item in value]
    return value


def format_optional(value: Any) -> str:
    if value is None:
        return ""
    return f"{float(value):.6g}"
