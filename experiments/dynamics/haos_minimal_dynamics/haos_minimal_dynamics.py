"""Stripped-back HAOS relational-address dynamics probe."""

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
CONTROL_MARGIN = 0.05


@dataclass(frozen=True)
class MinimalConfig:
    output_dir: Path = ROOT / "outputs"
    seed: int = 20260425
    steps: int = 120
    dt: float = 0.08
    diffusion: float = 0.15
    address_gain: float = 0.05
    perturbation_step: int = 20
    perturbation_scale: float = 0.85
    permutation_trials: int = 64
    max_nodes: int = 96


@dataclass(frozen=True)
class GraphData:
    adjacency: np.ndarray
    positions: np.ndarray
    source_label: str
    source_kind: str


@dataclass(frozen=True)
class MinimalRun:
    label: str
    states: np.ndarray
    rows: list[dict[str, float | int | str]]
    summary: dict[str, Any]


def run_minimal_probe(config: MinimalConfig) -> dict[str, Any]:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(config.seed)
    graph = load_geometry_graph(config)
    reference = build_reference_field(graph.positions)
    perturbation_nodes = choose_perturbation_nodes(graph.positions, rng)

    observed = run_single("observed", graph.adjacency, reference, perturbation_nodes, config, rng)
    controls = [
        run_single("node_shuffle_control", graph.adjacency, reference, perturbation_nodes, config, rng, initial_override=rng.permutation(reference)),
        run_single("edge_shuffle_control", edge_shuffle(graph.adjacency, rng), reference, perturbation_nodes, config, rng),
        run_single("topology_randomization_control", topology_randomize(graph.adjacency, rng), reference, perturbation_nodes, config, rng),
    ]
    status = classify(graph, observed, controls)
    write_outputs(config.output_dir, status, observed, controls)
    return status


def run_single(
    label: str,
    adjacency: np.ndarray,
    reference: np.ndarray,
    perturbation_nodes: np.ndarray,
    config: MinimalConfig,
    rng: np.random.Generator,
    initial_override: np.ndarray | None = None,
) -> MinimalRun:
    transition = build_transition(adjacency)
    laplacian = np.eye(adjacency.shape[0]) - transition
    target_address = local_address(reference, adjacency)
    initial = reference + 0.04 * rng.normal(size=reference.shape) if initial_override is None else initial_override
    states = np.zeros((config.steps + 1, reference.size), dtype=float)
    states[0] = normalize_field(initial, reference)

    for step in range(1, config.steps + 1):
        prior = states[step - 1].copy()
        if step == config.perturbation_step:
            prior[perturbation_nodes] += config.perturbation_scale * rng.normal(size=perturbation_nodes.size)
            prior = normalize_field(prior, reference)
        current_address = local_address(prior, adjacency)
        address_error = target_address - current_address
        update = prior - config.dt * config.diffusion * (laplacian @ prior)
        update += config.dt * config.address_gain * address_error / np.maximum(adjacency.sum(axis=1), 1.0)
        states[step] = normalize_field(update, reference)

    rows = compute_rows(label, states, reference, adjacency, config)
    summary = summarize_run(label, rows, states[-1], reference, adjacency, config)
    return MinimalRun(label, states, rows, summary)


def local_address(field: np.ndarray, adjacency: np.ndarray) -> np.ndarray:
    """Weighted neighbor-difference signature for each node."""

    degree = np.maximum(adjacency.sum(axis=1), EPS)
    return (adjacency @ field) / degree - field


def address_retention(field: np.ndarray, reference: np.ndarray, adjacency: np.ndarray) -> float:
    score = address_error_score(field, reference, adjacency)
    return float(np.clip(1.0 + score, 0.0, 1.0))


def address_error_score(field: np.ndarray, reference: np.ndarray, adjacency: np.ndarray) -> float:
    ref = local_address(reference, adjacency)
    cur = local_address(field, adjacency)
    distance = float(np.linalg.norm(cur - ref) / max(np.linalg.norm(ref), EPS))
    return float(-distance)


def address_specificity_test(
    field: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    trials: int,
    seed: int,
) -> dict[str, float | bool]:
    observed = address_error_score(field, reference, adjacency)
    rng = np.random.default_rng(seed)
    nulls = np.zeros(trials, dtype=float)
    for idx in range(trials):
        nulls[idx] = address_error_score(rng.permutation(field), reference, adjacency)
    mean_null = float(np.mean(nulls))
    std_null = float(np.std(nulls))
    z_score = float((observed - mean_null) / (std_null if std_null > 1.0e-8 else 1.0))
    p_value = float((np.count_nonzero(nulls >= observed) + 1) / (trials + 1))
    return {
        "address_specificity": observed,
        "address_retention": address_retention(field, reference, adjacency),
        "address_specificity_null_mean": mean_null,
        "address_specificity_p": p_value,
        "address_specificity_z": z_score,
        "address_specificity_pass": bool(p_value < 0.05 and z_score > 1.5),
    }


def compute_rows(
    label: str,
    states: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    config: MinimalConfig,
) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    prev = None
    for step, state in enumerate(states):
        recoverability = state_recoverability(state, reference)
        specificity = address_retention(state, reference, adjacency)
        delta = 0.0 if prev is None else recoverability - prev
        prev = recoverability
        rows.append(
            {
                "run": label,
                "step": step,
                "time": float(step * config.dt),
                "recoverability": recoverability,
                "delta_persistence": float(delta),
                "address_retention": specificity,
            }
        )
    return rows


def summarize_run(
    label: str,
    rows: list[dict[str, float | int | str]],
    final_state: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    config: MinimalConfig,
) -> dict[str, Any]:
    recoverability = np.array([float(row["recoverability"]) for row in rows])
    address_values = np.array([float(row["address_retention"]) for row in rows])
    start = min(config.perturbation_step, len(rows) - 1)
    recovery_window = recoverability[max(start, len(rows) // 2) :]
    specificity = address_specificity_test(final_state, reference, adjacency, config.permutation_trials, config.seed + stable_label_offset(label))
    return {
        "run": label,
        "recoverability_score": float(np.mean(recovery_window)),
        "final_recoverability": float(recoverability[-1]),
        "min_recoverability": float(np.min(recoverability)),
        "final_address_retention": float(address_values[-1]),
        "min_address_retention": float(np.min(address_values)),
        "recovered": bool(float(recoverability[-1]) >= RECOVERABILITY_FLOOR),
        **specificity,
    }


def classify(graph: GraphData, observed: MinimalRun, controls: list[MinimalRun]) -> dict[str, Any]:
    best_control = max(float(run.summary["recoverability_score"]) for run in controls)
    control_contrast = float(observed.summary["recoverability_score"]) >= best_control + CONTROL_MARGIN
    control_specificity_passes = sum(1 for run in controls if bool(run.summary["address_specificity_pass"]))
    observed_specific = bool(observed.summary["address_specificity_pass"])
    observed_recovered = bool(observed.summary["recovered"])

    if graph.source_kind == "OPEN_NO_DATA_SYNTHETIC":
        status = "OPEN_NO_DATA_SYNTHETIC"
        failure = "NO_HAOS_GRAPH_ARTIFACT_FOUND"
    elif observed_recovered and observed_specific and control_contrast and control_specificity_passes == 0:
        status = "PASS"
        failure = "RECOVERY_WITH_ADDRESS_SPECIFICITY_AND_CONTROL_CONTRAST"
    elif observed_recovered and observed_specific and control_specificity_passes > 0:
        status = "MARGINAL"
        failure = "ADDRESS_SPECIFICITY_CONTROL_MATCH"
    elif observed_recovered and observed_specific:
        status = "MARGINAL"
        failure = "ADDRESS_SIGNAL_WITHOUT_CONTROL_CONTRAST"
    elif not observed_recovered:
        status = "FAIL"
        failure = "RECOVERABILITY_COLLAPSE"
    else:
        status = "FAIL"
        failure = "ADDRESS_SPECIFICITY_FAILED"

    return {
        "bridge_status": status,
        "failure_mode": failure,
        "graph_source": graph.source_label,
        "graph_source_kind": graph.source_kind,
        "nodes": int(graph.adjacency.shape[0]),
        "edges": int(np.count_nonzero(np.triu(graph.adjacency, k=1))),
        "control_contrast": control_contrast,
        "control_specificity_pass_count": control_specificity_passes,
        "best_control_recoverability_score": best_control,
        "observed_summary": observed.summary,
        "control_summaries": {run.label: run.summary for run in controls},
    }


def write_outputs(output_dir: Path, status: dict[str, Any], observed: MinimalRun, controls: list[MinimalRun]) -> None:
    write_json(output_dir / "bridge_status.json", status)
    write_csv(output_dir / "minimal_timeseries.csv", [observed, *controls])
    write_report(output_dir / "minimal_dynamics_report.md", status, observed, controls)
    write_plots(output_dir, observed, controls)


def write_report(path: Path, status: dict[str, Any], observed: MinimalRun, controls: list[MinimalRun]) -> None:
    lines = [
        "# HAOS Minimal Dynamics Report",
        "",
        f"- bridge_status: {status['bridge_status']}",
        f"- failure_mode: {status['failure_mode']}",
        f"- graph_source: {status['graph_source']}",
        f"- control_contrast: {status['control_contrast']}",
        f"- control_specificity_pass_count: {status['control_specificity_pass_count']}",
        "",
        "| run | recoverability_score | final_recoverability | final_address_retention | specificity_z | specificity_p | specificity_pass |",
        "| --- | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for run in [observed, *controls]:
        summary = run.summary
        lines.append(
            "| {run} | {score:.6f} | {final:.6f} | {address:.6f} | {z:.6f} | {p:.6g} | {passed} |".format(
                run=run.label,
                score=float(summary["recoverability_score"]),
                final=float(summary["final_recoverability"]),
                address=float(summary["final_address_retention"]),
                z=float(summary["address_specificity_z"]),
                p=float(summary["address_specificity_p"]),
                passed=summary["address_specificity_pass"],
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_csv(path: Path, runs: list[MinimalRun]) -> None:
    fields = ["run", "step", "time", "recoverability", "delta_persistence", "address_retention"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for run in runs:
            writer.writerows(run.rows)


def write_plots(output_dir: Path, observed: MinimalRun, controls: list[MinimalRun]) -> None:
    import matplotlib.pyplot as plt

    for key, filename, ylabel in (
        ("recoverability", "recoverability.png", "recoverability"),
        ("address_retention", "address_specificity.png", "address retention"),
    ):
        plt.figure(figsize=(8, 4))
        for run in [observed, *controls]:
            plt.plot([float(row["time"]) for row in run.rows], [float(row[key]) for row in run.rows], label=run.label)
        plt.xlabel("time")
        plt.ylabel(ylabel)
        plt.title(f"Minimal dynamics {ylabel}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_dir / filename, dpi=160)
        plt.close()


def load_geometry_graph(config: MinimalConfig) -> GraphData:
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
    radius = float(payload.get("locality_radius_factor", 3.0)) * kernel_width
    rng = np.random.default_rng(int(payload.get("seed", 17)))
    centers = rng.uniform(0.15, 0.85, size=(clusters, dim))
    counts = np.full(clusters, n_nodes // clusters, dtype=int)
    counts[: n_nodes % clusters] += 1
    assignments = np.repeat(np.arange(clusters), counts)
    rng.shuffle(assignments)
    positions = np.clip(centers[assignments] + rng.normal(0.0, spread, size=(n_nodes, dim)), 0.0, 1.0)
    distances = pairwise_distances(positions)
    adjacency = np.exp(-(distances ** 2) / max(2.0 * kernel_width * kernel_width, EPS))
    adjacency[distances > radius] = 0.0
    np.fill_diagonal(adjacency, 0.0)
    return GraphData(normalize_adjacency(adjacency), positions, "geometry_emergence/default_geometry_probe_config.json", "haos_config_reconstruction")


def build_reference_field(positions: np.ndarray) -> np.ndarray:
    x = positions[:, 0]
    y = positions[:, 1] if positions.shape[1] > 1 else np.zeros_like(x)
    field = np.sin(2.0 * math.pi * x) + 0.55 * np.cos(2.0 * math.pi * y)
    return normalize_field(field, field)


def choose_perturbation_nodes(positions: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    center = positions[int(rng.integers(0, positions.shape[0]))]
    distances = np.linalg.norm(positions - center, axis=1)
    return np.argsort(distances)[: max(4, int(round(0.12 * positions.shape[0])))]


def state_recoverability(state: np.ndarray, reference: np.ndarray) -> float:
    distance = float(np.linalg.norm(state - reference) / max(math.sqrt(reference.size), EPS))
    return float(np.clip(1.0 - distance, 0.0, 1.0))


def build_transition(adjacency: np.ndarray) -> np.ndarray:
    row_sums = adjacency.sum(axis=1, keepdims=True)
    transition = np.divide(adjacency, np.maximum(row_sums, EPS), out=np.zeros_like(adjacency), where=row_sums > EPS)
    isolated = np.flatnonzero(row_sums[:, 0] <= EPS)
    transition[isolated, isolated] = 1.0
    return transition


def normalize_field(field: np.ndarray, reference: np.ndarray) -> np.ndarray:
    centered = np.asarray(field, dtype=float) - float(np.mean(field))
    target = np.asarray(reference, dtype=float) - float(np.mean(reference))
    return centered * (max(float(np.linalg.norm(target)), EPS) / max(float(np.linalg.norm(centered)), EPS))


def normalize_adjacency(adjacency: np.ndarray) -> np.ndarray:
    arr = 0.5 * (np.asarray(adjacency, dtype=float) + np.asarray(adjacency, dtype=float).T)
    np.fill_diagonal(arr, 0.0)
    positive = arr[arr > 0.0]
    return arr / max(float(np.median(positive)), EPS) if positive.size else arr


def edge_shuffle(adjacency: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    upper = adjacency[np.triu_indices(adjacency.shape[0], k=1)]
    shuffled = rng.permutation(upper)
    out = np.zeros_like(adjacency)
    out[np.triu_indices(adjacency.shape[0], k=1)] = shuffled
    return normalize_adjacency(out + out.T)


def topology_randomize(adjacency: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    n = adjacency.shape[0]
    edge_count = int(np.count_nonzero(np.triu(adjacency, k=1)))
    weights = adjacency[np.triu_indices(n, k=1)]
    weights = weights[weights > 0.0]
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    selected = rng.choice(len(pairs), size=min(edge_count, len(pairs)), replace=False)
    sampled = rng.choice(weights, size=len(selected), replace=True) if weights.size else np.ones(len(selected))
    out = np.zeros_like(adjacency)
    for pair_idx, weight in zip(selected, sampled):
        i, j = pairs[int(pair_idx)]
        out[i, j] = float(weight)
        out[j, i] = float(weight)
    return normalize_adjacency(out)


def synthetic_ring_graph(node_count: int) -> GraphData:
    angles = np.linspace(0.0, 2.0 * math.pi, node_count, endpoint=False)
    positions = np.column_stack([np.cos(angles), np.sin(angles)])
    adjacency = np.zeros((node_count, node_count), dtype=float)
    for idx in range(node_count):
        adjacency[idx, (idx - 1) % node_count] = 1.0
        adjacency[idx, (idx + 1) % node_count] = 1.0
    return GraphData(adjacency, positions, "synthetic_ring_fallback", "OPEN_NO_DATA_SYNTHETIC")


def pairwise_distances(points: np.ndarray) -> np.ndarray:
    delta = points[:, None, :] - points[None, :, :]
    return np.sqrt(np.sum(delta * delta, axis=2))


def stable_label_offset(label: str) -> int:
    return sum((idx + 1) * ord(char) for idx, char in enumerate(label))


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
