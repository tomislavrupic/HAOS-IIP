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
    invariant_gain: float = 0.025
    perturbation_step: int = 20
    perturbation_scale: float = 0.85
    permutation_trials: int = 64
    identity_bins: int = 4
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
        run_single("degree_preserving_rewire_control", degree_preserving_rewire(graph.adjacency, rng), reference, perturbation_nodes, config, rng),
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
    target_shell_variance = local_shell_variance(reference, adjacency)
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
        invariant_pull = shell_variance_restoration(prior, adjacency, target_shell_variance)
        update = prior - config.dt * config.diffusion * (laplacian @ prior)
        update += config.dt * config.address_gain * address_error / np.maximum(adjacency.sum(axis=1), 1.0)
        update += config.dt * config.invariant_gain * invariant_pull
        states[step] = normalize_field(update, reference)

    rows = compute_rows(label, states, reference, adjacency, config)
    summary = summarize_run(label, rows, states[-1], reference, adjacency, config)
    return MinimalRun(label, states, rows, summary)


def local_address(field: np.ndarray, adjacency: np.ndarray) -> np.ndarray:
    """Weighted neighbor-difference signature for each node."""

    degree = np.maximum(adjacency.sum(axis=1), EPS)
    return (adjacency @ field) / degree - field


def local_shell_variance(field: np.ndarray, adjacency: np.ndarray) -> np.ndarray:
    """Frozen branch-local contrast energy around each node.

    This is deliberately small-scope: for each node, measure the weighted
    variance of neighbor offsets relative to that node. It gives the dynamics a
    local invariant target that is not just "return to the reference value".
    """

    degree = np.maximum(adjacency.sum(axis=1), EPS)
    offsets = field[None, :] - field[:, None]
    return np.sum(adjacency * offsets * offsets, axis=1) / degree


def shell_variance_restoration(
    field: np.ndarray,
    adjacency: np.ndarray,
    target_shell_variance: np.ndarray,
) -> np.ndarray:
    """Return a bounded local pull toward the frozen shell-variance invariant."""

    neighbor_mean = field + local_address(field, adjacency)
    current = local_shell_variance(field, adjacency)
    error = target_shell_variance - current
    scale = np.maximum(np.abs(target_shell_variance) + np.abs(current), 1.0)
    direction = field - neighbor_mean
    return np.clip(error / scale, -1.0, 1.0) * direction


def address_retention(field: np.ndarray, reference: np.ndarray, adjacency: np.ndarray) -> float:
    score = address_error_score(field, reference, adjacency)
    return float(np.clip(1.0 + score, 0.0, 1.0))


def address_error_score(field: np.ndarray, reference: np.ndarray, adjacency: np.ndarray) -> float:
    ref = local_address(reference, adjacency)
    cur = local_address(field, adjacency)
    distance = float(np.linalg.norm(cur - ref) / max(np.linalg.norm(ref), EPS))
    return float(-distance)


def invariant_retention(field: np.ndarray, reference: np.ndarray, adjacency: np.ndarray) -> float:
    score = invariant_error_score(field, reference, adjacency)
    return float(np.clip(1.0 + score, 0.0, 1.0))


def invariant_error_score(field: np.ndarray, reference: np.ndarray, adjacency: np.ndarray) -> float:
    ref = local_shell_variance(reference, adjacency)
    cur = local_shell_variance(field, adjacency)
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


def invariant_specificity_test(
    field: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    trials: int,
    seed: int,
) -> dict[str, float | bool]:
    observed = invariant_error_score(field, reference, adjacency)
    rng = np.random.default_rng(seed)
    nulls = np.zeros(trials, dtype=float)
    for idx in range(trials):
        nulls[idx] = invariant_error_score(rng.permutation(field), reference, adjacency)
    mean_null = float(np.mean(nulls))
    std_null = float(np.std(nulls))
    z_score = float((observed - mean_null) / (std_null if std_null > 1.0e-8 else 1.0))
    p_value = float((np.count_nonzero(nulls >= observed) + 1) / (trials + 1))
    return {
        "invariant_specificity": observed,
        "invariant_retention": invariant_retention(field, reference, adjacency),
        "invariant_specificity_null_mean": mean_null,
        "invariant_specificity_p": p_value,
        "invariant_specificity_z": z_score,
        "invariant_specificity_pass": bool(p_value < 0.05 and z_score > 1.5),
    }


def branch_identity_specificity_test(
    field: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    trials: int,
    seed: int,
    bins: int,
) -> dict[str, float | bool | int]:
    """Specificity test with degree/shell-stat preserving node shuffles.

    The global permutation tests ask whether the final field keeps any relational
    structure. This stricter null asks whether the signal survives after
    preserving two cheap local summaries: weighted degree and frozen shell
    variance. Only values inside matched degree/shell buckets are permuted, so
    an observed win is harder to explain as generic density or contrast effects.
    """

    observed = branch_identity_score(field, reference, adjacency)
    rng = np.random.default_rng(seed)
    nulls = np.zeros(trials, dtype=float)
    for idx in range(trials):
        shuffled = degree_shell_stratified_permutation(field, reference, adjacency, rng, bins)
        nulls[idx] = branch_identity_score(shuffled, reference, adjacency)
    mean_null = float(np.mean(nulls))
    std_null = float(np.std(nulls))
    z_score = float((observed - mean_null) / (std_null if std_null > 1.0e-8 else 1.0))
    p_value = float((np.count_nonzero(nulls >= observed) + 1) / (trials + 1))
    return {
        "branch_identity_specificity": observed,
        "branch_identity_null_mean": mean_null,
        "branch_identity_p": p_value,
        "branch_identity_z": z_score,
        "branch_identity_bins": int(bins),
        "branch_identity_pass": bool(p_value < 0.05 and z_score > 1.5),
    }


def branch_identity_score(field: np.ndarray, reference: np.ndarray, adjacency: np.ndarray) -> float:
    """Combined address + invariant score used by the stricter identity null."""

    return float(0.5 * (address_error_score(field, reference, adjacency) + invariant_error_score(field, reference, adjacency)))


def degree_shell_stratified_permutation(
    field: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    rng: np.random.Generator,
    bins: int,
) -> np.ndarray:
    """Shuffle field values within matched weighted-degree/shell-variance buckets."""

    labels = degree_shell_strata(reference, adjacency, bins)
    out = np.asarray(field, dtype=float).copy()
    moved = False
    for label in np.unique(labels):
        idx = np.flatnonzero(labels == label)
        if idx.size < 2:
            continue
        permuted_idx = rng.permutation(idx)
        if np.array_equal(permuted_idx, idx):
            permuted_idx = np.roll(idx, 1)
        out[idx] = field[permuted_idx]
        moved = True
    return out if moved else rng.permutation(field)


def degree_shell_strata(reference: np.ndarray, adjacency: np.ndarray, bins: int) -> np.ndarray:
    safe_bins = max(2, int(bins))
    degree = adjacency.sum(axis=1)
    shell = local_shell_variance(reference, adjacency)
    return quantile_bins(degree, safe_bins) * safe_bins + quantile_bins(shell, safe_bins)


def quantile_bins(values: np.ndarray, bins: int) -> np.ndarray:
    if values.size == 0:
        return np.zeros(0, dtype=int)
    edges = np.quantile(values, np.linspace(0.0, 1.0, bins + 1)[1:-1])
    return np.digitize(values, edges, right=True).astype(int)


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
        address_score = address_retention(state, reference, adjacency)
        invariant_score = invariant_retention(state, reference, adjacency)
        delta = 0.0 if prev is None else recoverability - prev
        prev = recoverability
        rows.append(
            {
                "run": label,
                "step": step,
                "time": float(step * config.dt),
                "recoverability": recoverability,
                "delta_persistence": float(delta),
                "address_retention": address_score,
                "invariant_retention": invariant_score,
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
    invariant_values = np.array([float(row["invariant_retention"]) for row in rows])
    start = min(config.perturbation_step, len(rows) - 1)
    recovery_window = recoverability[max(start, len(rows) // 2) :]
    specificity = address_specificity_test(final_state, reference, adjacency, config.permutation_trials, config.seed + stable_label_offset(label))
    invariant_specificity = invariant_specificity_test(
        final_state,
        reference,
        adjacency,
        config.permutation_trials,
        config.seed + stable_label_offset(label) + 10_000,
    )
    branch_identity_specificity = branch_identity_specificity_test(
        final_state,
        reference,
        adjacency,
        config.permutation_trials,
        config.seed + stable_label_offset(label) + 20_000,
        config.identity_bins,
    )
    combined_specificity_pass = bool(specificity["address_specificity_pass"] and invariant_specificity["invariant_specificity_pass"])
    strict_specificity_pass = bool(combined_specificity_pass and branch_identity_specificity["branch_identity_pass"])
    return {
        "run": label,
        "recoverability_score": float(np.mean(recovery_window)),
        "final_recoverability": float(recoverability[-1]),
        "min_recoverability": float(np.min(recoverability)),
        "final_address_retention": float(address_values[-1]),
        "min_address_retention": float(np.min(address_values)),
        "final_invariant_retention": float(invariant_values[-1]),
        "min_invariant_retention": float(np.min(invariant_values)),
        "recovered": bool(float(recoverability[-1]) >= RECOVERABILITY_FLOOR),
        **specificity,
        **invariant_specificity,
        **branch_identity_specificity,
        "combined_specificity_pass": combined_specificity_pass,
        "strict_specificity_pass": strict_specificity_pass,
    }


def classify(graph: GraphData, observed: MinimalRun, controls: list[MinimalRun]) -> dict[str, Any]:
    best_control = max(float(run.summary["recoverability_score"]) for run in controls)
    control_contrast = float(observed.summary["recoverability_score"]) >= best_control + CONTROL_MARGIN
    control_specificity_passes = sum(1 for run in controls if bool(run.summary["combined_specificity_pass"]))
    control_strict_specificity_passes = sum(1 for run in controls if bool(run.summary["strict_specificity_pass"]))
    observed_specific = bool(observed.summary["strict_specificity_pass"])
    observed_recovered = bool(observed.summary["recovered"])

    if graph.source_kind == "OPEN_NO_DATA_SYNTHETIC":
        status = "OPEN_NO_DATA_SYNTHETIC"
        failure = "NO_HAOS_GRAPH_ARTIFACT_FOUND"
    elif observed_recovered and observed_specific and control_contrast and control_strict_specificity_passes == 0:
        status = "PASS"
        failure = "RECOVERY_WITH_STRICT_BRANCH_IDENTITY_SPECIFICITY"
    elif observed_recovered and observed_specific and control_strict_specificity_passes > 0:
        status = "MARGINAL"
        failure = "STRICT_SPECIFICITY_CONTROL_MATCH"
    elif observed_recovered and observed_specific:
        status = "MARGINAL"
        failure = "STRICT_SIGNAL_WITHOUT_CONTROL_CONTRAST"
    elif not observed_recovered:
        status = "FAIL"
        failure = "RECOVERABILITY_COLLAPSE"
    else:
        status = "FAIL"
        failure = "STRICT_BRANCH_IDENTITY_SPECIFICITY_FAILED"

    return {
        "bridge_status": status,
        "failure_mode": failure,
        "graph_source": graph.source_label,
        "graph_source_kind": graph.source_kind,
        "nodes": int(graph.adjacency.shape[0]),
        "edges": int(np.count_nonzero(np.triu(graph.adjacency, k=1))),
        "control_contrast": control_contrast,
        "control_combined_specificity_pass_count": control_specificity_passes,
        "control_strict_specificity_pass_count": control_strict_specificity_passes,
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
        f"- control_combined_specificity_pass_count: {status['control_combined_specificity_pass_count']}",
        f"- control_strict_specificity_pass_count: {status['control_strict_specificity_pass_count']}",
        "",
        "| run | recoverability_score | final_recoverability | address_retention | invariant_retention | address_z | invariant_z | branch_z | strict_pass |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for run in [observed, *controls]:
        summary = run.summary
        lines.append(
            "| {run} | {score:.6f} | {final:.6f} | {address:.6f} | {invariant:.6f} | {address_z:.6f} | {invariant_z:.6f} | {branch_z:.6f} | {passed} |".format(
                run=run.label,
                score=float(summary["recoverability_score"]),
                final=float(summary["final_recoverability"]),
                address=float(summary["final_address_retention"]),
                invariant=float(summary["final_invariant_retention"]),
                address_z=float(summary["address_specificity_z"]),
                invariant_z=float(summary["invariant_specificity_z"]),
                branch_z=float(summary["branch_identity_z"]),
                passed=summary["strict_specificity_pass"],
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_csv(path: Path, runs: list[MinimalRun]) -> None:
    fields = ["run", "step", "time", "recoverability", "delta_persistence", "address_retention", "invariant_retention"]
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
        ("invariant_retention", "invariant_retention.png", "invariant retention"),
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


def degree_preserving_rewire(adjacency: np.ndarray, rng: np.random.Generator, swap_factor: int = 8) -> np.ndarray:
    """Randomize topology with double-edge swaps while preserving unweighted degrees."""

    binary = np.asarray(adjacency > 0.0, dtype=bool)
    n = binary.shape[0]
    edge_i, edge_j = np.triu_indices(n, k=1)
    edges = [(int(i), int(j)) for i, j in zip(edge_i, edge_j) if binary[i, j]]
    edge_set = set(edges)
    if len(edges) < 2:
        return adjacency.copy()

    attempts = max(len(edges) * swap_factor * 4, 1)
    target_swaps = len(edges) * swap_factor
    swaps = 0
    for _ in range(attempts):
        if swaps >= target_swaps:
            break
        idx_a, idx_b = rng.choice(len(edges), size=2, replace=False)
        a, b = edges[int(idx_a)]
        c, d = edges[int(idx_b)]
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
        if new_a in edge_set or new_b in edge_set:
            continue
        old_a = edges[int(idx_a)]
        old_b = edges[int(idx_b)]
        edge_set.remove(old_a)
        edge_set.remove(old_b)
        edge_set.add(new_a)
        edge_set.add(new_b)
        edges[int(idx_a)] = new_a
        edges[int(idx_b)] = new_b
        swaps += 1

    weights = adjacency[np.triu_indices(n, k=1)]
    weights = weights[weights > 0.0]
    shuffled_weights = rng.permutation(weights)
    out = np.zeros_like(adjacency, dtype=float)
    for (i, j), weight in zip(edges, shuffled_weights):
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
