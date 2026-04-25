"""Kuramoto bridge for HAOS-IIP oscillator experiments.

This module is intentionally a sidecar: it reads HAOS-adjacent graph artifacts
when possible, runs deterministic Kuramoto scans, and writes JSON/Markdown-first
diagnostics. It does not import or modify frozen HAOS core APIs.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np


ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parent.parent.parent
STRONG_SYNC_THRESHOLD = 0.82
RECOVERABILITY_THRESHOLD = 0.72
SUSTAIN_STEPS = 2
EPS = 1.0e-12


@dataclass(frozen=True)
class KuramotoConfig:
    """Runtime configuration for the bridge scan."""

    explicit_graph_path: Path | None = None
    output_dir: Path = ROOT / "outputs"
    seed: int = 20260425
    max_nodes: int = 96
    k_min: float = 0.0
    k_max: float = 4.0
    k_count: int = 17
    steps: int = 900
    dt: float = 0.03
    omega_mean: float = 0.0
    omega_std: float = 1.0
    frequency_mode: str = "gaussian"
    frequency_file: Path | None = None
    higher_order: bool = False
    higher_order_strength: float = 0.15


@dataclass(frozen=True)
class GraphData:
    """Undirected weighted graph used by the Kuramoto model."""

    adjacency: np.ndarray
    source_label: str
    source_path: str | None
    source_kind: str

    @property
    def node_count(self) -> int:
        return int(self.adjacency.shape[0])

    @property
    def edge_count(self) -> int:
        return int(np.count_nonzero(np.triu(self.adjacency, k=1)))

    @property
    def mean_degree(self) -> float:
        return float(np.count_nonzero(self.adjacency, axis=1).mean())


@dataclass(frozen=True)
class SimulationResult:
    """One K-scan trajectory and derived time-series proxies."""

    label: str
    k_value: float
    times: np.ndarray
    phases: np.ndarray
    order_parameter: np.ndarray
    recoverability: np.ndarray
    delta_persistence: np.ndarray
    trajectory_identity_retention: np.ndarray
    k_star_time: float | None
    r_star_time: float | None
    early_detection: bool


@dataclass(frozen=True)
class ScanRun:
    """Collection of simulations for an observed or control condition."""

    label: str
    graph: GraphData
    omega: np.ndarray
    simulations: list[SimulationResult]
    summary: dict[str, Any]


@dataclass(frozen=True)
class BridgeOutcome:
    """Top-level bridge output returned to the launcher."""

    status: str
    graph: GraphData
    primary_summary: dict[str, Any]
    output_dir: Path


def run_bridge(config: KuramotoConfig) -> BridgeOutcome:
    """Run the full bridge: load graph, scan K, run controls, write outputs."""

    config.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(config.seed)
    graph = load_haos_graph(config, rng)
    omega = make_frequencies(config, graph.node_count, rng)
    k_values = np.linspace(config.k_min, config.k_max, config.k_count)

    observed = run_scan("observed", graph, omega, k_values, config, rng)
    controls = [
        run_scan("frequency_shuffle_control", graph, rng.permutation(omega), k_values, config, rng),
        run_scan("edge_shuffle_control", edge_shuffle_graph(graph, rng), omega, k_values, config, rng),
        run_scan("topology_randomization_control", randomize_topology_graph(graph, rng), omega, k_values, config, rng),
    ]

    status = classify_status(graph, observed, controls)
    write_outputs(config, status, graph, observed, controls)
    return BridgeOutcome(status=status["bridge_status"], graph=graph, primary_summary=observed.summary, output_dir=config.output_dir)


def load_haos_graph(config: KuramotoConfig, rng: np.random.Generator) -> GraphData:
    """Load an explicit or discoverable HAOS graph, with honest fallback labels."""

    candidates: list[Path] = []
    if config.explicit_graph_path is not None:
        candidates.append(config.explicit_graph_path)
    candidates.extend(discover_candidate_graph_paths())

    for path in candidates:
        try:
            graph = load_graph_path(path, config.max_nodes, rng)
        except Exception:
            continue
        if graph.node_count >= 3 and graph.edge_count > 0:
            return graph

    geometry_config = REPO_ROOT / "geometry_emergence" / "experiments" / "default_geometry_probe_config.json"
    if geometry_config.exists():
        return reconstruct_geometry_emergence_graph(geometry_config, config.max_nodes, rng)

    return synthetic_ring_graph(min(32, config.max_nodes))


def discover_candidate_graph_paths() -> list[Path]:
    """Return likely graph artifact paths without treating metrics as proof."""

    roots = [REPO_ROOT / name for name in ("data", "telemetry", "geometry_emergence", "experiments")]
    suffixes = {".json", ".csv", ".npz", ".npy"}
    paths: list[Path] = []
    for root in roots:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if path.is_file() and path.suffix.lower() in suffixes:
                paths.append(path)
    prioritized = sorted(
        paths,
        key=lambda p: (
            0 if any(word in p.name.lower() for word in ("graph", "laplacian", "adjacency", "kernel")) else 1,
            len(p.parts),
            p.name,
        ),
    )
    return prioritized[:250]


def load_graph_path(path: Path, max_nodes: int, rng: np.random.Generator) -> GraphData:
    suffix = path.suffix.lower()
    if suffix == ".json":
        data = json.loads(path.read_text(encoding="utf-8"))
        adjacency = graph_from_json_like(data)
    elif suffix == ".csv":
        adjacency = graph_from_csv(path)
    elif suffix == ".npy":
        adjacency = np.load(path)
    elif suffix == ".npz":
        archive = np.load(path)
        key = next((name for name in archive.files if name.lower() in {"adjacency", "adj", "a", "laplacian", "l"}), archive.files[0])
        adjacency = archive[key]
    else:
        raise ValueError(f"unsupported graph suffix: {suffix}")

    adjacency = normalize_adjacency(adjacency)
    adjacency = cap_graph_nodes(adjacency, max_nodes, rng)
    return GraphData(adjacency, source_label=path.name, source_path=relative_path(path), source_kind="loaded_artifact")


def graph_from_json_like(data: Any) -> np.ndarray:
    """Extract an adjacency matrix from common JSON graph shapes."""

    if isinstance(data, list):
        arr = np.array(data, dtype=float)
        if arr.ndim == 2:
            return arr
        raise ValueError("JSON list is not a 2D graph")

    if not isinstance(data, dict):
        raise ValueError("JSON graph must be a dict or 2D list")

    for key in ("adjacency", "adjacency_matrix", "adj", "A", "weights", "kernel"):
        if key in data:
            arr = np.array(data[key], dtype=float)
            if arr.ndim == 2:
                return arr

    for key in ("laplacian", "L", "graph_laplacian"):
        if key in data:
            laplacian = np.array(data[key], dtype=float)
            if laplacian.ndim == 2:
                adjacency = -laplacian.copy()
                np.fill_diagonal(adjacency, 0.0)
                adjacency[adjacency < 0.0] = 0.0
                return adjacency

    for key in ("edges", "edge_list", "links"):
        if key in data:
            return adjacency_from_edges(data[key], data.get("nodes"))

    graph_wrappers = ("graph", "operator", "laplacian", "adjacency", "kernel")
    for key, value in data.items():
        if isinstance(value, dict) and any(token in str(key).lower() for token in graph_wrappers):
            try:
                return graph_from_json_like(value)
            except ValueError:
                continue

    raise ValueError("no graph-like JSON field found")


def adjacency_from_edges(edges: Any, nodes: Any = None) -> np.ndarray:
    node_names: list[Any] = list(nodes) if isinstance(nodes, list) else []
    edge_rows: list[tuple[Any, Any, float]] = []
    for edge in edges:
        if isinstance(edge, dict):
            source = edge.get("source", edge.get("u", edge.get("i")))
            target = edge.get("target", edge.get("v", edge.get("j")))
            weight = float(edge.get("weight", edge.get("w", 1.0)))
        elif isinstance(edge, (list, tuple)) and len(edge) >= 2:
            source, target = edge[0], edge[1]
            weight = float(edge[2]) if len(edge) >= 3 else 1.0
        else:
            continue
        node_names.extend([source, target])
        edge_rows.append((source, target, weight))

    unique_nodes = list(dict.fromkeys(node_names))
    node_index = {node: idx for idx, node in enumerate(unique_nodes)}
    adjacency = np.zeros((len(unique_nodes), len(unique_nodes)), dtype=float)
    for source, target, weight in edge_rows:
        i = node_index[source]
        j = node_index[target]
        if i != j:
            adjacency[i, j] = max(adjacency[i, j], weight)
            adjacency[j, i] = max(adjacency[j, i], weight)
    return adjacency


def graph_from_csv(path: Path) -> np.ndarray:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames and {"source", "target"}.issubset({name.lower() for name in reader.fieldnames}):
            rows = list(reader)
            lower_names = {name.lower(): name for name in reader.fieldnames}
            edges = [
                {
                    "source": row[lower_names["source"]],
                    "target": row[lower_names["target"]],
                    "weight": row.get(lower_names.get("weight", ""), 1.0),
                }
                for row in rows
            ]
            return adjacency_from_edges(edges)

    numeric_rows: list[list[float]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader_plain = csv.reader(handle)
        for row in reader_plain:
            try:
                numeric_rows.append([float(cell) for cell in row])
            except ValueError:
                continue
    arr = np.array(numeric_rows, dtype=float)
    if arr.ndim == 2 and arr.shape[0] == arr.shape[1] and arr.shape[0] >= 3:
        return arr
    raise ValueError("CSV is neither edge list nor square adjacency")


def normalize_adjacency(adjacency: np.ndarray) -> np.ndarray:
    arr = np.asarray(adjacency, dtype=float)
    if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
        raise ValueError("adjacency must be square")
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    arr = 0.5 * (arr + arr.T)
    np.fill_diagonal(arr, 0.0)
    arr[arr < 0.0] = 0.0
    positive = arr[arr > 0.0]
    if positive.size == 0:
        raise ValueError("empty adjacency")
    scale = float(np.median(positive))
    if scale > EPS:
        arr = arr / scale
    return arr


def cap_graph_nodes(adjacency: np.ndarray, max_nodes: int, rng: np.random.Generator) -> np.ndarray:
    if adjacency.shape[0] <= max_nodes:
        return adjacency
    degrees = np.count_nonzero(adjacency, axis=1)
    keep = np.argsort(-degrees)[:max_nodes]
    keep = np.sort(keep)
    capped = adjacency[np.ix_(keep, keep)]
    if np.count_nonzero(capped) == 0:
        keep = np.sort(rng.choice(adjacency.shape[0], size=max_nodes, replace=False))
        capped = adjacency[np.ix_(keep, keep)]
    return capped


def reconstruct_geometry_emergence_graph(path: Path, max_nodes: int, rng: np.random.Generator) -> GraphData:
    """Rebuild the committed geometry-emergence graph from its public config."""

    config = json.loads(path.read_text(encoding="utf-8"))
    n_nodes = min(int(config.get("n_nodes", 72)), max_nodes)
    dim = int(config.get("embedding_dim", 2))
    clusters = int(config.get("cluster_count", 3))
    spread = float(config.get("cluster_spread", 0.075))
    width_values = config.get("kernel_widths", [0.13])
    kernel_width = float(width_values[len(width_values) // 2])
    k_nearest = int(config.get("k_nearest", 6))
    local_rng = np.random.default_rng(int(config.get("seed", 17)))

    centers = local_rng.uniform(0.15, 0.85, size=(clusters, dim))
    points = np.zeros((n_nodes, dim), dtype=float)
    for idx in range(n_nodes):
        center = centers[idx % clusters]
        points[idx] = center + local_rng.normal(0.0, spread, size=dim)

    distances = pairwise_distances(points)
    weights = np.exp(-(distances ** 2) / max(2.0 * kernel_width * kernel_width, EPS))
    np.fill_diagonal(weights, 0.0)
    adjacency = np.zeros_like(weights)
    for idx in range(n_nodes):
        neighbors = np.argsort(distances[idx])[1 : k_nearest + 1]
        adjacency[idx, neighbors] = weights[idx, neighbors]
    adjacency = normalize_adjacency(np.maximum(adjacency, adjacency.T))
    return GraphData(adjacency, "geometry_emergence/default_geometry_probe_config.json", relative_path(path), "haos_config_reconstruction")


def synthetic_ring_graph(node_count: int) -> GraphData:
    adjacency = np.zeros((node_count, node_count), dtype=float)
    for idx in range(node_count):
        adjacency[idx, (idx - 1) % node_count] = 1.0
        adjacency[idx, (idx + 1) % node_count] = 1.0
    return GraphData(adjacency, "synthetic_ring_fallback", None, "OPEN_NO_DATA_SYNTHETIC")


def make_frequencies(config: KuramotoConfig, node_count: int, rng: np.random.Generator) -> np.ndarray:
    if config.frequency_mode == "custom":
        if config.frequency_file is None:
            raise ValueError("--frequency-file is required for custom frequencies")
        values = read_frequency_file(config.frequency_file)
        if len(values) < node_count:
            raise ValueError("frequency file has fewer values than graph nodes")
        omega = np.array(values[:node_count], dtype=float)
    else:
        omega = rng.normal(config.omega_mean, config.omega_std, size=node_count)
    omega = omega - float(np.mean(omega))
    return omega


def read_frequency_file(path: Path) -> list[float]:
    values: list[float] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        for token in raw.replace(",", " ").split():
            try:
                values.append(float(token))
            except ValueError:
                continue
    return values


def run_scan(
    label: str,
    graph: GraphData,
    omega: np.ndarray,
    k_values: Iterable[float],
    config: KuramotoConfig,
    rng: np.random.Generator,
) -> ScanRun:
    initial_phase = rng.uniform(-math.pi, math.pi, size=graph.node_count)
    triangles = find_triangles(graph.adjacency, limit=5000) if config.higher_order else []
    simulations = [
        simulate_kuramoto(label, float(k), graph.adjacency, omega, initial_phase, triangles, config)
        for k in k_values
    ]
    summary = summarize_scan(simulations)
    return ScanRun(label=label, graph=graph, omega=omega, simulations=simulations, summary=summary)


def simulate_kuramoto(
    label: str,
    k_value: float,
    adjacency: np.ndarray,
    omega: np.ndarray,
    initial_phase: np.ndarray,
    triangles: list[tuple[int, int, int]],
    config: KuramotoConfig,
) -> SimulationResult:
    """Integrate graph Kuramoto dynamics with optional triangle coupling."""

    node_count = adjacency.shape[0]
    phases = np.zeros((config.steps + 1, node_count), dtype=float)
    phases[0] = initial_phase.copy()
    order = np.zeros(config.steps + 1, dtype=float)
    order[0] = kuramoto_order(phases[0])

    degree = np.maximum(adjacency.sum(axis=1), EPS)
    for step in range(1, config.steps + 1):
        theta = phases[step - 1]
        pairwise = np.sin(theta[None, :] - theta[:, None])
        coupling = k_value * (adjacency * pairwise).sum(axis=1) / degree
        if triangles:
            coupling += config.higher_order_strength * k_value * triangle_coupling(theta, triangles, node_count)
        phases[step] = wrap_phase(theta + config.dt * (omega + coupling))
        order[step] = kuramoto_order(phases[step])

    recoverability, delta, identity = compute_haos_proxies(phases, order)
    k_star = first_sustained_crossing(recoverability, RECOVERABILITY_THRESHOLD, mode="ge", sustain=SUSTAIN_STEPS)
    r_star = first_sustained_crossing(order, STRONG_SYNC_THRESHOLD, mode="ge", sustain=SUSTAIN_STEPS)
    k_star_time = float(k_star * config.dt) if k_star is not None else None
    r_star_time = float(r_star * config.dt) if r_star is not None else None
    early = k_star_time is not None and (r_star_time is None or k_star_time < r_star_time)
    return SimulationResult(
        label=label,
        k_value=k_value,
        times=np.arange(config.steps + 1, dtype=float) * config.dt,
        phases=phases,
        order_parameter=order,
        recoverability=recoverability,
        delta_persistence=delta,
        trajectory_identity_retention=identity,
        k_star_time=k_star_time,
        r_star_time=r_star_time,
        early_detection=early,
    )


def compute_haos_proxies(phases: np.ndarray, order: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute bounded HAOS-style proxies from phase trajectories."""

    node_count = phases.shape[1]
    initial_offsets = centered_phase_offsets(phases[0])
    recoverability = np.zeros(phases.shape[0], dtype=float)
    identity = np.zeros(phases.shape[0], dtype=float)

    for idx, theta in enumerate(phases):
        offsets = centered_phase_offsets(theta)
        retention = 1.0 - circular_rmse(offsets, initial_offsets) / math.pi
        identity[idx] = float(np.clip(retention, 0.0, 1.0))
        local_coherence = float(np.mean(np.cos(offsets))) if node_count else 0.0
        recoverability[idx] = float(np.clip(0.55 * order[idx] + 0.30 * identity[idx] + 0.15 * max(local_coherence, 0.0), 0.0, 1.0))

    delta = np.zeros_like(recoverability)
    delta[1:] = recoverability[1:] - recoverability[:-1]
    return recoverability, delta, identity


def centered_phase_offsets(theta: np.ndarray) -> np.ndarray:
    mean_angle = math.atan2(float(np.mean(np.sin(theta))), float(np.mean(np.cos(theta))))
    return wrap_phase(theta - mean_angle)


def circular_rmse(a: np.ndarray, b: np.ndarray) -> float:
    diff = wrap_phase(a - b)
    return float(np.sqrt(np.mean(diff * diff)))


def kuramoto_order(theta: np.ndarray) -> float:
    return float(abs(np.mean(np.exp(1j * theta))))


def wrap_phase(theta: np.ndarray) -> np.ndarray:
    return (theta + math.pi) % (2.0 * math.pi) - math.pi


def first_sustained_crossing(values: np.ndarray, threshold: float, mode: str, sustain: int) -> int | None:
    for idx in range(0, len(values) - sustain):
        window = values[idx : idx + sustain + 1]
        if mode == "ge" and bool(np.all(window >= threshold)):
            return idx
        if mode == "le" and bool(np.all(window <= threshold)):
            return idx
    return None


def summarize_scan(simulations: list[SimulationResult]) -> dict[str, Any]:
    best = max(simulations, key=lambda item: float(item.order_parameter[-1]))
    earliest_proxy = min((sim for sim in simulations if sim.k_star_time is not None), key=lambda item: item.k_star_time, default=None)
    earliest_r = min((sim for sim in simulations if sim.r_star_time is not None), key=lambda item: item.r_star_time, default=None)
    early_count = sum(1 for sim in simulations if sim.early_detection)
    return {
        "best_k": best.k_value,
        "best_final_R": float(best.order_parameter[-1]),
        "best_final_recoverability": float(best.recoverability[-1]),
        "k_star_time": earliest_proxy.k_star_time if earliest_proxy else None,
        "r_star_time": earliest_r.r_star_time if earliest_r else None,
        "k_star_K": earliest_proxy.k_value if earliest_proxy else None,
        "r_star_K": earliest_r.k_value if earliest_r else None,
        "early_detection": bool(early_count > 0 and earliest_proxy is not None and (earliest_r is None or earliest_proxy.k_star_time < earliest_r.r_star_time)),
        "early_detection_count": early_count,
        "mean_final_R": float(np.mean([sim.order_parameter[-1] for sim in simulations])),
        "mean_final_recoverability": float(np.mean([sim.recoverability[-1] for sim in simulations])),
        "min_delta_persistence": float(min(np.min(sim.delta_persistence) for sim in simulations)),
        "max_identity_retention": float(max(np.max(sim.trajectory_identity_retention) for sim in simulations)),
    }


def edge_shuffle_graph(graph: GraphData, rng: np.random.Generator) -> GraphData:
    upper = graph.adjacency[np.triu_indices(graph.node_count, k=1)]
    shuffled = rng.permutation(upper)
    adjacency = np.zeros_like(graph.adjacency)
    adjacency[np.triu_indices(graph.node_count, k=1)] = shuffled
    adjacency = adjacency + adjacency.T
    return GraphData(normalize_adjacency(adjacency), f"{graph.source_label}:edge_shuffle", graph.source_path, "control")


def randomize_topology_graph(graph: GraphData, rng: np.random.Generator) -> GraphData:
    n = graph.node_count
    edge_count = graph.edge_count
    positive = graph.adjacency[np.triu_indices(n, k=1)]
    weights = positive[positive > 0.0]
    adjacency = np.zeros((n, n), dtype=float)
    possible = [(i, j) for i in range(n) for j in range(i + 1, n)]
    chosen = rng.choice(len(possible), size=min(edge_count, len(possible)), replace=False)
    sampled_weights = rng.choice(weights, size=len(chosen), replace=True) if weights.size else np.ones(len(chosen))
    for pick, weight in zip(chosen, sampled_weights):
        i, j = possible[int(pick)]
        adjacency[i, j] = float(weight)
        adjacency[j, i] = float(weight)
    return GraphData(normalize_adjacency(adjacency), f"{graph.source_label}:topology_randomized", graph.source_path, "control")


def find_triangles(adjacency: np.ndarray, limit: int) -> list[tuple[int, int, int]]:
    binary = adjacency > 0.0
    triangles: list[tuple[int, int, int]] = []
    n = adjacency.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            if not binary[i, j]:
                continue
            common = np.flatnonzero(binary[i] & binary[j])
            for k in common:
                if k > j:
                    triangles.append((i, j, int(k)))
                    if len(triangles) >= limit:
                        return triangles
    return triangles


def triangle_coupling(theta: np.ndarray, triangles: list[tuple[int, int, int]], node_count: int) -> np.ndarray:
    values = np.zeros(node_count, dtype=float)
    counts = np.zeros(node_count, dtype=float)
    for i, j, k in triangles:
        values[i] += math.sin(theta[j] + theta[k] - 2.0 * theta[i])
        values[j] += math.sin(theta[i] + theta[k] - 2.0 * theta[j])
        values[k] += math.sin(theta[i] + theta[j] - 2.0 * theta[k])
        counts[[i, j, k]] += 1.0
    return values / np.maximum(counts, 1.0)


def classify_status(graph: GraphData, observed: ScanRun, controls: list[ScanRun]) -> dict[str, Any]:
    control_early_count = sum(1 for run in controls if bool(run.summary["early_detection"]))
    observed_early = bool(observed.summary["early_detection"])
    observed_sync = float(observed.summary["best_final_R"]) >= STRONG_SYNC_THRESHOLD
    control_contrast = control_early_count == 0

    if graph.source_kind == "OPEN_NO_DATA_SYNTHETIC":
        bridge_status = "OPEN_NO_DATA_SYNTHETIC"
        failure_mode = "NO_HAOS_GRAPH_ARTIFACT_FOUND"
    elif observed_early and control_contrast:
        bridge_status = "PASS"
        failure_mode = "OBSERVED_EARLY_DETECTION_WITH_CONTROL_CONTRAST"
    elif observed_sync and control_contrast:
        bridge_status = "MARGINAL"
        failure_mode = "SYNC_WITHOUT_HAOS_EARLY_DETECTION"
    elif observed_early and not control_contrast:
        bridge_status = "FAIL"
        failure_mode = "CONTROL_MATCHED_EARLY_DETECTION"
    else:
        bridge_status = "FAIL"
        failure_mode = "PROXIES_UNDERPERFORM_STANDARD_R"

    return {
        "bridge_status": bridge_status,
        "failure_mode": failure_mode,
        "graph_source": graph.source_label,
        "graph_source_path": graph.source_path,
        "graph_source_kind": graph.source_kind,
        "nodes": graph.node_count,
        "edges": graph.edge_count,
        "mean_degree": graph.mean_degree,
        "observed_early_detection": observed_early,
        "control_early_detection_count": control_early_count,
        "control_contrast": control_contrast,
        "strong_sync_threshold": STRONG_SYNC_THRESHOLD,
        "recoverability_threshold": RECOVERABILITY_THRESHOLD,
    }


def write_outputs(config: KuramotoConfig, status: dict[str, Any], graph: GraphData, observed: ScanRun, controls: list[ScanRun]) -> None:
    output_dir = config.output_dir
    write_json(output_dir / "bridge_status.json", {**status, "observed_summary": observed.summary, "control_summaries": {run.label: run.summary for run in controls}})
    write_probe_comparison(output_dir / "probe_comparison.md", status, observed, controls)
    write_probe_csv(output_dir / "probe_comparison.csv", observed, controls)
    write_failure_analysis(output_dir / "failure_analysis.md", status, observed, controls)
    write_plots(output_dir, observed, controls)


def write_probe_comparison(path: Path, status: dict[str, Any], observed: ScanRun, controls: list[ScanRun]) -> None:
    runs = [observed, *controls]
    lines = [
        "# Kuramoto Bridge Probe Comparison",
        "",
        f"- bridge_status: {status['bridge_status']}",
        f"- failure_mode: {status['failure_mode']}",
        f"- graph_source: {status['graph_source']}",
        f"- graph_source_kind: {status['graph_source_kind']}",
        "",
        "| run | best_K | best_final_R | best_final_recoverability | k_star_time | r_star_time | early_detection |",
        "| --- | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for run in runs:
        summary = run.summary
        lines.append(
            "| {label} | {best_k:.6g} | {best_r:.6f} | {best_rec:.6f} | {kst} | {rst} | {early} |".format(
                label=run.label,
                best_k=float(summary["best_k"]),
                best_r=float(summary["best_final_R"]),
                best_rec=float(summary["best_final_recoverability"]),
                kst=format_optional(summary["k_star_time"]),
                rst=format_optional(summary["r_star_time"]),
                early=summary["early_detection"],
            )
        )
    lines.extend(
        [
            "",
            "## Metric Notes",
            "",
            "- `R` is the standard Kuramoto order parameter.",
            "- `recoverability_score` combines order, trajectory identity retention, and local phase coherence.",
            "- `delta_persistence` is the time derivative of recoverability.",
            "- `k_star_time` is the first sustained recoverability crossing.",
            "- `early_detection` means `k_star_time` occurs before the standard `R` crossing.",
            "",
            "## Boundary",
            "",
            "This bridge compares diagnostics under perturbation and controls. It does not promote a HAOS-IIP physics claim by itself.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_probe_csv(path: Path, observed: ScanRun, controls: list[ScanRun]) -> None:
    fields = ["run", "best_K", "best_final_R", "best_final_recoverability", "k_star_time", "r_star_time", "early_detection"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for run in [observed, *controls]:
            summary = run.summary
            writer.writerow(
                {
                    "run": run.label,
                    "best_K": summary["best_k"],
                    "best_final_R": summary["best_final_R"],
                    "best_final_recoverability": summary["best_final_recoverability"],
                    "k_star_time": summary["k_star_time"],
                    "r_star_time": summary["r_star_time"],
                    "early_detection": summary["early_detection"],
                }
            )


def write_failure_analysis(path: Path, status: dict[str, Any], observed: ScanRun, controls: list[ScanRun]) -> None:
    lines = [
        "# Kuramoto Bridge Failure Analysis",
        "",
        f"- bridge_status: {status['bridge_status']}",
        f"- failure_mode: {status['failure_mode']}",
        f"- observed_early_detection: {observed.summary['early_detection']}",
        f"- control_early_detection_count: {status['control_early_detection_count']}",
        "",
    ]
    if status["bridge_status"] == "PASS":
        lines.append("The observed run shows HAOS-proxy early detection and the controls do not reproduce it.")
    elif status["bridge_status"] == "MARGINAL":
        lines.append("The oscillator field synchronizes, but the HAOS proxies do not clearly beat the classic order parameter.")
    elif status["bridge_status"] == "OPEN_NO_DATA_SYNTHETIC":
        lines.append("No HAOS graph artifact was available. The run is a deterministic plumbing check, not a bridge result.")
    elif status["failure_mode"] == "CONTROL_MATCHED_EARLY_DETECTION":
        lines.append("At least one negative control reproduced early detection, so the proxy is not specific enough yet.")
    else:
        lines.append("The current proxies underperform or tie the standard Kuramoto order parameter on this scan.")
    lines.extend(
        [
            "",
            "Next bounded refinements:",
            "",
            "- replace the generic recoverability blend with the exact HAOS Line E recoverability function when available;",
            "- test graph-specific frozen cochain-Laplacian operators instead of reconstructed adjacency where serialized operators exist;",
            "- widen the K scan only after control contrast survives the default window.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_plots(output_dir: Path, observed: ScanRun, controls: list[ScanRun]) -> None:
    import matplotlib.pyplot as plt

    best = max(observed.simulations, key=lambda item: item.order_parameter[-1])
    plt.figure(figsize=(8, 4))
    for sim in observed.simulations:
        plt.plot(sim.times, sim.order_parameter, alpha=0.25, linewidth=1)
    plt.plot(best.times, best.order_parameter, color="black", linewidth=2, label=f"best K={best.k_value:.3g}")
    plt.axhline(STRONG_SYNC_THRESHOLD, color="firebrick", linestyle="--", linewidth=1, label="R threshold")
    plt.xlabel("time")
    plt.ylabel("R")
    plt.title("Kuramoto synchronization curves")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "synchronization_curves.png", dpi=160)
    plt.close()

    plt.figure(figsize=(8, 4))
    plt.plot(best.times, best.order_parameter, label="R")
    plt.plot(best.times, best.recoverability, label="recoverability_score")
    plt.plot(best.times, best.trajectory_identity_retention, label="trajectory_identity_retention")
    plt.plot(best.times, best.delta_persistence, label="delta_persistence")
    plt.axhline(RECOVERABILITY_THRESHOLD, color="firebrick", linestyle="--", linewidth=1)
    plt.xlabel("time")
    plt.ylabel("proxy value")
    plt.title("HAOS proxy evolution")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "proxy_evolution.png", dpi=160)
    plt.close()

    labels = [observed.label, *[run.label for run in controls]]
    best_r = [observed.summary["best_final_R"], *[run.summary["best_final_R"] for run in controls]]
    best_recoverability = [
        observed.summary["best_final_recoverability"],
        *[run.summary["best_final_recoverability"] for run in controls],
    ]
    x = np.arange(len(labels))
    width = 0.38
    plt.figure(figsize=(9, 4))
    plt.bar(x - width / 2, best_r, width, label="best final R")
    plt.bar(x + width / 2, best_recoverability, width, label="best final recoverability")
    plt.xticks(x, labels, rotation=20, ha="right")
    plt.ylim(0, 1.05)
    plt.ylabel("score")
    plt.title("Observed bridge vs negative controls")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "control_comparison.png", dpi=160)
    plt.close()


def pairwise_distances(points: np.ndarray) -> np.ndarray:
    diff = points[:, None, :] - points[None, :, :]
    return np.sqrt(np.sum(diff * diff, axis=2))


def write_json(path: Path, data: dict[str, Any]) -> None:
    serializable = to_jsonable(data)
    path.write_text(json.dumps(serializable, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def to_jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "__dataclass_fields__"):
        return to_jsonable(asdict(value))
    if isinstance(value, dict):
        return {str(key): to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(item) for item in value]
    return value


def relative_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT.resolve()))
    except ValueError:
        return str(path.resolve())


def format_optional(value: Any) -> str:
    if value is None:
        return ""
    return f"{float(value):.6g}"
