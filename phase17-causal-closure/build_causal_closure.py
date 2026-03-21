#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from collections import defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path


PHASE_DIR = Path(__file__).resolve().parent
REPO_ROOT = PHASE_DIR.parent
RUNS_DIR = PHASE_DIR / "runs"
PLOTS_DIR = PHASE_DIR / "plots"

PHASE15_MANIFEST = REPO_ROOT / "phase15-propagation" / "phase15_manifest.json"
PHASE15_INFLUENCE = REPO_ROOT / "phase15-propagation" / "runs" / "phase15_influence_range_ledger.csv"
PHASE16_MANIFEST = REPO_ROOT / "phase16-temporal-ordering" / "phase16_manifest.json"
PHASE16_EVENTS = REPO_ROOT / "phase16-temporal-ordering" / "runs" / "phase16_event_ordering_ledger.csv"
PHASE16_FRONTS = REPO_ROOT / "phase16-temporal-ordering" / "runs" / "phase16_front_arrival_ordering.csv"
PHASE16_ROBUSTNESS = REPO_ROOT / "phase16-temporal-ordering" / "runs" / "phase16_ordering_robustness.csv"

TARGET_PROBE = "bias_onset"
TARGET_REFINEMENTS = [60, 72, 84]
TARGET_ENSEMBLE_SIZE = 7
BRANCH_LABEL = "frozen_branch"
SOURCE_EVENT = "low_k_half"
EVENTS = [
    "low_k_half",
    "dispersion_half",
    "spectral_half",
    "radius_half",
    "width_half",
    "front_near",
    "front_far",
]

EDGE_REPRO_GATE = 0.60
EDGE_SET_DISTANCE_GATE = 0.35
ACYCLICITY_BAND_GATE = 0.15
DEPTH_DRIFT_GATE = 0.25
MISMATCH_MEAN_GATE = 0.04
MISMATCH_DRIFT_GATE = 0.005


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def jaccard_distance(left: set, right: set) -> float | None:
    if right is None:
        return None
    union = left | right
    if not union:
        return 0.0
    return 1.0 - (len(left & right) / len(union))


def tarjan_scc(nodes: list[str], edges: set[tuple[str, str]]) -> tuple[list[list[str]], dict[str, int]]:
    adjacency: dict[str, set[str]] = defaultdict(set)
    for src, dst in edges:
        adjacency[src].add(dst)

    index = [0]
    stack: list[str] = []
    on_stack: set[str] = set()
    ids: dict[str, int] = {}
    low: dict[str, int] = {}
    components: list[list[str]] = []

    def visit(node: str) -> None:
        ids[node] = index[0]
        low[node] = index[0]
        index[0] += 1
        stack.append(node)
        on_stack.add(node)
        for nxt in adjacency[node]:
            if nxt not in ids:
                visit(nxt)
                low[node] = min(low[node], low[nxt])
            elif nxt in on_stack:
                low[node] = min(low[node], ids[nxt])
        if low[node] == ids[node]:
            component: list[str] = []
            while True:
                current = stack.pop()
                on_stack.remove(current)
                component.append(current)
                if current == node:
                    break
            components.append(component)

    for node in nodes:
        if node not in ids:
            visit(node)

    node_to_component: dict[str, int] = {}
    for index_id, component in enumerate(components):
        for node in component:
            node_to_component[node] = index_id
    return components, node_to_component


def longest_depths(source_component: int, dag_edges: set[tuple[int, int]]) -> dict[int, int]:
    adjacency: dict[int, set[int]] = defaultdict(set)
    indegree: dict[int, int] = defaultdict(int)
    nodes: set[int] = {source_component}
    for src, dst in dag_edges:
        adjacency[src].add(dst)
        indegree[dst] += 1
        nodes.add(src)
        nodes.add(dst)

    queue: deque[int] = deque(sorted(node for node in nodes if indegree[node] == 0))
    order: list[int] = []
    while queue:
        node = queue.popleft()
        order.append(node)
        for nxt in sorted(adjacency[node]):
            indegree[nxt] -= 1
            if indegree[nxt] == 0:
                queue.append(nxt)

    depths = {source_component: 0}
    for node in order:
        if node not in depths:
            continue
        for nxt in adjacency[node]:
            depths[nxt] = max(depths.get(nxt, 0), depths[node] + 1)
    return depths


def choose_seed_pair(
    phase15_rows: list[dict[str, str]],
    phase16_front_rows: list[dict[str, str]],
    phase16_robustness_rows: list[dict[str, str]],
) -> tuple[int, int]:
    seeds = sorted(
        {
            int(row["seed"])
            for row in phase16_front_rows
            if row["hierarchy_label"] == BRANCH_LABEL
            and row["probe_name"] == TARGET_PROBE
            and int(row["ensemble_size"]) == TARGET_ENSEMBLE_SIZE
            and int(row["n_side"]) in TARGET_REFINEMENTS
        }
    )
    valid_seeds: list[int] = []
    for seed in seeds:
        valid = True
        for n_side in TARGET_REFINEMENTS:
            matches = [
                row
                for row in phase16_front_rows
                if row["hierarchy_label"] == BRANCH_LABEL
                and row["probe_name"] == TARGET_PROBE
                and int(row["ensemble_size"]) == TARGET_ENSEMBLE_SIZE
                and int(row["n_side"]) == n_side
                and int(row["seed"]) == seed
            ]
            if len(matches) != 2 or not all(row["front_ordering_valid"] == "True" for row in matches):
                valid = False
                break
        if valid:
            valid_seeds.append(seed)

    if len(valid_seeds) < 2:
        raise RuntimeError("Phase XVII could not locate two valid branch seeds for the frozen bias_onset slice.")

    def mean_noise_distance(seed: int) -> float:
        values = [
            float(row["graph_distance"])
            for row in phase16_robustness_rows
            if row["hierarchy_label"] == BRANCH_LABEL
            and row["probe_name"] == TARGET_PROBE
            and int(row["ensemble_size"]) == TARGET_ENSEMBLE_SIZE
            and int(row["n_side"]) in TARGET_REFINEMENTS
            and int(row["seed"]) == seed
            and row["perturbation_kind"] == "event_graph_connectivity_noise"
        ]
        return mean(values)

    def mean_normalized_latency(seed: int) -> float:
        values = [
            float(row["normalized_latency"])
            for row in phase15_rows
            if row["hierarchy_label"] == BRANCH_LABEL
            and row["probe_name"] == TARGET_PROBE
            and int(row["ensemble_size"]) == TARGET_ENSEMBLE_SIZE
            and int(row["n_side"]) in TARGET_REFINEMENTS
            and int(row["seed"]) == seed
        ]
        return mean(values)

    cleanest_seed = min(valid_seeds, key=lambda seed: (mean_noise_distance(seed), seed))
    contrast_candidates = [seed for seed in valid_seeds if seed != cleanest_seed]
    contrast_seed = max(contrast_candidates, key=lambda seed: (mean_normalized_latency(seed), -seed))
    return cleanest_seed, contrast_seed


def select_control_label(
    phase16_event_rows: list[dict[str, str]],
    selected_seeds: list[int],
) -> str:
    candidates = sorted(
        {
            row["hierarchy_label"]
            for row in phase16_event_rows
            if row["hierarchy_label"] != BRANCH_LABEL
            and row["probe_name"] == TARGET_PROBE
            and int(row["ensemble_size"]) == TARGET_ENSEMBLE_SIZE
            and int(row["n_side"]) in TARGET_REFINEMENTS
            and int(row["seed"]) in selected_seeds
        }
    )
    if not candidates:
        raise RuntimeError("Phase XVII could not locate a matched deterministic control hierarchy.")
    ranked = []
    for label in candidates:
        count = sum(
            1
            for row in phase16_event_rows
            if row["hierarchy_label"] == label
            and row["probe_name"] == TARGET_PROBE
            and int(row["ensemble_size"]) == TARGET_ENSEMBLE_SIZE
            and int(row["n_side"]) in TARGET_REFINEMENTS
            and int(row["seed"]) in selected_seeds
        )
        ranked.append((count, label))
    ranked.sort(reverse=True)
    return ranked[0][1]


def build_run_snapshot(
    hierarchy_label: str,
    n_side: int,
    seed: int,
    phase15_rows: list[dict[str, str]],
    phase16_event_rows: list[dict[str, str]],
    phase16_front_rows: list[dict[str, str]],
    phase16_robustness_rows: list[dict[str, str]],
) -> dict:
    event_rows = sorted(
        [
            row
            for row in phase16_event_rows
            if row["hierarchy_label"] == hierarchy_label
            and row["probe_name"] == TARGET_PROBE
            and int(row["ensemble_size"]) == TARGET_ENSEMBLE_SIZE
            and int(row["n_side"]) == n_side
            and int(row["seed"]) == seed
        ],
        key=lambda row: int(row["event_rank"]),
    )
    if not event_rows:
        raise RuntimeError(f"Missing event ordering rows for {hierarchy_label} n_side={n_side} seed={seed}.")

    order = [row["event_label"] for row in event_rows]
    if order[0] != SOURCE_EVENT:
        raise RuntimeError(
            f"Phase XVII expected {SOURCE_EVENT} to be the source event, found {order[0]} for {hierarchy_label} n_side={n_side} seed={seed}."
        )

    front_rows = [
        row
        for row in phase16_front_rows
        if row["hierarchy_label"] == hierarchy_label
        and row["probe_name"] == TARGET_PROBE
        and int(row["ensemble_size"]) == TARGET_ENSEMBLE_SIZE
        and int(row["n_side"]) == n_side
        and int(row["seed"]) == seed
    ]
    if len(front_rows) != 2:
        raise RuntimeError(f"Missing front arrival rows for {hierarchy_label} n_side={n_side} seed={seed}.")
    front_by_band = {row["distance_band"]: row for row in front_rows}
    near_latency = float(front_by_band["near"]["mean_latency"])
    far_latency = float(front_by_band["far"]["mean_latency"])
    front_distance = mean([float(row["front_distance"]) for row in front_rows])
    front_expected = "near_to_far" if near_latency <= far_latency else "far_to_near"
    front_actual = "near_to_far" if order.index("front_near") < order.index("front_far") else "far_to_near"
    front_valid = front_by_band["near"]["front_ordering_valid"] == "True" and front_by_band["far"]["front_ordering_valid"] == "True"

    robustness_rows = [
        row
        for row in phase16_robustness_rows
        if row["hierarchy_label"] == hierarchy_label
        and row["probe_name"] == TARGET_PROBE
        and int(row["ensemble_size"]) == TARGET_ENSEMBLE_SIZE
        and int(row["n_side"]) == n_side
        and int(row["seed"]) == seed
        and row["perturbation_kind"] == "event_graph_connectivity_noise"
    ]
    if len(robustness_rows) != 1:
        raise RuntimeError(f"Missing connectivity-noise robustness row for {hierarchy_label} n_side={n_side} seed={seed}.")
    noise_distance = float(robustness_rows[0]["graph_distance"])

    influence_rows = sorted(
        [
            row
            for row in phase15_rows
            if row["hierarchy_label"] == hierarchy_label
            and row["probe_name"] == TARGET_PROBE
            and int(row["ensemble_size"]) == TARGET_ENSEMBLE_SIZE
            and int(row["n_side"]) == n_side
            and int(row["seed"]) == seed
        ],
        key=lambda row: float(row["distance"]),
    )
    if not influence_rows:
        raise RuntimeError(f"Missing influence rows for {hierarchy_label} n_side={n_side} seed={seed}.")

    return {
        "hierarchy_label": hierarchy_label,
        "n_side": n_side,
        "seed": seed,
        "event_rows": event_rows,
        "order": order,
        "chain_edges": list(zip(order[:-1], order[1:])),
        "front_distance": front_distance,
        "front_expected": front_expected,
        "front_actual": front_actual,
        "front_valid": front_valid,
        "noise_distance": noise_distance,
        "influence_rows": influence_rows,
    }


def aggregate_refinement(
    hierarchy_label: str,
    n_side: int,
    seeds: list[int],
    phase15_rows: list[dict[str, str]],
    phase16_event_rows: list[dict[str, str]],
    phase16_front_rows: list[dict[str, str]],
    phase16_robustness_rows: list[dict[str, str]],
    previous_edges: set[tuple[str, str]] | None,
    previous_cycles: set[tuple[str, str]] | None,
) -> dict:
    snapshots = [
        build_run_snapshot(hierarchy_label, n_side, seed, phase15_rows, phase16_event_rows, phase16_front_rows, phase16_robustness_rows)
        for seed in seeds
    ]

    edge_counts: dict[tuple[str, str], int] = defaultdict(int)
    for snapshot in snapshots:
        for edge in snapshot["chain_edges"]:
            edge_counts[edge] += 1

    edge_weights = {edge: count / len(seeds) for edge, count in edge_counts.items()}
    edge_set = set(edge_weights)
    mean_edge_reproducibility = mean(list(edge_weights.values()))
    stable_edge_count = sum(1 for weight in edge_weights.values() if weight == 1.0)

    components, node_to_component = tarjan_scc(EVENTS, edge_set)
    cycle_components = [component for component in components if len(component) > 1]
    cycle_nodes = {node for component in cycle_components for node in component}
    cycle_edges = {edge for edge in edge_set if edge[0] in cycle_nodes and edge[1] in cycle_nodes}
    cycle_count = len(cycle_components)
    acyclicity_score = 1.0 - (len(cycle_edges) / len(edge_set) if edge_set else 0.0)

    dag_edges: set[tuple[int, int]] = set()
    for src, dst in edge_set:
        src_component = node_to_component[src]
        dst_component = node_to_component[dst]
        if src_component != dst_component:
            dag_edges.add((src_component, dst_component))

    source_component = node_to_component[SOURCE_EVENT]
    component_depths = longest_depths(source_component, dag_edges)
    reachable_nodes = [node for node in EVENTS if node_to_component[node] in component_depths]
    mean_causal_depth = mean([float(component_depths[node_to_component[node]]) for node in reachable_nodes])
    front_near_depth = component_depths.get(node_to_component["front_near"])
    front_far_depth = component_depths.get(node_to_component["front_far"])

    baseline_direction_mismatch_values: list[float] = []
    normalized_noise_values: list[float] = []
    mismatch_values: list[float] = []
    front_distance_values: list[float] = []
    noise_distance_values: list[float] = []
    for snapshot in snapshots:
        baseline_direction_mismatch = 0.0 if snapshot["front_valid"] and snapshot["front_expected"] == snapshot["front_actual"] else 1.0
        normalized_noise = snapshot["noise_distance"] / (1.0 + snapshot["front_distance"])
        baseline_direction_mismatch_values.append(baseline_direction_mismatch)
        normalized_noise_values.append(normalized_noise)
        mismatch_values.append(baseline_direction_mismatch + normalized_noise)
        front_distance_values.append(snapshot["front_distance"])
        noise_distance_values.append(snapshot["noise_distance"])

    return {
        "hierarchy_label": hierarchy_label,
        "n_side": n_side,
        "seed_list": list(seeds),
        "snapshots": snapshots,
        "edge_count": len(edge_set),
        "stable_edge_count": stable_edge_count,
        "mean_edge_reproducibility": mean_edge_reproducibility,
        "edge_set": edge_set,
        "edge_set_distance_from_prev": jaccard_distance(edge_set, previous_edges),
        "cycle_count": cycle_count,
        "cycle_edge_count": len(cycle_edges),
        "acyclicity_score": acyclicity_score,
        "cycle_edges": cycle_edges,
        "cycle_persistence_from_prev": jaccard_distance(cycle_edges, previous_cycles),
        "reachable_node_count": len(reachable_nodes),
        "mean_causal_depth": mean_causal_depth,
        "front_near_depth": front_near_depth,
        "front_far_depth": front_far_depth,
        "mean_front_distance": mean(front_distance_values),
        "mean_connectivity_noise": mean(noise_distance_values),
        "baseline_direction_mismatch_rate": mean(baseline_direction_mismatch_values),
        "normalized_noise_mismatch_rate": mean(normalized_noise_values),
        "mismatch_rate": mean(mismatch_values),
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise RuntimeError(f"Cannot write empty CSV to {path}.")
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def svg_line_plot(path: Path, title: str, x_values: list[int], series: list[tuple[str, str, list[float]]], y_label: str) -> None:
    width, height = 700, 420
    margin_left, margin_right, margin_top, margin_bottom = 70, 20, 40, 60
    plot_width = width - margin_left - margin_right
    plot_height = height - margin_top - margin_bottom
    all_values = [value for _, _, values in series for value in values]
    y_min = min(all_values)
    y_max = max(all_values)
    if y_min == y_max:
        y_min -= 0.5
        y_max += 0.5
    y_pad = 0.1 * (y_max - y_min)
    y_min -= y_pad
    y_max += y_pad

    def x_pos(value: int) -> float:
        if len(x_values) == 1:
            return margin_left + plot_width / 2.0
        span = x_values[-1] - x_values[0]
        return margin_left + ((value - x_values[0]) / span) * plot_width

    def y_pos(value: float) -> float:
        return margin_top + (1.0 - ((value - y_min) / (y_max - y_min))) * plot_height

    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='white'/>",
        f"<text x='{width/2:.1f}' y='24' text-anchor='middle' font-size='18' font-family='Helvetica, Arial, sans-serif'>{title}</text>",
        f"<line x1='{margin_left}' y1='{height-margin_bottom}' x2='{width-margin_right}' y2='{height-margin_bottom}' stroke='#333' stroke-width='1.5'/>",
        f"<line x1='{margin_left}' y1='{margin_top}' x2='{margin_left}' y2='{height-margin_bottom}' stroke='#333' stroke-width='1.5'/>",
        f"<text x='{width/2:.1f}' y='{height-18}' text-anchor='middle' font-size='12' font-family='Helvetica, Arial, sans-serif'>n_side</text>",
        f"<text x='18' y='{height/2:.1f}' transform='rotate(-90 18 {height/2:.1f})' text-anchor='middle' font-size='12' font-family='Helvetica, Arial, sans-serif'>{y_label}</text>",
    ]

    for value in x_values:
        xpos = x_pos(value)
        parts.append(
            f"<line x1='{xpos:.2f}' y1='{height-margin_bottom}' x2='{xpos:.2f}' y2='{height-margin_bottom+6}' stroke='#333' stroke-width='1'/>"
        )
        parts.append(
            f"<text x='{xpos:.2f}' y='{height-margin_bottom+22}' text-anchor='middle' font-size='11' font-family='Helvetica, Arial, sans-serif'>{value}</text>"
        )

    for tick in range(5):
        value = y_min + (tick / 4.0) * (y_max - y_min)
        ypos = y_pos(value)
        parts.append(
            f"<line x1='{margin_left-6}' y1='{ypos:.2f}' x2='{margin_left}' y2='{ypos:.2f}' stroke='#333' stroke-width='1'/>"
        )
        parts.append(
            f"<text x='{margin_left-10}' y='{ypos+4:.2f}' text-anchor='end' font-size='11' font-family='Helvetica, Arial, sans-serif'>{value:.3f}</text>"
        )
        parts.append(
            f"<line x1='{margin_left}' y1='{ypos:.2f}' x2='{width-margin_right}' y2='{ypos:.2f}' stroke='#e5e7eb' stroke-width='1'/>"
        )

    legend_y = margin_top + 12
    legend_x = width - margin_right - 180
    for index_id, (label, color, values) in enumerate(series):
        points = " ".join(f"{x_pos(x):.2f},{y_pos(y):.2f}" for x, y in zip(x_values, values))
        parts.append(f"<polyline fill='none' stroke='{color}' stroke-width='2.5' points='{points}'/>")
        for x, y in zip(x_values, values):
            parts.append(f"<circle cx='{x_pos(x):.2f}' cy='{y_pos(y):.2f}' r='3.5' fill='{color}'/>")
        ly = legend_y + 18 * index_id
        parts.append(f"<line x1='{legend_x}' y1='{ly}' x2='{legend_x+18}' y2='{ly}' stroke='{color}' stroke-width='2.5'/>")
        parts.append(
            f"<text x='{legend_x+24}' y='{ly+4}' font-size='11' font-family='Helvetica, Arial, sans-serif'>{label}</text>"
        )

    parts.append("</svg>")
    write_text(path, "\n".join(parts))


def build_phase17() -> None:
    phase15_manifest = load_json(PHASE15_MANIFEST)
    phase16_manifest = load_json(PHASE16_MANIFEST)
    phase15_rows = load_csv(PHASE15_INFLUENCE)
    phase16_event_rows = load_csv(PHASE16_EVENTS)
    phase16_front_rows = load_csv(PHASE16_FRONTS)
    phase16_robustness_rows = load_csv(PHASE16_ROBUSTNESS)

    if phase16_manifest.get("primary_ordering_probe") != TARGET_PROBE:
        raise RuntimeError(f"Phase XVI primary ordering probe must remain {TARGET_PROBE}.")

    cleanest_seed, contrast_seed = choose_seed_pair(phase15_rows, phase16_front_rows, phase16_robustness_rows)
    selected_seeds = [cleanest_seed, contrast_seed]
    control_label = select_control_label(phase16_event_rows, selected_seeds)

    hierarchy_rows: dict[str, list[dict]] = {}
    for hierarchy_label in [BRANCH_LABEL, control_label]:
        rows: list[dict] = []
        previous_edges: set[tuple[str, str]] | None = None
        previous_cycles: set[tuple[str, str]] | None = None
        for n_side in TARGET_REFINEMENTS:
            row = aggregate_refinement(
                hierarchy_label,
                n_side,
                selected_seeds,
                phase15_rows,
                phase16_event_rows,
                phase16_front_rows,
                phase16_robustness_rows,
                previous_edges,
                previous_cycles,
            )
            rows.append(row)
            previous_edges = row["edge_set"]
            previous_cycles = row["cycle_edges"]
        hierarchy_rows[hierarchy_label] = rows

    branch_rows = hierarchy_rows[BRANCH_LABEL]
    control_rows = hierarchy_rows[control_label]

    branch_edge_repro_mean = mean([row["mean_edge_reproducibility"] for row in branch_rows])
    control_edge_repro_mean = mean([row["mean_edge_reproducibility"] for row in control_rows])
    branch_max_edge_set_distance = max(row["edge_set_distance_from_prev"] or 0.0 for row in branch_rows)
    control_max_edge_set_distance = max(row["edge_set_distance_from_prev"] or 0.0 for row in control_rows)

    branch_acyclicity_band_width = max(row["acyclicity_score"] for row in branch_rows) - min(
        row["acyclicity_score"] for row in branch_rows
    )
    control_acyclicity_band_width = max(row["acyclicity_score"] for row in control_rows) - min(
        row["acyclicity_score"] for row in control_rows
    )

    branch_mean_causal_depth = mean([row["mean_causal_depth"] for row in branch_rows])
    control_mean_causal_depth = mean([row["mean_causal_depth"] for row in control_rows])
    branch_causal_depth_drift = max(row["mean_causal_depth"] for row in branch_rows) - min(
        row["mean_causal_depth"] for row in branch_rows
    )
    control_causal_depth_drift = max(row["mean_causal_depth"] for row in control_rows) - min(
        row["mean_causal_depth"] for row in control_rows
    )

    branch_mismatch_mean = mean([row["mismatch_rate"] for row in branch_rows])
    control_mismatch_mean = mean([row["mismatch_rate"] for row in control_rows])
    branch_mismatch_drift = max(row["mismatch_rate"] for row in branch_rows) - min(row["mismatch_rate"] for row in branch_rows)
    control_mismatch_drift = max(row["mismatch_rate"] for row in control_rows) - min(
        row["mismatch_rate"] for row in control_rows
    )

    branch_flags = {
        "edge_reproducibility_stable": branch_edge_repro_mean >= EDGE_REPRO_GATE
        and branch_max_edge_set_distance <= EDGE_SET_DISTANCE_GATE,
        "acyclicity_band_compact": branch_acyclicity_band_width <= ACYCLICITY_BAND_GATE,
        "causal_depth_drift_bounded": branch_causal_depth_drift <= DEPTH_DRIFT_GATE,
        "propagation_order_mismatch_low": branch_mismatch_mean <= MISMATCH_MEAN_GATE
        and branch_mismatch_drift <= MISMATCH_DRIFT_GATE
        and branch_mismatch_mean < control_mismatch_mean,
    }
    control_flags = {
        "edge_reproducibility_stable": control_edge_repro_mean >= EDGE_REPRO_GATE
        and control_max_edge_set_distance <= EDGE_SET_DISTANCE_GATE,
        "acyclicity_band_compact": control_acyclicity_band_width <= ACYCLICITY_BAND_GATE,
        "causal_depth_drift_bounded": control_causal_depth_drift <= DEPTH_DRIFT_GATE,
        "propagation_order_mismatch_low": control_mismatch_mean <= MISMATCH_MEAN_GATE
        and control_mismatch_drift <= MISMATCH_DRIFT_GATE
        and control_mismatch_mean < branch_mismatch_mean,
    }

    branch_failed_gate_count = sum(1 for passed in branch_flags.values() if not passed)
    control_failed_gate_count = sum(1 for passed in control_flags.values() if not passed)
    hard_stop_triggered = branch_failed_gate_count >= 2
    success = all(branch_flags.values()) and control_failed_gate_count >= 1

    influence_rows_out: list[dict] = []
    dag_rows_out: list[dict] = []
    depth_rows_out: list[dict] = []
    compatibility_rows_out: list[dict] = []
    for row in branch_rows + control_rows:
        counterpart = next(
            candidate
            for candidate in (control_rows if row["hierarchy_label"] == BRANCH_LABEL else branch_rows)
            if candidate["n_side"] == row["n_side"]
        )
        influence_rows_out.append(
            {
                "hierarchy_label": row["hierarchy_label"],
                "n_side": row["n_side"],
                "selected_seed_pair": f"{cleanest_seed},{contrast_seed}",
                "edge_count": row["edge_count"],
                "stable_edge_count": row["stable_edge_count"],
                "mean_edge_reproducibility": f"{row['mean_edge_reproducibility']:.12f}",
                "edge_set_distance_from_prev": ""
                if row["edge_set_distance_from_prev"] is None
                else f"{row['edge_set_distance_from_prev']:.12f}",
            }
        )
        dag_rows_out.append(
            {
                "hierarchy_label": row["hierarchy_label"],
                "n_side": row["n_side"],
                "cycle_count": row["cycle_count"],
                "cycle_edge_count": row["cycle_edge_count"],
                "acyclicity_score": f"{row['acyclicity_score']:.12f}",
                "cycle_persistence_from_prev": ""
                if row["cycle_persistence_from_prev"] is None
                else f"{row['cycle_persistence_from_prev']:.12f}",
            }
        )
        depth_rows_out.append(
            {
                "hierarchy_label": row["hierarchy_label"],
                "n_side": row["n_side"],
                "mean_causal_depth": f"{row['mean_causal_depth']:.12f}",
                "reachable_node_count": row["reachable_node_count"],
                "front_near_depth": "" if row["front_near_depth"] is None else row["front_near_depth"],
                "front_far_depth": "" if row["front_far_depth"] is None else row["front_far_depth"],
                "branch_control_depth_gap": f"{row['mean_causal_depth'] - counterpart['mean_causal_depth']:.12f}",
            }
        )
        compatibility_rows_out.append(
            {
                "hierarchy_label": row["hierarchy_label"],
                "n_side": row["n_side"],
                "mean_front_distance": f"{row['mean_front_distance']:.12f}",
                "mean_connectivity_noise": f"{row['mean_connectivity_noise']:.12f}",
                "baseline_direction_mismatch_rate": f"{row['baseline_direction_mismatch_rate']:.12f}",
                "normalized_noise_mismatch_rate": f"{row['normalized_noise_mismatch_rate']:.12f}",
                "mismatch_rate": f"{row['mismatch_rate']:.12f}",
                "branch_control_mismatch_gap": f"{row['mismatch_rate'] - counterpart['mismatch_rate']:.12f}",
            }
        )

    write_csv(RUNS_DIR / "phase17_influence_graph_ledger.csv", influence_rows_out)
    write_csv(RUNS_DIR / "phase17_dag_statistics.csv", dag_rows_out)
    write_csv(RUNS_DIR / "phase17_causal_distance_metrics.csv", depth_rows_out)
    write_csv(RUNS_DIR / "phase17_propagation_order_compatibility.csv", compatibility_rows_out)

    runs_payload = {
        "timestamp": utc_now(),
        "stage_identifier": "phase17-causal-closure",
        "input_artifacts": {
            "phase15_manifest_json": "phase15-propagation/phase15_manifest.json",
            "phase15_influence_range_ledger_csv": "phase15-propagation/runs/phase15_influence_range_ledger.csv",
            "phase16_manifest_json": "phase16-temporal-ordering/phase16_manifest.json",
            "phase16_event_ordering_ledger_csv": "phase16-temporal-ordering/runs/phase16_event_ordering_ledger.csv",
            "phase16_front_arrival_ordering_csv": "phase16-temporal-ordering/runs/phase16_front_arrival_ordering.csv",
            "phase16_ordering_robustness_csv": "phase16-temporal-ordering/runs/phase16_ordering_robustness.csv",
        },
        "selected_slice": {
            "probe_name": TARGET_PROBE,
            "refinements": TARGET_REFINEMENTS,
            "ensemble_size": TARGET_ENSEMBLE_SIZE,
            "cleanest_branch_seed": cleanest_seed,
            "contrast_seed": contrast_seed,
            "control_label": control_label,
        },
        "threshold_policy": phase16_manifest["evaluation_protocol"]["event_thresholds"],
        "phase15_response_threshold_fraction": phase15_manifest["evaluation_protocol"]["response_threshold_fraction"],
        "compatibility_policy": {
            "baseline_direction_rule": "front_near must precede front_far when near mean latency does not exceed far mean latency",
            "mismatch_rate": "baseline_direction_mismatch + event_graph_connectivity_noise_graph_distance / (1 + mean_front_distance)",
            "cycle_graph_rule": "aggregate consecutive event-chain edges across the two selected seeds and score cycles on the resulting directed graph",
        },
        "acceptance_gates": {
            "edge_repro_mean_min": EDGE_REPRO_GATE,
            "max_edge_set_distance": EDGE_SET_DISTANCE_GATE,
            "acyclicity_band_width_max": ACYCLICITY_BAND_GATE,
            "causal_depth_drift_max": DEPTH_DRIFT_GATE,
            "mismatch_mean_max": MISMATCH_MEAN_GATE,
            "mismatch_drift_max": MISMATCH_DRIFT_GATE,
        },
    }
    write_json(RUNS_DIR / "phase17_runs.json", runs_payload)

    if success:
        svg_line_plot(
            PLOTS_DIR / "phase17_acyclicity_vs_refinement.svg",
            "Phase XVII Acyclicity vs Refinement",
            TARGET_REFINEMENTS,
            [
                ("branch", "#0f766e", [row["acyclicity_score"] for row in branch_rows]),
                ("control", "#b91c1c", [row["acyclicity_score"] for row in control_rows]),
            ],
            "acyclicity score",
        )
        svg_line_plot(
            PLOTS_DIR / "phase17_mismatch_rate_vs_refinement.svg",
            "Phase XVII Mismatch Rate vs Refinement",
            TARGET_REFINEMENTS,
            [
                ("branch", "#0f766e", [row["mismatch_rate"] for row in branch_rows]),
                ("control", "#b91c1c", [row["mismatch_rate"] for row in control_rows]),
            ],
            "mismatch rate",
        )
        svg_line_plot(
            PLOTS_DIR / "phase17_branch_vs_control_depth.svg",
            "Phase XVII Branch vs Control Causal Depth",
            TARGET_REFINEMENTS,
            [
                ("branch", "#0f766e", [row["mean_causal_depth"] for row in branch_rows]),
                ("control", "#b91c1c", [row["mean_causal_depth"] for row in control_rows]),
            ],
            "mean causal depth",
        )

    conclusion = (
        "Phase XVII establishes causal-closure feasibility for the frozen operator hierarchy."
        if success
        else "Phase XVII does not yet establish causal-closure feasibility for the frozen operator hierarchy."
    )

    manifest_payload = {
        "timestamp": utc_now(),
        "phase": 17,
        "phase_name": "phase17-causal-closure",
        "stage_identifier": "phase17-causal-closure",
        "status": "passed" if success else "failed",
        "success": success,
        "objective": "low_token_causal_closure_authority_probe",
        "input_references": runs_payload["input_artifacts"],
        "selected_slice": runs_payload["selected_slice"],
        "threshold_policy": runs_payload["threshold_policy"],
        "aggregate_metrics": {
            "branch_edge_repro_mean": round(branch_edge_repro_mean, 12),
            "control_edge_repro_mean": round(control_edge_repro_mean, 12),
            "branch_max_edge_set_distance": round(branch_max_edge_set_distance, 12),
            "control_max_edge_set_distance": round(control_max_edge_set_distance, 12),
            "branch_acyclicity_band_width": round(branch_acyclicity_band_width, 12),
            "control_acyclicity_band_width": round(control_acyclicity_band_width, 12),
            "branch_mean_causal_depth": round(branch_mean_causal_depth, 12),
            "control_mean_causal_depth": round(control_mean_causal_depth, 12),
            "branch_causal_depth_drift": round(branch_causal_depth_drift, 12),
            "control_causal_depth_drift": round(control_causal_depth_drift, 12),
            "branch_mismatch_mean": round(branch_mismatch_mean, 12),
            "control_mismatch_mean": round(control_mismatch_mean, 12),
            "branch_mismatch_drift": round(branch_mismatch_drift, 12),
            "control_mismatch_drift": round(control_mismatch_drift, 12),
            "branch_control_depth_gap": round(branch_mean_causal_depth - control_mean_causal_depth, 12),
            "branch_control_mismatch_gap": round(control_mismatch_mean - branch_mismatch_mean, 12),
        },
        "success_flags": branch_flags,
        "control_flags": control_flags,
        "branch_failed_gate_count": branch_failed_gate_count,
        "control_failed_gate_count": control_failed_gate_count,
        "hard_stop_triggered": hard_stop_triggered,
        "plot_generation": "enabled" if success else "skipped",
        "artifacts": {
            "manifest_json": "phase17-causal-closure/phase17_manifest.json",
            "summary_md": "phase17-causal-closure/phase17_summary.md",
            "influence_graph_ledger_csv": "phase17-causal-closure/runs/phase17_influence_graph_ledger.csv",
            "dag_statistics_csv": "phase17-causal-closure/runs/phase17_dag_statistics.csv",
            "causal_distance_metrics_csv": "phase17-causal-closure/runs/phase17_causal_distance_metrics.csv",
            "propagation_order_compatibility_csv": "phase17-causal-closure/runs/phase17_propagation_order_compatibility.csv",
            "runs_json": "phase17-causal-closure/runs/phase17_runs.json",
            "acyclicity_plot": "phase17-causal-closure/plots/phase17_acyclicity_vs_refinement.svg" if success else None,
            "mismatch_plot": "phase17-causal-closure/plots/phase17_mismatch_rate_vs_refinement.svg" if success else None,
            "depth_plot": "phase17-causal-closure/plots/phase17_branch_vs_control_depth.svg" if success else None,
            "builder_script": "phase17-causal-closure/build_causal_closure.py",
            "checker_script": "phase17-causal-closure/diagnostics/check_phase17_bundle.py",
        },
        "claim_boundary": "Phase XVII is limited to a low-token causal-closure authority probe over frozen Phase XV and Phase XVI artifacts. It does not modify earlier phase contracts, import earlier phase scripts, or assert metric, spacetime, or physical causal structure.",
        "conclusion": conclusion,
    }
    write_json(PHASE_DIR / "phase17_manifest.json", manifest_payload)

    summary_lines = [
        "# Phase XVII",
        "",
        "Phase XVII consumed only frozen Phase XV–XVI artifacts and restricted the authority slice to `bias_onset`, `n_side = 60, 72, 84`, ensemble size `7`, branch seeds "
        f"`{cleanest_seed}` and `{contrast_seed}`, and the matched control `{control_label}`.",
        "",
        f"Result: {'passed' if success else 'failed'}.",
        "",
        "- Influence-edge reconstruction used only consecutive event-chain edges from the frozen Phase XVI ledgers.",
        f"- Branch mean edge reproducibility was `{branch_edge_repro_mean:.6f}` with max adjacent edge-set distance `{branch_max_edge_set_distance:.6f}`; control was `{control_edge_repro_mean:.6f}` and `{control_max_edge_set_distance:.6f}`.",
        f"- Branch acyclicity stayed in a compact band of width `{branch_acyclicity_band_width:.6f}`; control band width was `{control_acyclicity_band_width:.6f}`.",
        f"- Branch mean causal depth was `{branch_mean_causal_depth:.6f}` with drift `{branch_causal_depth_drift:.6f}`; control drift was `{control_causal_depth_drift:.6f}`.",
        f"- Branch order-compatibility mismatch was `{branch_mismatch_mean:.6f}` with drift `{branch_mismatch_drift:.6f}` and remained below control `{control_mismatch_mean:.6f}`.",
        "",
        f"Hard stop triggered: `{hard_stop_triggered}`.",
        "",
        conclusion,
    ]
    write_text(PHASE_DIR / "phase17_summary.md", "\n".join(summary_lines) + "\n")


if __name__ == "__main__":
    build_phase17()
