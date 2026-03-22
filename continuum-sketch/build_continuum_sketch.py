#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
from collections import defaultdict, deque
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = Path(__file__).resolve().parent

TABLE_PATH = OUT_DIR / "continuum_sketch_table.csv"
PLOT_PATH = OUT_DIR / "continuum_sketch_plot_family.svg"
SUMMARY_PATH = OUT_DIR / "continuum_sketch_summary.md"

PHASE11_PERSISTENCE = ROOT / "phase11-protection" / "runs" / "phase11_persistence_scaling.csv"
PHASE15_SPEED = ROOT / "phase15-propagation" / "runs" / "phase15_effective_speed_ledger.csv"
PHASE16_EVENTS = ROOT / "phase16-temporal-ordering" / "runs" / "phase16_event_ordering_ledger.csv"
PHASE16_FRONTS = ROOT / "phase16-temporal-ordering" / "runs" / "phase16_front_arrival_ordering.csv"

LEVELS = [48, 60, 72, 84]
SEEDS = ["1302", "1303"]
ENSEMBLE_SIZE = "7"
PROBE_NAME = "bias_onset"
BRANCH = "frozen_branch"
CONTROL = "periodic_diagonal_augmented_control"
HIERARCHIES = [BRANCH, CONTROL]
TAU_MAX = 4.2
EVENT_LABELS = [
    "low_k_half",
    "dispersion_half",
    "spectral_half",
    "radius_half",
    "width_half",
    "front_near",
    "front_far",
]
SOURCE_EVENT = "low_k_half"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def mean(values) -> float:
    values = list(values)
    return sum(values) / len(values) if values else 0.0


def std(values: list[float]) -> float:
    if not values:
        return 0.0
    mu = mean(values)
    return math.sqrt(sum((value - mu) ** 2 for value in values) / len(values))


def coeff_var(values: list[float]) -> float:
    mu = mean(values)
    return 0.0 if abs(mu) < 1.0e-12 else std(values) / abs(mu)


def relative_span(values: list[float]) -> float:
    mu = mean(values)
    return 0.0 if abs(mu) < 1.0e-12 else (max(values) - min(values)) / abs(mu)


def linear_fit(xs: list[float], ys: list[float]) -> tuple[float, float, float]:
    x_mean = mean(xs)
    y_mean = mean(ys)
    numerator = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys))
    denominator = sum((x - x_mean) ** 2 for x in xs)
    slope = 0.0 if denominator == 0.0 else numerator / denominator
    intercept = y_mean - slope * x_mean
    predictions = [intercept + slope * x for x in xs]
    ss_res = sum((y - prediction) ** 2 for y, prediction in zip(ys, predictions))
    ss_tot = sum((y - y_mean) ** 2 for y in ys)
    r2 = 1.0 if ss_tot == 0.0 else 1.0 - ss_res / ss_tot
    return intercept, slope, r2


def log_fit(xs: list[float], ys: list[float]) -> tuple[float, float, float]:
    log_xs = [math.log(value) for value in xs]
    log_ys = [math.log(value) for value in ys]
    return linear_fit(log_xs, log_ys)


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


def least_squares_slope(xs: list[float], ys: list[float]) -> float:
    _, slope, _ = linear_fit(xs, ys)
    return slope


def shell_slope_by_level(
    event_rows: list[dict[str, str]],
    front_rows: list[dict[str, str]],
    hierarchy_label: str,
    n_side: int,
) -> tuple[float, float]:
    edge_counts: dict[tuple[str, str], int] = defaultdict(int)
    for seed in SEEDS:
        rows = sorted(
            [
                row
                for row in event_rows
                if row["hierarchy_label"] == hierarchy_label
                and row["n_side"] == str(n_side)
                and row["seed"] == seed
                and row["probe_name"] == PROBE_NAME
                and row["ensemble_size"] == ENSEMBLE_SIZE
            ],
            key=lambda row: int(row["event_rank"]),
        )
        order = [row["event_label"] for row in rows]
        for left, right in zip(order[:-1], order[1:]):
            edge_counts[(left, right)] += 1

    edge_set = set(edge_counts)
    _, node_to_component = tarjan_scc(EVENT_LABELS, edge_set)
    dag_edges: set[tuple[int, int]] = set()
    for src, dst in edge_set:
        src_component = node_to_component[src]
        dst_component = node_to_component[dst]
        if src_component != dst_component:
            dag_edges.add((src_component, dst_component))
    component_depths = longest_depths(node_to_component[SOURCE_EVENT], dag_edges)

    front_near_depth = float(component_depths[node_to_component["front_near"]])
    front_far_depth = float(component_depths[node_to_component["front_far"]])

    near_values = [
        float(row["mean_latency"])
        for row in front_rows
        if row["hierarchy_label"] == hierarchy_label
        and row["n_side"] == str(n_side)
        and row["seed"] in SEEDS
        and row["probe_name"] == PROBE_NAME
        and row["ensemble_size"] == ENSEMBLE_SIZE
        and row["distance_band"] == "near"
    ]
    far_values = [
        float(row["mean_latency"])
        for row in front_rows
        if row["hierarchy_label"] == hierarchy_label
        and row["n_side"] == str(n_side)
        and row["seed"] in SEEDS
        and row["probe_name"] == PROBE_NAME
        and row["ensemble_size"] == ENSEMBLE_SIZE
        and row["distance_band"] == "far"
    ]

    slope = least_squares_slope(
        [0.0, front_near_depth, front_far_depth],
        [0.0, mean(near_values), mean(far_values)],
    )
    return slope, mean(far_values) - mean(near_values)


def ordering_consistency_by_level(event_rows: list[dict[str, str]], hierarchy_label: str, n_side: int) -> float:
    event_times: dict[str, dict[str, float]] = defaultdict(dict)
    for row in event_rows:
        if (
            row["hierarchy_label"] == hierarchy_label
            and row["n_side"] == str(n_side)
            and row["seed"] in SEEDS
            and row["probe_name"] == PROBE_NAME
            and row["ensemble_size"] == ENSEMBLE_SIZE
        ):
            event_times[row["seed"]][row["event_label"]] = float(row["event_time"])

    left_seed, right_seed = SEEDS
    distance = mean(
        [
            abs(event_times[left_seed][label] - event_times[right_seed][label])
            for label in EVENT_LABELS
        ]
    )
    return 1.0 - distance / TAU_MAX


def generate_svg(
    path: Path,
    levels: list[int],
    series: dict[str, dict[str, list[float]]],
    titles: dict[str, str],
) -> None:
    width = 1180
    height = 820
    cols = 2
    rows = 3
    panel_width = 520
    panel_height = 210
    margin_x = 50
    margin_y = 50
    gap_x = 50
    gap_y = 45
    colors = {BRANCH: "#0f766e", CONTROL: "#b45309"}
    labels = {BRANCH: "branch", CONTROL: "control"}

    def panel_origin(index_id: int) -> tuple[float, float]:
        row = index_id // cols
        col = index_id % cols
        x0 = margin_x + col * (panel_width + gap_x)
        y0 = margin_y + row * (panel_height + gap_y)
        return x0, y0

    def x_position(x0: float, value: int) -> float:
        inner_left = x0 + 55
        inner_right = x0 + panel_width - 20
        span = levels[-1] - levels[0]
        return inner_left + ((value - levels[0]) / span) * (inner_right - inner_left)

    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='white'/>",
        "<text x='590' y='28' text-anchor='middle' font-size='22' font-family='Helvetica, Arial, sans-serif'>Continuum Sketch Observable Family</text>",
    ]

    observable_order = list(series.keys())
    for index_id, observable_name in enumerate(observable_order):
        x0, y0 = panel_origin(index_id)
        inner_left = x0 + 55
        inner_right = x0 + panel_width - 20
        inner_top = y0 + 35
        inner_bottom = y0 + panel_height - 35
        values = series[observable_name][BRANCH] + series[observable_name][CONTROL]
        y_min = min(values)
        y_max = max(values)
        if abs(y_max - y_min) < 1.0e-12:
            y_min -= 0.5
            y_max += 0.5
        pad = 0.10 * (y_max - y_min)
        y_min -= pad
        y_max += pad

        def y_position(value: float) -> float:
            return inner_bottom - ((value - y_min) / (y_max - y_min)) * (inner_bottom - inner_top)

        parts.append(f"<rect x='{x0}' y='{y0}' width='{panel_width}' height='{panel_height}' fill='none' stroke='#d1d5db' stroke-width='1'/>")
        parts.append(f"<text x='{x0 + panel_width / 2:.1f}' y='{y0 + 20}' text-anchor='middle' font-size='14' font-family='Helvetica, Arial, sans-serif'>{titles[observable_name]}</text>")
        parts.append(f"<line x1='{inner_left}' y1='{inner_bottom}' x2='{inner_right}' y2='{inner_bottom}' stroke='#374151' stroke-width='1.2'/>")
        parts.append(f"<line x1='{inner_left}' y1='{inner_top}' x2='{inner_left}' y2='{inner_bottom}' stroke='#374151' stroke-width='1.2'/>")

        for tick in range(4):
            value = y_min + (tick / 3.0) * (y_max - y_min)
            ypos = y_position(value)
            parts.append(f"<line x1='{inner_left}' y1='{ypos:.2f}' x2='{inner_right}' y2='{ypos:.2f}' stroke='#f3f4f6' stroke-width='1'/>")
            parts.append(f"<text x='{inner_left - 8}' y='{ypos + 4:.2f}' text-anchor='end' font-size='10' font-family='Helvetica, Arial, sans-serif'>{value:.3f}</text>")

        for level in levels:
            xpos = x_position(x0, level)
            parts.append(f"<line x1='{xpos:.2f}' y1='{inner_bottom}' x2='{xpos:.2f}' y2='{inner_bottom + 4}' stroke='#374151' stroke-width='1'/>")
            parts.append(f"<text x='{xpos:.2f}' y='{inner_bottom + 18}' text-anchor='middle' font-size='10' font-family='Helvetica, Arial, sans-serif'>{level}</text>")

        for hierarchy_label in HIERARCHIES:
            points = " ".join(
                f"{x_position(x0, level):.2f},{y_position(value):.2f}"
                for level, value in zip(levels, series[observable_name][hierarchy_label])
            )
            color = colors[hierarchy_label]
            parts.append(f"<polyline fill='none' stroke='{color}' stroke-width='2.2' points='{points}'/>")
            for level, value in zip(levels, series[observable_name][hierarchy_label]):
                parts.append(f"<circle cx='{x_position(x0, level):.2f}' cy='{y_position(value):.2f}' r='3' fill='{color}'/>")

        legend_y = inner_top + 12
        legend_x = inner_right - 110
        for offset, hierarchy_label in enumerate(HIERARCHIES):
            ly = legend_y + offset * 16
            color = colors[hierarchy_label]
            parts.append(f"<line x1='{legend_x}' y1='{ly}' x2='{legend_x + 18}' y2='{ly}' stroke='{color}' stroke-width='2.2'/>")
            parts.append(f"<text x='{legend_x + 24}' y='{ly + 4}' font-size='10' font-family='Helvetica, Arial, sans-serif'>{labels[hierarchy_label]}</text>")

    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    phase11_rows = read_csv(PHASE11_PERSISTENCE)
    phase15_speed_rows = read_csv(PHASE15_SPEED)
    phase16_event_rows = read_csv(PHASE16_EVENTS)
    phase16_front_rows = read_csv(PHASE16_FRONTS)

    table_rows: list[dict[str, object]] = []
    series: dict[str, dict[str, list[float]]] = {
        "dispersion_proxy": {BRANCH: [], CONTROL: []},
        "effective_speed_band": {BRANCH: [], CONTROL: []},
        "persistence_time": {BRANCH: [], CONTROL: []},
        "ordering_consistency_score": {BRANCH: [], CONTROL: []},
        "distance_surrogate_shell_slope": {BRANCH: [], CONTROL: []},
    }

    for hierarchy_label in HIERARCHIES:
        for n_side in LEVELS:
            h_value = 1.0 / float(n_side)

            persistence_row = next(
                row
                for row in phase11_rows
                if row["hierarchy_label"] == hierarchy_label and int(row["n_side"]) == n_side
            )
            spectral_gap = float(persistence_row["spectral_gap"])
            persistence_time = float(persistence_row["persistence_time_tau"])

            speed_values = [
                float(row["effective_speed"])
                for row in phase15_speed_rows
                if row["hierarchy_label"] == hierarchy_label
                and int(row["n_side"]) == n_side
                and row["probe_name"] == PROBE_NAME
                and row["ensemble_size"] == ENSEMBLE_SIZE
                and row["seed"] in SEEDS
            ]
            effective_speed = mean(speed_values)

            ordering_score = ordering_consistency_by_level(phase16_event_rows, hierarchy_label, n_side)
            shell_slope, shell_latency_gap = shell_slope_by_level(phase16_event_rows, phase16_front_rows, hierarchy_label, n_side)

            observable_values = {
                "dispersion_proxy": spectral_gap,
                "effective_speed_band": effective_speed,
                "persistence_time": persistence_time,
                "ordering_consistency_score": ordering_score,
                "distance_surrogate_shell_slope": shell_slope,
            }
            sources = {
                "dispersion_proxy": "phase11:spectral_gap",
                "effective_speed_band": "phase15:effective_speed",
                "persistence_time": "phase11:persistence_time_tau",
                "ordering_consistency_score": "phase16:event_time_distance_seed_pair",
                "distance_surrogate_shell_slope": "phase18:arrival_vs_depth_slope_formula",
            }

            for observable_name, value in observable_values.items():
                table_rows.append(
                    {
                        "hierarchy_label": hierarchy_label,
                        "n_side": n_side,
                        "h": f"{h_value:.12f}",
                        "seeds": ",".join(SEEDS),
                        "observable_name": observable_name,
                        "value": f"{value:.12f}",
                        "source_metric": sources[observable_name],
                        "notes": (
                            "primary_bias"
                            if observable_name in {"effective_speed_band", "ordering_consistency_score", "distance_surrogate_shell_slope"}
                            else "frozen"
                        ),
                    }
                )
                series[observable_name][hierarchy_label].append(value)

    with TABLE_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "hierarchy_label",
                "n_side",
                "h",
                "seeds",
                "observable_name",
                "value",
                "source_metric",
                "notes",
            ],
        )
        writer.writeheader()
        writer.writerows(table_rows)

    titles = {
        "dispersion_proxy": "Low-mode Dispersion Proxy",
        "effective_speed_band": "Effective Propagation Speed",
        "persistence_time": "Localized-mode Persistence",
        "ordering_consistency_score": "Temporal-ordering Consistency",
        "distance_surrogate_shell_slope": "Distance-surrogate Shell Slope",
    }
    generate_svg(PLOT_PATH, LEVELS, series, titles)

    h_values = [1.0 / float(level) for level in LEVELS]
    fit_notes: list[str] = []
    smooth_count = 0
    branch_control_separations: list[str] = []
    for observable_name in titles:
        branch_values = series[observable_name][BRANCH]
        control_values = series[observable_name][CONTROL]
        _, _, branch_linear_r2 = linear_fit(h_values, branch_values)
        _, branch_log_slope, branch_log_r2 = log_fit(h_values, branch_values)
        _, _, control_linear_r2 = linear_fit(h_values, control_values)
        _, control_log_slope, control_log_r2 = log_fit(h_values, control_values)
        branch_span = relative_span(branch_values)
        control_span = relative_span(control_values)
        if max(branch_linear_r2, branch_log_r2) >= 0.90:
            smooth_count += 1
        fit_notes.append(
            f"- `{observable_name}`: branch relative span `{branch_span:.6f}`, control relative span `{control_span:.6f}`, "
            f"branch best fit `R^2={max(branch_linear_r2, branch_log_r2):.6f}`, control best fit `R^2={max(control_linear_r2, control_log_r2):.6f}`, "
            f"branch log-slope `{branch_log_slope:.6f}`, control log-slope `{control_log_slope:.6f}`."
        )
        separation = mean(abs(left - right) for left, right in zip(branch_values, control_values))
        branch_control_separations.append(f"- `{observable_name}` mean branch/control gap: `{separation:.6f}`.")

    branch_speed_values = series["effective_speed_band"][BRANCH]
    control_speed_values = series["effective_speed_band"][CONTROL]
    branch_speed_cv = coeff_var(branch_speed_values)
    control_speed_cv = coeff_var(control_speed_values)
    branch_shell_values = series["distance_surrogate_shell_slope"][BRANCH]
    control_shell_values = series["distance_surrogate_shell_slope"][CONTROL]
    branch_shell_drift = max(
        abs(branch_shell_values[index] - branch_shell_values[index - 1]) for index in range(1, len(branch_shell_values))
    )
    control_shell_drift = max(
        abs(control_shell_values[index] - control_shell_values[index - 1]) for index in range(1, len(control_shell_values))
    )

    feasible = (
        smooth_count >= 2
        and branch_speed_cv <= 0.05
        and branch_shell_drift <= 0.10
        and control_speed_cv > branch_speed_cv
        and control_shell_drift > branch_shell_drift
    )

    lines = [
        "# Continuum Sketch Summary",
        "",
        "## Inputs",
        "",
        "- Refinement slice: `n_side = {48, 60, 72, 84}`",
        f"- Branch seeds: `{', '.join(SEEDS)}`",
        f"- Matched control hierarchy: `{CONTROL}`",
        "- Disturbance family for propagation / ordering / shell metrics: `bias_onset`",
        "",
        "## Observable Conventions",
        "",
        "- Low-mode dispersion proxy: frozen `spectral_gap` from the Phase XI persistence ledger.",
        "- Effective propagation speed band: mean `effective_speed` on the Phase XV `bias_onset` / ensemble-7 / seed-pair slice.",
        "- Temporal-ordering consistency score: `1 - event_time_distance / 4.2` on the Phase XVI primary event ledger for the same seed pair.",
        "- Distance-surrogate shell slope: the Phase XVIII least-squares arrival-vs-depth slope, extended to `n_side = 48` using the same artifact-only reconstruction rule.",
        "",
        "## Feasibility Checks",
        "",
        f"- Smooth branch trends identified: `{smooth_count}` of `5` observables.",
        f"- Branch propagation-speed coefficient of variation: `{branch_speed_cv:.6f}` vs control `{control_speed_cv:.6f}`.",
        f"- Branch shell-slope max adjacent drift: `{branch_shell_drift:.6f}` vs control `{control_shell_drift:.6f}`.",
        f"- Branch persistence-time band: `{min(series['persistence_time'][BRANCH]):.6f}` to `{max(series['persistence_time'][BRANCH]):.6f}`.",
        f"- Branch temporal-ordering score band: `{min(series['ordering_consistency_score'][BRANCH]):.6f}` to `{max(series['ordering_consistency_score'][BRANCH]):.6f}`.",
        "",
        "## Fit Notes",
        "",
        *fit_notes,
        "",
        "## Branch vs Control",
        "",
        *branch_control_separations,
        "",
    ]

    if feasible:
        lines.append(
            "The minimal scaling sketch indicates bounded refinement-ordered behavior of key frozen observables. This is consistent with a proto-continuum effective description within the tested branch regime. No continuum ontology or universality is asserted."
        )
    else:
        lines.append(
            "The minimal scaling sketch does not yet show stable refinement-ordered behavior sufficient for a proto-continuum effective description within the tested branch regime."
        )

    SUMMARY_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
