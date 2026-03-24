from __future__ import annotations

import argparse
import csv
import json
import shutil
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np

from telemetry import (
    SurvivalThresholds,
    acyclicity_score,
    causal_depths,
    front_arrival_order,
    order_compatibility,
    persistence_time,
    reconstruct_influence_edges,
)


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_DIR = ROOT / "examples"
SCENARIO_DIR = EXAMPLES_DIR / "scenarios"
OUTPUT_DIR = EXAMPLES_DIR / "output"
ROOT_OVERVIEW_PLOT = ROOT / "stability_demo_overview.svg"

THRESHOLDS = SurvivalThresholds(
    max_width_growth=0.10,
    min_concentration=0.88,
    max_participation_growth=0.20,
    min_overlap=0.90,
    min_recovery_score=0.90,
)

SCENARIO_ORDER = ("baseline", "perturbed", "fragmented")
TECHNICAL_METRICS = ("persistence", "ordering", "depth_drift", "distance_coherence")
PUBLIC_METRIC_ALIASES = {
    "persistence": "structural_retention",
    "ordering": "temporal_consistency",
    "depth_drift": "causal_deformation",
    "distance_coherence": "geometric_integrity",
}
PUBLIC_METRIC_TITLES = {
    "persistence": "Structural Retention",
    "ordering": "Temporal Consistency",
    "depth_drift": "Causal Deformation (lower is better)",
    "distance_coherence": "Geometric Integrity",
}
HUMAN_TABLE_HEADERS = (
    ("scenario", "scenario"),
    ("structural_retention", "retention"),
    ("temporal_consistency", "consistency"),
    ("causal_deformation", "deformation"),
    ("geometric_integrity", "integrity"),
    ("class", "class"),
)
CLASS_COLORS = {"stable": "#0f766e", "marginal": "#b45309", "unstable": "#b91c1c"}
SCAN_AXIS_ALIASES = {
    "noise": "noise",
    "n": "noise",
    "connectivity": "connectivity_drop",
    "connectivity_drop": "connectivity_drop",
    "drop": "connectivity_drop",
}
SCAN_TITLES = {
    "noise": "noise",
    "connectivity_drop": "connectivity drop",
}


def load_json(path: Path) -> dict[str, object]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def clone_json_object(data: dict[str, object]) -> dict[str, object]:
    return json.loads(json.dumps(data))


def clamp(value: float, lower: float = 0.0, upper: float = 1.0) -> float:
    return max(lower, min(upper, value))


def normalize_persistence(scenario: dict[str, object]) -> float:
    tau_grid = [float(value) for value in scenario["tau_grid"]]
    value = persistence_time(list(scenario["persistence_history"]), tau_grid, THRESHOLDS)
    return 0.0 if not tau_grid or tau_grid[-1] <= 0.0 else float(value / tau_grid[-1])


def arrival_windows(scenario: dict[str, object]) -> list[dict[str, int | None]]:
    threshold = float(scenario["ordering_threshold"])
    windows: list[dict[str, int | None]] = []
    for window in scenario["windows"]:
        probe_histories = {
            label: np.asarray(history, dtype=float)
            for label, history in window["probe_histories"].items()
        }
        windows.append(front_arrival_order(probe_histories, threshold))
    return windows


def completeness_fraction(arrival_ordering: dict[str, int | None], nodes: list[str]) -> float:
    arrived = sum(1 for node in nodes if arrival_ordering.get(node) is not None)
    return float(arrived / max(len(nodes), 1))


def normalized_depth_drift(
    reference_depths: dict[str, int],
    current_depths: dict[str, int],
    nodes: list[str],
) -> float:
    max_depth = max(reference_depths.values(), default=0)
    missing_depth = max_depth + 1
    diffs = [
        abs(current_depths.get(node, missing_depth) - reference_depths.get(node, missing_depth))
        for node in nodes
    ]
    return float(sum(diffs) / max(len(nodes) * max(missing_depth, 1), 1))


def shell_overlap_fraction(
    expected_shells: dict[str, int],
    current_depths: dict[str, int],
) -> float:
    max_shell = max(expected_shells.values(), default=0)
    missing_shell = max_shell + 1
    mismatches = sum(
        1
        for node, expected in expected_shells.items()
        if current_depths.get(node, missing_shell) != expected
    )
    return float(mismatches / max(len(expected_shells), 1))


def triangle_violation_rate(
    distance_proxy: dict[str, dict[str, float]],
    nodes: list[str],
    tolerance: float = 1.0e-9,
) -> float:
    triples = list(combinations(nodes, 3))
    if not triples:
        return 0.0

    violations = 0
    for left, middle, right in triples:
        distances = sorted(
            [
                float(distance_proxy[left][middle]),
                float(distance_proxy[left][right]),
                float(distance_proxy[middle][right]),
            ]
        )
        if distances[2] > distances[0] + distances[1] + tolerance:
            violations += 1
    return float(violations / len(triples))


def classify_row(
    persistence: float,
    ordering: float,
    depth_drift: float,
    distance_coherence: float,
) -> str:
    if persistence >= 0.90 and ordering >= 0.90 and depth_drift <= 0.05 and distance_coherence >= 0.90:
        return "stable"
    if persistence >= 0.65 and ordering >= 0.75 and depth_drift <= 0.15 and distance_coherence >= 0.60:
        return "marginal"
    return "unstable"


def evaluate_scenario(scenario: dict[str, object]) -> dict[str, object]:
    nodes = list(scenario["nodes"])
    source = str(scenario["source"])
    expected_shells = {key: int(value) for key, value in scenario["expected_shells"].items()}
    windows = arrival_windows(scenario)

    reference_arrival = windows[0]
    reference_edges = reconstruct_influence_edges(reference_arrival)
    reference_depths = causal_depths(reference_edges, source)

    ordering_scores: list[float] = []
    depth_drifts: list[float] = []
    shell_overlaps: list[float] = []

    for arrival_ordering in windows:
        current_edges = reconstruct_influence_edges(arrival_ordering)
        current_depths = causal_depths(current_edges, source)
        mismatch = order_compatibility(reference_edges, arrival_ordering)
        completeness = completeness_fraction(arrival_ordering, nodes)
        ordering_scores.append(
            float(acyclicity_score(current_edges, nodes) * (1.0 - mismatch) * completeness)
        )
        depth_drifts.append(normalized_depth_drift(reference_depths, current_depths, nodes))
        shell_overlaps.append(shell_overlap_fraction(expected_shells, current_depths))

    persistence = normalize_persistence(scenario)
    ordering = float(sum(ordering_scores) / max(len(ordering_scores), 1))
    depth_drift = float(sum(depth_drifts) / max(len(depth_drifts), 1))
    shell_overlap = float(sum(shell_overlaps) / max(len(shell_overlaps), 1))
    triangle_violation = triangle_violation_rate(scenario["distance_proxy"], nodes)
    distance_coherence = float(max(0.0, 1.0 - 0.5 * (shell_overlap + triangle_violation)))

    row: dict[str, object] = {
        "scenario": str(scenario["name"]),
        "description": str(scenario["description"]),
        "persistence": persistence,
        "ordering": ordering,
        "depth_drift": depth_drift,
        "distance_coherence": distance_coherence,
        "shell_overlap": shell_overlap,
        "triangle_violation_rate": triangle_violation,
        "class": classify_row(persistence, ordering, depth_drift, distance_coherence),
    }
    for technical_name, public_name in PUBLIC_METRIC_ALIASES.items():
        row[public_name] = row[technical_name]
    return row


def scenario_path_from_value(value: str) -> Path:
    candidate = Path(value)
    if candidate.exists():
        return candidate.resolve()

    if not value.endswith(".json"):
        named = SCENARIO_DIR / f"{value}.json"
        if named.exists():
            return named

    named = SCENARIO_DIR / value
    if named.exists():
        return named

    raise FileNotFoundError(f"Unknown scenario reference: {value}")


def load_scenarios(target: str | None) -> tuple[list[dict[str, object]], str]:
    if target is None or target == "all":
        return ([load_json(SCENARIO_DIR / f"{name}.json") for name in SCENARIO_ORDER], "default_bundle")

    path = scenario_path_from_value(target)
    return ([load_json(path)], path.stem)


def set_symmetric_distance(matrix: dict[str, dict[str, float]], left: str, right: str, value: float) -> None:
    matrix[left][right] = value
    matrix[right][left] = value


def apply_generator(
    scenario: dict[str, object],
    noise: float,
    connectivity_drop: float,
    cluster_split: bool,
) -> dict[str, object]:
    if noise <= 0.0 and connectivity_drop <= 0.0 and not cluster_split:
        return scenario

    generated = clone_json_object(scenario)
    modifier_bits: list[str] = []
    if noise > 0.0:
        modifier_bits.append(f"noise{noise:.2f}")
    if connectivity_drop > 0.0:
        modifier_bits.append(f"drop{connectivity_drop:.2f}")
    if cluster_split:
        modifier_bits.append("split")
    modifier_label = "_".join(bit.replace(".", "p") for bit in modifier_bits)

    generated["name"] = f"{generated['name']}_{modifier_label}"
    generated["description"] = (
        f"{generated['description']} Deterministic generator modifiers: {', '.join(modifier_bits)}."
    )

    history = list(generated["persistence_history"])
    for index, record in enumerate(history):
        progress = float(index + 1) / max(len(history), 1)
        split_weight = 0.12 * progress if cluster_split else 0.0
        record["width_growth"] = round(
            float(record["width_growth"]) + noise * 0.20 * progress + connectivity_drop * 0.40 * progress + split_weight,
            6,
        )
        record["concentration_retention"] = round(
            clamp(float(record["concentration_retention"]) - noise * 0.22 * progress - connectivity_drop * 0.35 * progress - split_weight),
            6,
        )
        record["participation_growth"] = round(
            float(record["participation_growth"]) + noise * 0.25 * progress + connectivity_drop * 0.30 * progress + split_weight,
            6,
        )
        record["overlap"] = round(
            clamp(float(record["overlap"]) - noise * 0.30 * progress - connectivity_drop * 0.35 * progress - 0.18 * progress * int(cluster_split)),
            6,
        )
        record["recovery_score"] = round(
            clamp(float(record["recovery_score"]) - noise * 0.34 * progress - connectivity_drop * 0.42 * progress - 0.20 * progress * int(cluster_split)),
            6,
        )

    expected_shells = {key: int(value) for key, value in generated["expected_shells"].items()}
    max_depth = max(expected_shells.values(), default=1)
    for window_index, window in enumerate(generated["windows"]):
        probe_histories = window["probe_histories"]
        for node, values in probe_histories.items():
            if node == generated["source"]:
                continue
            depth = expected_shells.get(node, 1)
            shift = int(round(connectivity_drop * depth * 2.0))
            if cluster_split and node in {"mid", "far"}:
                shift += 1
            attenuation = noise * (0.10 + 0.05 * depth) + connectivity_drop * (0.10 * depth)
            if cluster_split and node == "mid":
                attenuation += 0.08
            if cluster_split and node == "far":
                attenuation += 0.16

            adjusted: list[float] = []
            for sample_index, value in enumerate(values):
                progress = float(sample_index + 1) / max(len(values), 1)
                wobble = noise * 0.03 * ((window_index + sample_index + depth) % 3)
                degraded = clamp(float(value) - attenuation * progress - wobble)
                adjusted.append(round(degraded, 6))

            if shift > 0:
                adjusted = ([0.0] * shift + adjusted)[: len(adjusted)]

            if cluster_split and node == "mid" and adjusted:
                adjusted = ([0.0] + adjusted)[: len(adjusted)]

            if cluster_split and node == "far" and len(adjusted) >= 3:
                adjusted[1] = round(max(adjusted[1], clamp(0.74 - 0.05 * window_index)), 6)
                adjusted[2] = round(max(adjusted[2], clamp(0.86 - 0.04 * window_index)), 6)
            probe_histories[node] = adjusted

    distance_proxy = generated["distance_proxy"]
    base_scale = 1.0 + noise * 0.10 + connectivity_drop * 0.20
    nodes = list(generated["nodes"])
    for left, right in combinations(nodes, 2):
        current = float(distance_proxy[left][right]) * base_scale
        depth_gap = abs(expected_shells.get(left, 0) - expected_shells.get(right, 0))
        current += connectivity_drop * 0.20 * depth_gap / max(max_depth, 1)
        set_symmetric_distance(distance_proxy, left, right, round(current, 6))

    if connectivity_drop > 0.0:
        set_symmetric_distance(
            distance_proxy,
            "source",
            "far",
            round(float(distance_proxy["source"]["far"]) + 0.80 * connectivity_drop, 6),
        )
        set_symmetric_distance(
            distance_proxy,
            "near",
            "mid",
            round(max(0.2, float(distance_proxy["near"]["mid"]) - 0.30 * connectivity_drop), 6),
        )

    if cluster_split:
        set_symmetric_distance(distance_proxy, "source", "far", round(float(distance_proxy["source"]["far"]) + 0.90, 6))
        set_symmetric_distance(distance_proxy, "near", "far", round(max(0.2, float(distance_proxy["near"]["far"]) - 0.95), 6))
        set_symmetric_distance(distance_proxy, "mid", "far", round(float(distance_proxy["mid"]["far"]) + 0.35, 6))

    return generated


def output_base_name(bundle_label: str, noise: float, connectivity_drop: float, cluster_split: bool) -> str:
    if bundle_label == "default_bundle" and noise <= 0.0 and connectivity_drop <= 0.0 and not cluster_split:
        return "stability_demo"

    safe_label = bundle_label.replace(".", "_").replace("-", "_")
    parts = [f"stability_{safe_label}"]
    if noise > 0.0:
        parts.append(f"noise{noise:.2f}".replace(".", "p"))
    if connectivity_drop > 0.0:
        parts.append(f"drop{connectivity_drop:.2f}".replace(".", "p"))
    if cluster_split:
        parts.append("split")
    return "_".join(parts)


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "scenario",
        "description",
        "noise",
        "connectivity_drop",
        "cluster_split",
        "structural_retention",
        "temporal_consistency",
        "causal_deformation",
        "geometric_integrity",
        "persistence",
        "ordering",
        "depth_drift",
        "distance_coherence",
        "shell_overlap",
        "triangle_violation_rate",
        "class",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "scenario": row["scenario"],
                    "description": row["description"],
                    "noise": "" if "noise" not in row else f"{float(row['noise']):.3f}",
                    "connectivity_drop": "" if "connectivity_drop" not in row else f"{float(row['connectivity_drop']):.3f}",
                    "cluster_split": "" if "cluster_split" not in row else str(bool(row["cluster_split"])).lower(),
                    "structural_retention": f"{float(row['structural_retention']):.3f}",
                    "temporal_consistency": f"{float(row['temporal_consistency']):.3f}",
                    "causal_deformation": f"{float(row['causal_deformation']):.3f}",
                    "geometric_integrity": f"{float(row['geometric_integrity']):.3f}",
                    "persistence": f"{float(row['persistence']):.3f}",
                    "ordering": f"{float(row['ordering']):.3f}",
                    "depth_drift": f"{float(row['depth_drift']):.3f}",
                    "distance_coherence": f"{float(row['distance_coherence']):.3f}",
                    "shell_overlap": f"{float(row['shell_overlap']):.3f}",
                    "triangle_violation_rate": f"{float(row['triangle_violation_rate']):.3f}",
                    "class": row["class"],
                }
            )


def render_metric_plot(path: Path, rows: list[dict[str, object]]) -> None:
    width = 1100
    height = 760
    panel_width = 470
    panel_height = 250
    margin_x = 60
    margin_y = 70
    gap_x = 40
    gap_y = 45
    scenario_labels = [str(row["scenario"]) for row in rows]
    threshold_lines = {
        "persistence": (0.90, 0.65),
        "ordering": (0.90, 0.75),
        "depth_drift": (0.05, 0.15),
        "distance_coherence": (0.90, 0.60),
    }

    def panel_origin(index_id: int) -> tuple[float, float]:
        row = index_id // 2
        col = index_id % 2
        return (
            margin_x + col * (panel_width + gap_x),
            margin_y + row * (panel_height + gap_y),
        )

    def y_position(top: float, bottom: float, value: float) -> float:
        return bottom - value * (bottom - top)

    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='white'/>",
        "<text x='550' y='32' text-anchor='middle' font-size='22' font-family='Helvetica, Arial, sans-serif'>Trajectory Stability Demo</text>",
        "<text x='550' y='54' text-anchor='middle' font-size='12' fill='#4b5563' font-family='Helvetica, Arial, sans-serif'>Synthetic scenarios evaluated with the frozen telemetry layer plus a lightweight distance proxy</text>",
    ]

    for index_id, metric_name in enumerate(TECHNICAL_METRICS):
        x0, y0 = panel_origin(index_id)
        inner_left = x0 + 55
        inner_right = x0 + panel_width - 20
        inner_top = y0 + 35
        inner_bottom = y0 + panel_height - 35
        bar_width = 70
        bar_gap = 45
        stable_line, marginal_line = threshold_lines[metric_name]

        parts.append(f"<rect x='{x0}' y='{y0}' width='{panel_width}' height='{panel_height}' fill='none' stroke='#d1d5db' stroke-width='1'/>")
        parts.append(
            f"<text x='{x0 + panel_width / 2:.1f}' y='{y0 + 20}' text-anchor='middle' font-size='14' font-family='Helvetica, Arial, sans-serif'>{PUBLIC_METRIC_TITLES[metric_name]}</text>"
        )
        parts.append(f"<line x1='{inner_left}' y1='{inner_bottom}' x2='{inner_right}' y2='{inner_bottom}' stroke='#374151' stroke-width='1.2'/>")
        parts.append(f"<line x1='{inner_left}' y1='{inner_top}' x2='{inner_left}' y2='{inner_bottom}' stroke='#374151' stroke-width='1.2'/>")

        for tick in range(6):
            value = tick / 5.0
            ypos = y_position(inner_top, inner_bottom, value)
            parts.append(f"<line x1='{inner_left}' y1='{ypos:.2f}' x2='{inner_right}' y2='{ypos:.2f}' stroke='#f3f4f6' stroke-width='1'/>")
            parts.append(f"<text x='{inner_left - 8}' y='{ypos + 4:.2f}' text-anchor='end' font-size='10' font-family='Helvetica, Arial, sans-serif'>{value:.1f}</text>")

        for value, color, label in (
            (stable_line, "#0f766e", "stable line"),
            (marginal_line, "#b45309", "marginal line"),
        ):
            ypos = y_position(inner_top, inner_bottom, value)
            parts.append(f"<line x1='{inner_left}' y1='{ypos:.2f}' x2='{inner_right}' y2='{ypos:.2f}' stroke='{color}' stroke-dasharray='5 4' stroke-width='1.2'/>")
            parts.append(f"<text x='{inner_right - 4}' y='{ypos - 4:.2f}' text-anchor='end' font-size='9' fill='{color}' font-family='Helvetica, Arial, sans-serif'>{label}</text>")

        start_x = inner_left + 25
        for bar_index, row in enumerate(rows):
            xpos = start_x + bar_index * (bar_width + bar_gap)
            value = float(row[metric_name])
            color = CLASS_COLORS[str(row["class"])]
            bar_top = y_position(inner_top, inner_bottom, value)
            parts.append(f"<rect x='{xpos}' y='{bar_top:.2f}' width='{bar_width}' height='{inner_bottom - bar_top:.2f}' fill='{color}' opacity='0.88'/>")
            parts.append(f"<text x='{xpos + bar_width / 2:.2f}' y='{bar_top - 6:.2f}' text-anchor='middle' font-size='10' font-family='Helvetica, Arial, sans-serif'>{value:.2f}</text>")
            parts.append(f"<text x='{xpos + bar_width / 2:.2f}' y='{inner_bottom + 18}' text-anchor='middle' font-size='10' font-family='Helvetica, Arial, sans-serif'>{scenario_labels[bar_index]}</text>")

    legend_y = height - 28
    legend_x = 260
    for offset, label in enumerate(("stable", "marginal", "unstable")):
        lx = legend_x + offset * 180
        parts.append(f"<rect x='{lx}' y='{legend_y - 12}' width='16' height='16' fill='{CLASS_COLORS[label]}'/>")
        parts.append(f"<text x='{lx + 24}' y='{legend_y}' font-size='11' font-family='Helvetica, Arial, sans-serif'>{label}</text>")

    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def classification_score(label: str) -> int:
    return {"unstable": 0, "marginal": 1, "stable": 2}[label]


def render_scan_heatmap(
    path: Path,
    rows: list[dict[str, object]],
    noise_values: list[float],
    connectivity_values: list[float],
    scenario_label: str,
    cluster_split: bool,
) -> None:
    cell_w = 84
    cell_h = 52
    left_margin = 110
    top_margin = 88
    width = left_margin + len(connectivity_values) * cell_w + 60
    height = top_margin + len(noise_values) * cell_h + 110
    lookup = {
        (float(row["noise"]), float(row["connectivity_drop"])): row
        for row in rows
    }

    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='white'/>",
        "<text x='50%' y='30' text-anchor='middle' font-size='22' font-family='Helvetica, Arial, sans-serif'>Stability Scan Heatmap</text>",
        f"<text x='50%' y='52' text-anchor='middle' font-size='12' fill='#4b5563' font-family='Helvetica, Arial, sans-serif'>scenario={scenario_label} | cluster_split={'on' if cluster_split else 'off'} | cell text=classification</text>",
        f"<text x='50%' y='{height - 20}' text-anchor='middle' font-size='12' fill='#4b5563' font-family='Helvetica, Arial, sans-serif'>x axis: connectivity drop | y axis: noise</text>",
    ]

    for x_index, connectivity in enumerate(connectivity_values):
        xpos = left_margin + x_index * cell_w + cell_w / 2
        parts.append(f"<text x='{xpos:.2f}' y='74' text-anchor='middle' font-size='10' font-family='Helvetica, Arial, sans-serif'>{connectivity:.2f}</text>")

    for y_index, noise in enumerate(noise_values):
        ypos = top_margin + y_index * cell_h + cell_h / 2 + 4
        parts.append(f"<text x='94' y='{ypos:.2f}' text-anchor='end' font-size='10' font-family='Helvetica, Arial, sans-serif'>{noise:.2f}</text>")

    parts.append(f"<text x='22' y='{top_margin + (len(noise_values) * cell_h) / 2:.2f}' transform='rotate(-90 22 {top_margin + (len(noise_values) * cell_h) / 2:.2f})' text-anchor='middle' font-size='11' font-family='Helvetica, Arial, sans-serif'>noise</text>")
    parts.append(f"<text x='{left_margin + (len(connectivity_values) * cell_w) / 2:.2f}' y='96' text-anchor='middle' font-size='11' font-family='Helvetica, Arial, sans-serif'>connectivity drop</text>")

    for y_index, noise in enumerate(noise_values):
        for x_index, connectivity in enumerate(connectivity_values):
            row = lookup[(noise, connectivity)]
            xpos = left_margin + x_index * cell_w
            ypos = top_margin + y_index * cell_h
            color = CLASS_COLORS[str(row["class"])]
            parts.append(f"<rect x='{xpos}' y='{ypos}' width='{cell_w}' height='{cell_h}' fill='{color}' opacity='0.85' stroke='white' stroke-width='1'/>")
            parts.append(f"<text x='{xpos + cell_w / 2:.2f}' y='{ypos + 23:.2f}' text-anchor='middle' font-size='10' fill='white' font-family='Helvetica, Arial, sans-serif'>{row['class']}</text>")
            parts.append(f"<text x='{xpos + cell_w / 2:.2f}' y='{ypos + 38:.2f}' text-anchor='middle' font-size='9' fill='white' font-family='Helvetica, Arial, sans-serif'>score={classification_score(str(row['class']))}</text>")

    legend_x = left_margin
    legend_y = height - 52
    for offset, label in enumerate(("stable", "marginal", "unstable")):
        lx = legend_x + offset * 170
        parts.append(f"<rect x='{lx}' y='{legend_y}' width='16' height='16' fill='{CLASS_COLORS[label]}'/>")
        parts.append(f"<text x='{lx + 24}' y='{legend_y + 12}' font-size='11' font-family='Helvetica, Arial, sans-serif'>{label}</text>")

    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def format_table(rows: list[dict[str, object]]) -> str:
    body = [
        (
            str(row["scenario"]),
            f"{float(row['structural_retention']):.3f}",
            f"{float(row['temporal_consistency']):.3f}",
            f"{float(row['causal_deformation']):.3f}",
            f"{float(row['geometric_integrity']):.3f}",
            str(row["class"]),
        )
        for row in rows
    ]
    widths = [
        max(len(header), max(len(record[index]) for record in body))
        for index, (_, header) in enumerate(HUMAN_TABLE_HEADERS)
    ]

    lines = [
        "  ".join(HUMAN_TABLE_HEADERS[index][1].ljust(widths[index]) for index in range(len(HUMAN_TABLE_HEADERS))),
        "  ".join("-" * widths[index] for index in range(len(HUMAN_TABLE_HEADERS))),
    ]
    for record in body:
        lines.append("  ".join(record[index].ljust(widths[index]) for index in range(len(record))))
    return "\n".join(lines)


def parse_range_spec(spec: str) -> list[float]:
    parts = [segment.strip() for segment in spec.split(":")]
    if len(parts) != 3:
        raise ValueError(f"Range spec must be start:stop:step, got: {spec}")
    start, stop, step = (float(part) for part in parts)
    if step <= 0.0:
        raise ValueError(f"Range step must be positive, got: {spec}")
    values: list[float] = []
    current = start
    while current <= stop + 1.0e-12:
        values.append(round(current, 10))
        current += step
    return values


def parse_scan_specs(tokens: list[str]) -> dict[str, list[float]]:
    axes: dict[str, list[float]] = {}
    for token in tokens:
        if "=" not in token:
            raise ValueError(f"Scan token must have axis=value form, got: {token}")
        name, spec = token.split("=", 1)
        canonical_name = SCAN_AXIS_ALIASES.get(name.strip().lower())
        if canonical_name is None:
            valid = ", ".join(sorted(SCAN_AXIS_ALIASES))
            raise ValueError(f"Unknown scan axis '{name}'. Valid aliases: {valid}")
        axes[canonical_name] = parse_range_spec(spec.strip())
    return axes


def round_metric(value: float) -> float:
    return round(float(value), 6)


def row_json_payload(row: dict[str, object]) -> dict[str, object]:
    payload = {
        "scenario": row["scenario"],
        "classification": row["class"],
        "structural_retention": round_metric(float(row["structural_retention"])),
        "temporal_consistency": round_metric(float(row["temporal_consistency"])),
        "causal_deformation": round_metric(float(row["causal_deformation"])),
        "geometric_integrity": round_metric(float(row["geometric_integrity"])),
        "metrics": {
            "persistence": round_metric(float(row["persistence"])),
            "ordering": round_metric(float(row["ordering"])),
            "depth_drift": round_metric(float(row["depth_drift"])),
            "distance_coherence": round_metric(float(row["distance_coherence"])),
            "shell_overlap": round_metric(float(row["shell_overlap"])),
            "triangle_violation_rate": round_metric(float(row["triangle_violation_rate"])),
        },
    }
    if "description" in row:
        payload["description"] = row["description"]
    if "noise" in row:
        payload["noise"] = round_metric(float(row["noise"]))
    if "connectivity_drop" in row:
        payload["connectivity_drop"] = round_metric(float(row["connectivity_drop"]))
    if "cluster_split" in row:
        payload["cluster_split"] = bool(row["cluster_split"])
    return payload


def ordered_classification_counts(rows: list[dict[str, object]]) -> dict[str, int]:
    counts = Counter(str(row["class"]) for row in rows)
    return {label: int(counts.get(label, 0)) for label in ("stable", "marginal", "unstable")}


def single_run_json_payload(rows: list[dict[str, object]], csv_path: Path, plot_path: Path, overview_path: Path | None) -> dict[str, object]:
    if len(rows) == 1:
        payload = row_json_payload(rows[0])
        payload["artifacts"] = {"csv": str(csv_path), "plot": str(plot_path)}
        if overview_path is not None:
            payload["artifacts"]["overview"] = str(overview_path)
        return payload

    payload = {
        "mode": "stability_bundle",
        "classification_counts": ordered_classification_counts(rows),
        "scenarios": [row_json_payload(row) for row in rows],
        "artifacts": {"csv": str(csv_path), "plot": str(plot_path)},
    }
    if overview_path is not None:
        payload["artifacts"]["overview"] = str(overview_path)
    return payload


def scan_json_payload(
    rows: list[dict[str, object]],
    scenario_label: str,
    noise_values: list[float],
    connectivity_values: list[float],
    csv_path: Path,
    heatmap_path: Path,
    cluster_split: bool,
) -> dict[str, object]:
    return {
        "mode": "stability_scan",
        "scenario": scenario_label,
        "cluster_split": cluster_split,
        "scan_axes": {
            "noise": [round_metric(value) for value in noise_values],
            "connectivity_drop": [round_metric(value) for value in connectivity_values],
        },
        "classification_counts": ordered_classification_counts(rows),
        "results": [row_json_payload(row) for row in rows],
        "artifacts": {"csv": str(csv_path), "heatmap": str(heatmap_path)},
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Public HAOS-IIP demo utilities.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    stability = subparsers.add_parser(
        "stability",
        help="Classify synthetic propagation trajectories with the frozen telemetry layer.",
    )
    stability.add_argument(
        "scenario",
        nargs="?",
        default="all",
        help="Scenario name, scenario JSON path, or `all` for the default bundle.",
    )
    stability.add_argument(
        "--noise",
        type=float,
        default=0.0,
        help="Deterministic amplitude-style degradation applied to the selected synthetic scenario.",
    )
    stability.add_argument(
        "--connectivity-drop",
        type=float,
        default=0.0,
        help="Deterministic propagation-delay and distance distortion strength.",
    )
    stability.add_argument(
        "--cluster-split",
        action="store_true",
        help="Apply a deterministic cluster split that isolates the deeper shell structure.",
    )
    stability.add_argument(
        "--json",
        action="store_true",
        help="Emit machine-readable JSON instead of the human-readable summary.",
    )
    stability.add_argument(
        "--scan",
        nargs="+",
        metavar="AXIS=START:STOP:STEP",
        help="Run a deterministic parameter scan, for example: --scan noise=0.0:0.2:0.02 connectivity=0.0:0.4:0.05",
    )
    return parser


def run_single_stability(args: argparse.Namespace) -> int:
    noise = max(float(args.noise), 0.0)
    connectivity_drop = max(float(args.connectivity_drop), 0.0)

    scenarios, bundle_label = load_scenarios(args.scenario)
    rows = [
        evaluate_scenario(
            apply_generator(
                scenario,
                noise=noise,
                connectivity_drop=connectivity_drop,
                cluster_split=bool(args.cluster_split),
            )
        )
        for scenario in scenarios
    ]

    base_name = output_base_name(
        bundle_label,
        noise=noise,
        connectivity_drop=connectivity_drop,
        cluster_split=bool(args.cluster_split),
    )
    csv_path = OUTPUT_DIR / f"{base_name}_results.csv"
    plot_path = OUTPUT_DIR / f"{base_name}_metrics.svg"
    write_csv(csv_path, rows)
    render_metric_plot(plot_path, rows)

    overview_path: Path | None = None
    if bundle_label == "default_bundle" and noise <= 0.0 and connectivity_drop <= 0.0 and not args.cluster_split:
        shutil.copyfile(plot_path, ROOT_OVERVIEW_PLOT)
        overview_path = ROOT_OVERVIEW_PLOT

    if args.json:
        print(json.dumps(single_run_json_payload(rows, csv_path, plot_path, overview_path), indent=2))
        return 0

    print("Trajectory stability demo completed.")
    print(f"CSV: {csv_path}")
    print(f"Plot: {plot_path}")
    if overview_path is not None:
        print(f"Overview: {overview_path}")
    print()
    print(format_table(rows))
    print()
    print(
        "This demo illustrates how HAOS-IIP telemetry can be used as a structural-stability diagnostic for evolving "
        "graph-like systems. It is not a claim about physical spacetime emergence."
    )
    print(
        "The classification relies only on previously established persistence, ordering, causal-depth, and "
        "proto-distance diagnostics. No claims about continuum limits or physical ontology are implied."
    )
    return 0


def run_scan(args: argparse.Namespace) -> int:
    scan_axes = parse_scan_specs(list(args.scan or []))
    noise_values = scan_axes.get("noise", [max(float(args.noise), 0.0)])
    connectivity_values = scan_axes.get("connectivity_drop", [max(float(args.connectivity_drop), 0.0)])

    scenario_target = args.scenario if args.scenario != "all" else "baseline"
    scenarios, bundle_label = load_scenarios(scenario_target)
    if len(scenarios) != 1:
        raise ValueError("Scan mode requires exactly one scenario.")
    base_scenario = scenarios[0]

    rows: list[dict[str, object]] = []
    for noise in noise_values:
        for connectivity_drop in connectivity_values:
            evaluated = evaluate_scenario(
                apply_generator(
                    base_scenario,
                    noise=float(noise),
                    connectivity_drop=float(connectivity_drop),
                    cluster_split=bool(args.cluster_split),
                )
            )
            evaluated["noise"] = float(noise)
            evaluated["connectivity_drop"] = float(connectivity_drop)
            evaluated["cluster_split"] = bool(args.cluster_split)
            rows.append(evaluated)

    base_name = output_base_name(
        bundle_label,
        noise=0.0,
        connectivity_drop=0.0,
        cluster_split=bool(args.cluster_split),
    )
    csv_path = OUTPUT_DIR / f"{base_name}_scan_grid.csv"
    heatmap_path = OUTPUT_DIR / f"{base_name}_scan_heatmap.svg"
    write_csv(csv_path, rows)
    render_scan_heatmap(
        heatmap_path,
        rows,
        noise_values=noise_values,
        connectivity_values=connectivity_values,
        scenario_label=str(base_scenario["name"]),
        cluster_split=bool(args.cluster_split),
    )

    if args.json:
        print(
            json.dumps(
                scan_json_payload(
                    rows,
                    scenario_label=str(base_scenario["name"]),
                    noise_values=noise_values,
                    connectivity_values=connectivity_values,
                    csv_path=csv_path,
                    heatmap_path=heatmap_path,
                    cluster_split=bool(args.cluster_split),
                ),
                indent=2,
            )
        )
        return 0

    counts = Counter(str(row["class"]) for row in rows)
    print("Trajectory stability scan completed.")
    print(f"Grid CSV: {csv_path}")
    print(f"Heatmap: {heatmap_path}")
    print(
        "Classification counts: "
        + ", ".join(f"{label}={counts.get(label, 0)}" for label in ("stable", "marginal", "unstable"))
    )
    print(
        f"Scan axes: noise={noise_values[0]:.2f}:{noise_values[-1]:.2f}:{(noise_values[1] - noise_values[0]) if len(noise_values) > 1 else 0.0:.2f}, "
        f"connectivity_drop={connectivity_values[0]:.2f}:{connectivity_values[-1]:.2f}:{(connectivity_values[1] - connectivity_values[0]) if len(connectivity_values) > 1 else 0.0:.2f}"
    )
    return 0


def run_stability_command(args: argparse.Namespace) -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if args.scan:
        return run_scan(args)
    return run_single_stability(args)


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "stability":
        return run_stability_command(args)

    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
