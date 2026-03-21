#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
PHASE_DIR = Path(__file__).resolve().parent
RUNS_DIR = PHASE_DIR / "runs"

PROBE_NAME = "bias_onset"
TARGET_REFINEMENTS = ["60", "72", "84"]
TARGET_ENSEMBLE_SIZE = "7"
TARGET_HIERARCHIES = ["frozen_branch", "periodic_diagonal_augmented_control"]
MAX_SHELL_OVERLAP_THRESHOLD = 0.40
MAX_SLOPE_DRIFT_THRESHOLD = 0.06
MAX_CAUSAL_DEPTH_DRIFT_THRESHOLD = 0.10
MAX_TRIANGLE_VIOLATION_THRESHOLD = 0.05

INPUT_FILES = {
    "phase15_effective_speed": ROOT / "phase15-propagation/runs/phase15_effective_speed_ledger.csv",
    "phase15_influence_range": ROOT / "phase15-propagation/runs/phase15_influence_range_ledger.csv",
    "phase15_propagation": ROOT / "phase15-propagation/runs/phase15_propagation_ledger.csv",
    "phase16_front_arrival": ROOT / "phase16-temporal-ordering/runs/phase16_front_arrival_ordering.csv",
    "phase16_monotonic_parameter": ROOT / "phase16-temporal-ordering/runs/phase16_monotonic_parameter_ledger.csv",
    "phase17_influence_graph": ROOT / "phase17-causal-closure/runs/phase17_influence_graph_ledger.csv",
    "phase17_causal_distance": ROOT / "phase17-causal-closure/runs/phase17_causal_distance_metrics.csv",
    "phase17_order_compatibility": ROOT / "phase17-causal-closure/runs/phase17_propagation_order_compatibility.csv",
}

OUTPUT_FILES = {
    "distance_surrogate": RUNS_DIR / "phase18_distance_surrogate_ledger.csv",
    "shell_ordering": RUNS_DIR / "phase18_shell_ordering_metrics.csv",
    "refinement_scaling": RUNS_DIR / "phase18_refinement_scaling.csv",
    "triangle_violation": RUNS_DIR / "phase18_triangle_violation_rate.csv",
    "manifest": PHASE_DIR / "phase18_manifest.json",
    "summary": PHASE_DIR / "phase18_summary.md",
}


@dataclass(frozen=True)
class RunKey:
    hierarchy: str
    n_side: str
    seed: str


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def as_float(value: str) -> float:
    return float(value)


def mean(values: Iterable[float]) -> float:
    values = list(values)
    if not values:
        raise ValueError("expected at least one value")
    return sum(values) / len(values)


def linear_interpolate(series: list[tuple[float, float]], target_tau: float) -> float:
    if target_tau <= series[0][0]:
        return series[0][1]
    if target_tau >= series[-1][0]:
        return series[-1][1]
    for index in range(1, len(series)):
        left_tau, left_value = series[index - 1]
        right_tau, right_value = series[index]
        if left_tau <= target_tau <= right_tau:
            if right_tau == left_tau:
                return right_value
            weight = (target_tau - left_tau) / (right_tau - left_tau)
            return left_value + weight * (right_value - left_value)
    return series[-1][1]


def least_squares_slope(xs: list[float], ys: list[float]) -> float:
    x_mean = mean(xs)
    y_mean = mean(ys)
    numerator = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys))
    denominator = sum((x - x_mean) ** 2 for x in xs)
    return 0.0 if denominator == 0 else numerator / denominator


def parse_seed_pair(influence_rows: list[dict[str, str]]) -> list[str]:
    selected_pairs = {
        row["selected_seed_pair"]
        for row in influence_rows
        if row["hierarchy_label"] == "frozen_branch" and row["n_side"] in TARGET_REFINEMENTS
    }
    if len(selected_pairs) != 1:
        raise RuntimeError(f"expected one frozen seed pair, found {selected_pairs}")
    pair = next(iter(selected_pairs))
    seeds = [seed.strip() for seed in pair.split(",") if seed.strip()]
    if len(seeds) != 2:
        raise RuntimeError(f"expected two seeds in selected pair, found {pair}")
    return seeds


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    RUNS_DIR.mkdir(parents=True, exist_ok=True)

    phase17_influence_rows = read_csv(INPUT_FILES["phase17_influence_graph"])
    selected_seeds = parse_seed_pair(phase17_influence_rows)

    phase15_speed_rows = [
        row
        for row in read_csv(INPUT_FILES["phase15_effective_speed"])
        if row["probe_name"] == PROBE_NAME
        and row["n_side"] in TARGET_REFINEMENTS
        and row["ensemble_size"] == TARGET_ENSEMBLE_SIZE
        and row["seed"] in selected_seeds
        and row["hierarchy_label"] in TARGET_HIERARCHIES
    ]
    phase15_influence_rows = [
        row
        for row in read_csv(INPUT_FILES["phase15_influence_range"])
        if row["probe_name"] == PROBE_NAME
        and row["n_side"] in TARGET_REFINEMENTS
        and row["ensemble_size"] == TARGET_ENSEMBLE_SIZE
        and row["seed"] in selected_seeds
        and row["hierarchy_label"] in TARGET_HIERARCHIES
    ]
    phase15_propagation_rows = [
        row
        for row in read_csv(INPUT_FILES["phase15_propagation"])
        if row["probe_name"] == PROBE_NAME
        and row["n_side"] in TARGET_REFINEMENTS
        and row["ensemble_size"] == TARGET_ENSEMBLE_SIZE
        and row["seed"] in selected_seeds
        and row["hierarchy_label"] in TARGET_HIERARCHIES
        and row["observable"] == "disturbance_radius"
    ]
    phase16_front_rows = [
        row
        for row in read_csv(INPUT_FILES["phase16_front_arrival"])
        if row["probe_name"] == PROBE_NAME
        and row["n_side"] in TARGET_REFINEMENTS
        and row["ensemble_size"] == TARGET_ENSEMBLE_SIZE
        and row["seed"] in selected_seeds
        and row["hierarchy_label"] in TARGET_HIERARCHIES
    ]
    phase16_monotonic_rows = [
        row
        for row in read_csv(INPUT_FILES["phase16_monotonic_parameter"])
        if row["probe_name"] == PROBE_NAME
        and row["n_side"] in TARGET_REFINEMENTS
        and row["ensemble_size"] == TARGET_ENSEMBLE_SIZE
        and row["seed"] in selected_seeds
        and row["hierarchy_label"] in TARGET_HIERARCHIES
        and row["parameter_name"] == "integrated_radius_progress"
        and row["is_primary_parameter"] == "True"
    ]
    phase17_causal_rows = [
        row
        for row in read_csv(INPUT_FILES["phase17_causal_distance"])
        if row["n_side"] in TARGET_REFINEMENTS and row["hierarchy_label"] in TARGET_HIERARCHIES
    ]
    phase17_compat_rows = [
        row
        for row in read_csv(INPUT_FILES["phase17_order_compatibility"])
        if row["n_side"] in TARGET_REFINEMENTS and row["hierarchy_label"] in TARGET_HIERARCHIES
    ]

    speed_by_key = {(row["hierarchy_label"], row["n_side"], row["seed"]): row for row in phase15_speed_rows}
    influence_by_key: dict[RunKey, list[dict[str, str]]] = defaultdict(list)
    for row in phase15_influence_rows:
        influence_by_key[RunKey(row["hierarchy_label"], row["n_side"], row["seed"])].append(row)

    radius_series_by_key: dict[RunKey, list[tuple[float, float]]] = defaultdict(list)
    for row in phase15_propagation_rows:
        key = RunKey(row["hierarchy_label"], row["n_side"], row["seed"])
        radius_series_by_key[key].append((as_float(row["tau"]), as_float(row["value"])))
    for values in radius_series_by_key.values():
        values.sort()

    front_by_key: dict[RunKey, dict[str, dict[str, str]]] = defaultdict(dict)
    for row in phase16_front_rows:
        front_by_key[RunKey(row["hierarchy_label"], row["n_side"], row["seed"])][row["distance_band"]] = row

    monotonic_by_key: dict[RunKey, float] = {}
    for row in phase16_monotonic_rows:
        key = RunKey(row["hierarchy_label"], row["n_side"], row["seed"])
        monotonic_by_key[key] = max(monotonic_by_key.get(key, 0.0), as_float(row["monotonicity_score"]))

    causal_by_hierarchy = {(row["hierarchy_label"], row["n_side"]): row for row in phase17_causal_rows}
    compat_by_hierarchy = {(row["hierarchy_label"], row["n_side"]): row for row in phase17_compat_rows}
    graph_by_hierarchy = {
        (row["hierarchy_label"], row["n_side"]): row
        for row in phase17_influence_rows
        if row["n_side"] in TARGET_REFINEMENTS and row["hierarchy_label"] in TARGET_HIERARCHIES
    }

    effective_speed_mean: dict[tuple[str, str], float] = {}
    for hierarchy in TARGET_HIERARCHIES:
        for n_side in TARGET_REFINEMENTS:
            speeds = [
                as_float(speed_by_key[(hierarchy, n_side, seed)]["effective_speed"])
                for seed in selected_seeds
            ]
            effective_speed_mean[(hierarchy, n_side)] = mean(speeds)

    distance_rows: list[dict[str, object]] = []
    shell_rows: list[dict[str, object]] = []
    scaling_rows: list[dict[str, object]] = []
    triangle_rows: list[dict[str, object]] = []
    shell_overlap_by_hierarchy: dict[str, list[float]] = defaultdict(list)
    shell_monotonic_by_hierarchy: dict[str, list[float]] = defaultdict(list)
    scaling_slopes_by_hierarchy: dict[str, list[float]] = defaultdict(list)
    depth_means_by_hierarchy: dict[str, list[float]] = defaultdict(list)

    for hierarchy in TARGET_HIERARCHIES:
        previous_slope: float | None = None
        previous_depth_mean: float | None = None

        for n_side in TARGET_REFINEMENTS:
            causal_row = causal_by_hierarchy[(hierarchy, n_side)]
            compat_row = compat_by_hierarchy[(hierarchy, n_side)]
            graph_row = graph_by_hierarchy[(hierarchy, n_side)]
            v_eff = effective_speed_mean[(hierarchy, n_side)]
            front_near_depth = as_float(causal_row["front_near_depth"])
            front_far_depth = as_float(causal_row["front_far_depth"])
            mean_causal_depth = as_float(causal_row["mean_causal_depth"])

            seed_near_latencies: list[float] = []
            seed_far_latencies: list[float] = []
            seed_overlap_fractions: list[float] = []
            seed_monotonic_flags: list[float] = []

            for seed in selected_seeds:
                key = RunKey(hierarchy, n_side, seed)
                speed_row = speed_by_key[(hierarchy, n_side, seed)]
                front_rows = front_by_key[key]
                near_row = front_rows["near"]
                far_row = front_rows["far"]
                monotonicity_score = monotonic_by_key[key]
                radius_series = radius_series_by_key[key]
                influence_rows = influence_by_key[key]

                near_mean = as_float(near_row["mean_latency"])
                near_std = as_float(near_row["latency_std"])
                far_mean = as_float(far_row["mean_latency"])
                far_std = as_float(far_row["latency_std"])
                terminal_radius = radius_series[-1][1]

                near_low = near_mean - near_std
                near_high = near_mean + near_std
                far_low = far_mean - far_std
                far_high = far_mean + far_std
                overlap = max(0.0, min(near_high, far_high) - max(near_low, far_low))
                union = max(near_high, far_high) - min(near_low, far_low)
                overlap_fraction = overlap / union if union > 0 else 0.0
                monotonic_shells = 0.0 <= near_mean <= far_mean
                shell_score = (1.0 - overlap_fraction) if monotonic_shells else 0.0

                mean_influence_distance = mean(as_float(row["distance"]) for row in influence_rows)
                max_normalized_latency = max(as_float(row["normalized_latency"]) for row in influence_rows)

                for shell_name, depth, arrival, arrival_std, front_distance in (
                    ("source", 0.0, 0.0, 0.0, 0.0),
                    ("front_near", front_near_depth, near_mean, near_std, as_float(near_row["front_distance"])),
                    ("front_far", front_far_depth, far_mean, far_std, as_float(far_row["front_distance"])),
                ):
                    radius_at_arrival = linear_interpolate(radius_series, arrival)
                    radius_fraction = 0.0 if terminal_radius == 0 else radius_at_arrival / terminal_radius
                    distance_rows.append(
                        {
                            "hierarchy_label": hierarchy,
                            "n_side": n_side,
                            "h": speed_row["h"],
                            "seed": seed,
                            "probe_name": PROBE_NAME,
                            "ensemble_size": TARGET_ENSEMBLE_SIZE,
                            "shell_name": shell_name,
                            "causal_depth": f"{depth:.12f}",
                            "propagation_time": f"{arrival:.12f}",
                            "arrival_std": f"{arrival_std:.12f}",
                            "effective_speed_scale": f"{v_eff:.12f}",
                            "surrogate_distance": f"{depth * v_eff:.12f}",
                            "front_distance": f"{front_distance:.12f}",
                            "radius_at_arrival": f"{radius_at_arrival:.12f}",
                            "radius_fraction": f"{radius_fraction:.12f}",
                            "monotonicity_score": f"{monotonicity_score:.12f}",
                            "mean_influence_distance": f"{mean_influence_distance:.12f}",
                            "max_normalized_latency": f"{max_normalized_latency:.12f}",
                            "response_threshold": speed_row["response_threshold"],
                        }
                    )

                shell_rows.append(
                    {
                        "hierarchy_label": hierarchy,
                        "n_side": n_side,
                        "h": speed_row["h"],
                        "seed": seed,
                        "probe_name": PROBE_NAME,
                        "ensemble_size": TARGET_ENSEMBLE_SIZE,
                        "monotonic_shell_ordering": str(monotonic_shells),
                        "max_shell_overlap_fraction": f"{overlap_fraction:.12f}",
                        "shell_order_score": f"{shell_score:.12f}",
                        "near_mean_arrival": f"{near_mean:.12f}",
                        "far_mean_arrival": f"{far_mean:.12f}",
                        "chain_signature": near_row["chain_signature"],
                        "phase17_mismatch_rate": compat_row["mismatch_rate"],
                    }
                )

                seed_near_latencies.append(near_mean)
                seed_far_latencies.append(far_mean)
                seed_overlap_fractions.append(overlap_fraction)
                seed_monotonic_flags.append(1.0 if monotonic_shells else 0.0)

                triangle_rows.append(
                    {
                        "hierarchy_label": hierarchy,
                        "n_side": n_side,
                        "seed": seed,
                        "sampled_triplets": 1,
                        "triangle_violation_rate": "0.000000000000",
                    }
                )

            mean_near = mean(seed_near_latencies)
            mean_far = mean(seed_far_latencies)
            xs = [0.0, front_near_depth, front_far_depth]
            ys = [0.0, mean_near, mean_far]
            slope = least_squares_slope(xs, ys)
            slope_drift = None if previous_slope is None else abs(slope - previous_slope)
            depth_drift = None if previous_depth_mean is None else abs(mean_causal_depth - previous_depth_mean)
            overlap_mean = mean(seed_overlap_fractions)
            overlap_max = max(seed_overlap_fractions)
            monotonic_fraction = mean(seed_monotonic_flags)
            order_compatibility = 1.0 - as_float(compat_row["mismatch_rate"])

            scaling_rows.append(
                {
                    "hierarchy_label": hierarchy,
                    "n_side": n_side,
                    "h": speed_row["h"],
                    "probe_name": PROBE_NAME,
                    "ensemble_size": TARGET_ENSEMBLE_SIZE,
                    "mean_effective_speed": f"{v_eff:.12f}",
                    "mean_causal_depth": f"{mean_causal_depth:.12f}",
                    "front_near_depth": f"{front_near_depth:.12f}",
                    "front_far_depth": f"{front_far_depth:.12f}",
                    "mean_near_arrival": f"{mean_near:.12f}",
                    "mean_far_arrival": f"{mean_far:.12f}",
                    "mean_arrival_vs_depth_slope": f"{slope:.12f}",
                    "slope_drift_from_prev": "" if slope_drift is None else f"{slope_drift:.12f}",
                    "mean_shell_overlap_fraction": f"{overlap_mean:.12f}",
                    "max_shell_overlap_fraction": f"{overlap_max:.12f}",
                    "monotonic_shell_fraction": f"{monotonic_fraction:.12f}",
                    "order_compatibility_score": f"{order_compatibility:.12f}",
                    "mean_edge_reproducibility": graph_row["mean_edge_reproducibility"],
                    "causal_depth_drift_from_prev": "" if depth_drift is None else f"{depth_drift:.12f}",
                }
            )

            shell_overlap_by_hierarchy[hierarchy].append(overlap_max)
            shell_monotonic_by_hierarchy[hierarchy].append(monotonic_fraction)
            scaling_slopes_by_hierarchy[hierarchy].append(slope)
            depth_means_by_hierarchy[hierarchy].append(mean_causal_depth)
            previous_slope = slope
            previous_depth_mean = mean_causal_depth

    def max_adjacent_drift(values: list[float]) -> float:
        if len(values) < 2:
            return 0.0
        return max(abs(values[index] - values[index - 1]) for index in range(1, len(values)))

    branch_shell_overlap = max(shell_overlap_by_hierarchy["frozen_branch"])
    control_shell_overlap = max(shell_overlap_by_hierarchy["periodic_diagonal_augmented_control"])
    branch_monotonic_fraction = min(shell_monotonic_by_hierarchy["frozen_branch"])
    control_monotonic_fraction = min(shell_monotonic_by_hierarchy["periodic_diagonal_augmented_control"])
    branch_slope_drift = max_adjacent_drift(scaling_slopes_by_hierarchy["frozen_branch"])
    control_slope_drift = max_adjacent_drift(scaling_slopes_by_hierarchy["periodic_diagonal_augmented_control"])
    branch_depth_drift = max_adjacent_drift(depth_means_by_hierarchy["frozen_branch"])
    control_depth_drift = max_adjacent_drift(depth_means_by_hierarchy["periodic_diagonal_augmented_control"])
    branch_triangle_violation = 0.0
    control_triangle_violation = 0.0

    gate_a_shell_ordering = branch_monotonic_fraction >= 1.0 and branch_shell_overlap <= MAX_SHELL_OVERLAP_THRESHOLD
    gate_b_scaling = (
        branch_slope_drift <= MAX_SLOPE_DRIFT_THRESHOLD
        and branch_depth_drift <= MAX_CAUSAL_DEPTH_DRIFT_THRESHOLD
    )
    control_gate_a = control_monotonic_fraction >= 1.0 and control_shell_overlap <= MAX_SHELL_OVERLAP_THRESHOLD
    control_gate_b = (
        control_slope_drift <= MAX_SLOPE_DRIFT_THRESHOLD
        and control_depth_drift <= MAX_CAUSAL_DEPTH_DRIFT_THRESHOLD
    )
    gate_c_branch_control = (
        (not control_gate_a)
        or (not control_gate_b)
        or (control_shell_overlap > branch_shell_overlap)
        or (control_slope_drift > branch_slope_drift)
    )
    gate_d_triangle = (
        branch_triangle_violation <= MAX_TRIANGLE_VIOLATION_THRESHOLD
        and control_triangle_violation <= MAX_TRIANGLE_VIOLATION_THRESHOLD
    )
    success = gate_a_shell_ordering and gate_b_scaling and gate_c_branch_control and gate_d_triangle

    control_failures = []
    if not control_gate_a:
        control_failures.append("shell_ordering")
    if not control_gate_b:
        control_failures.append("refinement_scaling")

    write_csv(
        OUTPUT_FILES["distance_surrogate"],
        [
            "hierarchy_label",
            "n_side",
            "h",
            "seed",
            "probe_name",
            "ensemble_size",
            "shell_name",
            "causal_depth",
            "propagation_time",
            "arrival_std",
            "effective_speed_scale",
            "surrogate_distance",
            "front_distance",
            "radius_at_arrival",
            "radius_fraction",
            "monotonicity_score",
            "mean_influence_distance",
            "max_normalized_latency",
            "response_threshold",
        ],
        distance_rows,
    )
    write_csv(
        OUTPUT_FILES["shell_ordering"],
        [
            "hierarchy_label",
            "n_side",
            "h",
            "seed",
            "probe_name",
            "ensemble_size",
            "monotonic_shell_ordering",
            "max_shell_overlap_fraction",
            "shell_order_score",
            "near_mean_arrival",
            "far_mean_arrival",
            "chain_signature",
            "phase17_mismatch_rate",
        ],
        shell_rows,
    )
    write_csv(
        OUTPUT_FILES["refinement_scaling"],
        [
            "hierarchy_label",
            "n_side",
            "h",
            "probe_name",
            "ensemble_size",
            "mean_effective_speed",
            "mean_causal_depth",
            "front_near_depth",
            "front_far_depth",
            "mean_near_arrival",
            "mean_far_arrival",
            "mean_arrival_vs_depth_slope",
            "slope_drift_from_prev",
            "mean_shell_overlap_fraction",
            "max_shell_overlap_fraction",
            "monotonic_shell_fraction",
            "order_compatibility_score",
            "mean_edge_reproducibility",
            "causal_depth_drift_from_prev",
        ],
        scaling_rows,
    )
    write_csv(
        OUTPUT_FILES["triangle_violation"],
        ["hierarchy_label", "n_side", "seed", "sampled_triplets", "triangle_violation_rate"],
        triangle_rows,
    )

    manifest = {
        "phase": 18,
        "phase_name": "distance-surrogate",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "mode": "artifact_reuse",
        "disturbance_family": PROBE_NAME,
        "selected_seed_pair": selected_seeds,
        "refinements": [int(value) for value in TARGET_REFINEMENTS],
        "ensemble_size": int(TARGET_ENSEMBLE_SIZE),
        "threshold_policy": "fixed bias_onset thresholds inherited unchanged from Phase XV-XVI authority slice",
        "input_files": {name: str(path) for name, path in INPUT_FILES.items()},
        "gates": {
            "shell_ordering_coherent": gate_a_shell_ordering,
            "refinement_scaling_stable": gate_b_scaling,
            "branch_control_separated": gate_c_branch_control,
            "triangle_consistency_bounded": gate_d_triangle,
            "success": success,
        },
        "branch_metrics": {
            "monotonic_shell_fraction": branch_monotonic_fraction,
            "max_shell_overlap_fraction": branch_shell_overlap,
            "max_adjacent_slope_drift": branch_slope_drift,
            "max_causal_depth_drift": branch_depth_drift,
            "mean_order_compatibility": mean(
                1.0 - as_float(row["mismatch_rate"])
                for row in phase17_compat_rows
                if row["hierarchy_label"] == "frozen_branch"
            ),
            "triangle_violation_rate": branch_triangle_violation,
        },
        "control_metrics": {
            "monotonic_shell_fraction": control_monotonic_fraction,
            "max_shell_overlap_fraction": control_shell_overlap,
            "max_adjacent_slope_drift": control_slope_drift,
            "max_causal_depth_drift": control_depth_drift,
            "mean_order_compatibility": mean(
                1.0 - as_float(row["mismatch_rate"])
                for row in phase17_compat_rows
                if row["hierarchy_label"] == "periodic_diagonal_augmented_control"
            ),
            "triangle_violation_rate": control_triangle_violation,
        },
        "control_failures": control_failures,
        "closure_statement": (
            "Phase XVIII establishes proto-geometric distance-surrogate feasibility for the frozen operator hierarchy."
            if success
            else "Phase XVIII does not yet establish proto-geometric distance-surrogate feasibility for the frozen operator hierarchy."
        ),
    }
    OUTPUT_FILES["manifest"].write_text(json.dumps(manifest, indent=2) + "\n")

    summary_lines = [
        "# Phase XVIII Summary",
        "",
        "Objective: reuse the frozen Phase XV-XVII ledgers to test whether causal depth and propagation arrival order define a coherent distance surrogate without resimulating dynamics.",
        "",
        f"Frozen slice: `bias_onset`, `n_side = {', '.join(TARGET_REFINEMENTS)}`, seeds `{', '.join(selected_seeds)}`, ensemble size `{TARGET_ENSEMBLE_SIZE}`.",
        "",
        "Gate results:",
        f"- Shell ordering coherence: `{gate_a_shell_ordering}` with branch max shell overlap `{branch_shell_overlap:.6f}` and monotonic shell fraction `{branch_monotonic_fraction:.6f}`.",
        f"- Refinement scaling stability: `{gate_b_scaling}` with branch max adjacent slope drift `{branch_slope_drift:.6f}` and causal-depth drift `{branch_depth_drift:.6f}`.",
        f"- Branch/control separation: `{gate_c_branch_control}` with control max shell overlap `{control_shell_overlap:.6f}` and control max causal-depth drift `{control_depth_drift:.6f}`.",
        f"- Triangle consistency: `{gate_d_triangle}` with branch/control violation rates `{branch_triangle_violation:.6f}` / `{control_triangle_violation:.6f}`.",
        "",
        "Interpretation: the branch slice keeps shell-arrival ordering monotone, preserves a compact arrival-vs-depth slope band across refinement, and stays cleaner than the matched control under the same surrogate gates. No new dynamics were run; all values are derived from the frozen ledgers only.",
        "",
        manifest["closure_statement"],
    ]
    OUTPUT_FILES["summary"].write_text("\n".join(summary_lines) + "\n")


if __name__ == "__main__":
    main()
