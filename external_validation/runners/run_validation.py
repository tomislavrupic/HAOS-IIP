#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol


EPSILON = 1.0e-12
RESULT_FIELDS = (
    "trial_id",
    "family",
    "schedule",
    "lambda_value",
    "lambda_step",
    "curvature_scale",
    "noise_pct",
    "reshuffle_pct",
    "sector_weights",
    "collapse_order",
    "alpha_k_star_step",
    "beta_k_star_step",
    "gamma_k_star_step",
    "global_first_k_star_step",
    "global_final_delta_persistence",
    "global_final_safety_margin",
    "global_final_confidence",
    "min_adjacent_gap",
)


@dataclass(frozen=True)
class Node:
    node_id: str
    sector: str


@dataclass(frozen=True)
class Edge:
    left: str
    right: str
    strength: float
    distance: float
    curvature_penalty: float


@dataclass(frozen=True)
class WeightedEdge:
    left: str
    right: str
    strength: float
    distance: float
    curvature_penalty: float
    score: float


@dataclass(frozen=True)
class GraphSpec:
    name: str
    nodes: tuple[Node, ...]
    edges: tuple[Edge, ...]


@dataclass(frozen=True)
class Snapshot:
    edges: tuple[WeightedEdge, ...]
    step_index: int
    focus_sector: str | None = None


@dataclass(frozen=True)
class SensorReading:
    delta_persistence: float
    k_star: int | None
    safety_margin: float
    confidence: float


@dataclass(frozen=True)
class TrialSpec:
    trial_id: str
    family: str
    schedule: str
    lambda_value: float
    lambda_step: float
    sector_weights: dict[str, float]
    curvature_scale: float
    noise_pct: float
    reshuffle_pct: float
    seed: int


class SensorFactory(Protocol):
    def __call__(self, nodes_by_id: dict[str, str], full_edges: tuple[WeightedEdge, ...]) -> "ToyHaosSensor":
        ...


class ToyHaosSensor:
    def __init__(
        self,
        nodes_by_id: dict[str, str],
        full_edges: tuple[WeightedEdge, ...],
        focus_sector: str | None = None,
        collapse_threshold: float = 0.33,
        activation_persistence: float = 0.62,
        activation_confidence: float = 0.44,
        drop_threshold: float = 0.18,
    ) -> None:
        self.nodes_by_id = nodes_by_id
        self.full_edges = full_edges
        self.focus_sector = focus_sector
        self.collapse_threshold = collapse_threshold
        self.activation_persistence = activation_persistence
        self.activation_confidence = activation_confidence
        self.drop_threshold = drop_threshold
        self.previous_persistence = 0.0
        self.peak_persistence = 0.0
        self.activated = False
        self.sectors = tuple(sorted(set(nodes_by_id.values())))
        self.distance_scale = max((edge.distance for edge in full_edges), default=1.0)
        self.reference = {
            sector: self._sector_observables(full_edges, sector)
            for sector in self.sectors
        }

    def observe(self, snapshot: Snapshot) -> SensorReading:
        if self.focus_sector is None:
            persistence_values: list[float] = []
            confidence_values: list[float] = []
            for sector in self.sectors:
                persistence, confidence = self._sector_state(snapshot.edges, sector)
                persistence_values.append(persistence)
                confidence_values.append(confidence)
            if persistence_values:
                persistence = sum(persistence_values) / len(persistence_values)
                spread = max(persistence_values) - min(persistence_values)
                persistence = clamp(persistence - 0.12 * spread)
                safety_margin = min(persistence_values) - self.collapse_threshold
                confidence = sum(confidence_values) / len(confidence_values)
                if confidence >= self.activation_confidence and any(
                    value >= self.activation_persistence for value in persistence_values
                ):
                    self.activated = True
                if self.activated:
                    collapsed = [
                        sector
                        for sector, value in zip(self.sectors, persistence_values, strict=True)
                        if value < self.collapse_threshold
                    ]
                    k_star = len(collapsed) if collapsed else None
                else:
                    k_star = None
            else:
                persistence = 0.0
                confidence = 0.0
                safety_margin = -self.collapse_threshold
                k_star = None
        else:
            persistence, confidence = self._sector_state(snapshot.edges, self.focus_sector)
            safety_margin = persistence - self.collapse_threshold
            delta_persistence = persistence - self.previous_persistence
            if confidence >= self.activation_confidence and persistence >= self.activation_persistence:
                self.activated = True
            if (
                self.activated
                and (self.peak_persistence - persistence) >= self.drop_threshold
                and delta_persistence < 0.0
            ):
                k_star = snapshot.step_index
            else:
                k_star = None
            self.previous_persistence = persistence
            self.peak_persistence = max(self.peak_persistence, persistence)
            return SensorReading(
                delta_persistence=delta_persistence,
                k_star=k_star,
                safety_margin=safety_margin,
                confidence=confidence,
            )

        delta_persistence = persistence - self.previous_persistence
        self.previous_persistence = persistence
        self.peak_persistence = max(self.peak_persistence, persistence)
        return SensorReading(
            delta_persistence=delta_persistence,
            k_star=k_star,
            safety_margin=safety_margin,
            confidence=confidence,
        )

    def _sector_observables(self, edges: tuple[WeightedEdge, ...], sector: str) -> dict[str, float]:
        within = 0.0
        cross = 0.0
        long_cross = 0.0
        incident = 0.0
        for edge in edges:
            left_sector = self.nodes_by_id[edge.left]
            right_sector = self.nodes_by_id[edge.right]
            if sector not in (left_sector, right_sector):
                continue
            incident += edge.score
            if left_sector == sector and right_sector == sector:
                within += edge.score
                continue
            cross += edge.score
            long_cross += edge.score * normalized_distance(edge.distance, self.distance_scale)
        return {
            "within": within,
            "cross": cross,
            "long_cross": long_cross,
            "incident": incident,
        }

    def _sector_state(self, edges: tuple[WeightedEdge, ...], sector: str) -> tuple[float, float]:
        current = self._sector_observables(edges, sector)
        reference = self.reference[sector]
        within_ratio = safe_ratio(current["within"], reference["within"])
        cross_ratio = safe_ratio(current["cross"], reference["cross"])
        long_ratio = safe_ratio(current["long_cross"], reference["long_cross"])
        incident_ratio = safe_ratio(current["incident"], reference["incident"])
        fragility = clamp(reference["cross"] / max(reference["within"], EPSILON), 0.0, 2.0)

        persistence = clamp(
            0.06
            + 1.18 * within_ratio
            - (0.50 + 0.28 * fragility) * cross_ratio
            - (0.16 + 0.12 * fragility) * long_ratio
            - 0.10 * max(0.0, cross_ratio - within_ratio)
        )
        confidence = clamp(0.35 * within_ratio + 0.65 * incident_ratio)
        return persistence, confidence


def clamp(value: float, lower: float = 0.0, upper: float = 1.0) -> float:
    return max(lower, min(upper, value))


def safe_ratio(value: float, reference: float) -> float:
    if reference <= EPSILON:
        return 0.0
    return clamp(value / reference)


def normalized_distance(distance: float, distance_scale: float) -> float:
    return clamp(distance / max(distance_scale, EPSILON))


def canonical_pair(left: str, right: str) -> tuple[str, str]:
    return (left, right) if left < right else (right, left)


def load_graph(path: Path) -> GraphSpec:
    payload = json.loads(path.read_text(encoding="utf-8"))
    nodes = tuple(
        Node(
            node_id=str(node["id"]),
            sector=str(node["sector"]),
        )
        for node in payload["nodes"]
    )
    node_ids = {node.node_id for node in nodes}
    edges = []
    for edge in payload["edges"]:
        left = str(edge["source"])
        right = str(edge["target"])
        if left not in node_ids or right not in node_ids:
            raise ValueError(f"edge references unknown node: {left}-{right}")
        if left == right:
            raise ValueError(f"self-loop not allowed in input graph: {left}")
        edges.append(
            Edge(
                left=canonical_pair(left, right)[0],
                right=canonical_pair(left, right)[1],
                strength=float(edge["strength"]),
                distance=float(edge["distance"]),
                curvature_penalty=float(edge["curvature_penalty"]),
            )
        )

    return GraphSpec(
        name=str(payload.get("name", path.stem)),
        nodes=nodes,
        edges=tuple(edges),
    )


def build_nodes_by_id(graph: GraphSpec) -> dict[str, str]:
    return {node.node_id: node.sector for node in graph.nodes}


def perturb_edges(graph: GraphSpec, noise_pct: float, reshuffle_pct: float, seed: int) -> tuple[Edge, ...]:
    rng = random.Random(seed)
    edges = [
        Edge(
            left=edge.left,
            right=edge.right,
            strength=max(0.01, edge.strength * (1.0 + rng.uniform(-noise_pct, noise_pct))),
            distance=max(0.05, edge.distance),
            curvature_penalty=clamp(edge.curvature_penalty),
        )
        for edge in graph.edges
    ]

    if reshuffle_pct <= 0.0:
        return tuple(edges)

    count = max(2, round(len(edges) * reshuffle_pct))
    count += count % 2
    indexes = sorted(rng.sample(range(len(edges)), min(count, len(edges))))
    rotated_right_nodes = [edges[index].right for index in indexes[1:]] + [edges[indexes[0]].right]
    mutated = list(edges)
    seen_pairs = {canonical_pair(edge.left, edge.right) for edge in edges}

    for index, new_right in zip(indexes, rotated_right_nodes, strict=True):
        edge = mutated[index]
        if edge.left == new_right:
            continue
        new_pair = canonical_pair(edge.left, new_right)
        old_pair = canonical_pair(edge.left, edge.right)
        if new_pair in seen_pairs and new_pair != old_pair:
            continue
        seen_pairs.discard(old_pair)
        seen_pairs.add(new_pair)
        mutated[index] = Edge(
            left=new_pair[0],
            right=new_pair[1],
            strength=edge.strength,
            distance=edge.distance,
            curvature_penalty=edge.curvature_penalty,
        )
    return tuple(mutated)


def score_edge(
    edge: Edge,
    nodes_by_id: dict[str, str],
    lambda_value: float,
    sector_weights: dict[str, float],
    curvature_scale: float,
) -> float:
    left_sector = nodes_by_id[edge.left]
    right_sector = nodes_by_id[edge.right]
    coupling = math.sqrt(sector_weights[left_sector] * sector_weights[right_sector])
    curvature = edge.curvature_penalty ** curvature_scale
    locality = math.exp(-edge.distance / lambda_value)
    return edge.strength * locality * curvature * coupling


def build_weighted_edges(graph: GraphSpec, trial: TrialSpec) -> tuple[WeightedEdge, ...]:
    nodes_by_id = build_nodes_by_id(graph)
    perturbed_edges = perturb_edges(
        graph,
        noise_pct=trial.noise_pct,
        reshuffle_pct=trial.reshuffle_pct,
        seed=trial.seed,
    )
    weighted = [
        WeightedEdge(
            left=edge.left,
            right=edge.right,
            strength=edge.strength,
            distance=edge.distance,
            curvature_penalty=edge.curvature_penalty,
            score=score_edge(
                edge=edge,
                nodes_by_id=nodes_by_id,
                lambda_value=trial.lambda_value,
                sector_weights=trial.sector_weights,
                curvature_scale=trial.curvature_scale,
            ),
        )
        for edge in perturbed_edges
    ]
    weighted.sort(
        key=lambda edge: (
            -edge.score,
            edge.distance,
            edge.left,
            edge.right,
        )
    )
    return tuple(weighted)


def run_filtration(
    graph: GraphSpec,
    full_edges: tuple[WeightedEdge, ...],
    sensor_factory: SensorFactory,
) -> dict[str, object]:
    nodes_by_id = build_nodes_by_id(graph)
    sectors = tuple(sorted(set(nodes_by_id.values())))
    global_sensor = sensor_factory(nodes_by_id, full_edges)
    sector_sensors = {
        sector: sensor_factory(nodes_by_id, full_edges, focus_sector=sector)
        for sector in sectors
    }

    series: list[dict[str, object]] = []
    sector_crossings: dict[str, int | None] = {sector: None for sector in sectors}
    global_first_k_star_step: int | None = None

    included: list[WeightedEdge] = []
    for step_index, edge in enumerate(full_edges, start=1):
        included.append(edge)
        snapshot = Snapshot(edges=tuple(included), step_index=step_index)
        global_reading = global_sensor.observe(snapshot)
        if global_first_k_star_step is None and global_reading.k_star is not None:
            global_first_k_star_step = step_index

        sector_readings: dict[str, SensorReading] = {}
        for sector, sensor in sector_sensors.items():
            reading = sensor.observe(
                Snapshot(edges=tuple(included), step_index=step_index, focus_sector=sector)
            )
            sector_readings[sector] = reading
            if sector_crossings[sector] is None and reading.k_star is not None:
                sector_crossings[sector] = step_index

        series.append(
            {
                "step_index": step_index,
                "included_edge": f"{edge.left}-{edge.right}",
                "score": edge.score,
                "global": global_reading,
                "sectors": sector_readings,
            }
        )

    ordered = sorted(
        sectors,
        key=lambda sector: (
            float("inf") if sector_crossings[sector] is None else sector_crossings[sector],
            sector,
        ),
    )
    ordered_steps = [sector_crossings[sector] for sector in ordered]
    min_adjacent_gap = min_adjacent_step_gap(ordered_steps)
    final_global = series[-1]["global"] if series else SensorReading(0.0, None, 0.0, 0.0)

    return {
        "series": series,
        "sector_crossings": sector_crossings,
        "collapse_order": ordered,
        "global_first_k_star_step": global_first_k_star_step,
        "global_final": final_global,
        "min_adjacent_gap": min_adjacent_gap,
    }


def min_adjacent_step_gap(steps: list[int | None]) -> int:
    concrete = [step for step in steps if step is not None]
    if len(concrete) < 2:
        return 0
    return min(max(next_step - step, 0) for step, next_step in zip(concrete, concrete[1:]))


def collapse_order_string(sector_crossings: dict[str, int | None]) -> str:
    grouped: dict[int | None, list[str]] = {}
    for sector, step in sector_crossings.items():
        grouped.setdefault(step, []).append(sector)

    ordered_parts: list[str] = []
    for step in sorted(grouped, key=lambda value: (float("inf") if value is None else value, value)):
        sectors = sorted(grouped[step])
        if step is None:
            ordered_parts.append(f"uncollapsed({', '.join(sectors)})")
        else:
            ordered_parts.append(" = ".join(sectors))
    return " > ".join(ordered_parts)


def enumerate_trials() -> list[TrialSpec]:
    trials: list[TrialSpec] = []
    seed = 1303
    base_weights = {"alpha": 1.00, "beta": 1.00, "gamma": 1.00}
    locality_schedules = {
        "coarse": (0.34, 0.46, 0.58),
        "baseline": (0.34, 0.42, 0.50, 0.58),
        "fine": (0.34, 0.38, 0.42, 0.46, 0.50, 0.54, 0.58),
    }
    locality_perturbations = {
        "coarse": ((0.00, 0.00),),
        "fine": ((0.00, 0.00),),
        "baseline": ((0.00, 0.00), (0.05, 0.00), (0.10, 0.10)),
    }
    for schedule, lambda_values in locality_schedules.items():
        lambda_step = round(lambda_values[1] - lambda_values[0], 2) if len(lambda_values) > 1 else 0.0
        for lambda_value in lambda_values:
            for noise_pct, reshuffle_pct in locality_perturbations[schedule]:
                trials.append(
                    TrialSpec(
                        trial_id=f"locality_{schedule}_l{int(round(lambda_value * 100)):02d}_n{int(noise_pct * 100):02d}_r{int(reshuffle_pct * 100):02d}",
                        family="locality",
                        schedule=schedule,
                        lambda_value=lambda_value,
                        lambda_step=lambda_step,
                        sector_weights=dict(base_weights),
                        curvature_scale=1.0,
                        noise_pct=noise_pct,
                        reshuffle_pct=reshuffle_pct,
                        seed=seed,
                    )
                )
                seed += 1

    coupling_variants = {
        "baseline": {"alpha": 1.00, "beta": 1.00, "gamma": 1.00},
        "alpha_bias": {"alpha": 1.08, "beta": 0.96, "gamma": 1.00},
        "beta_bias": {"alpha": 0.96, "beta": 1.08, "gamma": 1.00},
        "gamma_bias": {"alpha": 1.00, "beta": 0.96, "gamma": 1.08},
    }
    for schedule, sector_weights in coupling_variants.items():
        for noise_pct, reshuffle_pct in ((0.00, 0.00), (0.05, 0.10)):
            trials.append(
                TrialSpec(
                    trial_id=f"couplings_{schedule}_n{int(noise_pct * 100):02d}_r{int(reshuffle_pct * 100):02d}",
                    family="couplings",
                    schedule=schedule,
                    lambda_value=0.42,
                    lambda_step=0.08,
                    sector_weights=sector_weights,
                    curvature_scale=1.0,
                    noise_pct=noise_pct,
                    reshuffle_pct=reshuffle_pct,
                    seed=seed,
                )
            )
            seed += 1

    for schedule, curvature_scale in (("light", 0.85), ("baseline", 1.0), ("strong", 1.15)):
        for noise_pct, reshuffle_pct in ((0.00, 0.00), (0.05, 0.10)):
            trials.append(
                TrialSpec(
                    trial_id=f"curvature_{schedule}_n{int(noise_pct * 100):02d}_r{int(reshuffle_pct * 100):02d}",
                    family="curvature",
                    schedule=schedule,
                    lambda_value=0.42,
                    lambda_step=0.08,
                    sector_weights=dict(base_weights),
                    curvature_scale=curvature_scale,
                    noise_pct=noise_pct,
                    reshuffle_pct=reshuffle_pct,
                    seed=seed,
                )
            )
            seed += 1
    return trials


def compare_pairwise(
    reference_steps: dict[str, int | None],
    candidate_steps: dict[str, int | None],
) -> tuple[int, int]:
    flips = 0
    total = 0
    sectors = sorted(reference_steps)
    for index, left in enumerate(sectors):
        for right in sectors[index + 1 :]:
            ref_left = float("inf") if reference_steps[left] is None else reference_steps[left]
            ref_right = float("inf") if reference_steps[right] is None else reference_steps[right]
            cand_left = float("inf") if candidate_steps[left] is None else candidate_steps[left]
            cand_right = float("inf") if candidate_steps[right] is None else candidate_steps[right]
            ref_sign = 0 if ref_left == ref_right else (-1 if ref_left < ref_right else 1)
            cand_sign = 0 if cand_left == cand_right else (-1 if cand_left < cand_right else 1)
            if cand_sign != ref_sign:
                flips += 1
            total += 1
    return flips, total


def gap_behavior(baseline_gap: int, observed_gaps: list[int]) -> str:
    if not observed_gaps:
        return "vanishes"
    zero_share = sum(1 for gap in observed_gaps if gap <= 0) / len(observed_gaps)
    median_gap = sorted(observed_gaps)[len(observed_gaps) // 2]
    if zero_share >= 0.25 or median_gap <= 0:
        return "vanishes"
    if median_gap >= max(baseline_gap, 1) and zero_share <= 0.10:
        return "persists"
    return "shrinks"


def invariance_assessment(
    dominant_share: float,
    flip_rate: float,
    gap_state: str,
) -> str:
    if dominant_share >= 0.80 and flip_rate <= 0.20 and gap_state != "vanishes":
        return "invariant structural signal"
    return "artifact of construction"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def write_summary(
    path: Path,
    title: str,
    baseline_row: dict[str, object],
    dominant_order: str,
    dominant_count: int,
    total_trials: int,
    flip_rate: float,
    gap_state: str,
    assessment: str,
) -> None:
    lines = [
        f"# {title}",
        "",
        f"- Collapse ordering by sector: `{baseline_row['collapse_order']}`",
        f"- Baseline k* steps: `alpha={baseline_row['alpha_k_star_step']}`, `beta={baseline_row['beta_k_star_step']}`, `gamma={baseline_row['gamma_k_star_step']}`",
        f"- Consistency across trials (flip rate): `{flip_rate:.3f}`",
        f"- Dominant ordering count: `{dominant_order}` in `{dominant_count}/{total_trials}` trials",
        f"- Gap behavior: `{gap_state}`",
        f"- Invariance assessment: `{assessment}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    validation_root = root / "external_validation"
    parser = argparse.ArgumentParser(description="Run the external HAOS validation runner on a graph artifact.")
    parser.add_argument(
        "--sensor",
        choices=("toy", "haos_iip"),
        default="toy",
        help="Sensor backend to use.",
    )
    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Optional output prefix. Defaults to 'toy' paths for the toy sensor and 'haos_iip' for the HAOS-IIP adapter.",
    )
    parser.add_argument(
        "--haos-iip-root",
        default=str(root),
        help="Path to the frozen HAOS-IIP root used by the HAOS-IIP sensor adapter.",
    )
    parser.add_argument(
        "--graph-file",
        default=str(validation_root / "data" / "toy_graph.json"),
        help="Path to the graph JSON artifact consumed by the runner.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = Path(__file__).resolve().parents[2]
    validation_root = root / "external_validation"
    results_dir = validation_root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    if args.sensor == "toy" and args.output_prefix is None:
        csv_path = results_dir / "toy_sweeps.csv"
        summary_path = results_dir / "toy_summaries.md"
        summary_title = "Toy Validation Summary"
    else:
        prefix = args.output_prefix or args.sensor
        csv_path = results_dir / f"{prefix}_sweeps.csv"
        summary_path = results_dir / f"{prefix}_summaries.md"
        summary_title = "HAOS-IIP Validation Summary" if args.sensor == "haos_iip" else "Toy Validation Summary"

    graph = load_graph(Path(args.graph_file))
    rows: list[dict[str, object]] = []
    order_counts: dict[str, int] = {}
    pairwise_flips = 0
    pairwise_total = 0
    gaps: list[int] = []
    baseline_reference: dict[str, object] | None = None
    baseline_steps: dict[str, int | None] | None = None

    if args.sensor == "haos_iip":
        try:
            from haos_iip_sensor_adapter import build_haos_iip_factory
        except ModuleNotFoundError as exc:
            raise RuntimeError(
                "The haos_iip sensor backend requires a Python runtime with numpy available."
            ) from exc
        factory = build_haos_iip_factory(Path(args.haos_iip_root))
    else:
        def factory(
            nodes_by_id: dict[str, str],
            full_edges: tuple[WeightedEdge, ...],
            focus_sector: str | None = None,
        ) -> ToyHaosSensor:
            return ToyHaosSensor(nodes_by_id=nodes_by_id, full_edges=full_edges, focus_sector=focus_sector)

    for trial in enumerate_trials():
        weighted_edges = build_weighted_edges(graph, trial)
        result = run_filtration(graph, weighted_edges, factory)
        sector_crossings = result["sector_crossings"]
        collapse_order = collapse_order_string(sector_crossings)
        global_final = result["global_final"]
        row = {
            "trial_id": trial.trial_id,
            "family": trial.family,
            "schedule": trial.schedule,
            "lambda_value": f"{trial.lambda_value:.2f}",
            "lambda_step": f"{trial.lambda_step:.2f}",
            "curvature_scale": f"{trial.curvature_scale:.2f}",
            "noise_pct": f"{trial.noise_pct:.2f}",
            "reshuffle_pct": f"{trial.reshuffle_pct:.2f}",
            "sector_weights": json.dumps(trial.sector_weights, sort_keys=True),
            "collapse_order": collapse_order,
            "alpha_k_star_step": sector_crossings["alpha"],
            "beta_k_star_step": sector_crossings["beta"],
            "gamma_k_star_step": sector_crossings["gamma"],
            "global_first_k_star_step": result["global_first_k_star_step"],
            "global_final_delta_persistence": f"{global_final.delta_persistence:.6f}",
            "global_final_safety_margin": f"{global_final.safety_margin:.6f}",
            "global_final_confidence": f"{global_final.confidence:.6f}",
            "min_adjacent_gap": result["min_adjacent_gap"],
        }
        rows.append(row)
        order_counts[collapse_order] = order_counts.get(collapse_order, 0) + 1
        gaps.append(int(result["min_adjacent_gap"]))

        if trial.trial_id == "locality_baseline_l42_n00_r00":
            baseline_reference = row
            baseline_steps = dict(sector_crossings)
        elif baseline_steps is not None:
            flips, total = compare_pairwise(baseline_steps, sector_crossings)
            pairwise_flips += flips
            pairwise_total += total

    if baseline_reference is None or baseline_steps is None:
        raise RuntimeError("baseline reference trial was not generated")

    dominant_order, dominant_count = max(order_counts.items(), key=lambda item: (item[1], item[0]))
    baseline_gap = int(baseline_reference["min_adjacent_gap"])
    flip_rate = float(pairwise_flips / max(pairwise_total, 1))
    gap_state = gap_behavior(baseline_gap, gaps)
    dominant_share = dominant_count / max(len(rows), 1)
    assessment = invariance_assessment(dominant_share, flip_rate, gap_state)

    write_csv(csv_path, rows)
    write_summary(
        summary_path,
        title=summary_title,
        baseline_row=baseline_reference,
        dominant_order=dominant_order,
        dominant_count=dominant_count,
        total_trials=len(rows),
        flip_rate=flip_rate,
        gap_state=gap_state,
        assessment=assessment,
    )


if __name__ == "__main__":
    main()
