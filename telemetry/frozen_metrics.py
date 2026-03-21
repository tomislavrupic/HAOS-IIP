from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass

import numpy as np


_EPS = 1.0e-12
_RECOVERY_TOLERANCES = {
    "width_growth": 0.08,
    "concentration_loss": 0.12,
    "participation_growth": 0.20,
    "overlap_loss": 0.08,
}


def _probabilities(state: np.ndarray) -> np.ndarray:
    values = np.abs(np.asarray(state, dtype=complex)) ** 2
    total = max(float(np.sum(values)), _EPS)
    return np.asarray(values / total, dtype=float)


def _safe_relative_change(new_value: float, old_value: float) -> float:
    return abs(float(new_value) - float(old_value)) / max(abs(float(old_value)), _EPS)


def overlap(reference: np.ndarray, state: np.ndarray) -> float:
    reference = np.asarray(reference, dtype=complex)
    state = np.asarray(state, dtype=complex)
    denom = max(float(np.linalg.norm(reference) * np.linalg.norm(state)), _EPS)
    return float(abs(np.vdot(reference, state)) / denom)


def localization_width(state: np.ndarray, coords: np.ndarray) -> float:
    probabilities = _probabilities(state)
    coords = np.asarray(coords, dtype=float)
    center = coords[int(np.argmax(probabilities))]
    delta = coords - center
    return float(np.sqrt(np.sum(probabilities * np.sum(delta * delta, axis=1))))


def concentration_retention(state: np.ndarray, mask: np.ndarray) -> float:
    probabilities = _probabilities(state)
    mask = np.asarray(mask, dtype=bool)
    return float(np.sum(probabilities[mask]))


def participation_ratio(state: np.ndarray) -> float:
    probabilities = _probabilities(state)
    return float(1.0 / max(float(np.sum(probabilities * probabilities)), _EPS))


def recovery_score(
    reference: np.ndarray,
    state: np.ndarray,
    coords: np.ndarray,
    mask: np.ndarray,
) -> float:
    ref_width = localization_width(reference, coords)
    state_width = localization_width(state, coords)
    ref_participation = participation_ratio(reference)
    state_participation = participation_ratio(state)
    concentration = concentration_retention(state, mask)
    overlap_value = overlap(reference, state)

    penalties = [
        min(_safe_relative_change(state_width, ref_width) / _RECOVERY_TOLERANCES["width_growth"], 1.0),
        min((1.0 - concentration) / _RECOVERY_TOLERANCES["concentration_loss"], 1.0),
        min(_safe_relative_change(state_participation, ref_participation) / _RECOVERY_TOLERANCES["participation_growth"], 1.0),
        min((1.0 - overlap_value) / _RECOVERY_TOLERANCES["overlap_loss"], 1.0),
    ]
    return float(max(0.0, 1.0 - float(np.mean(penalties))))


@dataclass(frozen=True)
class SurvivalThresholds:
    max_width_growth: float
    min_concentration: float
    max_participation_growth: float
    min_overlap: float
    min_recovery_score: float


def classify_single_mode(
    reference,
    state,
    coords,
    mask,
    thresholds: SurvivalThresholds,
) -> str:
    ref_width = localization_width(reference, coords)
    state_width = localization_width(state, coords)
    width_growth = 0.0 if ref_width <= _EPS else max((state_width / ref_width) - 1.0, 0.0)

    ref_participation = participation_ratio(reference)
    state_participation = participation_ratio(state)
    participation_growth = 0.0 if ref_participation <= _EPS else max((state_participation / ref_participation) - 1.0, 0.0)

    concentration = concentration_retention(state, mask)
    overlap_value = overlap(reference, state)
    score = recovery_score(reference, state, coords, mask)

    if (
        width_growth <= float(thresholds.max_width_growth)
        and concentration >= float(thresholds.min_concentration)
        and participation_growth <= float(thresholds.max_participation_growth)
        and overlap_value >= float(thresholds.min_overlap)
        and score >= float(thresholds.min_recovery_score)
    ):
        return "persistent"

    if (
        concentration >= 0.5 * float(thresholds.min_concentration)
        and overlap_value >= 0.5 * float(thresholds.min_overlap)
        and score >= 0.5 * float(thresholds.min_recovery_score)
    ):
        return "diffusive"

    return "unstable"


def persistence_time(
    history: list[dict],
    tau_grid: list[float],
    thresholds: SurvivalThresholds,
) -> float:
    if len(history) != len(tau_grid):
        raise ValueError("history and tau_grid must have the same length")

    last_tau = float(tau_grid[0])
    for record, tau in zip(history, tau_grid):
        tau_value = float(tau)
        if {"reference", "state", "coords", "mask"}.issubset(record):
            label = classify_single_mode(
                record["reference"],
                record["state"],
                record["coords"],
                record["mask"],
                thresholds,
            )
        else:
            width_growth = float(record["width_growth"])
            concentration = float(record["concentration_retention"])
            participation_growth = float(record["participation_growth"])
            overlap_value = float(record["overlap"])
            score = float(record["recovery_score"])
            if (
                width_growth <= float(thresholds.max_width_growth)
                and concentration >= float(thresholds.min_concentration)
                and participation_growth <= float(thresholds.max_participation_growth)
                and overlap_value >= float(thresholds.min_overlap)
                and score >= float(thresholds.min_recovery_score)
            ):
                label = "persistent"
            elif (
                concentration >= 0.5 * float(thresholds.min_concentration)
                and overlap_value >= 0.5 * float(thresholds.min_overlap)
                and score >= 0.5 * float(thresholds.min_recovery_score)
            ):
                label = "diffusive"
            else:
                label = "unstable"

        if label != "persistent":
            break
        last_tau = tau_value
    return float(last_tau)


def first_threshold_crossing(signal: np.ndarray, threshold: float) -> int | None:
    values = np.asarray(signal, dtype=float)
    for index, value in enumerate(values):
        if float(value) >= float(threshold) - _EPS:
            return index
    return None


def front_arrival_order(probe_histories: dict[str, np.ndarray], threshold: float) -> dict[str, int | None]:
    return {
        probe_name: first_threshold_crossing(np.asarray(history, dtype=float), threshold)
        for probe_name, history in sorted(probe_histories.items())
    }


def reconstruct_influence_edges(event_times: dict[str, int | None]) -> list[tuple[str, str]]:
    ordered = [
        (label, value)
        for label, value in sorted(event_times.items(), key=lambda item: (float("inf") if item[1] is None else item[1], item[0]))
        if value is not None
    ]
    return [(ordered[index][0], ordered[index + 1][0]) for index in range(len(ordered) - 1)]


def _tarjan_scc(nodes: list[str], edges: list[tuple[str, str]]) -> tuple[list[list[str]], dict[str, int]]:
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
    for component_index, component in enumerate(components):
        for node in component:
            node_to_component[node] = component_index
    return components, node_to_component


def acyclicity_score(edges: list[tuple[str, str]], nodes: list[str]) -> float:
    if not edges:
        return 1.0
    components, node_to_component = _tarjan_scc(list(nodes), edges)
    cycle_edges = 0
    for src, dst in edges:
        src_component = node_to_component.get(src)
        dst_component = node_to_component.get(dst)
        if src_component is None or dst_component is None:
            continue
        if src_component == dst_component and len(components[src_component]) > 1:
            cycle_edges += 1
    return float(1.0 - (cycle_edges / max(len(edges), 1)))


def causal_depths(edges: list[tuple[str, str]], source: str) -> dict[str, int]:
    adjacency: dict[str, set[str]] = defaultdict(set)
    for src, dst in edges:
        adjacency[src].add(dst)

    depths = {source: 0}
    queue: deque[str] = deque([source])
    while queue:
        node = queue.popleft()
        for nxt in sorted(adjacency[node]):
            if nxt not in depths:
                depths[nxt] = depths[node] + 1
                queue.append(nxt)
    return depths


def order_compatibility(edge_set, arrival_ordering) -> float:
    comparable = 0
    mismatches = 0
    for src, dst in edge_set:
        src_time = arrival_ordering.get(src)
        dst_time = arrival_ordering.get(dst)
        if src_time is None or dst_time is None:
            continue
        comparable += 1
        if int(src_time) > int(dst_time):
            mismatches += 1
    if comparable == 0:
        return 0.0
    return float(mismatches / comparable)
