from __future__ import annotations

import importlib
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np


EPSILON = 1.0e-12


@dataclass(frozen=True)
class HaosSensorReading:
    delta_persistence: float
    k_star: int | None
    safety_margin: float
    confidence: float


@dataclass(frozen=True)
class TelemetryBundle:
    thresholds: Any
    acyclicity_score: Any
    causal_depths: Any
    classify_single_mode: Any
    front_arrival_order: Any
    order_compatibility: Any
    overlap: Any
    persistence_time: Any
    reconstruct_influence_edges: Any
    recovery_score: Any


def clamp(value: float, lower: float = 0.0, upper: float = 1.0) -> float:
    return max(lower, min(upper, value))


@lru_cache(maxsize=4)
def load_telemetry_bundle(root_str: str) -> TelemetryBundle:
    root = Path(root_str).resolve()
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))

    telemetry = importlib.import_module("telemetry")
    demo = importlib.import_module("haos_iip.demo")
    return TelemetryBundle(
        thresholds=demo.THRESHOLDS,
        acyclicity_score=telemetry.acyclicity_score,
        causal_depths=telemetry.causal_depths,
        classify_single_mode=telemetry.classify_single_mode,
        front_arrival_order=telemetry.front_arrival_order,
        order_compatibility=telemetry.order_compatibility,
        overlap=telemetry.overlap,
        persistence_time=telemetry.persistence_time,
        reconstruct_influence_edges=telemetry.reconstruct_influence_edges,
        recovery_score=telemetry.recovery_score,
    )


class HaosIipSensor:
    def __init__(
        self,
        nodes_by_id: dict[str, str],
        full_edges: tuple[Any, ...],
        focus_sector: str | None = None,
        haos_iip_root: Path | None = None,
        ordering_threshold: float = 0.70,
    ) -> None:
        root = Path(haos_iip_root) if haos_iip_root is not None else Path(__file__).resolve().parents[2]
        self.bundle = load_telemetry_bundle(str(root))
        self.nodes_by_id = nodes_by_id
        self.full_edges = full_edges
        self.focus_sector = focus_sector
        self.ordering_threshold = ordering_threshold
        self.node_order = tuple(sorted(nodes_by_id))
        self.node_index = {node_id: index for index, node_id in enumerate(self.node_order)}
        self.coords = self._build_coords()
        self.reference_weights = self._build_reference_weights()
        self.reference_state = self._weights_to_state(self.reference_weights)
        self.mask = self._build_mask()
        self.source_node = self.node_order[int(np.argmax(self.reference_weights))]
        self.arrival_threshold = self._build_arrival_threshold()
        self.history_records: list[dict[str, np.ndarray]] = []
        self.probability_history: list[np.ndarray] = []
        self.previous_persistence = 0.0
        self.previous_label = "diffusive"
        self.activated = False

    def observe(self, snapshot: Any) -> HaosSensorReading:
        current_weights = self._build_current_weights(snapshot.edges)
        current_state = self._weights_to_state(current_weights)
        record = {
            "reference": self.reference_state,
            "state": current_state,
            "coords": self.coords,
            "mask": self.mask,
        }
        self.history_records.append(record)
        current_probabilities = self._probabilities(current_state)
        self.probability_history.append(current_probabilities)

        tau_grid = list(range(1, len(self.history_records) + 1))
        persistence_tau = float(
            self.bundle.persistence_time(self.history_records, tau_grid, self.bundle.thresholds)
        )
        persistence_fraction = persistence_tau / max(float(tau_grid[-1]), 1.0)
        delta_persistence = persistence_fraction - self.previous_persistence

        label = self.bundle.classify_single_mode(
            self.reference_state,
            current_state,
            self.coords,
            self.mask,
            self.bundle.thresholds,
        )
        overlap_value = float(self.bundle.overlap(self.reference_state, current_state))
        recovery = float(
            self.bundle.recovery_score(
                self.reference_state,
                current_state,
                self.coords,
                self.mask,
            )
        )
        ordering_score = self._ordering_score()

        safety_margin = min(
            recovery - float(self.bundle.thresholds.min_recovery_score),
            overlap_value - float(self.bundle.thresholds.min_overlap),
            ordering_score - self.ordering_threshold,
        )
        confidence = clamp(0.45 * recovery + 0.35 * overlap_value + 0.20 * ordering_score)

        if label == "persistent":
            self.activated = True

        if self.activated and self.previous_label == "persistent" and label != "persistent":
            k_star = int(snapshot.step_index)
        else:
            k_star = None

        self.previous_persistence = persistence_fraction
        self.previous_label = label
        return HaosSensorReading(
            delta_persistence=delta_persistence,
            k_star=k_star,
            safety_margin=safety_margin,
            confidence=confidence,
        )

    def _build_coords(self) -> np.ndarray:
        sectors = tuple(sorted(set(self.nodes_by_id.values())))
        sector_index = {sector: index for index, sector in enumerate(sectors)}
        counts: dict[str, int] = {sector: 0 for sector in sectors}
        coords = []
        for node_id in self.node_order:
            sector = self.nodes_by_id[node_id]
            index_in_sector = counts[sector]
            counts[sector] += 1
            x = float(2 * sector_index[sector])
            y = float(index_in_sector - 1)
            coords.append((x, y))
        return np.asarray(coords, dtype=float)

    def _build_mask(self) -> np.ndarray:
        if self.focus_sector is None:
            return np.asarray(self.reference_weights > EPSILON, dtype=bool)
        return np.asarray(
            [self.nodes_by_id[node_id] == self.focus_sector for node_id in self.node_order],
            dtype=bool,
        )

    def _build_arrival_threshold(self) -> float:
        reference_probabilities = self._probabilities(self.reference_state)
        masked = reference_probabilities[self.mask]
        if masked.size == 0:
            return 0.0
        return float(0.45 * np.max(masked))

    def _build_reference_weights(self) -> np.ndarray:
        weights = np.zeros(len(self.node_order), dtype=float)
        for edge in self.full_edges:
            left_sector = self.nodes_by_id[edge.left]
            right_sector = self.nodes_by_id[edge.right]
            if self.focus_sector is None:
                if left_sector != right_sector:
                    continue
            else:
                if left_sector != self.focus_sector or right_sector != self.focus_sector:
                    continue
            weights[self.node_index[edge.left]] += float(edge.score)
            weights[self.node_index[edge.right]] += float(edge.score)

        if float(np.sum(weights)) <= EPSILON:
            for index, node_id in enumerate(self.node_order):
                if self.focus_sector is None or self.nodes_by_id[node_id] == self.focus_sector:
                    weights[index] = 1.0
        return weights

    def _build_current_weights(self, edges: tuple[Any, ...]) -> np.ndarray:
        weights = np.zeros(len(self.node_order), dtype=float)
        for edge in edges:
            left_sector = self.nodes_by_id[edge.left]
            right_sector = self.nodes_by_id[edge.right]
            if self.focus_sector is None:
                include = True
            else:
                include = self.focus_sector in (left_sector, right_sector)
            if not include:
                continue
            weights[self.node_index[edge.left]] += float(edge.score)
            weights[self.node_index[edge.right]] += float(edge.score)
        if float(np.sum(weights)) <= EPSILON:
            weights += np.asarray(self.mask, dtype=float)
        return weights

    def _ordering_score(self) -> float:
        if len(self.probability_history) < 2:
            return 1.0

        probe_histories = {
            node_id: np.asarray(
                [history[self.node_index[node_id]] for history in self.probability_history],
                dtype=float,
            )
            for node_id in self.node_order
        }
        arrival_order = self.bundle.front_arrival_order(probe_histories, self.arrival_threshold)
        influence_edges = self.bundle.reconstruct_influence_edges(arrival_order)
        acyclic = float(self.bundle.acyclicity_score(influence_edges, list(self.node_order)))
        mismatch = float(self.bundle.order_compatibility(influence_edges, arrival_order))
        depths = self.bundle.causal_depths(influence_edges, self.source_node)
        coverage = len(depths) / max(len(self.node_order), 1)
        return clamp(acyclic * (1.0 - mismatch) * coverage)

    def _weights_to_state(self, weights: np.ndarray) -> np.ndarray:
        return np.sqrt(np.asarray(weights, dtype=float) + EPSILON)

    def _probabilities(self, state: np.ndarray) -> np.ndarray:
        values = np.abs(np.asarray(state, dtype=float)) ** 2
        total = max(float(np.sum(values)), EPSILON)
        return np.asarray(values / total, dtype=float)


def build_haos_iip_factory(haos_iip_root: Path):
    def factory(
        nodes_by_id: dict[str, str],
        full_edges: tuple[Any, ...],
        focus_sector: str | None = None,
    ) -> HaosIipSensor:
        return HaosIipSensor(
            nodes_by_id=nodes_by_id,
            full_edges=full_edges,
            focus_sector=focus_sector,
            haos_iip_root=haos_iip_root,
        )

    return factory
