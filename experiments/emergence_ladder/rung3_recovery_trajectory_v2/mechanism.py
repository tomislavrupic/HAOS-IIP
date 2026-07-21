from __future__ import annotations

from dataclasses import dataclass, fields

import numpy as np


PROHIBITED_MEMORY_FIELD_FRAGMENTS = ("reference", "target", "future", "score", "classification", "checkpoint")


@dataclass(frozen=True)
class RelationalMemory:
    incidence: np.ndarray
    orientation: np.ndarray
    valid_edges: np.ndarray
    max_degree: float

    @classmethod
    def encode(cls, incidence: np.ndarray, state: np.ndarray, minimum_edge_magnitude: float = 0.08) -> "RelationalMemory":
        differences = incidence @ state
        valid = np.abs(differences) >= float(minimum_edge_magnitude)
        orientation = np.sign(differences)
        orientation[~valid] = 0.0
        max_degree = float(np.max(np.sum(np.abs(incidence), axis=0)))
        return cls(incidence.copy(), orientation.astype(float), valid.astype(bool), max(max_degree, 1.0))


@dataclass(frozen=True)
class StepDiagnostics:
    correction_norm: float
    active_fraction: float
    available_fraction: float
    compatibility_error: float


@dataclass(frozen=True)
class SimulationTrace:
    states: tuple[np.ndarray, ...]
    corrections: tuple[float, ...]
    active_fractions: tuple[float, ...]
    compatibility_errors: tuple[float, ...]


def assert_memory_has_no_reference_fields() -> None:
    names = {field.name.lower() for field in fields(RelationalMemory)}
    for fragment in PROHIBITED_MEMORY_FIELD_FRAGMENTS:
        if any(fragment in name for name in names):
            raise ValueError(f"prohibited memory field detected: {fragment}")


def compatibility(memory: RelationalMemory, state: np.ndarray, margin: float) -> tuple[np.ndarray, np.ndarray]:
    oriented = memory.orientation * (memory.incidence @ state)
    gaps = np.maximum(0.0, float(margin) - oriented)
    gaps[~memory.valid_edges] = 0.0
    return oriented, gaps


def correction_vector(
    state: np.ndarray,
    memory: RelationalMemory,
    eta: float,
    margin: float,
    availability: np.ndarray,
) -> tuple[np.ndarray, StepDiagnostics]:
    _, gaps = compatibility(memory, state, margin)
    active = memory.valid_edges & availability & (gaps > 0.0)
    edge_force = memory.orientation * gaps * active
    correction = float(eta) * (memory.incidence.T @ edge_force) / memory.max_degree
    correction -= np.mean(correction)
    valid_count = max(int(np.sum(memory.valid_edges)), 1)
    diagnostics = StepDiagnostics(
        correction_norm=float(np.sqrt(np.mean(correction**2))),
        active_fraction=float(np.sum(active) / valid_count),
        available_fraction=float(np.sum(availability & memory.valid_edges) / valid_count),
        compatibility_error=float(np.sqrt(np.sum(gaps**2) / valid_count)),
    )
    return correction, diagnostics


def simulate_feedback(
    state0: np.ndarray,
    memory: RelationalMemory,
    eta: float,
    margin: float,
    steps: int,
    dropout_fraction: float,
    seed: int,
) -> SimulationTrace:
    assert_memory_has_no_reference_fields()
    rng = np.random.default_rng(seed)
    state = np.asarray(state0, dtype=float).copy()
    states = [state.copy()]
    corrections: list[float] = []
    active_fractions: list[float] = []
    errors: list[float] = []
    for _ in range(int(steps)):
        availability = rng.random(len(memory.orientation)) >= float(dropout_fraction)
        correction, diagnostics = correction_vector(state, memory, eta, margin, availability)
        state = state + correction
        if not np.isfinite(state).all():
            raise ValueError("non-finite restorative state")
        states.append(state.copy())
        corrections.append(diagnostics.correction_norm)
        active_fractions.append(diagnostics.active_fraction)
        errors.append(diagnostics.compatibility_error)
    return SimulationTrace(tuple(states), tuple(corrections), tuple(active_fractions), tuple(errors))
