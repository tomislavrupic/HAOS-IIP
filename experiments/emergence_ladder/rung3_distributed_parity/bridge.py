from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .alphabet import RelationalAlphabet


@dataclass(frozen=True)
class BridgeTrace:
    states: tuple[np.ndarray, ...]
    correction_norms: tuple[float, ...]
    active_relations: tuple[int, ...]


def correction_vector(
    state: np.ndarray,
    target_symbols: np.ndarray,
    alphabet: RelationalAlphabet,
    eta: float,
    margin: float,
    enabled_relations: np.ndarray | None = None,
) -> np.ndarray:
    desired = 2.0 * np.asarray(target_symbols, dtype=float) - 1.0
    differences = alphabet.incidence @ np.asarray(state, dtype=float)
    gaps = np.maximum(0.0, float(margin) - desired * differences)
    enabled = np.ones(len(desired), dtype=bool) if enabled_relations is None else np.asarray(enabled_relations, dtype=bool)
    forces = desired * gaps * enabled
    degree = max(float(np.max(np.sum(np.abs(alphabet.incidence), axis=0))), 1.0)
    correction = float(eta) * (alphabet.incidence.T @ forces) / degree
    correction -= np.mean(correction)
    return correction


def simulate_bridge(
    state0: np.ndarray,
    target_symbols: np.ndarray,
    alphabet: RelationalAlphabet,
    eta: float,
    margin: float,
    steps: int,
    enabled_relations: np.ndarray | None = None,
) -> BridgeTrace:
    state = np.asarray(state0, dtype=float).copy()
    states = [state.copy()]
    costs: list[float] = []
    active: list[int] = []
    for _ in range(int(steps)):
        before = alphabet.encode(state)
        correction = correction_vector(state, target_symbols, alphabet, eta, margin, enabled_relations)
        state = state + correction
        if not np.isfinite(state).all():
            raise ValueError("non-finite parity bridge state")
        states.append(state.copy())
        costs.append(float(np.sqrt(np.mean(correction**2))))
        active.append(int(np.sum(before != target_symbols)))
    return BridgeTrace(tuple(states), tuple(costs), tuple(active))

