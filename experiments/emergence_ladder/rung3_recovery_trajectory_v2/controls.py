from __future__ import annotations

import numpy as np

from .fixtures import GridFixture, rewired_incidence, spanning_tree_incidence
from .mechanism import RelationalMemory, SimulationTrace, correction_vector, simulate_feedback


ADMISSIBLE_CONTROLS = (
    "passive_relaxation",
    "operator_only_filtering",
    "topology_altered",
    "trivial_attractor",
    "identity_scrambled",
    "signal_blind",
    "parameter_budget_matched",
    "memory_ablation",
    "redundancy_ablation",
)
ORACLE_CONTROL = "oracle_upper_bound_non_admissible"


def static_trace(state: np.ndarray, steps: int) -> SimulationTrace:
    states = tuple(np.asarray(state, dtype=float).copy() for _ in range(steps + 1))
    return SimulationTrace(states, tuple(0.0 for _ in range(steps)), tuple(0.0 for _ in range(steps)), tuple(0.0 for _ in range(steps)))


def operator_filter_trace(state0: np.ndarray, incidence: np.ndarray, eta: float, margin: float, steps: int) -> SimulationTrace:
    state = state0.copy()
    degree = max(float(np.max(np.sum(np.abs(incidence), axis=0))), 1.0)
    states = [state.copy()]
    corrections = []
    for _ in range(steps):
        differences = incidence @ state
        active = np.abs(differences) > margin
        correction = -eta * (incidence.T @ (differences * active)) / degree
        correction -= np.mean(correction)
        state = state + correction
        states.append(state.copy())
        corrections.append(float(np.sqrt(np.mean(correction**2))))
    zeros = tuple(0.0 for _ in range(steps))
    return SimulationTrace(tuple(states), tuple(corrections), zeros, zeros)


def signal_blind_trace(
    state0: np.ndarray,
    memory: RelationalMemory,
    eta: float,
    margin: float,
    steps: int,
    dropout_fraction: float,
    seed: int,
) -> SimulationTrace:
    rng = np.random.default_rng(seed)
    state = state0.copy()
    states = [state.copy()]
    costs = []
    active = []
    errors = []
    for _ in range(steps):
        availability = rng.random(len(memory.orientation)) >= dropout_fraction
        aligned, diagnostics = correction_vector(state, memory, eta, margin, availability)
        direction = rng.normal(size=len(state))
        direction -= np.mean(direction)
        direction /= max(float(np.sqrt(np.mean(direction**2))), 1.0e-12)
        correction = direction * diagnostics.correction_norm
        state = state + correction
        states.append(state.copy())
        costs.append(diagnostics.correction_norm)
        active.append(diagnostics.active_fraction)
        errors.append(diagnostics.compatibility_error)
    return SimulationTrace(tuple(states), tuple(costs), tuple(active), tuple(errors))


def trivial_attractor_trace(state0: np.ndarray, eta: float, steps: int) -> SimulationTrace:
    state = state0.copy()
    mean = float(np.mean(state0))
    states = [state.copy()]
    corrections = []
    rho = min(0.25, max(0.02, eta / 5.0))
    for _ in range(steps):
        correction = rho * (mean - state)
        state = state + correction
        states.append(state.copy())
        corrections.append(float(np.sqrt(np.mean(correction**2))))
    zeros = tuple(0.0 for _ in range(steps))
    return SimulationTrace(tuple(states), tuple(corrections), zeros, zeros)


def oracle_trace(state0: np.ndarray, reference: np.ndarray, steps: int) -> SimulationTrace:
    state = state0.copy()
    states = [state.copy()]
    corrections = []
    for _ in range(steps):
        correction = 0.25 * (reference - state)
        state = state + correction
        states.append(state.copy())
        corrections.append(float(np.sqrt(np.mean(correction**2))))
    zeros = tuple(0.0 for _ in range(steps))
    return SimulationTrace(tuple(states), tuple(corrections), zeros, zeros)


def run_condition(
    condition: str,
    fixture: GridFixture,
    reference: np.ndarray,
    perturbed: np.ndarray,
    memory: RelationalMemory,
    eta: float,
    margin: float,
    steps: int,
    dropout_fraction: float,
    seed: int,
) -> SimulationTrace:
    if condition == "target":
        return simulate_feedback(perturbed, memory, eta, margin, steps, dropout_fraction, seed)
    if condition == "passive_relaxation":
        return static_trace(perturbed, steps)
    if condition in {"operator_only_filtering", "parameter_budget_matched"}:
        return operator_filter_trace(perturbed, fixture.incidence, eta, margin, steps)
    if condition == "topology_altered":
        incidence = rewired_incidence(fixture, seed + 1701)
        altered = RelationalMemory(incidence, memory.orientation.copy(), memory.valid_edges.copy(), memory.max_degree)
        return simulate_feedback(perturbed, altered, eta, margin, steps, dropout_fraction, seed)
    if condition == "trivial_attractor":
        return trivial_attractor_trace(perturbed, eta, steps)
    if condition == "identity_scrambled":
        rng = np.random.default_rng(seed + 1801)
        scrambled = RelationalMemory(
            memory.incidence.copy(),
            memory.orientation[rng.permutation(len(memory.orientation))],
            memory.valid_edges.copy(),
            memory.max_degree,
        )
        return simulate_feedback(perturbed, scrambled, eta, margin, steps, dropout_fraction, seed)
    if condition == "signal_blind":
        return signal_blind_trace(perturbed, memory, eta, margin, steps, dropout_fraction, seed + 100000)
    if condition == "memory_ablation":
        absent = RelationalMemory(memory.incidence.copy(), np.zeros_like(memory.orientation), np.zeros_like(memory.valid_edges), memory.max_degree)
        return simulate_feedback(perturbed, absent, eta, margin, steps, dropout_fraction, seed)
    if condition == "redundancy_ablation":
        incidence = spanning_tree_incidence(fixture)
        sparse = RelationalMemory.encode(incidence, reference)
        return simulate_feedback(perturbed, sparse, eta, margin, steps, dropout_fraction, seed)
    if condition == ORACLE_CONTROL:
        return oracle_trace(perturbed, reference, steps)
    raise ValueError(f"unknown condition: {condition}")
