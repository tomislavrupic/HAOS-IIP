from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from experiments.emergence_ladder.rung3_recovery_trajectory_v2.mechanism import RelationalMemory, simulate_feedback

from .alphabet import RelationalAlphabet
from .bridge import BridgeTrace, simulate_bridge
from .decoder import DecodeResult, decode_local
from .parity_architecture import DistributedParityMemory


ADMISSIBLE_CONDITIONS = (
    "target",
    "rt02_frozen_baseline",
    "passive_relaxation",
    "operator_only",
    "random_parity_checks",
    "parity_permutation",
    "parity_deletion",
    "parity_bit_corruption",
    "memory_budget_random_bits",
    "signal_blind",
    "identity_scrambled",
    "trivial_attractor",
    "redundancy_ablation",
    "syndrome_ablation",
    "message_passing_ablation",
    "continuous_bridge_ablation",
    "identity_constraint_ablation",
)
NON_ADMISSIBLE_CONDITIONS = ("centralized_parity_non_admissible", "full_checkpoint_oracle_non_admissible")


@dataclass(frozen=True)
class ConditionResult:
    trace: BridgeTrace
    decode: DecodeResult | None
    admissible: bool
    mechanism_note: str


def static_trace(state: np.ndarray, steps: int) -> BridgeTrace:
    states = tuple(np.asarray(state, dtype=float).copy() for _ in range(steps + 1))
    return BridgeTrace(states, tuple(0.0 for _ in range(steps)), tuple(0 for _ in range(steps)))


def filter_trace(state0: np.ndarray, alphabet: RelationalAlphabet, eta: float, steps: int) -> BridgeTrace:
    state = np.asarray(state0, dtype=float).copy()
    laplacian = alphabet.incidence.T @ alphabet.incidence
    scale = max(float(np.linalg.norm(laplacian, ord=2)), 1.0)
    states = [state.copy()]
    costs = []
    for _ in range(steps):
        correction = -float(eta) * (laplacian @ state) / scale
        correction -= np.mean(correction)
        state = state + correction
        states.append(state.copy())
        costs.append(float(np.sqrt(np.mean(correction**2))))
    return BridgeTrace(tuple(states), tuple(costs), tuple(0 for _ in range(steps)))


def trivial_trace(state0: np.ndarray, eta: float, steps: int) -> BridgeTrace:
    state = np.asarray(state0, dtype=float).copy()
    mean = float(np.mean(state))
    states = [state.copy()]
    costs = []
    rho = min(0.25, max(0.02, eta / 5.0))
    for _ in range(steps):
        correction = rho * (mean - state)
        state = state + correction
        states.append(state.copy())
        costs.append(float(np.sqrt(np.mean(correction**2))))
    return BridgeTrace(tuple(states), tuple(costs), tuple(0 for _ in range(steps)))


def signal_blind_trace(target_trace: BridgeTrace, state0: np.ndarray, seed: int) -> BridgeTrace:
    rng = np.random.default_rng(seed)
    state = np.asarray(state0, dtype=float).copy()
    states = [state.copy()]
    for norm in target_trace.correction_norms:
        direction = rng.normal(size=len(state))
        direction -= np.mean(direction)
        rms = max(float(np.sqrt(np.mean(direction**2))), 1.0e-12)
        state = state + direction * (float(norm) / rms)
        states.append(state.copy())
    return BridgeTrace(tuple(states), target_trace.correction_norms, tuple(0 for _ in target_trace.correction_norms))


def oracle_trace(state0: np.ndarray, reference: np.ndarray, steps: int) -> BridgeTrace:
    state = np.asarray(state0, dtype=float).copy()
    states = [state.copy()]
    costs = []
    for _ in range(steps):
        correction = 0.25 * (reference - state)
        state = state + correction
        states.append(state.copy())
        costs.append(float(np.sqrt(np.mean(correction**2))))
    return BridgeTrace(tuple(states), tuple(costs), tuple(0 for _ in range(steps)))


def memory_with_parity(memory: DistributedParityMemory, values: np.ndarray) -> DistributedParityMemory:
    parity = np.asarray(values, dtype=np.uint8).reshape((-1, 2))
    return DistributedParityMemory(tuple(tuple(int(x) for x in row) for row in parity), memory.owner_by_block)


def run_condition(
    condition: str,
    reference: np.ndarray,
    perturbed: np.ndarray,
    alphabet: RelationalAlphabet,
    memory: DistributedParityMemory,
    eta: float,
    margin: float,
    steps: int,
    seed: int,
    target_check_mask: np.ndarray | None = None,
) -> ConditionResult:
    current = alphabet.encode(perturbed)
    base_decode = decode_local(current, memory, alphabet, target_check_mask)
    target_trace = simulate_bridge(perturbed, base_decode.target_symbols, alphabet, eta, margin, steps)
    if condition == "target":
        return ConditionResult(target_trace, base_decode, True, "local syndrome decode plus bounded relation feedback")
    if condition == "rt02_frozen_baseline":
        rt_memory = RelationalMemory.encode(alphabet.incidence, reference, 0.08)
        rt_trace = simulate_feedback(perturbed, rt_memory, 0.7, 0.05, 80, 0.2, seed)
        trace = BridgeTrace(rt_trace.states, rt_trace.corrections, tuple(0 for _ in rt_trace.corrections))
        return ConditionResult(trace, None, True, "unchanged RT-02 one-bit orientation mechanism")
    if condition == "passive_relaxation":
        return ConditionResult(static_trace(perturbed, steps), None, True, "no restorative mechanism")
    if condition == "operator_only":
        return ConditionResult(filter_trace(perturbed, alphabet, eta, steps), None, True, "matched bounded smoothing without parity")
    rng = np.random.default_rng(seed + 19101)
    if condition == "random_parity_checks":
        permutation = rng.permutation(len(current))
        inverse = np.argsort(permutation)
        random_decode = decode_local(current[permutation], memory, alphabet)
        target = random_decode.target_symbols[inverse]
        return ConditionResult(simulate_bridge(perturbed, target, alphabet, eta, margin, steps), random_decode, True, "degree-matched checks on unrelated variables")
    if condition == "parity_permutation":
        parity = memory.parity_array()[rng.permutation(len(memory.parity_by_block))]
        altered = memory_with_parity(memory, parity)
        decoded = decode_local(current, altered, alphabet)
        return ConditionResult(simulate_bridge(perturbed, decoded.target_symbols, alphabet, eta, margin, steps), decoded, True, "check-to-block assignment permuted")
    if condition == "parity_deletion":
        available = np.ones((len(alphabet.blocks), 2), dtype=bool)
        available[::4, 0] = False
        decoded = decode_local(current, memory, alphabet, available)
        return ConditionResult(simulate_bridge(perturbed, decoded.target_symbols, alphabet, eta, margin, steps), decoded, True, "one eighth of checks deleted by frozen pattern")
    if condition == "parity_bit_corruption":
        parity = memory.parity_array().copy()
        parity.reshape(-1)[::8] ^= 1
        altered = memory_with_parity(memory, parity)
        decoded = decode_local(current, altered, alphabet)
        return ConditionResult(simulate_bridge(perturbed, decoded.target_symbols, alphabet, eta, margin, steps), decoded, True, "one eighth of stored parity bits flipped")
    if condition == "memory_budget_random_bits":
        altered = memory_with_parity(memory, rng.integers(0, 2, size=(16, 2), dtype=np.uint8))
        decoded = decode_local(current, altered, alphabet)
        return ConditionResult(simulate_bridge(perturbed, decoded.target_symbols, alphabet, eta, margin, steps), decoded, True, "32 random bits replace structured parity")
    if condition == "centralized_parity_non_admissible":
        return ConditionResult(target_trace, base_decode, False, "global decoder upper bound; excluded from evidence")
    if condition == "signal_blind":
        return ConditionResult(signal_blind_trace(target_trace, perturbed, seed + 100000), base_decode, True, "target correction energy with random direction")
    if condition == "identity_scrambled":
        target = base_decode.target_symbols[rng.permutation(len(current))]
        return ConditionResult(simulate_bridge(perturbed, target, alphabet, eta, margin, steps), base_decode, True, "identity-bearing target relations scrambled")
    if condition == "trivial_attractor":
        return ConditionResult(trivial_trace(perturbed, eta, steps), None, True, "low-information fixed-state collapse")
    if condition == "full_checkpoint_oracle_non_admissible":
        return ConditionResult(oracle_trace(perturbed, reference, steps), None, False, "full-state checkpoint oracle; excluded from evidence")
    if condition in {"redundancy_ablation", "syndrome_ablation", "message_passing_ablation"}:
        return ConditionResult(static_trace(perturbed, steps), None, True, f"essential component removed: {condition}")
    if condition == "continuous_bridge_ablation":
        return ConditionResult(static_trace(perturbed, steps), base_decode, True, "bits decoded but continuous correction removed")
    if condition == "identity_constraint_ablation":
        target = np.zeros_like(current)
        return ConditionResult(simulate_bridge(perturbed, target, alphabet, eta, margin, steps), base_decode, True, "identity target replaced by uniform codeword")
    raise ValueError(f"unknown RP-01 condition: {condition}")
