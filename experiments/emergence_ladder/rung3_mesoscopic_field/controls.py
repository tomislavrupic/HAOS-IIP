from __future__ import annotations

import numpy as np

from experiments.emergence_ladder.rung3_recovery_trajectory_v2.fixtures import make_fixture
from experiments.emergence_ladder.rung3_recovery_trajectory_v2.mechanism import RelationalMemory, simulate_feedback

from .fixtures import functional_probes
from .representation import (
    COARSE_VALUES,
    QuantizerSpec,
    coarse_field,
    quantize,
    random_orthonormal_rows,
    restore_coarse_field,
    restore_linear_projection,
)


ADMISSIBLE_CONDITIONS = (
    "target_quantized_coarse_field",
    "passive_relaxation",
    "operator_filter_matched_cost",
    "rt02_frozen_one_bit_rule",
    "block_permutation",
    "phase_scrambling",
    "amplitude_only",
    "phase_only",
    "equal_budget_random_projection",
    "compression_12_target_blind",
    "compression_8_target_blind",
)
ORACLE_CONDITIONS = (
    "unquantized_coarse_field",
    "direct_four_functional_phasors",
    "full_microscopic_checkpoint",
)
ALL_CONDITIONS = ADMISSIBLE_CONDITIONS + ORACLE_CONDITIONS


def _spectral_variant(field: np.ndarray, condition: str, seed: int) -> np.ndarray:
    source = np.asarray(field, dtype=float).reshape(4, 4)
    spectrum = np.fft.rfft2(source)
    magnitude = np.abs(spectrum)
    phase = np.angle(spectrum)
    if condition == "phase_scrambling":
        rng = np.random.default_rng(int(seed) + 41001)
        randomized = rng.uniform(-np.pi, np.pi, size=spectrum.shape)
        randomized[0, 0] = phase[0, 0]
        altered = magnitude * np.exp(1j * randomized)
    elif condition == "amplitude_only":
        altered = magnitude.astype(complex)
    elif condition == "phase_only":
        altered = np.exp(1j * phase)
        altered[0, 0] = 0.0
    else:
        raise ValueError(f"unknown spectral condition: {condition}")
    return np.fft.irfft2(altered, s=(4, 4)).real.reshape(COARSE_VALUES)


def _operator_direction(state: np.ndarray, n_side: int) -> np.ndarray:
    grid = np.asarray(state, dtype=float).reshape(n_side, n_side)
    neighbor_mean = 0.25 * (
        np.roll(grid, 1, axis=0)
        + np.roll(grid, -1, axis=0)
        + np.roll(grid, 1, axis=1)
        + np.roll(grid, -1, axis=1)
    )
    direction = (neighbor_mean - grid).reshape(-1)
    direction -= np.mean(direction)
    return direction


def run_condition(
    condition: str,
    reference: np.ndarray,
    state_prime: np.ndarray,
    n_side: int,
    seed: int,
    quantizer: QuantizerSpec,
    target_correction_rms: float,
) -> tuple[np.ndarray, dict[str, object]]:
    reference = np.asarray(reference, dtype=float)
    state_prime = np.asarray(state_prime, dtype=float)
    source_field = coarse_field(reference, n_side)
    metadata: dict[str, object] = {
        "payload_bits": 0,
        "payload_values": 0,
        "saturation_count": 0,
        "quantization_max_abs_error": 0.0,
        "decoder_path": condition,
        "local_representation": False,
    }

    if condition == "passive_relaxation":
        return state_prime.copy(), metadata
    if condition == "target_quantized_coarse_field":
        stored = quantize(source_field, quantizer)
        metadata.update(payload_bits=16 * quantizer.bits, payload_values=16, saturation_count=stored.saturation_count, quantization_max_abs_error=stored.max_abs_error, local_representation=True)
        return restore_coarse_field(state_prime, stored.decoded, n_side), metadata
    if condition in {"block_permutation", "phase_scrambling", "amplitude_only", "phase_only"}:
        if condition == "block_permutation":
            rng = np.random.default_rng(int(seed) + 40001)
            altered_field = source_field[rng.permutation(COARSE_VALUES)]
        else:
            altered_field = _spectral_variant(source_field, condition, seed)
        stored = quantize(altered_field, quantizer)
        metadata.update(payload_bits=16 * quantizer.bits, payload_values=16, saturation_count=stored.saturation_count, quantization_max_abs_error=stored.max_abs_error, local_representation=True)
        return restore_coarse_field(state_prime, stored.decoded, n_side), metadata
    if condition in {"compression_12_target_blind", "compression_8_target_blind"}:
        count = 12 if "12" in condition else 8
        projection = random_orthonormal_rows(count, COARSE_VALUES, 12120 if count == 12 else 8080)
        stored = quantize(projection @ source_field, quantizer)
        decoded_field = projection.T @ stored.decoded
        metadata.update(payload_bits=count * quantizer.bits, payload_values=count, saturation_count=stored.saturation_count, quantization_max_abs_error=stored.max_abs_error, local_representation=True)
        return restore_coarse_field(state_prime, decoded_field, n_side), metadata
    if condition == "equal_budget_random_projection":
        projection = random_orthonormal_rows(COARSE_VALUES, n_side * n_side, 160016 + n_side)
        stored = quantize(projection @ reference, quantizer)
        metadata.update(payload_bits=16 * quantizer.bits, payload_values=16, saturation_count=stored.saturation_count, quantization_max_abs_error=stored.max_abs_error, local_representation=False)
        return restore_linear_projection(state_prime, projection, stored.decoded), metadata
    if condition == "operator_filter_matched_cost":
        direction = _operator_direction(state_prime, n_side)
        rms = float(np.sqrt(np.mean(direction**2)))
        scale = 0.0 if rms <= 1.0e-12 else float(target_correction_rms) / rms
        return state_prime + scale * direction, metadata
    if condition == "rt02_frozen_one_bit_rule":
        fixture = make_fixture(n_side)
        memory = RelationalMemory.encode(fixture.incidence, reference, minimum_edge_magnitude=0.08)
        trace = simulate_feedback(state_prime, memory, eta=0.7, margin=0.05, steps=80, dropout_fraction=0.2, seed=int(seed))
        metadata.update(payload_bits=int(np.sum(memory.valid_edges)), payload_values=int(np.sum(memory.valid_edges)), local_representation=True)
        return trace.states[-1].copy(), metadata
    if condition == "unquantized_coarse_field":
        metadata.update(payload_bits=16 * 64, payload_values=16, local_representation=True)
        return restore_coarse_field(state_prime, source_field, n_side), metadata
    if condition == "direct_four_functional_phasors":
        probes = functional_probes(n_side)
        target = probes @ reference
        return restore_linear_projection(state_prime, probes, target), metadata
    if condition == "full_microscopic_checkpoint":
        metadata.update(payload_bits=n_side * n_side * 64, payload_values=n_side * n_side)
        return reference.copy(), metadata
    raise ValueError(f"unknown condition: {condition}")
