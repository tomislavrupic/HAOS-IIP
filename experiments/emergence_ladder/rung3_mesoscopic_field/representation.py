from __future__ import annotations

from dataclasses import dataclass

import numpy as np


COARSE_SIDE = 4
COARSE_VALUES = COARSE_SIDE * COARSE_SIDE


@dataclass(frozen=True)
class QuantizerSpec:
    bits: int = 8
    minimum: float = -3.0
    maximum: float = 3.0

    @property
    def code_count(self) -> int:
        return 1 << self.bits

    @property
    def step(self) -> float:
        return (self.maximum - self.minimum) / (self.code_count - 1)


@dataclass(frozen=True)
class QuantizedField:
    codes: np.ndarray
    decoded: np.ndarray
    saturation_count: int
    max_abs_error: float


def _validate_n_side(n_side: int) -> int:
    n = int(n_side)
    if n <= 0 or n % COARSE_SIDE != 0:
        raise ValueError("n_side must be a positive multiple of four")
    return n


def coarse_field(state: np.ndarray, n_side: int) -> np.ndarray:
    n = _validate_n_side(n_side)
    values = np.asarray(state, dtype=float)
    if values.shape != (n * n,):
        raise ValueError(f"state must have shape {(n * n,)}")
    cell = n // COARSE_SIDE
    return values.reshape(COARSE_SIDE, cell, COARSE_SIDE, cell).mean(axis=(1, 3)).reshape(COARSE_VALUES)


def expand_coarse(values: np.ndarray, n_side: int) -> np.ndarray:
    n = _validate_n_side(n_side)
    coarse = np.asarray(values, dtype=float).reshape(COARSE_SIDE, COARSE_SIDE)
    cell = n // COARSE_SIDE
    return np.repeat(np.repeat(coarse, cell, axis=0), cell, axis=1).reshape(n * n)


def coarse_matrix(n_side: int) -> np.ndarray:
    n = _validate_n_side(n_side)
    matrix = np.zeros((COARSE_VALUES, n * n), dtype=float)
    for index in range(COARSE_VALUES):
        basis = np.zeros(COARSE_VALUES, dtype=float)
        basis[index] = 1.0
        expanded = expand_coarse(basis, n)
        cell_count = int(np.sum(expanded))
        matrix[index] = expanded / cell_count
    return matrix


def quantize(values: np.ndarray, spec: QuantizerSpec) -> QuantizedField:
    array = np.asarray(values, dtype=float)
    if not np.isfinite(array).all():
        raise ValueError("non-finite quantizer input")
    if spec.bits <= 0 or spec.minimum >= spec.maximum:
        raise ValueError("invalid quantizer specification")
    saturation_count = int(np.sum((array < spec.minimum) | (array > spec.maximum)))
    clipped = np.clip(array, spec.minimum, spec.maximum)
    scaled = (clipped - spec.minimum) / (spec.maximum - spec.minimum) * (spec.code_count - 1)
    lower = np.floor(scaled)
    fraction = scaled - lower
    ties = np.isclose(fraction, 0.5, rtol=0.0, atol=1.0e-12)
    rounded = lower + (fraction > 0.5)
    rounded = np.where(ties, lower + np.mod(lower, 2.0), rounded)
    codes = rounded.astype(np.uint16)
    decoded = spec.minimum + codes.astype(float) * spec.step
    return QuantizedField(
        codes=codes,
        decoded=decoded,
        saturation_count=saturation_count,
        max_abs_error=float(np.max(np.abs(decoded - clipped))) if decoded.size else 0.0,
    )


def restore_coarse_field(state_prime: np.ndarray, target_coarse: np.ndarray, n_side: int) -> np.ndarray:
    current = coarse_field(state_prime, n_side)
    correction = expand_coarse(np.asarray(target_coarse, dtype=float) - current, n_side)
    restored = np.asarray(state_prime, dtype=float) + correction
    if not np.isfinite(restored).all():
        raise ValueError("non-finite restored state")
    return restored


def random_orthonormal_rows(row_count: int, column_count: int, seed: int) -> np.ndarray:
    if not 0 < row_count <= column_count:
        raise ValueError("row_count must be in 1..column_count")
    rng = np.random.default_rng(int(seed))
    raw = rng.normal(size=(column_count, row_count))
    q, _ = np.linalg.qr(raw, mode="reduced")
    return q.T


def restore_linear_projection(state_prime: np.ndarray, matrix: np.ndarray, target_values: np.ndarray) -> np.ndarray:
    state = np.asarray(state_prime, dtype=float)
    projection = np.asarray(matrix, dtype=float)
    target = np.asarray(target_values, dtype=float)
    gram = projection @ projection.T
    correction = projection.T @ np.linalg.solve(gram, target - projection @ state)
    restored = state + correction
    if not np.isfinite(restored).all():
        raise ValueError("non-finite projected state")
    return restored


def rank_and_nullity(n_side: int) -> tuple[int, int]:
    n = _validate_n_side(n_side)
    rank = int(np.linalg.matrix_rank(coarse_matrix(n)))
    return rank, n * n - rank
