from __future__ import annotations

import numpy as np

from .representation import COARSE_SIDE, coarse_field, expand_coarse


PERTURBATION_FAMILIES = (
    "coarse_cell_offsets",
    "smooth_physical_twist",
    "within_cell_balanced_gradient",
)


def coordinates(n_side: int) -> tuple[np.ndarray, np.ndarray]:
    axis = np.arange(int(n_side), dtype=float) / float(n_side)
    return np.meshgrid(axis, axis, indexing="ij")


def functional_probes(n_side: int) -> np.ndarray:
    xx, yy = coordinates(n_side)
    probes = np.vstack(
        [
            np.cos(2.0 * np.pi * xx).reshape(-1),
            np.sin(2.0 * np.pi * xx).reshape(-1),
            np.cos(2.0 * np.pi * yy).reshape(-1),
            np.sin(2.0 * np.pi * yy).reshape(-1),
        ]
    )
    probes /= np.linalg.norm(probes, axis=1, keepdims=True)
    return probes


def initial_state(n_side: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(int(seed))
    xx, yy = coordinates(n_side)
    coefficients = rng.normal(size=8)
    state = (
        coefficients[0] * np.cos(2.0 * np.pi * xx)
        + coefficients[1] * np.sin(2.0 * np.pi * xx)
        + coefficients[2] * np.cos(2.0 * np.pi * yy)
        + coefficients[3] * np.sin(2.0 * np.pi * yy)
        + 0.45 * coefficients[4] * np.cos(2.0 * np.pi * (xx + yy))
        + 0.45 * coefficients[5] * np.sin(2.0 * np.pi * (xx - yy))
        + 0.20 * coefficients[6] * np.cos(4.0 * np.pi * xx)
        + 0.20 * coefficients[7] * np.sin(4.0 * np.pi * yy)
    ).reshape(-1)
    state -= np.mean(state)
    state /= max(float(np.std(state)), 1.0e-12)
    return state


def _normalize(values: np.ndarray, magnitude: float) -> np.ndarray:
    centered = np.asarray(values, dtype=float) - float(np.mean(values))
    rms = float(np.sqrt(np.mean(centered**2)))
    if rms <= 1.0e-12:
        raise ValueError("degenerate perturbation")
    return float(magnitude) * centered / rms


def perturbation(n_side: int, family: str, magnitude: float, seed: int) -> np.ndarray:
    n = int(n_side)
    offset = PERTURBATION_FAMILIES.index(family) * 1009
    rng = np.random.default_rng(int(seed) + offset + n * 17)
    if family == "coarse_cell_offsets":
        selected = rng.choice(16, size=4, replace=False)
        coarse = np.zeros(16, dtype=float)
        signs = np.array([-1.0, -1.0, 1.0, 1.0], dtype=float)
        rng.shuffle(signs)
        coarse[selected] = signs
        values = expand_coarse(coarse, n)
    elif family == "smooth_physical_twist":
        xx, yy = coordinates(n)
        phase_x, phase_y = rng.uniform(0.0, 2.0 * np.pi, size=2)
        values = np.sin(2.0 * np.pi * xx + phase_x) + 0.65 * np.cos(2.0 * np.pi * yy + phase_y)
        values = values.reshape(-1)
    elif family == "within_cell_balanced_gradient":
        cell = n // COARSE_SIDE
        local = np.arange(cell, dtype=float) - 0.5 * (cell - 1)
        local_x, local_y = np.meshgrid(local, local, indexing="ij")
        values_2d = np.zeros((n, n), dtype=float)
        weights = rng.choice((-1.0, 1.0), size=(COARSE_SIDE, COARSE_SIDE, 2))
        for row in range(COARSE_SIDE):
            for col in range(COARSE_SIDE):
                patch = weights[row, col, 0] * local_x + weights[row, col, 1] * local_y
                patch -= np.mean(patch)
                values_2d[row * cell : (row + 1) * cell, col * cell : (col + 1) * cell] = patch
        values = values_2d.reshape(-1)
        if np.max(np.abs(coarse_field(values, n))) > 1.0e-12:
            raise ValueError("within-cell perturbation leaked into coarse means")
    else:
        raise ValueError(f"unknown perturbation family: {family}")
    return _normalize(values, magnitude)
