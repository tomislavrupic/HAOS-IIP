#!/usr/bin/env python3
"""HAOS-GEN V0.2: tuning versus sealed held-out hostile validation."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haos_gen import TuningConfig, ValidationConfig, validate_expansion_hostile
from telemetry.frozen_metrics import SurvivalThresholds


def path_laplacian(size: int) -> np.ndarray:
    adjacency = np.zeros((size, size), dtype=float)
    for index in range(size - 1):
        adjacency[index, index + 1] = adjacency[index + 1, index] = 1.0
    return np.diag(np.sum(adjacency, axis=1)) - adjacency


def edge_weight_family(
    operator: np.ndarray,
) -> tuple[tuple[np.ndarray, ...], tuple[float, ...]]:
    levels = (0.00, 0.03, 0.06, 0.10, 0.15, 0.22)
    size = operator.shape[0]
    pattern = np.zeros_like(operator)
    for index in range(size - 1):
        weight = (-1.0 if index % 2 else 1.0) * (1.0 + index / size)
        pattern[index, index + 1] = pattern[index + 1, index] = weight
    pattern += np.diag(np.linspace(-0.5, 0.5, size))
    return tuple(operator + level * pattern for level in levels), levels


def sealed_long_range_family(
    operator: np.ndarray,
) -> tuple[tuple[np.ndarray, ...], tuple[float, ...]]:
    """Structurally distinct family not consumed by candidate generation."""

    levels = (0.00, 0.025, 0.05, 0.085, 0.13, 0.19)
    size = operator.shape[0]
    pattern = np.zeros_like(operator)
    for index in range(size - 2):
        weight = (-1.0 if index % 3 == 0 else 1.0) * (0.6 + index / size)
        pattern[index, index + 2] = pattern[index + 2, index] = weight
    pattern += np.diag(np.cos(np.arange(size, dtype=float) * 0.7))
    return tuple(operator + level * pattern for level in levels), levels


def main() -> None:
    operator = path_laplacian(16)
    _, eigenvectors = np.linalg.eigh(operator)
    seed = (
        0.82 * eigenvectors[:, 3]
        + 0.38 * eigenvectors[:, 2]
        + 0.28 * eigenvectors[:, 5]
        + 0.22 * eigenvectors[:, 7]
    )
    seed /= np.linalg.norm(seed)
    tuning_operators, tuning_levels = edge_weight_family(operator)
    heldout_operators, heldout_levels = sealed_long_range_family(operator)
    thresholds = SurvivalThresholds(
        max_width_growth=0.35,
        min_concentration=0.55,
        max_participation_growth=0.40,
        min_overlap=0.72,
        min_recovery_score=0.45,
    )
    audit = validate_expansion_hostile(
        operator,
        seed,
        tuning_operators,
        tuning_levels,
        heldout_operators,
        heldout_levels,
        thresholds,
        TuningConfig(
            max_candidates=8,
            min_resonance=0.05,
            min_novelty=0.001,
            control_margin=0.08,
        ),
        ValidationConfig(
            min_heldout_mean_gain=0.01,
            min_k_star_displacement=1,
            hostile_control_margin=0.03,
        ),
    )
    print(json.dumps(audit.to_jsonable(), indent=2))


if __name__ == "__main__":
    main()
