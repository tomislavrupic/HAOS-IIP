#!/usr/bin/env python3
"""Minimal HAOS-GEN generation -> stress -> selection demonstration."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haos_gen import TuningConfig, expand_spectral_addresses
from telemetry.frozen_metrics import SurvivalThresholds


def path_laplacian(size: int) -> np.ndarray:
    adjacency = np.zeros((size, size), dtype=float)
    for index in range(size - 1):
        adjacency[index, index + 1] = 1.0
        adjacency[index + 1, index] = 1.0
    return np.diag(np.sum(adjacency, axis=1)) - adjacency


def perturbation_ladder(operator: np.ndarray) -> tuple[tuple[np.ndarray, ...], tuple[float, ...]]:
    levels = (0.00, 0.03, 0.06, 0.10, 0.15, 0.22)
    size = operator.shape[0]
    pattern = np.zeros_like(operator)
    for index in range(size - 1):
        weight = (-1.0 if index % 2 else 1.0) * (1.0 + index / size)
        pattern[index, index + 1] = weight
        pattern[index + 1, index] = weight
    pattern += np.diag(np.linspace(-0.5, 0.5, size))
    return tuple(operator + level * pattern for level in levels), levels


def main() -> None:
    operator = path_laplacian(12)
    _, eigenvectors = np.linalg.eigh(operator)
    # An admitted but imperfect address: mostly mode 3 with a small mode-2
    # component. Generation must decide whether either nearby coordinate is a
    # recoverable expansion rather than assuming both are improvements.
    seed = 0.995 * eigenvectors[:, 3] + np.sqrt(1.0 - 0.995**2) * eigenvectors[:, 2]
    perturbations, levels = perturbation_ladder(operator)
    thresholds = SurvivalThresholds(
        max_width_growth=0.35,
        min_concentration=0.55,
        max_participation_growth=0.40,
        min_overlap=0.72,
        min_recovery_score=0.45,
    )
    audit = expand_spectral_addresses(
        operator,
        seed,
        perturbations,
        levels,
        thresholds,
        TuningConfig(
            max_candidates=6,
            min_resonance=0.05,
            min_novelty=0.001,
            control_margin=0.18,
        ),
    )
    print(json.dumps(audit.to_jsonable(), indent=2))


if __name__ == "__main__":
    main()
