#!/usr/bin/env python3
"""Reproducible HAOS-GEN V0.3 proposal-operator tournament."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haos_gen import (
    Regime,
    TournamentConfig,
    TuningConfig,
    ValidationConfig,
    run_proposal_tournament,
)
from telemetry.frozen_metrics import SurvivalThresholds


def path_laplacian(size: int) -> np.ndarray:
    adjacency = np.zeros((size, size), dtype=float)
    for index in range(size - 1):
        adjacency[index, index + 1] = adjacency[index + 1, index] = 1.0
    return np.diag(np.sum(adjacency, axis=1)) - adjacency


def perturbation_family(
    operator: np.ndarray,
    levels: tuple[float, ...],
    family: str,
) -> tuple[np.ndarray, ...]:
    size = operator.shape[0]
    pattern = np.zeros_like(operator)
    if family == "edge":
        for index in range(size - 1):
            weight = (-1.0 if index % 2 else 1.0) * (1.0 + index / size)
            pattern[index, index + 1] = pattern[index + 1, index] = weight
        pattern += np.diag(np.linspace(-0.5, 0.5, size))
    elif family == "long_range":
        for index in range(size - 2):
            weight = (-1.0 if index % 3 == 0 else 1.0) * (0.6 + index / size)
            pattern[index, index + 2] = pattern[index + 2, index] = weight
        pattern += np.diag(np.cos(np.arange(size, dtype=float) * 0.7))
    elif family == "hard_disorder":
        for index in range(size - 1):
            weight = (1.0 if index % 4 in (1, 2) else -1.0) * (1.4 + index / size)
            pattern[index, index + 1] = pattern[index + 1, index] = weight
        for index in range(size - 3):
            pattern[index, index + 3] = pattern[index + 3, index] = 0.75 * (-1.0) ** index
        pattern += np.diag(1.5 * np.sin(np.arange(size, dtype=float) * 1.1))
    else:
        raise ValueError(f"unknown perturbation family: {family}")
    return tuple(operator + level * pattern for level in levels)


def admitted_seed_set(operator: np.ndarray) -> dict[str, np.ndarray]:
    """Five declared, normalized spectral addresses with broad mode support."""

    _, eigenvectors = np.linalg.eigh(operator)
    mode_indices = np.arange(1, operator.shape[0])
    seeds: dict[str, np.ndarray] = {}
    for seed_index in range(5):
        coefficients = (
            np.cos((seed_index + 1) * mode_indices * 0.37)
            + 0.55 * np.sin((seed_index + 2) * mode_indices * 0.23)
        )
        coefficients *= np.exp(-0.055 * mode_indices)
        vector = eigenvectors[:, 1:] @ coefficients
        seeds[f"seed_{seed_index:02d}"] = vector / np.linalg.norm(vector)
    return seeds


def main() -> None:
    operator = path_laplacian(18)
    tuning_levels = (0.00, 0.03, 0.06, 0.10, 0.15, 0.22)
    heldout_levels = (0.00, 0.025, 0.05, 0.085, 0.13, 0.19)
    hard_levels = (0.00, 0.05, 0.10, 0.17, 0.25, 0.36)
    regimes = (
        Regime(
            name="standard",
            tuning_operators=perturbation_family(operator, tuning_levels, "edge"),
            tuning_levels=tuning_levels,
            heldout_operators=perturbation_family(
                operator, heldout_levels, "long_range"
            ),
            heldout_levels=heldout_levels,
        ),
        Regime(
            name="hard",
            tuning_operators=perturbation_family(operator, tuning_levels, "edge"),
            tuning_levels=tuning_levels,
            heldout_operators=perturbation_family(
                operator, hard_levels, "hard_disorder"
            ),
            heldout_levels=hard_levels,
        ),
    )
    thresholds = SurvivalThresholds(
        max_width_growth=0.35,
        min_concentration=0.55,
        max_participation_growth=0.40,
        min_overlap=0.72,
        min_recovery_score=0.45,
    )
    audit = run_proposal_tournament(
        operator=operator,
        seed_addresses=admitted_seed_set(operator),
        regimes=regimes,
        thresholds=thresholds,
        tuning_config=TuningConfig(
            max_candidates=4,
            min_resonance=0.015,
            min_novelty=0.001,
            control_margin=0.03,
        ),
        validation_config=ValidationConfig(
            min_heldout_mean_gain=0.01,
            min_k_star_displacement=1,
            hostile_control_margin=0.03,
        ),
        config=TournamentConfig(
            proposal_budget=4,
            min_spectral_separation=0.18,
            multi_seed_k=2,
            candidate_seed_overlap_limit=0.94,
            walk_length=4,
            walk_extraction_radius=1,
            address_distance_threshold=0.08,
            subspace_overlap_threshold=0.92,
            minimum_positive_seeds=2,
            hard_regime_name="hard",
        ),
    )
    print(json.dumps(audit.to_jsonable(), indent=2))


if __name__ == "__main__":
    main()
