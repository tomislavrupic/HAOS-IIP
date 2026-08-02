#!/usr/bin/env python3
"""Canonical HAOS-GEN V0.4 Class-1 edge-intervention experiment."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haos_gen import (
    EdgeInterventionConfig,
    Regime,
    TuningConfig,
    ValidationConfig,
    run_edge_intervention_experiment,
)
from telemetry.frozen_metrics import SurvivalThresholds


def weighted_path_laplacian(size: int, realization: int) -> np.ndarray:
    adjacency = np.zeros((size, size), dtype=float)
    for index in range(size - 1):
        weight = 1.0 + 0.12 * np.sin((realization + 1) * (index + 1) * 0.61)
        adjacency[index, index + 1] = adjacency[index + 1, index] = weight
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
            weight = (-1.0 if index % 2 else 1.0) * (0.8 + index / size)
            pattern[index, index + 1] = pattern[index + 1, index] = weight
        pattern += np.diag(np.linspace(-0.4, 0.4, size))
    elif family == "long_range":
        for index in range(size - 2):
            weight = (-1.0 if index % 3 == 0 else 1.0) * (0.55 + index / size)
            pattern[index, index + 2] = pattern[index + 2, index] = weight
        pattern += np.diag(np.cos(np.arange(size, dtype=float) * 0.73))
    elif family == "hard":
        for index in range(size - 1):
            weight = (-1.0) ** index * (1.3 + index / size)
            pattern[index, index + 1] = pattern[index + 1, index] = weight
        for index in range(size - 3):
            pattern[index, index + 3] = pattern[index + 3, index] = 0.7 * (-1.0) ** index
        pattern += np.diag(1.4 * np.sin(np.arange(size, dtype=float) * 1.13))
    else:
        raise ValueError(f"unknown family: {family}")
    return tuple(operator + level * pattern for level in levels)


def admitted_seeds(operator: np.ndarray, realization: int) -> dict[str, np.ndarray]:
    _, eigenvectors = np.linalg.eigh(operator)
    indices = np.arange(1, operator.shape[0])
    seeds: dict[str, np.ndarray] = {}
    for seed_index in range(4):
        coefficients = np.cos(
            (seed_index + 1) * indices * (0.31 + 0.01 * realization)
        )
        coefficients += 0.45 * np.sin((seed_index + 2) * indices * 0.21)
        coefficients *= np.exp(-0.06 * indices)
        address = eigenvectors[:, 1:] @ coefficients
        seeds[f"seed_{seed_index:02d}"] = address / np.linalg.norm(address)
    return seeds


def main() -> None:
    size = 14
    tuning_levels = (0.00, 0.04, 0.08, 0.13, 0.20)
    standard_levels = (0.00, 0.035, 0.07, 0.12, 0.18)
    hard_levels = (0.00, 0.06, 0.12, 0.20, 0.30)
    substrates = {}
    regimes = {}
    for realization in range(3):
        substrate_id = f"weighted_path_{realization:02d}"
        operator = weighted_path_laplacian(size, realization)
        substrates[substrate_id] = (operator, admitted_seeds(operator, realization))
        regimes[substrate_id] = (
            Regime(
                name="standard",
                tuning_operators=perturbation_family(operator, tuning_levels, "edge"),
                tuning_levels=tuning_levels,
                heldout_operators=perturbation_family(
                    operator, standard_levels, "long_range"
                ),
                heldout_levels=standard_levels,
            ),
            Regime(
                name="hard",
                tuning_operators=perturbation_family(operator, tuning_levels, "edge"),
                tuning_levels=tuning_levels,
                heldout_operators=perturbation_family(operator, hard_levels, "hard"),
                heldout_levels=hard_levels,
            ),
        )
    audit = run_edge_intervention_experiment(
        substrates,
        regimes,
        SurvivalThresholds(
            max_width_growth=0.35,
            min_concentration=0.55,
            max_participation_growth=0.40,
            min_overlap=0.72,
            min_recovery_score=0.45,
        ),
        TuningConfig(
            max_candidates=3,
            min_resonance=0.015,
            min_novelty=0.001,
            control_margin=0.03,
        ),
        ValidationConfig(
            min_heldout_mean_gain=0.01,
            min_k_star_displacement=1,
            hostile_control_margin=0.03,
        ),
        EdgeInterventionConfig(
            edge_budget=2,
            max_weight_change=0.08,
            retention_tolerance=0.02,
            reversal_tolerance=1.0e-10,
            subspace_equivalence_threshold=0.92,
            random_trials=4,
            minimum_expansion_substrates=2,
            hard_regime_name="hard",
        ),
    )
    print(json.dumps(audit.to_jsonable(), indent=2))


if __name__ == "__main__":
    main()
