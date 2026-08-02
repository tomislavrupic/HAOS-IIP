from __future__ import annotations

import unittest

import numpy as np

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


def family(operator: np.ndarray, levels: tuple[float, ...], stride: int) -> tuple[np.ndarray, ...]:
    pattern = np.zeros_like(operator)
    for index in range(operator.shape[0] - stride):
        pattern[index, index + stride] = pattern[index + stride, index] = (-1.0) ** index
    pattern += np.diag(np.linspace(-0.3, 0.3, operator.shape[0]))
    return tuple(operator + level * pattern for level in levels)


class HaosGenV03Tests(unittest.TestCase):
    def test_tournament_enforces_budget_and_fail_closed_status(self) -> None:
        operator = path_laplacian(12)
        _, eigenvectors = np.linalg.eigh(operator)
        mode_indices = np.arange(1, 12)
        seeds = {}
        for seed_index in range(3):
            coefficients = np.cos((seed_index + 1) * mode_indices * 0.31)
            coefficients += 0.4 * np.sin((seed_index + 2) * mode_indices * 0.19)
            vector = eigenvectors[:, 1:] @ coefficients
            seeds[f"seed_{seed_index}"] = vector / np.linalg.norm(vector)
        levels = (0.0, 0.05, 0.10, 0.18)
        regimes = (
            Regime("standard", family(operator, levels, 1), levels, family(operator, levels, 2), levels),
            Regime("hard", family(operator, levels, 1), levels, family(operator, levels, 3), levels),
        )
        budget = 2
        audit = run_proposal_tournament(
            operator,
            seeds,
            regimes,
            SurvivalThresholds(0.5, 0.4, 0.6, 0.6, 0.35),
            TuningConfig(
                max_candidates=budget,
                min_resonance=0.01,
                min_novelty=0.001,
                control_margin=0.01,
            ),
            ValidationConfig(hostile_control_margin=0.01),
            TournamentConfig(
                proposal_budget=budget,
                min_spectral_separation=0.12,
                candidate_seed_overlap_limit=0.98,
                walk_length=3,
                address_distance_threshold=0.05,
                subspace_overlap_threshold=0.95,
            ),
        )

        expected_per_operator = len(seeds) * len(regimes) * budget
        for operator_name, summary in audit.aggregate.items():
            if operator_name in {"status", "judge"}:
                continue
            self.assertEqual(summary["proposal_count"], expected_per_operator)
        self.assertIn(
            audit.aggregate["status"],
            {
                "OPEN_NO_UNIQUE_DIRECTED_ADVANTAGE",
                "PASS_BOUNDED_UNIQUE_DIRECTED_ADVANTAGE",
            },
        )
        self.assertTrue(audit.contract["uniqueness"]["both_conditions_required"])

    def test_tournament_requires_three_seeds_and_hard_regime(self) -> None:
        operator = path_laplacian(6)
        levels = (0.0, 0.1)
        regime = Regime("standard", (operator, operator), levels, (operator, operator), levels)
        thresholds = SurvivalThresholds(0.5, 0.4, 0.6, 0.6, 0.35)
        with self.assertRaisesRegex(ValueError, "at least three"):
            run_proposal_tournament(
                operator,
                {"a": np.ones(6), "b": np.arange(1, 7)},
                [regime],
                thresholds,
                TuningConfig(),
                ValidationConfig(),
            )


if __name__ == "__main__":
    unittest.main()
