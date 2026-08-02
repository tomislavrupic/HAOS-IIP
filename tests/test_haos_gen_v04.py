from __future__ import annotations

import unittest

import numpy as np

from haos_gen import (
    EdgeInterventionConfig,
    Regime,
    TuningConfig,
    ValidationConfig,
    run_edge_intervention_experiment,
    select_directed_edge_updates,
    interpret_result_status,
)
from telemetry.frozen_metrics import SurvivalThresholds


def path_laplacian(size: int, scale: float = 1.0) -> np.ndarray:
    adjacency = np.zeros((size, size), dtype=float)
    for index in range(size - 1):
        weight = scale * (1.0 + 0.03 * index)
        adjacency[index, index + 1] = adjacency[index + 1, index] = weight
    return np.diag(np.sum(adjacency, axis=1)) - adjacency


def family(operator: np.ndarray, levels: tuple[float, ...], stride: int) -> tuple[np.ndarray, ...]:
    pattern = np.zeros_like(operator)
    for index in range(operator.shape[0] - stride):
        pattern[index, index + stride] = pattern[index + stride, index] = (-1.0) ** index
    return tuple(operator + level * pattern for level in levels)


def seeds(operator: np.ndarray) -> dict[str, np.ndarray]:
    _, eigenvectors = np.linalg.eigh(operator)
    indices = np.arange(1, operator.shape[0])
    result = {}
    for seed_index in range(3):
        coefficients = np.cos((seed_index + 1) * indices * 0.29)
        coefficients += 0.35 * np.sin((seed_index + 2) * indices * 0.17)
        vector = eigenvectors[:, 1:] @ coefficients
        result[f"seed_{seed_index}"] = vector / np.linalg.norm(vector)
    return result


class HaosGenV04Tests(unittest.TestCase):
    def test_clean_negative_is_not_execution_failure(self) -> None:
        result = interpret_result_status(
            "OPEN_NO_STRUCTURAL_INTERVENTION_ADVANTAGE"
        )
        self.assertEqual(result.execution_health, "VERIFIED")
        self.assertEqual(result.evidence_outcome, "NEGATIVE")
        self.assertEqual(result.display_tone, "AMBER")
        self.assertNotEqual(result.display_tone, "RED")
        paused = interpret_result_status(
            "VERIFIED_NEGATIVE_SYNTHETIC_GENERATIVE_LINE_PAUSED"
        )
        self.assertEqual(paused.execution_health, "VERIFIED")
        self.assertEqual(paused.research_state, "TERMINAL_BOUNDARY")
        self.assertEqual(paused.display_tone, "BLUE")

    def test_directed_selector_is_deterministic_and_budget_bounded(self) -> None:
        operator = path_laplacian(10)
        levels = (0.0, 0.05, 0.1)
        tuning = family(operator, levels, 1)
        thresholds = SurvivalThresholds(0.5, 0.4, 0.6, 0.6, 0.35)
        tuning_config = TuningConfig(max_candidates=2, min_resonance=0.01)
        config = EdgeInterventionConfig(edge_budget=2, max_weight_change=0.05)

        first = select_directed_edge_updates(
            operator, seeds(operator), tuning, levels, thresholds, tuning_config, config
        )
        second = select_directed_edge_updates(
            operator, seeds(operator), tuning, levels, thresholds, tuning_config, config
        )

        self.assertEqual(first, second)
        self.assertEqual(len(first), config.edge_budget)
        self.assertLessEqual(
            sum(abs(update.delta_weight) for update in first),
            config.edge_budget * config.max_weight_change + 1.0e-12,
        )
        self.assertEqual(len({(update.src, update.dst) for update in first}), len(first))

    def test_full_matrix_is_fail_closed_and_reversible(self) -> None:
        levels = (0.0, 0.06, 0.12)
        substrates = {}
        regimes = {}
        for index in range(3):
            substrate_id = f"graph_{index}"
            operator = path_laplacian(10, 1.0 + 0.04 * index)
            substrates[substrate_id] = (operator, seeds(operator))
            regimes[substrate_id] = (
                Regime("standard", family(operator, levels, 1), levels, family(operator, levels, 2), levels),
                Regime("hard", family(operator, levels, 1), levels, family(operator, levels, 3), levels),
            )
        audit = run_edge_intervention_experiment(
            substrates,
            regimes,
            SurvivalThresholds(0.5, 0.4, 0.6, 0.6, 0.35),
            TuningConfig(
                max_candidates=2,
                min_resonance=0.01,
                min_novelty=0.001,
                control_margin=0.01,
            ),
            ValidationConfig(hostile_control_margin=0.01),
            EdgeInterventionConfig(
                edge_budget=1,
                max_weight_change=0.04,
                random_trials=1,
                minimum_expansion_substrates=2,
            ),
        )

        self.assertEqual(len(audit.ledger), 6)
        self.assertTrue(all(row.edges_touched == 1 for row in audit.ledger))
        self.assertTrue(all(row.reversal_ok for row in audit.ledger))
        self.assertIn(
            audit.aggregate["status"],
            {
                "STRUCTURAL_EXPANSION_CONFIRMED",
                "OPEN_NO_STRUCTURAL_INTERVENTION_ADVANTAGE",
            },
        )
        self.assertFalse(audit.contract["selector_uses_heldout"])
        self.assertTrue(audit.contract["random_uses_identical_signed_magnitude_multiset"])

    def test_requires_three_substrates(self) -> None:
        operator = path_laplacian(6)
        levels = (0.0, 0.1)
        regime = Regime("hard", (operator, operator), levels, (operator, operator), levels)
        with self.assertRaisesRegex(ValueError, "at least three"):
            run_edge_intervention_experiment(
                {"one": (operator, seeds(operator))},
                {"one": (regime,)},
                SurvivalThresholds(0.5, 0.4, 0.6, 0.6, 0.35),
                TuningConfig(),
                ValidationConfig(),
            )


if __name__ == "__main__":
    unittest.main()
