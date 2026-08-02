from __future__ import annotations

import unittest
from hashlib import sha256
from pathlib import Path

import numpy as np

from haos_gen import (
    TuningConfig,
    ValidationConfig,
    expand_spectral_addresses,
    first_sustained_collapse,
    validate_expansion_hostile,
)
from telemetry.frozen_metrics import SurvivalThresholds


def path_laplacian(size: int) -> np.ndarray:
    adjacency = np.zeros((size, size), dtype=float)
    for index in range(size - 1):
        adjacency[index, index + 1] = adjacency[index + 1, index] = 1.0
    return np.diag(np.sum(adjacency, axis=1)) - adjacency


class HaosGenTests(unittest.TestCase):
    def test_v0_2_judge_is_byte_frozen_for_v0_3(self) -> None:
        root = Path(__file__).resolve().parents[1]
        expected = {
            root / "haos_gen" / "hostile_validation.py": (
                "7409b6b96cc213fcaf30237bb5686bcd146be0c0effd4ee89084890a2fa56d8c"
            ),
            root / "haos_gen" / "spectral_tuning.py": (
                "e69eb0123159fb2a94913e47cad6472a04f340157d384d7841f7a89ab33df063"
            ),
            root / "telemetry" / "frozen_metrics.py": (
                "daa7a759ad8f0ed67249f194a7e7e0753c78a4eb4c0a8b44d47e2e7a0d121f87"
            ),
        }
        for path, expected_digest in expected.items():
            with self.subTest(path=path.name):
                self.assertEqual(sha256(path.read_bytes()).hexdigest(), expected_digest)

    def test_first_sustained_collapse_requires_full_window(self) -> None:
        self.assertEqual(first_sustained_collapse([0.8, 0.4, 0.3], 0.5, 2), 1)
        self.assertIsNone(first_sustained_collapse([0.8, 0.4, 0.7], 0.5, 2))

    def test_directed_generation_is_deterministic_and_audited(self) -> None:
        operator = path_laplacian(10)
        seed = np.zeros(10)
        seed[2:5] = (0.4, 1.0, 0.6)
        diagonal = np.diag(np.linspace(-0.5, 0.5, 10))
        levels = (0.0, 0.05, 0.10, 0.20)
        perturbations = tuple(operator + level * diagonal for level in levels)
        thresholds = SurvivalThresholds(0.5, 0.4, 0.6, 0.6, 0.35)
        config = TuningConfig(max_candidates=4, min_resonance=0.05)

        first = expand_spectral_addresses(operator, seed, perturbations, levels, thresholds, config)
        second = expand_spectral_addresses(operator, seed, perturbations, levels, thresholds, config)

        self.assertEqual(first.to_jsonable(), second.to_jsonable())
        self.assertGreater(len(first.candidates), 0)
        self.assertTrue(all(candidate.null_profile.scores for candidate in first.candidates))
        self.assertEqual(first.contract["mutates_frozen_state"], False)

    def test_rejects_shape_mismatch(self) -> None:
        operator = path_laplacian(5)
        thresholds = SurvivalThresholds(0.5, 0.4, 0.6, 0.6, 0.35)
        with self.assertRaisesRegex(ValueError, "match baseline shape"):
            expand_spectral_addresses(
                operator,
                np.ones(5),
                [np.eye(4)],
                [0.0],
                thresholds,
            )

    def test_hostile_validation_is_equal_budget_and_deterministic(self) -> None:
        operator = path_laplacian(12)
        _, eigenvectors = np.linalg.eigh(operator)
        seed = 0.8 * eigenvectors[:, 2] + 0.5 * eigenvectors[:, 4]
        seed /= np.linalg.norm(seed)
        levels = (0.0, 0.05, 0.10, 0.18)
        diagonal = np.diag(np.linspace(-0.4, 0.4, 12))
        tuning = tuple(operator + level * diagonal for level in levels)
        long_range = np.zeros_like(operator)
        for index in range(10):
            long_range[index, index + 2] = long_range[index + 2, index] = (-1.0) ** index
        heldout = tuple(operator + level * long_range for level in levels)
        thresholds = SurvivalThresholds(0.5, 0.4, 0.6, 0.6, 0.35)
        tuning_config = TuningConfig(
            max_candidates=5,
            min_resonance=0.05,
            min_novelty=0.001,
        )
        validation_config = ValidationConfig(hostile_control_margin=0.02)

        first = validate_expansion_hostile(
            operator,
            seed,
            tuning,
            levels,
            heldout,
            levels,
            thresholds,
            tuning_config,
            validation_config,
        )
        second = validate_expansion_hostile(
            operator,
            seed,
            tuning,
            levels,
            heldout,
            levels,
            thresholds,
            tuning_config,
            validation_config,
        )

        self.assertEqual(first.to_jsonable(), second.to_jsonable())
        self.assertEqual(
            first.yield_summary["proposed_count"],
            first.yield_summary["random_budget"],
        )
        self.assertTrue(first.contract["amplification_required"])
        self.assertFalse(first.contract["heldout_family_used_for_proposals"])
        self.assertEqual(len(first.contract["null_families"]), 3)
        directed_ids = {candidate.candidate_id for candidate in first.candidates}
        random_ids = {
            baseline.baseline_id.replace("random_eigenmode", "spectral_mode")
            for baseline in first.random_baseline
        }
        self.assertTrue(directed_ids.isdisjoint(random_ids))
        self.assertIn(
            first.yield_summary["probe_status"],
            {
                "FAIL_NO_HOSTILE_SURVIVORS",
                "OPEN_NO_DIRECTED_YIELD_ADVANTAGE",
                "PASS_BOUNDED_DIRECTED_YIELD_ADVANTAGE",
            },
        )

    def test_hostile_validation_requires_heldout_family(self) -> None:
        operator = path_laplacian(5)
        thresholds = SurvivalThresholds(0.5, 0.4, 0.6, 0.6, 0.35)
        with self.assertRaisesRegex(ValueError, "heldout_operators"):
            validate_expansion_hostile(
                operator,
                np.ones(5),
                [operator],
                [0.0],
                [],
                [],
                thresholds,
            )


if __name__ == "__main__":
    unittest.main()
