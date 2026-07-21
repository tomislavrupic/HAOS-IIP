from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from experiments.emergence_ladder.rung3_recovery_trajectory.run_recovery_trajectory import (
    classify,
    orthogonal_perturbation,
    run,
)


class RecoveryTrajectoryTests(unittest.TestCase):
    def test_perturbation_is_orthogonal(self) -> None:
        reference = np.array([1.0, 0.0, 0.0], dtype=complex)
        _, overlap = orthogonal_perturbation(reference, 0.5, 3301)
        self.assertLess(overlap, 1.0e-10)

    def test_failed_gate_is_negative_result(self) -> None:
        self.assertEqual(classify({"a": True, "b": False}), "NEGATIVE_RESULT")
        self.assertEqual(classify({"a": True}, instrument_valid=False), "INSTRUMENT_INVALID")

    def test_result_is_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            self.assertEqual(run(Path(first))["result_hash"], run(Path(second))["result_hash"])

    def test_passive_control_has_zero_gain(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            result = run(Path(directory))
            self.assertAlmostEqual(result["means"]["passive_recovery_gain"], 0.0, places=12)

    def test_all_controls_execute(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            result = run(Path(directory))
            self.assertEqual(
                {row["control"] for row in result["controls"]},
                {"passive_relaxation", "operator_only_filtering", "topology_altered_cycle_phase", "trivial_attractor"},
            )
            self.assertTrue(all(row["mechanism_valid"] for row in result["controls"]))


if __name__ == "__main__":
    unittest.main()
