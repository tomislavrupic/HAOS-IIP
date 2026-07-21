from __future__ import annotations

import unittest

import numpy as np

from experiments.emergence_ladder.rung3_recovery_trajectory_v2.fixtures import initial_state, make_fixture, perturbation
from experiments.emergence_ladder.rung3_recovery_trajectory_v2.mechanism import (
    RelationalMemory,
    assert_memory_has_no_reference_fields,
    compatibility,
    simulate_feedback,
)


class RestorativeMechanismTests(unittest.TestCase):
    def setUp(self) -> None:
        self.fixture = make_fixture(8)
        self.reference = initial_state(self.fixture, 5001)
        self.memory = RelationalMemory.encode(self.fixture.incidence, self.reference, 0.08)

    def test_memory_contains_only_one_bit_relational_identity(self) -> None:
        self.assertTrue(set(np.unique(self.memory.orientation)) <= {-1.0, 0.0, 1.0})
        self.assertFalse(hasattr(self.memory, "reference"))
        self.assertFalse(hasattr(self.memory, "checkpoint"))
        assert_memory_has_no_reference_fields()

    def test_correction_preserves_global_mean(self) -> None:
        perturbed = self.reference + perturbation(self.fixture, "localized_block_shock", 0.5, 5001)
        trace = simulate_feedback(perturbed, self.memory, 0.7, 0.03, 10, 0.2, 5001)
        self.assertAlmostEqual(float(np.mean(trace.states[0])), float(np.mean(trace.states[-1])), places=12)

    def test_constraint_error_decreases_on_development_fixture(self) -> None:
        perturbed = self.reference + perturbation(self.fixture, "sparse_node_corruption", 0.6, 5002)
        initial_error = compatibility(self.memory, perturbed, 0.03)[1]
        trace = simulate_feedback(perturbed, self.memory, 0.7, 0.03, 40, 0.2, 5002)
        final_error = compatibility(self.memory, trace.states[-1], 0.03)[1]
        self.assertLess(float(np.mean(final_error)), float(np.mean(initial_error)))

    def test_no_future_or_score_enters_simulation(self) -> None:
        import inspect

        parameters = set(inspect.signature(simulate_feedback).parameters)
        for forbidden in ("reference", "future", "target", "score", "classification"):
            self.assertFalse(any(forbidden in name for name in parameters))


if __name__ == "__main__":
    unittest.main()
