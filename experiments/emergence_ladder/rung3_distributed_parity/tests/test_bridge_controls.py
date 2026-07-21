from __future__ import annotations

import unittest

import numpy as np

from experiments.emergence_ladder.rung3_distributed_parity.alphabet import RelationalAlphabet
from experiments.emergence_ladder.rung3_distributed_parity.controls import run_condition
from experiments.emergence_ladder.rung3_distributed_parity.parity_architecture import DistributedParityMemory
from experiments.emergence_ladder.rung3_recovery_trajectory_v2.fixtures import initial_state, make_fixture, perturbation


class BridgeControlTests(unittest.TestCase):
    def setUp(self) -> None:
        self.fixture = make_fixture(8)
        self.alphabet = RelationalAlphabet.build(8)
        self.reference = initial_state(self.fixture, 7002)
        self.memory = DistributedParityMemory.encode(self.alphabet.encode(self.reference), self.alphabet)
        self.perturbed = self.reference + perturbation(self.fixture, "sparse_node_corruption", 0.5, 7002)

    def test_target_is_deterministic_and_finite(self) -> None:
        left = run_condition("target", self.reference, self.perturbed, self.alphabet, self.memory, 0.4, 0.02, 12, 77)
        right = run_condition("target", self.reference, self.perturbed, self.alphabet, self.memory, 0.4, 0.02, 12, 77)
        np.testing.assert_allclose(left.trace.states[-1], right.trace.states[-1])
        self.assertTrue(np.isfinite(left.trace.states[-1]).all())

    def test_passive_has_zero_correction_cost(self) -> None:
        result = run_condition("passive_relaxation", self.reference, self.perturbed, self.alphabet, self.memory, 0.4, 0.02, 12, 77)
        self.assertEqual(sum(result.trace.correction_norms), 0.0)

    def test_signal_blind_matches_target_step_energy(self) -> None:
        target = run_condition("target", self.reference, self.perturbed, self.alphabet, self.memory, 0.4, 0.02, 12, 77)
        blind = run_condition("signal_blind", self.reference, self.perturbed, self.alphabet, self.memory, 0.4, 0.02, 12, 77)
        np.testing.assert_allclose(target.trace.correction_norms, blind.trace.correction_norms)
        self.assertFalse(np.allclose(target.trace.states[-1], blind.trace.states[-1]))

    def test_oracle_is_non_admissible(self) -> None:
        result = run_condition("full_checkpoint_oracle_non_admissible", self.reference, self.perturbed, self.alphabet, self.memory, 0.4, 0.02, 12, 77)
        self.assertFalse(result.admissible)


if __name__ == "__main__":
    unittest.main()
