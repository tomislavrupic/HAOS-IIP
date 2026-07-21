from __future__ import annotations

import unittest

import numpy as np

from experiments.emergence_ladder.rung3_recovery_trajectory_v2.controls import ADMISSIBLE_CONTROLS, ORACLE_CONTROL, run_condition
from experiments.emergence_ladder.rung3_recovery_trajectory_v2.fixtures import initial_state, make_fixture, perturbation
from experiments.emergence_ladder.rung3_recovery_trajectory_v2.mechanism import RelationalMemory


class RestorativeControlTests(unittest.TestCase):
    def setUp(self) -> None:
        self.fixture = make_fixture(8)
        self.reference = initial_state(self.fixture, 5003)
        self.perturbed = self.reference + perturbation(self.fixture, "relational_twist", 0.5, 5003)
        self.memory = RelationalMemory.encode(self.fixture.incidence, self.reference, 0.08)

    def run_control(self, name: str):
        return run_condition(name, self.fixture, self.reference, self.perturbed, self.memory, 0.7, 0.03, 10, 0.2, 5003)

    def test_every_declared_control_executes(self) -> None:
        for name in ADMISSIBLE_CONTROLS + (ORACLE_CONTROL,):
            trace = self.run_control(name)
            self.assertEqual(len(trace.states), 11)
            self.assertTrue(np.isfinite(trace.states[-1]).all())

    def test_passive_has_zero_corrective_cost(self) -> None:
        self.assertEqual(sum(self.run_control("passive_relaxation").corrections), 0.0)

    def test_signal_blind_spends_nonzero_corrective_energy(self) -> None:
        self.assertGreater(sum(self.run_control("signal_blind").corrections), 0.0)

    def test_oracle_is_not_an_admissible_control(self) -> None:
        self.assertNotIn(ORACLE_CONTROL, ADMISSIBLE_CONTROLS)


if __name__ == "__main__":
    unittest.main()
