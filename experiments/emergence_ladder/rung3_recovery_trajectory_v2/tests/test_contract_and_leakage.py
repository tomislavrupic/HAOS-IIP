from __future__ import annotations

import json
import unittest

from experiments.emergence_ladder.rung3_recovery_trajectory_v2.run_experiment import (
    CONTRACT_PATH,
    leakage_guard,
    select_parameters,
)


class ContractAndLeakageTests(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))

    def test_contract_is_frozen_and_partitions_are_disjoint(self) -> None:
        self.assertEqual(self.contract["status"], "FROZEN")
        self.assertTrue(all(leakage_guard(self.contract).values()))

    def test_parameter_selection_accepts_calibration_rows_only(self) -> None:
        import inspect

        self.assertEqual(list(inspect.signature(select_parameters).parameters), ["rows_by_parameter"])

    def test_final_partition_is_not_in_calibration_or_validation(self) -> None:
        final = set(self.contract["final_evaluation_partition"]["seeds"])
        calibration = set(self.contract["calibration_partition"]["seeds"])
        validation = set(self.contract["validation_partition"]["seeds"])
        self.assertFalse(final & calibration)
        self.assertFalse(final & validation)

    def test_full_state_checkpoint_is_explicitly_prohibited(self) -> None:
        prohibited = " ".join(self.contract["prohibited_information"]).lower()
        self.assertIn("full pre-intervention node-state checkpoint", prohibited)
        self.assertIn("unperturbed future trajectory", prohibited)


if __name__ == "__main__":
    unittest.main()
