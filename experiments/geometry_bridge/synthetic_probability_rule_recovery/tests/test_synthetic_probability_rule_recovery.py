from __future__ import annotations

import csv
import importlib.util
import json
import sys
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_synthetic_probability_rule_recovery.py"
SPEC = importlib.util.spec_from_file_location("synthetic_probability_rule_recovery", MODULE_PATH)
assert SPEC and SPEC.loader
probability = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = probability
SPEC.loader.exec_module(probability)


class SyntheticProbabilityRuleRecoveryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = probability.evaluate()

    def test_contract_written(self) -> None:
        contract = json.loads(probability.CONTRACT_PATH.read_text(encoding="utf-8"))
        self.assertEqual(contract["bridge_id"], "GEO-03")
        self.assertIn("synthetic calibration only", contract["claim_boundary"])

    def test_observations_and_predictions_written(self) -> None:
        with probability.OBSERVATIONS_PATH.open("r", encoding="utf-8", newline="") as handle:
            obs = list(csv.DictReader(handle))
        with probability.PREDICTIONS_PATH.open("r", encoding="utf-8", newline="") as handle:
            preds = list(csv.DictReader(handle))
        self.assertGreater(len(obs), 0)
        self.assertGreater(len(preds), 0)
        self.assertEqual({row["transformation"] for row in obs}, set(probability.TRANSFORMATIONS))

    def test_result_is_open_not_physical_claim(self) -> None:
        self.assertIn("PROBABILITY_RULE_OPEN", self.result["labels"])
        self.assertIn("MIXED_OPEN", self.result["labels"])


if __name__ == "__main__":
    unittest.main()

