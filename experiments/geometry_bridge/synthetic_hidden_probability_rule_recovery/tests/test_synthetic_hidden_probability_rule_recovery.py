from __future__ import annotations

import csv
import importlib.util
import json
import sys
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_synthetic_hidden_probability_rule_recovery.py"
SPEC = importlib.util.spec_from_file_location("synthetic_hidden_probability_rule_recovery", MODULE_PATH)
assert SPEC and SPEC.loader
probability = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = probability
SPEC.loader.exec_module(probability)


class SyntheticHiddenProbabilityRuleRecoveryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = probability.evaluate()

    def test_contract_written(self) -> None:
        contract = json.loads(probability.CONTRACT_PATH.read_text(encoding="utf-8"))
        self.assertEqual(contract["bridge_id"], "GEO-P1-01")
        self.assertIn("synthetic calibration only", contract["claim_boundary"])
        self.assertIn("probability semantics", contract["purpose"].lower())

    def test_observations_and_predictions_written(self) -> None:
        with probability.OBSERVATIONS_PATH.open("r", encoding="utf-8", newline="") as handle:
            obs = list(csv.DictReader(handle))
        with probability.PREDICTIONS_PATH.open("r", encoding="utf-8", newline="") as handle:
            preds = list(csv.DictReader(handle))
        self.assertGreater(len(obs), 0)
        self.assertGreater(len(preds), 0)
        self.assertEqual({row["transformation"] for row in obs}, set(probability.TRANSFORMATIONS))

    def test_result_is_valid_or_open_not_physical_claim(self) -> None:
        self.assertTrue(
            {"PROBABILITY_RULE_VALID", "PROBABILITY_RULE_OPEN"}.intersection(self.result["labels"])
        )

    def test_holdout_metrics_present(self) -> None:
        metrics = self.result["holdout_metrics"]
        self.assertIn("accuracy", metrics)
        self.assertIn("brier_score", metrics)
        self.assertIn("log_loss", metrics)
        self.assertIn("calibration_error", metrics)

    def test_deterministic_reproducibility(self) -> None:
        first = probability.evaluate()
        second = probability.evaluate()
        self.assertEqual(first["result_hash"], second["result_hash"])


if __name__ == "__main__":
    unittest.main()
