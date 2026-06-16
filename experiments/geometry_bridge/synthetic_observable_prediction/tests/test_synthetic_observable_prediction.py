from __future__ import annotations

import csv
import importlib.util
import json
import sys
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_synthetic_observable_prediction.py"
SPEC = importlib.util.spec_from_file_location("synthetic_observable_prediction", MODULE_PATH)
assert SPEC and SPEC.loader
observable = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = observable
SPEC.loader.exec_module(observable)


class SyntheticObservablePredictionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = observable.evaluate()

    def test_contract_written(self) -> None:
        contract = json.loads(observable.CONTRACT_PATH.read_text(encoding="utf-8"))
        self.assertEqual(contract["bridge_id"], "GEO-04")
        self.assertIn("synthetic calibration only", contract["claim_boundary"])

    def test_observations_written(self) -> None:
        with observable.OBSERVATIONS_PATH.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertGreater(len(rows), 0)
        self.assertEqual({row["transformation"] for row in rows}, set(observable.TRANSFORMATIONS))

    def test_result_open(self) -> None:
        self.assertIn("OBSERVABLE_PREDICTION_OPEN", self.result["labels"])
        self.assertIn("MIXED_OPEN", self.result["labels"])
        holdout = self.result["holdout_metrics"]
        self.assertGreater(holdout["pairwise_accuracy"], holdout["pairwise_null_accuracy"])
        self.assertEqual(holdout["pairwise_accuracy"], 1.0)


if __name__ == "__main__":
    unittest.main()
