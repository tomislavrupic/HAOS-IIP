from __future__ import annotations

import importlib.util
import json
import sys
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_synthetic_hidden_transformation_recovery.py"
SPEC = importlib.util.spec_from_file_location("synthetic_hidden_transformation_recovery", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class SyntheticHiddenTransformationRecoveryTests(unittest.TestCase):
    def test_contract_exists(self) -> None:
        contract = json.loads(MODULE.CONTRACT_PATH.read_text(encoding="utf-8"))
        self.assertEqual(contract["bridge_id"], "GEO-T1-01")
        self.assertIn("Recover identity", contract["purpose"])

    def test_result_has_metrics(self) -> None:
        result = MODULE.evaluate()
        self.assertIn("holdout_metrics", result)
        self.assertTrue(result["holdout_metrics"]["identity_pass"])
        self.assertTrue(result["holdout_metrics"]["inverse_pass"])
        self.assertTrue(result["holdout_metrics"]["composition_pass"])
        self.assertTrue(result["holdout_metrics"]["relation_pass"])
        self.assertTrue(result["holdout_pass"])
        self.assertEqual(result["labels"], ["BENCHMARK_VALID"])

    def test_outputs_written(self) -> None:
        MODULE.evaluate()
        for filename in (
            "precommitment_contract.json",
            "hidden_transformation_source_manifest.json",
            "hidden_transformation_observations.csv",
            "hidden_transformation_predictions.csv",
            "hidden_transformation_report.md",
            "hidden_transformation_result.json",
        ):
            self.assertTrue((MODULE.ROOT / filename).exists(), filename)

    def test_deterministic_reproducibility(self) -> None:
        first = MODULE.evaluate()
        second = MODULE.evaluate()
        self.assertEqual(first["result_hash"], second["result_hash"])

    def test_control_results_exist(self) -> None:
        MODULE.evaluate()
        self.assertTrue((MODULE.ROOT / "control_results.csv").exists())


if __name__ == "__main__":
    unittest.main()
