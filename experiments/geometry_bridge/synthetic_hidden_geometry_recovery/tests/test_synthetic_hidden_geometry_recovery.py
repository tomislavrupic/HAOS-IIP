from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_synthetic_hidden_geometry_recovery.py"
SPEC = importlib.util.spec_from_file_location("synthetic_hidden_geometry_recovery", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class SyntheticHiddenGeometryRecoveryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = MODULE.HiddenGeometryConfig()
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmpdir.name)
        cls.result = MODULE.evaluate()

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmpdir.cleanup()

    def test_contract_exists(self) -> None:
        contract = json.loads((MODULE.CONTRACT_PATH).read_text(encoding="utf-8"))
        self.assertEqual(contract["bridge_id"], "GEO-HIDDEN-01")
        self.assertIn("Recover intrinsic distance", contract["purpose"])

    def test_result_has_metrics(self) -> None:
        self.assertIn("holdout_metrics", self.result)
        self.assertIn("distance_spearman", self.result["holdout_metrics"])
        self.assertIn("orientation_accuracy", self.result["holdout_metrics"])
        self.assertIn("transform_accuracy", self.result["holdout_metrics"])
        self.assertIn("fiedler_sign_stability", self.result["holdout_metrics"])
        self.assertIn("cheeger_conductance", self.result["holdout_metrics"])
        self.assertIn("fiedler_sym_accuracy", self.result["holdout_metrics"])
        self.assertIn("fiedler_rw_accuracy", self.result["holdout_metrics"])
        self.assertIn("relation_accuracy", self.result["holdout_metrics"])
        self.assertTrue(self.result["holdout_metrics"]["distance_pass"])
        self.assertTrue(self.result["holdout_metrics"]["orientation_pass"])
        self.assertTrue(self.result["holdout_metrics"]["relation_pass"])
        self.assertFalse(self.result["holdout_metrics"]["transform_pass"])
        self.assertIn("fiedler_transform_pass", self.result["holdout_metrics"])
        self.assertIn("fiedler_sym_pass", self.result["holdout_metrics"])
        self.assertIn("fiedler_rw_pass", self.result["holdout_metrics"])
        self.assertIn("fiedler_sign_stability_pass", self.result["holdout_metrics"])
        self.assertIn("cheeger_pass", self.result["holdout_metrics"])
        self.assertFalse(self.result["holdout_pass"])
        self.assertEqual(
            self.result["labels"],
            ["BENCHMARK_OPEN", "MIXED_OPEN"],
        )

    def test_outputs_written(self) -> None:
        for filename in (
            "precommitment_contract.json",
            "hidden_geometry_source_manifest.json",
            "hidden_geometry_observations.csv",
            "hidden_geometry_predictions.csv",
            "hidden_geometry_report.md",
            "hidden_geometry_result.json",
        ):
            self.assertTrue((MODULE.ROOT / filename).exists(), filename)

    def test_deterministic_reproducibility(self) -> None:
        with tempfile.TemporaryDirectory() as other:
            rerun = MODULE.evaluate()
        self.assertEqual(self.result["result_hash"], rerun["result_hash"])

    def test_control_results_exist(self) -> None:
        self.assertTrue((MODULE.ROOT / "control_results.csv").exists())
        control_text = (MODULE.ROOT / "control_results.csv").read_text(encoding="utf-8")
        self.assertIn("edge_rewiring_control", control_text)
        self.assertIn("bridge_removal_control", control_text)
        self.assertIn("fiedler_variant_gap", (MODULE.ROOT / "hidden_geometry_observations.csv").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
