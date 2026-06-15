from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_synthetic_calibration.py"
SPEC = importlib.util.spec_from_file_location("synthetic_calibration", MODULE_PATH)
assert SPEC and SPEC.loader
calibration = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = calibration
SPEC.loader.exec_module(calibration)


class SyntheticCalibrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = calibration.CalibrationConfig()
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmpdir.name)
        cls.result = calibration.run_calibration(cls.config, cls.output_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmpdir.cleanup()

    def test_precommitment_contract_declares_non_claims(self) -> None:
        contract = json.loads((self.output_dir / "precommitment_contract.json").read_text(encoding="utf-8"))
        self.assertEqual(contract["name"], "Synthetic Relational Geometry Calibration v0.1")
        self.assertIn("not a physical Bell experiment", contract["non_claims"])
        self.assertFalse(contract["chsh_policy"]["computed"])
        self.assertFalse(contract["chsh_policy"]["authorized"])

    def test_known_positive_passes(self) -> None:
        condition = self.result["condition_results"]["known_positive"]
        self.assertEqual(condition["observed_status"], "PASS")
        self.assertTrue(condition["semantic"]["semantic_pass"])
        self.assertTrue(condition["refinement"]["refinement_pass"])

    def test_destructive_controls_fail_expected_surfaces(self) -> None:
        semantic = self.result["condition_results"]["semantic_permutation_control"]
        topology = self.result["condition_results"]["topology_destroyed_control"]
        refinement = self.result["condition_results"]["refinement_broken_control"]
        self.assertEqual(semantic["observed_status"], "FAIL_SEMANTIC")
        self.assertIn(topology["observed_status"], {"FAIL_SEMANTIC", "FAIL_BOTH"})
        self.assertIn(refinement["observed_status"], {"FAIL_REFINEMENT", "FAIL_BOTH"})

    def test_ambiguous_case_remains_open(self) -> None:
        ambiguous = self.result["condition_results"]["ambiguous_partial_preservation"]
        self.assertEqual(ambiguous["observed_status"], "OPEN")
        self.assertFalse(ambiguous["semantic"]["semantic_pass"])
        self.assertTrue(ambiguous["refinement"]["refinement_pass"])

    def test_calibration_passes_without_chsh(self) -> None:
        self.assertTrue(self.result["calibration_pass"])
        self.assertIsNone(self.result["chsh_score"])
        self.assertFalse(self.result["chsh_scoring_computed"])
        self.assertIn("CALIBRATION_PASS", self.result["labels"])
        self.assertIn("HAOS_BELL_DERIVATION_NOT_ESTABLISHED", self.result["labels"])

    def test_geometry_matrix_has_all_conditions(self) -> None:
        with (self.output_dir / "calibration_geometry_matrix.csv").open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual({row["condition"] for row in rows}, set(calibration.CONDITIONS))
        expected_rows = len(calibration.CONDITIONS) * len(calibration.SEEDS) * len(calibration.REFINEMENT_LEVELS) * len(calibration.SETTINGS) ** 2
        self.assertEqual(len(rows), expected_rows)

    def test_deterministic_reproducibility(self) -> None:
        with tempfile.TemporaryDirectory() as other:
            rerun = calibration.run_calibration(self.config, Path(other))
        self.assertEqual(self.result["result_hash"], rerun["result_hash"])

    def test_required_outputs_exist(self) -> None:
        for filename in (
            "precommitment_contract.json",
            "synthetic_source_manifest.json",
            "semantic_relation_matrix.csv",
            "calibration_geometry_matrix.csv",
            "calibration_control_results.csv",
            "calibration_refinement_persistence.csv",
            "synthetic_calibration_report.md",
            "synthetic_calibration_result.json",
        ):
            self.assertTrue((self.output_dir / filename).exists(), filename)


if __name__ == "__main__":
    unittest.main()
