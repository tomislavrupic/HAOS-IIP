from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_b3_geometry_audit.py"
SPEC = importlib.util.spec_from_file_location("b3_geometry_audit", MODULE_PATH)
assert SPEC and SPEC.loader
geometry = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = geometry
SPEC.loader.exec_module(geometry)


class B3GeometryAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = geometry.GeometryConfig(source_limit=6)
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmpdir.name)
        cls.result = geometry.run_audit(cls.config, cls.output_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmpdir.cleanup()

    def test_precommitment_contract_is_explicit(self) -> None:
        contract = json.loads((self.output_dir / "b3_2_precommitment_contract.json").read_text(encoding="utf-8"))
        self.assertEqual(contract["name"], "B3.2_PRECOMMITMENT_CONTRACT")
        self.assertIn("relational_invariant", contract)
        self.assertIn("pairing_operator", contract)
        self.assertIn("null_operator_contracts", contract)
        self.assertFalse(contract["target_leakage_guard"]["detected"])

    def test_no_chsh_score_is_computed(self) -> None:
        self.assertIsNone(self.result["chsh_score"])
        self.assertFalse(self.result["chsh_scoring_computed"])
        self.assertIn("CHSH_SCORING_NOT_AUTHORIZED", self.result["gate_results"]["labels"])

    def test_geometry_detects_sensitivity_and_sign_change(self) -> None:
        gates = self.result["gate_results"]["gates"]
        self.assertTrue(gates["G1_RELATIONAL_SENSITIVITY"]["passed"])
        self.assertTrue(gates["G2_SIGN_STRUCTURE"]["passed"])
        stats = gates["G2_SIGN_STRUCTURE"]["stats"]
        self.assertGreater(stats["positive_count"], 0)
        self.assertGreater(stats["negative_count"], 0)

    def test_covariance_diagnostics_pass(self) -> None:
        gates = self.result["gate_results"]["gates"]
        self.assertTrue(gates["G3_COVARIANCE"]["passed"])
        self.assertEqual(gates["G3_COVARIANCE"]["max_error"], 0.0)

    def test_holdout_transfer_is_measured_before_scoreboard(self) -> None:
        gates = self.result["gate_results"]["gates"]
        self.assertTrue(gates["G4_HOLDOUT_TRANSFER"]["passed"])
        self.assertFalse(gates["G6_CHSH_SCOREBOARD"]["computed"])

    def test_null_failure_blocks_chsh_authorization(self) -> None:
        gates = self.result["gate_results"]["gates"]
        self.assertFalse(gates["G5_NULL_REJECTION"]["passed"])
        self.assertFalse(gates["G6_CHSH_SCOREBOARD"]["passed"])
        self.assertFalse(self.result["gate_results"]["chsh_scoring_authorized"])

    def test_required_null_controls_are_present(self) -> None:
        with (self.output_dir / "b3_2_control_results.csv").open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        control_ids = {row["control_id"] for row in rows}
        self.assertEqual(set(geometry.CONTROL_CONTRACTS), control_ids)
        failing = {row["control_id"] for row in rows if row["status"] == "NULL_REJECTION_FAIL"}
        self.assertTrue(failing)

    def test_geometry_leakage_guard(self) -> None:
        detected, hits = geometry.detect_geometry_leakage(["value = -cos(theta_a - theta_b)\n"])
        self.assertTrue(detected)
        self.assertTrue(hits)
        clean, clean_hits = geometry.detect_geometry_leakage(geometry.generator_sources())
        self.assertFalse(clean)
        self.assertEqual(clean_hits, [])

    def test_deterministic_reproducibility(self) -> None:
        with tempfile.TemporaryDirectory() as other:
            rerun = geometry.run_audit(self.config, Path(other))
        self.assertEqual(self.result["result_hash"], rerun["result_hash"])

    def test_required_outputs_exist(self) -> None:
        for filename in (
            "b3_2_precommitment_contract.json",
            "b3_2_geometry_matrix.csv",
            "b3_2_covariance_diagnostics.csv",
            "b3_2_control_results.csv",
            "b3_2_geometry_report.md",
            "b3_2_result.json",
        ):
            self.assertTrue((self.output_dir / filename).exists(), filename)


if __name__ == "__main__":
    unittest.main()
