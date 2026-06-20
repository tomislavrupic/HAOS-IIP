from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_gaussian_prime_norm_lift.py"
SPEC = importlib.util.spec_from_file_location("gaussian_prime_norm_lift", MODULE_PATH)
assert SPEC and SPEC.loader
bridge = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = bridge
SPEC.loader.exec_module(bridge)


class GaussianPrimeNormLiftTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = bridge.GaussianPrimeNormLiftConfig()
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmpdir.name)
        cls.result = bridge.run_bridge(cls.config, cls.output_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmpdir.cleanup()

    def test_gaussian_prime_classifier_known_cases(self) -> None:
        self.assertEqual(bridge.gaussian_prime_class(1, 1), "ramified_norm2")
        self.assertEqual(bridge.gaussian_prime_class(3, 0), "inert_axis_prime")
        self.assertEqual(bridge.gaussian_prime_class(2, 1), "split_norm_prime")
        self.assertEqual(bridge.gaussian_prime_class(2, 0), "composite_or_nonprime")
        self.assertEqual(bridge.gaussian_prime_class(1, 0), "unit")

    def test_contract_declares_non_physics_ceiling(self) -> None:
        contract = json.loads((self.output_dir / "precommitment_contract.json").read_text(encoding="utf-8"))
        self.assertEqual(contract["claim_ceiling"], "SYNTHETIC_ARITHMETIC_GEOMETRY_CALIBRATION")
        self.assertIn("not a physical bridge", contract["non_claims"])
        self.assertIn("not Monster moonshine", contract["non_claims"])

    def test_rotation_invariance_and_self_recovery(self) -> None:
        known = self.result["condition_results"]["known_positive"]
        rotated = self.result["condition_results"]["rotation_invariance_control"]
        self.assertEqual(known["observed_status"], "PASS")
        self.assertEqual(rotated["observed_status"], "PASS")
        self.assertLessEqual(max(known["distances"].values()), self.config.invariance_max)
        self.assertLessEqual(max(rotated["distances"].values()), self.config.invariance_max)

    def test_destructive_controls_move_intended_components(self) -> None:
        class_shuffle = self.result["condition_results"]["class_shuffle_control"]
        norm_shuffle = self.result["condition_results"]["norm_shuffle_control"]
        topology = self.result["condition_results"]["topology_destroyed_control"]
        self.assertEqual(class_shuffle["observed_status"], "PASS")
        self.assertGreaterEqual(class_shuffle["distances"]["shell"], self.config.class_shuffle_shell_min)
        self.assertEqual(norm_shuffle["observed_status"], "PASS")
        self.assertGreaterEqual(norm_shuffle["distances"]["norm_lift"], self.config.norm_shuffle_lift_min)
        self.assertGreaterEqual(norm_shuffle["distances"]["spectral"], self.config.norm_shuffle_spectral_min)
        self.assertEqual(topology["observed_status"], "PASS")
        self.assertGreaterEqual(topology["distances"]["spectral"], self.config.topology_spectral_min)
        self.assertGreaterEqual(topology["distances"]["cochain"], self.config.topology_cochain_min)

    def test_result_remains_claim_gated(self) -> None:
        self.assertEqual(self.result["status"], "PASS")
        self.assertEqual(self.result["classification"], "SYNTHETIC_ARITHMETIC_GEOMETRY_CALIBRATION")
        self.assertIn("PHYSICAL_BRIDGE_NOT_ESTABLISHED", self.result["labels"])
        self.assertIn("MONSTER_MOONSHINE_NOT_TESTED", self.result["labels"])

    def test_outputs_exist_and_are_reproducible(self) -> None:
        for filename in (
            "precommitment_contract.json",
            "source_manifest.json",
            "gaussian_prime_nodes.csv",
            "component_scores.csv",
            "control_results.csv",
            "gaussian_prime_norm_lift_report.md",
            "gaussian_prime_norm_lift_result.json",
        ):
            self.assertTrue((self.output_dir / filename).exists(), filename)
        with tempfile.TemporaryDirectory() as other:
            rerun = bridge.run_bridge(self.config, Path(other))
        self.assertEqual(self.result["result_hash"], rerun["result_hash"])

    def test_control_csv_contains_all_conditions(self) -> None:
        with (self.output_dir / "control_results.csv").open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual({row["condition"] for row in rows}, set(bridge.CONDITIONS))
        self.assertTrue(all(row["observed_status"] == "PASS" for row in rows))


if __name__ == "__main__":
    unittest.main()

