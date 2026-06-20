from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_moonshine_betti_bridge.py"
SPEC = importlib.util.spec_from_file_location("moonshine_betti_bridge", MODULE_PATH)
assert SPEC and SPEC.loader
bridge = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = bridge
SPEC.loader.exec_module(bridge)


class MoonshineBettiBridgeBoundaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = bridge.MoonshineBettiBridgeConfig()
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmpdir.name)
        cls.result = bridge.run_bridge(cls.config, cls.output_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmpdir.cleanup()

    def test_manifest_declares_shared_support_not_derivation(self) -> None:
        manifest = json.loads((self.output_dir / "moonshine_betti_bridge_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["shared_support"], "15 supersingular primes")
        self.assertEqual(manifest["moonshine_channel"], "arithmetic support diagnostic")
        self.assertEqual(manifest["betti_channel"], "Gaussian-prime threshold graph topology")
        self.assertEqual(manifest["bridge_type"], "shared-support diagnostic coupling")
        self.assertEqual(manifest["allowed_phrase"], "shared-support perturbation covariance")
        self.assertEqual(manifest["dangerous_phrase"], "Moonshine-Betti mechanism")
        for non_claim in (
            "Moonshine proof",
            "Monster action",
            "Moonshine module construction",
            "physical bridge",
            "continuum limit",
            "quantum claim",
            "gravity claim",
        ):
            self.assertIn(non_claim, manifest["not_claimed"])

    def test_covariance_rows_are_narrow_and_predeclared(self) -> None:
        rows = self.result["covariance_rows"]
        self.assertEqual([row["bridge_case"] for row in rows], ["self_recovery", "shared_support_replacement"])
        self.assertEqual(rows[0]["observed_status"], "PASS")
        self.assertEqual(rows[0]["moonshine_total_distance"], "0")
        self.assertEqual(rows[0]["betti_vector_distance"], "0")
        self.assertEqual(rows[1]["observed_status"], "PASS")
        self.assertGreaterEqual(float(rows[1]["moonshine_total_distance"]), self.config.moonshine_move_min)
        betti_moved = float(rows[1]["betti_vector_distance"]) >= self.config.betti_vector_move_min or float(
            rows[1]["betti_edge_jaccard_distance"]
        ) >= self.config.betti_edge_move_min
        self.assertTrue(betti_moved)

    def test_result_remains_claim_gated(self) -> None:
        self.assertEqual(self.result["status"], "PASS")
        self.assertEqual(self.result["classification"], "SHARED_SUPPORT_DIAGNOSTIC_COUPLING")
        self.assertIn("SHARED_SUPPORT_PERTURBATION_COVARIANCE_REPORTED", self.result["labels"])
        self.assertIn("MOONSHINE_BETTI_MECHANISM_NOT_ESTABLISHED", self.result["labels"])
        self.assertIn("MOONSHINE_PROOF_NOT_ESTABLISHED", self.result["labels"])
        self.assertIn("PHYSICAL_BRIDGE_NOT_ESTABLISHED", self.result["labels"])
        self.assertIn("Shared-support perturbation covariance only", self.result["claim_ceiling"])

    def test_unmatched_controls_are_visible(self) -> None:
        self.assertEqual(
            self.result["unmatched_controls"]["moonshine_only"],
            ["exponent_shuffle_control", "decomposition_break_control", "gaussian_class_swap_control"],
        )
        self.assertEqual(
            self.result["unmatched_controls"]["betti_only"],
            [
                "d4_rotate_90",
                "d4_reflect_real",
                "isomorphism_relabel",
                "threshold_mutation_control",
                "topology_destroyed_control",
            ],
        )

    def test_outputs_exist_and_are_reproducible(self) -> None:
        for filename in (
            "moonshine_betti_bridge_manifest.json",
            "moonshine_betti_bridge_covariance.csv",
            "moonshine_betti_bridge_report.md",
            "moonshine_betti_bridge_result.json",
        ):
            self.assertTrue((self.output_dir / filename).exists(), filename)
        with tempfile.TemporaryDirectory() as other:
            rerun = bridge.run_bridge(self.config, Path(other))
        self.assertEqual(self.result["result_hash"], rerun["result_hash"])

    def test_covariance_csv_contains_required_fields(self) -> None:
        with (self.output_dir / "moonshine_betti_bridge_covariance.csv").open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[1]["moonshine_condition"], "nonsupersingular_prime_control")
        self.assertEqual(rows[1]["betti_condition"], "support_replacement_control")
        self.assertEqual(rows[1]["observed_status"], "PASS")


if __name__ == "__main__":
    unittest.main()
