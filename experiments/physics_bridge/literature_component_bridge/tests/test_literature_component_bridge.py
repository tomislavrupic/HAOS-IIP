from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_literature_component_bridge.py"
SPEC = importlib.util.spec_from_file_location("literature_component_bridge", MODULE_PATH)
assert SPEC and SPEC.loader
bridge = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = bridge
SPEC.loader.exec_module(bridge)


class LiteratureComponentBridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = bridge.BridgeComponentConfig()
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmpdir.name)
        cls.result = bridge.run_bridge(cls.config, cls.output_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmpdir.cleanup()

    def test_contract_declares_claim_ceiling(self) -> None:
        contract = json.loads((self.output_dir / "precommitment_contract.json").read_text(encoding="utf-8"))
        self.assertEqual(contract["claim_ceiling"], "OPERATIONAL_MAPPING")
        self.assertIn("not a physical bridge", contract["non_claims"])
        self.assertTrue(all(status == "HEURISTIC" for status in contract["mapping_status"].values()))

    def test_self_and_label_permutation_do_not_separate(self) -> None:
        known = self.result["condition_results"]["known_positive"]
        label = self.result["condition_results"]["label_permutation_control"]
        self.assertEqual(known["observed_status"], "PASS")
        self.assertEqual(label["observed_status"], "PASS")
        self.assertLessEqual(max(known["distances"].values()), self.config.label_invariance_max)
        self.assertLessEqual(max(label["distances"].values()), self.config.label_invariance_max)

    def test_destructive_controls_move_intended_components(self) -> None:
        weight = self.result["condition_results"]["weight_shuffle_control"]
        topology = self.result["condition_results"]["topology_destroyed_control"]
        hodge = self.result["condition_results"]["hodge_triangle_removed_control"]
        self.assertEqual(weight["observed_status"], "PASS")
        self.assertGreaterEqual(weight["distances"]["transport_kernel"], self.config.transport_distance_min)
        self.assertEqual(topology["observed_status"], "PASS")
        self.assertGreaterEqual(topology["distances"]["spectral"], self.config.spectral_distance_min)
        self.assertGreaterEqual(topology["distances"]["curvature"], self.config.curvature_distance_min)
        self.assertEqual(hodge["observed_status"], "PASS")
        self.assertGreaterEqual(hodge["distances"]["hodge"], self.config.hodge_distance_min)
        self.assertEqual(hodge["distances"]["spectral"], 0.0)
        self.assertEqual(hodge["distances"]["curvature"], 0.0)

    def test_sinkhorn_divergence_removes_self_bias(self) -> None:
        known_transport = self.result["condition_results"]["known_positive"]["components"]["transport_kernel"]
        self.assertGreater(known_transport["sinkhorn_distance"], 0.0)
        self.assertEqual(known_transport["sinkhorn_divergence"], 0.0)
        self.assertEqual(self.result["condition_results"]["known_positive"]["distances"]["transport_kernel"], 0.0)

    def test_result_remains_claim_gated(self) -> None:
        self.assertEqual(self.result["status"], "PASS")
        self.assertEqual(self.result["classification"], "OPERATIONAL_MAPPING")
        self.assertIn("PHYSICAL_BRIDGE_NOT_ESTABLISHED", self.result["labels"])
        self.assertIn("CLAIM_GATED_OPERATIONAL_MAPPING", self.result["labels"])

    def test_outputs_exist_and_are_reproducible(self) -> None:
        for filename in (
            "precommitment_contract.json",
            "source_manifest.json",
            "component_scores.csv",
            "control_results.csv",
            "literature_component_bridge_report.md",
            "literature_component_bridge_result.json",
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
