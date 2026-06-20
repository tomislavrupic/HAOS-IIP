from __future__ import annotations

import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_betti_vector.py"
SPEC = importlib.util.spec_from_file_location("betti_vector", MODULE_PATH)
assert SPEC and SPEC.loader
betti_vector = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = betti_vector
SPEC.loader.exec_module(betti_vector)


class BettiVectorRecoverabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = betti_vector.BettiVectorConfig()
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmpdir.name)
        cls.result = betti_vector.run_betti_vector(cls.config, cls.output_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmpdir.cleanup()

    def test_reference_betti_vector_matches_declared_graph(self) -> None:
        self.assertEqual(
            self.result["reference_signature"],
            {"Betti_0": 4, "Betti_1": 22, "nodes": 15, "edges": 33},
        )
        self.assertEqual(self.result["reference_component_size_distribution"], [12, 1, 1, 1])

    def test_betti1_uses_graph_cycle_count_formula(self) -> None:
        reference = self.result["reference_signature"]
        self.assertEqual(reference["Betti_1"], reference["edges"] - reference["nodes"] + reference["Betti_0"])
        for condition_result in self.result["condition_results"].values():
            signature = condition_result["signature"]
            self.assertEqual(signature["Betti_1"], signature["edges"] - signature["nodes"] + signature["Betti_0"])

    def test_stable_controls_preserve_betti_vector_and_edges(self) -> None:
        for condition in ("known_positive", "d4_rotate_90", "d4_reflect_real", "isomorphism_relabel"):
            result = self.result["condition_results"][condition]
            self.assertEqual(result["observed_status"], "PASS", condition)
            self.assertLessEqual(result["betti_vector_distance"], self.config.invariance_max)
            self.assertLessEqual(result["edge_jaccard_distance"], self.config.invariance_max)
            self.assertEqual(result["signature"]["Betti_0"], 4)
            self.assertEqual(result["signature"]["Betti_1"], 22)

    def test_destructive_controls_move_betti_vector_or_edges(self) -> None:
        for condition in ("threshold_mutation_control", "topology_destroyed_control", "support_replacement_control"):
            result = self.result["condition_results"][condition]
            self.assertEqual(result["observed_status"], "PASS", condition)
            moved = (
                result["betti_vector_distance"] >= self.config.destructive_betti_vector_min
                or result["edge_jaccard_distance"] >= self.config.destructive_edge_min
            )
            self.assertTrue(moved, condition)

    def test_outputs_remain_claim_gated(self) -> None:
        self.assertEqual(self.result["status"], "PASS")
        self.assertEqual(self.result["classification"], "TOPOLOGICAL_DIAGNOSTIC_SIDECAR")
        self.assertIn("BETTI_VECTOR_BUILT", self.result["labels"])
        self.assertIn("BETTI_1_CYCLE_COUNT_AVAILABLE", self.result["labels"])
        self.assertIn("LEAN_THEOREM_NOT_INCLUDED", self.result["labels"])
        self.assertIn("PHYSICAL_BRIDGE_NOT_ESTABLISHED", self.result["labels"])
        self.assertIn("CLAIM_GATED_GRAPH_TOPOLOGY", self.result["labels"])
        self.assertIn("Graph-topology diagnostic only", self.result["claim_ceiling"])

    def test_outputs_exist_and_are_reproducible(self) -> None:
        for filename in (
            "betti_vector_control_results.csv",
            "betti_vector_report.md",
            "betti_vector_result.json",
        ):
            self.assertTrue((self.output_dir / filename).exists(), filename)
        with tempfile.TemporaryDirectory() as other:
            rerun = betti_vector.run_betti_vector(self.config, Path(other))
        self.assertEqual(self.result["result_hash"], rerun["result_hash"])

    def test_control_csv_contains_all_conditions(self) -> None:
        with (self.output_dir / "betti_vector_control_results.csv").open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual({row["condition"] for row in rows}, set(betti_vector.CONDITIONS))
        self.assertTrue(all(row["observed_status"] == "PASS" for row in rows))
        self.assertEqual(rows[0]["betti0"], "4")
        self.assertEqual(rows[0]["betti1"], "22")


if __name__ == "__main__":
    unittest.main()
