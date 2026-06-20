from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_betti_component_count.py"
SPEC = importlib.util.spec_from_file_location("betti_component_count", MODULE_PATH)
assert SPEC and SPEC.loader
betti = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = betti
SPEC.loader.exec_module(betti)


class BettiComponentRecoverabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = betti.BettiComponentConfig()
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmpdir.name)
        cls.result = betti.run_betti(cls.config, cls.output_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmpdir.cleanup()

    def test_graph_spec_declares_claim_ceiling(self) -> None:
        spec = json.loads((betti.SPEC_PATH).read_text(encoding="utf-8"))
        self.assertEqual(spec["claim_ceiling"], "TOPOLOGICAL_DIAGNOSTIC_SIDECAR")
        self.assertEqual(spec["primary_invariant"], "Betti_0 / componentCount")
        self.assertIn("not a Lean-certified theorem inside this repository", spec["non_claims"])

    def test_representatives_cover_supersingular_primes(self) -> None:
        nodes = betti.base_nodes()
        self.assertEqual([node["prime"] for node in nodes], list(betti.SUPERSINGULAR_PRIMES))
        self.assertEqual(len(nodes), 15)
        self.assertEqual(betti.representative_for_prime(2), (1, 1))
        self.assertEqual(betti.representative_for_prime(3), (3, 0))
        self.assertEqual(betti.representative_for_prime(5), (1, 2))

    def test_stable_controls_preserve_betti_and_edges(self) -> None:
        for condition in ("known_positive", "d4_rotate_90", "d4_reflect_real", "isomorphism_relabel"):
            result = self.result["condition_results"][condition]
            self.assertEqual(result["observed_status"], "PASS", condition)
            self.assertLessEqual(result["betti0_distance"], self.config.invariance_max)
            self.assertLessEqual(result["edge_jaccard_distance"], self.config.invariance_max)

    def test_destructive_controls_move_betti_or_edges(self) -> None:
        for condition in ("threshold_mutation_control", "topology_destroyed_control", "support_replacement_control"):
            result = self.result["condition_results"][condition]
            self.assertEqual(result["observed_status"], "PASS", condition)
            moved = (
                result["betti0_distance"] >= self.config.destructive_betti_min
                or result["edge_jaccard_distance"] >= self.config.destructive_edge_min
            )
            self.assertTrue(moved, condition)

    def test_result_remains_claim_gated(self) -> None:
        self.assertEqual(self.result["status"], "PASS")
        self.assertEqual(self.result["classification"], "TOPOLOGICAL_DIAGNOSTIC_SIDECAR")
        self.assertIn("LEAN_THEOREM_NOT_INCLUDED", self.result["labels"])
        self.assertIn("PHYSICAL_BRIDGE_NOT_ESTABLISHED", self.result["labels"])

    def test_outputs_exist_and_are_reproducible(self) -> None:
        for filename in (
            "betti_relation_graph_nodes.csv",
            "betti_relation_graph_edges.csv",
            "betti_control_results.csv",
            "betti_component_count_report.md",
            "betti_component_count_result.json",
        ):
            self.assertTrue((self.output_dir / filename).exists(), filename)
        with tempfile.TemporaryDirectory() as other:
            rerun = betti.run_betti(self.config, Path(other))
        self.assertEqual(self.result["result_hash"], rerun["result_hash"])

    def test_control_csv_contains_all_conditions(self) -> None:
        with (self.output_dir / "betti_control_results.csv").open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual({row["condition"] for row in rows}, set(betti.CONDITIONS))
        self.assertTrue(all(row["observed_status"] == "PASS" for row in rows))


if __name__ == "__main__":
    unittest.main()

