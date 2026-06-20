from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_formal_lean_targets.py"
SPEC = importlib.util.spec_from_file_location("formal_lean_targets", MODULE_PATH)
assert SPEC and SPEC.loader
formal_targets = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = formal_targets
SPEC.loader.exec_module(formal_targets)


class FormalLeanTargetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = formal_targets.FormalLeanTargetConfig()
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmpdir.name)
        cls.result = formal_targets.run_formal_targets(cls.config, cls.output_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmpdir.cleanup()

    def test_manifest_is_open_target_scaffold_not_proof(self) -> None:
        manifest = json.loads((self.output_dir / "formal_lean_targets_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["status"], "OPEN")
        self.assertEqual(manifest["classification"], "FORMAL_TARGET_SCAFFOLD_ONLY")
        self.assertFalse(manifest["config"]["lean_project_present"])
        self.assertFalse(manifest["config"]["lean_check_run"])
        self.assertIn("LEAN_THEOREM_NOT_INCLUDED", manifest["labels"])
        self.assertIn("LEAN_CHECK_NOT_RUN", manifest["labels"])
        self.assertIn("not Lean-certified", manifest["non_claims"])

    def test_required_definitions_are_declared_missing(self) -> None:
        names = {definition["name"]: definition for definition in self.result["required_definitions"]}
        for name in (
            "GaussianPrimeNode",
            "RelationGraph",
            "D4Action",
            "edgeRelation",
            "componentCount",
            "bettiOne",
        ):
            self.assertIn(name, names)
            self.assertEqual(names[name]["local_lean_status"], "MISSING")

    def test_ladder_order_starts_with_graph_isomorphism_before_d4(self) -> None:
        ladder = self.result["theorem_ladder"]
        self.assertEqual([item["id"] for item in ladder], ["L1", "L2", "L3", "L4", "L5", "L6"])
        self.assertEqual(ladder[1]["name"], "graph_iso_preserves_component_count")
        self.assertEqual(ladder[3]["name"], "d4_preserves_component_count")
        self.assertTrue(all(item["status"] == "TARGET_ONLY_NOT_LEAN_CHECKED" for item in ladder))

    def test_executable_evidence_is_betti_vector_only(self) -> None:
        evidence = self.result["executable_evidence"]
        self.assertEqual(evidence["betti_vector_result_hash"], "betti_vector_1d5373fa671fd513fe2d5ac0")
        self.assertEqual(evidence["reference_signature"], {"Betti_0": 4, "Betti_1": 22, "edges": 33, "nodes": 15})
        self.assertTrue(all(status == "PASS" for status in evidence["stable_condition_status"].values()))
        self.assertTrue(all(distance == 0 for distance in evidence["stable_condition_betti_vector_distance"].values()))

    def test_target_lean_file_contains_no_fake_proof_tokens(self) -> None:
        target = (self.output_dir / "lean_graph_invariance_targets.lean").read_text(encoding="utf-8")
        self.assertIn("TARGET_ONLY_NOT_LEAN_CHECKED", target)
        self.assertIn("theorem graph_iso_preserves_component_count", target)
        self.assertNotIn("axiom ", target)
        self.assertNotIn("sorry", target)
        self.assertNotIn("by", target)

    def test_outputs_exist_and_are_reproducible(self) -> None:
        for filename in (
            "formal_lean_targets_manifest.json",
            "formal_lean_targets_report.md",
            "lean_graph_invariance_targets.lean",
        ):
            self.assertTrue((self.output_dir / filename).exists(), filename)
        with tempfile.TemporaryDirectory() as other:
            rerun = formal_targets.run_formal_targets(self.config, Path(other))
        self.assertEqual(self.result["result_hash"], rerun["result_hash"])


if __name__ == "__main__":
    unittest.main()
