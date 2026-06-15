from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_b3_semantic_refinement_audit.py"
SPEC = importlib.util.spec_from_file_location("b3_semantic_refinement_audit", MODULE_PATH)
assert SPEC and SPEC.loader
audit = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = audit
SPEC.loader.exec_module(audit)


class B3SemanticRefinementAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = audit.SemanticRefinementConfig(source_limit=99)
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmpdir.name)
        cls.result = audit.run_audit(cls.config, cls.output_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmpdir.cleanup()

    def test_precommitment_freezes_inputs_and_forbids_chsh(self) -> None:
        contract = json.loads((self.output_dir / "b3_2_2_precommitment_contract.json").read_text(encoding="utf-8"))
        self.assertEqual(contract["name"], "B3.2.2_SEMANTIC_REFINEMENT_SPECIFICITY_PRECOMMITMENT")
        self.assertEqual(contract["frozen_inputs"]["J_lambda"], "chain_orientation_operator from B3.2")
        self.assertFalse(contract["chsh_policy"]["computed"])
        self.assertFalse(contract["chsh_policy"]["authorized"])
        self.assertIn("thresholds moved after inspection", contract["hard_invalidation_conditions"])

    def test_chsh_is_not_computed(self) -> None:
        self.assertIsNone(self.result["chsh_score"])
        self.assertFalse(self.result["chsh_scoring_computed"])
        self.assertFalse(self.result["chsh_scoring_authorized"])
        self.assertIn("CHSH_SCORING_NOT_AUTHORIZED", self.result["labels"])

    def test_semantic_gate_fails_despite_label_degradation(self) -> None:
        with (self.output_dir / "b3_2_2_semantic_ordering.csv").open("r", encoding="utf-8", newline="") as handle:
            rows = {row["control_id"]: row for row in csv.DictReader(handle)}
        target = rows["target_geometry"]
        label = rows["label_permuted_settings"]
        self.assertEqual(target["status"], "SETTING_SEMANTICS_NOT_ESTABLISHED")
        self.assertEqual(label["status"], "LABEL_CONTROL_DEGRADED")
        self.assertGreater(float(target["semantic_composite"]), float(label["semantic_composite"]))
        self.assertFalse(self.result["S1_semantic_ordering_passed"])

    def test_refinement_gate_fails_when_control_degradation_is_insufficient(self) -> None:
        with (self.output_dir / "b3_2_2_refinement_persistence.csv").open("r", encoding="utf-8", newline="") as handle:
            rows = {row["hierarchy_label"]: row for row in csv.DictReader(handle)}
        target = rows["frozen_branch"]
        control = rows["periodic_diagonal_augmented_control"]
        self.assertEqual(target["status"], "REFINEMENT_SPECIFICITY_NOT_ESTABLISHED")
        self.assertEqual(control["status"], "REFINEMENT_CONTROL_NOT_DEGRADED")
        self.assertEqual(int(target["degradation_points_vs_control"]), 3)
        self.assertFalse(self.result["S2_refinement_persistence_passed"])

    def test_dual_failure_returns_generic_relational_geometry(self) -> None:
        self.assertIn("SETTING_SEMANTICS_NOT_ESTABLISHED", self.result["labels"])
        self.assertIn("REFINEMENT_SPECIFICITY_NOT_ESTABLISHED", self.result["labels"])
        self.assertIn("GENERIC_RELATIONAL_GEOMETRY", self.result["labels"])
        self.assertNotIn("HAOS_SPECIFIC_RELATIONAL_STRUCTURE_CANDIDATE", self.result["labels"])

    def test_semantic_edges_include_target_and_label_controls(self) -> None:
        with (self.output_dir / "b3_2_2_semantic_edges.csv").open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), len(audit.LEFT_SETTINGS) * len(audit.RIGHT_SETTINGS) * 2)
        self.assertEqual({row["control_id"] for row in rows}, {"target_geometry", "label_permuted_settings"})

    def test_deterministic_reproducibility(self) -> None:
        with tempfile.TemporaryDirectory() as other:
            rerun = audit.run_audit(self.config, Path(other))
        self.assertEqual(self.result["result_hash"], rerun["result_hash"])

    def test_required_outputs_exist(self) -> None:
        for filename in (
            "b3_2_2_precommitment_contract.json",
            "b3_2_2_semantic_ordering.csv",
            "b3_2_2_semantic_edges.csv",
            "b3_2_2_refinement_persistence.csv",
            "b3_2_2_refinement_pairs.csv",
            "b3_2_2_report.md",
            "b3_2_2_result.json",
        ):
            self.assertTrue((self.output_dir / filename).exists(), filename)


if __name__ == "__main__":
    unittest.main()
