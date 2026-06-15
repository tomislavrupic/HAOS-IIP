from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_b3_null_failure_localization.py"
SPEC = importlib.util.spec_from_file_location("b3_null_failure_localization", MODULE_PATH)
assert SPEC and SPEC.loader
localization = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = localization
SPEC.loader.exec_module(localization)


class B3NullFailureLocalizationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = localization.NullLocalizationConfig(source_limit=6)
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmpdir.name)
        cls.result = localization.run_localization(cls.config, cls.output_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmpdir.cleanup()

    def test_precommitment_freezes_b3_2_inputs(self) -> None:
        contract = json.loads((self.output_dir / "b3_2_1_precommitment_contract.json").read_text(encoding="utf-8"))
        self.assertEqual(contract["name"], "B3.2.1_NULL_FAILURE_LOCALIZATION_PRECOMMITMENT")
        self.assertIn("do not alter J_lambda", contract["frozen_inputs"]["prohibited_changes"])
        self.assertFalse(contract["chsh_policy"]["computed"])
        self.assertFalse(contract["chsh_policy"]["authorized"])

    def test_chsh_remains_absent(self) -> None:
        self.assertIsNone(self.result["chsh_score"])
        self.assertFalse(self.result["chsh_scoring_computed"])
        self.assertFalse(self.result["chsh_scoring_authorized"])
        self.assertIn("CHSH_SCORING_NOT_AUTHORIZED", self.result["labels"])

    def test_failing_nulls_are_localized(self) -> None:
        self.assertIn("label_permuted_settings", self.result["failing_nulls"])
        self.assertIn("refinement_broken_control", self.result["failing_nulls"])
        self.assertIn("SETTING_SEMANTICS_NOT_ESTABLISHED", self.result["labels"])
        self.assertIn("REFINEMENT_SPECIFICITY_NOT_ESTABLISHED", self.result["labels"])

    def test_distribution_comparison_contains_required_metrics(self) -> None:
        with (self.output_dir / "b3_2_1_distribution_comparison.csv").open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        target = next(row for row in rows if row["control_id"] == "target_geometry")
        self.assertIn("pair_order_signature", target)
        self.assertIn("distribution_overlap_vs_target", target)
        self.assertEqual(float(target["pair_order_spearman_vs_target"]), 1.0)

    def test_matrix_comparison_distinguishes_random_orthogonal_from_refinement_failure(self) -> None:
        with (self.output_dir / "b3_2_1_matrix_comparison.csv").open("r", encoding="utf-8", newline="") as handle:
            rows = {row["control_id"]: row for row in csv.DictReader(handle)}
        self.assertGreater(float(rows["random_orthogonal_pairing"]["scale_aligned_distance_to_target"]), 0.75)
        self.assertLess(float(rows["refinement_broken_control"]["scale_aligned_distance_to_target"]), 0.75)

    def test_operator_dependence_variants_are_present(self) -> None:
        with (self.output_dir / "b3_2_1_operator_dependence.csv").open("r", encoding="utf-8", newline="") as handle:
            variants = {row["operator_variant"]: row for row in csv.DictReader(handle)}
        self.assertEqual(set(localization.OPERATOR_VARIANTS), set(variants))
        self.assertEqual(variants["target_J"]["status"], "REFERENCE_OPERATOR")
        self.assertIn(variants["degree_matched_J"]["status"], {"OPERATOR_GEOMETRY_PERSISTS", "OPERATOR_GEOMETRY_DEGRADED"})

    def test_null_invariant_ledger_records_preserved_destroyed_fields(self) -> None:
        ledger = json.loads((self.output_dir / "b3_2_1_null_invariant_ledger.json").read_text(encoding="utf-8"))
        self.assertEqual(ledger["frozen_status"]["chsh_scoring"], "NOT_AUTHORIZED")
        self.assertEqual(ledger["controls"]["label_permuted_settings"]["setting_semantics"], "destroyed")
        self.assertEqual(ledger["controls"]["refinement_broken_control"]["refinement_continuity"], "destroyed")

    def test_deterministic_reproducibility(self) -> None:
        with tempfile.TemporaryDirectory() as other:
            rerun = localization.run_localization(self.config, Path(other))
        self.assertEqual(self.result["result_hash"], rerun["result_hash"])

    def test_required_outputs_exist(self) -> None:
        for filename in (
            "b3_2_1_precommitment_contract.json",
            "b3_2_1_distribution_comparison.csv",
            "b3_2_1_matrix_comparison.csv",
            "b3_2_1_topology_dependence.csv",
            "b3_2_1_operator_dependence.csv",
            "b3_2_1_pair_ordering.csv",
            "b3_2_1_null_invariant_ledger.json",
            "b3_2_1_report.md",
            "b3_2_1_result.json",
        ):
            self.assertTrue((self.output_dir / filename).exists(), filename)


if __name__ == "__main__":
    unittest.main()
