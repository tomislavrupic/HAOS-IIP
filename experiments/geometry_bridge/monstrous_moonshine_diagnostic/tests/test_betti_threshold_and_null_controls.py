from __future__ import annotations

import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "threshold_sweep_betti_stability.py"
SPEC = importlib.util.spec_from_file_location("threshold_sweep_betti_stability_combined", MODULE_PATH)
assert SPEC and SPEC.loader
combined = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = combined
SPEC.loader.exec_module(combined)


class BettiThresholdAndNullControlsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = combined.BettiSweepConfig()
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmpdir.name)
        cls.result = combined.run_threshold_sweep(cls.config, cls.output_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmpdir.cleanup()

    def test_requested_output_names_exist(self) -> None:
        for filename in (
            "betti_threshold_sweep.csv",
            "betti_threshold_sweep_report.md",
            "betti_null_ensemble_results.csv",
            "betti_null_ensemble_report.md",
        ):
            self.assertTrue((self.output_dir / filename).exists(), filename)

    def test_required_minimum_labels_are_present(self) -> None:
        labels = set(self.result["labels"])
        for label in (
            "BETTI_THRESHOLD_SWEEP_BUILT",
            "BETTI_NULL_ENSEMBLE_BUILT",
            "BETTI_ROBUSTNESS_OPEN",
            "FALSE_POSITIVE_RATE_REPORTED",
            "LEAN_THEOREM_NOT_INCLUDED",
            "PHYSICAL_BRIDGE_NOT_ESTABLISHED",
        ):
            self.assertIn(label, labels)

    def test_threshold_sweep_records_required_fields(self) -> None:
        with (self.output_dir / "betti_threshold_sweep.csv").open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 21)
        self.assertEqual(rows[0]["threshold"], "2.0")
        self.assertEqual(rows[-1]["threshold"], "12.0")
        for field in (
            "threshold",
            "edge_count",
            "betti0",
            "component_size_distribution",
            "largest_component_size",
            "edge_jaccard_distance_to_reference",
            "status_region",
        ):
            self.assertIn(field, rows[0])

    def test_null_ensemble_reports_damage_not_flattery(self) -> None:
        with (self.output_dir / "betti_null_ensemble_results.csv").open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        false_positives = [row for row in rows if row["false_positive"] == "true"]
        self.assertEqual(len(rows), 100)
        self.assertEqual(len(false_positives), 11)
        self.assertAlmostEqual(self.result["null_summary"]["false_positive_rate"], 0.11)
        self.assertIn("NULL_FALSE_POSITIVES_DETECTED", self.result["labels"])


if __name__ == "__main__":
    unittest.main()

