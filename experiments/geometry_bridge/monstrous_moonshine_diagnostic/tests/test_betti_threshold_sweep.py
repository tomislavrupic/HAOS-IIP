from __future__ import annotations

import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "threshold_sweep_betti_stability.py"
SPEC = importlib.util.spec_from_file_location("threshold_sweep_betti_stability", MODULE_PATH)
assert SPEC and SPEC.loader
sweep = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = sweep
SPEC.loader.exec_module(sweep)


class BettiThresholdSweepTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = sweep.BettiSweepConfig()
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmpdir.name)
        cls.result = sweep.run_threshold_sweep(cls.config, cls.output_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmpdir.cleanup()

    def test_threshold_grid_and_reference_row(self) -> None:
        with (self.output_dir / "betti_threshold_sweep.csv").open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 21)
        reference = next(row for row in rows if row["threshold"] == "8.0")
        self.assertEqual(reference["betti0"], "4")
        self.assertEqual(reference["edge_count"], "33")
        self.assertEqual(reference["status_region"], "REFERENCE_THRESHOLD")

    def test_stability_bands_are_reported(self) -> None:
        summary = self.result["sweep_summary"]
        self.assertEqual(summary["reference_betti0"], 4)
        self.assertEqual(summary["reference_edge_count"], 33)
        self.assertEqual(summary["betti0_stability_band"], [8.0, 11.5])
        self.assertEqual(summary["edge_exact_band"], [8.0, 8.5])
        self.assertEqual(summary["edge_neighborhood_band"], [8.0, 10.0])
        self.assertIn("BETTI0_ROBUSTNESS_BAND_DETECTED", self.result["labels"])

    def test_null_ensemble_false_positive_rate_is_visible(self) -> None:
        null_summary = self.result["null_summary"]
        self.assertEqual(null_summary["null_seed_count"], 100)
        self.assertEqual(null_summary["false_positive_count"], 11)
        self.assertAlmostEqual(null_summary["false_positive_rate"], 0.11)
        self.assertIn("NULL_FALSE_POSITIVES_DETECTED", self.result["labels"])

    def test_null_csv_contains_100_deterministic_impostors(self) -> None:
        with (self.output_dir / "betti_null_ensemble.csv").open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 100)
        self.assertEqual(rows[0]["seed"], "0")
        self.assertEqual(rows[-1]["seed"], "99")
        self.assertEqual(sum(1 for row in rows if row["false_positive"] == "true"), 11)

    def test_outputs_exist_and_are_reproducible(self) -> None:
        for filename in (
            "betti_threshold_sweep.csv",
            "betti_null_ensemble.csv",
            "betti_threshold_sweep_report.md",
            "betti_threshold_sweep_result.json",
        ):
            self.assertTrue((self.output_dir / filename).exists(), filename)
        with tempfile.TemporaryDirectory() as other:
            rerun = sweep.run_threshold_sweep(self.config, Path(other))
        self.assertEqual(self.result["result_hash"], rerun["result_hash"])

    def test_result_remains_claim_gated(self) -> None:
        self.assertEqual(self.result["status"], "PASS")
        self.assertEqual(self.result["classification"], "TOPOLOGICAL_DIAGNOSTIC_SIDECAR")
        self.assertIn("LEAN_THEOREM_NOT_INCLUDED", self.result["labels"])
        self.assertIn("PHYSICAL_BRIDGE_NOT_ESTABLISHED", self.result["labels"])


if __name__ == "__main__":
    unittest.main()

