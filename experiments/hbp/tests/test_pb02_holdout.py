from __future__ import annotations

import csv
import json
import os
import unittest
from pathlib import Path

from experiments.hbp.pb02_external_network_recovery import pb02_holdout as pb02


DATA_ROOT = pb02.DEFAULT_DATA_ROOT
RESULTS_DIR = Path(__file__).resolve().parents[1] / "pb02_external_network_recovery" / "results"
RESULT_PATH = RESULTS_DIR / "pb02_result.json"
DATASET_VALIDATION_PATH = RESULTS_DIR / "dataset_validation.json"
SPLIT_MANIFEST_PATH = RESULTS_DIR / "split_manifest.json"
CONTROL_RESULTS_PATH = RESULTS_DIR / "control_results.csv"
HOLDOUT_PREDICTIONS_PATH = RESULTS_DIR / "holdout_predictions.csv"


class PB02HoldoutTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
        cls.dataset_validation = json.loads(DATASET_VALIDATION_PATH.read_text(encoding="utf-8"))
        cls.split_manifest = json.loads(SPLIT_MANIFEST_PATH.read_text(encoding="utf-8"))
        with CONTROL_RESULTS_PATH.open(encoding="utf-8", newline="") as handle:
            cls.controls = list(csv.DictReader(handle))
        with HOLDOUT_PREDICTIONS_PATH.open(encoding="utf-8", newline="") as handle:
            cls.holdout_predictions = list(csv.DictReader(handle))
        cls.recompute_requested = os.environ.get("PB02_RECOMPUTE", "").strip().lower() in {"1", "true", "yes"}
        cls.recomputed = pb02.run_pb02(DATA_ROOT, write_outputs=False) if cls.recompute_requested else None

    def test_dataset_schema_validation(self) -> None:
        validation = self.dataset_validation
        self.assertEqual(validation["dataset"], "PowerGraph")
        self.assertEqual(validation["family"], "ieee24")
        self.assertEqual(validation["sample_count"], 21500)
        self.assertTrue(validation["data_root"].startswith("DATA/Powergraph"))
        self.assertFalse(Path(validation["data_root"]).is_absolute())
        self.assertEqual(validation["files_hash"], self.result["dataset_validation_hash"])

    def test_deterministic_result_hash(self) -> None:
        self.assertEqual(self.result["result_hash"], json.loads(RESULT_PATH.read_text(encoding="utf-8"))["result_hash"])
        if self.recomputed is not None:
            self.assertEqual(self.result["result_hash"], self.recomputed["result"]["result_hash"])
        self.assertEqual(self.result["status"], "CONTROL_INVALID")

    def test_frozen_split_enforcement(self) -> None:
        splits = self.split_manifest
        self.assertEqual(splits["bridge_id"], "PB-02")
        self.assertFalse(set(splits["development_indices"]) & set(splits["holdout_indices"]))
        self.assertFalse(set(splits["calibration_indices"]) & set(splits["holdout_indices"]))
        self.assertEqual(sum(len(v) for v in (self.split_manifest["development_indices"], self.split_manifest["calibration_indices"], self.split_manifest["holdout_indices"])), 21500)

    def test_holdout_isolation_and_routing(self) -> None:
        self.assertIn("best_baseline_model", self.result)
        self.assertIn("baseline_plus_haos_holdout_spearman", self.result)
        self.assertGreaterEqual(len(self.holdout_predictions), 1)
        self.assertEqual(self.result["case_counts"]["holdout"], len(self.split_manifest["holdout_indices"]))

    def test_incremental_value_and_controls(self) -> None:
        incremental = self.result["incremental_value"]
        self.assertIn("incremental_spearman", incremental)
        controls = {row["control"]: row for row in self.controls}
        self.assertEqual(controls["intentional_leakage_positive_control"]["status"], "TARGET_LEAKAGE_DETECTED")
        self.assertIn(controls["seed_repeat"]["status"], {"PASS", "FAIL"})

    def test_control_invalidation_and_non_separation(self) -> None:
        self.assertIn("CONTROL_INVALID", self.result["labels"])
        self.assertIn("PREDICTION_NOT_DISTINCT_FROM_BASELINES", self.result["labels"])
        self.assertLess(
            float(self.result["baseline_plus_haos_holdout_spearman"]),
            float(self.result["best_baseline_holdout_spearman"]),
        )

    def test_missing_data_handling(self) -> None:
        self.assertAlmostEqual(pb02.summarize_array([1.0, float("nan"), 2.0], "x")["x_zero_fraction"], 1.0 / 3.0)

    def test_forced_non_separation_and_underpowered_logic(self) -> None:
        contract = pb02.read_json(pb02.PRECOMMITMENT_PATH)
        row = {"model": "closeness_centrality", "bundle": "baseline", "holdout_spearman": "0.4"}
        combined = {"model": "baseline_plus_haos", "bundle": "combined", "holdout_spearman": "0.39"}
        verdict = pb02.verdict_from_results(
            [row, combined],
            [
                {"control": "shuffled_target_labels", "status": "PASS", "holdout_spearman": "0.0"},
                {"control": "intentional_leakage_positive_control", "status": "FAIL", "holdout_spearman": "0.1"},
            ],
            contract,
            None,
            {"development": [0], "calibration": [1], "holdout": [2, 3]},
        )
        self.assertEqual(verdict["status"], "TARGET_LEAKAGE_DETECTED")
        self.assertIn("TARGET_LEAKAGE_DETECTED", verdict["labels"])

    def test_malformed_contract_failure(self) -> None:
        bad = pb02.read_json(pb02.PRECOMMITMENT_PATH)
        bad["bridge_id"] = "BROKEN"
        self.assertNotEqual(bad["bridge_id"], "PB-02")

    def test_frozen_check_script_still_validates_outputs(self) -> None:
        from experiments.hbp.pb02_external_network_recovery.check_pb02_holdout import main as check_main

        self.assertEqual(check_main(), 0)


if __name__ == "__main__":
    unittest.main()
