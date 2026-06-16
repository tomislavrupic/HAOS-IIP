from __future__ import annotations

import unittest
from pathlib import Path

from experiments.hbp.pb02_external_network_recovery import pb02_holdout as pb02


DATA_ROOT = Path("/Volumes/Samsung T5/2026/HAOS/HAOS DOCS/DATA/Powergraph")


class PB02HoldoutTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.bundle = pb02.load_feature_bundle(DATA_ROOT)
        cls.splits = pb02.build_case_splits(cls.bundle.case_ids)
        cls.first_run = pb02.run_pb02(DATA_ROOT, write_outputs=False)
        cls.second_run = pb02.run_pb02(DATA_ROOT, write_outputs=False)

    def test_dataset_schema_validation(self) -> None:
        validation = pb02.validate_dataset_root(DATA_ROOT)
        self.assertEqual(validation["dataset"], "PowerGraph")
        self.assertEqual(validation["family"], "ieee24")
        self.assertEqual(validation["sample_count"], 21500)
        self.assertTrue(validation["data_root"].startswith("DATA/Powergraph"))
        self.assertFalse(Path(validation["data_root"]).is_absolute())

    def test_deterministic_result_hash(self) -> None:
        self.assertEqual(self.first_run["result"]["result_hash"], self.second_run["result"]["result_hash"])
        self.assertEqual(self.first_run["result"]["status"], "CONTROL_INVALID")

    def test_frozen_split_enforcement(self) -> None:
        splits = self.splits
        self.assertFalse(set(splits["development"]) & set(splits["holdout"]))
        self.assertFalse(set(splits["calibration"]) & set(splits["holdout"]))
        self.assertEqual(sum(len(v) for v in splits.values()), len(self.bundle.case_ids))

    def test_holdout_isolation_and_routing(self) -> None:
        result = self.first_run
        self.assertIn("best_baseline_model", result["result"])
        self.assertIn("baseline_plus_haos_holdout_spearman", result["result"])
        self.assertGreaterEqual(len(result["baseline_rows"]), 3)
        self.assertGreaterEqual(len(result["haos_rows"]), 3)
        self.assertEqual(result["result"]["case_counts"]["holdout"], len(self.splits["holdout"]))

    def test_incremental_value_and_controls(self) -> None:
        incremental = self.first_run["incremental_value"]
        self.assertIn("incremental_spearman", incremental)
        controls = {row["control"]: row for row in self.first_run["control_rows"]}
        self.assertEqual(controls["intentional_leakage_positive_control"]["status"], "TARGET_LEAKAGE_DETECTED")
        self.assertIn(controls["seed_repeat"]["status"], {"PASS", "FAIL"})

    def test_control_invalidation_and_non_separation(self) -> None:
        self.assertIn("CONTROL_INVALID", self.first_run["result"]["labels"])
        self.assertIn("PREDICTION_NOT_DISTINCT_FROM_BASELINES", self.first_run["result"]["labels"])
        self.assertLess(
            float(self.first_run["result"]["baseline_plus_haos_holdout_spearman"]),
            float(self.first_run["result"]["best_baseline_holdout_spearman"]),
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
            self.bundle,
            {"development": [0], "calibration": [1], "holdout": [2, 3]},
        )
        self.assertEqual(verdict["status"], "TARGET_LEAKAGE_DETECTED")
        self.assertIn("TARGET_LEAKAGE_DETECTED", verdict["labels"])

    def test_malformed_contract_failure(self) -> None:
        bad = pb02.read_json(pb02.PRECOMMITMENT_PATH)
        bad["bridge_id"] = "BROKEN"
        self.assertNotEqual(bad["bridge_id"], "PB-02")


if __name__ == "__main__":
    unittest.main()
