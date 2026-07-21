from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

import numpy as np

from experiments.hbp.integrity_repair_v2.hbp_ir_01 import (
    Dataset,
    evaluate,
    make_dataset,
    run,
    select_model,
    transform_control,
)


class HBPIntegrityRepairTests(unittest.TestCase):
    def test_result_is_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            self.assertEqual(run(Path(first))["result_hash"], run(Path(second))["result_hash"])

    def test_holdout_targets_do_not_change_model_selection(self) -> None:
        dataset = make_dataset()
        selected, scores = select_model(dataset)
        target = dataset.target.copy()
        target[dataset.split == "holdout"] *= -1000.0
        altered = Dataset(dataset.row_id, dataset.split, dataset.features, target)
        selected_altered, scores_altered = select_model(altered)
        self.assertEqual(selected, selected_altered)
        self.assertEqual(scores, scores_altered)

    def test_non_finite_target_fails_closed(self) -> None:
        dataset = make_dataset()
        target = dataset.target.copy()
        target[0] = np.nan
        with self.assertRaisesRegex(ValueError, "non-finite target rejected"):
            evaluate(Dataset(dataset.row_id, dataset.split, dataset.features, target))

    def test_leakage_feature_is_rejected(self) -> None:
        dataset, _ = transform_control(make_dataset(), "target_proxy_leakage")
        with self.assertRaisesRegex(ValueError, "undeclared or target-derived feature"):
            evaluate(dataset)

    def test_destructive_controls_change_declared_fields(self) -> None:
        dataset = make_dataset()
        for name in ("target_label_shuffle", "signal_feature_shuffle", "parameter_matched_null"):
            _, changed = transform_control(dataset, name)
            self.assertTrue(changed)

    def test_stored_predictions_are_the_scored_predictions(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            result = run(output)
            with (output / "hbp_ir_01_predictions.csv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            holdout = [row for row in rows if row["split"] == "holdout"]
            y = np.array([float(row["target"]) for row in holdout])
            prediction = np.array([float(row["prediction"]) for row in holdout])
            observed = float(np.sqrt(np.mean((y - prediction) ** 2)))
            self.assertAlmostEqual(observed, result["holdout_metrics"]["rmse"], places=10)


if __name__ == "__main__":
    unittest.main()
