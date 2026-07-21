from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from experiments.emergence_ladder.rung4_operational_closure.run_operational_closure import (
    CONTRACT,
    REPO,
    bootstrap_ci,
    classify,
    load_runs,
    run,
    shuffled_training_runs,
)


class OperationalClosureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        import json

        cls.contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
        cls.runs = load_runs(REPO / cls.contract["source"]["artifact"])

    def test_source_has_expected_independent_run_keys(self) -> None:
        self.assertEqual(len(self.runs), 216)
        self.assertEqual({run.n_side for run in self.runs}, {48, 60, 72, 84})
        self.assertEqual(len({run.run_id for run in self.runs}), len(self.runs))

    def test_holdout_level_is_not_used_for_construction(self) -> None:
        construction = set(self.contract["source"]["construction_levels"])
        self.assertNotIn(self.contract["source"]["holdout_level"], construction)

    def test_transition_shuffle_changes_training_chains(self) -> None:
        construction = [run for run in self.runs if run.n_side in {48, 60}]
        shuffled = shuffled_training_runs(construction, 4402)
        self.assertTrue(any(before.chain != after.chain for before, after in zip(construction, shuffled)))

    def test_bootstrap_is_deterministic(self) -> None:
        values = [-0.1, 0.0, 0.2, 0.3]
        self.assertEqual(bootstrap_ci(values, 4401, 100), bootstrap_ci(values, 4401, 100))

    def test_failed_gate_returns_negative_result(self) -> None:
        self.assertEqual(classify({"a": True, "b": False}), "NEGATIVE_RESULT")
        self.assertEqual(classify({"a": True}, instrument_valid=False), "INSTRUMENT_INVALID")

    def test_result_and_predictions_are_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            first_path = Path(first)
            second_path = Path(second)
            first_result = run(first_path)
            second_result = run(second_path)
            self.assertEqual(first_result["result_hash"], second_result["result_hash"])
            self.assertEqual(
                (first_path / "operational_closure_predictions.csv").read_bytes(),
                (second_path / "operational_closure_predictions.csv").read_bytes(),
            )

    def test_prediction_rows_do_not_include_future_context(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            run(output)
            with (output / "operational_closure_predictions.csv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertNotIn("future_event", rows[0])
            self.assertIn("current_event", rows[0])
            self.assertIn("target_event", rows[0])


if __name__ == "__main__":
    unittest.main()
