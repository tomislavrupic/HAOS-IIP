from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_synthetic_transformation_semantics.py"
SPEC = importlib.util.spec_from_file_location("synthetic_transformation_semantics", MODULE_PATH)
assert SPEC and SPEC.loader
semantics = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = semantics
SPEC.loader.exec_module(semantics)


class SyntheticTransformationSemanticsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmpdir.name)
        # The script writes into its own directory; run by import-side evaluation is avoided.
        cls.result = semantics.evaluate()

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmpdir.cleanup()

    def test_contract_written(self) -> None:
        self.assertTrue((semantics.CONTRACT_PATH).exists())
        contract = json.loads(semantics.CONTRACT_PATH.read_text(encoding="utf-8"))
        self.assertEqual(contract["bridge_id"], "GEO-02")
        self.assertIn("synthetic calibration only", contract["claim_boundary"])

    def test_results_written(self) -> None:
        with semantics.RESULTS_PATH.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertGreater(len(rows), 0)
        self.assertEqual({row["transformation"] for row in rows}, set(semantics.TRANSFORMATIONS))

    def test_holdout_and_controls_present(self) -> None:
        self.assertIn("TRANSFORMATION_SEMANTICS_OPEN", self.result["labels"])
        self.assertIn("MIXED_OPEN", self.result["labels"])
        self.assertIn("leakage_positive_control", self.result["control_summary"])


if __name__ == "__main__":
    unittest.main()

