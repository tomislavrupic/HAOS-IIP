from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_synthetic_geometry_recovery.py"
SPEC = importlib.util.spec_from_file_location("synthetic_geometry_recovery", MODULE_PATH)
assert SPEC and SPEC.loader
geometry = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = geometry
SPEC.loader.exec_module(geometry)


class SyntheticGeometryRecoveryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = geometry.GeometryConfig()
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmpdir.name)
        cls.result = geometry.run_benchmark(cls.config, cls.output_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmpdir.cleanup()

    def test_precommitment_exists(self) -> None:
        contract = json.loads((self.output_dir / "precommitment_contract.json").read_text(encoding="utf-8"))
        self.assertEqual(contract["bridge_id"], "GEO-01")
        self.assertEqual(contract["classification"], "OPERATIONAL_MAPPING_PARTIAL")
        self.assertIn("synthetic calibration only", contract["claim_boundary"])
        self.assertEqual(contract["holdout_split"], "spiral")

    def test_outputs_exist(self) -> None:
        for filename in (
            "precommitment_contract.json",
            "geometry_source_manifest.json",
            "geometry_matrix.csv",
            "control_results.csv",
            "geometry_recovery_report.md",
            "geometry_recovery_result.json",
        ):
            self.assertTrue((self.output_dir / filename).exists(), filename)

    def test_holdout_is_blindly_frozen(self) -> None:
        manifest = json.loads((self.output_dir / "geometry_source_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["holdout_families"], ["spiral"])
        self.assertEqual(manifest["development_families"], ["ring"])

    def test_geometry_matrix_has_expected_rows(self) -> None:
        with (self.output_dir / "geometry_matrix.csv").open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        expected_rows = (len(geometry.DEVELOPMENT_FAMILIES) + len(geometry.CALIBRATION_FAMILIES) + len(geometry.HOLDOUT_FAMILIES)) * len(geometry.SEEDS) * (geometry.GeometryConfig().n_nodes * (geometry.GeometryConfig().n_nodes - 1) // 2)
        self.assertEqual(len(rows), expected_rows)
        self.assertEqual({row["split"] for row in rows}, {"all"})

    def test_deterministic_reproducibility(self) -> None:
        self.assertIn("result_hash", self.result)
        with tempfile.TemporaryDirectory() as other:
            rerun = geometry.run_benchmark(self.config, Path(other))
        self.assertEqual(self.result["result_hash"], rerun["result_hash"])

    def test_holdout_result_present(self) -> None:
        self.assertIn("holdout_metrics", self.result)
        self.assertIn("haos", self.result["holdout_metrics"])
        self.assertIn("best_baseline_spearman", self.result["holdout_metrics"])

    def test_control_results_exist(self) -> None:
        with (self.output_dir / "control_results.csv").open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertGreater(len(rows), 0)
        self.assertIn("leakage_positive_control", {row["control_name"] for row in rows})


if __name__ == "__main__":
    unittest.main()
