from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DASHBOARD_PATH = ROOT / "run_geometry_bridge_recoverability_dashboard.py"
SPEC = importlib.util.spec_from_file_location("geometry_bridge_dashboard", DASHBOARD_PATH)
assert SPEC and SPEC.loader
dashboard = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = dashboard
SPEC.loader.exec_module(dashboard)


class FailureConditionsAndDashboardTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = dashboard.GeometryBridgeDashboardConfig()
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmpdir.name)
        cls.result = dashboard.run_dashboard(cls.config, cls.output_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmpdir.cleanup()

    def test_failure_ledger_contains_required_failure_modes(self) -> None:
        ledger = (ROOT / "failure_conditions.md").read_text(encoding="utf-8")
        for phrase in (
            "Threshold Failure",
            "Null-Impostor Failure",
            "D4 Symmetry Failure",
            "Destructive-Control Failure",
            "Moonshine-Betti Decoupling Failure",
            "Formalization Failure",
            "Source-Provenance Failure",
            "Physics Claim Failure",
        ):
            self.assertIn(phrase, ledger)
        self.assertIn("Moonshine-Betti mechanism", ledger)
        self.assertIn("physical bridge established", ledger)

    def test_dashboard_reports_pass_with_fragility_not_global_pass(self) -> None:
        self.assertEqual(self.result["status"], "PASS_WITH_FRAGILITY")
        self.assertEqual(self.result["classification"], "COMPOSITE_DIAGNOSTIC_DASHBOARD_ONLY")
        self.assertIn("PASS_WITH_FRAGILITY", self.result["labels"])
        self.assertIn("GLOBAL_CLAIM_OPEN", self.result["labels"])
        self.assertIn("PHYSICAL_BRIDGE_NOT_ESTABLISHED", self.result["labels"])

    def test_dashboard_channels_include_open_and_partial_items(self) -> None:
        channels = {item["name"]: item for item in self.result["channels"]}
        self.assertEqual(channels["Gaussian-prime norm-lift support"]["status"], "PARTIAL")
        self.assertEqual(channels["Null ensemble rarity"]["status"], "OPEN")
        self.assertEqual(channels["Formal Lean targets"]["status"], "OPEN")
        self.assertEqual(channels["Claim-boundary checks"]["status"], "PASS")
        self.assertEqual(channels["Betti_0 / Betti_1 diagnostic"]["evidence"]["reference_signature"]["Betti_1"], 22)

    def test_dashboard_best_case_is_bounded(self) -> None:
        self.assertEqual(
            self.result["best_case_classification"],
            [
                "LOCAL_SIGNAL_STABLE",
                "MATHEMATICALLY_INTERESTING",
                "ENGINEERING_FEASIBLE",
                "FORMALLY_PARTIAL",
                "GLOBAL_CLAIM_OPEN",
            ],
        )
        self.assertIn("no Moonshine proof", self.result["claim_ceiling"])
        self.assertIn("physical bridge", self.result["claim_ceiling"])

    def test_outputs_exist_and_are_reproducible(self) -> None:
        for filename in (
            "geometry_bridge_recoverability_report.md",
            "geometry_bridge_recoverability_result.json",
        ):
            self.assertTrue((self.output_dir / filename).exists(), filename)
        with tempfile.TemporaryDirectory() as other:
            rerun = dashboard.run_dashboard(self.config, Path(other))
        self.assertEqual(self.result["result_hash"], rerun["result_hash"])

    def test_dashboard_json_round_trips(self) -> None:
        payload = json.loads((self.output_dir / "geometry_bridge_recoverability_result.json").read_text(encoding="utf-8"))
        self.assertEqual(payload["result_hash"], self.result["result_hash"])
        self.assertEqual(payload["status"], "PASS_WITH_FRAGILITY")


if __name__ == "__main__":
    unittest.main()
