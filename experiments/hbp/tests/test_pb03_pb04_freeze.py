from __future__ import annotations

import json
import subprocess
import unittest
from pathlib import Path

from experiments.hbp.pb03_temporal_recovery.loader import build_loader_summary as build_pb03_loader_summary
from experiments.hbp.pb04_cross_domain_transfer.loader import build_loader_summary as build_pb04_loader_summary


REPO_ROOT = Path(__file__).resolve().parents[3]
PB03 = REPO_ROOT / "experiments/hbp/pb03_temporal_recovery"
PB04 = REPO_ROOT / "experiments/hbp/pb04_cross_domain_transfer"


def run_script(path: Path) -> dict[str, object]:
    completed = subprocess.run(
        ["uv", "run", "python", str(path)],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    return json.loads(completed.stdout)


class PB03PB04FreezeTests(unittest.TestCase):
    def test_pb03_loader_reads_frozen_powergraph_cascades(self) -> None:
        summary = build_pb03_loader_summary()
        self.assertEqual(summary["bridge_id"], "PB-03-TEMPORAL-RECOVERY")
        self.assertEqual(summary["selection_status"], "frozen")
        self.assertEqual(summary["family"], "ieee24")
        self.assertEqual(summary["file_count"], 8)
        self.assertEqual(summary["file_summaries"]["of_reg.mat"]["shape"], [1, 21500])

    def test_pb04_loader_reads_proxy_source_and_target(self) -> None:
        summary = build_pb04_loader_summary()
        self.assertEqual(summary["bridge_id"], "PB-04-CROSS-DOMAIN-TRANSFER")
        self.assertEqual(summary["selection_status"], "frozen")
        self.assertTrue(summary["source_exists"])
        self.assertTrue(summary["target_exists"])
        self.assertGreaterEqual(summary["source_file_count"], 1)
        self.assertGreaterEqual(summary["target_file_count"], 1)

    def test_pb03_runner_stub_reports_frozen_artifacts(self) -> None:
        payload = run_script(PB03 / "run_pb03_temporal_recovery.py")
        self.assertIn(payload["status"], {"PREDICTION_NOT_DISTINCT_FROM_BASELINES", "CONTROL_INVALID", "EMPIRICAL_BRIDGE_CANDIDATE"})
        result = json.loads((PB03 / "results" / "pb03_result.json").read_text(encoding="utf-8"))
        self.assertEqual(result["bridge_id"], "PB-03-TEMPORAL-RECOVERY")
        self.assertEqual(result["case_counts"]["holdout"], 24)
        self.assertEqual(len(result["controls"]), 4)
        self.assertTrue((PB03 / "results" / "pb03_report.md").exists())
        self.assertIn(result["status"], {"PREDICTION_NOT_DISTINCT_FROM_BASELINES", "CONTROL_INVALID", "EMPIRICAL_BRIDGE_CANDIDATE"})

    def test_pb04_runner_stub_reports_frozen_artifacts(self) -> None:
        payload = run_script(PB04 / "run_pb04_cross_domain_transfer.py")
        self.assertIn(payload["status"], {"PREDICTION_NOT_DISTINCT_FROM_BASELINES", "CONTROL_INVALID", "EMPIRICAL_BRIDGE_CANDIDATE"})
        result = json.loads((PB04 / "results" / "pb04_result.json").read_text(encoding="utf-8"))
        self.assertEqual(result["bridge_id"], "PB-04-CROSS-DOMAIN-TRANSFER")
        self.assertEqual(result["source_case_count"], 31)
        self.assertEqual(len(result["controls"]), 4)
        self.assertTrue((PB04 / "results" / "pb04_report.md").exists())
        self.assertIn(result["status"], {"PREDICTION_NOT_DISTINCT_FROM_BASELINES", "CONTROL_INVALID", "EMPIRICAL_BRIDGE_CANDIDATE"})


if __name__ == "__main__":
    unittest.main()
