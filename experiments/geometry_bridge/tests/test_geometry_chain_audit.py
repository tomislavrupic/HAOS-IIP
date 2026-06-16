from __future__ import annotations

import importlib.util
import json
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_geometry_chain_audit.py"
SPEC = importlib.util.spec_from_file_location("geometry_chain_audit", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class GeometryChainAuditTests(unittest.TestCase):
    def test_contract_has_frozen_open_boundary(self) -> None:
        contract = json.loads((Path(__file__).resolve().parents[1] / "pregeometry_mechanism_contract.json").read_text())
        self.assertEqual(contract["bridge_id"], "GEO-05")
        self.assertEqual(contract["status"], "OPEN")
        self.assertIn("geometry", contract["boundary_statement"])
        self.assertIn("not a Bell mechanism", contract["boundary_statement"])

    def test_audit_result_marks_chain_open(self) -> None:
        result = MODULE.evaluate()
        self.assertEqual(result["bridge_id"], "GEO-05")
        self.assertEqual(result["status"], "OPEN")
        chain = {item["name"]: item for item in result["chain"]}
        self.assertEqual(chain["intrinsic_geometry"]["status"], "OPEN")
        self.assertIn(chain["observable_prediction"]["status"], {"OPEN", "PARTIAL", "VALID"})
        self.assertTrue((Path(__file__).resolve().parents[1] / "geometry_chain_audit_result.json").exists())
        self.assertTrue((Path(__file__).resolve().parents[1] / "geometry_chain_audit_report.md").exists())

    def test_objective_completion_audit_stays_incomplete(self) -> None:
        audit = json.loads((Path(__file__).resolve().parents[1] / "goal_completion_audit.json").read_text())
        self.assertEqual(audit["status"], "INCOMPLETE")
        statuses = {item["requirement"]: item["status"] for item in audit["requirements"]}
        self.assertEqual(statuses["HAOS-IIP current state is pregeometric / proto-geometric"], "SUPPORTED")
        self.assertEqual(statuses["Geometry requirement is likely necessary"], "SUPPORTED")
        self.assertEqual(statuses["Geometry alone sufficient for Bell is NO"], "SUPPORTED")
        self.assertEqual(statuses["Primary missing chain: geometry -> transformation semantics -> probability rule -> observable prediction"], "PARTIAL")
        self.assertEqual(statuses["Bell failures explained by missing geometry"], "PARTIAL")
        self.assertEqual(statuses["Bell failures proved to be caused by missing geometry"], "NOT_ESTABLISHED")

    def test_transfer_ledger_stays_partial_open(self) -> None:
        ledger = json.loads((Path(__file__).resolve().parents[1] / "geometry_transfer_ledger.json").read_text())
        self.assertEqual(ledger["bridge_id"], "GEO-05")
        self.assertEqual(ledger["status"], "OPEN")
        steps = {(item["from"], item["to"]): item for item in ledger["steps"]}
        self.assertEqual(steps[("geometry", "transformation_semantics")]["status"], "OPEN")
        self.assertEqual(steps[("transformation_semantics", "probability_rule")]["status"], "PARTIAL")
        self.assertEqual(steps[("probability_rule", "observable_prediction")]["status"], "PARTIAL")
        self.assertEqual(steps[("observable_prediction", "bell_causality_from_missing_geometry")]["status"], "NOT_ESTABLISHED")

    def test_audit_artifacts_are_consistent(self) -> None:
        root = Path(__file__).resolve().parents[1]
        transfer = json.loads((root / "geometry_transfer_consistency.json").read_text())
        gap = json.loads((root / "geometry_gap_map.json").read_text())
        ledger = json.loads((root / "geometry_transfer_ledger.json").read_text())
        summary = json.loads((root / "geometry_transfer_summary.json").read_text())
        self.assertEqual(transfer["status"], "OPEN")
        self.assertEqual(gap["status"], "OPEN")
        self.assertEqual(ledger["status"], "OPEN")
        self.assertEqual(summary["status"], "OPEN")
        self.assertEqual(transfer["chain"][-1]["status"], "NOT_ESTABLISHED")
        self.assertEqual(gap["layers"][-1]["status"], "NOT_ESTABLISHED")
        self.assertEqual(ledger["steps"][-1]["status"], "NOT_ESTABLISHED")
        self.assertIn("not sufficient", summary["conclusion"])


if __name__ == "__main__":
    unittest.main()
