from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_monstrous_moonshine_diagnostic.py"
SPEC = importlib.util.spec_from_file_location("monstrous_moonshine_diagnostic", MODULE_PATH)
assert SPEC and SPEC.loader
diagnostic = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = diagnostic
SPEC.loader.exec_module(diagnostic)


class MonstrousMoonshineDiagnosticTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = diagnostic.MoonshineDiagnosticConfig()
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmpdir.name)
        cls.result = diagnostic.run_diagnostic(cls.config, cls.output_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmpdir.cleanup()

    def test_factorization_support_for_embedded_irreps(self) -> None:
        support = set(diagnostic.SUPERSINGULAR_PRIMES)
        for name, value in diagnostic.IRREP_DIMS.items():
            if value == 1:
                continue
            self.assertLessEqual(set(diagnostic.factor_integer(value)), support, name)

    def test_decomposition_witnesses_are_exact_in_reference(self) -> None:
        components = self.result["condition_results"]["known_positive"]["components"]
        self.assertEqual(components["decomposition"]["residual_abs_sum"], 0)
        self.assertEqual(components["decomposition"]["exact_witness_count"], len(diagnostic.J_DECOMPOSITION_WITNESSES))

    def test_contract_declares_claim_ceiling(self) -> None:
        contract = json.loads((self.output_dir / "precommitment_contract.json").read_text(encoding="utf-8"))
        self.assertEqual(contract["claim_ceiling"], "ARITHMETIC_SUPPORT_DIAGNOSTIC")
        self.assertIn("not a proof of Monstrous Moonshine", contract["non_claims"])
        self.assertIn("not a physical bridge", contract["non_claims"])

    def test_self_recovery_and_controls(self) -> None:
        known = self.result["condition_results"]["known_positive"]
        support_control = self.result["condition_results"]["nonsupersingular_prime_control"]
        exponent_control = self.result["condition_results"]["exponent_shuffle_control"]
        decomposition_control = self.result["condition_results"]["decomposition_break_control"]
        gaussian_control = self.result["condition_results"]["gaussian_class_swap_control"]
        self.assertEqual(known["observed_status"], "PASS")
        self.assertLessEqual(max(known["distances"].values()), self.config.invariance_max)
        self.assertEqual(support_control["observed_status"], "PASS")
        self.assertGreaterEqual(support_control["distances"]["support"], self.config.support_control_min)
        self.assertEqual(exponent_control["observed_status"], "PASS")
        self.assertGreaterEqual(exponent_control["distances"]["exponent"], self.config.exponent_control_min)
        self.assertEqual(decomposition_control["observed_status"], "PASS")
        self.assertGreaterEqual(decomposition_control["distances"]["decomposition"], self.config.decomposition_control_min)
        self.assertEqual(gaussian_control["observed_status"], "PASS")
        self.assertGreaterEqual(gaussian_control["distances"]["gaussian_class"], self.config.gaussian_class_control_min)

    def test_result_remains_claim_gated(self) -> None:
        self.assertEqual(self.result["status"], "PASS")
        self.assertEqual(self.result["classification"], "ARITHMETIC_SUPPORT_DIAGNOSTIC")
        self.assertIn("MOONSHINE_PROOF_NOT_ESTABLISHED", self.result["labels"])
        self.assertIn("PHYSICAL_BRIDGE_NOT_ESTABLISHED", self.result["labels"])

    def test_outputs_exist_and_are_reproducible(self) -> None:
        for filename in (
            "precommitment_contract.json",
            "source_manifest.json",
            "supersingular_prime_table.csv",
            "moonshine_witnesses.csv",
            "component_scores.csv",
            "control_results.csv",
            "monstrous_moonshine_diagnostic_report.md",
            "monstrous_moonshine_diagnostic_result.json",
        ):
            self.assertTrue((self.output_dir / filename).exists(), filename)
        with tempfile.TemporaryDirectory() as other:
            rerun = diagnostic.run_diagnostic(self.config, Path(other))
        self.assertEqual(self.result["result_hash"], rerun["result_hash"])

    def test_control_csv_contains_all_conditions(self) -> None:
        with (self.output_dir / "control_results.csv").open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual({row["condition"] for row in rows}, set(diagnostic.CONDITIONS))
        self.assertTrue(all(row["observed_status"] == "PASS" for row in rows))


if __name__ == "__main__":
    unittest.main()

