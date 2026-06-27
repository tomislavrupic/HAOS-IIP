from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "arithmetic_source_validation.py"
SPEC = importlib.util.spec_from_file_location("arithmetic_source_validation", MODULE_PATH)
assert SPEC and SPEC.loader
source_validation = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = source_validation
SPEC.loader.exec_module(source_validation)


class ArithmeticSourceValidationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = source_validation.SourceValidationConfig()
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmpdir.name)
        cls.result = source_validation.run_source_validation(cls.config, cls.output_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmpdir.cleanup()

    def test_pinned_manifest_entries_have_required_provenance_fields(self) -> None:
        manifest = json.loads((self.output_dir / "source_manifest_pinned.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["status"], "PASS")
        self.assertEqual(manifest["classification"], "PINNED_SOURCE_VALIDATION_ONLY")
        self.assertEqual(len(manifest["source_entries"]), 5)
        for entry in manifest["source_entries"]:
            for field in ("claim", "source", "retrieval_date", "hash_or_version", "used_in"):
                self.assertIn(field, entry)
                self.assertTrue(entry[field])
            self.assertEqual(entry["retrieval_date"], "2026-06-20")

    def test_validates_all_required_constant_families(self) -> None:
        claims = {validation["claim"]: validation for validation in self.result["validations"]}
        for claim in (
            "15 supersingular primes",
            "Monster order prime exponents",
            "j-function witness coefficients and sums",
            "Monster irrep dimensions used by witness subset",
            "Gaussian-prime representatives",
        ):
            self.assertIn(claim, claims)
            self.assertEqual(claims[claim]["status"], "PASS")
        self.assertEqual(claims["15 supersingular primes"]["observed"], list(source_validation.PINNED_SUPERSINGULAR_PRIMES))

    def test_j_witnesses_are_exact_against_pinned_coefficients(self) -> None:
        witness_validation = next(
            validation for validation in self.result["validations"] if validation["claim"] == "j-function witness coefficients and sums"
        )
        self.assertTrue(all(row["status"] == "PASS" for row in witness_validation["rows"]))
        self.assertTrue(all(row["residual"] == 0 for row in witness_validation["rows"]))

    def test_irrep_factor_support_is_visible(self) -> None:
        irrep_validation = next(
            validation for validation in self.result["validations"] if validation["claim"] == "Monster irrep dimensions used by witness subset"
        )
        nontrivial = [row for row in irrep_validation["rows"] if row["dimension"] > 1]
        self.assertEqual(len(nontrivial), 5)
        self.assertTrue(all(row["factor_support_subset"] for row in nontrivial))

    def test_result_remains_claim_gated(self) -> None:
        self.assertEqual(self.result["status"], "PASS")
        self.assertEqual(self.result["classification"], "PINNED_SOURCE_VALIDATION_ONLY")
        self.assertIn("SOURCE_MANIFEST_PINNED_BUILT", self.result["labels"])
        self.assertIn("MOONSHINE_PROOF_NOT_ESTABLISHED", self.result["labels"])
        self.assertIn("PHYSICAL_BRIDGE_NOT_ESTABLISHED", self.result["labels"])
        self.assertIn("not a proof of Monstrous Moonshine", self.result["non_claims"])

    def test_outputs_exist_and_are_reproducible(self) -> None:
        for filename in (
            "source_manifest_pinned.json",
            "source_validation_report.md",
            "arithmetic_source_validation_result.json",
        ):
            self.assertTrue((self.output_dir / filename).exists(), filename)
        with tempfile.TemporaryDirectory() as other:
            rerun = source_validation.run_source_validation(self.config, Path(other))
        self.assertEqual(self.result["result_hash"], rerun["result_hash"])
        self.assertEqual(self.result["manifest_hash"], rerun["manifest_hash"])


if __name__ == "__main__":
    unittest.main()
