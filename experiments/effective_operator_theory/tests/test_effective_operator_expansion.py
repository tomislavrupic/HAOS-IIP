from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from experiments.effective_operator_theory.run_effective_operator_expansion import (
    EffectiveOperatorConfig,
    run_effective_operator_expansion,
)


class EffectiveOperatorExpansionTests(unittest.TestCase):
    def run_bundle(self) -> tuple[dict, Path]:
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        output_dir = Path(temp.name)
        result = run_effective_operator_expansion(EffectiveOperatorConfig(), output_dir)
        return result, output_dir

    def test_reproducible_frozen_hash(self) -> None:
        first, _ = self.run_bundle()
        second, _ = self.run_bundle()
        self.assertEqual(first["result_hash"], second["result_hash"])
        self.assertEqual(first["result_hash"], "effective_operator_a740ce933dd1dbd479931847")

    def test_claim_ceiling_blocks_physics_upgrade(self) -> None:
        result, output_dir = self.run_bundle()
        contract = json.loads((output_dir / "precommitment_contract.json").read_text(encoding="utf-8"))
        for nonclaim in [
            "EFT quantum gravity",
            "physical field theory",
            "spacetime derivation",
            "Lorentz invariance",
            "empirical physics",
            "ontology",
        ]:
            self.assertIn(nonclaim, contract["not_claimed"])
        self.assertIn("EFT_QG_NOT_DERIVED", result["labels"])
        self.assertIn("PHYSICAL_FIELD_THEORY_NOT_ESTABLISHED", result["labels"])
        self.assertIn("PHYSICAL_BRIDGE_NOT_ESTABLISHED", result["labels"])

    def test_allowed_terms_are_synthetic_and_bounded(self) -> None:
        _, output_dir = self.run_bundle()
        allowed = json.loads((output_dir / "allowed_terms.json").read_text(encoding="utf-8"))
        names = {term["name"] for term in allowed["terms"]}
        self.assertEqual(names, {"identity", "laplacian", "laplacian_squared"})
        self.assertEqual(allowed["mapping_status"]["physical_field"], "UNAVAILABLE")
        self.assertEqual(allowed["mapping_status"]["metric"], "UNAVAILABLE")
        self.assertEqual(allowed["mapping_status"]["stress_energy"], "UNAVAILABLE")

    def test_effective_operator_scaffold_passes_only_synthetic_gates(self) -> None:
        result, _ = self.run_bundle()
        fit = result["fit"]
        config = result["config"]
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["classification"], "SYNTHETIC_EFFECTIVE_OPERATOR_SCAFFOLD")
        self.assertGreater(fit["leading_laplacian_coeff"], config["min_leading_coefficient"])
        self.assertLessEqual(fit["correction_ratio"], config["max_correction_ratio"])
        self.assertGreaterEqual(fit["r2"], config["min_fit_r2"])
        self.assertLessEqual(
            fit["max_long_wavelength_relative_residual"],
            config["max_long_mode_relative_residual"],
        )

    def test_controls_detect_missing_coupling_and_unstable_sign(self) -> None:
        result, output_dir = self.run_bundle()
        controls = {row["control"]: row for row in result["controls"]}
        self.assertEqual(controls["identity_operator_control"]["observed_status"], "PASS")
        self.assertEqual(controls["unstable_sign_control"]["observed_status"], "PASS")
        self.assertLess(float(controls["unstable_sign_control"]["observed_leading_coefficient"]), 0.0)

        with (output_dir / "control_results.csv").open(encoding="utf-8") as handle:
            csv_rows = list(csv.DictReader(handle))
        self.assertEqual({row["observed_status"] for row in csv_rows}, {"PASS"})

    def test_sweep_contains_frozen_modes_and_cutoff_classes(self) -> None:
        _, output_dir = self.run_bundle()
        with (output_dir / "coefficient_sweep.csv").open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual([int(row["mode"]) for row in rows], list(range(1, 13)))
        long_modes = [int(row["mode"]) for row in rows if row["scale_class"] == "long_wavelength"]
        self.assertEqual(long_modes, [1, 2, 3, 4])
        self.assertTrue(all(float(row["relative_residual"]) >= 0.0 for row in rows))


if __name__ == "__main__":
    unittest.main()
