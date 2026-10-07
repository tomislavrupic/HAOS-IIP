from __future__ import annotations

import copy
import unittest

import numpy as np

from experiments.emergence_ladder.rung3_mesoscopic_field.check_bundle import BundleError, validate_claim_language, validate_contract, validate_rows_semantics
from experiments.emergence_ladder.rung3_mesoscopic_field.controls import ALL_CONDITIONS, run_condition
from experiments.emergence_ladder.rung3_mesoscopic_field.fixtures import initial_state, perturbation
from experiments.emergence_ladder.rung3_mesoscopic_field.representation import QuantizerSpec, coarse_field, coarse_matrix, quantize, rank_and_nullity, restore_coarse_field
from experiments.emergence_ladder.rung3_mesoscopic_field.run_experiment import aggregate_rows, execute_rows, load_contract


class Rung3MesoscopicFieldTests(unittest.TestCase):
    def test_coarse_operator_rank_gram_and_nullity(self) -> None:
        for n_side in (8, 12, 16):
            with self.subTest(n_side=n_side):
                matrix = coarse_matrix(n_side)
                cell_values = (n_side // 4) ** 2
                self.assertEqual(np.linalg.matrix_rank(matrix), 16)
                self.assertTrue(np.allclose(matrix @ matrix.T, np.eye(16) / cell_values))
                self.assertEqual(rank_and_nullity(n_side), (16, n_side * n_side - 16))

    def test_decoder_restores_stored_quantized_field_exactly(self) -> None:
        for n_side in (8, 12, 16):
            with self.subTest(n_side=n_side):
                spec = QuantizerSpec()
                source = initial_state(n_side, 77)
                disrupted = source + perturbation(n_side, "smooth_physical_twist", 0.7, 81)
                stored = quantize(coarse_field(source, n_side), spec)
                restored = restore_coarse_field(disrupted, stored.decoded, n_side)
                self.assertLessEqual(np.max(np.abs(coarse_field(restored, n_side) - stored.decoded)), 1.0e-12)
                self.assertFalse(np.allclose(restored, source))

    def test_quantizer_is_half_even_fixed_range_and_rejects_nonfinite(self) -> None:
        spec = QuantizerSpec()
        values = spec.minimum + np.array([10.5, 11.5]) * spec.step
        self.assertEqual(quantize(values, spec).codes.tolist(), [10, 12])
        saturated = quantize(np.array([-4.0, 4.0]), spec)
        self.assertEqual(saturated.codes.tolist(), [0, 255])
        self.assertEqual(saturated.saturation_count, 2)
        with self.assertRaisesRegex(ValueError, "non-finite"):
            quantize(np.array([np.nan]), spec)

    def test_within_cell_perturbation_is_in_coarse_nullspace(self) -> None:
        for n_side in (8, 12, 16):
            delta = perturbation(n_side, "within_cell_balanced_gradient", 1.05, 88)
            self.assertLessEqual(np.max(np.abs(coarse_field(delta, n_side))), 1.0e-12)

    def test_all_conditions_are_deterministic_and_finite_on_smoke_fixture(self) -> None:
        n_side = 8
        spec = QuantizerSpec()
        source = initial_state(n_side, 9301)
        disrupted = source + perturbation(n_side, "coarse_cell_offsets", 0.35, 9301)
        target, _ = run_condition("target_quantized_coarse_field", source, disrupted, n_side, 9301, spec, 0.0)
        target_rms = float(np.sqrt(np.mean((target - disrupted) ** 2)))
        for condition in ALL_CONDITIONS:
            first, first_metadata = run_condition(condition, source, disrupted, n_side, 9301, spec, target_rms)
            second, second_metadata = run_condition(condition, source, disrupted, n_side, 9301, spec, target_rms)
            self.assertTrue(np.isfinite(first).all())
            self.assertTrue(np.array_equal(first, second))
            self.assertEqual(first_metadata, second_metadata)

    def test_frozen_contract_excludes_target_leakage_and_adaptive_state(self) -> None:
        contract = load_contract()
        validate_contract(contract)
        self.assertEqual(contract["information_accounting"]["payload_bits"], 128)
        self.assertFalse(contract["quantizer"]["adaptive_range"])
        self.assertFalse(contract["quantizer"]["per_sample_normalization"])
        self.assertNotIn("direct_four_functional_phasors", contract["conditions"]["admissible"])

    def test_semantic_checker_rejects_dropped_duplicate_and_mutated_rows(self) -> None:
        contract = load_contract()
        rows = execute_rows(contract, "smoke")
        validate_rows_semantics(rows, contract, "smoke")
        with self.assertRaisesRegex(BundleError, "schedule"):
            validate_rows_semantics(rows[:-1], contract, "smoke")
        with self.assertRaisesRegex(BundleError, "duplicate"):
            validate_rows_semantics(rows + [copy.deepcopy(rows[0])], contract, "smoke")
        mutated = copy.deepcopy(rows)
        next(row for row in mutated if row["condition"] == "target_quantized_coarse_field")["payload_bits"] = 127
        with self.assertRaisesRegex(BundleError, "payload"):
            validate_rows_semantics(mutated, contract, "smoke")
        mutated = copy.deepcopy(rows)
        next(row for row in mutated if row["condition"] == "target_quantized_coarse_field")["stored_quantized_field_error"] = 1.0
        with self.assertRaisesRegex(BundleError, "exactly"):
            validate_rows_semantics(mutated, contract, "smoke")

    def test_claim_checker_rejects_affirmative_overclaim(self) -> None:
        contract = load_contract()
        valid = "## Open\nfinite levels\n## Not claimed\nnone\nNO_HIGHER_RUNG_PROMOTION\n"
        validate_claim_language(valid, contract)
        with self.assertRaisesRegex(BundleError, "forbidden"):
            validate_claim_language(valid + "We establish universality.\n", contract)

    def test_failed_recovery_rate_without_scale_drop_is_partial_not_unstable(self) -> None:
        contract = load_contract()
        rows = []
        for n_side in (8, 12, 16):
            for seed in contract["final_evaluation"]["seeds"]:
                for family in contract["final_evaluation"]["perturbation_families"]:
                    for magnitude in contract["final_evaluation"]["magnitudes"]:
                        for condition in contract["conditions"]["admissible"] + contract["conditions"]["non_admissible_oracles"]:
                            recovery = 0.8 if condition == "target_quantized_coarse_field" else 0.0
                            recovered = condition == "target_quantized_coarse_field" and family != "within_cell_balanced_gradient"
                            rows.append({
                                "n_side": n_side, "seed": seed, "perturbation_family": family, "magnitude": magnitude,
                                "condition": condition, "function_applicable": True, "functional_recovery": recovery,
                                "recovered": recovered, "coarse_field_recovery": 1.0, "intervention_rms": 0.0,
                                "intervention_density": 0.0, "saturation_count": 0, "admissible": condition in contract["conditions"]["admissible"],
                                "payload_bits": 128 if condition not in {"compression_12_target_blind", "compression_8_target_blind"} else (96 if "12" in condition else 64),
                                "coarse_rank": 16, "microscopic_nullity": n_side * n_side - 16, "state_dimension": n_side * n_side,
                                "variance_ratio": 1.0, "stored_quantized_field_error": 0.0,
                                "edge_sign_identity": 0.5, "state_recovery_gain": 0.0,
                            })
        aggregate, _, _ = aggregate_rows(rows, contract)
        self.assertTrue(aggregate["gates"]["fixed_budget_scale_stable"])
        self.assertFalse(aggregate["gates"]["functional_recovery_rate_by_level"])
        self.assertEqual(aggregate["classification"], "PARTIAL_RECOVERY_ONLY")


if __name__ == "__main__":
    unittest.main()
