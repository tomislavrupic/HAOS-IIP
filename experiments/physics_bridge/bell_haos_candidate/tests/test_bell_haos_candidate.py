from __future__ import annotations

import importlib.util
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_bell_haos_candidate.py"
SPEC = importlib.util.spec_from_file_location("bell_haos_candidate", MODULE_PATH)
assert SPEC and SPEC.loader
bell = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = bell
SPEC.loader.exec_module(bell)


class BellHaosCandidateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = bell.BellCandidateConfig(trials_per_setting_pair=64, seeds=(7201, 7202))
        cls.tmpdir = tempfile.TemporaryDirectory()
        cls.output_dir = Path(cls.tmpdir.name)
        cls.result = bell.run_bridge(cls.config, cls.output_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.tmpdir.cleanup()

    def test_local_code_path_isolation(self) -> None:
        signature = str(bell.inspect.signature(bell.local_sanity_outcome))
        self.assertIn("local_setting", signature)
        self.assertNotIn("remote", signature)
        signature_b2 = str(bell.inspect.signature(bell.haos_local_recoverability_outcome))
        self.assertIn("local_setting", signature_b2)
        self.assertNotIn("remote", signature_b2)

    def test_deterministic_reproducibility(self) -> None:
        with tempfile.TemporaryDirectory() as other:
            rerun = bell.run_bridge(self.config, Path(other))
        self.assertEqual(self.result["result_hash"], rerun["result_hash"])

    def test_binary_outcome_validity(self) -> None:
        trials_path = self.output_dir / "candidate_trials.csv"
        with trials_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for index, row in enumerate(reader):
                self.assertIn(row["alice_outcome"], {"-1", "1"})
                self.assertIn(row["bob_outcome"], {"-1", "1"})
                if index > 200:
                    break

    def test_correlation_and_chsh_calculation(self) -> None:
        config = bell.BellCandidateConfig(trials_per_setting_pair=1, seeds=(1,))
        source = bell.synthetic_source_state(1, 0, config)
        trials = []
        products = {"E00": 1, "E01": 1, "E10": 1, "E11": -1}
        for alice, bob in config.setting_pairs:
            key = bell.pair_key(alice, bob)
            outcome = bell.BellOutcome(1, products[key], 0.0, 0.0)
            trials.append(
                bell.BellTrial(
                    "manual",
                    "manual",
                    "test",
                    1,
                    0,
                    bell.pair_label(alice, bob),
                    alice,
                    bob,
                    source.source_id,
                    source.source_bucket,
                    source.source_hash,
                    outcome,
                    "manual",
                )
            )
        table = bell.correlation_table("manual", "manual", "test", "test", trials, config)
        self.assertEqual(table.correlations["E00"], 1.0)
        self.assertEqual(table.S, 4.0)

    def test_no_signalling_metric(self) -> None:
        result = next(item for item in self.result["candidate_results"] if item["candidate_id"] == "local_sanity_candidate")
        self.assertLessEqual(result["no_signalling_max_deviation"], self.config.no_signalling_tolerance)

    def test_source_setting_independence(self) -> None:
        result = next(item for item in self.result["candidate_results"] if item["candidate_id"] == "haos_local_recoverability_candidate")
        self.assertLessEqual(result["setting_source_mutual_information_bits"], self.config.setting_independence_mi_tolerance_bits)

    def test_post_selection_detection(self) -> None:
        result = next(item for item in self.result["candidate_results"] if item["candidate_id"] == "post_selection_positive_control")
        self.assertIn("POST_SELECTION_DETECTED", result["verdict_labels"])
        self.assertGreater(result["correlation_table"]["rejected_trials"], 0)

    def test_target_leakage_guard(self) -> None:
        detected, hits = bell.detect_target_leakage(["value = 2 * sqrt(2)\ncorrelation = -cos(theta_a - theta_b)"])
        self.assertTrue(detected)
        self.assertTrue(hits)
        clean, clean_hits = bell.detect_target_leakage(bell.generator_sources([bell.local_sanity_outcome]))
        self.assertFalse(clean)
        self.assertEqual(clean_hits, [])
        b3_clean, b3_hits = bell.detect_target_leakage(
            bell.generator_sources(
                [
                    bell.b3_source_features,
                    bell.b3_orientation_score,
                    bell.b3_joint_closure_costs,
                    bell.b3_sample_product,
                    bell.b3_joint_closure_trials,
                ]
            )
        )
        self.assertFalse(b3_clean)
        self.assertEqual(b3_hits, [])

    def test_b3_precommitment_authorizes_only_explicit_invoice(self) -> None:
        contract = json.loads((self.output_dir / "b3_precommitment_contract.json").read_text(encoding="utf-8"))
        self.assertTrue(contract["implementation_authorized"])
        self.assertEqual(contract["authorization_status"], "IMPLEMENTATION_AUTHORIZED_WITH_EXPLICIT_INVOICE")
        self.assertTrue(contract["joint_probability_model"]["provided"])
        self.assertIn("FACTORIZATION_RELAXED", contract["bell_assumption_invoice"]["selected_invoice"])
        self.assertIn("OUTCOME_INDEPENDENCE_RELAXED", contract["bell_assumption_invoice"]["selected_invoice"])

    def test_b3_joint_probability_is_nonfactorizable(self) -> None:
        source = bell.read_phase18_sources(self.config)[0]
        costs = bell.b3_joint_closure_costs(source, "a0", "b0", self.config)
        p_plus = costs["p_product_plus"]
        self.assertNotAlmostEqual(p_plus, 0.5)
        joint_pp = 0.5 * p_plus
        product_of_marginals = 0.25
        self.assertNotAlmostEqual(joint_pp, product_of_marginals)

    def test_b3_diagnostics_are_recorded(self) -> None:
        result = next(item for item in self.result["candidate_results"] if item["candidate_id"] == "haos_joint_closure_candidate")
        self.assertEqual(result["stage"], "B3")
        self.assertIn(result["assumption_ledger"]["setting_independence"], {"SETTING_INDEPENDENCE_PASS", "SETTING_INDEPENDENCE_FAIL"})
        self.assertEqual(result["assumption_ledger"]["factorization"], "FACTORIZATION_RELAXED_BY_EXPLICIT_JOINT_CLOSURE_COST")
        self.assertTrue(result["assumption_ledger"]["joint_setting_access_declared"])
        self.assertFalse(result["assumption_ledger"]["remote_setting_access"])
        self.assertTrue((self.output_dir / "b3_joint_cost_audit.csv").exists())

    def test_finite_sample_local_overshoot_classified_correctly(self) -> None:
        table = bell.BellCorrelationTable(
            candidate_id="local_sanity_candidate",
            run_id="fake",
            stage="B1",
            sample_role="candidate",
            correlations={"E00": 0.51, "E01": 0.51, "E10": 0.51, "E11": -0.50},
            standard_errors={"E00": 0.05, "E01": 0.05, "E10": 0.05, "E11": 0.05},
            trials_by_pair={"E00": 10, "E01": 10, "E10": 10, "E11": 10},
            alice_marginals_by_pair={},
            bob_marginals_by_pair={},
            S=2.03,
            ci_low=1.84,
            ci_high=2.22,
            retained_trials=40,
            rejected_trials=0,
            retention_rate=1.0,
        )
        labels = bell.classify_result(
            "local_sanity_candidate",
            table,
            "NO_SIGNALLING_PASS",
            "SETTING_INDEPENDENCE_PASS",
            False,
            False,
            False,
            bell.BellCandidateConfig(),
        )
        self.assertIn("CHSH_VIOLATION_NOT_DETECTED", labels)
        self.assertNotIn("CHSH_VIOLATION_DETECTED", labels)

    def test_intentionally_invalid_controls_detected(self) -> None:
        leakage = next(item for item in self.result["candidate_results"] if item["candidate_id"] == "source_setting_leakage_positive_control")
        post = next(item for item in self.result["candidate_results"] if item["candidate_id"] == "post_selection_positive_control")
        self.assertIn("SETTING_INDEPENDENCE_FAIL", leakage["verdict_labels"])
        self.assertIn("POST_SELECTION_DETECTED", post["verdict_labels"])

    def test_frozen_reference_sidecar_unchanged(self) -> None:
        before = bell.reference_sidecar_hashes()
        with tempfile.TemporaryDirectory() as other:
            bell.run_bridge(self.config, Path(other))
        after = bell.reference_sidecar_hashes()
        self.assertEqual(before, after)

    def test_required_outputs_exist(self) -> None:
        for filename in (
            "precommitment_contract.json",
            "assumption_ledger.json",
            "candidate_trials.csv",
            "candidate_correlations.json",
            "b3_precommitment_contract.json",
            "b3_precommitment_report.md",
            "b3_joint_cost_audit.csv",
            "no_signalling_diagnostics.csv",
            "setting_independence_diagnostics.csv",
            "control_results.csv",
            "bell_candidate_report.md",
            "bell_candidate_result.json",
        ):
            self.assertTrue((self.output_dir / filename).exists(), filename)
        result = json.loads((self.output_dir / "bell_candidate_result.json").read_text(encoding="utf-8"))
        self.assertIn("CST_BELL_STATUS_OPEN", result["labels"])


if __name__ == "__main__":
    unittest.main()
