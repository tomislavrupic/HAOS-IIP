from __future__ import annotations

from pathlib import Path
import json
import unittest

from experiments.cst.cst_metrics import COMPONENT_NAMES
from experiments.cst.diagnostics.discriminative_instrument_repair import (
    DEFAULT_CONTRACT_PATH,
    PHASE_2_2_CONTROL_LABELS,
    evaluate_component_value_map,
    phase_2_2_labels,
    run_phase_2_2_repair,
)


class Phase22RepairTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.repair = run_phase_2_2_repair()

    def test_precommitment_contract_exists_independently(self) -> None:
        self.assertTrue(DEFAULT_CONTRACT_PATH.exists())
        contract = json.loads(DEFAULT_CONTRACT_PATH.read_text(encoding="utf-8"))
        self.assertTrue(contract["frozen_before_repaired_benchmark_execution"])
        self.assertEqual(contract["authorization"], "GO_REPAIR_ONLY")

    def test_provenance_hash_and_fine_signature_are_separate(self) -> None:
        for signature in self.repair["fine_signatures"]:
            self.assertNotEqual(signature["run_instance_id"], signature["coarse_signature_id"])
            self.assertNotEqual(signature["run_instance_id"], signature["fine_signature_record_id"])
            self.assertEqual(signature["coarse_signature_policy"], "grouping_only_not_equivalence")
            self.assertEqual(signature["fine_signature_hash_policy"], "names_representation_only_not_equivalence")

    def test_fine_signature_is_inspectable(self) -> None:
        signature = self.repair["fine_signatures"][0]
        self.assertIn("equivalence_payload", signature)
        self.assertEqual(
            set(signature["equivalence_payload"]),
            {
                "invariant_summary",
                "structural_graph_summary",
                "spectral_summary",
                "causal_summary",
                "temporal_sequence_summary",
                "stable_address_summary",
            },
        )
        for section in signature["equivalence_payload"].values():
            self.assertIn("values", section)
            self.assertIn("semantic_justification", section)
            self.assertTrue(section["contributes_to_equivalence"])

    def test_seed_independent_equivalence(self) -> None:
        target_rows = [
            row for row in self.repair["equivalence"]
            if row["mode"] == "strict" and row["pair_type"] == "target_target"
        ]
        self.assertGreater(len(target_rows), 0)
        self.assertTrue(all(row["equivalent"] for row in target_rows))

    def test_meaningful_non_equivalence_for_destructive_controls(self) -> None:
        rows = [
            row for row in self.repair["equivalence"]
            if row["mode"] == "strict" and row["pair_type"].startswith("target_control:phase_2_2")
        ]
        self.assertGreater(len(rows), 0)
        self.assertTrue(all(not row["equivalent"] for row in rows))

    def test_no_control_label_leakage_into_equivalence_payload(self) -> None:
        for signature in self.repair["fine_signatures"]:
            payload = json.dumps(signature["equivalence_payload"], sort_keys=True)
            self.assertNotIn("control_group_label", payload)
            self.assertNotIn("phase_2_2_control", payload)

    def test_unavailable_component_policy(self) -> None:
        values = {component: 0.0 for component in COMPONENT_NAMES}
        availability = {component: "available" for component in COMPONENT_NAMES}
        values["d_time"] = None
        availability["d_time"] = "unavailable"
        strict = evaluate_component_value_map(values, availability, self.repair["contract"], mode="strict")
        exploratory = evaluate_component_value_map(values, availability, self.repair["contract"], mode="exploratory")
        self.assertFalse(strict["equivalent"])
        self.assertTrue(exploratory["equivalent"])

    def test_strict_vs_exploratory_equivalence(self) -> None:
        values = {component: 0.0 for component in COMPONENT_NAMES}
        availability = {component: "available" for component in COMPONENT_NAMES}
        values["d_causal"] = 0.15
        strict = evaluate_component_value_map(values, availability, self.repair["contract"], mode="strict")
        exploratory = evaluate_component_value_map(values, availability, self.repair["contract"], mode="exploratory")
        self.assertFalse(strict["equivalent"])
        self.assertTrue(exploratory["equivalent"])

    def test_component_specific_control_behavior(self) -> None:
        affected_rows = [
            row for row in self.repair["validation_rows"]
            if row["expectation"] == "affected_should_increase"
        ]
        self.assertGreater(len(affected_rows), 0)
        self.assertTrue(all(row["passed"] for row in affected_rows))
        for label in PHASE_2_2_CONTROL_LABELS:
            self.assertTrue(any(row["control_label"] == label for row in affected_rows), label)

    def test_control_invalidation_conditions_and_forced_invalid_result(self) -> None:
        for contract in self.repair["contract"]["component_control_contracts"].values():
            self.assertIn("invalidation_condition", contract)
        fake_labels = phase_2_2_labels(
            equivalence=[
                {"mode": "strict", "pair_type": "target_target", "equivalent": True},
                {"mode": "strict", "pair_type": "known_positive:seed_repeat_control", "equivalent": True},
                {"mode": "strict", "pair_type": "target_control:phase_2_2_invariant_control", "equivalent": False},
            ],
            validation_rows=[{"expectation": "affected_should_increase", "passed": False}],
            collisions={"unexplained_fine_collision_count": 0},
            power={"power_gate_pass": True},
        )
        self.assertIn("CONTROL INVALID", fake_labels)
        self.assertIn("INSTRUMENT INVALID", fake_labels)

    def test_fixed_tolerance_determinism(self) -> None:
        repeat = run_phase_2_2_repair(write_outputs=False)
        self.assertEqual(self.repair["result"]["result_hash"], repeat["result"]["result_hash"])

    def test_threshold_boundary_sensitivity(self) -> None:
        tolerance = self.repair["contract"]["equivalence"]["strict"]["tolerances"]["d_inv"]
        values = {component: 0.0 for component in COMPONENT_NAMES}
        availability = {component: "available" for component in COMPONENT_NAMES}
        values["d_inv"] = tolerance - 0.01
        result = evaluate_component_value_map(values, availability, self.repair["contract"], mode="strict")
        self.assertTrue(result["boundary_sensitive"])
        self.assertIn("d_inv", result["boundary_hits"])

    def test_underpowered_verdict(self) -> None:
        self.assertIn("UNDERPOWERED", self.repair["result"]["labels"])
        self.assertEqual(self.repair["power"]["status"], "UNDERPOWERED")

    def test_scalar_disabled_official_verdict(self) -> None:
        self.assertFalse(self.repair["result"]["official_scalar_enabled"])

    def test_forced_target_not_distinct_result(self) -> None:
        self.assertIn("TARGET NOT DISTINCT", self.repair["result"]["labels"])
        self.assertNotIn("TARGET DISTINCT", self.repair["result"]["labels"])

    def test_required_outputs_exist(self) -> None:
        for path in self.repair["result"]["outputs"].values():
            self.assertTrue(Path(path).exists(), path)


if __name__ == "__main__":
    unittest.main()
