from __future__ import annotations

import unittest

from experiments.cst.run_cst_benchmark import run_cst_benchmark


REQUIRED_CONTROLS = {
    "shuffled_structure_control",
    "randomized_edge_control",
    "degraded_signal_control",
    "label_permutation_control",
    "generic_graph_operator_control",
    "perturbation_free_control",
    "seed_repeat_control",
    "parameter_matched_null_control",
    "periodic_diagonal_augmented_control",
}


class CSTBenchmarkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.benchmark = run_cst_benchmark(write_outputs=False)

    def test_benchmark_is_deterministic(self) -> None:
        repeat = run_cst_benchmark(write_outputs=False)
        self.assertEqual(self.benchmark["result"].result_hash, repeat["result"].result_hash)
        self.assertEqual(self.benchmark["canonical_payload"], repeat["canonical_payload"])

    def test_branch_signature_id_is_not_run_instance_id(self) -> None:
        targets = [observation for observation in self.benchmark["observations"] if observation.observation_kind == "target"]
        self.assertGreater(len(targets), 0)
        for observation in targets:
            self.assertNotEqual(
                observation.provenance.run_instance_id,
                observation.closure_signature.branch_signature_id,
            )
        target_signature_ids = {observation.closure_signature.branch_signature_id for observation in targets}
        self.assertEqual(len(target_signature_ids), 1)

    def test_required_controls_are_reported(self) -> None:
        labels = {row["control_label"] for row in self.benchmark["control_distributions"]}
        self.assertTrue(REQUIRED_CONTROLS.issubset(labels))
        unavailable = {
            row["control_label"]
            for row in self.benchmark["control_distributions"]
            if row["availability"] == "unavailable"
        }
        self.assertEqual(unavailable, {"perturbation_free_control"})

    def test_toy_slice_can_fail_when_control_family_beats_target(self) -> None:
        result = self.benchmark["result"]
        self.assertEqual(result.verdict, "FAIL")
        joined_reasons = "\n".join(result.reasons)
        self.assertIn("negative control family shuffled_structure_control", joined_reasons)
        self.assertIn("negative control family randomized_edge_control", joined_reasons)
        self.assertIn("negative control family periodic_diagonal_augmented_control", joined_reasons)

    def test_optional_recoverability_scalar_can_be_disabled(self) -> None:
        scalarless = run_cst_benchmark(write_outputs=False, scalar_enabled_override=False)
        vector = scalarless["result"].recoverability_vector
        self.assertIsNone(vector.optional_scalar)
        self.assertIsNone(vector.scalar_formula)
        self.assertEqual(vector.authoritative_mode, "vector")
        self.assertGreater(len(scalarless["target_distances"]), 0)
        self.assertIsNotNone(scalarless["target_distances"][0].scalar_distance)


if __name__ == "__main__":
    unittest.main()
