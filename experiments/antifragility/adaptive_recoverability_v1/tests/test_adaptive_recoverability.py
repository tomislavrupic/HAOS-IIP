from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from experiments.antifragility.adaptive_recoverability_v1.model import (
    adaptive_backup_fraction,
    apply_stress,
    disjoint_primary_locations,
    graph_metrics,
    interaction_budget,
    operator_identity,
    ring_with_backups,
)
from experiments.antifragility.adaptive_recoverability_v1.run_experiment import (
    CONTRACT_PATH,
    expected_row_count,
    generate_rows,
    run,
)


class AdaptiveRecoverabilityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))

    def test_contract_is_frozen_and_seed_partitions_are_disjoint(self) -> None:
        self.assertEqual(self.contract["status"], "FROZEN")
        final = set(self.contract["final_evaluation"]["seeds"])
        smoke = set(self.contract["smoke_partition"]["seeds"])
        self.assertFalse(final & smoke)

    def test_adaptation_is_monotone_and_capped(self) -> None:
        baseline = 0.05
        values = [adaptive_backup_fraction(baseline, severity, 0.22, 0.27) for severity in (0.0, 0.4, 0.7, 1.0)]
        self.assertEqual(values, sorted(values))
        self.assertEqual(values[0], baseline)
        self.assertLessEqual(values[-1], 0.27)

    def test_weight_budget_and_identity_are_preserved(self) -> None:
        baseline = ring_with_backups(24, 0.05)
        adapted = ring_with_backups(24, 0.27)
        self.assertAlmostEqual(interaction_budget(baseline), interaction_budget(adapted), places=12)
        self.assertGreaterEqual(operator_identity(baseline, adapted), 0.95)

    def test_stress_does_not_mutate_its_input(self) -> None:
        graph = ring_with_backups(24, 0.05)
        before = graph.copy()
        stressed = apply_stress(graph, "contiguous_primary_damage", 3, 3, 0.7)
        np.testing.assert_array_equal(graph, before)
        self.assertFalse(np.array_equal(stressed, graph))

    def test_disconnected_control_has_zero_algebraic_connectivity(self) -> None:
        graph = ring_with_backups(24, 0.0)
        disconnected = apply_stress(graph, "distributed_primary_damage", 0, 3, 1.0)
        metrics = graph_metrics(disconnected)
        self.assertAlmostEqual(metrics.algebraic_connectivity, 0.0, places=12)
        self.assertGreaterEqual(metrics.effective_conductance, 0.0)

    def test_disjoint_location_guard(self) -> None:
        self.assertTrue(disjoint_primary_locations(24, 0, 12, 3, "displaced_primary_block"))
        self.assertFalse(disjoint_primary_locations(24, 0, 1, 3, "displaced_primary_block"))
        self.assertTrue(disjoint_primary_locations(24, 0, 0, 3, "backup_channel_damage"))

    def test_smoke_schedule_and_outputs(self) -> None:
        rows = generate_rows(self.contract, "smoke")
        self.assertEqual(len(rows), expected_row_count(self.contract, "smoke"))
        self.assertTrue(all(row["partition"] == "smoke" for row in rows))
        with tempfile.TemporaryDirectory() as temp_dir:
            result = run("smoke", Path(temp_dir))
            self.assertEqual(result["candidate_id"], "HAOS-AF-AR-01")
            self.assertTrue((Path(temp_dir) / "results" / "per_run_results.jsonl").is_file())


if __name__ == "__main__":
    unittest.main()
