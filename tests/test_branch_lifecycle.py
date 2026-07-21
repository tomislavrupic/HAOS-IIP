from __future__ import annotations

import copy
import unittest

from scripts.check_branch_lifecycle import (
    LifecycleValidationError,
    ROOT,
    load_registry,
    render_summary,
    validate_registry,
)


class BranchLifecycleTests(unittest.TestCase):
    def setUp(self) -> None:
        self.registry = load_registry()

    def test_authoritative_registry_is_valid(self) -> None:
        validate_registry(self.registry, ROOT)

    def test_summary_is_deterministic(self) -> None:
        self.assertEqual(render_summary(self.registry), render_summary(copy.deepcopy(self.registry)))

    def test_closed_branch_cannot_authorize_expansion(self) -> None:
        mutated = copy.deepcopy(self.registry)
        row = next(item for item in mutated["branches"] if item["branch_id"] == "CST-V0.2.2")
        row["expansion_authorized"] = True
        with self.assertRaisesRegex(LifecycleValidationError, "inactive branch cannot authorize expansion"):
            validate_registry(mutated, ROOT)

    def test_quarantined_branch_cannot_authorize_implementation(self) -> None:
        mutated = copy.deepcopy(self.registry)
        row = next(item for item in mutated["branches"] if item["branch_id"] == "HBP-PB02-V1")
        row["implementation_authorized"] = True
        with self.assertRaisesRegex(LifecycleValidationError, "inactive branch cannot authorize implementation"):
            validate_registry(mutated, ROOT)

    def test_active_candidate_requires_precommitment_gate(self) -> None:
        mutated = copy.deepcopy(self.registry)
        active_id = mutated["active_priority_order"][0]
        row = next(item for item in mutated["branches"] if item["branch_id"] == active_id)
        row["next_precommitment"] = {}
        with self.assertRaisesRegex(LifecycleValidationError, "precommitment status is invalid"):
            validate_registry(mutated, ROOT)

    def test_unfrozen_precommitment_cannot_authorize_implementation(self) -> None:
        mutated = copy.deepcopy(self.registry)
        active_id = mutated["active_priority_order"][0]
        row = next(item for item in mutated["branches"] if item["branch_id"] == active_id)
        row["implementation_authorized"] = True
        with self.assertRaisesRegex(LifecycleValidationError, "implementation requires a frozen precommitment"):
            validate_registry(mutated, ROOT)

    def test_active_priority_must_cover_exact_active_set(self) -> None:
        mutated = copy.deepcopy(self.registry)
        mutated["active_priority_order"] = mutated["active_priority_order"][:-1]
        with self.assertRaisesRegex(LifecycleValidationError, "every and only active branch"):
            validate_registry(mutated, ROOT)

    def test_prior_negative_results_remain_preserved(self) -> None:
        rows = {row["branch_id"]: row for row in self.registry["branches"]}
        for branch_id in ("CST-V0.2.2", "BELL-HAOS-B0-B3.2.2", "HBP-PB01-V1", "HBP-PB02-V1", "HBP-PB03-V1", "HBP-PB04-V1"):
            self.assertTrue(rows[branch_id]["preserved_artifacts"])

    def test_active_queue_does_not_authorize_physics_claims(self) -> None:
        rows = {row["branch_id"]: row for row in self.registry["branches"]}
        for branch_id in self.registry["active_priority_order"]:
            ceiling = rows[branch_id]["claim_ceiling"].lower()
            self.assertTrue("no " in ceiling or "only" in ceiling)


if __name__ == "__main__":
    unittest.main()
