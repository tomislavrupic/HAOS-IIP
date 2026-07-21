from __future__ import annotations

import inspect
import json
import copy
import unittest

from experiments.emergence_ladder.rung3_distributed_parity.bridge import correction_vector, simulate_bridge
from experiments.emergence_ladder.rung3_distributed_parity.decoder import decode_local
from experiments.emergence_ladder.rung3_distributed_parity.run_experiment import CONTRACT_PATH, leakage_guard, memory_proof


class LocalityLeakageTests(unittest.TestCase):
    def test_mechanism_signatures_exclude_hidden_reference(self) -> None:
        forbidden = ("reference", "future", "function", "score", "label", "checkpoint")
        for function in (decode_local, correction_vector, simulate_bridge):
            names = inspect.signature(function).parameters
            self.assertFalse(any(fragment in name.lower() for fragment in forbidden for name in names))

    def test_contract_seed_partitions_and_guards(self) -> None:
        contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
        self.assertTrue(all(leakage_guard(contract).values()))

    def test_empirical_memory_non_uniqueness(self) -> None:
        proof = memory_proof()
        self.assertTrue(proof["empirical_non_uniqueness_pass"])
        self.assertFalse(proof["continuous_state_reconstruction_possible"])

    def test_function_label_leakage_positive_control_is_detected(self) -> None:
        contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
        contaminated = copy.deepcopy(contract)
        contaminated["architecture"]["function_label"] = "recovered"
        self.assertFalse(leakage_guard(contaminated)["function_label_absent_from_contract_mechanism"])


if __name__ == "__main__":
    unittest.main()
