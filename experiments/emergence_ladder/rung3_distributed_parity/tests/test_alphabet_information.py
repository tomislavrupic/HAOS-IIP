from __future__ import annotations

import unittest

import numpy as np

from experiments.emergence_ladder.rung3_distributed_parity.alphabet import RelationalAlphabet
from experiments.emergence_ladder.rung3_distributed_parity.parity_architecture import (
    DistributedParityMemory,
    assert_memory_schema,
    information_accounting,
)
from experiments.emergence_ladder.rung3_recovery_trajectory_v2.fixtures import initial_state, make_fixture


class AlphabetInformationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.fixture = make_fixture(8)
        self.alphabet = RelationalAlphabet.build(8)
        self.state = initial_state(self.fixture, 7001)

    def test_alphabet_shape_and_locality(self) -> None:
        self.assertEqual(self.alphabet.incidence.shape, (48, 64))
        self.assertEqual(len(self.alphabet.blocks), 16)
        self.assertTrue(all(np.count_nonzero(row) == 2 for row in self.alphabet.incidence))
        self.assertTrue(all(len(self.alphabet.local_nodes(block)) == 4 for block in range(16)))

    def test_distinct_continuous_states_share_memory(self) -> None:
        second = 1.75 * self.state + 3.25
        first_symbols = self.alphabet.encode(self.state)
        second_symbols = self.alphabet.encode(second)
        first_memory = DistributedParityMemory.encode(first_symbols, self.alphabet)
        second_memory = DistributedParityMemory.encode(second_symbols, self.alphabet)
        self.assertFalse(np.allclose(self.state, second))
        np.testing.assert_array_equal(first_memory.parity_array(), second_memory.parity_array())

    def test_information_budget_is_not_checkpoint(self) -> None:
        accounting = information_accounting()
        self.assertEqual(accounting["stored_parity_bits"], 32)
        self.assertEqual(accounting["full_checkpoint_bits"], 4096)
        self.assertEqual(accounting["global_binary_nullity"], 16)
        self.assertFalse(accounting["continuous_state_reconstruction_possible"])
        assert_memory_schema()


if __name__ == "__main__":
    unittest.main()

