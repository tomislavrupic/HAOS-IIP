from __future__ import annotations

import itertools
import unittest

import numpy as np

from experiments.emergence_ladder.rung3_distributed_parity.alphabet import RelationalAlphabet
from experiments.emergence_ladder.rung3_distributed_parity.decoder import decode_local
from experiments.emergence_ladder.rung3_distributed_parity.parity_architecture import (
    DistributedParityMemory,
    global_check_matrix,
)


class CodeDistanceDecoderTests(unittest.TestCase):
    def setUp(self) -> None:
        self.alphabet = RelationalAlphabet.build(8)
        self.symbols = np.asarray(([0, 1, 0] * 16), dtype=np.uint8)
        self.memory = DistributedParityMemory.encode(self.symbols, self.alphabet)

    def test_minimum_distance_is_three(self) -> None:
        matrix = global_check_matrix(1)
        codewords = []
        for bits in itertools.product((0, 1), repeat=3):
            value = np.asarray(bits, dtype=np.uint8)
            if not np.any((matrix @ value) % 2):
                codewords.append(value)
        nonzero_weights = [int(np.sum(value)) for value in codewords if np.any(value)]
        self.assertEqual(min(nonzero_weights), 3)

    def test_every_single_error_is_corrected(self) -> None:
        for index in range(len(self.symbols)):
            observed = self.symbols.copy()
            observed[index] ^= 1
            result = decode_local(observed, self.memory, self.alphabet)
            np.testing.assert_array_equal(result.target_symbols, self.symbols)
            self.assertTrue(result.syndrome_eliminated)
            self.assertLessEqual(max(result.visible_bits_per_block), 5)

    def test_two_errors_can_miscorrect_to_wrong_codeword(self) -> None:
        observed = self.symbols.copy()
        observed[0] ^= 1
        observed[1] ^= 1
        result = decode_local(observed, self.memory, self.alphabet)
        self.assertTrue(result.syndrome_eliminated)
        self.assertFalse(np.array_equal(result.target_symbols, self.symbols))

    def test_missing_check_never_falls_back_globally(self) -> None:
        observed = self.symbols.copy()
        observed[0] ^= 1
        available = np.ones((16, 2), dtype=bool)
        available[0, 0] = False
        result = decode_local(observed, self.memory, self.alphabet, available)
        self.assertEqual(result.block_status[0], "unresolved_missing_check")
        self.assertEqual(result.visible_bits_per_block[0], 4)
        self.assertFalse(np.array_equal(result.target_symbols, self.symbols))


if __name__ == "__main__":
    unittest.main()

