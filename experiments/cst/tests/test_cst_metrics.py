from __future__ import annotations

import unittest

from experiments.cst.cst_metrics import forced_fail_probe, jaccard_distance, sequence_rank_distance


class CSTMetricTests(unittest.TestCase):
    def test_forced_fail_probe_returns_fail(self) -> None:
        self.assertTrue(forced_fail_probe())

    def test_jaccard_distance_bounds(self) -> None:
        self.assertEqual(jaccard_distance(["a"], ["a"]), 0.0)
        self.assertEqual(jaccard_distance(["a"], ["b"]), 1.0)
        self.assertGreater(jaccard_distance(["a", "b"], ["b", "c"]), 0.0)
        self.assertLessEqual(jaccard_distance(["a", "b"], ["b", "c"]), 1.0)

    def test_sequence_rank_distance_bounds(self) -> None:
        self.assertEqual(sequence_rank_distance(["a", "b"], ["a", "b"]), 0.0)
        self.assertGreater(sequence_rank_distance(["a", "b"], ["b", "a"]), 0.0)
        self.assertLessEqual(sequence_rank_distance(["a", "b"], ["b", "a"]), 1.0)


if __name__ == "__main__":
    unittest.main()
