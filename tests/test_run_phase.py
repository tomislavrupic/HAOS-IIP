from __future__ import annotations

import unittest

import run_phase


class RunPhaseResolutionTests(unittest.TestCase):
    def test_registered_checker_is_selected(self) -> None:
        target = run_phase.resolve_target("3", check=True)
        self.assertEqual(target, run_phase.PHASE_CHECKERS["3"])

    def test_missing_checker_never_falls_back_to_builder(self) -> None:
        with self.assertRaisesRegex(ValueError, "Refusing to run its builder"):
            run_phase.resolve_target("6", check=True)

    def test_builder_mode_remains_available(self) -> None:
        target = run_phase.resolve_target("6", check=False)
        self.assertEqual(target, run_phase.PHASE_COMMANDS["6"][0])


if __name__ == "__main__":
    unittest.main()
