from __future__ import annotations

from pathlib import Path
import unittest

from experiments.cst.diagnostics.failure_localization import run_failure_localization


REQUIRED_OUTPUT_NAMES = {
    "failure_localization_report",
    "component_separation",
    "ablation_matrix",
    "signature_collision_report",
    "control_integrity_report",
    "timing_ablation",
    "normalization_audit",
    "seed_power_analysis",
    "failure_localization_result",
}


class FailureLocalizationTests(unittest.TestCase):
    def test_failure_localization_outputs_and_labels(self) -> None:
        result = run_failure_localization()
        self.assertTrue({"TARGET NOT DISTINCT", "METRIC INSUFFICIENT", "CONTROL INVALID", "UNDERPOWERED"}.issubset(result["labels"]))
        self.assertEqual(set(result["outputs"]), REQUIRED_OUTPUT_NAMES)
        for path in result["outputs"].values():
            self.assertTrue(Path(path).exists(), path)


if __name__ == "__main__":
    unittest.main()
