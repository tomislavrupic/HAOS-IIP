#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from experiments.emergence_ladder.rung3_recovery_trajectory_v2.run_experiment import (
    CONTRACT_PATH,
    FINAL_DIR,
    ROOT,
    file_hash,
    stable_hash,
)


REQUIRED = [
    ROOT / "MECHANISM_SELECTION.md",
    CONTRACT_PATH,
    ROOT / "calibration/parameter_selection.json",
    ROOT / "calibration/calibration_runs.csv",
    ROOT / "calibration/derived_thresholds.json",
    ROOT / "calibration/validation_runs.csv",
    FINAL_DIR / "seed_registry.json",
    FINAL_DIR / "run_level_results.csv",
    FINAL_DIR / "control_results.csv",
    FINAL_DIR / "aggregate_result.json",
    FINAL_DIR / "uncertainty.json",
    FINAL_DIR / "CLASSIFICATION.md",
    ROOT / "THEORY_INTERPRETATION.md",
]


def check() -> dict[str, object]:
    missing = [str(path.relative_to(ROOT)) for path in REQUIRED if not path.exists()]
    if missing:
        raise ValueError(f"missing RT-02 artifacts: {missing}")
    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    result = json.loads((FINAL_DIR / "aggregate_result.json").read_text(encoding="utf-8"))
    comparable = dict(result)
    observed_hash = comparable.pop("result_hash")
    if observed_hash != stable_hash("el_r3_rt_02", comparable):
        raise ValueError("aggregate result hash mismatch")
    if result["contract_hash"] != stable_hash("el_r3_rt_02_contract", contract):
        raise ValueError("contract hash mismatch")
    expected_sources = {
        "mechanism": file_hash(ROOT / "mechanism.py"),
        "controls": file_hash(ROOT / "controls.py"),
        "runner": file_hash(ROOT / "run_experiment.py"),
    }
    if result["source_hashes"] != expected_sources:
        raise ValueError("source hash mismatch")
    if result["classification"] not in {
        "RUNG_3_SUPPORTED",
        "PARTIAL_RECOVERY_ONLY",
        "MECHANISM_NEGATIVE",
        "CONTROL_INVALID",
        "REFERENCE_LEAKAGE",
        "IDENTITY_FAILURE",
        "TRIVIAL_ATTRACTOR",
        "UNSTABLE_RESTORATION",
        "INCONCLUSIVE",
    }:
        raise ValueError("invalid classification")
    if not result["leakage_guard"]["memory_stores_only_orientation_bits"]:
        raise ValueError("memory leakage guard failed")
    return {"status": "ok", "classification": result["classification"], "result_hash": result["result_hash"]}


def main() -> int:
    print(json.dumps(check(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
