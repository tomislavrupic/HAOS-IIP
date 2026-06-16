#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from experiments.hbp.hbp_validation import stable_hash


HERE = Path(__file__).resolve().parent
RESULTS_DIR = HERE / "results"
RESULT_PATH = RESULTS_DIR / "pb02_result.json"
REQUIRED_OUTPUTS = [
    RESULTS_DIR / "dataset_validation.json",
    RESULTS_DIR / "split_manifest.json",
    RESULTS_DIR / "baseline_results.csv",
    RESULTS_DIR / "haos_results.csv",
    RESULTS_DIR / "incremental_value.csv",
    RESULTS_DIR / "control_results.csv",
    RESULTS_DIR / "holdout_predictions.csv",
    RESULTS_DIR / "uncertainty_report.json",
    RESULTS_DIR / "pb02_result.json",
    RESULTS_DIR / "pb02_report.md",
]


def fail(message: str) -> None:
    raise SystemExit(f"PB-02 holdout check failed: {message}")


def main() -> int:
    missing = [str(path.relative_to(HERE)) for path in REQUIRED_OUTPUTS if not path.exists()]
    if missing:
        fail(f"missing outputs: {missing}")
    payload = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    comparable = dict(payload)
    comparable.pop("result_hash", None)
    if payload.get("result_hash") != stable_hash(comparable, "pb02_result_"):
        fail("result hash mismatch")
    summary = {
        "status": payload.get("status"),
        "result_hash": payload.get("result_hash"),
        "controls": len(payload.get("controls", [])),
        "holdout_n": payload.get("case_counts", {}).get("holdout"),
        "outputs": [str(path.relative_to(HERE)) for path in REQUIRED_OUTPUTS],
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
