#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
PRECOMMITMENT_PATH = HERE / "precommitment_contract.json"
DATASET_MANIFEST_PATH = HERE / "dataset_manifest.json"
EXECUTION_CONTRACT_PATH = HERE / "execution_contract.json"
DEVELOPMENT_SPLIT_PATH = HERE / "development_split.json"
CALIBRATION_SPLIT_PATH = HERE / "calibration_split.json"
HOLDOUT_SPLIT_PATH = HERE / "holdout_split.json"


def fail(message: str) -> None:
    raise SystemExit(f"PB-02 bundle check failed: {message}")


def main() -> int:
    if not PRECOMMITMENT_PATH.exists():
        fail("missing precommitment_contract.json")
    if not DATASET_MANIFEST_PATH.exists():
        fail("missing dataset_manifest.json")
    if not EXECUTION_CONTRACT_PATH.exists():
        fail("missing execution_contract.json")
    if not DEVELOPMENT_SPLIT_PATH.exists():
        fail("missing development_split.json")
    if not CALIBRATION_SPLIT_PATH.exists():
        fail("missing calibration_split.json")
    if not HOLDOUT_SPLIT_PATH.exists():
        fail("missing holdout_split.json")

    precommitment = json.loads(PRECOMMITMENT_PATH.read_text(encoding="utf-8"))
    manifest = json.loads(DATASET_MANIFEST_PATH.read_text(encoding="utf-8"))
    execution = json.loads(EXECUTION_CONTRACT_PATH.read_text(encoding="utf-8"))
    required = {
        "bridge_id": "PB-02",
        "status": "PRECOMMITMENT_DRAFT",
        "dataset": "PowerGraph",
        "holdout_policy": "untouched holdout, no feature changes after inspection",
        "execution_mode": "frozen_external_benchmark_candidate",
    }
    checks = {
        "bridge_id": precommitment.get("bridge_id"),
        "status": precommitment.get("status"),
        "dataset": manifest.get("dataset"),
        "holdout_policy": manifest.get("holdout_policy"),
        "execution_mode": execution.get("execution_mode"),
    }
    mismatches = {key: {"expected": expected, "observed": checks[key]} for key, expected in required.items() if checks.get(key) != expected}
    if mismatches:
        fail(f"mismatched required fields: {json.dumps(mismatches, sort_keys=True)}")

    status = {
        "status": "ok",
        "bridge_id": precommitment["bridge_id"],
        "dataset": manifest["dataset"],
        "dataset_manifest": str(DATASET_MANIFEST_PATH.relative_to(HERE)),
        "precommitment": str(PRECOMMITMENT_PATH.relative_to(HERE)),
        "execution_contract": str(EXECUTION_CONTRACT_PATH.relative_to(HERE)),
        "development_split": str(DEVELOPMENT_SPLIT_PATH.relative_to(HERE)),
        "calibration_split": str(CALIBRATION_SPLIT_PATH.relative_to(HERE)),
        "holdout_split": str(HOLDOUT_SPLIT_PATH.relative_to(HERE)),
        "claim_boundary": precommitment["claim_boundary"],
        "result": "PRECOMMITMENT_DRAFT_ONLY",
    }
    print(json.dumps(status, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
