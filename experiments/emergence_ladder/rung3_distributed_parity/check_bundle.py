#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from experiments.emergence_ladder.rung3_distributed_parity.run_experiment import (
    CONTRACT_PATH,
    FINAL_DIR,
    ROOT,
    file_hash,
    stable_hash,
)


REQUIRED = (
    ROOT / "README.md",
    ROOT / "RELATIONAL_ALPHABET.md",
    ROOT / "PARITY_ARCHITECTURE_SELECTION.md",
    ROOT / "FUNCTIONAL_TARGET.md",
    ROOT / "THEORY_INTERPRETATION.md",
    ROOT / "CLASSIFICATION.md",
    CONTRACT_PATH,
    ROOT / "information_accounting.json",
    ROOT / "calibration/calibration_runs.csv",
    ROOT / "calibration/parameter_selection.json",
    ROOT / "calibration/derived_thresholds.json",
    ROOT / "validation/validation_runs.csv",
    ROOT / "validation/validation_result.json",
    ROOT / "validation/uncertainty.json",
    ROOT / "validation/control_results.csv",
    ROOT / "validation/radius_results.csv",
    ROOT / "validation/perturbation_results.csv",
)


def check() -> dict[str, object]:
    missing = [str(path.relative_to(ROOT)) for path in REQUIRED if not path.exists()]
    if missing:
        raise ValueError(f"missing RP-01 artifacts: {missing}")
    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    result = json.loads((ROOT / "validation/validation_result.json").read_text(encoding="utf-8"))
    comparable = dict(result)
    observed_hash = comparable.pop("result_hash")
    if observed_hash != stable_hash("el_r3_rp_01", comparable):
        raise ValueError("RP-01 result hash mismatch")
    if result["contract_hash"] != stable_hash("el_r3_rp_01_contract", contract):
        raise ValueError("RP-01 contract hash mismatch")
    expected_sources = {
        name: file_hash(ROOT / filename) for name, filename in {
            "alphabet": "alphabet.py", "parity_architecture": "parity_architecture.py",
            "decoder": "decoder.py", "bridge": "bridge.py", "controls": "controls.py", "runner": "run_experiment.py",
        }.items()
    }
    if result["source_hashes"] != expected_sources:
        raise ValueError("RP-01 source hash mismatch")
    if result["classification"] != "VALIDATION_GATE_FAILED":
        raise ValueError("unexpected RP-01 classification")
    if result["final_evaluation_authorized"] or result["final_seed_count_consumed"] != 0:
        raise ValueError("final evaluation access policy violated")
    if FINAL_DIR.exists() and any(FINAL_DIR.iterdir()):
        raise ValueError("final artifacts exist after failed validation gate")
    with (ROOT / "validation/validation_runs.csv").open(newline="", encoding="utf-8") as handle:
        seeds = {int(row["seed"]) for row in csv.DictReader(handle)}
    if seeds & set(contract["final_evaluation"]["seeds"]):
        raise ValueError("final seed appeared in validation artifacts")
    if result["historical_rt02_result_hash"] != "el_r3_rt_02_cade4e057a54d119e6af143c":
        raise ValueError("frozen RT-02 baseline hash changed")
    if not all(result["validation_gate"]["leakage_guard"].values()):
        raise ValueError("leakage guard failed")
    return {
        "status": "ok",
        "classification": result["classification"],
        "result_hash": observed_hash,
        "final_seeds_consumed": 0,
    }


def main() -> int:
    print(json.dumps(check(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

