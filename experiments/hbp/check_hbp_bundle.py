#!/usr/bin/env python3
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from experiments.hbp.hbp_validation import stable_hash
from experiments.hbp.run_hbp_registry import HBP_RESULT, REGISTRY_JSON, build_registry


REQUIRED_OUTPUTS = [
    Path("experiments/hbp/results/hbp_bridge_registry.json"),
    Path("experiments/hbp/results/hbp_bridge_registry.csv"),
    Path("experiments/hbp/results/hbp_result.json"),
    Path("experiments/hbp/results/hbp_report.md"),
    Path("experiments/hbp/pb01_network_recovery/results/precommitment_contract.json"),
    Path("experiments/hbp/pb01_network_recovery/results/pb01_result.json"),
    Path("experiments/hbp/pb01_network_recovery/results/pb01_predictions.csv"),
    Path("experiments/hbp/pb01_network_recovery/results/pb01_baseline_metrics.csv"),
    Path("experiments/hbp/pb01_network_recovery/results/pb01_control_results.csv"),
    Path("experiments/hbp/pb02_external_network_recovery/README.md"),
    Path("experiments/hbp/pb02_external_network_recovery/observation_map.md"),
    Path("experiments/hbp/pb02_external_network_recovery/precommitment_contract.json"),
    Path("experiments/hbp/pb02_external_network_recovery/dataset_manifest.json"),
    Path("experiments/hbp/pb02_external_network_recovery/run_pb02_external_network_recovery.py"),
    Path("experiments/hbp/pb03_temporal_recovery/README.md"),
    Path("experiments/hbp/pb03_temporal_recovery/precommitment_contract.json"),
    Path("experiments/hbp/pb04_cross_domain_transfer/README.md"),
    Path("experiments/hbp/pb04_cross_domain_transfer/precommitment_contract.json"),
]


ROOT = Path(__file__).resolve().parents[2]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def fail(message: str) -> None:
    raise SystemExit(f"HBP bundle check failed: {message}")


def main() -> int:
    missing = [str(path) for path in REQUIRED_OUTPUTS if not (ROOT / path).exists()]
    if missing:
        fail(f"missing outputs: {missing}")

    disk_result = read_json(HBP_RESULT)
    comparable = dict(disk_result)
    comparable.pop("result_hash", None)
    if disk_result.get("result_hash") != stable_hash(comparable, "hbp_result_"):
        fail("hbp_result hash mismatch")

    memory = build_registry(write_outputs=False)
    if memory["hbp_result"]["result_hash"] != disk_result.get("result_hash"):
        fail("fresh in-memory HBP result differs from disk")

    registry = read_json(REGISTRY_JSON)
    if not registry.get("registry_rows"):
        fail("registry has no rows")
    if any(row["bridge_id"] == "BELL-BRIDGE-V0-3-TERMINAL" and row["effective_classification"] != "FORMAL_CORRESPONDENCE" for row in registry["registry_rows"]):
        fail("terminal Bell bridge was promoted")
    if any(row["bridge_id"] == "SYNTHETIC-RELATIONAL-CALIBRATION" and row["effective_classification"] not in {"OPERATIONAL_MAPPING", "FORMAL_CORRESPONDENCE"} for row in registry["registry_rows"]):
        fail("synthetic calibration was promoted too far")

    print(
        json.dumps(
            {
                "status": "ok",
                "bridge_count": disk_result["bridge_count"],
                "pb01_status": disk_result["pb01_status"],
                "result_hash": disk_result["result_hash"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
