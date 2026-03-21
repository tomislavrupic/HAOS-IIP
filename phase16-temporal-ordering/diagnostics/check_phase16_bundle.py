#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PHASE_ROOT = ROOT / "phase16-temporal-ordering"
MANIFEST_PATH = PHASE_ROOT / "phase16_manifest.json"


def main() -> None:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    required = [
        "phase16-temporal-ordering/runs/phase16_event_ordering_ledger.csv",
        "phase16-temporal-ordering/runs/phase16_monotonic_parameter_ledger.csv",
        "phase16-temporal-ordering/runs/phase16_front_arrival_ordering.csv",
        "phase16-temporal-ordering/runs/phase16_ordering_robustness.csv",
        "phase16-temporal-ordering/runs/phase16_runs.json",
        "phase16-temporal-ordering/phase16_summary.md",
        "phase16-temporal-ordering/phase16_manifest.json",
    ]
    success = bool(manifest.get("success"))
    for rel_path in required:
        success = success and (ROOT / rel_path).exists()
    print(json.dumps({"success": success, "manifest": str(MANIFEST_PATH)}, indent=2))


if __name__ == "__main__":
    main()
