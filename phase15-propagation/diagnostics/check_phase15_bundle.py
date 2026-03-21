#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PHASE_ROOT = ROOT / "phase15-propagation"
MANIFEST_PATH = PHASE_ROOT / "phase15_manifest.json"


def main() -> None:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    required = [
        "phase15-propagation/runs/phase15_propagation_ledger.csv",
        "phase15-propagation/runs/phase15_effective_speed_ledger.csv",
        "phase15-propagation/runs/phase15_influence_range_ledger.csv",
        "phase15-propagation/runs/phase15_transport_descriptor_ledger.csv",
        "phase15-propagation/runs/phase15_runs.json",
        "phase15-propagation/phase15_summary.md",
        "phase15-propagation/phase15_manifest.json",
    ]
    success = bool(manifest.get("success"))
    for rel_path in required:
        success = success and (ROOT / rel_path).exists()
    print(json.dumps({"success": success, "manifest": str(MANIFEST_PATH)}, indent=2))


if __name__ == "__main__":
    main()
