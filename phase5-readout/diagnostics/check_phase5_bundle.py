#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haos_core import read_json, relpath


def build_report(bundle: dict) -> dict:
    summary = bundle["summary"]
    return {
        "phase": bundle["phase_name"],
        "bundle": bundle["artifacts"]["latest_json"],
        "input_mode": bundle["input_mode"],
        "recovery_ready": summary["recovery_ready"],
        "phase3_freeze_valid": summary["phase3_freeze_valid"],
        "phase4_freeze_valid": summary["phase4_freeze_valid"],
        "combined_status_histogram": summary["combined_status_histogram"],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Check a Phase V readout bundle.")
    parser.add_argument(
        "--bundle",
        type=Path,
        default=ROOT / "phase5-readout" / "runs" / "phase5_readout_bundle_latest.json",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    bundle = read_json(args.bundle)
    if "artifacts" not in bundle:
        bundle["artifacts"] = {"latest_json": relpath(args.bundle)}
    print(json.dumps(build_report(bundle), indent=2))


if __name__ == "__main__":
    main()
