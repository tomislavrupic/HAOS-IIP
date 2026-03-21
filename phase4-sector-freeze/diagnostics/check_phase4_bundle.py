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
        "phaseIV_freeze_valid": summary["phaseIV_freeze_valid"],
        "selector_threshold_fixed": summary["selector_threshold_fixed"],
        "first_selector_failure_level": summary["first_selector_failure_level"],
        "portable_claims": sum(1 for row in bundle["claim_ledger"] if row["claim_id"].startswith("P4_POS")),
        "negative_claims": sum(1 for row in bundle["claim_ledger"] if row["claim_id"].startswith("P4_NEG")),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Check a Phase IV sector freeze bundle.")
    parser.add_argument(
        "--bundle",
        type=Path,
        default=ROOT / "phase4-sector-freeze" / "runs" / "phase4_sector_freeze_bundle_latest.json",
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
