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


def build_report(bundle: dict[str, object]) -> dict[str, object]:
    audit_name = bundle.get("audit_name")
    if isinstance(audit_name, str):
        summary = dict(bundle.get("summary", {}))
        aggregate = dict(bundle.get("aggregate", {}))
        success = dict(bundle.get("success", {}))
        controls = dict(bundle.get("controls", {}))
        stable_control = dict(controls.get("stable_frozen_sector", {}))
        stable_aggregate = dict(stable_control.get("aggregate", {}))
        return {
            "phase": bundle.get("phase_name", "phase5-readout"),
            "audit": audit_name,
            "dominant_scaling_class": summary.get(
                "stable_dominant_scaling_class",
                aggregate.get("dominant_scaling_class", stable_aggregate.get("dominant_scaling_class")),
            ),
            "class_consistency_fraction": summary.get(
                "stable_class_consistency_fraction",
                aggregate.get("class_consistency_fraction", stable_aggregate.get("class_consistency_fraction")),
            ),
            "overall_classifier_stability_fraction": summary.get(
                "policy_classifier_stability_fraction",
                aggregate.get("overall_classifier_stability_fraction"),
            ),
            "overall_pass": summary.get("overall_pass", success.get("overall_pass")),
        }
    summary = dict(bundle["summary"])
    artifacts = dict(bundle.get("artifacts", {"latest_json": relpath(ROOT / "phase5-readout" / "runs" / "phase5_readout_bundle_latest.json")}))
    return {
        "phase": bundle["phase_name"],
        "bundle": artifacts["latest_json"],
        "input_mode": bundle["input_mode"],
        "control_class": bundle["control_class"],
        "sector_identifier": summary["sector_identifier"],
        "trial_count": summary["trial_count"],
        "mean_recovery_score": summary["mean_recovery_score"],
        "nonzero_recovery_fraction": summary["nonzero_recovery_fraction"],
        "scaling_class": summary["scaling_class"],
        "runtime_seconds": summary["runtime_seconds"],
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
