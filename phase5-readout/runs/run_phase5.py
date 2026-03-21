#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from haos_core import read_json, relpath, write_json, write_timestamped_json


def unwrap_phase_payload(payload: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if "summary" in payload and isinstance(payload["summary"], dict):
        summary = payload["summary"]
        claims = list(payload.get("claim_ledger", summary.get("claim_ledger", [])))
        return summary, claims
    return payload, list(payload.get("claim_ledger", []))


def resolve_candidate_path(explicit_path: Path | None, candidates: list[str]) -> Path:
    if explicit_path is not None:
        return explicit_path
    for candidate in candidates:
        candidate_path = ROOT / candidate
        if candidate_path.exists():
            return candidate_path
    raise FileNotFoundError(f"no frozen input found in candidates: {candidates}")


def normalize_phase3(summary: dict[str, Any], claims: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "experiment": summary.get("experiment", "phase3_unknown"),
        "freeze_valid": bool(summary.get("phaseIII_freeze_valid", summary.get("phase3_freeze_valid", False))),
        "minimal_read": summary.get("phaseIII_minimal_read", "phase3 read unavailable"),
        "claim_ledger": claims,
    }


def normalize_phase4(summary: dict[str, Any], claims: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "experiment": summary.get("experiment", "phase4_unknown"),
        "freeze_valid": bool(summary.get("phaseIV_freeze_valid", summary.get("phase4_freeze_valid", False))),
        "minimal_read": summary.get("phaseIV_minimal_read", "phase4 read unavailable"),
        "selector_threshold_fixed": summary.get("selector_threshold_fixed"),
        "extension_status_counts": dict(summary.get("extension_status_counts", {})),
        "claim_ledger": claims,
    }


def load_dummy_inputs(config: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], dict[str, str]]:
    dummy = dict(config["dummy_inputs"])
    phase3 = normalize_phase3(dummy["phase3"], list(dummy["phase3"]["claim_ledger"]))
    phase4 = normalize_phase4(dummy["phase4"], list(dummy["phase4"]["claim_ledger"]))
    refs = {
        "phase3_input": "config:dummy_inputs.phase3",
        "phase4_input": "config:dummy_inputs.phase4",
    }
    return phase3, phase4, refs


def load_frozen_inputs(
    config: dict[str, Any],
    phase3_input: Path | None,
    phase4_input: Path | None,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, str]]:
    frozen = dict(config["frozen_inputs"])
    phase3_path = resolve_candidate_path(phase3_input, list(frozen["phase3_candidates"]))
    phase4_path = resolve_candidate_path(phase4_input, list(frozen["phase4_candidates"]))
    phase3_payload = read_json(phase3_path)
    phase4_payload = read_json(phase4_path)
    phase3_summary, phase3_claims = unwrap_phase_payload(phase3_payload)
    phase4_summary, phase4_claims = unwrap_phase_payload(phase4_payload)
    refs = {
        "phase3_input": relpath(phase3_path),
        "phase4_input": relpath(phase4_path),
    }
    return normalize_phase3(phase3_summary, phase3_claims), normalize_phase4(phase4_summary, phase4_claims), refs


def build_recovery_histogram(phase3: dict[str, Any], phase4: dict[str, Any]) -> dict[str, Any]:
    phase3_counts = Counter(str(row.get("status", "unknown")) for row in phase3["claim_ledger"])
    phase4_counts = Counter(str(row.get("status", "unknown")) for row in phase4["claim_ledger"])
    combined = Counter(phase3_counts) + Counter(phase4_counts)
    return {
        "phase3": dict(sorted(phase3_counts.items())),
        "phase4": dict(sorted(phase4_counts.items())),
        "combined": dict(sorted(combined.items())),
    }


def build_summary(mode: str, phase3: dict[str, Any], phase4: dict[str, Any], histogram: dict[str, Any]) -> dict[str, Any]:
    return {
        "experiment": "5_readout_scaffold",
        "description": "Phase V scaffold over frozen Phase III and Phase IV outputs.",
        "input_mode": mode,
        "phase3_freeze_valid": phase3["freeze_valid"],
        "phase4_freeze_valid": phase4["freeze_valid"],
        "recovery_ready": bool(phase3["freeze_valid"] and phase4["freeze_valid"]),
        "phase3_minimal_read": phase3["minimal_read"],
        "phase4_minimal_read": phase4["minimal_read"],
        "selector_threshold_fixed": phase4["selector_threshold_fixed"],
        "extension_status_counts": phase4["extension_status_counts"],
        "combined_status_histogram": histogram["combined"],
    }


def run_phase5(
    config_path: Path,
    output_dir: Path,
    mode: str,
    phase3_input: Path | None,
    phase4_input: Path | None,
) -> dict[str, Any]:
    config = read_json(config_path)
    if mode == "dummy":
        phase3, phase4, input_refs = load_dummy_inputs(config)
    else:
        phase3, phase4, input_refs = load_frozen_inputs(config, phase3_input, phase4_input)

    histogram = build_recovery_histogram(phase3, phase4)
    histogram_payload = {
        "phase": 5,
        "phase_name": "phase5-readout",
        "input_mode": mode,
        "inputs": input_refs,
        "histogram": histogram,
    }
    histogram_stamped, histogram_latest = write_timestamped_json(output_dir, "phase5_recovery_histogram", histogram_payload)
    histogram_payload["artifacts"] = {
        "stamped_json": relpath(histogram_stamped),
        "latest_json": relpath(histogram_latest),
    }
    write_json(histogram_stamped, histogram_payload)
    write_json(histogram_latest, histogram_payload)

    bundle = {
        "phase": 5,
        "phase_name": "phase5-readout",
        "source_config": relpath(config_path),
        "input_mode": mode,
        "inputs": input_refs,
        "summary": build_summary(mode, phase3, phase4, histogram),
        "readout": {
            "phase3": phase3,
            "phase4": phase4,
        },
    }
    bundle_stamped, bundle_latest = write_timestamped_json(output_dir, "phase5_readout_bundle", bundle)
    bundle["artifacts"] = {
        "stamped_json": relpath(bundle_stamped),
        "latest_json": relpath(bundle_latest),
        "recovery_histogram_stamped_json": relpath(histogram_stamped),
        "recovery_histogram_latest_json": relpath(histogram_latest),
    }
    write_json(bundle_stamped, bundle)
    write_json(bundle_latest, bundle)
    return bundle


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the Phase V readout scaffold.")
    parser.add_argument(
        "--config",
        type=Path,
        default=ROOT / "phase5-readout" / "configs" / "phase5_readout_scaffold.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "phase5-readout" / "runs",
    )
    parser.add_argument(
        "--mode",
        choices=["dummy", "frozen"],
        default="dummy",
    )
    parser.add_argument(
        "--phase3-input",
        type=Path,
        default=None,
        help="Optional explicit Phase III input bundle/JSON path for frozen mode.",
    )
    parser.add_argument(
        "--phase4-input",
        type=Path,
        default=None,
        help="Optional explicit Phase IV input bundle/JSON path for frozen mode.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print the final bundle as JSON.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    bundle = run_phase5(args.config, args.output_dir, args.mode, args.phase3_input, args.phase4_input)
    if args.json:
        print(json.dumps(bundle, indent=2))
        return
    print(f"Phase V recovery ready: {bundle['summary']['recovery_ready']}")
    print(f"Latest bundle: {bundle['artifacts']['latest_json']}")


if __name__ == "__main__":
    main()
