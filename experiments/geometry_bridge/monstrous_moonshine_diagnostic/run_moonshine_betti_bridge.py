#!/usr/bin/env python3
"""Shared-support bridge manifest for Moonshine and Betti diagnostics.

This runner compares declared perturbation responses from two already bounded
sidecars:

- Moonshine arithmetic support diagnostic.
- Gaussian-prime Betti vector diagnostic.

The bridge type is shared-support diagnostic coupling. It is not a Moonshine
mechanism, Monster action, formal proof, or physical bridge.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from run_betti_vector import BettiVectorConfig, run_betti_vector  # noqa: E402
from run_monstrous_moonshine_diagnostic import (  # noqa: E402
    MoonshineDiagnosticConfig,
    run_diagnostic,
    stable_hash,
    write_csv,
    write_json,
)


MANIFEST_PATH = ROOT / "moonshine_betti_bridge_manifest.json"
COVARIANCE_PATH = ROOT / "moonshine_betti_bridge_covariance.csv"
REPORT_PATH = ROOT / "moonshine_betti_bridge_report.md"
RESULT_PATH = ROOT / "moonshine_betti_bridge_result.json"

COVARIANCE_FIELDS = [
    "bridge_case",
    "moonshine_condition",
    "betti_condition",
    "moonshine_total_distance",
    "moonshine_support_distance",
    "moonshine_exponent_distance",
    "moonshine_decomposition_distance",
    "moonshine_gaussian_class_distance",
    "betti_vector_distance",
    "betti_edge_jaccard_distance",
    "expected_response",
    "observed_status",
]


@dataclass(frozen=True)
class MoonshineBettiBridgeConfig:
    version: str = "moonshine-betti-bridge-v0.1"
    invariance_max: float = 1.0e-12
    moonshine_move_min: float = 0.10
    betti_vector_move_min: float = 1.0
    betti_edge_move_min: float = 0.10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Moonshine-Betti shared-support bridge diagnostic.")
    parser.add_argument("--output-dir", type=Path, default=ROOT)
    return parser.parse_args()


def manifest(config: MoonshineBettiBridgeConfig) -> dict[str, Any]:
    return {
        "bridge_id": "GEO-MM-BETTI-01",
        "version": config.version,
        "shared_support": "15 supersingular primes",
        "moonshine_channel": "arithmetic support diagnostic",
        "betti_channel": "Gaussian-prime threshold graph topology",
        "bridge_type": "shared-support diagnostic coupling",
        "allowed_phrase": "shared-support perturbation covariance",
        "dangerous_phrase": "Moonshine-Betti mechanism",
        "matched_perturbations": {
            "self_recovery": {
                "moonshine_condition": "known_positive",
                "betti_condition": "known_positive",
                "expected_response": "both channels recover with zero distance",
            },
            "shared_support_replacement": {
                "moonshine_condition": "nonsupersingular_prime_control",
                "betti_condition": "support_replacement_control",
                "expected_response": "support replacement moves both arithmetic support telemetry and Betti vector or edge signature",
            },
        },
        "not_claimed": [
            "Moonshine proof",
            "Monster action",
            "Moonshine module construction",
            "physical bridge",
            "continuum limit",
            "quantum claim",
            "gravity claim",
        ],
        "failure_conditions": [
            "self-recovery moves either channel",
            "shared-support replacement moves only one channel",
            "report language upgrades shared support into a derivation",
        ],
        "thresholds": asdict(config),
        "status": "PRECOMMITTED_SHARED_SUPPORT_DIAGNOSTIC_COUPLING",
    }


def moonshine_total_distance(distances: dict[str, float]) -> float:
    return float(sum(float(value) for value in distances.values()))


def row_for_case(
    case_name: str,
    case: dict[str, str],
    moonshine_result: dict[str, Any],
    betti_result: dict[str, Any],
    config: MoonshineBettiBridgeConfig,
) -> dict[str, Any]:
    moonshine = moonshine_result["condition_results"][case["moonshine_condition"]]
    betti = betti_result["condition_results"][case["betti_condition"]]
    moonshine_distances = {name: float(value) for name, value in moonshine["distances"].items()}
    moonshine_total = moonshine_total_distance(moonshine_distances)
    betti_vector_distance = float(betti["betti_vector_distance"])
    betti_edge_distance = float(betti["edge_jaccard_distance"])

    if case_name == "self_recovery":
        status = (
            "PASS"
            if moonshine_total <= config.invariance_max
            and betti_vector_distance <= config.invariance_max
            and betti_edge_distance <= config.invariance_max
            else "FAIL_SELF_RECOVERY"
        )
    elif case_name == "shared_support_replacement":
        moonshine_moves = moonshine_total >= config.moonshine_move_min
        betti_moves = betti_vector_distance >= config.betti_vector_move_min or betti_edge_distance >= config.betti_edge_move_min
        status = "PASS" if moonshine_moves and betti_moves else "FAIL_SHARED_SUPPORT_COVARIANCE"
    else:
        raise ValueError(f"unknown bridge case {case_name}")

    return {
        "bridge_case": case_name,
        "moonshine_condition": case["moonshine_condition"],
        "betti_condition": case["betti_condition"],
        "moonshine_total_distance": f"{moonshine_total:.12g}",
        "moonshine_support_distance": f"{moonshine_distances['support']:.12g}",
        "moonshine_exponent_distance": f"{moonshine_distances['exponent']:.12g}",
        "moonshine_decomposition_distance": f"{moonshine_distances['decomposition']:.12g}",
        "moonshine_gaussian_class_distance": f"{moonshine_distances['gaussian_class']:.12g}",
        "betti_vector_distance": f"{betti_vector_distance:.12g}",
        "betti_edge_jaccard_distance": f"{betti_edge_distance:.12g}",
        "expected_response": case["expected_response"],
        "observed_status": status,
    }


def run_bridge(config: MoonshineBettiBridgeConfig, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    moonshine_result = run_diagnostic(MoonshineDiagnosticConfig(), output_dir)
    betti_result = run_betti_vector(BettiVectorConfig(), output_dir)
    bridge_manifest = manifest(config)
    rows = [
        row_for_case(case_name, case, moonshine_result, betti_result, config)
        for case_name, case in bridge_manifest["matched_perturbations"].items()
    ]

    covariance_pass = all(row["observed_status"] == "PASS" for row in rows)
    labels = [
        "MOONSHINE_BETTI_BRIDGE_MANIFEST_BUILT",
        "SHARED_SUPPORT_DIAGNOSTIC_COUPLING_BUILT",
        "SHARED_SUPPORT_PERTURBATION_COVARIANCE_REPORTED",
        "MOONSHINE_BETTI_MECHANISM_NOT_ESTABLISHED",
        "MOONSHINE_PROOF_NOT_ESTABLISHED",
        "PHYSICAL_BRIDGE_NOT_ESTABLISHED",
        "CLAIM_GATED_SHARED_SUPPORT_DIAGNOSTIC",
    ]
    labels.append("SHARED_SUPPORT_COVARIANCE_PASS" if covariance_pass else "SHARED_SUPPORT_COVARIANCE_OPEN")

    result = {
        "bridge_id": "moonshine_betti_bridge_v0_1",
        "status": "PASS" if covariance_pass else "MIXED_OPEN",
        "classification": "SHARED_SUPPORT_DIAGNOSTIC_COUPLING",
        "labels": labels,
        "manifest": bridge_manifest,
        "input_hashes": {
            "moonshine_result_hash": moonshine_result["result_hash"],
            "betti_vector_result_hash": betti_result["result_hash"],
        },
        "covariance_rows": rows,
        "unmatched_controls": {
            "moonshine_only": [
                "exponent_shuffle_control",
                "decomposition_break_control",
                "gaussian_class_swap_control",
            ],
            "betti_only": [
                "d4_rotate_90",
                "d4_reflect_real",
                "isomorphism_relabel",
                "threshold_mutation_control",
                "topology_destroyed_control",
            ],
        },
        "claim_ceiling": "Shared-support perturbation covariance only; no Moonshine-Betti mechanism, Moonshine proof, or physical bridge claim.",
    }
    result["result_hash"] = stable_hash("moonshine_betti_bridge", result)

    write_json(output_dir / MANIFEST_PATH.name, bridge_manifest)
    write_csv(output_dir / COVARIANCE_PATH.name, rows, COVARIANCE_FIELDS)
    write_json(output_dir / RESULT_PATH.name, result)
    write_report(output_dir / REPORT_PATH.name, result)
    return result


def write_report(path: Path, result: dict[str, Any]) -> None:
    lines = [
        "# Moonshine-Betti Bridge Report",
        "",
        "## Status",
        "",
        f"- Result: `{result['status']}`",
        f"- Classification: `{result['classification']}`",
        f"- Result hash: `{result['result_hash']}`",
        f"- Moonshine input hash: `{result['input_hashes']['moonshine_result_hash']}`",
        f"- Betti vector input hash: `{result['input_hashes']['betti_vector_result_hash']}`",
        "",
        "## Bridge Type",
        "",
        "`shared-support diagnostic coupling`",
        "",
        "Useful phrase: `shared-support perturbation covariance`.",
        "",
        "Dangerous phrase: `Moonshine-Betti mechanism`.",
        "",
        "## Labels",
        "",
        *[f"- `{label}`" for label in result["labels"]],
        "",
        "## Matched Perturbation Response",
        "",
        "| Case | Moonshine condition | Betti condition | Moonshine total | Betti vector | Betti edge | Status |",
        "| --- | --- | --- | ---: | ---: | ---: | --- |",
    ]
    for row in result["covariance_rows"]:
        lines.append(
            "| {bridge_case} | {moonshine_condition} | {betti_condition} | {moonshine_total_distance} | {betti_vector_distance} | {betti_edge_jaccard_distance} | `{observed_status}` |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "This bridge records whether a declared shared-support perturbation moves both bounded diagnostic channels.",
            "It does not derive the Betti graph from Moonshine data, prove Monstrous Moonshine, construct a Monster action, or establish a physical bridge.",
            "",
            "Unmatched channel-specific controls remain visible in `moonshine_betti_bridge_result.json` and must not be compressed into a mechanism claim.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    result = run_bridge(MoonshineBettiBridgeConfig(), args.output_dir)
    print(json.dumps({"status": result["status"], "result_hash": result["result_hash"]}, sort_keys=True))


if __name__ == "__main__":
    main()
