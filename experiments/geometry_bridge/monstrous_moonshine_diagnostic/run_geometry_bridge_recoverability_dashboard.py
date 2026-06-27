#!/usr/bin/env python3
"""Composite recoverability dashboard for GEO-MM-01.

This runner aggregates the bounded Moonshine arithmetic, Betti, robustness,
null, source-validation, bridge-covariance, and formal-target channels. It is a
dashboard over already declared diagnostics, not a new mechanism.
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

from arithmetic_source_validation import SourceValidationConfig, run_source_validation  # noqa: E402
from run_betti_vector import BettiVectorConfig, run_betti_vector  # noqa: E402
from run_formal_lean_targets import FormalLeanTargetConfig, run_formal_targets  # noqa: E402
from run_monstrous_moonshine_diagnostic import MoonshineDiagnosticConfig, run_diagnostic, stable_hash, write_json  # noqa: E402
from run_moonshine_betti_bridge import MoonshineBettiBridgeConfig, run_bridge  # noqa: E402
from threshold_sweep_betti_stability import BettiSweepConfig, run_threshold_sweep  # noqa: E402


REPORT_PATH = ROOT / "geometry_bridge_recoverability_report.md"
RESULT_PATH = ROOT / "geometry_bridge_recoverability_result.json"
FAILURE_LEDGER_PATH = ROOT / "failure_conditions.md"


@dataclass(frozen=True)
class GeometryBridgeDashboardConfig:
    version: str = "geometry-bridge-recoverability-dashboard-v0.1"
    null_false_positive_watch: float = 0.05
    claim_ceiling: str = "COMPOSITE_DIAGNOSTIC_DASHBOARD_ONLY"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run GEO-MM-01 composite recoverability dashboard.")
    parser.add_argument("--output-dir", type=Path, default=ROOT)
    return parser.parse_args()


def channel(name: str, status: str, classification: str, evidence: dict[str, Any]) -> dict[str, Any]:
    return {
        "name": name,
        "status": status,
        "classification": classification,
        "evidence": evidence,
    }


def claim_boundary_pass(results: list[dict[str, Any]]) -> bool:
    forbidden = {
        "MOONSHINE_PROOF_ESTABLISHED",
        "PHYSICAL_BRIDGE_ESTABLISHED",
        "LEAN_THEOREM_ESTABLISHED",
        "MONSTER_ACTION_ESTABLISHED",
    }
    for result in results:
        labels = set(result.get("labels", []))
        if labels & forbidden:
            return False
    return True


def run_dashboard(config: GeometryBridgeDashboardConfig, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    moonshine = run_diagnostic(MoonshineDiagnosticConfig(), output_dir)
    betti = run_betti_vector(BettiVectorConfig(), output_dir)
    sweep = run_threshold_sweep(BettiSweepConfig(), output_dir)
    bridge = run_bridge(MoonshineBettiBridgeConfig(), output_dir)
    source_validation = run_source_validation(SourceValidationConfig(), output_dir)
    formal_targets = run_formal_targets(FormalLeanTargetConfig(), output_dir)

    null_rate = float(sweep["null_summary"]["false_positive_rate"])
    null_status = "OPEN" if null_rate > config.null_false_positive_watch else "PASS"
    formal_status = "OPEN" if formal_targets["status"] == "OPEN" else formal_targets["status"]
    boundary_pass = claim_boundary_pass([moonshine, betti, sweep, bridge, source_validation, formal_targets])

    channels = [
        channel(
            "Moonshine arithmetic diagnostic",
            moonshine["status"],
            moonshine["classification"],
            {"result_hash": moonshine["result_hash"]},
        ),
        channel(
            "Gaussian-prime norm-lift support",
            "PARTIAL",
            "REPRESENTATIVE_INPUT_TO_BETTI_GRAPH",
            {
                "reason": "Gaussian representatives and residue classes are validated as inputs, but no standalone norm-lift recovery runner exists in this bundle.",
                "source_validation_hash": source_validation["result_hash"],
            },
        ),
        channel(
            "Betti_0 / Betti_1 diagnostic",
            betti["status"],
            betti["classification"],
            {"result_hash": betti["result_hash"], "reference_signature": betti["reference_signature"]},
        ),
        channel(
            "Threshold robustness",
            "PASS_WITH_FRAGILITY",
            "LOCAL_ROBUSTNESS_BAND",
            {
                "result_hash": sweep["result_hash"],
                "betti0_stability_band": sweep["sweep_summary"]["betti0_stability_band"],
                "edge_exact_band": sweep["sweep_summary"]["edge_exact_band"],
            },
        ),
        channel(
            "Null ensemble rarity",
            null_status,
            "FALSE_POSITIVE_RATE_ESTIMATE",
            {
                "false_positive_rate": null_rate,
                "watch_threshold": config.null_false_positive_watch,
                "result_hash": sweep["result_hash"],
            },
        ),
        channel(
            "Cross-channel perturbation covariance",
            bridge["status"],
            bridge["classification"],
            {"result_hash": bridge["result_hash"], "covariance_rows": bridge["covariance_rows"]},
        ),
        channel(
            "Pinned source validation",
            source_validation["status"],
            source_validation["classification"],
            {"result_hash": source_validation["result_hash"], "manifest_hash": source_validation["manifest_hash"]},
        ),
        channel(
            "Formal Lean targets",
            formal_status,
            formal_targets["classification"],
            {
                "result_hash": formal_targets["result_hash"],
                "lean_project_present": formal_targets["config"]["lean_project_present"],
                "lean_check_run": formal_targets["config"]["lean_check_run"],
            },
        ),
        channel(
            "Failure ledger",
            "PASS" if FAILURE_LEDGER_PATH.exists() else "FAIL",
            "RECOVERABILITY_FAILURE_CONTRACT",
            {"path": FAILURE_LEDGER_PATH.name},
        ),
        channel(
            "Claim-boundary checks",
            "PASS" if boundary_pass else "FAIL",
            "NO_FORBIDDEN_UPGRADE_LABELS",
            {"forbidden_upgrade_labels_absent": boundary_pass},
        ),
    ]

    if any(item["status"] == "FAIL" for item in channels):
        status = "FAIL"
    elif any(item["status"] in {"OPEN", "PARTIAL", "PASS_WITH_FRAGILITY"} for item in channels):
        status = "PASS_WITH_FRAGILITY"
    else:
        status = "PASS"

    labels = [
        "GEOMETRY_BRIDGE_RECOVERABILITY_DASHBOARD_BUILT",
        "MOONSHINE_ARITHMETIC_CHANNEL_INCLUDED",
        "GAUSSIAN_PRIME_SUPPORT_CHANNEL_PARTIAL",
        "BETTI_VECTOR_CHANNEL_INCLUDED",
        "THRESHOLD_ROBUSTNESS_INCLUDED",
        "NULL_RARITY_INCLUDED",
        "CROSS_CHANNEL_COVARIANCE_INCLUDED",
        "FORMAL_TARGET_STATUS_INCLUDED",
        "CLAIM_BOUNDARY_CHECKS_PASS" if boundary_pass else "CLAIM_BOUNDARY_CHECKS_FAIL",
        "GLOBAL_CLAIM_OPEN",
        "PHYSICAL_BRIDGE_NOT_ESTABLISHED",
    ]
    if status == "PASS_WITH_FRAGILITY":
        labels.append("PASS_WITH_FRAGILITY")

    result = {
        "bridge_id": "geometry_bridge_recoverability_dashboard_v0_1",
        "status": status,
        "classification": "COMPOSITE_DIAGNOSTIC_DASHBOARD_ONLY",
        "labels": labels,
        "config": asdict(config),
        "channels": channels,
        "best_case_classification": [
            "LOCAL_SIGNAL_STABLE",
            "MATHEMATICALLY_INTERESTING",
            "ENGINEERING_FEASIBLE",
            "FORMALLY_PARTIAL",
            "GLOBAL_CLAIM_OPEN",
        ],
        "claim_ceiling": "Composite recoverability dashboard only; no Moonshine proof, Monster action, physical bridge, continuum, quantum, or gravity claim.",
    }
    result["result_hash"] = stable_hash("geometry_bridge_dashboard", result)

    write_json(output_dir / RESULT_PATH.name, result)
    write_report(output_dir / REPORT_PATH.name, result)
    return result


def write_report(path: Path, result: dict[str, Any]) -> None:
    lines = [
        "# Geometry Bridge Recoverability Dashboard",
        "",
        "## Status",
        "",
        f"- Result: `{result['status']}`",
        f"- Classification: `{result['classification']}`",
        f"- Result hash: `{result['result_hash']}`",
        "",
        "## Labels",
        "",
        *[f"- `{label}`" for label in result["labels"]],
        "",
        "## Channels",
        "",
        "| Channel | Status | Classification |",
        "| --- | --- | --- |",
    ]
    for item in result["channels"]:
        lines.append(f"| {item['name']} | `{item['status']}` | `{item['classification']}` |")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "`PASS_WITH_FRAGILITY` means the local executable diagnostics recover declared structure, but null rarity, standalone norm-lift recovery, and local Lean certification remain incomplete or open.",
            "",
            "## Boundary",
            "",
            "This dashboard aggregates bounded diagnostics only. It does not prove Monstrous Moonshine, recover a Monster action, establish a physical bridge, or derive continuum, quantum, or gravity structure.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    result = run_dashboard(GeometryBridgeDashboardConfig(), args.output_dir)
    print(json.dumps({"status": result["status"], "result_hash": result["result_hash"]}, sort_keys=True))


if __name__ == "__main__":
    main()
