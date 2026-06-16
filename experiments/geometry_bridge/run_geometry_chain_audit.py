#!/usr/bin/env python3
"""Frozen audit for the geometry -> semantics -> probability -> observable chain.

This verifier does not fit or tune anything. It only reads the frozen GEO-01
through GEO-04 artifacts and reports where the chain is established, partial,
or still open.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
REPO_ROOT = Path(__file__).resolve().parents[2]

GEO01 = ROOT / "synthetic_intrinsic_geometry_recovery" / "geometry_recovery_result.json"
GEO02 = ROOT / "synthetic_transformation_semantics_recovery" / "semantics_result.json"
GEO03 = ROOT / "synthetic_probability_rule_recovery" / "probability_rule_result.json"
GEO_P1 = ROOT / "synthetic_hidden_probability_rule_recovery" / "probability_rule_result.json"
GEO04 = ROOT / "synthetic_observable_prediction" / "observable_prediction_result.json"

REPORT_PATH = ROOT / "geometry_chain_audit_report.md"
RESULT_PATH = ROOT / "geometry_chain_audit_result.json"

def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def repo_rel(path: Path) -> str:
    resolved = path.resolve()
    root = str(REPO_ROOT)
    text = str(resolved)
    if text.startswith(root + "/"):
        return text[len(root) + 1 :]
    return str(path)


def evaluate() -> dict[str, Any]:
    geo01 = read_json(GEO01)
    geo02 = read_json(GEO02)
    geo03 = read_json(GEO03)
    geo_p1 = read_json(GEO_P1) if GEO_P1.exists() else None
    geo04 = read_json(GEO04)

    probability_result = geo_p1 or geo03
    probability_name = "hidden_probability_rule" if geo_p1 else "probability_rule"

    statuses = [
        {
            "name": "intrinsic_geometry",
            "status": "VALID" if geo01.get("holdout_pass") else "OPEN",
            "synthetic_status": "VALID" if geo01.get("holdout_pass") else "OPEN",
            "external_status": "OPEN",
            "unresolved_criterion": "holdout transfer does not outperform the best baseline on the spiral family",
            "evidence": f"holdout_pass={geo01.get('holdout_pass')} / best_baseline_spearman={geo01.get('holdout_metrics', {}).get('best_baseline_spearman')}",
        },
        {
            "name": "transformation_semantics",
            "status": "VALID" if geo02.get("labels") else "OPEN",
            "synthetic_status": "VALID" if geo02.get("labels") else "OPEN",
            "external_status": "OPEN",
            "evidence": "synthetic semantics layer present; frozen control summary exists but no Bell relevance",
        },
        {
            "name": probability_name,
            "status": "VALID" if probability_result.get("holdout_pass") else "PARTIAL",
            "synthetic_status": "VALID" if probability_result.get("holdout_pass") else "PARTIAL",
            "external_status": "OPEN",
            "evidence": f"holdout_accuracy={probability_result.get('holdout_metrics', {}).get('accuracy')} / null_accuracy={probability_result.get('holdout_metrics', {}).get('null_accuracy')} / brier={probability_result.get('holdout_metrics', {}).get('brier_score')}",
        },
        {
            "name": "observable_prediction",
            "status": "VALID"
            if geo04.get("holdout_metrics", {}).get("accuracy", 0.0) > geo04.get("holdout_metrics", {}).get("null_accuracy", 0.0)
            and geo04.get("holdout_metrics", {}).get("pairwise_accuracy", 0.0) > geo04.get("holdout_metrics", {}).get("pairwise_null_accuracy", 0.0)
            else "PARTIAL"
            if geo04.get("holdout_metrics", {}).get("pairwise_accuracy", 0.0) > geo04.get("holdout_metrics", {}).get("pairwise_null_accuracy", 0.0)
            else "OPEN",
            "synthetic_status": "VALID" if geo04.get("holdout_metrics", {}).get("accuracy", 0.0) > geo04.get("holdout_metrics", {}).get("null_accuracy", 0.0) else "PARTIAL",
            "external_status": "OPEN",
            "evidence": f"pairwise_accuracy={geo04.get('holdout_metrics', {}).get('pairwise_accuracy')} / pairwise_null_accuracy={geo04.get('holdout_metrics', {}).get('pairwise_null_accuracy')} / accuracy={geo04.get('holdout_metrics', {}).get('accuracy')} / null_accuracy={geo04.get('holdout_metrics', {}).get('null_accuracy')}",
        },
    ]

    result = {
        "bridge_id": "GEO-05",
        "status": "OPEN",
        "chain": [
            item for item in statuses
        ],
        "labels": ["PREGEOMETRIC_CHAIN_OPEN", "MIXED_OPEN"],
        "status_axes": {
            "synthetic_calibration": "CLOSED",
            "external_semantics": "OPEN",
        },
        "boundary": "Synthetic calibration only; no Bell scoring; no physical mechanism claim.",
        "source_artifacts": [repo_rel(GEO01), repo_rel(GEO02), repo_rel(GEO03), repo_rel(GEO_P1) if GEO_P1.exists() else None, repo_rel(GEO04)],
    }
    result["source_artifacts"] = [path for path in result["source_artifacts"] if path is not None]
    RESULT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    REPORT_PATH.write_text(
        "\n".join(
            [
                "# Geometry Chain Audit",
                "",
                f"- bridge id: {result['bridge_id']}",
                f"- status: {result['status']}",
                "",
                "## Chain Readout",
                *[
                    (
                        f"- {item['name']}: synthetic={item.get('synthetic_status', item['status'])} / "
                        f"external={item.get('external_status', 'OPEN')}"
                        + (
                            f" / unresolved={item['unresolved_criterion']}"
                            if item.get("unresolved_criterion")
                            else ""
                        )
                        + f" ({item['evidence']})"
                    )
                    for item in result["chain"]
                ],
                "",
                "## Status Axes",
                f"- synthetic_calibration: {result['status_axes']['synthetic_calibration']}",
                f"- external_semantics: {result['status_axes']['external_semantics']}",
                "",
                "## Boundary",
                result["boundary"],
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.parse_args()
    evaluate()


if __name__ == "__main__":
    main()
