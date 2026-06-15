#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Any

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from experiments.hbp.hbp_validation import assess_contracts, contract_from_dict, stable_hash
from experiments.hbp.pb01_network_recovery.run_pb01_network_recovery import run_pb01


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
CONTRACTS_PATH = HERE / "bridge_contracts.json"
RESULTS_DIR = HERE / "results"
REGISTRY_JSON = RESULTS_DIR / "hbp_bridge_registry.json"
REGISTRY_CSV = RESULTS_DIR / "hbp_bridge_registry.csv"
HBP_RESULT = RESULTS_DIR / "hbp_result.json"
HBP_REPORT = RESULTS_DIR / "hbp_report.md"


def repo_rel(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def artifact_status(paths: list[str]) -> dict[str, bool]:
    return {path: (REPO_ROOT / path).exists() for path in paths}


def build_registry(write_outputs: bool = True) -> dict[str, Any]:
    pb01_payload = run_pb01(write_outputs=write_outputs)["result"]
    raw = read_json(CONTRACTS_PATH)
    contracts = [contract_from_dict(item) for item in raw["contracts"]]
    assessments = assess_contracts(contracts)
    assessment_by_id = {item.bridge_id: item for item in assessments}

    rows: list[dict[str, Any]] = []
    for contract in contracts:
        assessment = assessment_by_id[contract.bridge_id]
        artifacts = artifact_status(contract.provenance.source_artifacts)
        rows.append(
            {
                "bridge_id": contract.bridge_id,
                "domain": contract.domain,
                "requested_classification": assessment.requested_classification,
                "effective_classification": assessment.effective_classification,
                "verdict_authorized": assessment.verdict_authorized,
                "labels": ";".join(assessment.labels),
                "missing_required_fields": ";".join(assessment.missing_required_fields),
                "ceiling_applied": assessment.ceiling_applied,
                "source_artifacts_exist": all(artifacts.values()) if artifacts else False,
                "assessment_hash": assessment.result_hash,
            }
        )

    counts: dict[str, int] = {}
    for row in rows:
        key = str(row["effective_classification"])
        counts[key] = counts.get(key, 0) + 1

    registry_payload = {
        "schema_version": raw["schema_version"],
        "scope": raw["scope"],
        "contracts_hash": stable_hash(raw, "hbp_contracts_"),
        "registry_rows": rows,
        "assessments": [item.to_dict() for item in assessments],
        "counts_by_effective_classification": dict(sorted(counts.items())),
        "pb01_result_hash": pb01_payload["result_hash"],
        "labels": ["REGISTRY_VALID", "PHYSICAL_MECHANISM_NOT_ESTABLISHED"],
    }
    registry_payload["registry_hash"] = stable_hash(registry_payload, "hbp_registry_")

    hbp_result = {
        "status": "REGISTRY_VALID",
        "registry_hash": registry_payload["registry_hash"],
        "contracts_hash": registry_payload["contracts_hash"],
        "bridge_count": len(rows),
        "counts_by_effective_classification": registry_payload["counts_by_effective_classification"],
        "pb01_status": pb01_payload["verdict"]["status"],
        "pb01_labels": pb01_payload["verdict"]["labels"],
        "non_claims": [
            "no empirical physical validation",
            "no physical mechanism claim",
            "no Bell, hydrogen, or semiconductor derivation",
            "synthetic PB-01 calibration is not external reality correspondence",
        ],
        "outputs": {
            "registry_json": repo_rel(REGISTRY_JSON),
            "registry_csv": repo_rel(REGISTRY_CSV),
            "hbp_result": repo_rel(HBP_RESULT),
            "hbp_report": repo_rel(HBP_REPORT),
        },
    }
    hbp_result["result_hash"] = stable_hash(hbp_result, "hbp_result_")

    if write_outputs:
        write_json(REGISTRY_JSON, registry_payload)
        write_csv(
            REGISTRY_CSV,
            rows,
            [
                "bridge_id",
                "domain",
                "requested_classification",
                "effective_classification",
                "verdict_authorized",
                "labels",
                "missing_required_fields",
                "ceiling_applied",
                "source_artifacts_exist",
                "assessment_hash",
            ],
        )
        write_json(HBP_RESULT, hbp_result)
        write_report(hbp_result, rows)

    return {
        "registry_payload": registry_payload,
        "hbp_result": hbp_result,
        "pb01_result": pb01_payload,
        "contracts": contracts,
        "assessments": assessments,
    }


def write_report(result: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    lines = [
        "# HBP Registry Report",
        "",
        f"Status: `{result['status']}`",
        "",
        "This registry classifies bridge evidence. It does not upgrade frozen HAOS-IIP claims.",
        "",
        "## Classification Counts",
        "",
    ]
    for key, value in result["counts_by_effective_classification"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## PB-01", ""])
    lines.append(f"- status: `{result['pb01_status']}`")
    lines.extend(f"- `{label}`" for label in result["pb01_labels"])
    lines.extend(["", "## Registered Bridges", "", "| bridge_id | effective classification | labels |", "| --- | --- | --- |"])
    for row in rows:
        lines.append(f"| {row['bridge_id']} | {row['effective_classification']} | {row['labels']} |")
    lines.extend(["", "## Non-Claims", ""])
    lines.extend(f"- {item}" for item in result["non_claims"])
    HBP_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    payload = build_registry(write_outputs=True)["hbp_result"]
    print(json.dumps({"status": payload["status"], "result_hash": payload["result_hash"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
