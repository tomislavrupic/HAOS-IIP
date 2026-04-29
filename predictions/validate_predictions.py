"""Validate the HAOS-IIP prediction ledger.

This script intentionally uses only the Python standard library. It checks that
prediction records are falsifiable, scoped, and explicit about non-claims.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
LEDGER_PATH = ROOT / "prediction_ledger.json"
OUTPUT_DIR = ROOT / "outputs"

ALLOWED_LEVELS = {
    "P0_toy_internal",
    "P1_cross_system",
    "P2_external_dataset",
    "P3_field_facing",
    "P4_physics_candidate",
}

ALLOWED_STATUSES = {
    "candidate_unvalidated",
    "toy_supported_not_external",
    "external_test_ready",
    "externally_supported",
    "failed",
    "retired",
}

REQUIRED_TOP_LEVEL = [
    "schema_version",
    "scope",
    "prediction_levels",
    "records",
]

REQUIRED_RECORD_FIELDS = [
    "id",
    "title",
    "prediction_level",
    "status",
    "domain",
    "claim_type",
    "prediction_statement",
    "novelty_boundary",
    "system_scope",
    "perturbation",
    "observables",
    "pass_criteria",
    "fail_criteria",
    "falsifiers",
    "current_evidence",
    "next_external_test",
    "limitations",
    "non_claims",
]

LIST_FIELDS = {
    "system_scope",
    "observables",
    "pass_criteria",
    "fail_criteria",
    "falsifiers",
    "current_evidence",
    "limitations",
    "non_claims",
}


def _is_non_empty(value: Any) -> bool:
    if isinstance(value, str):
        return bool(value.strip())
    if isinstance(value, list):
        return len(value) > 0
    if isinstance(value, dict):
        return len(value) > 0
    return value is not None


def load_ledger(path: Path = LEDGER_PATH) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_top_level(ledger: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    for field in REQUIRED_TOP_LEVEL:
        if field not in ledger or not _is_non_empty(ledger[field]):
            errors.append(f"missing top-level field: {field}")
    if not isinstance(ledger.get("records"), list):
        errors.append("top-level records must be a list")
    return errors


def validate_record(record: dict[str, Any], seen_ids: set[str]) -> list[str]:
    errors: list[str] = []
    record_id = str(record.get("id", "<missing id>"))

    for field in REQUIRED_RECORD_FIELDS:
        if field not in record or not _is_non_empty(record[field]):
            errors.append(f"{record_id}: missing or empty field: {field}")

    for field in LIST_FIELDS:
        if field in record and not isinstance(record[field], list):
            errors.append(f"{record_id}: field must be a list: {field}")

    if record.get("prediction_level") not in ALLOWED_LEVELS:
        errors.append(f"{record_id}: invalid prediction_level: {record.get('prediction_level')}")

    if record.get("status") not in ALLOWED_STATUSES:
        errors.append(f"{record_id}: invalid status: {record.get('status')}")

    if record_id in seen_ids:
        errors.append(f"{record_id}: duplicate prediction id")
    seen_ids.add(record_id)

    if len(record.get("pass_criteria", [])) < 2:
        errors.append(f"{record_id}: at least two pass criteria are required")

    if len(record.get("fail_criteria", [])) < 2:
        errors.append(f"{record_id}: at least two fail criteria are required")

    if len(record.get("falsifiers", [])) < 1:
        errors.append(f"{record_id}: at least one falsifier is required")

    if len(record.get("non_claims", [])) < 2:
        errors.append(f"{record_id}: at least two non-claims are required")

    return errors


def validate_ledger(ledger: dict[str, Any]) -> list[str]:
    errors = validate_top_level(ledger)
    seen_ids: set[str] = set()
    for record in ledger.get("records", []):
        if not isinstance(record, dict):
            errors.append("record entry must be an object")
            continue
        errors.extend(validate_record(record, seen_ids))
    return errors


def summarize(ledger: dict[str, Any], validation_errors: list[str]) -> dict[str, Any]:
    records = ledger.get("records", [])
    by_status: dict[str, int] = {}
    by_level: dict[str, int] = {}
    for record in records:
        status = record.get("status", "unknown")
        level = record.get("prediction_level", "unknown")
        by_status[status] = by_status.get(status, 0) + 1
        by_level[level] = by_level.get(level, 0) + 1

    return {
        "schema_version": ledger.get("schema_version"),
        "prediction_records": len(records),
        "validation_passed": not validation_errors,
        "error_count": len(validation_errors),
        "errors": validation_errors,
        "by_status": dict(sorted(by_status.items())),
        "by_level": dict(sorted(by_level.items())),
        "records": [
            {
                "id": record.get("id"),
                "title": record.get("title"),
                "prediction_level": record.get("prediction_level"),
                "status": record.get("status"),
                "claim_type": record.get("claim_type"),
            }
            for record in records
        ],
    }


def write_summary_files(summary: dict[str, Any]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    json_path = OUTPUT_DIR / "prediction_registry_summary.json"
    with json_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    md_path = OUTPUT_DIR / "prediction_registry_summary.md"
    lines = [
        "# Prediction Registry Summary",
        "",
        f"- prediction_records: {summary['prediction_records']}",
        f"- validation_passed: {summary['validation_passed']}",
        f"- error_count: {summary['error_count']}",
        "",
        "## Records",
        "",
        "| id | level | status | claim_type | title |",
        "| --- | --- | --- | --- | --- |",
    ]

    for record in summary["records"]:
        lines.append(
            "| {id} | {prediction_level} | {status} | {claim_type} | {title} |".format(
                **record
            )
        )

    if summary["errors"]:
        lines.extend(["", "## Errors", ""])
        lines.extend(f"- {error}" for error in summary["errors"])

    with md_path.open("w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
        handle.write("\n")


def main() -> int:
    ledger = load_ledger()
    errors = validate_ledger(ledger)
    summary = summarize(ledger, errors)
    write_summary_files(summary)

    print(f"prediction_records: {summary['prediction_records']}")
    print(f"validation_passed: {summary['validation_passed']}")
    print(f"outputs_written: {OUTPUT_DIR.relative_to(ROOT.parent)}/")

    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())

