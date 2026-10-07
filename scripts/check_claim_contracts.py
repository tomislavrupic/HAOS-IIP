#!/usr/bin/env python3
"""Validate machine-readable public claim contracts without third-party packages."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_DIR = ROOT / "docs" / "claim_governance"
SCHEMA_PATH = CONTRACT_DIR / "claim_contract.schema.json"

ENUMS = {
    "claim_class": {"A", "B", "C", "D", "E", "F", "G", "H"},
    "overall_status": {"PASS", "OPEN", "FAIL"},
    "bridge_status": {
        "INTERNAL_DIAGNOSTIC_ONLY", "NUMERICALLY_SUGGESTIVE", "ANALOGY_ONLY",
        "BRIDGE_HYPOTHESIS", "EXTERNALLY_TESTABLE_CANDIDATE",
        "EXTERNALLY_SUPPORTED", "FAILED_UNDER_CURRENT_TESTS",
    },
    "evidence_scope": {
        "INTERNAL_NUMERICAL", "STRUCTURAL", "SCALING", "EXTERNAL_PROXY",
        "PHYSICAL_CORRESPONDENCE",
    },
}

REQUIRED_ROOT = {
    "contract_version", "claim_id", "claim_title", "claim_class",
    "overall_status", "bridge_status", "claim_text", "evidence_scope",
    "selection_history", "perturbation_domain", "telemetry", "controls",
    "reproduction", "claim_boundary",
}

REQUIRED_NESTED = {
    "selection_history": {
        "disclosure_status", "selection_rule", "alternatives_considered",
        "unreported_exploration_known_absent",
    },
    "perturbation_domain": {"included", "omitted", "applicability_statement"},
    "telemetry": {"primary", "alternatives_tested", "sensitivity_status"},
    "controls": {"tested", "omitted_adversarial_controls", "specificity_status"},
    "reproduction": {
        "artifact_status", "execution_status", "independent_replication_status",
        "commands", "artifacts",
    },
    "claim_boundary": {
        "allowed_language", "not_claimed", "downgrade_condition",
        "kill_condition", "physical_correspondence_status",
    },
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _validate_schema(value: Any, schema: dict[str, Any], location: str) -> None:
    """Validate the JSON-Schema subset used by the committed claim schema."""
    if "const" in schema:
        _require(value == schema["const"], f"{location}: expected {schema['const']!r}")
    if "enum" in schema:
        _require(value in schema["enum"], f"{location}: invalid value {value!r}")

    expected_type = schema.get("type")
    type_checks = {
        "object": lambda item: isinstance(item, dict),
        "array": lambda item: isinstance(item, list),
        "string": lambda item: isinstance(item, str),
        "boolean": lambda item: isinstance(item, bool),
    }
    if expected_type in type_checks:
        _require(type_checks[expected_type](value), f"{location}: expected {expected_type}")

    if isinstance(value, str) and "minLength" in schema:
        _require(len(value) >= schema["minLength"], f"{location}: string is too short")
    if isinstance(value, list):
        if "minItems" in schema:
            _require(len(value) >= schema["minItems"], f"{location}: too few items")
        if "items" in schema:
            for index, item in enumerate(value):
                _validate_schema(item, schema["items"], f"{location}[{index}]")
    if isinstance(value, dict):
        required = set(schema.get("required", []))
        missing = required - value.keys()
        _require(not missing, f"{location}: missing fields: {sorted(missing)}")
        properties = schema.get("properties", {})
        if schema.get("additionalProperties") is False:
            unexpected = value.keys() - properties.keys()
            _require(not unexpected, f"{location}: unexpected fields: {sorted(unexpected)}")
        for key, child in value.items():
            if key in properties:
                _validate_schema(child, properties[key], f"{location}.{key}")


def validate_contract(data: dict[str, Any], source: Path) -> None:
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    _validate_schema(data, schema, source.name)

    missing = REQUIRED_ROOT - data.keys()
    _require(not missing, f"{source.name}: missing root fields: {sorted(missing)}")
    _require(data["contract_version"] == "1.0", f"{source.name}: unsupported contract_version")

    for field, allowed in ENUMS.items():
        _require(data[field] in allowed, f"{source.name}: invalid {field}: {data[field]!r}")

    for section, required in REQUIRED_NESTED.items():
        value = data[section]
        _require(isinstance(value, dict), f"{source.name}: {section} must be an object")
        section_missing = required - value.keys()
        _require(not section_missing, f"{source.name}: {section} missing: {sorted(section_missing)}")

    _require(data["selection_history"]["disclosure_status"] in {"COMPLETE", "PARTIAL", "UNKNOWN"}, f"{source.name}: invalid disclosure_status")
    _require(data["telemetry"]["sensitivity_status"] in {"PASS", "OPEN", "FAIL"}, f"{source.name}: invalid telemetry sensitivity status")
    _require(data["controls"]["specificity_status"] in {"PASS", "OPEN", "FAIL"}, f"{source.name}: invalid control specificity status")
    _require(bool(data["perturbation_domain"]["included"]), f"{source.name}: perturbation domain cannot be empty")
    _require(bool(data["telemetry"]["primary"]), f"{source.name}: primary telemetry cannot be empty")
    _require(bool(data["controls"]["tested"]), f"{source.name}: at least one tested control is required")

    for control in data["controls"]["tested"]:
        missing_control = {"name", "matched_properties", "result", "status"} - control.keys()
        _require(not missing_control, f"{source.name}: control missing: {sorted(missing_control)}")
        _require(control["status"] in {"PASS", "OPEN", "FAIL"}, f"{source.name}: invalid control status")

    for artifact in data["reproduction"]["artifacts"]:
        artifact_path = ROOT / artifact
        _require(artifact_path.is_file(), f"{source.name}: missing artifact: {artifact}")

    boundary = data["claim_boundary"]
    replication = data["reproduction"]["independent_replication_status"]
    correspondence = boundary["physical_correspondence_status"]
    if data["bridge_status"] == "EXTERNALLY_SUPPORTED":
        _require(replication == "INDEPENDENTLY_REPLICATED", f"{source.name}: external support requires independent replication")
    if data["evidence_scope"] == "PHYSICAL_CORRESPONDENCE" and data["overall_status"] == "PASS":
        _require(correspondence == "SUPPORTED", f"{source.name}: physical PASS requires supported correspondence")
        _require(replication == "INDEPENDENTLY_REPLICATED", f"{source.name}: physical PASS requires independent replication")
    if data["selection_history"]["disclosure_status"] != "COMPLETE":
        _require(data["overall_status"] != "PASS" or data["evidence_scope"] == "INTERNAL_NUMERICAL", f"{source.name}: incomplete selection history caps non-internal claims at OPEN")


def main() -> int:
    _require(SCHEMA_PATH.is_file(), "claim contract schema is missing")
    contract_paths = sorted(CONTRACT_DIR.glob("*.claim.json"))
    _require(bool(contract_paths), "no claim contracts found")
    for path in contract_paths:
        validate_contract(json.loads(path.read_text(encoding="utf-8")), path)
    print(f"validated {len(contract_paths)} claim contract(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
