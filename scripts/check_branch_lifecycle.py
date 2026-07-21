#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = ROOT / "docs/branch_governance/branch_lifecycle_registry.json"
SUMMARY_PATH = ROOT / "docs/branch_governance/branch_lifecycle_summary.md"

ACTIVE_STATUSES = {"ACTIVE_CANDIDATE", "ACTIVE_INSTRUMENT_REPAIR"}
INACTIVE_STATUSES = {
    "FROZEN_REFERENCE",
    "QUARANTINED_INVALID",
    "SPECULATIVE_ONLY",
    "SUPPORTING_CALIBRATION",
    "TERMINAL_NEGATIVE",
}
ALLOWED_STATUSES = ACTIVE_STATUSES | INACTIVE_STATUSES
REQUIRED_CLOSED = {
    "CST-V0.2.2": "TERMINAL_NEGATIVE",
    "BELL-HAOS-B0-B3.2.2": "TERMINAL_NEGATIVE",
    "GEO-INTRINSIC-SPIRAL-V1": "TERMINAL_NEGATIVE",
    "GEO-HIDDEN-01-V1": "TERMINAL_NEGATIVE",
    "HBP-PB01-V1": "QUARANTINED_INVALID",
    "HBP-PB02-V1": "QUARANTINED_INVALID",
    "HBP-PB03-V1": "QUARANTINED_INVALID",
    "HBP-PB04-V1": "QUARANTINED_INVALID",
    "EL-R4-OC-01": "TERMINAL_NEGATIVE",
    "EL-R3-RT-01": "TERMINAL_NEGATIVE",
    "EL-R3-RT-02-RELATIONAL-FEEDBACK": "TERMINAL_NEGATIVE",
    "EL-R3-RP-01-DISTRIBUTED-PARITY": "TERMINAL_NEGATIVE",
}


class LifecycleValidationError(ValueError):
    pass


def load_registry(path: Path = REGISTRY_PATH) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise LifecycleValidationError(message)


def validate_registry(data: dict[str, Any], root: Path = ROOT) -> None:
    _require(data.get("schema_version") == "1.0", "unsupported schema version")
    _require(data.get("authority"), "authority boundary is missing")
    _require(data.get("reopen_policy", {}).get("required_conditions"), "reopen policy is missing")

    definitions = data.get("status_definitions", {})
    _require(set(definitions) == ALLOWED_STATUSES, "status definitions do not match the allowed lifecycle states")

    branches = data.get("branches")
    _require(isinstance(branches, list) and branches, "branches must be a non-empty list")

    ids: set[str] = set()
    active_ids: set[str] = set()
    rows: dict[str, dict[str, Any]] = {}
    for branch in branches:
        branch_id = branch.get("branch_id")
        _require(isinstance(branch_id, str) and branch_id, "branch_id is required")
        _require(branch_id not in ids, f"duplicate branch_id: {branch_id}")
        ids.add(branch_id)
        rows[branch_id] = branch

        status = branch.get("lifecycle_status")
        _require(status in ALLOWED_STATUSES, f"{branch_id}: invalid lifecycle_status {status!r}")
        _require(isinstance(branch.get("expansion_authorized"), bool), f"{branch_id}: expansion_authorized must be boolean")
        _require(isinstance(branch.get("implementation_authorized"), bool), f"{branch_id}: implementation_authorized must be boolean")
        _require(branch.get("result_status"), f"{branch_id}: result_status is required")
        _require(branch.get("reason"), f"{branch_id}: reason is required")
        _require(branch.get("claim_ceiling"), f"{branch_id}: claim_ceiling is required")
        _require(branch.get("terminal_labels"), f"{branch_id}: terminal_labels are required")

        branch_path = root / branch.get("path", "")
        _require(branch_path.exists(), f"{branch_id}: branch path does not exist: {branch_path}")
        artifacts = branch.get("preserved_artifacts")
        _require(isinstance(artifacts, list) and artifacts, f"{branch_id}: preserved_artifacts are required")
        for artifact in artifacts:
            _require((root / artifact).exists(), f"{branch_id}: preserved artifact does not exist: {artifact}")

        if status in INACTIVE_STATUSES:
            _require(not branch["expansion_authorized"], f"{branch_id}: inactive branch cannot authorize expansion")
            _require(not branch["implementation_authorized"], f"{branch_id}: inactive branch cannot authorize implementation")
        else:
            active_ids.add(branch_id)
            _require(branch["expansion_authorized"], f"{branch_id}: active branch must authorize bounded expansion")
            _require(branch.get("evidence_basis"), f"{branch_id}: active branch needs an evidence basis")
            _require(branch.get("stop_conditions"), f"{branch_id}: active branch needs stop conditions")
            _require(branch.get("success_criterion"), f"{branch_id}: active branch needs a success criterion")
            precommitment = branch.get("next_precommitment", {})
            _require(precommitment.get("status") in {"REQUIRED", "FROZEN"}, f"{branch_id}: precommitment status is invalid")
            _require(precommitment.get("path"), f"{branch_id}: precommitment path is required")
            if branch["implementation_authorized"]:
                _require(precommitment["status"] == "FROZEN", f"{branch_id}: implementation requires a frozen precommitment")
                _require((root / precommitment["path"]).exists(), f"{branch_id}: frozen precommitment is missing")

    priority = data.get("active_priority_order")
    _require(isinstance(priority, list), "active_priority_order must be a list")
    _require(len(priority) == len(set(priority)), "active_priority_order contains duplicates")
    _require(set(priority) == active_ids, "active_priority_order must contain every and only active branch")

    for branch_id, expected_status in REQUIRED_CLOSED.items():
        _require(branch_id in rows, f"required closed branch is missing: {branch_id}")
        _require(rows[branch_id]["lifecycle_status"] == expected_status, f"{branch_id}: expected {expected_status}")
        _require(not rows[branch_id]["expansion_authorized"], f"{branch_id}: closed branch was reopened")


def render_summary(data: dict[str, Any]) -> str:
    rows = {row["branch_id"]: row for row in data["branches"]}
    lines = [
        "# HAOS-IIP Branch Lifecycle Summary",
        "",
        "> Generated from `branch_lifecycle_registry.json`. Edit the JSON authority, then regenerate this view.",
        "",
        f"Audit date: `{data['audit_date']}`",
        "",
        "This registry closes exhausted mechanisms without deleting their evidence. A closed branch may only be followed by a new candidate with a new identifier, a new precommitment, independent motivation or evidence, and the prior negative result preserved.",
        "",
        "## Active Queue",
        "",
        "| Priority | Candidate | Lifecycle | Implementation | Next gate |",
        "| ---: | --- | --- | --- | --- |",
    ]
    for index, branch_id in enumerate(data["active_priority_order"], start=1):
        row = rows[branch_id]
        gate = row["next_precommitment"]
        implementation = "authorized" if row["implementation_authorized"] else "blocked"
        lines.append(
            f"| {index} | `{branch_id}` | `{row['lifecycle_status']}` | `{implementation}` | "
            f"`{gate['status']}`: `{gate['path']}` |"
        )

    lines.extend(
        [
            "",
            "The queue is ranked by directness, falsifiability, control quality, holdout readiness, and dependency readiness. It is not ranked by which target is easiest to make pass.",
            "",
            "## Closed And Quarantined",
            "",
            "| Branch | Lifecycle | Frozen reading | Why work stops |",
            "| --- | --- | --- | --- |",
        ]
    )
    closed = [row for row in data["branches"] if row["lifecycle_status"] in {"TERMINAL_NEGATIVE", "QUARANTINED_INVALID"}]
    for row in sorted(closed, key=lambda item: item["branch_id"]):
        lines.append(
            f"| `{row['branch_id']}` | `{row['lifecycle_status']}` | {row['result_status']} | {row['reason']} |"
        )

    lines.extend(
        [
            "",
            "## Retained But Inactive",
            "",
            "| Branch | Lifecycle | Role |",
            "| --- | --- | --- |",
        ]
    )
    retained = [row for row in data["branches"] if row["lifecycle_status"] in {"FROZEN_REFERENCE", "SUPPORTING_CALIBRATION", "SPECULATIVE_ONLY"}]
    for row in sorted(retained, key=lambda item: item["branch_id"]):
        lines.append(f"| `{row['branch_id']}` | `{row['lifecycle_status']}` | {row['reason']} |")

    lines.extend(
        [
            "",
            "## Enforcement",
            "",
            "- Terminal and quarantined branches cannot authorize expansion or implementation.",
            "- Active candidates cannot authorize implementation until their new precommitment exists and is marked `FROZEN`.",
            "- Reference and calibration bundles remain reusable but cannot promote a mechanism claim.",
            "- Reopening requires a new candidate ID, new precommitment, new mechanism or independent evidence, preserved prior results, and no threshold or post-hoc field movement.",
            "- The registry governs development attention only. It does not rewrite any scientific result artifact.",
            "",
            "Validation command:",
            "",
            "```bash",
            "uv run python scripts/check_branch_lifecycle.py",
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate HAOS-IIP branch lifecycle governance.")
    parser.add_argument("--write-summary", action="store_true", help="Regenerate the human-readable summary from JSON.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    data = load_registry()
    validate_registry(data)
    rendered = render_summary(data)
    if args.write_summary:
        SUMMARY_PATH.write_text(rendered, encoding="utf-8")
    else:
        _require(SUMMARY_PATH.exists(), "generated summary is missing")
        _require(SUMMARY_PATH.read_text(encoding="utf-8") == rendered, "generated summary is stale")
    counts: dict[str, int] = {}
    for row in data["branches"]:
        status = row["lifecycle_status"]
        counts[status] = counts.get(status, 0) + 1
    print(json.dumps({"status": "ok", "active_priority_order": data["active_priority_order"], "counts": counts}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
