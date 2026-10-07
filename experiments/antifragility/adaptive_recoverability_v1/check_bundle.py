#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from experiments.antifragility.adaptive_recoverability_v1.run_experiment import (  # type: ignore
        CONDITIONS,
        aggregate_rows,
        expected_row_count,
        sha256_file,
        stable_hash,
    )
else:
    from .run_experiment import CONDITIONS, aggregate_rows, expected_row_count, sha256_file, stable_hash


ROOT = Path(__file__).resolve().parent


class BundleError(ValueError):
    pass


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def validate_contract(contract: dict[str, Any]) -> None:
    failures: list[str] = []
    if contract.get("status") != "FROZEN":
        failures.append("contract status is not FROZEN")
    if contract.get("candidate_id") != "HAOS-AF-AR-01":
        failures.append("candidate id changed")
    if set(contract["conditions"]["admissible"]) != set(CONDITIONS[:-1]):
        failures.append("admissible control set changed")
    if contract["conditions"]["noncausal_reference_only"] != ["static_hardened_reference"]:
        failures.append("static reference boundary changed")
    if set(contract["final_evaluation"]["seeds"]) & set(contract["smoke_partition"]["seeds"]):
        failures.append("smoke and final seeds overlap")
    if contract["mechanism"]["update_depends_on"] != "exposure severity only":
        failures.append("mechanism gained undeclared information")
    if failures:
        raise BundleError("; ".join(failures))


def row_identity(row: dict[str, Any]) -> tuple[int, int, float, str, float, str]:
    return (
        int(row["node_count"]),
        int(row["seed"]),
        float(row["exposure_severity"]),
        str(row["stress_family"]),
        float(row["test_severity"]),
        str(row["condition"]),
    )


def validate_rows(rows: list[dict[str, Any]], contract: dict[str, Any]) -> None:
    failures: list[str] = []
    identities = [row_identity(row) for row in rows]
    if len(rows) != expected_row_count(contract, "final"):
        failures.append("row count differs from frozen final schedule")
    if len(identities) != len(set(identities)):
        failures.append("duplicate scheduled row")
    expected_conditions = set(CONDITIONS)
    if {str(row["condition"]) for row in rows} != expected_conditions:
        failures.append("condition set differs from frozen schedule")
    if any(row.get("partition") != "final" for row in rows):
        failures.append("non-final partition row in final bundle")
    if any(bool(row.get("location_memory")) for row in rows):
        failures.append("location memory entered an experiment row")
    if any(not bool(row.get("held_out_disjoint")) for row in rows):
        failures.append("exposure and held-out locations overlap")
    if any(bool(row["admissible"]) for row in rows if row["condition"] == "static_hardened_reference"):
        failures.append("static reference became admissible causal evidence")
    if any(not bool(row["admissible"]) for row in rows if row["condition"] != "static_hardened_reference"):
        failures.append("admissible control mislabeled")
    if failures:
        raise BundleError("; ".join(failures))


def validate_claim_language(summary: str) -> None:
    required = ["finite toy experiment", "## Open", "## Not Claimed", "NO_HIGHER_RUNG_PROMOTION"]
    missing = [phrase for phrase in required if phrase not in summary]
    forbidden = [
        "we demonstrate general antifragility",
        "we establish a physical law",
        "we demonstrate biological adaptation",
        "we derive new physics",
    ]
    found = [phrase for phrase in forbidden if phrase in summary.lower()]
    if missing or found:
        raise BundleError(f"claim language invalid; missing={missing}; forbidden={found}")


def canonical_aggregate(result: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in result.items() if key not in {"contract_sha256", "source_hashes", "result_hash"}}


def check_bundle(root: Path = ROOT) -> dict[str, Any]:
    required = [
        "precommitment_contract.json",
        "results/per_run_results.jsonl",
        "results/per_family_summary.csv",
        "final/uncertainty.json",
        "final/aggregate_result.json",
        "final/result_manifest.json",
        "final/summary.md",
    ]
    missing = [relative for relative in required if not (root / relative).is_file()]
    if missing:
        raise BundleError(f"missing required artifacts: {missing}")
    contract = read_json(root / "precommitment_contract.json")
    validate_contract(contract)
    rows = read_jsonl(root / "results/per_run_results.jsonl")
    validate_rows(rows, contract)
    summary = (root / "final" / "summary.md").read_text(encoding="utf-8")
    validate_claim_language(summary)
    stored = read_json(root / "final" / "aggregate_result.json")
    manifest = read_json(root / "final" / "result_manifest.json")
    recomputed, uncertainty, _ = aggregate_rows(rows, contract)
    if canonical_aggregate(stored) != recomputed:
        raise BundleError("stored aggregate differs from deterministic recomputation")
    if read_json(root / "final" / "uncertainty.json") != uncertainty:
        raise BundleError("stored uncertainty differs from deterministic recomputation")
    contract_hash = sha256_file(root / "precommitment_contract.json")
    if stored.get("contract_sha256") != contract_hash or manifest.get("contract_sha256") != contract_hash:
        raise BundleError("contract hash mismatch")
    source_files = {
        "checker": root / "check_bundle.py",
        "model": root / "model.py",
        "runner": root / "run_experiment.py",
    }
    source_hashes = {name: sha256_file(path) for name, path in source_files.items()}
    if stored.get("source_hashes") != source_hashes or manifest.get("source_hashes") != source_hashes:
        raise BundleError("source hash mismatch")
    without_hash = dict(stored)
    stored_hash = without_hash.pop("result_hash", None)
    expected_hash = stable_hash("haos_af_ar_01", without_hash)
    if stored_hash != expected_hash or manifest.get("result_hash") != expected_hash:
        raise BundleError("result hash mismatch")
    if set(manifest["final_seeds"]) & set(manifest["smoke_seeds"]):
        raise BundleError("manifest seed partitions overlap")
    return {
        "candidate_id": stored["candidate_id"],
        "classification": stored["classification"],
        "result_hash": expected_hash,
        "rows_checked": len(rows),
        "status": "PASS",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate HAOS-AF-AR-01 final bundle.")
    parser.add_argument("--root", type=Path, default=ROOT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    print(json.dumps(check_bundle(args.root.resolve()), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
