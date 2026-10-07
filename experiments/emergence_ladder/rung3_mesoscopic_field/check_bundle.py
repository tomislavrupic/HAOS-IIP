#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from experiments.emergence_ladder.rung3_mesoscopic_field.controls import (  # type: ignore
        ADMISSIBLE_CONDITIONS,
        ALL_CONDITIONS,
        ORACLE_CONDITIONS,
    )
    from experiments.emergence_ladder.rung3_mesoscopic_field.fixtures import PERTURBATION_FAMILIES  # type: ignore
    from experiments.emergence_ladder.rung3_mesoscopic_field.run_experiment import (  # type: ignore
        aggregate_rows,
        stable_hash,
    )
else:
    from .controls import ADMISSIBLE_CONDITIONS, ALL_CONDITIONS, ORACLE_CONDITIONS
    from .fixtures import PERTURBATION_FAMILIES
    from .run_experiment import aggregate_rows, stable_hash


ROOT = Path(__file__).resolve().parent
EXPECTED_PARENT_HASHES = {
    "experiments/emergence_ladder/rung3_recovery_trajectory_v2/final/aggregate_result.json": "7f88faba8018f75e9700de240d2b142f5a05dedec7a28e8e8874ad36bc196962",
    "experiments/emergence_ladder/rung3_distributed_parity/validation/validation_result.json": "0cbf2dd8ada3baa67f8dc498f5929d9b8a4a3f7700f234b58967d7d3f6d0fbe6",
    "telemetry/frozen_metrics.py": "daa7a759ad8f0ed67249f194a7e7e0753c78a4eb4c0a8b44d47e2e7a0d121f87",
}


class BundleError(ValueError):
    pass


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def expected_keys(contract: dict[str, Any], partition: str) -> set[tuple[int, int, str, float, str]]:
    if partition == "final":
        levels = [int(item["n_side"]) for item in contract["scale_levels"]]
        seeds = [int(value) for value in contract["final_evaluation"]["seeds"]]
        magnitudes = [float(value) for value in contract["final_evaluation"]["magnitudes"]]
    elif partition == "smoke":
        levels = [8]
        seeds = [int(value) for value in contract["smoke_partition"]["seeds"]]
        magnitudes = [float(contract["final_evaluation"]["magnitudes"][0])]
    else:
        raise BundleError(f"unknown partition {partition}")
    return {
        (level, seed, family, magnitude, condition)
        for level in levels
        for seed in seeds
        for family in PERTURBATION_FAMILIES
        for magnitude in magnitudes
        for condition in ALL_CONDITIONS
    }


def validate_contract(contract: dict[str, Any]) -> None:
    failures: list[str] = []
    if contract.get("status") != "FROZEN":
        failures.append("contract status is not FROZEN")
    if contract.get("candidate_id") != "EL-R3-MESOSCOPIC-FIELD-01":
        failures.append("candidate id changed")
    quantizer = contract.get("quantizer", {})
    frozen_quantizer = {
        "adaptive_range": False,
        "bits_per_value": 8,
        "code_count": 256,
        "fixed_max": 3.0,
        "fixed_min": -3.0,
        "non_finite_input": "abort run",
        "overflow": "saturate to fixed endpoints and report saturation count",
        "per_sample_normalization": False,
        "rounding": "round half to even",
        "signed": True,
    }
    for key, expected in frozen_quantizer.items():
        if quantizer.get(key) != expected:
            failures.append(f"quantizer {key} changed")
    accounting = contract.get("information_accounting", {})
    if accounting.get("payload_values") != 16 or accounting.get("payload_bits") != 128:
        failures.append("target information budget is not exactly 16 values / 128 bits")
    if accounting.get("addressing_bits_per_sample") != 0:
        failures.append("per-sample addressing bits are nonzero")
    if list(contract.get("conditions", {}).get("admissible", [])) != list(ADMISSIBLE_CONDITIONS):
        failures.append("admissible condition list changed")
    if list(contract.get("conditions", {}).get("non_admissible_oracles", [])) != list(ORACLE_CONDITIONS):
        failures.append("oracle condition list changed")
    if "direct_four_functional_phasors" in contract.get("conditions", {}).get("admissible", []):
        failures.append("direct functional phasors became admissible")
    if set(contract.get("final_evaluation", {}).get("seeds", [])) & set(contract.get("smoke_partition", {}).get("seeds", [])):
        failures.append("smoke and final seeds overlap")
    if failures:
        raise BundleError("; ".join(failures))


def validate_rows_semantics(rows: list[dict[str, Any]], contract: dict[str, Any], partition: str) -> None:
    expected = expected_keys(contract, partition)
    observed = [
        (int(row["n_side"]), int(row["seed"]), str(row["perturbation_family"]), float(row["magnitude"]), str(row["condition"]))
        for row in rows
    ]
    failures: list[str] = []
    if len(observed) != len(set(observed)):
        failures.append("duplicate scheduled rows")
    if set(observed) != expected:
        failures.append("row schedule differs from frozen Cartesian product")
    final_seeds = set(int(value) for value in contract["final_evaluation"]["seeds"])
    smoke_seeds = set(int(value) for value in contract["smoke_partition"]["seeds"])
    allowed = final_seeds if partition == "final" else smoke_seeds
    if any(int(row["seed"]) not in allowed for row in rows):
        failures.append("row uses a seed outside its partition")
    for row in rows:
        condition = str(row["condition"])
        n_side = int(row["n_side"])
        if str(row.get("partition")) != partition:
            failures.append("row partition label mismatch")
            break
        if int(row.get("coarse_rank", -1)) != 16 or int(row.get("microscopic_nullity", -1)) != n_side * n_side - 16:
            failures.append("rank or nullity mismatch")
            break
        if int(row.get("target_payload_bits", -1)) != 128:
            failures.append("target budget annotation mismatch")
            break
        expected_bits = {
            "target_quantized_coarse_field": 128,
            "block_permutation": 128,
            "phase_scrambling": 128,
            "amplitude_only": 128,
            "phase_only": 128,
            "equal_budget_random_projection": 128,
            "compression_12_target_blind": 96,
            "compression_8_target_blind": 64,
        }.get(condition)
        if expected_bits is not None and int(row.get("payload_bits", -1)) != expected_bits:
            failures.append(f"payload mismatch for {condition}")
            break
        if condition in ADMISSIBLE_CONDITIONS and not bool(row.get("admissible")):
            failures.append(f"admissible condition mislabeled: {condition}")
            break
        if condition in ORACLE_CONDITIONS and bool(row.get("admissible")):
            failures.append(f"oracle condition mislabeled admissible: {condition}")
            break
        if condition == "target_quantized_coarse_field":
            tolerance = float(contract["decision_rules"]["field_exact_tolerance"])
            if float(row.get("stored_quantized_field_error", float("inf"))) > tolerance:
                failures.append("target did not restore stored quantized field exactly")
                break
            if not bool(row.get("local_representation")):
                failures.append("target representation lost locality label")
                break
    if failures:
        raise BundleError("; ".join(failures))


def validate_claim_language(summary: str, contract: dict[str, Any]) -> None:
    required = ["## Open", "## Not claimed", "NO_HIGHER_RUNG_PROMOTION"]
    missing = [phrase for phrase in required if phrase not in summary]
    affirmative_forbidden = (
        "we demonstrate a continuum limit",
        "we establish universality",
        "we derive spacetime",
        "we demonstrate general emergence",
        "we establish a new physical law",
    )
    lowered = summary.lower()
    found = [phrase for phrase in affirmative_forbidden if phrase in lowered]
    if missing or found:
        raise BundleError(f"claim language invalid; missing={missing}; forbidden={found}")
    if not contract["claim_ceiling"].lower().startswith("bounded functional recovery"):
        raise BundleError("claim ceiling is not bounded")


def _canonical_aggregate(result: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in result.items() if key not in {"contract_sha256", "source_hashes", "result_hash"}}


def check_bundle(root: Path = ROOT) -> dict[str, Any]:
    required = [
        "precommitment_contract.json", "results/per_run_results.jsonl", "results/per_level_summary.csv",
        "results/information_budget.csv", "results/control_validity.json", "final/uncertainty.json",
        "final/aggregate_result.json", "final/result_manifest.json", "final/summary.md",
        "audit/hostile_audit.md", "plots/functional_recovery_vs_refinement.svg", "plots/fixed_budget_vs_nullity.svg",
    ]
    missing = [name for name in required if not (root / name).is_file()]
    if missing:
        raise BundleError(f"missing required artifacts: {missing}")
    contract = read_json(root / "precommitment_contract.json")
    validate_contract(contract)
    rows = read_jsonl(root / "results/per_run_results.jsonl")
    validate_rows_semantics(rows, contract, "final")
    aggregate = read_json(root / "final/aggregate_result.json")
    manifest = read_json(root / "final/result_manifest.json")
    validate_claim_language((root / "final/summary.md").read_text(encoding="utf-8"), contract)
    if aggregate.get("classification") == "FIXED_BUDGET_NOT_SCALE_STABLE" and aggregate.get("gates", {}).get("fixed_budget_scale_stable") is True:
        raise BundleError("classification contradicts the frozen scale-stability gate")

    recomputed, uncertainty, level_rows = aggregate_rows(rows, contract)
    if _canonical_aggregate(aggregate) != recomputed:
        raise BundleError("stored aggregate differs from deterministic row recomputation")
    if read_json(root / "final/uncertainty.json") != uncertainty:
        raise BundleError("stored uncertainty differs from deterministic recomputation")
    if read_json(root / "results/control_validity.json") != aggregate["control_assertions"]:
        raise BundleError("control validity artifact differs from aggregate")
    with (root / "results/per_level_summary.csv").open(encoding="utf-8", newline="") as handle:
        csv_rows = list(csv.DictReader(handle))
    if len(csv_rows) != len(level_rows):
        raise BundleError("per-level CSV row count mismatch")
    for csv_row, level_row in zip(csv_rows, level_rows, strict=True):
        if int(csv_row["n_side"]) != int(level_row["n_side"]):
            raise BundleError("per-level CSV ordering mismatch")
        for metric in ("functional_recovery_median", "functional_recovery_rate", "applicable_fraction"):
            if abs(float(csv_row[metric]) - float(level_row[metric])) > 1.0e-12:
                raise BundleError(f"per-level CSV metric mismatch: {metric}")

    contract_hash = sha256_file(root / "precommitment_contract.json")
    if aggregate.get("contract_sha256") != contract_hash or manifest.get("contract_sha256") != contract_hash:
        raise BundleError("contract hash mismatch")
    source_map = {"representation": "representation.py", "fixtures": "fixtures.py", "controls": "controls.py", "runner": "run_experiment.py", "checker": "check_bundle.py"}
    source_hashes = {name: sha256_file(root / filename) for name, filename in source_map.items()}
    if aggregate.get("source_hashes") != source_hashes or manifest.get("source_hashes") != source_hashes:
        raise BundleError("source hash mismatch")
    if manifest.get("parent_artifact_hashes") != EXPECTED_PARENT_HASHES:
        raise BundleError("frozen parent artifact hash mismatch")
    if set(manifest.get("final_seeds", [])) & set(manifest.get("smoke_seeds", [])):
        raise BundleError("manifest smoke/final seed overlap")
    for relative, expected_hash in manifest.get("artifact_hashes", {}).items():
        path = root / relative
        if not path.is_file() or sha256_file(path) != expected_hash:
            raise BundleError(f"artifact hash mismatch: {relative}")
    result_without_hash = dict(aggregate)
    stored_result_hash = result_without_hash.pop("result_hash", None)
    expected_result_hash = stable_hash("el_r3_mesoscopic_field_01", result_without_hash)
    if stored_result_hash != expected_result_hash or manifest.get("result_hash") != expected_result_hash:
        raise BundleError("result hash mismatch")
    return {"candidate_id": aggregate["candidate_id"], "classification": aggregate["classification"], "result_hash": expected_result_hash, "rows_checked": len(rows), "status": "PASS"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Semantically validate EL-R3-MESOSCOPIC-FIELD-01.")
    parser.add_argument("--root", type=Path, default=ROOT)
    return parser.parse_args()


def main() -> int:
    try:
        result = check_bundle(parse_args().root.resolve())
    except (BundleError, OSError, KeyError, TypeError, json.JSONDecodeError) as error:
        print(json.dumps({"status": "FAIL", "error": str(error)}, indent=2, sort_keys=True))
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
