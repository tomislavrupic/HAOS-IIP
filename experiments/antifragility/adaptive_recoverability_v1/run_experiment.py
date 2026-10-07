#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Iterable

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from experiments.antifragility.adaptive_recoverability_v1.model import (  # type: ignore
        adaptive_backup_fraction,
        apply_stress,
        disjoint_primary_locations,
        graph_metrics,
        interaction_budget,
        operator_identity,
        primary_capacity_retention,
        random_backup_modifiers,
        recovery_capacity,
        ring_with_backups,
    )
else:
    from .model import (
        adaptive_backup_fraction,
        apply_stress,
        disjoint_primary_locations,
        graph_metrics,
        interaction_budget,
        operator_identity,
        primary_capacity_retention,
        random_backup_modifiers,
        recovery_capacity,
        ring_with_backups,
    )


ROOT = Path(__file__).resolve().parent
CONTRACT_PATH = ROOT / "precommitment_contract.json"
CONDITIONS = (
    "target_adaptive",
    "passive_no_update",
    "matched_random_reallocation",
    "reversed_update",
    "sham_no_stress",
    "static_hardened_reference",
)


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def stable_hash(prefix: str, value: Any) -> str:
    return f"{prefix}_{hashlib.sha256(canonical_json(value).encode('utf-8')).hexdigest()[:24]}"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_contract() -> dict[str, Any]:
    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    if contract.get("status") != "FROZEN":
        raise ValueError("precommitment contract must be FROZEN")
    return contract


def condition_graph(
    condition: str,
    node_count: int,
    seed: int,
    exposure_severity: float,
    contract: dict[str, Any],
) -> tuple[np.ndarray, float, float, bool]:
    baseline = float(contract["graph_family"]["baseline_backup_fraction"])
    gain = float(contract["mechanism"]["adaptation_gain"])
    cap = float(contract["mechanism"]["backup_fraction_cap"])
    target_fraction = adaptive_backup_fraction(baseline, exposure_severity, gain, cap)
    if condition == "target_adaptive":
        return ring_with_backups(node_count, target_fraction), target_fraction, exposure_severity, True
    if condition == "passive_no_update":
        return ring_with_backups(node_count, baseline), baseline, exposure_severity, True
    if condition == "matched_random_reallocation":
        modifiers = random_backup_modifiers(node_count, seed + int(round(1000 * exposure_severity)))
        return ring_with_backups(node_count, target_fraction, modifiers), target_fraction, exposure_severity, True
    if condition == "reversed_update":
        reversed_fraction = max(0.0, baseline - gain * exposure_severity)
        return ring_with_backups(node_count, reversed_fraction), reversed_fraction, exposure_severity, True
    if condition == "sham_no_stress":
        return ring_with_backups(node_count, baseline), baseline, 0.0, True
    if condition == "static_hardened_reference":
        return ring_with_backups(node_count, cap), cap, 0.0, False
    raise ValueError(f"unknown condition: {condition}")


def held_out_start(node_count: int, exposure_start: int, family: str, block_size: int, seed: int) -> int:
    candidates = [int((exposure_start + node_count // 2 + offset) % node_count) for offset in range(-3, 4)]
    rng = np.random.default_rng(seed)
    rng.shuffle(candidates)
    for candidate in candidates:
        if disjoint_primary_locations(node_count, exposure_start, candidate, block_size, family):
            return candidate
    raise ValueError("could not construct a disjoint held-out location")


def generate_rows(contract: dict[str, Any], partition: str) -> list[dict[str, Any]]:
    if partition == "final":
        seeds = [int(value) for value in contract["final_evaluation"]["seeds"]]
        node_counts = [int(value) for value in contract["graph_family"]["node_counts"]]
        exposure_severities = [float(value) for value in contract["exposure"]["severities"]]
        families = list(contract["final_evaluation"]["stress_families"])
        test_severities = [float(value) for value in contract["final_evaluation"]["test_severities"]]
    elif partition == "smoke":
        seeds = [int(value) for value in contract["smoke_partition"]["seeds"]]
        node_counts = [int(contract["graph_family"]["node_counts"][0])]
        exposure_severities = [float(contract["exposure"]["severities"][0])]
        families = list(contract["final_evaluation"]["stress_families"])
        test_severities = [float(contract["final_evaluation"]["test_severities"][0])]
    else:
        raise ValueError(f"unknown partition: {partition}")

    baseline_fraction = float(contract["graph_family"]["baseline_backup_fraction"])
    rows: list[dict[str, Any]] = []
    for node_count in node_counts:
        block_size = int(round(node_count * float(contract["graph_family"]["stress_block_fraction"])))
        baseline_graph = ring_with_backups(node_count, baseline_fraction)
        baseline_metrics = graph_metrics(baseline_graph)
        baseline_budget = interaction_budget(baseline_graph)
        for seed in seeds:
            rng = np.random.default_rng(seed + 17 * node_count)
            exposure_start = int(rng.integers(0, node_count))
            for exposure_severity in exposure_severities:
                for family_index, family in enumerate(families):
                    test_start = held_out_start(
                        node_count,
                        exposure_start,
                        family,
                        block_size,
                        seed + 101 * family_index + int(100 * exposure_severity),
                    )
                    disjoint = disjoint_primary_locations(node_count, exposure_start, test_start, block_size, family)
                    for test_severity in test_severities:
                        for condition in CONDITIONS:
                            adapted_graph, backup_fraction, observed_severity, admissible = condition_graph(
                                condition, node_count, seed, exposure_severity, contract
                            )
                            stressed_graph = apply_stress(
                                adapted_graph, family, test_start, block_size, test_severity
                            )
                            metrics = graph_metrics(stressed_graph)
                            rows.append(
                                {
                                    "admissible": admissible,
                                    "algebraic_connectivity": metrics.algebraic_connectivity,
                                    "backup_fraction": backup_fraction,
                                    "budget_error": abs(interaction_budget(adapted_graph) - baseline_budget),
                                    "condition": condition,
                                    "effective_conductance": metrics.effective_conductance,
                                    "exposure_severity": exposure_severity,
                                    "exposure_start": exposure_start,
                                    "held_out_disjoint": disjoint,
                                    "location_memory": False,
                                    "node_count": node_count,
                                    "observed_exposure_severity": observed_severity,
                                    "operator_identity": operator_identity(baseline_graph, adapted_graph),
                                    "partition": partition,
                                    "primary_capacity_retention": primary_capacity_retention(
                                        backup_fraction, baseline_fraction
                                    ),
                                    "recovery_capacity": recovery_capacity(metrics, baseline_metrics),
                                    "seed": seed,
                                    "stress_family": family,
                                    "test_severity": test_severity,
                                    "test_start": test_start,
                                    "weight_budget": interaction_budget(adapted_graph),
                                }
                            )
    return rows


def row_key(row: dict[str, Any]) -> tuple[int, int, float, str, float]:
    return (
        int(row["node_count"]),
        int(row["seed"]),
        float(row["exposure_severity"]),
        str(row["stress_family"]),
        float(row["test_severity"]),
    )


def paired_differences(rows: list[dict[str, Any]], left: str, right: str) -> list[tuple[int, float]]:
    indexed = {(row_key(row), str(row["condition"])): row for row in rows}
    differences: list[tuple[int, float]] = []
    keys = sorted({row_key(row) for row in rows})
    for key in keys:
        left_row = indexed[(key, left)]
        right_row = indexed[(key, right)]
        differences.append((int(key[1]), float(left_row["recovery_capacity"]) - float(right_row["recovery_capacity"])))
    return differences


def seed_block_interval(
    paired: list[tuple[int, float]], confidence: float, resamples: int, seed: int
) -> dict[str, float]:
    by_seed: dict[int, list[float]] = {}
    for block, value in paired:
        by_seed.setdefault(block, []).append(value)
    block_values = np.asarray([float(np.mean(by_seed[key])) for key in sorted(by_seed)], dtype=float)
    rng = np.random.default_rng(seed)
    boot = np.empty(resamples, dtype=float)
    for index in range(resamples):
        boot[index] = float(np.mean(rng.choice(block_values, size=len(block_values), replace=True)))
    tail = (1.0 - confidence) / 2.0
    return {
        "estimate": float(np.mean(block_values)),
        "lower": float(np.quantile(boot, tail)),
        "upper": float(np.quantile(boot, 1.0 - tail)),
    }


def collapse_index(rows: Iterable[dict[str, Any]], threshold: float) -> float | None:
    ordered = sorted(rows, key=lambda row: float(row["test_severity"]))
    for row in ordered:
        if float(row["recovery_capacity"]) < threshold:
            return float(row["test_severity"])
    return None


def aggregate_rows(rows: list[dict[str, Any]], contract: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    confidence = float(contract["uncertainty"]["confidence"])
    resamples = int(contract["uncertainty"]["resamples"])
    uncertainty_seed = int(contract["uncertainty"]["seed"])
    target_passive = paired_differences(rows, "target_adaptive", "passive_no_update")
    target_random = paired_differences(rows, "target_adaptive", "matched_random_reallocation")
    target_reversed = paired_differences(rows, "target_adaptive", "reversed_update")
    sham_passive = paired_differences(rows, "sham_no_stress", "passive_no_update")
    intervals = {
        "target_minus_passive": seed_block_interval(target_passive, confidence, resamples, uncertainty_seed),
        "target_minus_random": seed_block_interval(target_random, confidence, resamples, uncertainty_seed + 1),
        "target_minus_reversed": seed_block_interval(target_reversed, confidence, resamples, uncertainty_seed + 2),
        "sham_minus_passive": seed_block_interval(sham_passive, confidence, resamples, uncertainty_seed + 3),
    }

    target_index = {row_key(row): row for row in rows if row["condition"] == "target_adaptive"}
    passive_index = {row_key(row): row for row in rows if row["condition"] == "passive_no_update"}
    target_gains = np.asarray(
        [float(target_index[key]["recovery_capacity"]) - float(passive_index[key]["recovery_capacity"]) for key in sorted(target_index)],
        dtype=float,
    )
    family_summary: list[dict[str, Any]] = []
    for family in contract["final_evaluation"]["stress_families"]:
        family_values = [
            float(target_index[key]["recovery_capacity"]) - float(passive_index[key]["recovery_capacity"])
            for key in sorted(target_index)
            if key[3] == family
        ]
        family_summary.append(
            {
                "stress_family": family,
                "median_gain": float(np.median(family_values)),
                "minimum_gain": float(np.min(family_values)),
                "rows": len(family_values),
            }
        )
    dose_summary: list[dict[str, Any]] = []
    observed_exposure_severities = sorted({float(key[2]) for key in target_index})
    for severity in observed_exposure_severities:
        values = [
            float(target_index[key]["recovery_capacity"]) - float(passive_index[key]["recovery_capacity"])
            for key in sorted(target_index)
            if key[2] == float(severity)
        ]
        dose_summary.append({"exposure_severity": float(severity), "median_gain": float(np.median(values))})

    threshold = float(contract["decision_rules"]["collapse_capacity_threshold"])
    collapse_groups: list[dict[str, Any]] = []
    group_keys = sorted({(key[0], key[1], key[2], key[3]) for key in target_index})
    for group_key in group_keys:
        target_group = [
            row for key, row in target_index.items() if (key[0], key[1], key[2], key[3]) == group_key
        ]
        passive_group = [
            row for key, row in passive_index.items() if (key[0], key[1], key[2], key[3]) == group_key
        ]
        target_k = collapse_index(target_group, threshold)
        passive_k = collapse_index(passive_group, threshold)
        sentinel = 1.0 + max(float(value) for value in contract["final_evaluation"]["test_severities"])
        target_order = sentinel if target_k is None else target_k
        passive_order = sentinel if passive_k is None else passive_k
        collapse_groups.append(
            {
                "node_count": group_key[0],
                "seed": group_key[1],
                "exposure_severity": group_key[2],
                "stress_family": group_key[3],
                "target_k_star": target_k,
                "passive_k_star": passive_k,
                "not_earlier": target_order >= passive_order,
                "shifted_later": target_order > passive_order,
            }
        )

    rules = contract["decision_rules"]
    target_rows = list(target_index.values())
    all_admissible_rows = [row for row in rows if bool(row["admissible"])]
    primary_groups = [row for row in collapse_groups if row["stress_family"] != "backup_channel_damage"]
    gates = {
        "budget_invariant": max(float(row["budget_error"]) for row in all_admissible_rows) <= float(rules["budget_error_max"]),
        "collapse_never_earlier": all(bool(row["not_earlier"]) for row in collapse_groups),
        "collapse_shift_present": any(bool(row["shifted_later"]) for row in primary_groups),
        "dose_response_nondecreasing": all(
            dose_summary[index]["median_gain"] <= dose_summary[index + 1]["median_gain"] + 1.0e-12
            for index in range(len(dose_summary) - 1)
        ),
        "family_generalization": all(
            float(row["median_gain"]) >= float(rules["minimum_family_median_gain"]) for row in family_summary
        ),
        "held_out_locations_disjoint": all(bool(row["held_out_disjoint"]) for row in target_rows),
        "identity_preserved": min(float(row["operator_identity"]) for row in target_rows) >= float(rules["identity_cosine_min"]),
        "no_target_regression": float(np.min(target_gains)) >= -1.0e-12,
        "primary_capacity_preserved": min(float(row["primary_capacity_retention"]) for row in target_rows)
        >= float(rules["minimum_primary_capacity_retention"]),
        "sham_inert": abs(float(intervals["sham_minus_passive"]["estimate"])) <= float(rules["sham_gain_abs_max"]),
        "tail_gain": float(np.quantile(target_gains, 0.10)) >= float(rules["minimum_tail_gain_p10"]),
        "target_beats_passive": float(intervals["target_minus_passive"]["lower"]) >= float(rules["minimum_target_gain_ci_lower"]),
        "target_beats_random": float(intervals["target_minus_random"]["lower"]) > float(rules["minimum_target_minus_random_ci_lower"]),
        "target_beats_reversed": float(intervals["target_minus_reversed"]["lower"]) > float(rules["minimum_target_minus_reversed_ci_lower"]),
    }
    if all(gates.values()):
        classification = "BOUNDED_TOY_ANTIFRAGILITY_SUPPORTED"
    elif not gates["identity_preserved"] or not gates["budget_invariant"] or not gates["primary_capacity_preserved"]:
        classification = "IDENTITY_OR_BUDGET_GATE_FAILED"
    elif not gates["family_generalization"] or not gates["tail_gain"] or not gates["held_out_locations_disjoint"]:
        classification = "HELD_OUT_GENERALIZATION_FAILED"
    elif not gates["target_beats_passive"] or not gates["target_beats_random"] or not gates["target_beats_reversed"]:
        classification = "ADAPTIVE_GAIN_NOT_CONTROL_DISTINCT"
    else:
        classification = "ANTIFRAGILITY_GATE_FAILED"

    static_rows = [row for row in rows if row["condition"] == "static_hardened_reference"]
    aggregate = {
        "candidate_id": contract["candidate_id"],
        "claim_ceiling": contract["claim_ceiling"],
        "classification": classification,
        "collapse_groups": collapse_groups,
        "dose_summary": dose_summary,
        "family_summary": family_summary,
        "gates": gates,
        "labels": [classification, "FINITE_TOY_GRAPH_ONLY", "NO_HIGHER_RUNG_PROMOTION"],
        "run_counts": {
            "admissible": len(all_admissible_rows),
            "static_reference": len(static_rows),
            "target": len(target_rows),
            "total": len(rows),
        },
        "status": "PASS" if all(gates.values()) else "FAIL",
        "static_reference_disclosure": "A statically hardened graph can reach the same endpoint without exposure; it is a noncausal ceiling reference, not evidence for the adaptive mechanism.",
        "target_gain": {
            "median": float(np.median(target_gains)),
            "minimum": float(np.min(target_gains)),
            "p10": float(np.quantile(target_gains, 0.10)),
            "maximum": float(np.max(target_gains)),
        },
        "target_minima": {
            "operator_identity": min(float(row["operator_identity"]) for row in target_rows),
            "primary_capacity_retention": min(float(row["primary_capacity_retention"]) for row in target_rows),
        },
    }
    return aggregate, intervals, family_summary


def expected_row_count(contract: dict[str, Any], partition: str) -> int:
    if partition == "final":
        return (
            len(contract["graph_family"]["node_counts"])
            * len(contract["final_evaluation"]["seeds"])
            * len(contract["exposure"]["severities"])
            * len(contract["final_evaluation"]["stress_families"])
            * len(contract["final_evaluation"]["test_severities"])
            * len(CONDITIONS)
        )
    return (
        len(contract["smoke_partition"]["seeds"])
        * len(contract["final_evaluation"]["stress_families"])
        * len(CONDITIONS)
    )


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(json.dumps(row, sort_keys=True) + "\n" for row in rows), encoding="utf-8")


def write_family_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["stress_family", "median_gain", "minimum_gain", "rows"])
        writer.writeheader()
        writer.writerows(rows)


def summary_markdown(aggregate: dict[str, Any], uncertainty: dict[str, Any]) -> str:
    demonstrated = [name for name, passed in aggregate["gates"].items() if passed]
    failed = [name for name, passed in aggregate["gates"].items() if not passed]
    lines = [
        "# HAOS-AF-AR-01 Final Summary",
        "",
        f"Classification: `{aggregate['classification']}`",
        "",
        "This finite toy experiment tests an operational antifragility criterion. It does not derive general antifragility from HAOS-IIP.",
        "",
        "## Demonstrated",
        "",
    ]
    lines.extend(f"- `{name}`" for name in demonstrated)
    lines.extend(["", "## Failed", ""])
    lines.extend([f"- `{name}`" for name in failed] or ["- None under the frozen finite schedule."])
    lines.extend(
        [
            "",
            "## Core Result",
            "",
            f"- median held-out recovery-capacity gain: `{aggregate['target_gain']['median']:.12g}`",
            f"- tenth-percentile gain: `{aggregate['target_gain']['p10']:.12g}`",
            f"- paired target-minus-passive CI: `[{uncertainty['target_minus_passive']['lower']:.12g}, {uncertainty['target_minus_passive']['upper']:.12g}]`",
            f"- paired target-minus-random CI: `[{uncertainty['target_minus_random']['lower']:.12g}, {uncertainty['target_minus_random']['upper']:.12g}]`",
            f"- minimum operator-identity cosine: `{aggregate['target_minima']['operator_identity']:.12g}`",
            f"- minimum primary-capacity retention: `{aggregate['target_minima']['primary_capacity_retention']:.12g}`",
            "",
            "## Static Reference Boundary",
            "",
            aggregate["static_reference_disclosure"],
            "The experiment therefore attributes a bounded stress-triggered transition, not a unique or optimal design principle.",
            "",
            "## Open",
            "",
            "- Other graph families, adaptive laws, resource budgets, functions, and stress distributions.",
            "- Durability beyond the declared single-update held-out test schedule.",
            "- Independent implementation and replication.",
            "- Whether the criterion remains discriminative in empirical systems.",
            "",
            "## Not Claimed",
            "",
            "- No general antifragility, physical law, biological adaptation, economic advantage, agency, universality, continuum limit, or new physics.",
            "- `NO_HIGHER_RUNG_PROMOTION`: this sidecar does not alter the emergence ladder or any frozen parent result.",
            "",
        ]
    )
    return "\n".join(lines)


def run(partition: str, output_root: Path) -> dict[str, Any]:
    contract = load_contract()
    final_result_path = output_root / "final" / "aggregate_result.json"
    if partition == "final" and final_result_path.exists():
        raise FileExistsError(f"final partition is single-use at this output root: {final_result_path}")
    rows = generate_rows(contract, partition)
    expected = expected_row_count(contract, partition)
    if len(rows) != expected:
        raise ValueError(f"row count {len(rows)} differs from frozen schedule {expected}")
    aggregate, uncertainty, family_summary = aggregate_rows(rows, contract)
    aggregate["contract_sha256"] = sha256_file(CONTRACT_PATH)
    source_files = {
        "checker": ROOT / "check_bundle.py",
        "model": ROOT / "model.py",
        "runner": ROOT / "run_experiment.py",
    }
    aggregate["source_hashes"] = {name: sha256_file(path) for name, path in source_files.items()}
    without_hash = dict(aggregate)
    aggregate["result_hash"] = stable_hash("haos_af_ar_01", without_hash)

    write_rows(output_root / "results" / "per_run_results.jsonl", rows)
    write_family_csv(output_root / "results" / "per_family_summary.csv", family_summary)
    write_json(output_root / "final" / "uncertainty.json", uncertainty)
    write_json(output_root / "final" / "aggregate_result.json", aggregate)
    summary = summary_markdown(aggregate, uncertainty)
    (output_root / "final").mkdir(parents=True, exist_ok=True)
    (output_root / "final" / "summary.md").write_text(summary, encoding="utf-8")
    manifest = {
        "candidate_id": contract["candidate_id"],
        "contract_sha256": aggregate["contract_sha256"],
        "final_seeds": contract["final_evaluation"]["seeds"],
        "partition": partition,
        "result_hash": aggregate["result_hash"],
        "smoke_seeds": contract["smoke_partition"]["seeds"],
        "source_hashes": aggregate["source_hashes"],
    }
    write_json(output_root / "final" / "result_manifest.json", manifest)
    return aggregate


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run HAOS-AF-AR-01 bounded adaptive-recoverability experiment.")
    parser.add_argument("--partition", choices=("smoke", "final"), required=True)
    parser.add_argument("--output-root", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output_root is None:
        if args.partition == "final":
            output_root = ROOT
        else:
            raise ValueError("smoke runs require --output-root outside the frozen experiment folder")
    else:
        output_root = args.output_root.resolve()
    result = run(args.partition, output_root)
    print(json.dumps({key: result[key] for key in ("candidate_id", "classification", "result_hash", "status")}, indent=2))


if __name__ == "__main__":
    main()
