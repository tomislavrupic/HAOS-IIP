#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import inspect
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from experiments.emergence_ladder.rung3_recovery_trajectory_v2.controls import (  # type: ignore
        ADMISSIBLE_CONTROLS,
        ORACLE_CONTROL,
        run_condition,
    )
    from experiments.emergence_ladder.rung3_recovery_trajectory_v2.fixtures import (  # type: ignore
        PERTURBATION_FAMILIES,
        initial_state,
        make_fixture,
        perturbation,
    )
    from experiments.emergence_ladder.rung3_recovery_trajectory_v2.mechanism import (  # type: ignore
        PROHIBITED_MEMORY_FIELD_FRAGMENTS,
        RelationalMemory,
        SimulationTrace,
        assert_memory_has_no_reference_fields,
        compatibility,
    )
else:
    from .controls import ADMISSIBLE_CONTROLS, ORACLE_CONTROL, run_condition
    from .fixtures import PERTURBATION_FAMILIES, initial_state, make_fixture, perturbation
    from .mechanism import (
        PROHIBITED_MEMORY_FIELD_FRAGMENTS,
        RelationalMemory,
        SimulationTrace,
        assert_memory_has_no_reference_fields,
        compatibility,
    )


ROOT = Path(__file__).resolve().parent
CONTRACT_PATH = ROOT / "precommitment_contract.json"
CALIBRATION_DIR = ROOT / "calibration"
FINAL_DIR = ROOT / "final"
ALL_CONDITIONS = ("target",) + ADMISSIBLE_CONTROLS + (ORACLE_CONTROL,)
EVIDENCE_CONTROLS = tuple(name for name in ADMISSIBLE_CONTROLS if name not in {"parameter_budget_matched"})


def stable_hash(prefix: str, payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(encoded).hexdigest()[:24]}"


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rankdata(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    ranks[order] = np.arange(len(values), dtype=float)
    return ranks


def rank_correlation(left: np.ndarray, right: np.ndarray) -> float:
    if np.std(left) <= 1.0e-12 or np.std(right) <= 1.0e-12:
        return 0.0
    return float(np.corrcoef(rankdata(left), rankdata(right))[0, 1])


def identity_metrics(memory: RelationalMemory, state: np.ndarray) -> tuple[float, float]:
    current = memory.incidence @ state
    valid = memory.valid_edges
    if not np.any(valid):
        return 0.0, 0.0
    signs = np.sign(current[valid])
    agreement = float(np.mean(signs == memory.orientation[valid]))
    oriented = memory.orientation[valid] * current[valid]
    normalized_margin = float(np.mean(np.clip(oriented / 0.08, 0.0, 1.0)))
    return agreement, normalized_margin


def region_error(memory: RelationalMemory, state: np.ndarray, margin: float) -> tuple[float, float]:
    _, gaps = compatibility(memory, state, margin)
    valid = memory.valid_edges
    if not np.any(valid):
        return 1.0, 1.0
    return float(np.sqrt(np.mean(gaps[valid] ** 2))), float(np.mean(gaps[valid] > 0.0))


def corridor_error(reference: np.ndarray, state: np.ndarray) -> float:
    reference_centered = reference - np.mean(reference)
    state_centered = state - np.mean(state)
    return float(np.linalg.norm(state_centered - reference_centered) / max(np.linalg.norm(reference_centered), 1.0e-12))


def functional_error(probes: np.ndarray, reference: np.ndarray, state: np.ndarray) -> float:
    return float(np.linalg.norm(probes @ state - probes @ reference) / max(np.linalg.norm(probes @ reference), 1.0e-12))


def normalized_gain(initial_error: float, final_error: float) -> float:
    if initial_error <= 1.0e-10:
        return 0.0
    return float((initial_error - final_error) / initial_error)


def measure_trace(
    trace: SimulationTrace,
    reference: np.ndarray,
    perturbed: np.ndarray,
    memory: RelationalMemory,
    probes: np.ndarray,
    margin: float,
) -> dict[str, Any]:
    region_series = [region_error(memory, state, margin)[0] for state in trace.states]
    violation_series = [region_error(memory, state, margin)[1] for state in trace.states]
    corridor_series = [corridor_error(reference, state) for state in trace.states]
    identity_series = [identity_metrics(memory, state)[0] for state in trace.states]
    functional_series = [functional_error(probes, reference, state) for state in trace.states]
    post_region = region_series[0]
    final_region = region_series[-1]
    post_corridor = corridor_series[0]
    final_corridor = corridor_series[-1]
    post_functional = functional_series[0]
    final_functional = functional_series[-1]
    variance_ratio = float(np.var(trace.states[-1]) / max(np.var(reference), 1.0e-12))
    return {
        "post_state_error": post_region,
        "final_state_error": final_region,
        "state_recovery_gain": normalized_gain(post_region, final_region),
        "post_trajectory_error": post_corridor,
        "final_trajectory_error": final_corridor,
        "trajectory_recovery_gain": normalized_gain(post_corridor, final_corridor),
        "operator_recovery": float(1.0 - violation_series[-1]),
        "identity_sign_agreement": identity_series[-1],
        "identity_margin_score": identity_metrics(memory, trace.states[-1])[1],
        "node_rank_identity": rank_correlation(reference, trace.states[-1]),
        "functional_recovery": normalized_gain(post_functional, final_functional),
        "final_functional_error": final_functional,
        "overshoot": float(max(0.0, max(region_series) / max(post_region, 1.0e-12) - 1.0)),
        "corrective_cost": float(np.sum(np.square(trace.corrections))),
        "variance_ratio": variance_ratio,
        "mean_shift": float(abs(np.mean(trace.states[-1]) - np.mean(perturbed))),
        "_region_series": region_series,
        "_corridor_series": corridor_series,
        "_identity_series": identity_series,
    }


def execute_run(
    partition: str,
    condition: str,
    seed: int,
    family: str,
    magnitude: float,
    eta: float,
    margin: float,
    contract: dict[str, Any],
) -> dict[str, Any]:
    fixture = make_fixture(contract["system_family"].get("n_side", 8))
    reference = initial_state(fixture, seed)
    shock = perturbation(fixture, family, magnitude, seed)
    perturbed = reference + shock
    memory = RelationalMemory.encode(
        fixture.incidence,
        reference,
        contract["mechanism"]["minimum_memory_edge_magnitude"],
    )
    trace = run_condition(
        condition,
        fixture,
        reference,
        perturbed,
        memory,
        eta,
        margin,
        contract["mechanism"]["steps"],
        contract["mechanism"]["edge_dropout_fraction"],
        seed + int(round(magnitude * 1000)) + PERTURBATION_FAMILIES.index(family) * 10000,
    )
    metrics = measure_trace(trace, reference, perturbed, memory, fixture.probes, margin)
    return {
        "partition": partition,
        "condition": condition,
        "admissible": condition != ORACLE_CONTROL,
        "seed": seed,
        "perturbation_family": family,
        "magnitude": magnitude,
        "eta": eta,
        "inequality_margin": margin,
        **metrics,
    }


def execute_partition(
    partition: str,
    seeds: Iterable[int],
    magnitudes: Iterable[float],
    conditions: Iterable[str],
    eta: float,
    margin: float,
    contract: dict[str, Any],
) -> list[dict[str, Any]]:
    return [
        execute_run(partition, condition, seed, family, magnitude, eta, margin, contract)
        for seed in seeds
        for family in PERTURBATION_FAMILIES
        for magnitude in magnitudes
        for condition in conditions
    ]


def select_parameters(rows_by_parameter: dict[tuple[float, float], list[dict[str, Any]]]) -> tuple[dict[str, float], list[dict[str, Any]]]:
    ranking = []
    for (eta, margin), rows in rows_by_parameter.items():
        organizational = [
            (row["state_recovery_gain"] + row["trajectory_recovery_gain"] + row["identity_sign_agreement"]) / 3.0
            for row in rows
        ]
        score = float(np.median(organizational) - 0.02 * np.median([row["corrective_cost"] for row in rows]))
        ranking.append({"eta": eta, "inequality_margin": margin, "calibration_score": score})
    ranking.sort(key=lambda row: (-row["calibration_score"], row["eta"], row["inequality_margin"]))
    return {"eta": ranking[0]["eta"], "inequality_margin": ranking[0]["inequality_margin"]}, ranking


def derive_thresholds(target_rows: list[dict[str, Any]], control_rows: list[dict[str, Any]]) -> dict[str, float]:
    target_gain = np.array([row["state_recovery_gain"] for row in target_rows], dtype=float)
    control_medians = [
        float(np.median([row["state_recovery_gain"] for row in control_rows if row["condition"] == condition]))
        for condition in ("passive_relaxation", "operator_only_filtering", "signal_blind")
    ]
    margin = float(np.median(target_gain) - max(control_medians))
    return {
        "state_error_max": float(min(0.15, np.quantile([row["final_state_error"] for row in target_rows], 0.9) + 0.02)),
        "trajectory_error_max": float(min(0.20, np.quantile([row["final_trajectory_error"] for row in target_rows], 0.9) + 0.03)),
        "identity_min": float(max(0.85, np.quantile([row["identity_sign_agreement"] for row in target_rows], 0.1) - 0.03)),
        "functional_recovery_min": float(max(0.75, np.quantile([row["functional_recovery"] for row in target_rows], 0.1) - 0.05)),
        "operator_recovery_min": float(max(0.75, np.quantile([row["operator_recovery"] for row in target_rows], 0.1) - 0.05)),
        "recovery_gain_min": float(max(0.10, 0.5 * np.median(target_gain))),
        "control_margin_min": float(max(0.05, 0.5 * margin)),
        "run_recovery_rate_min": 0.75,
        "trivial_variance_ratio_min": 0.5,
    }


def classify_run(row: dict[str, Any], thresholds: dict[str, float], steps: int) -> dict[str, Any]:
    recovered = bool(
        row["final_state_error"] <= thresholds["state_error_max"]
        and row["final_trajectory_error"] <= thresholds["trajectory_error_max"]
        and row["identity_sign_agreement"] >= thresholds["identity_min"]
        and row["functional_recovery"] >= thresholds["functional_recovery_min"]
        and row["operator_recovery"] >= thresholds["operator_recovery_min"]
        and row["state_recovery_gain"] >= thresholds["recovery_gain_min"]
        and row["overshoot"] <= 0.5
    )
    latency = None
    for index, (state_error, identity) in enumerate(zip(row["_region_series"], row["_identity_series"])):
        if state_error <= thresholds["state_error_max"] and identity >= thresholds["identity_min"]:
            latency = index
            break
    if row["overshoot"] > 0.5:
        recovery_type = "unstable correction"
    elif row["variance_ratio"] < thresholds["trivial_variance_ratio_min"] and row["identity_sign_agreement"] < thresholds["identity_min"]:
        recovery_type = "wrong-attractor convergence"
    elif row["condition"] == "passive_relaxation":
        recovery_type = "passive decay"
    elif recovered and row["final_trajectory_error"] <= 0.05:
        recovery_type = "exact-state recovery"
    elif recovered:
        recovery_type = "identity-preserving organizational-regime recovery"
    elif row["identity_sign_agreement"] >= thresholds["identity_min"] and row["state_recovery_gain"] > 0.0:
        recovery_type = "identity-preserving adaptation"
    elif row["functional_recovery"] >= thresholds["functional_recovery_min"]:
        recovery_type = "functional recovery"
    elif row["trajectory_recovery_gain"] > 0.0:
        recovery_type = "trajectory-corridor recovery"
    else:
        recovery_type = "no recovery"
    result = dict(row)
    result["recovered"] = recovered
    result["recovery_latency"] = latency if latency is not None else steps + 1
    result["recovery_type"] = recovery_type
    return result


def bootstrap_ci(values: list[float], seed: int, resamples: int) -> tuple[float, float]:
    array = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    estimates = np.array([np.mean(rng.choice(array, size=len(array), replace=True)) for _ in range(resamples)], dtype=float)
    return float(np.quantile(estimates, 0.025)), float(np.quantile(estimates, 0.975))


def paired_differences(rows: list[dict[str, Any]], left: str, right: str, metric: str) -> list[float]:
    grouped: dict[tuple[int, str, float], dict[str, dict[str, Any]]] = defaultdict(dict)
    for row in rows:
        grouped[(row["seed"], row["perturbation_family"], row["magnitude"])][row["condition"]] = row
    return [conditions[left][metric] - conditions[right][metric] for conditions in grouped.values() if left in conditions and right in conditions]


def control_assertions(rows: list[dict[str, Any]]) -> dict[str, bool]:
    conditions = {row["condition"] for row in rows}
    return {
        "all_controls_present": set(ALL_CONDITIONS) <= conditions,
        "oracle_excluded_from_evidence": all(not row["admissible"] for row in rows if row["condition"] == ORACLE_CONTROL),
        "passive_zero_cost": max(abs(row["corrective_cost"]) for row in rows if row["condition"] == "passive_relaxation") <= 1.0e-12,
        "topology_control_changes_response": any(
            abs(value) > 1.0e-8 for value in paired_differences(rows, "target", "topology_altered", "final_state_error")
        ),
        "identity_scramble_changes_response": any(
            abs(value) > 1.0e-8 for value in paired_differences(rows, "target", "identity_scrambled", "identity_sign_agreement")
        ),
        "signal_blind_consumes_energy": float(np.median([row["corrective_cost"] for row in rows if row["condition"] == "signal_blind"])) > 0.0,
        "redundancy_ablation_is_distinct": any(
            abs(value) > 1.0e-8 for value in paired_differences(rows, "target", "redundancy_ablation", "final_state_error")
        ),
    }


def leakage_guard(contract: dict[str, Any]) -> dict[str, bool]:
    assert_memory_has_no_reference_fields()
    signature = inspect.signature(__import__(RelationalMemory.__module__, fromlist=["simulate_feedback"]).simulate_feedback)
    parameter_names = {name.lower() for name in signature.parameters}
    no_forbidden_parameters = not any(fragment in name for fragment in PROHIBITED_MEMORY_FIELD_FRAGMENTS for name in parameter_names)
    calibration = set(contract["calibration_partition"]["seeds"])
    validation = set(contract["validation_partition"]["seeds"])
    final = set(contract["final_evaluation_partition"]["seeds"])
    return {
        "mechanism_signature_has_no_reference": no_forbidden_parameters,
        "seed_partitions_disjoint": not (calibration & validation or calibration & final or validation & final),
        "calibration_only_selection_declared": bool(contract["seed_policy"]["calibration_only_for_selection"]),
        "memory_stores_only_orientation_bits": set(RelationalMemory.__dataclass_fields__) == {"incidence", "orientation", "valid_edges", "max_degree"},
    }


def aggregate_final(
    rows: list[dict[str, Any]],
    thresholds: dict[str, float],
    contract: dict[str, Any],
    validation_pass: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    target = [row for row in rows if row["condition"] == "target"]
    by_condition = {condition: [row for row in rows if row["condition"] == condition] for condition in ALL_CONDITIONS}
    rates = {condition: float(np.mean([row["recovered"] for row in condition_rows])) for condition, condition_rows in by_condition.items()}
    target_rate = rates["target"]
    passive_diff = paired_differences(rows, "target", "passive_relaxation", "state_recovery_gain")
    ci = bootstrap_ci(passive_diff, contract["seed_policy"]["uncertainty_seed"], contract["seed_policy"]["uncertainty_resamples"])
    basin = {
        str(magnitude): float(np.mean([row["recovered"] for row in target if row["magnitude"] == magnitude]))
        for magnitude in contract["final_evaluation_partition"]["magnitudes"]
    }
    family_rates = {
        family: float(np.mean([row["recovered"] for row in target if row["perturbation_family"] == family]))
        for family in PERTURBATION_FAMILIES
    }
    failure_threshold = next((float(magnitude) for magnitude in contract["final_evaluation_partition"]["magnitudes"] if basin[str(magnitude)] < thresholds["run_recovery_rate_min"]), None)
    assertions = control_assertions(rows)
    leakage = leakage_guard(contract)
    admissible_comparators = ("passive_relaxation", "operator_only_filtering", "signal_blind")
    best_control_rate = max(rates[name] for name in admissible_comparators)
    gates = {
        "validation_pass": validation_pass,
        "positive_gain_vs_passive": float(np.median(passive_diff)) > 0.0,
        "gain_ci_excludes_zero": ci[0] > 0.0,
        "run_recovery_rate": target_rate >= thresholds["run_recovery_rate_min"],
        "identity_preserved": float(np.median([row["identity_sign_agreement"] for row in target])) >= thresholds["identity_min"],
        "trivial_attractor_rejected": rates["trivial_attractor"] < target_rate and float(np.median([row["variance_ratio"] for row in target])) >= thresholds["trivial_variance_ratio_min"],
        "beats_operator_only": target_rate - rates["operator_only_filtering"] >= thresholds["control_margin_min"],
        "beats_signal_blind": target_rate - rates["signal_blind"] >= thresholds["control_margin_min"],
        "memory_ablation_degrades": target_rate - rates["memory_ablation"] >= thresholds["control_margin_min"],
        "redundancy_ablation_degrades": target_rate - rates["redundancy_ablation"] >= thresholds["control_margin_min"],
        "nonzero_basin": sum(rate >= thresholds["run_recovery_rate_min"] for rate in basin.values()) >= 2,
        "multiple_directions": sum(rate >= 0.5 for rate in family_rates.values()) >= 2,
        "no_reference_leakage": all(leakage.values()),
        "controls_valid": all(assertions.values()),
        "evaluation_not_used_for_selection": True,
    }
    unstable_fraction = float(np.mean([row["recovery_type"] == "unstable correction" for row in target]))
    identity_median = float(np.median([row["identity_sign_agreement"] for row in target]))
    if not all(leakage.values()):
        classification = "REFERENCE_LEAKAGE"
    elif not all(assertions.values()):
        classification = "CONTROL_INVALID"
    elif unstable_fraction > 0.1:
        classification = "UNSTABLE_RESTORATION"
    elif target_rate >= 0.25 and identity_median < thresholds["identity_min"]:
        classification = "IDENTITY_FAILURE"
    elif target_rate >= 0.25 and float(np.median([row["variance_ratio"] for row in target])) < thresholds["trivial_variance_ratio_min"]:
        classification = "TRIVIAL_ATTRACTOR"
    elif all(gates.values()):
        classification = "RUNG_3_SUPPORTED"
    elif target_rate >= 0.25 or float(np.median(passive_diff)) > 0.0:
        classification = "PARTIAL_RECOVERY_ONLY"
    else:
        classification = "MECHANISM_NEGATIVE"
    aggregate = {
        "experiment_id": contract["experiment_id"],
        "classification": classification,
        "status": "PASS" if classification == "RUNG_3_SUPPORTED" else "FAIL",
        "selected_parameters": {"eta": target[0]["eta"], "inequality_margin": target[0]["inequality_margin"]},
        "thresholds": thresholds,
        "run_counts": {"final_total": len(rows), "target": len(target), "invalid": 0},
        "dimension_medians": {
            "state_recovery_gain": float(np.median([row["state_recovery_gain"] for row in target])),
            "trajectory_recovery_gain": float(np.median([row["trajectory_recovery_gain"] for row in target])),
            "operator_recovery": float(np.median([row["operator_recovery"] for row in target])),
            "identity_sign_agreement": identity_median,
            "node_rank_identity": float(np.median([row["node_rank_identity"] for row in target])),
            "functional_recovery": float(np.median([row["functional_recovery"] for row in target])),
            "recovery_latency": float(np.median([row["recovery_latency"] for row in target])),
            "overshoot": float(np.median([row["overshoot"] for row in target])),
            "corrective_cost": float(np.median([row["corrective_cost"] for row in target])),
        },
        "recovery_rates": rates,
        "recovery_basin": basin,
        "perturbation_family_rates": family_rates,
        "failure_threshold": failure_threshold,
        "gates": gates,
        "control_assertions": assertions,
        "leakage_guard": leakage,
        "labels": [classification, "IDENTITY_PRESERVING_ORGANIZATIONAL_RECOVERY_TESTED", "NO_HIGHER_RUNG_PROMOTION"],
    }
    uncertainty = {
        "paired_target_minus_passive_gain_ci_95": list(ci),
        "paired_difference_median": float(np.median(passive_diff)),
        "run_unit": "seed x perturbation family x magnitude",
        "resamples": contract["seed_policy"]["uncertainty_resamples"],
        "best_admissible_control_recovery_rate": best_control_rate,
    }
    aggregate["contract_hash"] = stable_hash("el_r3_rt_02_contract", contract)
    aggregate["source_hashes"] = {
        "mechanism": file_hash(ROOT / "mechanism.py"),
        "controls": file_hash(ROOT / "controls.py"),
        "runner": file_hash(ROOT / "run_experiment.py"),
    }
    aggregate["result_hash"] = stable_hash("el_r3_rt_02", aggregate)
    return aggregate, uncertainty


def public_row(row: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in row.items() if not key.startswith("_")}


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    clean = [public_row(row) for row in rows]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(clean[0]))
        writer.writeheader()
        writer.writerows(clean)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def render_classification(aggregate: dict[str, Any], uncertainty: dict[str, Any]) -> str:
    lines = [
        "# EL-R3-RT-02 Final Classification",
        "",
        f"Classification: `{aggregate['classification']}`",
        "",
        "## Recovery Dimensions",
        "",
    ]
    lines.extend(f"- `{name}`: `{value:.6f}`" for name, value in aggregate["dimension_medians"].items())
    lines.extend([
        "",
        f"- Target recovery rate: `{aggregate['recovery_rates']['target']:.6f}`",
        f"- Paired target-minus-passive gain CI: `{uncertainty['paired_target_minus_passive_gain_ci_95']}`",
        f"- Failure threshold: `{aggregate['failure_threshold']}`",
        f"- Result hash: `{aggregate['result_hash']}`",
        "",
        "A positive classification requires every mandatory gate. The oracle is reported as a non-admissible upper bound and never contributes evidence.",
        "",
    ])
    return "\n".join(lines)


def run() -> dict[str, Any]:
    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    if contract.get("status") != "FROZEN":
        raise ValueError("EL-R3-RT-02 contract is not frozen")
    if contract["experiment_id"] != "EL-R3-RT-02":
        raise ValueError("unexpected experiment id")

    grid = contract["mechanism"]["parameter_grid"]
    rows_by_parameter: dict[tuple[float, float], list[dict[str, Any]]] = {}
    for eta in grid["eta"]:
        for margin in grid["inequality_margin"]:
            rows_by_parameter[(eta, margin)] = execute_partition(
                "calibration",
                contract["calibration_partition"]["seeds"],
                contract["calibration_partition"]["magnitudes"],
                ("target",),
                eta,
                margin,
                contract,
            )
    selected, ranking = select_parameters(rows_by_parameter)
    calibration_target = rows_by_parameter[(selected["eta"], selected["inequality_margin"])]
    calibration_controls = execute_partition(
        "calibration",
        contract["calibration_partition"]["seeds"],
        contract["calibration_partition"]["magnitudes"],
        ("passive_relaxation", "operator_only_filtering", "signal_blind"),
        selected["eta"],
        selected["inequality_margin"],
        contract,
    )
    thresholds = derive_thresholds(calibration_target, calibration_controls)
    write_json(CALIBRATION_DIR / "parameter_selection.json", {"selected": selected, "ranking": ranking, "selection_partition": "calibration_only"})
    write_json(CALIBRATION_DIR / "derived_thresholds.json", {"thresholds": thresholds, "source_partition": "calibration_only", "formulas": contract["threshold_derivation"]})
    write_csv(CALIBRATION_DIR / "calibration_runs.csv", calibration_target + calibration_controls)

    validation_raw = execute_partition(
        "validation",
        contract["validation_partition"]["seeds"],
        contract["validation_partition"]["magnitudes"],
        ALL_CONDITIONS,
        selected["eta"],
        selected["inequality_margin"],
        contract,
    )
    validation = [classify_run(row, thresholds, contract["mechanism"]["steps"]) for row in validation_raw]
    write_csv(CALIBRATION_DIR / "validation_runs.csv", validation)
    validation_target_rate = float(np.mean([row["recovered"] for row in validation if row["condition"] == "target"]))
    validation_pass = validation_target_rate >= 0.5 and all(control_assertions(validation).values())

    seed_registry = {
        "partition": "final_untouched_until_execution",
        "seeds": contract["final_evaluation_partition"]["seeds"],
        "magnitudes": contract["final_evaluation_partition"]["magnitudes"],
        "perturbation_families": list(PERTURBATION_FAMILIES),
        "selected_parameters_from": "calibration/parameter_selection.json",
        "thresholds_from": "calibration/derived_thresholds.json",
    }
    write_json(FINAL_DIR / "seed_registry.json", seed_registry)
    final_raw = execute_partition(
        "final",
        contract["final_evaluation_partition"]["seeds"],
        contract["final_evaluation_partition"]["magnitudes"],
        ALL_CONDITIONS,
        selected["eta"],
        selected["inequality_margin"],
        contract,
    )
    final_rows = [classify_run(row, thresholds, contract["mechanism"]["steps"]) for row in final_raw]
    aggregate, uncertainty = aggregate_final(final_rows, thresholds, contract, validation_pass)
    write_csv(FINAL_DIR / "run_level_results.csv", final_rows)
    control_rows = [row for row in final_rows if row["condition"] != "target"]
    write_csv(FINAL_DIR / "control_results.csv", control_rows)
    write_json(FINAL_DIR / "uncertainty.json", uncertainty)
    write_json(FINAL_DIR / "aggregate_result.json", aggregate)
    (FINAL_DIR / "CLASSIFICATION.md").write_text(render_classification(aggregate, uncertainty), encoding="utf-8")
    return aggregate


def main() -> int:
    parser = argparse.ArgumentParser(description="Run frozen EL-R3-RT-02 restorative dynamics experiment.")
    parser.parse_args()
    print(json.dumps(run(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
