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
    from experiments.emergence_ladder.rung3_distributed_parity.alphabet import RelationalAlphabet, symbol_agreement
    from experiments.emergence_ladder.rung3_distributed_parity.bridge import correction_vector, simulate_bridge
    from experiments.emergence_ladder.rung3_distributed_parity.controls import (
        ADMISSIBLE_CONDITIONS,
        NON_ADMISSIBLE_CONDITIONS,
        run_condition,
    )
    from experiments.emergence_ladder.rung3_distributed_parity.decoder import decode_local
    from experiments.emergence_ladder.rung3_distributed_parity.parity_architecture import (
        DistributedParityMemory,
        assert_memory_schema,
        information_accounting,
        syndrome,
    )
    from experiments.emergence_ladder.rung3_recovery_trajectory_v2.fixtures import (
        initial_state,
        make_fixture,
        perturbation,
    )
else:
    from .alphabet import RelationalAlphabet, symbol_agreement
    from .bridge import correction_vector, simulate_bridge
    from .controls import ADMISSIBLE_CONDITIONS, NON_ADMISSIBLE_CONDITIONS, run_condition
    from .decoder import decode_local
    from .parity_architecture import DistributedParityMemory, assert_memory_schema, information_accounting, syndrome
    from ..rung3_recovery_trajectory_v2.fixtures import initial_state, make_fixture, perturbation


ROOT = Path(__file__).resolve().parent
CONTRACT_PATH = ROOT / "precommitment_contract.json"
CALIBRATION_DIR = ROOT / "calibration"
VALIDATION_DIR = ROOT / "validation"
FINAL_DIR = ROOT / "final"
PRIMARY_STATE_FAMILIES = (
    "localized_block_shock",
    "sparse_node_corruption",
    "relational_twist",
    "clustered_symbol_corruption",
    "dispersed_equal_weight_corruption",
    "topology_preserving_identity_disruption",
)
MECHANISM_DAMAGE_FAMILIES = ("parity_memory_corruption", "parity_check_deletion")
ALL_FAMILIES = PRIMARY_STATE_FAMILIES + MECHANISM_DAMAGE_FAMILIES
ALL_CONDITIONS = ADMISSIBLE_CONDITIONS + NON_ADMISSIBLE_CONDITIONS
CALIBRATION_CONTROLS = ("target", "rt02_frozen_baseline", "passive_relaxation", "operator_only", "random_parity_checks", "memory_budget_random_bits")


def stable_hash(prefix: str, payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(encoded).hexdigest()[:24]}"


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_contract() -> dict[str, Any]:
    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    if contract.get("status") != "FROZEN_BEFORE_EXECUTION":
        raise ValueError("RP-01 contract is not frozen")
    return contract


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    keys = sorted({key for row in rows for key in row if not key.startswith("_")})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in keys})


def normalized(values: np.ndarray, magnitude: float) -> np.ndarray:
    result = np.asarray(values, dtype=float).copy()
    result -= np.mean(result)
    rms = max(float(np.sqrt(np.mean(result**2))), 1.0e-12)
    return float(magnitude) * result / rms


def scenario_state(reference: np.ndarray, fixture: Any, family: str, magnitude: float, seed: int) -> np.ndarray:
    if family in {"localized_block_shock", "sparse_node_corruption", "relational_twist"}:
        return reference + perturbation(fixture, family, magnitude, seed)
    rng = np.random.default_rng(seed + 33011 + ALL_FAMILIES.index(family) * 101)
    n = fixture.n_side
    shock = np.zeros_like(reference)
    if family == "clustered_symbol_corruption":
        i0 = int(rng.integers(0, n - 2))
        j0 = int(rng.integers(0, n - 2))
        for di in range(3):
            for dj in range(3):
                shock[(i0 + di) * n + j0 + dj] = rng.choice((-1.0, 1.0))
    elif family == "dispersed_equal_weight_corruption":
        chosen = np.asarray([(2 * k * n + (3 * k + 1)) % (n * n) for k in range(9)], dtype=int)
        shock[chosen] = rng.choice((-1.0, 1.0), size=len(chosen))
    elif family == "topology_preserving_identity_disruption":
        shifted = np.roll(reference.reshape((n, n)), shift=1, axis=1).reshape(-1)
        shock = shifted - reference
    elif family in MECHANISM_DAMAGE_FAMILIES:
        return reference + perturbation(fixture, "localized_block_shock", magnitude, seed)
    else:
        raise ValueError(f"unknown perturbation family: {family}")
    return reference + normalized(shock, magnitude)


def scenario_memory(
    memory: DistributedParityMemory,
    family: str,
) -> tuple[DistributedParityMemory, np.ndarray | None]:
    if family == "parity_memory_corruption":
        parity = memory.parity_array().copy()
        parity.reshape(-1)[::8] ^= 1
        altered = DistributedParityMemory(
            tuple(tuple(int(value) for value in row) for row in parity),
            memory.owner_by_block,
        )
        return altered, None
    if family == "parity_check_deletion":
        available = np.ones((len(memory.parity_by_block), 2), dtype=bool)
        available[::4, 0] = False
        return memory, available
    return memory, None


def rank_correlation(left: np.ndarray, right: np.ndarray) -> float:
    left_order = np.argsort(np.argsort(left, kind="mergesort"), kind="mergesort").astype(float)
    right_order = np.argsort(np.argsort(right, kind="mergesort"), kind="mergesort").astype(float)
    if np.std(left_order) <= 1.0e-12 or np.std(right_order) <= 1.0e-12:
        return 0.0
    return float(np.corrcoef(left_order, right_order)[0, 1])


def normalized_gain(initial_error: float, final_error: float) -> float:
    if initial_error <= 1.0e-10:
        return 0.0
    return float((initial_error - final_error) / initial_error)


def region_error(state: np.ndarray, expected_symbols: np.ndarray, alphabet: RelationalAlphabet, margin: float) -> float:
    desired = 2.0 * np.asarray(expected_symbols, dtype=float) - 1.0
    gaps = np.maximum(0.0, float(margin) - desired * (alphabet.incidence @ state))
    return float(np.sqrt(np.mean(gaps**2)))


def corridor_error(reference: np.ndarray, state: np.ndarray) -> float:
    centered_reference = reference - np.mean(reference)
    centered_state = state - np.mean(state)
    return float(np.linalg.norm(centered_state - centered_reference) / max(np.linalg.norm(centered_reference), 1.0e-12))


def functional_error(probes: np.ndarray, reference: np.ndarray, state: np.ndarray) -> float:
    return float(np.linalg.norm(probes @ state - probes @ reference) / max(np.linalg.norm(probes @ reference), 1.0e-12))


def measure_run(
    result: Any,
    reference: np.ndarray,
    perturbed: np.ndarray,
    expected_symbols: np.ndarray,
    alphabet: RelationalAlphabet,
    probes: np.ndarray,
    margin: float,
) -> dict[str, Any]:
    final_state = result.trace.states[-1]
    post_symbols = alphabet.encode(perturbed)
    final_symbols = alphabet.encode(final_state)
    post_region = region_error(perturbed, expected_symbols, alphabet, margin)
    final_region = region_error(final_state, expected_symbols, alphabet, margin)
    post_corridor = corridor_error(reference, perturbed)
    final_corridor = corridor_error(reference, final_state)
    post_function = functional_error(probes, reference, perturbed)
    final_function = functional_error(probes, reference, final_state)
    identity = symbol_agreement(expected_symbols, final_symbols)
    final_syndrome_weight = None
    parity_decode_success = None
    decoded_relation_accuracy = None
    decode_status = "not_applicable"
    wrong_codeword = False
    if result.decode is not None:
        final_syndrome_weight = int(np.sum(result.decode.syndrome_after))
        parity_decode_success = bool(np.array_equal(result.decode.target_symbols, expected_symbols))
        decoded_relation_accuracy = symbol_agreement(expected_symbols, result.decode.target_symbols)
        decode_status = "eliminated" if result.decode.syndrome_eliminated else "unresolved"
        wrong_codeword = bool(result.decode.syndrome_eliminated and not parity_decode_success)
    variance_ratio = float(np.var(final_state) / max(np.var(reference), 1.0e-12))
    return {
        "syndrome_weight_before": int(np.sum(result.decode.syndrome_before)) if result.decode is not None else None,
        "syndrome_weight_after_decode": final_syndrome_weight,
        "syndrome_reduction": (
            int(np.sum(result.decode.syndrome_before)) - int(np.sum(result.decode.syndrome_after))
            if result.decode is not None
            else None
        ),
        "parity_decoding_success": parity_decode_success,
        "decode_status": decode_status,
        "decoded_relation_accuracy": decoded_relation_accuracy,
        "relational_symbol_restoration": identity,
        "post_state_region_error": post_region,
        "final_state_region_error": final_region,
        "state_region_recovery_gain": normalized_gain(post_region, final_region),
        "post_trajectory_corridor_error": post_corridor,
        "final_trajectory_corridor_error": final_corridor,
        "trajectory_corridor_recovery_gain": normalized_gain(post_corridor, final_corridor),
        "operator_restoration": float(max(0.0, 1.0 - final_region / max(float(margin), 1.0e-12))),
        "identity_preservation": identity,
        "node_rank_identity": rank_correlation(reference, final_state),
        "post_functional_error": post_function,
        "final_functional_error": final_function,
        "functional_restoration": normalized_gain(post_function, final_function),
        "corrective_cost": float(np.sum(np.square(result.trace.correction_norms))),
        "correction_peak": float(max(result.trace.correction_norms, default=0.0)),
        "recovery_latency": next((i for i, state in enumerate(result.trace.states) if symbol_agreement(expected_symbols, alphabet.encode(state)) >= 0.90), len(result.trace.states)),
        "overshoot": float(max(0.0, max(corridor_error(reference, state) for state in result.trace.states) / max(post_corridor, 1.0e-12) - 1.0)),
        "variance_ratio": variance_ratio,
        "trivial_attractor": bool(variance_ratio < 0.5),
        "wrong_codeword_convergence": wrong_codeword,
        "post_recovery_persistence": symbol_agreement(final_symbols, alphabet.encode(final_state.copy())),
        "post_symbol_error_weight": int(np.sum(post_symbols ^ expected_symbols)),
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
    fixture = make_fixture(contract["system"]["n_side"])
    alphabet = RelationalAlphabet.build(contract["system"]["n_side"])
    reference = initial_state(fixture, seed)
    expected_symbols = alphabet.encode(reference)
    base_memory = DistributedParityMemory.encode(expected_symbols, alphabet)
    perturbed = scenario_state(reference, fixture, family, magnitude, seed)
    observed_symbols = alphabet.encode(perturbed)
    block_errors = alphabet.block_error_counts(expected_symbols, observed_symbols)
    memory, check_mask = scenario_memory(base_memory, family)
    result = run_condition(
        condition,
        reference,
        perturbed,
        alphabet,
        memory,
        eta,
        margin,
        contract["system"]["correction_steps"],
        seed + int(round(magnitude * 1000)) + ALL_FAMILIES.index(family) * 10000,
        check_mask if condition == "target" else None,
    )
    metrics = measure_run(result, reference, perturbed, expected_symbols, alphabet, fixture.probes, margin)
    max_block_errors = int(np.max(block_errors))
    if max_block_errors <= 1:
        radius_class = "within_radius"
    elif max_block_errors == 2:
        radius_class = "near_above_boundary"
    else:
        radius_class = "above_radius"
    return {
        "partition": partition,
        "condition": condition,
        "admissible": result.admissible,
        "mechanism_note": result.mechanism_note,
        "seed": seed,
        "perturbation_family": family,
        "magnitude": magnitude,
        "eta": eta,
        "bridge_margin": margin,
        "relational_error_weight": int(np.sum(block_errors)),
        "maximum_block_error_weight": max_block_errors,
        "blocks_with_errors": int(np.sum(block_errors > 0)),
        "radius_class": radius_class,
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
    families: Iterable[str] = ALL_FAMILIES,
) -> list[dict[str, Any]]:
    return [
        execute_run(partition, condition, seed, family, magnitude, eta, margin, contract)
        for seed in seeds
        for family in families
        for magnitude in magnitudes
        for condition in conditions
    ]


def select_parameters(rows_by_parameter: dict[tuple[float, float], list[dict[str, Any]]]) -> tuple[dict[str, float], list[dict[str, Any]]]:
    ranking = []
    for (eta, margin), rows in rows_by_parameter.items():
        target = [row for row in rows if row["condition"] == "target"]
        score = float(
            np.median([row["functional_restoration"] for row in target])
            + np.median([row["identity_preservation"] for row in target])
            - 0.01 * np.median([row["corrective_cost"] for row in target])
        )
        ranking.append({"eta": eta, "bridge_margin": margin, "calibration_score": score})
    ranking.sort(key=lambda row: (-row["calibration_score"], row["eta"], row["bridge_margin"]))
    return {"eta": ranking[0]["eta"], "bridge_margin": ranking[0]["bridge_margin"]}, ranking


def derive_thresholds(rows: list[dict[str, Any]], contract: dict[str, Any]) -> dict[str, float]:
    target = [row for row in rows if row["condition"] == "target"]
    function_min = float(max(0.75, np.quantile([row["functional_restoration"] for row in target], 0.1) - 0.05))
    identity_min = float(max(0.90, np.quantile([row["identity_preservation"] for row in target], 0.1) - 0.03))
    provisional = {
        "functional_restoration_min": function_min,
        "identity_min": identity_min,
        "state_region_gain_min": float(contract["threshold_derivation"]["state_region_gain_min"]),
        "trajectory_corridor_gain_min": float(contract["threshold_derivation"]["trajectory_corridor_gain_min"]),
        "trivial_variance_ratio_min": float(contract["threshold_derivation"]["trivial_variance_ratio_min"]),
        "recovery_rate_min": float(contract["threshold_derivation"]["recovery_rate_min"]),
    }
    classified = [classify_row(row, provisional) for row in rows]
    rates = condition_rates(classified)
    best_control = max(rates.get(name, 0.0) for name in ("passive_relaxation", "rt02_frozen_baseline", "random_parity_checks", "memory_budget_random_bits"))
    provisional["superiority_margin_min"] = float(max(0.05, 0.5 * (rates.get("target", 0.0) - best_control)))
    return provisional


def classify_row(row: dict[str, Any], thresholds: dict[str, float]) -> dict[str, Any]:
    recovered = bool(
        row["functional_restoration"] >= thresholds["functional_restoration_min"]
        and row["identity_preservation"] >= thresholds["identity_min"]
        and row["state_region_recovery_gain"] >= thresholds["state_region_gain_min"]
        and row["trajectory_corridor_recovery_gain"] >= thresholds["trajectory_corridor_gain_min"]
        and row["variance_ratio"] >= thresholds["trivial_variance_ratio_min"]
        and not row["wrong_codeword_convergence"]
    )
    result = dict(row)
    result["recovered"] = recovered
    if row["trivial_attractor"]:
        result["recovery_type"] = "trivial-attractor convergence"
    elif row["wrong_codeword_convergence"]:
        result["recovery_type"] = "wrong-codeword convergence"
    elif recovered:
        result["recovery_type"] = "identity-preserving organizational-regime recovery"
    elif row["functional_restoration"] >= thresholds["functional_restoration_min"]:
        result["recovery_type"] = "functional recovery without full regime recovery"
    elif row["identity_preservation"] >= thresholds["identity_min"] and row["state_region_recovery_gain"] > 0:
        result["recovery_type"] = "relational-code recovery only"
    else:
        result["recovery_type"] = "no recovery"
    return result


def condition_rates(rows: list[dict[str, Any]]) -> dict[str, float]:
    grouped: dict[str, list[bool]] = defaultdict(list)
    for row in rows:
        grouped[row["condition"]].append(bool(row.get("recovered", False)))
    return {name: float(np.mean(values)) for name, values in grouped.items()}


def paired_metric(rows: list[dict[str, Any]], left: str, right: str, metric: str) -> list[float]:
    grouped: dict[tuple[int, str, float], dict[str, dict[str, Any]]] = defaultdict(dict)
    for row in rows:
        grouped[(row["seed"], row["perturbation_family"], row["magnitude"])][row["condition"]] = row
    return [pair[left][metric] - pair[right][metric] for pair in grouped.values() if left in pair and right in pair]


def bootstrap_ci(values: list[float], seed: int, resamples: int) -> tuple[float, float]:
    array = np.asarray(values, dtype=float)
    if len(array) == 0:
        return 0.0, 0.0
    rng = np.random.default_rng(seed)
    estimates = np.asarray([np.mean(rng.choice(array, len(array), replace=True)) for _ in range(resamples)])
    return float(np.quantile(estimates, 0.025)), float(np.quantile(estimates, 0.975))


def leakage_guard(contract: dict[str, Any]) -> dict[str, bool]:
    assert_memory_schema()
    final = set(contract["final_evaluation"]["seeds"])
    calibration = set(contract["calibration"]["seeds"])
    validation = set(contract["validation"]["seeds"])
    signatures = [inspect.signature(function).parameters for function in (decode_local, correction_vector, simulate_bridge)]
    forbidden = ("reference", "future", "function", "score", "label", "checkpoint")
    return {
        "seed_partitions_disjoint": not (final & calibration or final & validation or calibration & validation),
        "mechanism_signatures_exclude_prohibited_information": all(
            not any(fragment in name.lower() for fragment in forbidden for name in signature) for signature in signatures
        ),
        "memory_schema_excludes_continuous_state": set(DistributedParityMemory.__dataclass_fields__) == {
            "parity_by_block", "owner_by_block", "architecture_id"
        },
        "function_label_absent_from_contract_mechanism": "function" not in json.dumps(contract["architecture"]).lower(),
        "calibration_only_parameter_selection": "Calibration-only" in contract["parameter_selection"],
    }


def memory_proof(seed: int = 7001) -> dict[str, Any]:
    fixture = make_fixture(8)
    alphabet = RelationalAlphabet.build(8)
    first = initial_state(fixture, seed)
    second = 1.75 * first + 3.25
    z_first = alphabet.encode(first)
    z_second = alphabet.encode(second)
    p_first = DistributedParityMemory.encode(z_first, alphabet).parity_array()
    p_second = DistributedParityMemory.encode(z_second, alphabet).parity_array()
    accounting = information_accounting()
    accounting.update(
        {
            "distinct_continuous_states": bool(not np.allclose(first, second)),
            "same_primary_symbols_for_distinct_states": bool(np.array_equal(z_first, z_second)),
            "same_parity_for_distinct_states": bool(np.array_equal(p_first, p_second)),
            "empirical_non_uniqueness_pass": bool(not np.allclose(first, second) and np.array_equal(p_first, p_second)),
        }
    )
    return accounting


def validation_assessment(rows: list[dict[str, Any]], thresholds: dict[str, float], contract: dict[str, Any]) -> dict[str, Any]:
    primary = [row for row in rows if row["perturbation_family"] in PRIMARY_STATE_FAMILIES]
    rates = condition_rates(primary)
    target = [row for row in primary if row["condition"] == "target"]
    within = [row for row in target if row["radius_class"] == "within_radius"]
    above = [row for row in target if row["radius_class"] != "within_radius"]
    within_rate = float(np.mean([row["recovered"] for row in within])) if within else 0.0
    above_rate = float(np.mean([row["recovered"] for row in above])) if above else 0.0
    target_rate = rates.get("target", 0.0)
    margin = thresholds["superiority_margin_min"]
    deletion = [row for row in rows if row["condition"] == "parity_deletion" and row["perturbation_family"] in PRIMARY_STATE_FAMILIES]
    deletion_rate = float(np.mean([row["recovered"] for row in deletion])) if deletion else 0.0
    memory = memory_proof()
    leakage = leakage_guard(contract)
    ci_passive = bootstrap_ci(
        paired_metric(primary, "target", "passive_relaxation", "functional_restoration"),
        contract["uncertainty"]["seed"],
        contract["uncertainty"]["resamples"],
    )
    ci_rt02 = bootstrap_ci(
        paired_metric(primary, "target", "rt02_frozen_baseline", "functional_restoration"),
        contract["uncertainty"]["seed"] + 1,
        contract["uncertainty"]["resamples"],
    )
    gates = {
        "nonzero_functional_recovery": any(row["functional_restoration"] > 0.0 for row in target),
        "target_exceeds_passive": target_rate - rates.get("passive_relaxation", 0.0) >= margin and ci_passive[0] > 0.0,
        "target_exceeds_frozen_rt02": target_rate - rates.get("rt02_frozen_baseline", 0.0) >= margin and ci_rt02[0] > 0.0,
        "target_exceeds_random_parity_and_random_memory": all(
            target_rate - rates.get(name, 0.0) >= margin for name in ("random_parity_checks", "memory_budget_random_bits")
        ),
        "identity_threshold": float(np.median([row["identity_preservation"] for row in target])) >= thresholds["identity_min"],
        "trivial_attractor_rejected": target_rate > rates.get("trivial_attractor", 0.0) and float(np.median([row["variance_ratio"] for row in target])) >= thresholds["trivial_variance_ratio_min"],
        "no_leakage": all(leakage.values()),
        "parity_deletion_degrades": target_rate - deletion_rate >= margin,
        "beyond_radius_degrades": bool(within and above and within_rate > above_rate),
        "memory_insufficient_for_continuous_reconstruction": bool(memory["empirical_non_uniqueness_pass"] and not memory["continuous_state_reconstruction_possible"]),
    }
    return {
        "validation_pass": bool(all(gates.values())),
        "gates": gates,
        "rates": rates,
        "within_radius_recovery_rate": within_rate,
        "above_radius_recovery_rate": above_rate,
        "within_radius_count": len(within),
        "above_radius_count": len(above),
        "functional_difference_ci_vs_passive": list(ci_passive),
        "functional_difference_ci_vs_rt02": list(ci_rt02),
        "leakage_guard": leakage,
        "information_accounting": memory,
    }


def summarize_dimensions(rows: list[dict[str, Any]]) -> dict[str, float]:
    target = [row for row in rows if row["condition"] == "target" and row["perturbation_family"] in PRIMARY_STATE_FAMILIES]
    fields = (
        "parity_decoding_success", "syndrome_reduction", "relational_symbol_restoration",
        "state_region_recovery_gain", "trajectory_corridor_recovery_gain", "operator_restoration",
        "identity_preservation", "functional_restoration", "recovery_latency", "corrective_cost",
        "wrong_codeword_convergence", "trivial_attractor", "post_recovery_persistence",
    )
    return {
        field: float(np.median([float(row[field]) for row in target if row[field] is not None]))
        for field in fields
    }


def run() -> dict[str, Any]:
    contract = load_contract()
    memory = memory_proof()
    write_json(ROOT / "information_accounting.json", memory)

    rows_by_parameter: dict[tuple[float, float], list[dict[str, Any]]] = {}
    calibration_rows: list[dict[str, Any]] = []
    for eta in contract["calibration"]["bridge_eta"]:
        for margin in contract["calibration"]["bridge_margin"]:
            rows = execute_partition(
                "calibration", contract["calibration"]["seeds"], contract["calibration"]["magnitudes"],
                CALIBRATION_CONTROLS, eta, margin, contract, PRIMARY_STATE_FAMILIES,
            )
            rows_by_parameter[(eta, margin)] = rows
            calibration_rows.extend(rows)
    selected, ranking = select_parameters(rows_by_parameter)
    selected_rows = rows_by_parameter[(selected["eta"], selected["bridge_margin"])]
    thresholds = derive_thresholds(selected_rows, contract)
    write_csv(CALIBRATION_DIR / "calibration_runs.csv", calibration_rows)
    write_json(CALIBRATION_DIR / "parameter_selection.json", {"selected": selected, "ranking": ranking})
    write_json(CALIBRATION_DIR / "derived_thresholds.json", thresholds)

    validation_raw = execute_partition(
        "validation", contract["validation"]["seeds"], contract["validation"]["magnitudes"],
        ALL_CONDITIONS, selected["eta"], selected["bridge_margin"], contract,
    )
    validation_rows = [classify_row(row, thresholds) for row in validation_raw]
    assessment = validation_assessment(validation_rows, thresholds, contract)
    write_csv(VALIDATION_DIR / "validation_runs.csv", validation_rows)
    write_json(VALIDATION_DIR / "uncertainty.json", {
        "method": contract["uncertainty"]["method"],
        "confidence": contract["uncertainty"]["confidence"],
        "functional_difference_ci_vs_passive": assessment["functional_difference_ci_vs_passive"],
        "functional_difference_ci_vs_rt02": assessment["functional_difference_ci_vs_rt02"],
    })

    if not assessment["validation_pass"]:
        classification = "VALIDATION_GATE_FAILED"
        final_accessed = False
        final_rows: list[dict[str, Any]] = []
    else:
        final_accessed = True
        final_raw = execute_partition(
            "final", contract["final_evaluation"]["seeds"], contract["final_evaluation"]["magnitudes"],
            ALL_CONDITIONS, selected["eta"], selected["bridge_margin"], contract,
        )
        final_rows = [classify_row(row, thresholds) for row in final_raw]
        write_csv(FINAL_DIR / "run_level_results.csv", final_rows)
        final_assessment = validation_assessment(final_rows, thresholds, contract)
        classification = "RUNG_3_SUPPORTED" if final_assessment["validation_pass"] else "DISTRIBUTED_FUNCTIONAL_RECOVERY_PARTIAL"
        assessment = final_assessment

    result = {
        "experiment_id": contract["experiment_id"],
        "classification": classification,
        "validation_gate": assessment,
        "final_evaluation_authorized": final_accessed,
        "final_seed_count_consumed": len(contract["final_evaluation"]["seeds"]) if final_accessed else 0,
        "selected_parameters": selected,
        "thresholds": thresholds,
        "run_counts": {
            "calibration": len(calibration_rows),
            "validation": len(validation_rows),
            "final": len(final_rows),
            "invalid": 0,
        },
        "dimension_medians": summarize_dimensions(final_rows if final_rows else validation_rows),
        "radius_results": {
            "within": assessment["within_radius_recovery_rate"],
            "above": assessment["above_radius_recovery_rate"],
        },
        "contract_hash": stable_hash("el_r3_rp_01_contract", contract),
        "source_hashes": {
            name: file_hash(ROOT / filename) for name, filename in {
                "alphabet": "alphabet.py", "parity_architecture": "parity_architecture.py",
                "decoder": "decoder.py", "bridge": "bridge.py", "controls": "controls.py", "runner": "run_experiment.py",
            }.items()
        },
        "historical_rt02_result_hash": "el_r3_rt_02_cade4e057a54d119e6af143c",
        "labels": [classification, "NO_RUNG_4_PROMOTION" if classification != "RUNG_3_SUPPORTED" else "RUNG_4_RETEST_MAY_BE_CONSIDERED"],
    }
    result["result_hash"] = stable_hash("el_r3_rp_01", result)
    write_json(VALIDATION_DIR / "validation_result.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Run EL-R3-RP-01 distributed parity recovery")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    if args.check_only:
        from experiments.emergence_ladder.rung3_distributed_parity.check_bundle import check
        print(json.dumps(check(), indent=2, sort_keys=True))
    else:
        print(json.dumps(run(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

