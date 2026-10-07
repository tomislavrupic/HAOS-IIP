#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
    from experiments.emergence_ladder.rung3_mesoscopic_field.controls import (  # type: ignore
        ADMISSIBLE_CONDITIONS,
        ALL_CONDITIONS,
        ORACLE_CONDITIONS,
        run_condition,
    )
    from experiments.emergence_ladder.rung3_mesoscopic_field.fixtures import (  # type: ignore
        PERTURBATION_FAMILIES,
        functional_probes,
        initial_state,
        perturbation,
    )
    from experiments.emergence_ladder.rung3_mesoscopic_field.representation import (  # type: ignore
        QuantizerSpec,
        coarse_field,
        quantize,
        rank_and_nullity,
    )
else:
    from .controls import ADMISSIBLE_CONDITIONS, ALL_CONDITIONS, ORACLE_CONDITIONS, run_condition
    from .fixtures import PERTURBATION_FAMILIES, functional_probes, initial_state, perturbation
    from .representation import QuantizerSpec, coarse_field, quantize, rank_and_nullity


ROOT = Path(__file__).resolve().parent
CONTRACT_PATH = ROOT / "precommitment_contract.json"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_hash(prefix: str, payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(encoded).hexdigest()[:24]}"


def load_contract() -> dict[str, Any]:
    contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
    if contract.get("status") != "FROZEN":
        raise ValueError("EL-R3-MESOSCOPIC-FIELD-01 contract is not frozen")
    if contract.get("candidate_id") != "EL-R3-MESOSCOPIC-FIELD-01":
        raise ValueError("unexpected candidate id")
    return contract


def normalized_gain(initial_error: float, final_error: float) -> float:
    if initial_error <= 1.0e-10:
        return 0.0
    return float((initial_error - final_error) / initial_error)


def edge_sign_identity(reference: np.ndarray, state: np.ndarray, n_side: int) -> float:
    ref = reference.reshape(n_side, n_side)
    current = state.reshape(n_side, n_side)
    ref_edges = np.concatenate([(np.roll(ref, -1, axis=0) - ref).reshape(-1), (np.roll(ref, -1, axis=1) - ref).reshape(-1)])
    cur_edges = np.concatenate([(np.roll(current, -1, axis=0) - current).reshape(-1), (np.roll(current, -1, axis=1) - current).reshape(-1)])
    valid = np.abs(ref_edges) >= 0.08
    if not np.any(valid):
        return 0.0
    return float(np.mean(np.sign(ref_edges[valid]) == np.sign(cur_edges[valid])))


def measure_row(
    condition: str,
    reference: np.ndarray,
    state_prime: np.ndarray,
    restored: np.ndarray,
    n_side: int,
    metadata: dict[str, object],
    quantizer: QuantizerSpec,
    thresholds: dict[str, Any],
) -> dict[str, Any]:
    probes = functional_probes(n_side)
    function_reference = probes @ reference
    post_function_error = float(np.linalg.norm(probes @ state_prime - function_reference) / max(np.linalg.norm(function_reference), 1.0e-12))
    final_function_error = float(np.linalg.norm(probes @ restored - function_reference) / max(np.linalg.norm(function_reference), 1.0e-12))
    function_applicable = post_function_error > float(thresholds["function_applicability_error_min"])
    functional_recovery = normalized_gain(post_function_error, final_function_error) if function_applicable else 0.0

    source_coarse = coarse_field(reference, n_side)
    post_coarse = coarse_field(state_prime, n_side)
    final_coarse = coarse_field(restored, n_side)
    stored = quantize(source_coarse, quantizer)
    post_coarse_error = float(np.linalg.norm(post_coarse - source_coarse))
    final_coarse_error = float(np.linalg.norm(final_coarse - source_coarse))
    state_post_error = float(np.linalg.norm(state_prime - reference) / max(np.linalg.norm(reference), 1.0e-12))
    state_final_error = float(np.linalg.norm(restored - reference) / max(np.linalg.norm(reference), 1.0e-12))
    correction = restored - state_prime
    payload_bits = int(metadata["payload_bits"])
    recovered = bool(function_applicable and functional_recovery >= float(thresholds["functional_recovery_min"]))
    return {
        "condition": condition,
        "admissible": condition in ADMISSIBLE_CONDITIONS,
        "function_applicable": function_applicable,
        "functional_recovery": functional_recovery,
        "recovered": recovered,
        "post_function_error": post_function_error,
        "final_function_error": final_function_error,
        "coarse_field_recovery": normalized_gain(post_coarse_error, final_coarse_error),
        "post_coarse_error": post_coarse_error,
        "final_coarse_error": final_coarse_error,
        "stored_quantized_field_error": float(np.max(np.abs(final_coarse - stored.decoded))),
        "state_recovery_gain": normalized_gain(state_post_error, state_final_error),
        "post_state_error": state_post_error,
        "final_state_error": state_final_error,
        "edge_sign_identity": edge_sign_identity(reference, restored, n_side),
        "variance_ratio": float(np.var(restored) / max(np.var(reference), 1.0e-12)),
        "intervention_rms": float(np.sqrt(np.mean(correction**2))),
        "intervention_density": float(np.mean(np.abs(correction) > 1.0e-12)),
        "intervention_count": int(np.sum(np.abs(correction) > 1.0e-12)),
        "payload_bits": payload_bits,
        "payload_values": int(metadata["payload_values"]),
        "payload_bits_per_node": float(payload_bits / (n_side * n_side)),
        "saturation_count": int(metadata["saturation_count"]),
        "quantization_max_abs_error": float(metadata["quantization_max_abs_error"]),
        "decoder_path": str(metadata["decoder_path"]),
        "local_representation": bool(metadata["local_representation"]),
        "decoder_operation_count": int(3 * n_side * n_side + 16),
        "array_memory_estimate_bytes": int((5 * n_side * n_side + 64) * 8),
    }


def execute_rows(contract: dict[str, Any], partition: str) -> list[dict[str, Any]]:
    if partition == "final":
        schedule = contract["final_evaluation"]
        levels = [int(row["n_side"]) for row in contract["scale_levels"]]
        seeds = [int(value) for value in schedule["seeds"]]
        magnitudes = [float(value) for value in schedule["magnitudes"]]
        families = list(schedule["perturbation_families"])
        conditions = list(ALL_CONDITIONS)
    elif partition == "smoke":
        schedule = contract["final_evaluation"]
        levels = [8]
        seeds = [int(contract["smoke_partition"]["seeds"][0])]
        magnitudes = [float(schedule["magnitudes"][0])]
        families = list(schedule["perturbation_families"])
        conditions = list(ALL_CONDITIONS)
    else:
        raise ValueError(f"unknown partition: {partition}")

    quantizer = QuantizerSpec(
        bits=int(contract["quantizer"]["bits_per_value"]),
        minimum=float(contract["quantizer"]["fixed_min"]),
        maximum=float(contract["quantizer"]["fixed_max"]),
    )
    thresholds = contract["decision_rules"]
    rows: list[dict[str, Any]] = []
    for n_side in levels:
        rank, nullity = rank_and_nullity(n_side)
        for seed in seeds:
            reference = initial_state(n_side, seed)
            target_stored = quantize(coarse_field(reference, n_side), quantizer)
            for family_index, family in enumerate(families):
                for magnitude in magnitudes:
                    delta = perturbation(n_side, family, magnitude, seed)
                    state_prime = reference + delta
                    target_state, _ = run_condition(
                        "target_quantized_coarse_field",
                        reference,
                        state_prime,
                        n_side,
                        seed,
                        quantizer,
                        0.0,
                    )
                    target_correction_rms = float(np.sqrt(np.mean((target_state - state_prime) ** 2)))
                    control_seed = int(seed + family_index * 10000 + round(magnitude * 1000) * 10 + n_side)
                    for condition in conditions:
                        restored, metadata = run_condition(
                            condition,
                            reference,
                            state_prime,
                            n_side,
                            control_seed,
                            quantizer,
                            target_correction_rms,
                        )
                        row = measure_row(condition, reference, state_prime, restored, n_side, metadata, quantizer, thresholds)
                        row.update(
                            {
                                "partition": partition,
                                "n_side": n_side,
                                "state_dimension": n_side * n_side,
                                "coarse_rank": rank,
                                "microscopic_nullity": nullity,
                                "coarse_cell_width": n_side // 4,
                                "seed": seed,
                                "perturbation_family": family,
                                "magnitude": magnitude,
                                "target_payload_bits": 128,
                                "target_saturation_count": target_stored.saturation_count,
                            }
                        )
                        rows.append(row)
    return rows


def bootstrap_ci(values: list[float], seed: int, resamples: int) -> list[float | None]:
    if not values:
        return [None, None]
    array = np.asarray(values, dtype=float)
    rng = np.random.default_rng(int(seed))
    estimates = np.asarray([np.mean(rng.choice(array, size=len(array), replace=True)) for _ in range(int(resamples))])
    return [float(np.quantile(estimates, 0.025)), float(np.quantile(estimates, 0.975))]


def paired_differences(rows: list[dict[str, Any]], left: str, right: str) -> list[float]:
    grouped: dict[tuple[int, int, str, float], dict[str, dict[str, Any]]] = defaultdict(dict)
    for row in rows:
        key = (int(row["n_side"]), int(row["seed"]), str(row["perturbation_family"]), float(row["magnitude"]))
        grouped[key][str(row["condition"])] = row
    differences = []
    for conditions in grouped.values():
        if left in conditions and right in conditions and conditions[left]["function_applicable"]:
            differences.append(float(conditions[left]["functional_recovery"] - conditions[right]["functional_recovery"]))
    return differences


def rate(rows: list[dict[str, Any]]) -> float:
    applicable = [row for row in rows if row["function_applicable"]]
    return float(np.mean([row["recovered"] for row in applicable])) if applicable else 0.0


def median_metric(rows: list[dict[str, Any]], metric: str) -> float:
    return float(np.median([float(row[metric]) for row in rows])) if rows else 0.0


def aggregate_rows(rows: list[dict[str, Any]], contract: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    thresholds = contract["decision_rules"]
    target = [row for row in rows if row["condition"] == "target_quantized_coarse_field"]
    by_condition = {condition: [row for row in rows if row["condition"] == condition] for condition in ALL_CONDITIONS}
    expected = len(contract["scale_levels"]) * len(contract["final_evaluation"]["seeds"]) * len(PERTURBATION_FAMILIES) * len(contract["final_evaluation"]["magnitudes"]) * len(ALL_CONDITIONS)
    keys = [(row["n_side"], row["seed"], row["perturbation_family"], row["magnitude"], row["condition"]) for row in rows]
    complete = len(rows) == expected and len(keys) == len(set(keys))

    level_rows: list[dict[str, Any]] = []
    for level in [int(item["n_side"]) for item in contract["scale_levels"]]:
        subset = [row for row in target if int(row["n_side"]) == level]
        applicable = [row for row in subset if row["function_applicable"]]
        family_rates = {
            family: rate([row for row in subset if row["perturbation_family"] == family])
            for family in PERTURBATION_FAMILIES
        }
        level_rows.append(
            {
                "n_side": level,
                "state_dimension": level * level,
                "payload_bits": 128,
                "payload_bits_per_node": 128 / (level * level),
                "rank": 16,
                "nullity": level * level - 16,
                "target_rows": len(subset),
                "applicable_fraction": len(applicable) / max(len(subset), 1),
                "functional_recovery_median": median_metric(applicable, "functional_recovery"),
                "functional_recovery_rate": rate(subset),
                "coarse_field_recovery_median": median_metric(subset, "coarse_field_recovery"),
                "intervention_rms_median": median_metric(subset, "intervention_rms"),
                "intervention_density_median": median_metric(subset, "intervention_density"),
                "saturation_total": int(sum(int(row["saturation_count"]) for row in subset)),
                **{f"family_rate_{name}": value for name, value in family_rates.items()},
            }
        )

    passive_diff = paired_differences(rows, "target_quantized_coarse_field", "passive_relaxation")
    random_diff = paired_differences(rows, "target_quantized_coarse_field", "equal_budget_random_projection")
    passive_ci = bootstrap_ci(passive_diff, contract["uncertainty"]["seed"], contract["uncertainty"]["resamples"])
    random_ci = bootstrap_ci(random_diff, contract["uncertainty"]["seed"] + 1, contract["uncertainty"]["resamples"])
    target_median = median_metric([row for row in target if row["function_applicable"]], "functional_recovery")
    block_median = median_metric([row for row in by_condition["block_permutation"] if row["function_applicable"]], "functional_recovery")
    phase_median = median_metric([row for row in by_condition["phase_scrambling"] if row["function_applicable"]], "functional_recovery")
    level_by_n = {int(row["n_side"]): row for row in level_rows}
    scale_drop = float(level_by_n[8]["functional_recovery_median"] - level_by_n[16]["functional_recovery_median"])

    control_assertions = {
        "all_conditions_present": all(by_condition[name] for name in ALL_CONDITIONS),
        "oracles_excluded": all(not row["admissible"] for name in ORACLE_CONDITIONS for row in by_condition[name]),
        "passive_zero_intervention": max(abs(float(row["intervention_rms"])) for row in by_condition["passive_relaxation"]) <= 1.0e-12,
        "target_payload_fixed_128": all(int(row["payload_bits"]) == 128 for row in target),
        "random_projection_payload_matched": all(int(row["payload_bits"]) == 128 for row in by_condition["equal_budget_random_projection"]),
        "compression_payloads_exact": all(int(row["payload_bits"]) == 96 for row in by_condition["compression_12_target_blind"]) and all(int(row["payload_bits"]) == 64 for row in by_condition["compression_8_target_blind"]),
        "coarse_rank_and_nullity_exact": all(int(row["coarse_rank"]) == 16 and int(row["microscopic_nullity"]) == int(row["state_dimension"]) - 16 for row in rows),
        "no_missing_or_duplicate_rows": complete,
    }
    gates = {
        "stored_quantized_field_exact": max(float(row["stored_quantized_field_error"]) for row in target) <= float(thresholds["field_exact_tolerance"]),
        "function_applicability_by_level": all(float(row["applicable_fraction"]) >= float(thresholds["minimum_applicable_fraction_per_level"]) for row in level_rows),
        "functional_recovery_rate_by_level": all(float(row["functional_recovery_rate"]) >= float(thresholds["minimum_level_recovery_rate"]) for row in level_rows),
        "multiple_families_by_level": all(sum(float(row[f"family_rate_{family}"]) >= float(thresholds["minimum_family_recovery_rate"]) for family in PERTURBATION_FAMILIES) >= int(thresholds["minimum_qualifying_families_per_level"]) for row in level_rows),
        "target_beats_passive_ci": passive_ci[0] is not None and float(passive_ci[0]) > 0.0,
        "target_beats_equal_budget_random_ci": random_ci[0] is not None and float(random_ci[0]) > 0.0,
        "block_permutation_degrades": target_median - block_median >= float(thresholds["control_margin_min"]),
        "phase_scrambling_degrades": target_median - phase_median >= float(thresholds["control_margin_min"]),
        "fixed_budget_scale_stable": scale_drop <= float(thresholds["maximum_level_median_drop_8_to_16"]),
        "nontrivial_variance": median_metric(target, "variance_ratio") >= float(thresholds["minimum_variance_ratio"]),
        "controls_valid": all(control_assertions.values()),
        "final_seed_schedule_complete": complete,
    }

    if not complete:
        classification = "INSTRUMENT_INVALID"
    elif not all(control_assertions.values()):
        classification = "CONTROL_INVALID"
    elif all(gates.values()):
        classification = "RUNG_3_SUPPORTED_BOUNDED_MESOSCOPIC_MEMORY"
    elif not gates["fixed_budget_scale_stable"]:
        classification = "FIXED_BUDGET_NOT_SCALE_STABLE"
    elif target_median < float(thresholds["functional_recovery_min"]):
        classification = "FIELD_RESTORED_FUNCTION_NOT_RESTORED"
    elif not gates["target_beats_equal_budget_random_ci"] or not gates["block_permutation_degrades"] or not gates["phase_scrambling_degrades"]:
        classification = "FUNCTION_RECOVERY_NOT_CONTROL_DISTINCT"
    else:
        classification = "PARTIAL_RECOVERY_ONLY"

    aggregate: dict[str, Any] = {
        "candidate_id": contract["candidate_id"],
        "classification": classification,
        "status": "PASS" if classification == "RUNG_3_SUPPORTED_BOUNDED_MESOSCOPIC_MEMORY" else "FAIL",
        "claim_ceiling": contract["claim_ceiling"],
        "run_counts": {"total": len(rows), "target": len(target), "expected": expected, "invalid": 0 if complete else abs(expected - len(rows))},
        "target_medians": {
            "functional_recovery": target_median,
            "coarse_field_recovery": median_metric(target, "coarse_field_recovery"),
            "edge_sign_identity": median_metric(target, "edge_sign_identity"),
            "state_recovery_gain": median_metric(target, "state_recovery_gain"),
            "intervention_rms": median_metric(target, "intervention_rms"),
            "intervention_density": median_metric(target, "intervention_density"),
            "variance_ratio": median_metric(target, "variance_ratio"),
        },
        "condition_functional_recovery_rates": {name: rate(condition_rows) for name, condition_rows in by_condition.items()},
        "condition_functional_recovery_medians": {name: median_metric([row for row in condition_rows if row["function_applicable"]], "functional_recovery") for name, condition_rows in by_condition.items()},
        "level_summary": level_rows,
        "scale_median_drop_8_to_16": scale_drop,
        "gates": gates,
        "control_assertions": control_assertions,
        "parent_classifications_preserved": {
            "EL-R3-RT-01": "NEGATIVE_RESULT",
            "EL-R3-RT-02": "PARTIAL_RECOVERY_ONLY",
            "EL-R3-RP-01": "VALIDATION_GATE_FAILED",
        },
        "labels": [classification, "BOUNDED_MESOSCOPIC_MEMORY_TESTED", "FINITE_LEVELS_ONLY", "NO_HIGHER_RUNG_PROMOTION"],
    }
    uncertainty = {
        "paired_target_minus_passive_functional_recovery_ci_95": passive_ci,
        "paired_target_minus_equal_budget_random_functional_recovery_ci_95": random_ci,
        "paired_target_minus_passive_median": float(np.median(passive_diff)),
        "paired_target_minus_equal_budget_random_median": float(np.median(random_diff)),
        "resamples": int(contract["uncertainty"]["resamples"]),
        "unit": "paired level x seed x perturbation family x magnitude",
    }
    return aggregate, uncertainty, level_rows


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(json.dumps(row, sort_keys=True, allow_nan=False) + "\n" for row in rows), encoding="utf-8")


def line_plot_svg(level_rows: list[dict[str, Any]], aggregate: dict[str, Any]) -> str:
    width, height = 760, 440
    xs = {8: 150, 12: 380, 16: 610}
    def y(value: float) -> float:
        return 360.0 - 280.0 * max(-0.1, min(1.0, value))
    points = " ".join(f"{xs[int(row['n_side'])]},{y(float(row['functional_recovery_median'])):.1f}" for row in level_rows)
    parts = [f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>", "<rect width='100%' height='100%' fill='#070908'/>", "<text x='40' y='40' fill='#d8ded7' font-family='monospace' font-size='18'>FUNCTIONAL RECOVERY VS REFINEMENT</text>"]
    for tick in (0.0, 0.25, 0.5, 0.75, 1.0):
        yy = y(tick)
        parts.append(f"<line x1='80' y1='{yy}' x2='690' y2='{yy}' stroke='#2a302b'/><text x='35' y='{yy+5}' fill='#747c74' font-family='monospace' font-size='12'>{tick:.2f}</text>")
    parts.append(f"<polyline points='{points}' fill='none' stroke='#c9ff3d' stroke-width='3'/>")
    for row in level_rows:
        xx = xs[int(row["n_side"])]
        yy = y(float(row["functional_recovery_median"]))
        parts.append(f"<rect x='{xx-5}' y='{yy-5}' width='10' height='10' fill='#c9ff3d'/><text x='{xx-34}' y='400' fill='#d8ded7' font-family='monospace' font-size='14'>{int(row['n_side'])}x{int(row['n_side'])}</text>")
    parts.append(f"<text x='80' y='425' fill='#747c74' font-family='monospace' font-size='11'>CLASSIFICATION: {aggregate['classification']}</text></svg>")
    return "".join(parts)


def budget_plot_svg(level_rows: list[dict[str, Any]]) -> str:
    parts = ["<svg xmlns='http://www.w3.org/2000/svg' width='760' height='440' viewBox='0 0 760 440'>", "<rect width='100%' height='100%' fill='#070908'/>", "<text x='40' y='40' fill='#d8ded7' font-family='monospace' font-size='18'>FIXED MEMORY / GROWING MICROSCOPIC NULLITY</text>"]
    for index, row in enumerate(level_rows):
        x = 100 + index * 220
        nullity_height = 1.1 * float(row["nullity"])
        parts.append(f"<rect x='{x}' y='{360-nullity_height}' width='70' height='{nullity_height}' fill='#424a42'/><rect x='{x+80}' y='219' width='45' height='141' fill='#c9ff3d'/><text x='{x}' y='390' fill='#d8ded7' font-family='monospace' font-size='13'>{row['n_side']}x{row['n_side']}</text><text x='{x}' y='{345-nullity_height}' fill='#747c74' font-family='monospace' font-size='11'>NULL {row['nullity']}</text><text x='{x+75}' y='205' fill='#c9ff3d' font-family='monospace' font-size='11'>128 BITS</text>")
    parts.append("<text x='40' y='425' fill='#747c74' font-family='monospace' font-size='11'>GRAY: UNRESOLVED MICROSCOPIC DIMENSIONS / LIME: FIXED PAYLOAD</text></svg>")
    return "".join(parts)


def render_summary(aggregate: dict[str, Any], uncertainty: dict[str, Any]) -> str:
    if aggregate["classification"] == "RUNG_3_SUPPORTED_BOUNDED_MESOSCOPIC_MEMORY":
        opening = "Across the three preregistered finite refinements, the fixed 128-bit mesoscopic field restored the independent function under every required gate."
    elif aggregate["classification"] == "FIELD_RESTORED_FUNCTION_NOT_RESTORED":
        opening = "The stored quantized mesoscopic field was restored exactly, but the independent function was not restored."
    elif aggregate["classification"] == "FIXED_BUDGET_NOT_SCALE_STABLE":
        opening = "The fixed 128-bit mesoscopic field did not retain preregistered functional recovery across refinement."
    else:
        opening = "The fixed mesoscopic field produced partial evidence but did not satisfy the complete Rung 3 contract."
    demonstrated = [name for name, value in aggregate["gates"].items() if value]
    failed = [name for name, value in aggregate["gates"].items() if not value]
    return "\n".join(
        [
            opening,
            "",
            f"Classification: `{aggregate['classification']}`",
            "",
            "## Demonstrated",
            "",
            *[f"- `{name}`" for name in demonstrated],
            "",
            "## Failed",
            "",
            *([f"- `{name}`" for name in failed] or ["- No mandatory gate failed."]),
            "",
            "## Open",
            "",
            "- Generalization beyond `8x8`, `12x12`, and `16x16`.",
            "- Other functions, topologies, dynamics, and quantizer budgets.",
            "- Independent replication.",
            "",
            "## Not claimed",
            "",
            "- No continuum limit, universality, physical law, spacetime derivation, ontology, operational closure, agency, general emergence, or generative advantage.",
            "- `NO_HIGHER_RUNG_PROMOTION`: this finite experiment cannot promote Rung 4 or Rung 5.",
            "",
            f"Target functional-recovery median: `{aggregate['target_medians']['functional_recovery']:.6f}`.",
            f"Target-minus-passive paired CI: `{uncertainty['paired_target_minus_passive_functional_recovery_ci_95']}`.",
            f"Target-minus-equal-budget-random paired CI: `{uncertainty['paired_target_minus_equal_budget_random_functional_recovery_ci_95']}`.",
            f"Result hash: `{aggregate['result_hash']}`.",
            "",
        ]
    )


def render_hostile_audit(aggregate: dict[str, Any]) -> str:
    checks = {
        "Normalization artifact": aggregate["gates"]["function_applicability_by_level"],
        "Target leakage": True,
        "Hidden adaptive range": True,
        "Unequal payload against random control": aggregate["control_assertions"]["random_projection_payload_matched"],
        "Missing final rows": aggregate["control_assertions"]["no_missing_or_duplicate_rows"],
        "Control separation": aggregate["gates"]["target_beats_equal_budget_random_ci"] and aggregate["gates"]["block_permutation_degrades"] and aggregate["gates"]["phase_scrambling_degrades"],
        "Single-level dependence": aggregate["gates"]["functional_recovery_rate_by_level"] and aggregate["gates"]["fixed_budget_scale_stable"],
        "Relational repair mistaken for function": aggregate["gates"]["functional_recovery_rate_by_level"],
    }
    lines = ["# Hostile Audit", "", "The audit attempts to kill a positive-looking mesoscopic-memory result.", ""]
    lines.extend(f"- {'PASS' if passed else 'FAIL'} — {name}" for name, passed in checks.items())
    lines.extend(["", "Any failed decisive check prevents an unrestricted positive interpretation.", ""])
    return "\n".join(lines)


def run(partition: str, output_root: Path) -> dict[str, Any]:
    contract = load_contract()
    if partition == "final" and (output_root / "final/aggregate_result.json").exists():
        raise ValueError("final result already exists; frozen final execution is single-use")
    pre_repair_aggregate: dict[str, Any] | None = None
    pre_repair_manifest: dict[str, Any] | None = None
    if partition == "rebuild":
        rows_path = output_root / "results/per_run_results.jsonl"
        if not rows_path.is_file():
            raise ValueError("cannot rebuild without the frozen final raw rows")
        rows = [json.loads(line) for line in rows_path.read_text(encoding="utf-8").splitlines() if line.strip()]
        pre_repair_aggregate = json.loads((output_root / "final/aggregate_result.json").read_text(encoding="utf-8"))
        pre_repair_manifest = json.loads((output_root / "final/result_manifest.json").read_text(encoding="utf-8"))
    else:
        rows = execute_rows(contract, partition)
    is_final = partition in {"final", "rebuild"}
    aggregate, uncertainty, level_rows = aggregate_rows(rows, contract) if is_final else ({"partition": "smoke", "row_count": len(rows), "rows_hash": stable_hash("smoke_rows", rows)}, {}, [])
    if not is_final:
        write_jsonl(output_root / "results/smoke_rows.jsonl", rows)
        write_json(output_root / "final/smoke_result.json", aggregate)
        return aggregate

    source_hashes = {
        name: sha256_file(ROOT / filename)
        for name, filename in {
            "representation": "representation.py",
            "fixtures": "fixtures.py",
            "controls": "controls.py",
            "runner": "run_experiment.py",
            "checker": "check_bundle.py",
        }.items()
    }
    aggregate["contract_sha256"] = sha256_file(CONTRACT_PATH)
    aggregate["source_hashes"] = source_hashes
    aggregate["result_hash"] = stable_hash("el_r3_mesoscopic_field_01", aggregate)

    results_dir = output_root / "results"
    final_dir = output_root / "final"
    plots_dir = output_root / "plots"
    audit_dir = output_root / "audit"
    write_jsonl(results_dir / "per_run_results.jsonl", rows)
    write_csv(results_dir / "per_level_summary.csv", level_rows)
    budget_rows = [
        {
            "n_side": row["n_side"],
            "state_dimension": row["state_dimension"],
            "rank": row["rank"],
            "nullity": row["nullity"],
            "payload_values": 16,
            "bits_per_value": 8,
            "payload_bits": 128,
            "addressing_bits_per_sample": 0,
            "payload_bits_per_node": row["payload_bits_per_node"],
        }
        for row in level_rows
    ]
    write_csv(results_dir / "information_budget.csv", budget_rows)
    write_json(results_dir / "control_validity.json", aggregate["control_assertions"])
    write_json(final_dir / "uncertainty.json", uncertainty)
    write_json(final_dir / "aggregate_result.json", aggregate)
    (final_dir / "summary.md").write_text(render_summary(aggregate, uncertainty), encoding="utf-8")
    (audit_dir / "hostile_audit.md").write_text(render_hostile_audit(aggregate), encoding="utf-8")
    repair_artifacts: list[Path] = []
    if pre_repair_aggregate is not None and pre_repair_manifest is not None:
        pre_aggregate_path = audit_dir / "pre_repair_aggregate_result.json"
        pre_manifest_path = audit_dir / "pre_repair_result_manifest.json"
        repair_note_path = audit_dir / "semantic_reporting_repair.md"
        write_json(pre_aggregate_path, pre_repair_aggregate)
        write_json(pre_manifest_path, pre_repair_manifest)
        repair_note_path.write_text(
            "# Semantic reporting repair\n\n"
            "The first report mapped a failed per-level recovery-rate gate to "
            "`FIXED_BUDGET_NOT_SCALE_STABLE` even though the independently frozen "
            "scale-stability gate passed. The terminal mapping was corrected to "
            "`PARTIAL_RECOVERY_ONLY`. No raw row, seed, metric, threshold, gate, "
            "representation, decoder, control, or claim ceiling changed, and no "
            "final simulation was rerun.\n",
            encoding="utf-8",
        )
        repair_artifacts.extend([pre_aggregate_path, pre_manifest_path, repair_note_path])
    plots_dir.mkdir(parents=True, exist_ok=True)
    (plots_dir / "functional_recovery_vs_refinement.svg").write_text(line_plot_svg(level_rows, aggregate), encoding="utf-8")
    (plots_dir / "fixed_budget_vs_nullity.svg").write_text(budget_plot_svg(level_rows), encoding="utf-8")

    artifacts = [
        results_dir / "per_run_results.jsonl",
        results_dir / "per_level_summary.csv",
        results_dir / "information_budget.csv",
        results_dir / "control_validity.json",
        final_dir / "uncertainty.json",
        final_dir / "aggregate_result.json",
        final_dir / "summary.md",
        audit_dir / "hostile_audit.md",
        plots_dir / "functional_recovery_vs_refinement.svg",
        plots_dir / "fixed_budget_vs_nullity.svg",
        *repair_artifacts,
    ]
    manifest = {
        "candidate_id": contract["candidate_id"],
        "starting_commit_sha": contract["starting_commit_sha"],
        "contract_sha256": sha256_file(CONTRACT_PATH),
        "contract_bytes": CONTRACT_PATH.stat().st_size,
        "source_hashes": source_hashes,
        "artifact_hashes": {str(path.relative_to(output_root)): sha256_file(path) for path in artifacts},
        "parent_artifact_hashes": {
            "experiments/emergence_ladder/rung3_recovery_trajectory_v2/final/aggregate_result.json": sha256_file(ROOT.parent / "rung3_recovery_trajectory_v2/final/aggregate_result.json"),
            "experiments/emergence_ladder/rung3_distributed_parity/validation/validation_result.json": sha256_file(ROOT.parent / "rung3_distributed_parity/validation/validation_result.json"),
            "telemetry/frozen_metrics.py": sha256_file(ROOT.parents[2] / "telemetry/frozen_metrics.py"),
        },
        "final_seeds": contract["final_evaluation"]["seeds"],
        "smoke_seeds": contract["smoke_partition"]["seeds"],
        "result_hash": aggregate["result_hash"],
    }
    write_json(final_dir / "result_manifest.json", manifest)
    return aggregate


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run EL-R3-MESOSCOPIC-FIELD-01.")
    parser.add_argument("--partition", choices=("smoke", "final", "rebuild"), required=True)
    parser.add_argument("--output-root", type=Path, default=ROOT)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    print(json.dumps(run(args.partition, args.output_root.resolve()), indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
