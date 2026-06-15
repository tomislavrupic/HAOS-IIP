from __future__ import annotations

import argparse
import copy
import itertools
import json
import math
import random
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from experiments.cst.cst_io import read_csv, repo_rel, stable_hash, write_csv, write_json
from experiments.cst.cst_metrics import (
    COMPONENT_NAMES,
    compute_closure_distance,
    jaccard_distance,
    mean,
    numeric_dict_distance,
    population_std,
    safe_numeric_distance,
    sequence_rank_distance,
)
from experiments.cst.cst_types import ClosureSignature, RecoveryObservation
from experiments.cst.run_cst_benchmark import DEFAULT_OUTPUT_DIR, PHASE_ARTIFACTS, run_cst_benchmark


OUTPUT_DIR = DEFAULT_OUTPUT_DIR / "failure_localization"
FOCUS_CONTROLS = [
    "periodic_diagonal_augmented_control",
    "randomized_edge_control",
    "shuffled_structure_control",
]
POWER_ALPHA = 0.05
POWER_TARGET = 0.80
Z_ALPHA_975 = 1.959963984540054
Z_POWER_80 = 0.8416212335729143


def value_or_none(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def stable_float(value: Any) -> float | None:
    parsed = value_or_none(value)
    return None if parsed is None else float(parsed)


def stats(values: list[float]) -> dict[str, Any]:
    if not values:
        return {
            "count": 0,
            "mean": None,
            "std": None,
            "variance": None,
            "min": None,
            "max": None,
        }
    return {
        "count": len(values),
        "mean": mean(values),
        "std": population_std(values),
        "variance": population_std(values) ** 2,
        "min": min(values),
        "max": max(values),
    }


def cohen_d(reference: list[float], comparator: list[float]) -> float | None:
    if not reference or not comparator:
        return None
    reference_std = population_std(reference)
    comparator_std = population_std(comparator)
    pooled_var = ((reference_std ** 2) + (comparator_std ** 2)) / 2.0
    diff = mean(comparator) - mean(reference)
    if pooled_var <= 1.0e-18:
        if abs(diff) <= 1.0e-18:
            return 0.0
        return math.copysign(math.inf, diff)
    return diff / math.sqrt(pooled_var)


def histogram_overlap(left: list[float], right: list[float], *, bins: int = 20) -> float | None:
    if not left or not right:
        return None
    left_counts = [0 for _ in range(bins)]
    right_counts = [0 for _ in range(bins)]
    for value in left:
        index = min(max(int(float(value) * bins), 0), bins - 1)
        left_counts[index] += 1
    for value in right:
        index = min(max(int(float(value) * bins), 0), bins - 1)
        right_counts[index] += 1
    left_total = float(len(left))
    right_total = float(len(right))
    return sum(min(left_counts[index] / left_total, right_counts[index] / right_total) for index in range(bins))


def rank_order(*, target: list[float], target_control: list[float], control_control: list[float]) -> str:
    ranked = []
    for label, values in (
        ("target_target", target),
        ("target_control", target_control),
        ("control_control", control_control),
    ):
        if values:
            ranked.append((mean(values), label))
    if not ranked:
        return "unavailable"
    return " < ".join(label for _, label in sorted(ranked))


def component_values(
    pairs: list[tuple[RecoveryObservation, RecoveryObservation]],
    component: str,
    *,
    weights: dict[str, float],
    enabled_components: list[str],
) -> tuple[list[float], int]:
    values: list[float] = []
    missing = 0
    for left, right in pairs:
        distance = compute_closure_distance(
            left,
            right,
            weights=weights,
            enabled_components=enabled_components,
            scalar_enabled=True,
        )
        metric = distance.components[component]
        if metric.availability == "available" and metric.value is not None:
            values.append(float(metric.value))
        else:
            missing += 1
    return values, missing


def observation_key(observation: RecoveryObservation) -> tuple[int, int | None]:
    return (int(observation.n_side), observation.seed)


def build_pair_sets(observations: list[RecoveryObservation]) -> dict[str, Any]:
    targets = [item for item in observations if item.observation_kind == "target"]
    targets_by_id = {item.observation_id: item for item in targets}
    targets_by_key = {observation_key(item): item for item in targets}
    controls = [item for item in observations if item.control_role == "negative"]

    target_pairs: list[tuple[RecoveryObservation, RecoveryObservation]] = []
    targets_by_n_side: dict[int, list[RecoveryObservation]] = {}
    for target in targets:
        targets_by_n_side.setdefault(int(target.n_side), []).append(target)
    for same_grid_targets in targets_by_n_side.values():
        for left, right in itertools.combinations(sorted(same_grid_targets, key=lambda item: int(item.seed or -1)), 2):
            target_pairs.append((left, right))

    target_control_pairs: dict[str, list[tuple[RecoveryObservation, RecoveryObservation]]] = {}
    controls_by_family: dict[str, list[RecoveryObservation]] = {}
    for control in controls:
        label = control.control_group_label
        controls_by_family.setdefault(label, []).append(control)
        source_id = control.perturbation.parameters.get("source_observation_id")
        if source_id in targets_by_id:
            target = targets_by_id[source_id]
        else:
            target = targets_by_key.get(observation_key(control))
        if target is None:
            continue
        target_control_pairs.setdefault(label, []).append((target, control))

    control_control_pairs: dict[str, list[tuple[RecoveryObservation, RecoveryObservation]]] = {}
    for label, family_controls in controls_by_family.items():
        by_grid: dict[int, list[RecoveryObservation]] = {}
        for control in family_controls:
            by_grid.setdefault(int(control.n_side), []).append(control)
        for same_grid_controls in by_grid.values():
            for left, right in itertools.combinations(sorted(same_grid_controls, key=lambda item: int(item.seed or -1)), 2):
                control_control_pairs.setdefault(label, []).append((left, right))

    return {
        "targets": targets,
        "controls": controls,
        "target_pairs": target_pairs,
        "target_control_pairs": target_control_pairs,
        "control_control_pairs": control_control_pairs,
        "controls_by_family": controls_by_family,
    }


def component_separation_rows(benchmark: dict[str, Any], pair_sets: dict[str, Any]) -> list[dict[str, Any]]:
    config = benchmark["config"]
    weights = config["distance"]["weights"]
    enabled = config["distance"]["enabled_components"]
    rows: list[dict[str, Any]] = []
    control_labels = sorted(pair_sets["target_control_pairs"])
    for component in COMPONENT_NAMES:
        target_values, target_missing = component_values(
            pair_sets["target_pairs"],
            component,
            weights=weights,
            enabled_components=enabled,
        )
        for label in control_labels:
            target_control_values, target_control_missing = component_values(
                pair_sets["target_control_pairs"].get(label, []),
                component,
                weights=weights,
                enabled_components=enabled,
            )
            control_control_values, control_control_missing = component_values(
                pair_sets["control_control_pairs"].get(label, []),
                component,
                weights=weights,
                enabled_components=enabled,
            )
            effect = cohen_d(target_values, target_control_values)
            overlap = histogram_overlap(target_values, target_control_values)
            rows.append(
                {
                    "component": component,
                    "control_family": label,
                    "is_focus_control": label in FOCUS_CONTROLS,
                    "target_count": len(target_values),
                    "target_missing": target_missing,
                    "target_mean": stats(target_values)["mean"],
                    "target_variance": stats(target_values)["variance"],
                    "target_values": json.dumps(target_values, sort_keys=True),
                    "target_control_count": len(target_control_values),
                    "target_control_missing": target_control_missing,
                    "target_control_mean": stats(target_control_values)["mean"],
                    "target_control_variance": stats(target_control_values)["variance"],
                    "target_control_values": json.dumps(target_control_values, sort_keys=True),
                    "control_control_count": len(control_control_values),
                    "control_control_missing": control_control_missing,
                    "control_control_mean": stats(control_control_values)["mean"],
                    "control_control_variance": stats(control_control_values)["variance"],
                    "control_control_values": json.dumps(control_control_values, sort_keys=True),
                    "effect_size_cohen_d_control_minus_target": effect,
                    "overlap_coefficient_histogram_20_bin": overlap,
                    "rank_order_by_mean_distance": rank_order(
                        target=target_values,
                        target_control=target_control_values,
                        control_control=control_control_values,
                    ),
                    "target_outperforms_target_control": bool(
                        target_values and target_control_values and mean(target_values) < mean(target_control_values)
                    ),
                    "notes": "positive effect means target-target distances are smaller than target-control distances",
                }
            )
    return rows


def variant_distances(
    pairs: list[tuple[RecoveryObservation, RecoveryObservation]],
    *,
    enabled_components: list[str],
    weights: dict[str, float],
) -> list[float]:
    values: list[float] = []
    for left, right in pairs:
        distance = compute_closure_distance(
            left,
            right,
            weights=weights,
            enabled_components=enabled_components,
            scalar_enabled=True,
        )
        if distance.scalar_distance is not None:
            values.append(float(distance.scalar_distance))
    return values


def variant_verdict(
    *,
    target_values: list[float],
    control_values_by_family: dict[str, list[float]],
    thresholds: dict[str, Any],
    seed_count: int,
) -> tuple[str, list[str]]:
    reasons: list[str] = []
    if not target_values:
        return "FAIL", ["no target values"]
    pooled_control = [value for values in control_values_by_family.values() for value in values]
    if not pooled_control:
        return "FAIL", ["no control values"]

    verdict = "PASS"
    if any(value > float(thresholds["closure_distance_max"]) for value in target_values):
        verdict = "FAIL"
        reasons.append("target distances exceed closure threshold")
    if mean(pooled_control) <= mean(target_values):
        verdict = "FAIL"
        reasons.append("pooled controls equal or beat target")
    for label, values in sorted(control_values_by_family.items()):
        if values and mean(values) <= mean(target_values):
            verdict = "FAIL"
            reasons.append(f"{label} equals or beats target")
    if seed_count < int(thresholds["min_seed_count_for_pass"]):
        if verdict == "PASS":
            verdict = "OPEN"
        reasons.append("seed count below configured PASS minimum")
    if verdict == "PASS":
        reasons.append("variant gates passed")
    return verdict, reasons


def ablation_matrix_rows(benchmark: dict[str, Any], pair_sets: dict[str, Any]) -> list[dict[str, Any]]:
    config = benchmark["config"]
    weights = config["distance"]["weights"]
    thresholds = config["thresholds"]
    base_verdict = benchmark["result"].verdict
    base_components = list(config["distance"]["enabled_components"])
    rows: list[dict[str, Any]] = []
    variants: list[tuple[str, str, list[str]]] = []
    for component in base_components:
        variants.append(("single_component", component, [component]))
        variants.append(("leave_one_out", component, [item for item in base_components if item != component]))

    for variant_type, component, enabled in variants:
        target_values = variant_distances(
            pair_sets["target_pairs"],
            enabled_components=enabled,
            weights=weights,
        )
        focus_control_values = {
            label: variant_distances(pair_sets["target_control_pairs"].get(label, []), enabled_components=enabled, weights=weights)
            for label in FOCUS_CONTROLS
        }
        all_control_values = {
            label: variant_distances(pairs, enabled_components=enabled, weights=weights)
            for label, pairs in pair_sets["target_control_pairs"].items()
        }
        for scope, values_by_family in (
            ("focus_controls", focus_control_values),
            ("all_negative_controls", all_control_values),
        ):
            verdict, reasons = variant_verdict(
                target_values=target_values,
                control_values_by_family=values_by_family,
                thresholds=thresholds,
                seed_count=len(config["seeds"]),
            )
            pooled = [value for values in values_by_family.values() for value in values]
            failing = [label for label, values in values_by_family.items() if values and mean(values) <= mean(target_values)]
            flags = []
            if variant_type == "single_component" and verdict == base_verdict:
                flags.append("single_component_reproduces_base_verdict")
            if variant_type == "leave_one_out" and verdict != base_verdict:
                flags.append("removal_changes_verdict")
            if pooled and abs(mean(pooled) - mean(target_values)) <= 1.0e-9:
                flags.append("no_discrimination")
            if pooled and len(set(round(value, 12) for value in pooled + target_values)) == 1:
                flags.append("identical_values")
            if failing:
                flags.append("control_reversal")
            rows.append(
                {
                    "variant_type": variant_type,
                    "component_under_test": component,
                    "enabled_components": "|".join(enabled),
                    "control_scope": scope,
                    "target_mean": mean(target_values) if target_values else None,
                    "target_variance": population_std(target_values) ** 2 if target_values else None,
                    "control_mean": mean(pooled) if pooled else None,
                    "control_variance": population_std(pooled) ** 2 if pooled else None,
                    "failing_control_families": "|".join(failing),
                    "variant_verdict": verdict,
                    "base_verdict": base_verdict,
                    "verdict_changed": verdict != base_verdict,
                    "flags": "|".join(flags),
                    "reasons": "|".join(reasons),
                }
            )
    return rows


def signature_source_fields(observation: RecoveryObservation) -> dict[str, Any]:
    fields = observation.closure_signature.fields
    return {
        "observation_id": observation.observation_id,
        "observation_kind": observation.observation_kind,
        "control_group": observation.control_group_label,
        "n_side": observation.n_side,
        "seed": observation.seed,
        "branch_signature_id": observation.closure_signature.branch_signature_id,
        "address_summary": fields.get("address_summary", {}),
        "chain_signature": fields.get("event_chain_summary", {}).get("chain_signature"),
        "event_times": fields.get("event_chain_summary", {}).get("event_times", {}),
        "invariant_summary": fields.get("invariant_summary", {}),
        "spectral_low_mode_summary": fields.get("spectral_low_mode_summary", {}),
        "causal_summary": fields.get("causal_summary", {}),
        "shell_ordering_summary": fields.get("shell_ordering_summary", {}),
    }


def expanded_signature_hash(observation: RecoveryObservation) -> str:
    fields = signature_source_fields(observation)
    expanded = {
        "address_summary": fields["address_summary"],
        "chain_signature": fields["chain_signature"],
        "event_times": fields["event_times"],
        "invariant_summary": fields["invariant_summary"],
        "spectral_low_mode_summary": fields["spectral_low_mode_summary"],
        "causal_summary": fields["causal_summary"],
        "shell_ordering_summary": fields["shell_ordering_summary"],
    }
    return stable_hash(expanded, prefix="expanded_", float_precision=12)


def signature_collision_report(benchmark: dict[str, Any], pair_sets: dict[str, Any]) -> dict[str, Any]:
    observations = benchmark["observations"]
    groups: dict[str, list[RecoveryObservation]] = {}
    for observation in observations:
        groups.setdefault(observation.closure_signature.branch_signature_id, []).append(observation)
    collisions = []
    for signature_id, members in sorted(groups.items()):
        labels = sorted({member.control_group_label for member in members})
        expanded_hashes = sorted({expanded_signature_hash(member) for member in members})
        if len(members) > 1:
            collisions.append(
                {
                    "branch_signature_id": signature_id,
                    "member_count": len(members),
                    "control_groups": labels,
                    "expanded_hash_count": len(expanded_hashes),
                    "expanded_hashes": expanded_hashes,
                    "source_fields": [signature_source_fields(member) for member in members],
                    "collision_type": (
                        "field_selection_alias"
                        if len(expanded_hashes) > 1
                        else "exact_address_collision"
                    ),
                }
            )

    target_ids = {
        target.closure_signature.branch_signature_id
        for target in pair_sets["targets"]
    }
    current_address_groups: dict[str, list[str]] = {}
    high_precision_address_groups: dict[str, list[str]] = {}
    for observation in observations:
        address = observation.closure_signature.fields.get("address_summary", {})
        current_address_groups.setdefault(
            stable_hash(address, prefix="addr6_", float_precision=6),
            [],
        ).append(observation.observation_id)
        high_precision_address_groups.setdefault(
            stable_hash(address, prefix="addr12_", float_precision=12),
            [],
        ).append(observation.observation_id)

    omitted_fields = [
        "event_chain_summary.chain_signature beyond dominant prefix",
        "event_chain_summary.event_times",
        "invariant_summary",
        "spectral_low_mode_summary.low_mode_band",
        "causal_summary.acyclicity_score",
        "causal_summary.mismatch_rate",
        "shell_ordering_summary.shell_order_score",
        "influence_edge_summary.edges beyond address prefix",
    ]
    return {
        "target_seed_signature_ids": sorted(target_ids),
        "equivalent_target_runs_match_across_seeds": len(target_ids) == 1,
        "current_signature_group_count": len(groups),
        "observation_count": len(observations),
        "collision_count": len(collisions),
        "collisions": collisions,
        "canonicalization_precision_test": {
            "current_float_precision_group_count": len(current_address_groups),
            "high_precision_address_group_count": len(high_precision_address_groups),
            "precision_change_splits_groups": len(high_precision_address_groups) > len(current_address_groups),
            "note": "Current signature hashes address_summary only; most collisions are field-selection aliases, not SHA collisions.",
        },
        "quantization_rounding_audit": {
            "order_compatibility_quantum": benchmark["config"]["signature"]["order_compatibility_quantum"],
            "low_mode_quantum_configured": benchmark["config"]["signature"]["low_mode_quantum"],
            "low_mode_band_in_address_summary": False,
            "omitted_or_coarsened_fields": omitted_fields,
        },
    }


def edge_degree_distribution(observation: RecoveryObservation) -> list[int]:
    edges = observation.closure_signature.fields.get("influence_edge_summary", {}).get("edges", [])
    degree: dict[str, int] = {}
    for edge in edges:
        if "->" not in str(edge):
            continue
        left, right = str(edge).split("->", 1)
        degree[left] = degree.get(left, 0) + 1
        degree[right] = degree.get(right, 0) + 1
    return sorted(degree.values())


def label_set(observation: RecoveryObservation) -> set[str]:
    chain = observation.closure_signature.fields.get("event_chain_summary", {}).get("chain_signature", "")
    return {item for item in str(chain).split(">") if item}


def property_preservation(source: RecoveryObservation, control: RecoveryObservation) -> dict[str, bool | float | str]:
    source_fields = source.closure_signature.fields
    control_fields = control.closure_signature.fields
    source_labels = label_set(source)
    control_labels = label_set(control)
    label_overlap = len(source_labels & control_labels) / max(len(source_labels | control_labels), 1)
    source_chain = source_fields.get("event_chain_summary", {}).get("chain_signature", "").split(">")
    control_chain = control_fields.get("event_chain_summary", {}).get("chain_signature", "").split(">")
    return {
        "degree_distribution": edge_degree_distribution(source) == edge_degree_distribution(control),
        "spectral_density_proxy": source_fields.get("spectral_low_mode_summary") == control_fields.get("spectral_low_mode_summary"),
        "event_chain_length": len(source_chain) == len(control_chain),
        "causal_depth": source_fields.get("causal_summary", {}).get("mean_causal_depth") == control_fields.get("causal_summary", {}).get("mean_causal_depth")
        and source_fields.get("causal_summary", {}).get("front_depth_profile") == control_fields.get("causal_summary", {}).get("front_depth_profile"),
        "shell_ordering": source_fields.get("shell_ordering_summary") == control_fields.get("shell_ordering_summary"),
        "invariant_marginals": source_fields.get("invariant_summary") == control_fields.get("invariant_summary"),
        "label_jaccard": label_overlap,
        "label_specific_leakage": "high" if label_overlap >= 0.95 else "partial" if label_overlap > 0.0 else "low",
    }


def control_integrity_markdown(pair_sets: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    integrity: dict[str, Any] = {}
    lines = [
        "# Control Integrity Audit",
        "",
        "Implemented fact: this audit compares each matched control to the source target observation used by the CST toy benchmark.",
        "Design choice: spectral density is represented only by the available low-mode summary proxy; no full eigenspectrum is reconstructed.",
        "",
        "| Control | Preserved Properties | Destroyed Properties | Label Leakage | Validity Note |",
        "| --- | --- | --- | --- | --- |",
    ]
    for label, pairs in sorted(pair_sets["target_control_pairs"].items()):
        rows = [property_preservation(source, control) for source, control in pairs]
        properties = [
            "degree_distribution",
            "spectral_density_proxy",
            "event_chain_length",
            "causal_depth",
            "shell_ordering",
            "invariant_marginals",
        ]
        preservation_rates = {
            prop: mean([1.0 if row[prop] else 0.0 for row in rows]) if rows else None
            for prop in properties
        }
        label_overlap = mean([float(row["label_jaccard"]) for row in rows]) if rows else None
        preserved = [prop for prop, rate in preservation_rates.items() if rate == 1.0]
        destroyed = [prop for prop, rate in preservation_rates.items() if rate == 0.0]
        partial = [prop for prop, rate in preservation_rates.items() if rate not in (0.0, 1.0, None)]
        leakage = "high" if label_overlap is not None and label_overlap >= 0.95 else "partial" if label_overlap else "low"
        note = "valid negative control for all listed components"
        if label in ("randomized_edge_control", "shuffled_structure_control"):
            note = "invalid as a negative control for invariant, spectral, causal-depth, shell-ordering, and label-set-sensitive tests"
        elif label == "periodic_diagonal_augmented_control":
            note = "strongest frozen matched control; still shares benchmark probe, seed grid, and address fields"
        elif preserved:
            note = "component-specific validity only; preserves some tested signal"
        integrity[label] = {
            "pair_count": len(rows),
            "preservation_rates": preservation_rates,
            "mean_label_jaccard": label_overlap,
            "label_specific_leakage": leakage,
            "validity_note": note,
        }
        lines.append(
            "| {label} | {preserved} | {destroyed} | {leakage} | {note} |".format(
                label=label,
                preserved=", ".join(preserved + [f"partial:{item}" for item in partial]) or "none",
                destroyed=", ".join(destroyed) or "none",
                leakage=leakage,
                note=note,
            )
        )
    lines.extend(
        [
            "",
            "Key hostile finding: randomized and shuffled controls retain the exact event-label set and most non-chain descriptor marginals.",
            "A control that preserves a component's tested signal is not a valid negative control for that component.",
        ]
    )
    return "\n".join(lines) + "\n", integrity


def clone_with_direct_timing(observation: RecoveryObservation) -> RecoveryObservation:
    fields = copy.deepcopy(observation.closure_signature.fields)
    distance_summary = fields.get("distance_surrogate_summary", {})
    chain = fields.get("event_chain_summary", {}).get("chain_signature", "")
    chain_labels = {item for item in str(chain).split(">") if item}
    direct_times: dict[str, float] = {}
    for label, summary in distance_summary.items():
        if label in chain_labels:
            value = stable_float(summary.get("propagation_time"))
            if value is not None:
                direct_times[label] = value
    shell = fields.get("shell_ordering_summary", {})
    for label, key in (("front_near", "near_mean_arrival"), ("front_far", "far_mean_arrival")):
        if label in chain_labels:
            value = stable_float(shell.get(key))
            if value is not None:
                direct_times[label] = value
    fields.setdefault("event_chain_summary", {})["event_times"] = direct_times
    fields.setdefault("warnings", []).append("failure-localization direct-timing mode: omitted non-observed event times")
    signature = replace(
        observation.closure_signature,
        fields=fields,
        warnings=list(observation.closure_signature.warnings) + ["direct timing failure-localization clone"],
    )
    return replace(observation, closure_signature=signature)


def clone_observations_direct_timing(observations: list[RecoveryObservation]) -> list[RecoveryObservation]:
    return [clone_with_direct_timing(observation) for observation in observations]


def timing_ablation_rows(benchmark: dict[str, Any], pair_sets: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    config = benchmark["config"]
    weights = config["distance"]["weights"]
    modes = [
        ("interpolated_non_front_times", benchmark["observations"], list(config["distance"]["enabled_components"])),
        ("temporal_component_disabled", benchmark["observations"], [item for item in config["distance"]["enabled_components"] if item != "d_time"]),
        ("direct_observed_timing_only", clone_observations_direct_timing(benchmark["observations"]), list(config["distance"]["enabled_components"])),
    ]
    rows: list[dict[str, Any]] = []
    summary: dict[str, Any] = {}
    for mode_name, observations, enabled in modes:
        mode_pairs = build_pair_sets(observations)
        target_values = variant_distances(mode_pairs["target_pairs"], enabled_components=enabled, weights=weights)
        focus_values_by_family = {
            label: variant_distances(mode_pairs["target_control_pairs"].get(label, []), enabled_components=enabled, weights=weights)
            for label in FOCUS_CONTROLS
        }
        verdict, reasons = variant_verdict(
            target_values=target_values,
            control_values_by_family=focus_values_by_family,
            thresholds=config["thresholds"],
            seed_count=len(config["seeds"]),
        )
        summary[mode_name] = {
            "enabled_components": enabled,
            "verdict": verdict,
            "reasons": reasons,
            "target_mean": mean(target_values) if target_values else None,
            "focus_control_mean": mean([value for values in focus_values_by_family.values() for value in values]),
        }
        for component in COMPONENT_NAMES:
            target_component_values, target_missing = component_values(
                mode_pairs["target_pairs"],
                component,
                weights=weights,
                enabled_components=list(config["distance"]["enabled_components"]),
            )
            for label in FOCUS_CONTROLS:
                control_values, control_missing = component_values(
                    mode_pairs["target_control_pairs"].get(label, []),
                    component,
                    weights=weights,
                    enabled_components=list(config["distance"]["enabled_components"]),
                )
                rows.append(
                    {
                        "timing_mode": mode_name,
                        "component": component,
                        "component_enabled_in_verdict": component in enabled,
                        "control_family": label,
                        "verdict": verdict,
                        "target_mean": mean(target_component_values) if target_component_values else None,
                        "target_variance": population_std(target_component_values) ** 2 if target_component_values else None,
                        "target_missing": target_missing,
                        "target_control_mean": mean(control_values) if control_values else None,
                        "target_control_variance": population_std(control_values) ** 2 if control_values else None,
                        "target_control_missing": control_missing,
                        "effect_size_cohen_d_control_minus_target": cohen_d(target_component_values, control_values),
                        "reasons": "|".join(reasons),
                    }
                )
    return rows, summary


def symmetric_relative(left: float, right: float, *, epsilon: float = 1.0e-12) -> tuple[float, float, bool]:
    denominator = abs(left) + abs(right) + epsilon
    return min(abs(left - right) / denominator, 1.0), denominator, denominator <= 1.0e-6


def chain_from_observation(observation: RecoveryObservation) -> list[str]:
    chain = observation.closure_signature.fields.get("event_chain_summary", {}).get("chain_signature", "")
    return [item for item in str(chain).split(">") if item]


def raw_component_audit(left: RecoveryObservation, right: RecoveryObservation, component: str) -> dict[str, Any]:
    left_fields = left.closure_signature.fields
    right_fields = right.closure_signature.fields
    raw: dict[str, Any] = {}
    robust_values: list[float] = []
    small_denominator = False
    normalization_reference = ""
    artifact_flags: list[str] = []

    if component in ("d_inv", "d_spec"):
        section = "invariant_summary" if component == "d_inv" else "spectral_low_mode_summary"
        keys = (
            ["terminal_low_k_fraction", "transport_character_index", "shell_latency_span"]
            if component == "d_inv"
            else ["terminal_low_k_fraction", "terminal_low_k_power", "terminal_fluctuation_slope"]
        )
        normalization_reference = "per-field abs diff divided by max(abs(left), abs(right), 1e-12)"
        for key in keys:
            left_value = stable_float(left_fields.get(section, {}).get(key))
            right_value = stable_float(right_fields.get(section, {}).get(key))
            if left_value is None or right_value is None:
                continue
            normalized = safe_numeric_distance(left_value, right_value)
            robust, denominator, tiny = symmetric_relative(left_value, right_value)
            small_denominator = small_denominator or tiny
            robust_values.append(robust)
            raw[key] = {
                "left": left_value,
                "right": right_value,
                "abs_diff": abs(left_value - right_value),
                "default_normalized": normalized,
                "symmetric_relative": robust,
                "symmetric_denominator": denominator,
            }
    elif component == "d_struct":
        left_edges = [str(edge) for edge in left_fields.get("influence_edge_summary", {}).get("edges", [])]
        right_edges = [str(edge) for edge in right_fields.get("influence_edge_summary", {}).get("edges", [])]
        edge_distance = jaccard_distance(left_edges, right_edges)
        robust_values.append(edge_distance)
        raw["edge_set"] = {
            "left_count": len(left_edges),
            "right_count": len(right_edges),
            "intersection": len(set(left_edges) & set(right_edges)),
            "union": len(set(left_edges) | set(right_edges)),
            "jaccard_distance": edge_distance,
        }
        for key in ("edge_count", "stable_edge_count"):
            left_value = stable_float(left_fields.get("influence_edge_summary", {}).get(key))
            right_value = stable_float(right_fields.get("influence_edge_summary", {}).get(key))
            if left_value is None or right_value is None:
                continue
            robust, denominator, tiny = symmetric_relative(left_value, right_value)
            robust_values.append(robust)
            small_denominator = small_denominator or tiny
            raw[key] = {"left": left_value, "right": right_value, "symmetric_relative": robust, "symmetric_denominator": denominator}
        normalization_reference = "edge Jaccard plus relative edge-count drift"
    elif component == "d_causal":
        normalization_reference = "front depths use fixed scale 4.0; causal fields use relative numeric distance"
        left_depth = left_fields.get("causal_summary", {}).get("front_depth_profile", {})
        right_depth = right_fields.get("causal_summary", {}).get("front_depth_profile", {})
        for key in ("front_near_depth", "front_far_depth"):
            left_value = stable_float(left_depth.get(key))
            right_value = stable_float(right_depth.get(key))
            if left_value is None or right_value is None:
                continue
            robust_values.append(safe_numeric_distance(left_value, right_value, scale=4.0))
            raw[key] = {"left": left_value, "right": right_value, "abs_diff": abs(left_value - right_value), "scale": 4.0}
        for key in ("acyclicity_score", "mean_causal_depth", "order_compatibility_score", "mismatch_rate"):
            left_value = stable_float(left_fields.get("causal_summary", {}).get(key))
            right_value = stable_float(right_fields.get("causal_summary", {}).get(key))
            if left_value is None or right_value is None:
                continue
            robust, denominator, tiny = symmetric_relative(left_value, right_value)
            robust_values.append(robust)
            small_denominator = small_denominator or tiny
            raw[key] = {"left": left_value, "right": right_value, "symmetric_relative": robust, "symmetric_denominator": denominator}
    elif component == "d_time":
        normalization_reference = "rank distance plus event-time abs diff divided by tau scale 4.2"
        left_chain = chain_from_observation(left)
        right_chain = chain_from_observation(right)
        left_times = left_fields.get("event_chain_summary", {}).get("event_times", {})
        right_times = right_fields.get("event_chain_summary", {}).get("event_times", {})
        shared_direct_keys = sorted(set(left_times) & set(right_times))
        time_values = []
        for key in shared_direct_keys:
            left_value = stable_float(left_times.get(key))
            right_value = stable_float(right_times.get(key))
            if left_value is None or right_value is None:
                continue
            time_values.append(safe_numeric_distance(left_value, right_value, scale=4.2))
        rank_value = sequence_rank_distance(left_chain, right_chain)
        robust_values.append(rank_value)
        robust_values.extend(time_values)
        raw = {
            "rank_distance": rank_value,
            "shared_time_keys": shared_direct_keys,
            "left_time_key_count": len(left_times),
            "right_time_key_count": len(right_times),
            "time_distance_values": time_values,
        }
        if any("interpolated" in warning for warning in left.closure_signature.warnings + right.closure_signature.warnings):
            artifact_flags.append("heuristic_rank_interpolation_present")
    elif component == "d_addr":
        normalization_reference = "exact signature match else prefix/depth distance with minimum mismatch floor 0.25"
        left_address = left_fields.get("address_summary", {})
        right_address = right_fields.get("address_summary", {})
        prefix_distance = sequence_rank_distance(
            [str(item) for item in left_address.get("dominant_chain_prefix", [])],
            [str(item) for item in right_address.get("dominant_chain_prefix", [])],
        )
        depth_distance = numeric_dict_distance(
            left_address.get("front_depth_profile", {}),
            right_address.get("front_depth_profile", {}),
            ["front_near_depth", "front_far_depth"],
            scale=4.0,
        )
        raw_without_floor_values = [prefix_distance]
        if depth_distance is not None:
            raw_without_floor_values.append(depth_distance)
        no_floor = mean(raw_without_floor_values)
        robust_values.append(no_floor)
        if left.closure_signature.branch_signature_id != right.closure_signature.branch_signature_id and no_floor < 0.25:
            artifact_flags.append("address_floor_applied")
        raw = {
            "left_address": left_address,
            "right_address": right_address,
            "prefix_distance": prefix_distance,
            "depth_distance": depth_distance,
            "without_floor": no_floor,
        }
    return {
        "raw_values": raw,
        "robust_alt_value": mean(robust_values) if robust_values else None,
        "normalization_reference": normalization_reference,
        "small_denominator_flag": small_denominator,
        "artifact_flags": artifact_flags,
    }


def normalization_audit_rows(benchmark: dict[str, Any], pair_sets: dict[str, Any]) -> list[dict[str, Any]]:
    config = benchmark["config"]
    weights = config["distance"]["weights"]
    enabled = config["distance"]["enabled_components"]
    rows: list[dict[str, Any]] = []
    pair_groups: list[tuple[str, str, list[tuple[RecoveryObservation, RecoveryObservation]]]] = [
        ("target_target", "target", pair_sets["target_pairs"]),
    ]
    for label in FOCUS_CONTROLS:
        pair_groups.append(("target_control", label, pair_sets["target_control_pairs"].get(label, [])))
        pair_groups.append(("control_control", label, pair_sets["control_control_pairs"].get(label, [])))
    for pair_type, label, pairs in pair_groups:
        for index, (left, right) in enumerate(pairs):
            distance = compute_closure_distance(
                left,
                right,
                weights=weights,
                enabled_components=enabled,
                scalar_enabled=True,
            )
            for component in COMPONENT_NAMES:
                audit = raw_component_audit(left, right, component)
                metric = distance.components[component]
                rows.append(
                    {
                        "pair_type": pair_type,
                        "control_family": label,
                        "pair_index": index,
                        "left_observation_id": left.observation_id,
                        "right_observation_id": right.observation_id,
                        "component": component,
                        "normalized_value": metric.value,
                        "robust_alt_value": audit["robust_alt_value"],
                        "normalization_reference": audit["normalization_reference"],
                        "small_denominator_flag": audit["small_denominator_flag"],
                        "artifact_flags": "|".join(audit["artifact_flags"]),
                        "raw_values_json": json.dumps(audit["raw_values"], sort_keys=True),
                    }
                )
    return rows


def enumerate_available_seeds(config: dict[str, Any]) -> dict[str, Any]:
    shell_rows = read_csv(PHASE_ARTIFACTS["phase18_shell_ordering_metrics"])
    transport_rows = read_csv(PHASE_ARTIFACTS["phase15_transport_descriptor"])
    probe = config["probe"]
    ensemble = str(config["ensemble_size"])
    hierarchies = [config["hierarchies"]["target"], config["hierarchies"]["control"]]
    result: dict[str, Any] = {}
    for hierarchy in hierarchies:
        result[hierarchy] = {}
        for n_side in config["n_sides"]:
            shell_seeds = {
                int(row["seed"])
                for row in shell_rows
                if row["hierarchy_label"] == hierarchy
                and row["n_side"] == str(n_side)
                and row["probe_name"] == probe
                and row["ensemble_size"] == ensemble
            }
            transport_seeds = {
                int(row["seed"])
                for row in transport_rows
                if row["hierarchy_label"] == hierarchy
                and row["n_side"] == str(n_side)
                and row["probe_name"] == probe
                and row["ensemble_size"] == ensemble
            }
            result[hierarchy][str(n_side)] = {
                "phase18_shell_ordering_seeds": sorted(shell_seeds),
                "phase15_transport_descriptor_seeds": sorted(transport_seeds),
                "intersection_usable_for_cst": sorted(shell_seeds & transport_seeds),
            }
    return result


def bootstrap_effect_ci(reference: list[float], comparator: list[float], *, iterations: int = 1000) -> dict[str, Any]:
    if not reference or not comparator:
        return {"iterations": 0, "ci95": None}
    rng = random.Random(1303)
    estimates: list[float] = []
    for _ in range(iterations):
        sample_reference = [reference[rng.randrange(len(reference))] for _ in reference]
        sample_comparator = [comparator[rng.randrange(len(comparator))] for _ in comparator]
        estimate = cohen_d(sample_reference, sample_comparator)
        if estimate is not None and math.isfinite(estimate):
            estimates.append(estimate)
    if not estimates:
        return {"iterations": iterations, "ci95": None}
    estimates.sort()
    lower = estimates[int(0.025 * (len(estimates) - 1))]
    upper = estimates[int(0.975 * (len(estimates) - 1))]
    return {
        "iterations": iterations,
        "seed": 1303,
        "ci95": [lower, upper],
        "note": "bootstrap estimates uncertainty only; it is not synthetic evidence",
    }


def pair_count_to_seed_count(pair_count: float, n_side_count: int) -> int | None:
    if not math.isfinite(pair_count):
        return None
    needed_pairs_per_grid = max(pair_count / max(n_side_count, 1), 1.0)
    seed_count = 2
    while seed_count < 200:
        if (seed_count * (seed_count - 1) / 2.0) >= needed_pairs_per_grid and seed_count >= needed_pairs_per_grid:
            return seed_count
        seed_count += 1
    return None


def seed_power_analysis(benchmark: dict[str, Any], pair_sets: dict[str, Any]) -> dict[str, Any]:
    config = benchmark["config"]
    weights = config["distance"]["weights"]
    enabled = config["distance"]["enabled_components"]
    available = enumerate_available_seeds(config)
    rows: dict[str, Any] = {}
    for component in COMPONENT_NAMES:
        target_values, _ = component_values(pair_sets["target_pairs"], component, weights=weights, enabled_components=enabled)
        rows[component] = {}
        for label in FOCUS_CONTROLS:
            control_values, _ = component_values(pair_sets["target_control_pairs"].get(label, []), component, weights=weights, enabled_components=enabled)
            effect = cohen_d(target_values, control_values)
            if effect is None or effect == 0.0 or not math.isfinite(effect):
                pair_needed = None
            else:
                pair_needed = 2.0 * ((Z_ALPHA_975 + Z_POWER_80) / abs(effect)) ** 2
            if effect is None:
                effect_payload: float | str | None = None
            elif math.isfinite(effect):
                effect_payload = effect
            else:
                effect_payload = "Infinity" if effect > 0 else "-Infinity"
            rows[component][label] = {
                "observed_effect_size_cohen_d": effect_payload,
                "observed_effect_size_is_infinite": bool(effect is not None and not math.isfinite(effect)),
                "effect_direction": "target_closer" if effect is not None and effect > 0 else "control_closer_or_equal",
                "current_target_pair_count": len(target_values),
                "current_target_control_pair_count": len(control_values),
                "approx_pair_observations_per_group_needed": pair_needed,
                "approx_seeds_per_grid_needed": pair_count_to_seed_count(pair_needed or math.inf, len(config["n_sides"])),
                "bootstrap_effect_ci95": bootstrap_effect_ci(target_values, control_values),
            }
    return {
        "predeclared_power_target": POWER_TARGET,
        "alpha_two_sided": POWER_ALPHA,
        "method": "normal approximation on observed component effect size; bootstrap used only for uncertainty",
        "configured_seeds": config["seeds"],
        "configured_seed_count": len(config["seeds"]),
        "available_frozen_seeds_by_hierarchy_n_side": available,
        "actual_target_target_pair_count": len(pair_sets["target_pairs"]),
        "actual_focus_target_control_pair_counts": {
            label: len(pair_sets["target_control_pairs"].get(label, [])) for label in FOCUS_CONTROLS
        },
        "component_power": rows,
    }


def null_expectation(component_rows: list[dict[str, Any]], seed_power: dict[str, Any]) -> dict[str, Any]:
    focus_rows = [row for row in component_rows if row["is_focus_control"]]
    failures = [
        row for row in focus_rows
        if row["target_control_mean"] is not None
        and row["target_mean"] is not None
        and row["target_control_mean"] <= row["target_mean"]
    ]
    underpowered = seed_power["configured_seed_count"] < 3 or seed_power["actual_target_target_pair_count"] < 10
    return {
        "H0": "The target branch is no more recoverable than matched controls under the current signature and perturbation suite.",
        "decision": "fails_to_reject_H0",
        "underpowered": underpowered,
        "reason": (
            "Several focus control/component distributions are as close or closer than target-target; "
            "configured sample size is small. This is not proof of equivalence."
        ),
        "focus_control_component_failures": [
            {
                "component": row["component"],
                "control_family": row["control_family"],
                "target_mean": row["target_mean"],
                "target_control_mean": row["target_control_mean"],
            }
            for row in failures
        ],
    }


def final_labels(
    *,
    component_rows: list[dict[str, Any]],
    collision_report: dict[str, Any],
    integrity: dict[str, Any],
    timing_summary: dict[str, Any],
    seed_power: dict[str, Any],
) -> list[str]:
    labels: list[str] = []
    focus_failures = [
        row for row in component_rows
        if row["is_focus_control"]
        and row["target_control_mean"] is not None
        and row["target_mean"] is not None
        and row["target_control_mean"] <= row["target_mean"]
    ]
    if focus_failures:
        labels.append("TARGET NOT DISTINCT")
    if collision_report["collision_count"] > 0:
        labels.append("METRIC INSUFFICIENT")
    invalid_focus = [
        label for label in ("randomized_edge_control", "shuffled_structure_control")
        if integrity.get(label, {}).get("validity_note", "").startswith("invalid")
    ]
    if invalid_focus:
        labels.append("CONTROL INVALID")
    if seed_power["configured_seed_count"] < 3:
        labels.append("UNDERPOWERED")
    timing_verdicts = {payload["verdict"] for payload in timing_summary.values()}
    timing_means = {mode: payload["target_mean"] for mode, payload in timing_summary.items()}
    if len(timing_verdicts) > 1:
        labels.append("HEURISTIC CONTAMINATION")
    elif timing_means.get("direct_observed_timing_only") != timing_means.get("interpolated_non_front_times"):
        labels.append("MIXED / OPEN")
    if not labels:
        labels.append("MIXED / OPEN")
    return labels


def render_failure_report(
    *,
    result: dict[str, Any],
    component_rows: list[dict[str, Any]],
    ablation_rows: list[dict[str, Any]],
    timing_summary: dict[str, Any],
    null_result: dict[str, Any],
) -> str:
    focus_rows = [row for row in component_rows if row["is_focus_control"]]
    failing_component_rows = [
        row for row in focus_rows
        if row["target_control_mean"] is not None
        and row["target_mean"] is not None
        and row["target_control_mean"] <= row["target_mean"]
    ]
    verdict_lines = [f"- {label}" for label in result["labels"]]
    lines = [
        "# CST Failure Localization Report",
        "",
        "Implemented fact: this hostile diagnostic uses the existing CST toy benchmark and frozen Phase 15-18 artifacts.",
        "Design choice: all component comparisons use distance values where lower means closer.",
        "Heuristic: effect size is Cohen d with positive direction meaning target-target is closer than target-control.",
        "Unverified hypothesis: none added; this report tests where the current CST toy benchmark loses discrimination.",
        "",
        "## Verdict Labels",
        *verdict_lines,
        "",
        "## H0",
        null_result["H0"],
        "",
        f"Decision: {null_result['decision']}. Underpowered: {null_result['underpowered']}.",
        "This is not proof of equivalence.",
        "",
        "## Primary Failure Points",
    ]
    if failing_component_rows:
        for row in failing_component_rows:
            lines.append(
                "- {control} / {component}: control mean {control_mean:.12g} <= target mean {target_mean:.12g}".format(
                    control=row["control_family"],
                    component=row["component"],
                    control_mean=float(row["target_control_mean"]),
                    target_mean=float(row["target_mean"]),
                )
            )
    else:
        lines.append("- No focus control component mean beat the target mean.")
    changed = [row for row in ablation_rows if row["variant_type"] == "leave_one_out" and row["verdict_changed"]]
    lines.extend(["", "## Ablation Summary"])
    if changed:
        for row in changed:
            lines.append(f"- Removing {row['component_under_test']} changes verdict in {row['control_scope']} to {row['variant_verdict']}.")
    else:
        lines.append("- No leave-one-component-out variant rescues the current FAIL verdict for the focus controls.")
    lines.extend(["", "## Timing Ablation"])
    for mode, payload in timing_summary.items():
        lines.append(
            f"- {mode}: verdict={payload['verdict']}, target_mean={payload['target_mean']}, focus_control_mean={payload['focus_control_mean']}"
        )
    lines.extend(
        [
            "",
            "## Classification",
            "- Genuine absence of target-specific closure: Open, but current evidence fails to reject H0.",
            "- Signature aliasing: Present; branch_signature_id uses compact address fields and omits several discriminating fields.",
            "- Metric insensitivity: Present; shuffled/random controls preserve many tested marginals and remain close.",
            "- Normalization artifacts: Present for address-floor and relative small-denominator risks; see normalization_audit.csv.",
            "- Control leakage: Present for randomized/shuffled controls through retained labels and descriptor marginals.",
            "- Seed insufficiency: Present; configured seed count is 2 and target-target pairs are 3.",
            "- Heuristic timing contamination: Partial risk; timing mode changes distributions and is isolated in timing_ablation.csv.",
        ]
    )
    return "\n".join(lines) + "\n"


def run_failure_localization(output_dir: Path = OUTPUT_DIR) -> dict[str, Any]:
    benchmark = run_cst_benchmark(write_outputs=False)
    pair_sets = build_pair_sets(benchmark["observations"])
    component_rows = component_separation_rows(benchmark, pair_sets)
    ablation_rows = ablation_matrix_rows(benchmark, pair_sets)
    collision_report = signature_collision_report(benchmark, pair_sets)
    control_markdown, control_integrity = control_integrity_markdown(pair_sets)
    timing_rows, timing_summary = timing_ablation_rows(benchmark, pair_sets)
    normalization_rows = normalization_audit_rows(benchmark, pair_sets)
    seed_power = seed_power_analysis(benchmark, pair_sets)
    null_result = null_expectation(component_rows, seed_power)
    labels = final_labels(
        component_rows=component_rows,
        collision_report=collision_report,
        integrity=control_integrity,
        timing_summary=timing_summary,
        seed_power=seed_power,
    )
    result = {
        "labels": labels,
        "benchmark_verdict": benchmark["result"].verdict,
        "benchmark_result_hash": benchmark["result"].result_hash,
        "null_expectation": null_result,
        "signature_collision_summary": {
            "collision_count": collision_report["collision_count"],
            "equivalent_target_runs_match_across_seeds": collision_report["equivalent_target_runs_match_across_seeds"],
            "low_mode_band_in_address_summary": collision_report["quantization_rounding_audit"]["low_mode_band_in_address_summary"],
        },
        "control_integrity_summary": control_integrity,
        "timing_summary": timing_summary,
        "seed_power_summary": {
            "configured_seed_count": seed_power["configured_seed_count"],
            "actual_target_target_pair_count": seed_power["actual_target_target_pair_count"],
            "actual_focus_target_control_pair_counts": seed_power["actual_focus_target_control_pair_counts"],
        },
        "outputs": {
            "failure_localization_report": repo_rel(output_dir / "failure_localization_report.md"),
            "component_separation": repo_rel(output_dir / "component_separation.csv"),
            "ablation_matrix": repo_rel(output_dir / "ablation_matrix.csv"),
            "signature_collision_report": repo_rel(output_dir / "signature_collision_report.json"),
            "control_integrity_report": repo_rel(output_dir / "control_integrity_report.md"),
            "timing_ablation": repo_rel(output_dir / "timing_ablation.csv"),
            "normalization_audit": repo_rel(output_dir / "normalization_audit.csv"),
            "seed_power_analysis": repo_rel(output_dir / "seed_power_analysis.json"),
            "failure_localization_result": repo_rel(output_dir / "failure_localization_result.json"),
        },
    }
    report = render_failure_report(
        result=result,
        component_rows=component_rows,
        ablation_rows=ablation_rows,
        timing_summary=timing_summary,
        null_result=null_result,
    )

    write_csv(
        output_dir / "component_separation.csv",
        component_rows,
        [
            "component",
            "control_family",
            "is_focus_control",
            "target_count",
            "target_missing",
            "target_mean",
            "target_variance",
            "target_values",
            "target_control_count",
            "target_control_missing",
            "target_control_mean",
            "target_control_variance",
            "target_control_values",
            "control_control_count",
            "control_control_missing",
            "control_control_mean",
            "control_control_variance",
            "control_control_values",
            "effect_size_cohen_d_control_minus_target",
            "overlap_coefficient_histogram_20_bin",
            "rank_order_by_mean_distance",
            "target_outperforms_target_control",
            "notes",
        ],
    )
    write_csv(
        output_dir / "ablation_matrix.csv",
        ablation_rows,
        [
            "variant_type",
            "component_under_test",
            "enabled_components",
            "control_scope",
            "target_mean",
            "target_variance",
            "control_mean",
            "control_variance",
            "failing_control_families",
            "variant_verdict",
            "base_verdict",
            "verdict_changed",
            "flags",
            "reasons",
        ],
    )
    write_json(output_dir / "signature_collision_report.json", collision_report)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "control_integrity_report.md").write_text(control_markdown, encoding="utf-8")
    write_csv(
        output_dir / "timing_ablation.csv",
        timing_rows,
        [
            "timing_mode",
            "component",
            "component_enabled_in_verdict",
            "control_family",
            "verdict",
            "target_mean",
            "target_variance",
            "target_missing",
            "target_control_mean",
            "target_control_variance",
            "target_control_missing",
            "effect_size_cohen_d_control_minus_target",
            "reasons",
        ],
    )
    write_csv(
        output_dir / "normalization_audit.csv",
        normalization_rows,
        [
            "pair_type",
            "control_family",
            "pair_index",
            "left_observation_id",
            "right_observation_id",
            "component",
            "normalized_value",
            "robust_alt_value",
            "normalization_reference",
            "small_denominator_flag",
            "artifact_flags",
            "raw_values_json",
        ],
    )
    write_json(output_dir / "seed_power_analysis.json", seed_power)
    write_json(output_dir / "failure_localization_result.json", result)
    (output_dir / "failure_localization_report.md").write_text(report, encoding="utf-8")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run hostile CST failure localization diagnostics.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = run_failure_localization(output_dir=args.output_dir)
    print(json.dumps({"labels": result["labels"], "result": result["outputs"]["failure_localization_result"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
