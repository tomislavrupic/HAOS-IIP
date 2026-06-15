from __future__ import annotations

import argparse
import copy
import itertools
import json
import math
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from experiments.cst.cst_controls import edges_from_chain
from experiments.cst.cst_io import read_json, repo_rel, stable_hash, write_csv, write_json
from experiments.cst.cst_metrics import (
    COMPONENT_NAMES,
    branch_signature_id,
    compute_closure_distance,
    mean,
    population_std,
)
from experiments.cst.cst_types import ClosureSignature, PerturbationSpec, RecoveryObservation
from experiments.cst.diagnostics.failure_localization import (
    build_pair_sets,
    cohen_d,
    component_values,
    enumerate_available_seeds,
    variant_distances,
)
from experiments.cst.run_cst_benchmark import DEFAULT_OUTPUT_DIR, run_cst_benchmark


HERE = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT_PATH = HERE / "configs" / "phase_2_2_precommitment_contract.json"
DEFAULT_REPAIR_OUTPUT_DIR = DEFAULT_OUTPUT_DIR / "phase_2_2"
FOCUS_CONTROLS = [
    "periodic_diagonal_augmented_control",
    "randomized_edge_control",
    "shuffled_structure_control",
]
PHASE_2_2_CONTROL_LABELS = [
    "phase_2_2_invariant_control",
    "phase_2_2_structural_control",
    "phase_2_2_spectral_control",
    "phase_2_2_causal_control",
    "phase_2_2_temporal_control",
    "phase_2_2_address_control",
]


def chain_labels(fields: dict[str, Any]) -> list[str]:
    chain = fields.get("event_chain_summary", {}).get("chain_signature", "")
    return [item for item in str(chain).split(">") if item]


def update_chain(fields: dict[str, Any], labels: list[str]) -> None:
    fields.setdefault("event_chain_summary", {})["chain_signature"] = ">".join(labels)
    fields["event_chain_summary"]["dominant_chain_prefix"] = labels[:2]
    fields.setdefault("influence_edge_summary", {})["edges"] = edges_from_chain(labels)
    fields.setdefault("address_summary", {})["dominant_chain_prefix"] = labels[:2]


def rotated(values: list[Any], amount: int = 1) -> list[Any]:
    if not values:
        return []
    amount = amount % len(values)
    return values[amount:] + values[:amount]


def transform_invariant(fields: dict[str, Any]) -> dict[str, Any]:
    fields = copy.deepcopy(fields)
    summary = fields.setdefault("invariant_summary", {})
    keys = [
        "terminal_low_k_fraction",
        "terminal_low_k_power",
        "transport_character_index",
        "shell_latency_span",
        "shell_order_score",
    ]
    values = [summary.get(key) for key in keys]
    for key, value in zip(keys, rotated(values, 2)):
        summary[key] = value
    fields.setdefault("phase_2_2_transform", {})["contract"] = "invariant_relationship_rotation"
    return fields


def transform_structural(fields: dict[str, Any]) -> dict[str, Any]:
    fields = copy.deepcopy(fields)
    labels = chain_labels(fields)
    if len(labels) > 2:
        labels = labels[::2] + labels[1::2]
    update_chain(fields, labels)
    fields.setdefault("phase_2_2_transform", {})["contract"] = "degree_count_preserving_topology_scramble"
    return fields


def transform_spectral(fields: dict[str, Any]) -> dict[str, Any]:
    fields = copy.deepcopy(fields)
    spectral = fields.setdefault("spectral_low_mode_summary", {})
    for key in ("terminal_low_k_fraction", "terminal_low_k_power"):
        value = spectral.get(key)
        if isinstance(value, (int, float)):
            spectral[key] = max(0.0, min(1.0, 1.0 - float(value)))
    slope = spectral.get("terminal_fluctuation_slope")
    if isinstance(slope, (int, float)):
        spectral["terminal_fluctuation_slope"] = -float(slope)
    spectral["low_mode_band"] = "phase_2_2_spectral_destroyed"
    fields.setdefault("phase_2_2_transform", {})["contract"] = "low_mode_alignment_complement"
    return fields


def transform_causal(fields: dict[str, Any]) -> dict[str, Any]:
    fields = copy.deepcopy(fields)
    causal = fields.setdefault("causal_summary", {})
    depth = causal.setdefault("front_depth_profile", {})
    if isinstance(depth.get("front_near_depth"), (int, float)):
        depth["front_near_depth"] = float(depth["front_near_depth"]) + 1.0
    if isinstance(depth.get("front_far_depth"), (int, float)):
        depth["front_far_depth"] = max(0.0, float(depth["front_far_depth"]) - 1.0)
    if isinstance(causal.get("acyclicity_score"), (int, float)):
        causal["acyclicity_score"] = max(0.0, min(1.0, 1.0 - float(causal["acyclicity_score"])))
    if isinstance(causal.get("order_compatibility_score"), (int, float)):
        causal["order_compatibility_score"] = max(0.0, 1.0 - float(causal["order_compatibility_score"]))
    if isinstance(causal.get("mismatch_rate"), (int, float)):
        causal["mismatch_rate"] = min(1.0, float(causal["mismatch_rate"]) + 0.25)
    fields.setdefault("address_summary", {})["front_depth_profile"] = copy.deepcopy(depth)
    fields.setdefault("phase_2_2_transform", {})["contract"] = "causal_dependency_perturbation"
    return fields


def transform_temporal(fields: dict[str, Any]) -> dict[str, Any]:
    fields = copy.deepcopy(fields)
    labels = list(reversed(chain_labels(fields)))
    old_times = fields.get("event_chain_summary", {}).get("event_times", {})
    values = list(old_times.values())
    new_times = {label: value for label, value in zip(labels, values)}
    update_chain(fields, labels)
    fields.setdefault("event_chain_summary", {})["event_times"] = new_times
    fields.setdefault("phase_2_2_transform", {})["contract"] = "sequence_reversal_preserving_time_marginals"
    return fields


def transform_address(fields: dict[str, Any]) -> dict[str, Any]:
    fields = copy.deepcopy(fields)
    labels = chain_labels(fields)
    address = fields.setdefault("address_summary", {})
    address["dominant_chain_prefix"] = list(reversed(labels[-2:])) if len(labels) >= 2 else labels
    depth = copy.deepcopy(address.get("front_depth_profile", {}))
    if isinstance(depth.get("front_near_depth"), (int, float)):
        depth["front_near_depth"] = float(depth["front_near_depth"]) + 2.0
    if isinstance(depth.get("front_far_depth"), (int, float)):
        depth["front_far_depth"] = float(depth["front_far_depth"]) + 2.0
    address["front_depth_profile"] = depth
    fields.setdefault("phase_2_2_transform", {})["contract"] = "stable_address_correspondence_break"
    return fields


CONTROL_TRANSFORMS = {
    "phase_2_2_invariant_control": transform_invariant,
    "phase_2_2_structural_control": transform_structural,
    "phase_2_2_spectral_control": transform_spectral,
    "phase_2_2_causal_control": transform_causal,
    "phase_2_2_temporal_control": transform_temporal,
    "phase_2_2_address_control": transform_address,
}


def clone_phase_2_2_control(
    source: RecoveryObservation,
    *,
    label: str,
    contract_hash: str,
    signature_version: str,
) -> RecoveryObservation:
    fields = CONTROL_TRANSFORMS[label](source.closure_signature.fields)
    signature = ClosureSignature(
        branch_signature_id=branch_signature_id(fields, float_precision=6),
        signature_version=signature_version,
        operational_equivalence_note="Phase 2.2 repaired identity uses distance equivalence; this hash is grouping only.",
        fields=fields,
        unavailable_components={},
        source_artifacts=source.closure_signature.source_artifacts,
        warnings=list(source.closure_signature.warnings) + [f"phase_2_2_control={label}"],
    )
    provenance_payload = {
        "source_run_instance_id": source.provenance.run_instance_id,
        "source_observation_id": source.observation_id,
        "phase_2_2_control": label,
        "contract_hash": contract_hash,
    }
    provenance = replace(
        source.provenance,
        run_instance_id=stable_hash(provenance_payload, prefix="run_"),
        control_group=label,
    )
    observation_id = stable_hash(
        {
            "run_instance_id": provenance.run_instance_id,
            "coarse_signature_id": signature.branch_signature_id,
            "phase_2_2_control": label,
        },
        prefix="obs_",
    )
    return replace(
        source,
        observation_id=observation_id,
        observation_kind="phase_2_2_control",
        control_role="negative",
        control_group_label=label,
        perturbation=PerturbationSpec(
            name=label,
            family="phase_2_2_component_specific_control",
            parameters={"source_observation_id": source.observation_id},
            source="phase_2_2_precommitment_contract",
        ),
        provenance=provenance,
        closure_signature=signature,
        notes=[f"Phase 2.2 control derived from {source.observation_id} under precommitted contract"],
    )


def build_phase_2_2_observations(
    benchmark: dict[str, Any],
    *,
    contract_hash: str,
    signature_version: str,
) -> list[RecoveryObservation]:
    observations = list(benchmark["observations"])
    targets = [item for item in observations if item.observation_kind == "target"]
    for target in targets:
        for label in PHASE_2_2_CONTROL_LABELS:
            observations.append(
                clone_phase_2_2_control(
                    target,
                    label=label,
                    contract_hash=contract_hash,
                    signature_version=signature_version,
                )
            )
    return observations


def degree_sequence(edges: list[str]) -> list[int]:
    degree: dict[str, int] = {}
    for edge in edges:
        if "->" not in edge:
            continue
        left, right = edge.split("->", 1)
        degree[left] = degree.get(left, 0) + 1
        degree[right] = degree.get(right, 0) + 1
    return sorted(degree.values())


def direct_and_interpolated_time_keys(observation: RecoveryObservation) -> tuple[list[str], list[str]]:
    fields = observation.closure_signature.fields
    time_keys = set(fields.get("event_chain_summary", {}).get("event_times", {}))
    direct_keys = set(fields.get("distance_surrogate_summary", {}))
    direct_keys.update({"front_near", "front_far"})
    direct = sorted(time_keys & direct_keys)
    return direct, sorted(time_keys - set(direct))


def field_record(
    *,
    section: str,
    values: Any,
    schema: dict[str, Any],
) -> dict[str, Any]:
    meta = schema[section]
    return {
        "values": values,
        "semantic_justification": meta["semantic_justification"],
        "source_artifact": meta["source_artifact"],
        "extraction_method": meta["extraction_method"],
        "normalization": meta["normalization"],
        "availability": meta["availability"],
        "contributes_to_equivalence": bool(meta["contributes_to_equivalence"]),
        "contributes_to_reporting_only": not bool(meta["contributes_to_equivalence"]),
    }


def fine_signature_for(observation: RecoveryObservation, contract: dict[str, Any]) -> dict[str, Any]:
    fields = observation.closure_signature.fields
    schema = contract["fine_signature_schema"]
    chain = chain_labels(fields)
    edges = [str(edge) for edge in fields.get("influence_edge_summary", {}).get("edges", [])]
    direct_keys, interpolated_keys = direct_and_interpolated_time_keys(observation)
    equivalence_payload = {
        "invariant_summary": fields.get("invariant_summary", {}),
        "structural_graph_summary": {
            "transition_edges": edges,
            "degree_sequence": degree_sequence(edges),
            "edge_count": fields.get("influence_edge_summary", {}).get("edge_count"),
            "stable_edge_count": fields.get("influence_edge_summary", {}).get("stable_edge_count"),
            "mean_edge_reproducibility": fields.get("influence_edge_summary", {}).get("mean_edge_reproducibility"),
            "chain_length": len(chain),
        },
        "spectral_summary": fields.get("spectral_low_mode_summary", {}),
        "causal_summary": fields.get("causal_summary", {}),
        "temporal_sequence_summary": {
            "chain_signature": fields.get("event_chain_summary", {}).get("chain_signature"),
            "transition_pairs": edges_from_chain(chain),
            "event_times": fields.get("event_chain_summary", {}).get("event_times", {}),
            "direct_time_keys": direct_keys,
            "interpolated_time_keys": interpolated_keys,
            "timing_source_quality": "partial_interpolated" if interpolated_keys else "direct_observed",
        },
        "stable_address_summary": fields.get("address_summary", {}),
    }
    reporting_payload = {
        "shell_distance_surrogate_summary": {
            **fields.get("shell_ordering_summary", {}),
            "per_shell_distance_rows": fields.get("distance_surrogate_summary", {}),
        },
        "provenance_metadata": {
            "observation_id": observation.observation_id,
            "run_instance_id": observation.provenance.run_instance_id,
            "n_side": observation.n_side,
            "seed": observation.seed,
            "hierarchy_label": observation.hierarchy_label,
            "control_group_label": observation.control_group_label,
            "source_artifact_hashes": observation.provenance.source_artifact_hashes,
        },
    }
    fine_signature = {
        "observation_id": observation.observation_id,
        "observation_kind": observation.observation_kind,
        "control_group_label": observation.control_group_label,
        "run_instance_id": observation.provenance.run_instance_id,
        "coarse_signature_id": observation.closure_signature.branch_signature_id,
        "coarse_signature_policy": "grouping_only_not_equivalence",
        "fine_signature_record_id": stable_hash(equivalence_payload, prefix="fine_", float_precision=12),
        "fine_signature_hash_policy": "names_representation_only_not_equivalence",
        "equivalence_payload": {
            name: field_record(section=name, values=values, schema=schema)
            for name, values in equivalence_payload.items()
        },
        "reporting_payload": {
            name: field_record(section=name, values=values, schema=schema)
            for name, values in reporting_payload.items()
        },
    }
    return fine_signature


def evaluate_component_value_map(
    component_values_map: dict[str, float | None],
    availability_map: dict[str, str],
    contract: dict[str, Any],
    *,
    mode: str,
) -> dict[str, Any]:
    mode_contract = contract["equivalence"][mode]
    tolerances = mode_contract["tolerances"]
    rows = []
    available_count = 0
    equivalent = True
    failure_reasons = []
    boundary_epsilon = float(contract["equivalence"]["threshold_boundary_sensitivity"]["epsilon"])
    boundary_hits = []
    for component in COMPONENT_NAMES:
        value = component_values_map.get(component)
        availability = availability_map.get(component, "unavailable")
        tolerance = float(tolerances[component])
        if availability != "available" or value is None:
            passed = False if mode_contract["unavailable_component_policy"] == "fail" else None
            if mode_contract["unavailable_component_policy"] == "fail":
                equivalent = False
                failure_reasons.append(f"{component} unavailable")
        else:
            available_count += 1
            passed = float(value) <= tolerance
            if not passed:
                equivalent = False
                failure_reasons.append(f"{component}={value:.12g} > {tolerance:.12g}")
            if abs(float(value) - tolerance) <= boundary_epsilon:
                boundary_hits.append(component)
        rows.append(
            {
                "component": component,
                "value": value,
                "availability": availability,
                "tolerance": tolerance,
                "passed": passed,
                "within_boundary_epsilon": component in boundary_hits,
            }
        )
    if available_count < int(mode_contract["minimum_required_available_components"]):
        equivalent = False
        failure_reasons.append(
            f"available components {available_count} < required {mode_contract['minimum_required_available_components']}"
        )
    return {
        "mode": mode,
        "equivalent": bool(equivalent),
        "available_component_count": available_count,
        "component_results": rows,
        "failure_reasons": failure_reasons,
        "boundary_hits": boundary_hits,
        "boundary_sensitive": bool(boundary_hits),
    }


def evaluate_observation_equivalence(
    left: RecoveryObservation,
    right: RecoveryObservation,
    contract: dict[str, Any],
    *,
    mode: str,
) -> dict[str, Any]:
    weights = {component: 1.0 for component in COMPONENT_NAMES}
    distance = compute_closure_distance(
        left,
        right,
        weights=weights,
        enabled_components=list(COMPONENT_NAMES),
        scalar_enabled=False,
    )
    value_map = {
        component: distance.components[component].value
        for component in COMPONENT_NAMES
    }
    availability_map = {
        component: distance.components[component].availability
        for component in COMPONENT_NAMES
    }
    result = evaluate_component_value_map(value_map, availability_map, contract, mode=mode)
    result.update(
        {
            "left_observation_id": left.observation_id,
            "right_observation_id": right.observation_id,
            "left_group": left.control_group_label,
            "right_group": right.control_group_label,
            "left_kind": left.observation_kind,
            "right_kind": right.observation_kind,
            "left_coarse_signature_id": left.closure_signature.branch_signature_id,
            "right_coarse_signature_id": right.closure_signature.branch_signature_id,
            "hash_equality_used_for_equivalence": False,
        }
    )
    return result


def collect_known_positive_pairs(observations: list[RecoveryObservation]) -> list[tuple[RecoveryObservation, RecoveryObservation]]:
    targets_by_id = {
        observation.observation_id: observation
        for observation in observations
        if observation.observation_kind == "target"
    }
    pairs: list[tuple[RecoveryObservation, RecoveryObservation]] = []
    for observation in observations:
        if observation.control_group_label != "seed_repeat_control":
            continue
        source_id = observation.perturbation.parameters.get("source_observation_id")
        source = targets_by_id.get(source_id)
        if source is not None:
            pairs.append((source, observation))
    return pairs


def equivalence_pairs(pair_sets: dict[str, Any]) -> list[tuple[str, RecoveryObservation, RecoveryObservation]]:
    pairs: list[tuple[str, RecoveryObservation, RecoveryObservation]] = []
    for left, right in pair_sets["target_pairs"]:
        pairs.append(("target_target", left, right))
    for label in FOCUS_CONTROLS + PHASE_2_2_CONTROL_LABELS:
        for left, right in pair_sets["target_control_pairs"].get(label, []):
            pairs.append((f"target_control:{label}", left, right))
    for left, right in pair_sets.get("known_positive_pairs", []):
        pairs.append(("known_positive:seed_repeat_control", left, right))
    return pairs


def equivalence_rows(pair_sets: dict[str, Any], contract: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for pair_type, left, right in equivalence_pairs(pair_sets):
        for mode in ("strict", "exploratory"):
            result = evaluate_observation_equivalence(left, right, contract, mode=mode)
            rows.append(
                {
                    "pair_type": pair_type,
                    "mode": mode,
                    "left_observation_id": left.observation_id,
                    "right_observation_id": right.observation_id,
                    "left_group": left.control_group_label,
                    "right_group": right.control_group_label,
                    "equivalent": result["equivalent"],
                    "available_component_count": result["available_component_count"],
                    "failure_reasons": "|".join(result["failure_reasons"]),
                    "boundary_hits": "|".join(result["boundary_hits"]),
                    "boundary_sensitive": result["boundary_sensitive"],
                    "hash_equality_used_for_equivalence": False,
                    "component_results_json": json.dumps(result["component_results"], sort_keys=True),
                }
            )
    return rows


def strict_results_for(rows: list[dict[str, Any]], prefix: str) -> list[dict[str, Any]]:
    return [row for row in rows if row["mode"] == "strict" and str(row["pair_type"]).startswith(prefix)]


def component_control_validation_rows(
    pair_sets: dict[str, Any],
    contract: dict[str, Any],
) -> list[dict[str, Any]]:
    validation = contract["component_control_validation"]
    min_shift = float(validation["intended_metric_min_shift_over_target_mean"])
    max_unaffected_shift = float(validation["unaffected_metric_max_shift_over_target_mean"])
    rows: list[dict[str, Any]] = []
    target_means: dict[str, float] = {}
    for component in COMPONENT_NAMES:
        values, _ = component_values(
            pair_sets["target_pairs"],
            component,
            weights={name: 1.0 for name in COMPONENT_NAMES},
            enabled_components=list(COMPONENT_NAMES),
        )
        target_means[component] = mean(values) if values else 0.0
    for label in PHASE_2_2_CONTROL_LABELS:
        contract_row = contract["component_control_contracts"][label]
        affected = set(contract_row["expected_affected_metrics"])
        unaffected = set(contract_row["expected_unaffected_metrics"])
        for component in COMPONENT_NAMES:
            values, missing = component_values(
                pair_sets["target_control_pairs"].get(label, []),
                component,
                weights={name: 1.0 for name in COMPONENT_NAMES},
                enabled_components=list(COMPONENT_NAMES),
            )
            control_mean = mean(values) if values else None
            target_mean = target_means[component]
            shift = None if control_mean is None else float(control_mean) - float(target_mean)
            source_response = control_mean
            expectation = "not_specified"
            passed = True
            invalidation_reason = ""
            if component in affected:
                expectation = "affected_should_increase"
                passed = bool(source_response is not None and source_response >= min_shift)
                if not passed:
                    invalidation_reason = "intended metric did not increase by predeclared minimum"
            elif component in unaffected:
                expectation = "unaffected_should_remain_stable"
                passed = bool(source_response is not None and source_response <= max_unaffected_shift)
                if not passed:
                    invalidation_reason = "unaffected metric shifted beyond predeclared tolerance"
            rows.append(
                {
                    "control_label": label,
                    "control_class": contract_row["control_class"],
                    "component": component,
                    "expectation": expectation,
                    "target_target_mean": target_mean,
                    "target_control_mean": control_mean,
                    "shift_control_minus_target": shift,
                    "control_distance_from_source": source_response,
                    "missing_count": missing,
                    "passed": passed,
                    "invalidation_condition": contract_row["invalidation_condition"],
                    "invalidation_reason": invalidation_reason,
                    "expected_affected_metrics": "|".join(contract_row["expected_affected_metrics"]),
                    "expected_unaffected_metrics": "|".join(contract_row["expected_unaffected_metrics"]),
                }
            )
    return rows


def signature_collision_v2(fine_signatures: list[dict[str, Any]]) -> dict[str, Any]:
    coarse_groups: dict[str, list[dict[str, Any]]] = {}
    fine_groups: dict[str, list[dict[str, Any]]] = {}
    for signature in fine_signatures:
        coarse_groups.setdefault(signature["coarse_signature_id"], []).append(signature)
        fine_groups.setdefault(signature["fine_signature_record_id"], []).append(signature)
    coarse_collisions = [
        {
            "coarse_signature_id": key,
            "member_count": len(members),
            "control_groups": sorted({item["control_group_label"] for item in members}),
            "policy": "coarse hash is grouping only; not an equivalence decision",
        }
        for key, members in sorted(coarse_groups.items())
        if len(members) > 1
    ]
    allowed_duplicate_groups = []
    unexplained = []
    for key, members in sorted(fine_groups.items()):
        if len(members) <= 1:
            continue
        groups = sorted({item["control_group_label"] for item in members})
        payload = {
            "fine_signature_record_id": key,
            "member_count": len(members),
            "control_groups": groups,
            "observation_ids": [item["observation_id"] for item in members],
        }
        if set(groups).issubset({"target", "seed_repeat_control"}):
            allowed_duplicate_groups.append({**payload, "reason": "known positive replay may duplicate the fine structure"})
        else:
            unexplained.append(payload)
    return {
        "coarse_collision_count": len(coarse_collisions),
        "coarse_collisions": coarse_collisions,
        "fine_duplicate_count": len(allowed_duplicate_groups) + len(unexplained),
        "allowed_fine_duplicate_groups": allowed_duplicate_groups,
        "unexplained_fine_collision_count": len(unexplained),
        "unexplained_fine_collisions": unexplained,
        "hash_policy": "fine_signature_record_id names the inspectable record; equivalence remains distance-based",
    }


def power_gate(contract: dict[str, Any], benchmark: dict[str, Any], pair_sets: dict[str, Any]) -> dict[str, Any]:
    config = benchmark["config"]
    available = enumerate_available_seeds(config)
    usable_seed_sets: list[set[int]] = []
    for hierarchy_data in available.values():
        for n_side_data in hierarchy_data.values():
            usable_seed_sets.append(set(n_side_data["intersection_usable_for_cst"]))
    intersection = sorted(set.intersection(*usable_seed_sets)) if usable_seed_sets else []
    gate = contract["power_gate"]
    target_pair_count = len(pair_sets["target_pairs"])
    seed_count_pass = len(intersection) >= int(gate["minimum_independent_seed_count"])
    pair_count_pass = target_pair_count >= int(gate["minimum_target_target_pair_count"])
    return {
        "available_frozen_seeds_by_hierarchy_n_side": available,
        "valid_independent_seed_intersection": intersection,
        "actual_independent_seed_count": len(intersection),
        "minimum_independent_seed_count": gate["minimum_independent_seed_count"],
        "actual_target_target_pair_count": target_pair_count,
        "minimum_target_target_pair_count": gate["minimum_target_target_pair_count"],
        "alpha_two_sided": gate["alpha_two_sided"],
        "power_target": gate["power_target"],
        "bootstrap_policy": gate["bootstrap_policy"],
        "seed_count_pass": seed_count_pass,
        "target_pair_count_pass": pair_count_pass,
        "power_gate_pass": bool(seed_count_pass and pair_count_pass),
        "status": "PASS" if seed_count_pass and pair_count_pass else "UNDERPOWERED",
    }


def phase_2_2_labels(
    *,
    equivalence: list[dict[str, Any]],
    validation_rows: list[dict[str, Any]],
    collisions: dict[str, Any],
    power: dict[str, Any],
) -> list[str]:
    labels: list[str] = []
    strict_target = strict_results_for(equivalence, "target_target")
    strict_seed_repeat = strict_results_for(equivalence, "known_positive:seed_repeat_control")
    strict_focus = []
    for label in FOCUS_CONTROLS:
        strict_focus.extend(strict_results_for(equivalence, f"target_control:{label}"))
    strict_destructive = []
    for label in PHASE_2_2_CONTROL_LABELS:
        strict_destructive.extend(strict_results_for(equivalence, f"target_control:{label}"))

    seed_repeat_ok = bool(strict_seed_repeat) and all(row["equivalent"] for row in strict_seed_repeat)
    target_equiv_ok = bool(strict_target) and all(row["equivalent"] for row in strict_target)
    destructive_non_equiv = bool(strict_destructive) and all(not row["equivalent"] for row in strict_destructive)
    validation_failed = [row for row in validation_rows if row["expectation"] != "not_specified" and not row["passed"]]
    validation_has_success = any(
        row["expectation"] == "affected_should_increase" and row["passed"]
        for row in validation_rows
    )
    if seed_repeat_ok and destructive_non_equiv and validation_has_success:
        labels.append("INSTRUMENT PARTIALLY DISCRIMINATIVE")
    if seed_repeat_ok and destructive_non_equiv and not validation_failed and power["power_gate_pass"] and collisions["unexplained_fine_collision_count"] == 0:
        labels.append("INSTRUMENT DISCRIMINATIVE")
    if validation_failed or collisions["unexplained_fine_collision_count"] > 0:
        labels.append("INSTRUMENT INVALID")
    if validation_failed:
        labels.append("CONTROL INVALID")
    if collisions["unexplained_fine_collision_count"] > 0:
        labels.append("SIGNATURE INSUFFICIENT")
    if not power["power_gate_pass"]:
        labels.append("UNDERPOWERED")

    focus_non_equiv = bool(strict_focus) and all(not row["equivalent"] for row in strict_focus)
    if power["power_gate_pass"] and target_equiv_ok and focus_non_equiv and not validation_failed:
        labels.append("TARGET DISTINCT")
    else:
        labels.append("TARGET NOT DISTINCT")
    if "INSTRUMENT DISCRIMINATIVE" not in labels or "UNDERPOWERED" in labels or "CONTROL INVALID" in labels:
        labels.append("MIXED / OPEN")
    return list(dict.fromkeys(labels))


def component_separation_summary(pair_sets: dict[str, Any]) -> dict[str, Any]:
    weights = {component: 1.0 for component in COMPONENT_NAMES}
    enabled = list(COMPONENT_NAMES)
    summary: dict[str, Any] = {}
    for component in COMPONENT_NAMES:
        target_values, _ = component_values(pair_sets["target_pairs"], component, weights=weights, enabled_components=enabled)
        summary[component] = {}
        for label in FOCUS_CONTROLS:
            control_values, _ = component_values(pair_sets["target_control_pairs"].get(label, []), component, weights=weights, enabled_components=enabled)
            effect = cohen_d(target_values, control_values)
            if effect is None:
                effect_payload: float | str | None = None
            elif math.isfinite(effect):
                effect_payload = effect
            else:
                effect_payload = "Infinity" if effect > 0 else "-Infinity"
            summary[component][label] = {
                "target_mean": mean(target_values) if target_values else None,
                "target_control_mean": mean(control_values) if control_values else None,
                "cohen_d_control_minus_target": effect_payload,
                "cohen_d_is_infinite": bool(effect is not None and not math.isfinite(effect)),
                "target_outperforms": bool(target_values and control_values and mean(target_values) < mean(control_values)),
            }
    return summary


def result_hash_payload(payload: dict[str, Any]) -> str:
    comparable = copy.deepcopy(payload)
    comparable.pop("result_hash", None)
    return stable_hash(comparable, prefix="phase_2_2_")


def render_report(result: dict[str, Any], power: dict[str, Any], validation_rows: list[dict[str, Any]], collisions: dict[str, Any]) -> str:
    invalid_rows = [row for row in validation_rows if row["expectation"] != "not_specified" and not row["passed"]]
    lines = [
        "# CST Phase 2.2 Discriminative Instrument Repair",
        "",
        "Implemented fact: this repair layer consumes the precommitted contract before generating repaired outputs.",
        "Design choice: hashes name representations only; strict component-distance equivalence determines identity.",
        "Heuristic: component-specific controls are deterministic transforms over existing CST signatures, not new frozen physics runs.",
        "Unverified hypothesis: none added; target victory is not a success criterion.",
        "",
        "## Verdict Labels",
    ]
    for label in result["labels"]:
        lines.append(f"- {label}")
    lines.extend(
        [
            "",
            "## Power Gate",
            f"- status: {power['status']}",
            f"- actual independent seeds: {power['actual_independent_seed_count']}",
            f"- required independent seeds: {power['minimum_independent_seed_count']}",
            f"- actual target-target pairs: {power['actual_target_target_pair_count']}",
            f"- required target-target pairs: {power['minimum_target_target_pair_count']}",
            "",
            "## Control Validation",
        ]
    )
    if invalid_rows:
        for row in invalid_rows[:20]:
            lines.append(
                f"- {row['control_label']} / {row['component']}: {row['invalidation_reason']}"
            )
    else:
        lines.append("- all predeclared component-control expectations passed")
    lines.extend(
        [
            "",
            "## Signature Collision V2",
            f"- coarse collision count: {collisions['coarse_collision_count']} (grouping-only)",
            f"- unexplained fine collision count: {collisions['unexplained_fine_collision_count']}",
            "",
            "## Threshold Boundary Sensitivity",
            f"- official strict boundary-sensitive pair count: {result['threshold_boundary_sensitivity']['official_strict_boundary_sensitive_pair_count']}",
            f"- epsilon: {result['threshold_boundary_sensitivity']['epsilon']}",
            "",
            "## Limitations",
            "- This does not alter frozen HAOS-IIP phases or haos_core APIs.",
            "- This does not add fields because they favor the target; fields come from the static contract.",
            "- The official verdict does not use scalar compression or exploratory weight families.",
        ]
    )
    return "\n".join(lines) + "\n"


def run_phase_2_2_repair(
    *,
    contract_path: Path = DEFAULT_CONTRACT_PATH,
    output_dir: Path = DEFAULT_REPAIR_OUTPUT_DIR,
    write_outputs: bool = True,
) -> dict[str, Any]:
    contract = read_json(contract_path)
    contract_hash = stable_hash(contract, prefix="contract_")
    benchmark = run_cst_benchmark(write_outputs=False, scalar_enabled_override=False)
    observations = build_phase_2_2_observations(
        benchmark,
        contract_hash=contract_hash,
        signature_version=contract["contract_version"],
    )
    pair_sets = build_pair_sets(observations)
    pair_sets["known_positive_pairs"] = collect_known_positive_pairs(observations)
    fine_signatures = [fine_signature_for(observation, contract) for observation in observations]
    equivalence = equivalence_rows(pair_sets, contract)
    validation_rows = component_control_validation_rows(pair_sets, contract)
    collisions = signature_collision_v2(fine_signatures)
    power = power_gate(contract, benchmark, pair_sets)
    labels = phase_2_2_labels(
        equivalence=equivalence,
        validation_rows=validation_rows,
        collisions=collisions,
        power=power,
    )
    separation = component_separation_summary(pair_sets)
    result = {
        "contract_hash": contract_hash,
        "benchmark_result_hash": benchmark["result"].result_hash,
        "official_mode": contract["verdict_logic"]["official_mode"],
        "official_scalar_enabled": False,
        "labels": labels,
        "power_gate_status": power["status"],
        "component_control_invalid_count": len([
            row for row in validation_rows if row["expectation"] != "not_specified" and not row["passed"]
        ]),
        "signature_collision_v2": {
            "coarse_collision_count": collisions["coarse_collision_count"],
            "unexplained_fine_collision_count": collisions["unexplained_fine_collision_count"],
        },
        "target_equivalence": {
            "strict_target_target_all_equivalent": all(
                row["equivalent"] for row in strict_results_for(equivalence, "target_target")
            ),
            "strict_focus_controls_all_non_equivalent": all(
                not row["equivalent"]
                for label in FOCUS_CONTROLS
                for row in strict_results_for(equivalence, f"target_control:{label}")
            ),
        },
        "threshold_boundary_sensitivity": {
            "epsilon": contract["equivalence"]["threshold_boundary_sensitivity"]["epsilon"],
            "boundary_sensitive_pair_count": sum(1 for row in equivalence if row["boundary_sensitive"]),
            "official_strict_boundary_sensitive_pair_count": sum(
                1 for row in equivalence if row["mode"] == "strict" and row["boundary_sensitive"]
            ),
        },
        "component_separation_summary": separation,
        "outputs": {
            "precommitment_contract": repo_rel(output_dir / "precommitment_contract.json"),
            "signature_schema": repo_rel(output_dir / "signature_schema.json"),
            "fine_signatures": repo_rel(output_dir / "fine_signatures.json"),
            "equivalence_tolerances": repo_rel(output_dir / "equivalence_tolerances.json"),
            "equivalence_results": repo_rel(output_dir / "equivalence_results.csv"),
            "control_contracts": repo_rel(output_dir / "control_contracts.json"),
            "component_control_validation": repo_rel(output_dir / "component_control_validation.csv"),
            "signature_collision_v2": repo_rel(output_dir / "signature_collision_v2.json"),
            "power_gate": repo_rel(output_dir / "power_gate.json"),
            "discriminative_validity_report": repo_rel(output_dir / "discriminative_validity_report.md"),
            "phase_2_2_result": repo_rel(output_dir / "phase_2_2_result.json"),
        },
    }
    result["result_hash"] = result_hash_payload(result)
    report = render_report(result, power, validation_rows, collisions)

    if write_outputs:
        write_json(output_dir / "precommitment_contract.json", contract)
        write_json(output_dir / "signature_schema.json", contract["fine_signature_schema"])
        write_json(
            output_dir / "fine_signatures.json",
            {
                "contract_hash": contract_hash,
                "hash_policy": "fine_signature_record_id names the representation only",
                "fine_signatures": fine_signatures,
            },
        )
        write_json(output_dir / "equivalence_tolerances.json", contract["equivalence"])
        write_csv(
            output_dir / "equivalence_results.csv",
            equivalence,
            [
                "pair_type",
                "mode",
                "left_observation_id",
                "right_observation_id",
                "left_group",
                "right_group",
                "equivalent",
                "available_component_count",
                "failure_reasons",
                "boundary_hits",
                "boundary_sensitive",
                "hash_equality_used_for_equivalence",
                "component_results_json",
            ],
        )
        write_json(output_dir / "control_contracts.json", contract["component_control_contracts"])
        write_csv(
            output_dir / "component_control_validation.csv",
            validation_rows,
            [
                "control_label",
                "control_class",
                "component",
                "expectation",
                "target_target_mean",
                "target_control_mean",
                "shift_control_minus_target",
                "control_distance_from_source",
                "missing_count",
                "passed",
                "invalidation_condition",
                "invalidation_reason",
                "expected_affected_metrics",
                "expected_unaffected_metrics",
            ],
        )
        write_json(output_dir / "signature_collision_v2.json", collisions)
        write_json(output_dir / "power_gate.json", power)
        (output_dir / "discriminative_validity_report.md").write_text(report, encoding="utf-8")
        write_json(output_dir / "phase_2_2_result.json", result)
    return {
        "contract": contract,
        "contract_hash": contract_hash,
        "benchmark": benchmark,
        "observations": observations,
        "pair_sets": pair_sets,
        "fine_signatures": fine_signatures,
        "equivalence": equivalence,
        "validation_rows": validation_rows,
        "collisions": collisions,
        "power": power,
        "result": result,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run CST Phase 2.2 discriminative instrument repair.")
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_REPAIR_OUTPUT_DIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    run = run_phase_2_2_repair(contract_path=args.contract, output_dir=args.output_dir, write_outputs=True)
    print(
        json.dumps(
            {
                "labels": run["result"]["labels"],
                "result_hash": run["result"]["result_hash"],
                "result": run["result"]["outputs"]["phase_2_2_result"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
