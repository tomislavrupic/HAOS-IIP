from __future__ import annotations

import statistics
from typing import Any

from experiments.cst.cst_io import stable_hash
from experiments.cst.cst_types import (
    ClosureDistanceComponents,
    ClosureSignature,
    ComponentDistance,
    RecoverabilityVector,
    RecoveryObservation,
)


COMPONENT_NAMES = ("d_inv", "d_struct", "d_spec", "d_causal", "d_time", "d_addr")


def as_float(value: Any, default: float | None = None) -> float | None:
    if value in (None, ""):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def population_std(values: list[float]) -> float:
    if not values:
        return 0.0
    if len(values) == 1:
        return 0.0
    return statistics.pstdev(values)


def safe_numeric_distance(left: float, right: float, *, scale: float | None = None) -> float:
    denominator = max(abs(left), abs(right), 1.0e-12) if scale is None else max(float(scale), 1.0e-12)
    return min(abs(left - right) / denominator, 1.0)


def numeric_dict_distance(
    left: dict[str, Any],
    right: dict[str, Any],
    keys: list[str],
    *,
    scale: float | None = None,
) -> float | None:
    distances: list[float] = []
    for key in keys:
        left_value = as_float(left.get(key))
        right_value = as_float(right.get(key))
        if left_value is None or right_value is None:
            continue
        distances.append(safe_numeric_distance(left_value, right_value, scale=scale))
    if not distances:
        return None
    return mean(distances)


def jaccard_distance(left_values: list[str], right_values: list[str]) -> float:
    left = set(left_values)
    right = set(right_values)
    union = left | right
    if not union:
        return 0.0
    return 1.0 - (len(left & right) / len(union))


def sequence_rank_distance(left_values: list[str], right_values: list[str]) -> float:
    if not left_values or not right_values:
        return 0.0 if left_values == right_values else 1.0
    labels = sorted(set(left_values) | set(right_values))
    max_rank = max(len(labels) - 1, 1)
    left_positions = {label: index for index, label in enumerate(left_values)}
    right_positions = {label: index for index, label in enumerate(right_values)}
    distances = []
    for label in labels:
        left_rank = left_positions.get(label, max_rank)
        right_rank = right_positions.get(label, max_rank)
        distances.append(abs(left_rank - right_rank) / max_rank)
    return min(mean(distances), 1.0)


def unavailable_component(name: str, source_fields: list[str], reason: str) -> ComponentDistance:
    return ComponentDistance(
        name=name,
        value=None,
        availability="unavailable",
        normalization_method="not-computed",
        source_fields=source_fields,
        warnings=[reason],
    )


def available_component(
    name: str,
    value: float,
    source_fields: list[str],
    method: str,
    warnings: list[str] | None = None,
) -> ComponentDistance:
    return ComponentDistance(
        name=name,
        value=max(0.0, min(float(value), 1.0)),
        availability="available",
        normalization_method=method,
        source_fields=source_fields,
        warnings=list(warnings or []),
    )


def component_d_inv(left: ClosureSignature, right: ClosureSignature) -> ComponentDistance:
    left_summary = left.fields.get("invariant_summary")
    right_summary = right.fields.get("invariant_summary")
    source = ["invariant_summary"]
    if not isinstance(left_summary, dict) or not isinstance(right_summary, dict):
        return unavailable_component("d_inv", source, "invariant_summary missing from one or both signatures")
    value = numeric_dict_distance(
        left_summary,
        right_summary,
        ["terminal_low_k_fraction", "transport_character_index", "shell_latency_span"],
    )
    if value is None:
        return unavailable_component("d_inv", source, "no comparable numeric invariant fields")
    return available_component("d_inv", value, source, "mean bounded relative numeric difference")


def component_d_struct(left: ClosureSignature, right: ClosureSignature) -> ComponentDistance:
    left_summary = left.fields.get("influence_edge_summary")
    right_summary = right.fields.get("influence_edge_summary")
    source = ["influence_edge_summary.edges", "influence_edge_summary.edge_count"]
    if not isinstance(left_summary, dict) or not isinstance(right_summary, dict):
        return unavailable_component("d_struct", source, "influence_edge_summary missing from one or both signatures")
    left_edges = [str(edge) for edge in left_summary.get("edges", [])]
    right_edges = [str(edge) for edge in right_summary.get("edges", [])]
    edge_distance = jaccard_distance(left_edges, right_edges)
    count_distance = numeric_dict_distance(left_summary, right_summary, ["edge_count", "stable_edge_count"])
    values = [edge_distance]
    if count_distance is not None:
        values.append(count_distance)
    return available_component("d_struct", mean(values), source, "edge-set Jaccard plus bounded edge-count drift")


def component_d_spec(left: ClosureSignature, right: ClosureSignature) -> ComponentDistance:
    left_summary = left.fields.get("spectral_low_mode_summary")
    right_summary = right.fields.get("spectral_low_mode_summary")
    source = ["spectral_low_mode_summary"]
    if not isinstance(left_summary, dict) or not isinstance(right_summary, dict):
        return unavailable_component("d_spec", source, "spectral_low_mode_summary missing from one or both signatures")
    value = numeric_dict_distance(
        left_summary,
        right_summary,
        ["terminal_low_k_fraction", "terminal_low_k_power", "terminal_fluctuation_slope"],
    )
    if value is None:
        return unavailable_component("d_spec", source, "no comparable spectral or low-mode fields")
    return available_component("d_spec", value, source, "mean bounded relative low-mode descriptor difference")


def component_d_causal(left: ClosureSignature, right: ClosureSignature) -> ComponentDistance:
    left_summary = left.fields.get("causal_summary")
    right_summary = right.fields.get("causal_summary")
    source = ["causal_summary"]
    if not isinstance(left_summary, dict) or not isinstance(right_summary, dict):
        return unavailable_component("d_causal", source, "causal_summary missing from one or both signatures")
    depth_left = left_summary.get("front_depth_profile", {})
    depth_right = right_summary.get("front_depth_profile", {})
    depth_distance = numeric_dict_distance(depth_left, depth_right, ["front_near_depth", "front_far_depth"], scale=4.0)
    causal_distance = numeric_dict_distance(
        left_summary,
        right_summary,
        ["acyclicity_score", "mean_causal_depth", "order_compatibility_score", "mismatch_rate"],
    )
    values = [value for value in [depth_distance, causal_distance] if value is not None]
    if not values:
        return unavailable_component("d_causal", source, "no comparable causal fields")
    return available_component("d_causal", mean(values), source, "bounded causal-depth, acyclicity, and order-compatibility differences")


def component_d_time(left: ClosureSignature, right: ClosureSignature) -> ComponentDistance:
    left_summary = left.fields.get("event_chain_summary")
    right_summary = right.fields.get("event_chain_summary")
    source = ["event_chain_summary.chain_signature", "event_chain_summary.event_times"]
    if not isinstance(left_summary, dict) or not isinstance(right_summary, dict):
        return unavailable_component("d_time", source, "event_chain_summary missing from one or both signatures")
    left_chain = str(left_summary.get("chain_signature", "")).split(">") if left_summary.get("chain_signature") else []
    right_chain = str(right_summary.get("chain_signature", "")).split(">") if right_summary.get("chain_signature") else []
    rank_distance = sequence_rank_distance(left_chain, right_chain)
    left_times = left_summary.get("event_times", {})
    right_times = right_summary.get("event_times", {})
    time_distance = numeric_dict_distance(left_times, right_times, sorted(set(left_times) | set(right_times)), scale=4.2)
    values = [rank_distance]
    if time_distance is not None:
        values.append(time_distance)
    return available_component("d_time", mean(values), source, "event-chain rank distance plus tau-grid normalized event-time drift")


def component_d_addr(left: ClosureSignature, right: ClosureSignature) -> ComponentDistance:
    source = ["branch_signature_id", "address_summary"]
    address_left = left.fields.get("address_summary")
    address_right = right.fields.get("address_summary")
    if address_left is None or address_right is None:
        return unavailable_component("d_addr", source, "address_summary missing from one or both signatures")
    if left.branch_signature_id == right.branch_signature_id:
        return available_component("d_addr", 0.0, source, "exact operational branch-signature match")
    prefix_left = [str(value) for value in address_left.get("dominant_chain_prefix", [])] if isinstance(address_left, dict) else []
    prefix_right = [str(value) for value in address_right.get("dominant_chain_prefix", [])] if isinstance(address_right, dict) else []
    prefix_distance = sequence_rank_distance(prefix_left, prefix_right)
    depth_distance = None
    if isinstance(address_left, dict) and isinstance(address_right, dict):
        depth_distance = numeric_dict_distance(
            address_left.get("front_depth_profile", {}),
            address_right.get("front_depth_profile", {}),
            ["front_near_depth", "front_far_depth"],
            scale=4.0,
        )
    values = [prefix_distance]
    if depth_distance is not None:
        values.append(depth_distance)
    return available_component("d_addr", max(mean(values), 0.25), source, "operational address mismatch from dominant prefix/depth fields")


def compute_components(left: ClosureSignature, right: ClosureSignature) -> dict[str, ComponentDistance]:
    return {
        "d_inv": component_d_inv(left, right),
        "d_struct": component_d_struct(left, right),
        "d_spec": component_d_spec(left, right),
        "d_causal": component_d_causal(left, right),
        "d_time": component_d_time(left, right),
        "d_addr": component_d_addr(left, right),
    }


def normalized_weights(
    components: dict[str, ComponentDistance],
    weights: dict[str, float],
    enabled: list[str],
) -> tuple[dict[str, float], list[str], list[str]]:
    active: dict[str, float] = {}
    disabled: list[str] = []
    warnings: list[str] = []
    for name in COMPONENT_NAMES:
        component = components[name]
        if name not in enabled:
            disabled.append(name)
            continue
        if component.availability != "available" or component.value is None:
            disabled.append(name)
            warnings.append(f"{name} unavailable and excluded from scalar normalization")
            continue
        raw_weight = float(weights.get(name, 0.0))
        if raw_weight <= 0.0:
            disabled.append(name)
            warnings.append(f"{name} disabled by non-positive weight")
            continue
        active[name] = raw_weight
    total = sum(active.values())
    if total <= 0.0:
        return {}, disabled, warnings + ["no available enabled components for scalar distance"]
    return {name: value / total for name, value in active.items()}, disabled, warnings


def compute_closure_distance(
    left: RecoveryObservation,
    right: RecoveryObservation,
    *,
    weights: dict[str, float],
    enabled_components: list[str],
    scalar_enabled: bool,
) -> ClosureDistanceComponents:
    components = compute_components(left.closure_signature, right.closure_signature)
    used_weights, disabled, warnings = normalized_weights(components, weights, enabled_components)
    scalar = None
    if scalar_enabled and used_weights:
        scalar = sum(float(components[name].value) * weight for name, weight in used_weights.items())
    return ClosureDistanceComponents(
        left_observation_id=left.observation_id,
        right_observation_id=right.observation_id,
        left_signature_id=left.closure_signature.branch_signature_id,
        right_signature_id=right.closure_signature.branch_signature_id,
        components=components,
        scalar_distance=scalar,
        scalar_enabled=bool(scalar_enabled),
        weights_used=used_weights,
        disabled_components=disabled,
        warnings=warnings,
    )


def component_dominance(distance: ClosureDistanceComponents) -> float:
    if distance.scalar_distance is None or distance.scalar_distance <= 1.0e-12:
        return 0.0
    shares = []
    for name, weight in distance.weights_used.items():
        component = distance.components[name]
        if component.value is None:
            continue
        shares.append((component.value * weight) / distance.scalar_distance)
    return max(shares) if shares else 0.0


def scalar_values(distances: list[ClosureDistanceComponents]) -> list[float]:
    return [float(item.scalar_distance) for item in distances if item.scalar_distance is not None]


def compute_weight_sensitivity(
    target_pairs: list[tuple[RecoveryObservation, RecoveryObservation]],
    control_pairs: list[tuple[RecoveryObservation, RecoveryObservation]],
    *,
    weight_families: dict[str, dict[str, float]],
    enabled_components: list[str],
) -> dict[str, Any]:
    rows: dict[str, Any] = {}
    for family_name, weights in weight_families.items():
        target_distances = [
            compute_closure_distance(
                left,
                right,
                weights=weights,
                enabled_components=enabled_components,
                scalar_enabled=True,
            )
            for left, right in target_pairs
        ]
        control_distances = [
            compute_closure_distance(
                left,
                right,
                weights=weights,
                enabled_components=enabled_components,
                scalar_enabled=True,
            )
            for left, right in control_pairs
        ]
        target_scalars = scalar_values(target_distances)
        control_scalars = scalar_values(control_distances)
        rows[family_name] = {
            "target_mean": mean(target_scalars),
            "control_mean": mean(control_scalars),
            "target_std": population_std(target_scalars),
            "control_std": population_std(control_scalars),
            "target_better_than_control": bool(target_scalars and control_scalars and mean(target_scalars) < mean(control_scalars)),
        }
    return rows


def ablation_balance(
    target_pairs: list[tuple[RecoveryObservation, RecoveryObservation]],
    control_pairs: list[tuple[RecoveryObservation, RecoveryObservation]],
    *,
    base_weights: dict[str, float],
    enabled_components: list[str],
) -> tuple[float, dict[str, Any]]:
    outcomes: dict[str, Any] = {}
    checks = []
    for disabled_component in enabled_components:
        ablated_enabled = [name for name in enabled_components if name != disabled_component]
        target_distances = [
            compute_closure_distance(
                left,
                right,
                weights=base_weights,
                enabled_components=ablated_enabled,
                scalar_enabled=True,
            )
            for left, right in target_pairs
        ]
        control_distances = [
            compute_closure_distance(
                left,
                right,
                weights=base_weights,
                enabled_components=ablated_enabled,
                scalar_enabled=True,
            )
            for left, right in control_pairs
        ]
        target_mean = mean(scalar_values(target_distances))
        control_mean = mean(scalar_values(control_distances))
        survived = bool(target_mean < control_mean)
        checks.append(survived)
        outcomes[disabled_component] = {
            "target_mean": target_mean,
            "control_mean": control_mean,
            "target_better_than_control": survived,
        }
    return (sum(1 for item in checks if item) / len(checks) if checks else 0.0), outcomes


def selectivity_payload(
    *,
    recovery_rate: float,
    branch_identity_persistence: float,
    closure_fidelity: float,
    control_mean_distance: float,
    target_distances: list[float],
    epsilons: list[float],
    default_epsilon: float,
) -> dict[str, Any]:
    recoverable_signal_strength = float(recovery_rate * branch_identity_persistence * closure_fidelity)
    noise_floor = max(0.0, 1.0 - float(control_mean_distance))
    sensitivity = []
    for epsilon in epsilons:
        sensitivity.append(
            {
                "epsilon": float(epsilon),
                "generalized_cst_selectivity": recoverable_signal_strength / max(noise_floor, float(epsilon)),
            }
        )
    return {
        "label": "generalized CST selectivity",
        "not_physical_q": True,
        "recoverable_signal_strength": recoverable_signal_strength,
        "noise_floor": noise_floor,
        "epsilon": float(default_epsilon),
        "generalized_cst_selectivity": recoverable_signal_strength / max(noise_floor, float(default_epsilon)),
        "signal_variance": population_std([1.0 - value for value in target_distances]),
        "epsilon_sensitivity": sensitivity,
    }


def recoverability_vector(
    *,
    target_distances: list[ClosureDistanceComponents],
    target_signature_matches: list[bool],
    control_distances: list[ClosureDistanceComponents],
    thresholds: dict[str, Any],
    ablation_score: float,
    selectivity: dict[str, Any] | None,
    scalar_enabled: bool,
) -> RecoverabilityVector:
    target_scalars = scalar_values(target_distances)
    control_scalars = scalar_values(control_distances)
    closure_threshold = float(thresholds["closure_distance_max"])
    variance_threshold = float(thresholds["target_distance_std_max"])
    recovery_rate = mean([1.0 if value <= closure_threshold else 0.0 for value in target_scalars])
    identity_persistence = mean([1.0 if item else 0.0 for item in target_signature_matches])
    closure_fidelity = max(0.0, 1.0 - mean(target_scalars))
    control_margin = mean(control_scalars) - mean(target_scalars)
    variance_value = population_std(target_scalars)
    variance_bound = max(0.0, 1.0 - min(variance_value / max(variance_threshold, 1.0e-12), 1.0))
    optional_scalar = None
    formula = None
    warnings: list[str] = []
    if scalar_enabled:
        optional_scalar = mean([
            recovery_rate,
            identity_persistence,
            closure_fidelity,
            max(0.0, min(control_margin, 1.0)),
            variance_bound,
            ablation_score,
        ])
        formula = "mean(recovery_rate, branch_identity_persistence, closure_fidelity, clipped_control_margin, variance_bound, ablation_balance)"
        warnings.append("scalar compression is heuristic; vector fields are authoritative")
    return RecoverabilityVector(
        recovery_rate=recovery_rate,
        branch_identity_persistence=identity_persistence,
        closure_fidelity=closure_fidelity,
        control_margin=control_margin,
        variance_bound=variance_bound,
        ablation_balance=ablation_score,
        optional_scalar=optional_scalar,
        scalar_formula=formula,
        scalar_warnings=warnings,
        generalized_cst_selectivity=selectivity,
    )


def verdict_from_metrics(
    *,
    vector: RecoverabilityVector,
    target_distances: list[ClosureDistanceComponents],
    control_distances: list[ClosureDistanceComponents],
    target_seed_count: int,
    thresholds: dict[str, Any],
    weight_sensitivity: dict[str, Any],
    control_family_stats: dict[str, Any] | None = None,
    unavailable_controls: dict[str, Any],
    forced_fail_probe_passed: bool,
) -> tuple[str, list[str], list[str]]:
    reasons: list[str] = []
    warnings: list[str] = []
    verdict = "PASS"

    target_scalars = scalar_values(target_distances)
    control_scalars = scalar_values(control_distances)
    if not target_scalars:
        return "FAIL", ["no target closure distances were computable"], warnings
    if not control_scalars:
        return "FAIL", ["no negative-control closure distances were computable"], warnings
    if not forced_fail_probe_passed:
        return "FAIL", ["forced-fail probe did not produce FAIL"], warnings

    if vector.recovery_rate < 1.0:
        verdict = "FAIL"
        reasons.append("independent target runs exceed the closure-distance threshold")
    if vector.branch_identity_persistence < float(thresholds["branch_identity_persistence_min"]):
        verdict = "FAIL"
        reasons.append("independent target runs do not recover the same operational branch_signature_id")
    if mean(control_scalars) <= mean(target_scalars):
        verdict = "FAIL"
        reasons.append("negative controls perform equally or better than target under equal-weight scalar distance")
    elif vector.control_margin < float(thresholds["control_margin_min"]):
        verdict = "OPEN"
        reasons.append("target/control margin is positive but below configured threshold")

    for label, stats in sorted((control_family_stats or {}).items()):
        if stats.get("availability") != "available" or stats.get("role") != "negative":
            continue
        family_mean = stats.get("mean_distance")
        if family_mean is None:
            continue
        if float(family_mean) <= mean(target_scalars):
            verdict = "FAIL"
            reasons.append(f"negative control family {label} performs equally or better than target mean")

    target_std = population_std(target_scalars)
    if target_std > float(thresholds["target_distance_std_max"]):
        verdict = "FAIL"
        reasons.append("target closure-distance variance exceeds configured threshold")

    max_dominance = max([component_dominance(item) for item in target_distances + control_distances] or [0.0])
    if max_dominance > float(thresholds["max_component_dominance"]):
        verdict = "FAIL"
        reasons.append("one component dominates scalar distance beyond configured limit")

    if vector.ablation_balance < float(thresholds["ablation_balance_min"]):
        verdict = "FAIL"
        reasons.append("ablation reverses target/control ordering for at least one component")

    if any(not row.get("target_better_than_control") for row in weight_sensitivity.values()):
        verdict = "FAIL"
        reasons.append("weight-sensitivity family reverses or erases target/control ordering")

    if target_seed_count < int(thresholds["min_seed_count_for_pass"]):
        if verdict == "PASS":
            verdict = "OPEN"
        reasons.append("available seed count is below configured minimum for PASS")

    if unavailable_controls:
        if verdict == "PASS":
            verdict = "OPEN"
        warnings.append("one or more requested controls are unavailable in the frozen artifact slice")

    if verdict == "PASS":
        reasons.append("all configured CST toy gates passed")
    return verdict, reasons, warnings


def forced_fail_probe() -> bool:
    vector = RecoverabilityVector(
        recovery_rate=0.0,
        branch_identity_persistence=0.0,
        closure_fidelity=0.0,
        control_margin=-1.0,
        variance_bound=0.0,
        ablation_balance=0.0,
    )
    verdict, _, _ = verdict_from_metrics(
        vector=vector,
        target_distances=[],
        control_distances=[],
        target_seed_count=0,
        thresholds={
            "closure_distance_max": 0.0,
            "branch_identity_persistence_min": 1.0,
            "control_margin_min": 1.0,
            "target_distance_std_max": 0.0,
            "max_component_dominance": 0.0,
            "ablation_balance_min": 1.0,
            "min_seed_count_for_pass": 99,
        },
        weight_sensitivity={"equal": {"target_better_than_control": False}},
        control_family_stats={},
        unavailable_controls={},
        forced_fail_probe_passed=True,
    )
    return verdict == "FAIL"


def branch_signature_id(signature_fields: dict[str, Any], *, float_precision: int = 6) -> str:
    address_summary = signature_fields.get("address_summary", {})
    return stable_hash(address_summary, prefix="bsig_", float_precision=float_precision)
