from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path
from typing import Any

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from experiments.cst.cst_controls import NEGATIVE_CONTROL_LABELS, build_derived_control
from experiments.cst.cst_io import (
    REPO_ROOT,
    artifact_hashes,
    canonical_json,
    code_version,
    read_csv,
    read_json,
    repo_rel,
    stable_hash,
    write_csv,
    write_json,
)
from experiments.cst.cst_metrics import (
    ablation_balance,
    as_float,
    branch_signature_id,
    compute_closure_distance,
    compute_weight_sensitivity,
    forced_fail_probe,
    mean,
    population_std,
    recoverability_vector,
    scalar_values,
    selectivity_payload,
    verdict_from_metrics,
)
from experiments.cst.cst_types import (
    CSTBenchmarkResult,
    CSTRunProvenance,
    ClosureDistanceComponents,
    ClosureSignature,
    PerturbationSpec,
    RecoveryObservation,
)


HERE = Path(__file__).resolve().parent
DEFAULT_CONFIG_PATH = HERE / "configs" / "cst_toy_config.json"
DEFAULT_OUTPUT_DIR = HERE / "runs"

PHASE_ARTIFACTS = {
    "phase15_effective_speed": REPO_ROOT / "phase15-propagation" / "runs" / "phase15_effective_speed_ledger.csv",
    "phase15_influence_range": REPO_ROOT / "phase15-propagation" / "runs" / "phase15_influence_range_ledger.csv",
    "phase15_propagation": REPO_ROOT / "phase15-propagation" / "runs" / "phase15_propagation_ledger.csv",
    "phase15_transport_descriptor": REPO_ROOT / "phase15-propagation" / "runs" / "phase15_transport_descriptor_ledger.csv",
    "phase16_event_ordering": REPO_ROOT / "phase16-temporal-ordering" / "runs" / "phase16_event_ordering_ledger.csv",
    "phase16_front_arrival_ordering": REPO_ROOT / "phase16-temporal-ordering" / "runs" / "phase16_front_arrival_ordering.csv",
    "phase16_ordering_robustness": REPO_ROOT / "phase16-temporal-ordering" / "runs" / "phase16_ordering_robustness.csv",
    "phase17_causal_distance_metrics": REPO_ROOT / "phase17-causal-closure" / "runs" / "phase17_causal_distance_metrics.csv",
    "phase17_dag_statistics": REPO_ROOT / "phase17-causal-closure" / "runs" / "phase17_dag_statistics.csv",
    "phase17_influence_graph": REPO_ROOT / "phase17-causal-closure" / "runs" / "phase17_influence_graph_ledger.csv",
    "phase17_propagation_order_compatibility": REPO_ROOT / "phase17-causal-closure" / "runs" / "phase17_propagation_order_compatibility.csv",
    "phase18_distance_surrogate": REPO_ROOT / "phase18-distance-surrogate" / "runs" / "phase18_distance_surrogate_ledger.csv",
    "phase18_refinement_scaling": REPO_ROOT / "phase18-distance-surrogate" / "runs" / "phase18_refinement_scaling.csv",
    "phase18_shell_ordering_metrics": REPO_ROOT / "phase18-distance-surrogate" / "runs" / "phase18_shell_ordering_metrics.csv",
    "phase18_triangle_violation_rate": REPO_ROOT / "phase18-distance-surrogate" / "runs" / "phase18_triangle_violation_rate.csv",
}


def as_bool(value: Any) -> bool:
    return str(value).strip().lower() == "true"


def as_int(value: Any) -> int:
    return int(str(value).strip())


def quantize(value: float | None, quantum: float, *, digits: int = 6) -> float | None:
    if value is None:
        return None
    if quantum <= 0.0:
        return round(float(value), digits)
    return round(round(float(value) / quantum) * quantum, digits)


def row_float(row: dict[str, str] | None, key: str, default: float | None = None) -> float | None:
    if row is None:
        return default
    return as_float(row.get(key), default)


def chain_labels(chain_signature: str) -> list[str]:
    return [item for item in chain_signature.split(">") if item]


def edges_from_labels(labels: list[str]) -> list[str]:
    return [f"{labels[index]}->{labels[index + 1]}" for index in range(len(labels) - 1)]


def event_times_from_chain(
    labels: list[str],
    shell_row: dict[str, str],
    distance_rows: list[dict[str, str]],
    *,
    float_precision: int,
) -> dict[str, float]:
    explicit: dict[str, float] = {}
    for row in distance_rows:
        shell_name = row.get("shell_name")
        value = as_float(row.get("propagation_time"))
        if shell_name and value is not None:
            explicit[str(shell_name)] = value
    near = row_float(shell_row, "near_mean_arrival", 0.0) or 0.0
    far = row_float(shell_row, "far_mean_arrival", 0.0) or 0.0
    explicit.setdefault("front_near", near)
    explicit.setdefault("front_far", far)
    recovery_time = max([near, far] + list(explicit.values()) + [0.0])
    if not labels:
        return {}
    denominator = max(len(labels) - 1, 1)
    event_times: dict[str, float] = {}
    for index, label in enumerate(labels):
        if label in explicit:
            event_times[label] = round(explicit[label], float_precision)
        else:
            event_times[label] = round((index / denominator) * recovery_time, float_precision)
    return event_times


def index_ledgers(ledgers: dict[str, list[dict[str, str]]]) -> dict[str, Any]:
    indexes: dict[str, Any] = {}

    def seed_key(row: dict[str, str]) -> tuple[str, int, int, int, str]:
        return (
            row["hierarchy_label"],
            as_int(row["n_side"]),
            as_int(row["ensemble_size"]),
            as_int(row["seed"]),
            row["probe_name"],
        )

    def summary_key(row: dict[str, str]) -> tuple[str, int]:
        return (row["hierarchy_label"], as_int(row["n_side"]))

    def scaling_key(row: dict[str, str]) -> tuple[str, int, int, str]:
        return (
            row["hierarchy_label"],
            as_int(row["n_side"]),
            as_int(row["ensemble_size"]),
            row["probe_name"],
        )

    for name in ("phase15_transport_descriptor", "phase18_shell_ordering_metrics"):
        indexes[name] = {seed_key(row): row for row in ledgers[name]}
    indexes["phase18_refinement_scaling"] = {scaling_key(row): row for row in ledgers["phase18_refinement_scaling"]}
    for name in (
        "phase17_influence_graph",
        "phase17_dag_statistics",
        "phase17_causal_distance_metrics",
        "phase17_propagation_order_compatibility",
    ):
        indexes[name] = {summary_key(row): row for row in ledgers[name]}

    distance_groups: dict[tuple[str, int, int, int, str], list[dict[str, str]]] = {}
    for row in ledgers["phase18_distance_surrogate"]:
        distance_groups.setdefault(seed_key(row), []).append(row)
    indexes["phase18_distance_surrogate"] = distance_groups
    return indexes


def source_artifact_paths() -> dict[str, str]:
    return {name: repo_rel(path) for name, path in sorted(PHASE_ARTIFACTS.items())}


def make_provenance(
    *,
    config: dict[str, Any],
    config_hash: str,
    source_hashes: dict[str, str],
    hierarchy_label: str,
    n_side: int,
    seed: int,
    control_group: str,
) -> CSTRunProvenance:
    payload = {
        "source_artifact_hashes": source_hashes,
        "config_hash": config_hash,
        "config_version": config["config_version"],
        "candidate_label": config["candidate_label"],
        "hierarchy_label": hierarchy_label,
        "n_side": n_side,
        "ensemble_size": config["ensemble_size"],
        "seed": seed,
        "probe": config["probe"],
        "perturbation": config["probe"],
        "control_group": control_group,
    }
    return CSTRunProvenance(
        run_instance_id=stable_hash(payload, prefix="run_"),
        source_artifacts=source_artifact_paths(),
        source_artifact_hashes=source_hashes,
        config_version=config["config_version"],
        config_hash=config_hash,
        code_version=code_version(),
        candidate_label=config["candidate_label"],
        hierarchy_label=hierarchy_label,
        n_side=n_side,
        ensemble_size=int(config["ensemble_size"]),
        seed=seed,
        probe=config["probe"],
        perturbation=config["probe"],
        control_group=control_group,
    )


def extract_closure_signature(
    *,
    hierarchy_label: str,
    n_side: int,
    seed: int,
    config: dict[str, Any],
    indexes: dict[str, Any],
    source_artifacts: list[str],
) -> ClosureSignature:
    probe = config["probe"]
    ensemble_size = int(config["ensemble_size"])
    signature_config = config["signature"]
    float_precision = int(signature_config["float_precision"])
    key = (hierarchy_label, int(n_side), ensemble_size, int(seed), probe)
    summary_key = (hierarchy_label, int(n_side))
    scaling_key = (hierarchy_label, int(n_side), ensemble_size, probe)

    try:
        shell_row = indexes["phase18_shell_ordering_metrics"][key]
        transport_row = indexes["phase15_transport_descriptor"][key]
        refinement_row = indexes["phase18_refinement_scaling"][scaling_key]
        influence_row = indexes["phase17_influence_graph"][summary_key]
        dag_row = indexes["phase17_dag_statistics"][summary_key]
        causal_row = indexes["phase17_causal_distance_metrics"][summary_key]
        compatibility_row = indexes["phase17_propagation_order_compatibility"][summary_key]
    except KeyError as exc:
        raise KeyError(f"missing frozen CST input for {key}: {exc}") from exc

    distance_rows = indexes["phase18_distance_surrogate"].get(key, [])
    chain_signature = shell_row["chain_signature"]
    labels = chain_labels(chain_signature)
    prefix_len = int(signature_config["dominant_prefix_length"])
    event_times = event_times_from_chain(labels, shell_row, distance_rows, float_precision=float_precision)

    front_depth_profile = {
        "front_near_depth": row_float(causal_row, "front_near_depth"),
        "front_far_depth": row_float(causal_row, "front_far_depth"),
    }
    order_compatibility = row_float(refinement_row, "order_compatibility_score")
    low_mode_fraction = row_float(transport_row, "terminal_low_k_fraction")

    distance_summary = {
        row["shell_name"]: {
            "causal_depth": row_float(row, "causal_depth"),
            "propagation_time": row_float(row, "propagation_time"),
            "surrogate_distance": row_float(row, "surrogate_distance"),
            "front_distance": row_float(row, "front_distance"),
            "radius_fraction": row_float(row, "radius_fraction"),
            "monotonicity_score": row_float(row, "monotonicity_score"),
        }
        for row in distance_rows
    }

    fields: dict[str, Any] = {
        "candidate_branch": {
            "stable_identifier_rule": "branch_signature_id = sha256(address_summary)",
            "state_or_feature_representation": "phase18 chain signature plus phase15 low-mode and phase17 causal summaries",
            "control_group_label": hierarchy_label,
        },
        "parent_or_predecessor_relation": {
            "artifact_chain": [
                "phase15_transport_descriptor",
                "phase17_causal_closure",
                "phase18_distance_surrogate",
            ],
            "relation_note": "same hierarchy/n_side/probe with independent seeds; no dynamic parent inferred",
        },
        "perturbation_history": [
            {
                "name": probe,
                "family": "frozen_phase18_probe",
                "source": "phase18_shell_ordering_metrics",
            }
        ],
        "recovery_history": {
            "recovery_time": max(event_times.values()) if event_times else None,
            "recovery_time_source": "phase18 arrival fields with rank interpolation for non-front events",
            "chain_signature": chain_signature,
        },
        "event_chain_summary": {
            "chain_signature": chain_signature,
            "dominant_chain_prefix": labels[:prefix_len],
            "event_times": event_times,
        },
        "influence_edge_summary": {
            "edges": edges_from_labels(labels),
            "edge_count": row_float(influence_row, "edge_count"),
            "stable_edge_count": row_float(influence_row, "stable_edge_count"),
            "mean_edge_reproducibility": row_float(influence_row, "mean_edge_reproducibility"),
            "edge_set_distance_from_prev": row_float(influence_row, "edge_set_distance_from_prev"),
        },
        "invariant_summary": {
            "terminal_low_k_fraction": low_mode_fraction,
            "terminal_low_k_power": row_float(transport_row, "terminal_low_k_power"),
            "transport_character_index": row_float(transport_row, "transport_character_index"),
            "shell_latency_span": row_float(transport_row, "shell_latency_span"),
            "shell_order_score": row_float(shell_row, "shell_order_score"),
        },
        "spectral_low_mode_summary": {
            "terminal_low_k_fraction": low_mode_fraction,
            "terminal_low_k_power": row_float(transport_row, "terminal_low_k_power"),
            "terminal_fluctuation_slope": row_float(transport_row, "terminal_fluctuation_slope"),
            "low_mode_band": quantize(
                low_mode_fraction,
                float(signature_config["low_mode_quantum"]),
                digits=float_precision,
            ),
        },
        "causal_summary": {
            "front_depth_profile": front_depth_profile,
            "acyclicity_score": row_float(dag_row, "acyclicity_score"),
            "cycle_count": row_float(dag_row, "cycle_count"),
            "cycle_edge_count": row_float(dag_row, "cycle_edge_count"),
            "mean_causal_depth": row_float(causal_row, "mean_causal_depth"),
            "reachable_node_count": row_float(causal_row, "reachable_node_count"),
            "order_compatibility_score": order_compatibility,
            "mismatch_rate": row_float(compatibility_row, "mismatch_rate"),
        },
        "shell_ordering_summary": {
            "monotonic_shell_ordering": as_bool(shell_row.get("monotonic_shell_ordering")),
            "max_shell_overlap_fraction": row_float(shell_row, "max_shell_overlap_fraction"),
            "shell_order_score": row_float(shell_row, "shell_order_score"),
            "near_mean_arrival": row_float(shell_row, "near_mean_arrival"),
            "far_mean_arrival": row_float(shell_row, "far_mean_arrival"),
            "phase17_mismatch_rate": row_float(shell_row, "phase17_mismatch_rate"),
        },
        "distance_surrogate_summary": distance_summary,
        "address_summary": {
            "dominant_chain_prefix": labels[:prefix_len],
            "front_depth_profile": front_depth_profile,
            "monotonic_shell_ordering": as_bool(shell_row.get("monotonic_shell_ordering")),
            "order_compatibility_band": quantize(
                order_compatibility,
                float(signature_config["order_compatibility_quantum"]),
                digits=float_precision,
            ),
        },
        "warnings": [
            "event_times for non-front chain elements are rank-interpolated from phase18 recovery time"
        ],
    }
    signature_id = branch_signature_id(fields, float_precision=float_precision)
    return ClosureSignature(
        branch_signature_id=signature_id,
        signature_version=signature_config["signature_version"],
        operational_equivalence_note=(
            "Experimental CST branch identity: stable hash of a compact recovered-structure address. "
            "It does not assert fundamental identity."
        ),
        fields=fields,
        unavailable_components={},
        source_artifacts=source_artifacts,
        warnings=list(fields["warnings"]),
    )


def make_observation(
    *,
    hierarchy_label: str,
    n_side: int,
    seed: int,
    observation_kind: str,
    control_role: str,
    control_group_label: str,
    config: dict[str, Any],
    indexes: dict[str, Any],
    config_hash: str,
    source_hashes: dict[str, str],
    source_artifacts: list[str],
) -> RecoveryObservation:
    signature = extract_closure_signature(
        hierarchy_label=hierarchy_label,
        n_side=n_side,
        seed=seed,
        config=config,
        indexes=indexes,
        source_artifacts=source_artifacts,
    )
    provenance = make_provenance(
        config=config,
        config_hash=config_hash,
        source_hashes=source_hashes,
        hierarchy_label=hierarchy_label,
        n_side=n_side,
        seed=seed,
        control_group=control_group_label,
    )
    observation_id = stable_hash(
        {
            "run_instance_id": provenance.run_instance_id,
            "branch_signature_id": signature.branch_signature_id,
            "observation_kind": observation_kind,
        },
        prefix="obs_",
    )
    return RecoveryObservation(
        observation_id=observation_id,
        observation_kind=observation_kind,
        control_role=control_role,
        control_group_label=control_group_label,
        hierarchy_label=hierarchy_label,
        candidate_label=config["candidate_label"],
        n_side=int(n_side),
        ensemble_size=int(config["ensemble_size"]),
        seed=int(seed),
        probe=config["probe"],
        perturbation=PerturbationSpec(
            name=config["probe"],
            family="frozen_phase18_probe",
            parameters={"n_side": int(n_side), "seed": int(seed), "ensemble_size": int(config["ensemble_size"])},
            source="phase18_shell_ordering_metrics",
        ),
        provenance=provenance,
        closure_signature=signature,
        notes=["repository-native CST toy observation; no external physics claim"],
    )


def summary_stats(values: list[float]) -> dict[str, float | int | None]:
    if not values:
        return {
            "count": 0,
            "mean": None,
            "std": None,
            "min": None,
            "max": None,
        }
    return {
        "count": len(values),
        "mean": mean(values),
        "std": population_std(values),
        "min": min(values),
        "max": max(values),
    }


def distance_row(distance: ClosureDistanceComponents, observations_by_id: dict[str, RecoveryObservation]) -> dict[str, Any]:
    left = observations_by_id[distance.left_observation_id]
    right = observations_by_id[distance.right_observation_id]
    return {
        "left_observation_id": distance.left_observation_id,
        "right_observation_id": distance.right_observation_id,
        "left_kind": left.observation_kind,
        "right_kind": right.observation_kind,
        "left_control_group": left.control_group_label,
        "right_control_group": right.control_group_label,
        "left_signature_id": distance.left_signature_id,
        "right_signature_id": distance.right_signature_id,
        "scalar_distance": distance.scalar_distance,
        "scalar_enabled": distance.scalar_enabled,
        "disabled_components": "|".join(distance.disabled_components),
        "warnings": "|".join(distance.warnings),
    }


def component_rows(distance: ClosureDistanceComponents) -> list[dict[str, Any]]:
    rows = []
    for name in sorted(distance.components):
        component = distance.components[name]
        rows.append(
            {
                "left_observation_id": distance.left_observation_id,
                "right_observation_id": distance.right_observation_id,
                "component": name,
                "value": component.value,
                "availability": component.availability,
                "normalization_method": component.normalization_method,
                "source_fields": "|".join(component.source_fields),
                "warnings": "|".join(component.warnings),
            }
        )
    return rows


def render_cst_report(
    *,
    result: CSTBenchmarkResult,
    target_distances: list[ClosureDistanceComponents],
    control_distances: list[ClosureDistanceComponents],
    unavailable_controls: dict[str, Any],
    output_paths: dict[str, str],
) -> str:
    target_scalars = scalar_values(target_distances)
    control_scalars = scalar_values(control_distances)
    vector = result.recoverability_vector
    lines = [
        "# CST Toy Benchmark Report",
        "",
        "Status: " + result.verdict,
        "",
        "Implemented fact: this report is generated from frozen HAOS-IIP ledgers only.",
        "Design choice: branch identity is an operational hash of the compact address summary.",
        "Heuristic: non-front event times are rank-interpolated from Phase 18 recovery-time fields.",
        "Analogy: generalized CST selectivity is an engineering signal/noise ratio, not physical resonance Q.",
        "Unverified hypothesis: repeated recovery of this address might indicate closure stability in this toy slice.",
        "",
        "## Reasons",
    ]
    for reason in result.reasons:
        lines.append(f"- {reason}")
    if result.warnings:
        lines.append("")
        lines.append("## Warnings")
        for warning in result.warnings:
            lines.append(f"- {warning}")
    lines.extend(
        [
            "",
            "## Recoverability Vector",
            f"- recovery_rate: {vector.recovery_rate:.12g}",
            f"- branch_identity_persistence: {vector.branch_identity_persistence:.12g}",
            f"- closure_fidelity: {vector.closure_fidelity:.12g}",
            f"- control_margin: {vector.control_margin:.12g}",
            f"- variance_bound: {vector.variance_bound:.12g}",
            f"- ablation_balance: {vector.ablation_balance:.12g}",
            f"- optional_scalar: {vector.optional_scalar}",
            "",
            "## Distance Summary",
            f"- target_scalar_mean: {mean(target_scalars):.12g}",
            f"- target_scalar_std: {population_std(target_scalars):.12g}",
            f"- control_scalar_mean: {mean(control_scalars):.12g}",
            f"- control_scalar_std: {population_std(control_scalars):.12g}",
            "",
            "## Unavailable Controls",
        ]
    )
    if unavailable_controls:
        for label, payload in sorted(unavailable_controls.items()):
            lines.append(f"- {label}: {payload.get('reason')}")
    else:
        lines.append("- none")
    lines.extend(["", "## Output Files"])
    for label, path in sorted(output_paths.items()):
        lines.append(f"- {label}: {path}")
    lines.extend(
        [
            "",
            "## Limitations",
            "- This benchmark does not test Bell correlations, spectra, semiconductors, biology, consciousness, or cosmology.",
            "- A PASS/OPEN/FAIL verdict is scoped only to this frozen toy slice and configured thresholds.",
            "- Scalar compression is optional and heuristic; component metrics and the vector are the auditable result.",
        ]
    )
    return "\n".join(lines) + "\n"


def run_cst_benchmark(
    *,
    config_path: Path = DEFAULT_CONFIG_PATH,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    write_outputs: bool = True,
    scalar_enabled_override: bool | None = None,
) -> dict[str, Any]:
    config = read_json(config_path)
    if scalar_enabled_override is not None:
        config = dict(config)
        config["scalar_compression"] = dict(config["scalar_compression"])
        config["scalar_compression"]["enabled"] = bool(scalar_enabled_override)
    config_hash = stable_hash(config, prefix="cfg_")
    source_hashes = artifact_hashes(PHASE_ARTIFACTS)
    ledgers = {name: read_csv(path) for name, path in PHASE_ARTIFACTS.items()}
    indexes = index_ledgers(ledgers)
    source_artifacts = [repo_rel(path) for _, path in sorted(PHASE_ARTIFACTS.items())]

    target_label = config["hierarchies"]["target"]
    periodic_label = config["hierarchies"]["control"]
    observations: list[RecoveryObservation] = []
    target_observations: dict[tuple[int, int], RecoveryObservation] = {}
    periodic_observations: dict[tuple[int, int], RecoveryObservation] = {}

    for n_side in config["n_sides"]:
        for seed in config["seeds"]:
            target = make_observation(
                hierarchy_label=target_label,
                n_side=int(n_side),
                seed=int(seed),
                observation_kind="target",
                control_role="target",
                control_group_label="target",
                config=config,
                indexes=indexes,
                config_hash=config_hash,
                source_hashes=source_hashes,
                source_artifacts=source_artifacts,
            )
            target_observations[(int(n_side), int(seed))] = target
            observations.append(target)

            periodic = make_observation(
                hierarchy_label=periodic_label,
                n_side=int(n_side),
                seed=int(seed),
                observation_kind="frozen_control",
                control_role="negative",
                control_group_label=periodic_label,
                config=config,
                indexes=indexes,
                config_hash=config_hash,
                source_hashes=source_hashes,
                source_artifacts=source_artifacts,
            )
            periodic_observations[(int(n_side), int(seed))] = periodic
            observations.append(periodic)

    controls_by_source: dict[str, list[RecoveryObservation]] = {}
    unavailable_controls: dict[str, Any] = {}
    signature_version = config["signature"]["signature_version"]
    float_precision = int(config["signature"]["float_precision"])
    for control in config["controls"]:
        label = control["label"]
        role = control["role"]
        source = control["source"]
        if role == "unavailable" or source == "not_present_in_frozen_slice":
            unavailable_controls[label] = {
                "status": "Missing",
                "reason": "not present in frozen Phase 15-18 artifact slice",
            }
            continue
        if source != "derived_signature_transform":
            continue
        for target in target_observations.values():
            derived = build_derived_control(
                target,
                control_label=label,
                control_role=role,
                config_hash=config_hash,
                signature_version=signature_version,
                float_precision=float_precision,
            )
            controls_by_source.setdefault(target.observation_id, []).append(derived)
            observations.append(derived)

    target_pairs: list[tuple[RecoveryObservation, RecoveryObservation]] = []
    for n_side in config["n_sides"]:
        seeds = [int(seed) for seed in config["seeds"]]
        for left_seed, right_seed in itertools.combinations(seeds, 2):
            target_pairs.append((target_observations[(int(n_side), left_seed)], target_observations[(int(n_side), right_seed)]))

    negative_pairs: list[tuple[RecoveryObservation, RecoveryObservation]] = []
    invariance_pairs: list[tuple[RecoveryObservation, RecoveryObservation]] = []
    for key, target in target_observations.items():
        periodic = periodic_observations.get(key)
        if periodic is None:
            unavailable_controls[periodic_label] = {
                "status": "Missing",
                "reason": "matched periodic control row missing",
            }
        else:
            negative_pairs.append((target, periodic))
        for control in controls_by_source.get(target.observation_id, []):
            if control.control_group_label in NEGATIVE_CONTROL_LABELS:
                negative_pairs.append((target, control))
            else:
                invariance_pairs.append((target, control))

    weights = config["distance"]["weights"]
    enabled_components = config["distance"]["enabled_components"]
    scalar_enabled = bool(config["distance"]["scalar_enabled"])
    target_distances = [
        compute_closure_distance(
            left,
            right,
            weights=weights,
            enabled_components=enabled_components,
            scalar_enabled=scalar_enabled,
        )
        for left, right in target_pairs
    ]
    control_distances = [
        compute_closure_distance(
            left,
            right,
            weights=weights,
            enabled_components=enabled_components,
            scalar_enabled=scalar_enabled,
        )
        for left, right in negative_pairs
    ]
    invariance_distances = [
        compute_closure_distance(
            left,
            right,
            weights=weights,
            enabled_components=enabled_components,
            scalar_enabled=scalar_enabled,
        )
        for left, right in invariance_pairs
    ]
    target_signature_matches = [
        left.closure_signature.branch_signature_id == right.closure_signature.branch_signature_id
        for left, right in target_pairs
    ]
    negative_control_family_stats: dict[str, Any] = {}
    for label in sorted({right.control_group_label for _, right in negative_pairs}):
        relevant = [
            distance
            for distance, (_, right) in zip(control_distances, negative_pairs)
            if right.control_group_label == label
        ]
        stats = summary_stats(scalar_values(relevant))
        negative_control_family_stats[label] = {
            "availability": "available",
            "role": next(right.control_role for _, right in negative_pairs if right.control_group_label == label),
            "count": stats["count"],
            "mean_distance": stats["mean"],
            "std_distance": stats["std"],
            "min_distance": stats["min"],
            "max_distance": stats["max"],
        }

    weight_sensitivity = compute_weight_sensitivity(
        target_pairs,
        negative_pairs,
        weight_families=config["distance"]["weight_families"],
        enabled_components=enabled_components,
    )
    ablation_score, ablation_details = ablation_balance(
        target_pairs,
        negative_pairs,
        base_weights=weights,
        enabled_components=enabled_components,
    )
    target_scalars = scalar_values(target_distances)
    control_scalars = scalar_values(control_distances)
    selectivity = None
    if config["selectivity"]["enabled"]:
        provisional_recovery_rate = mean([
            1.0 if value <= float(config["thresholds"]["closure_distance_max"]) else 0.0
            for value in target_scalars
        ])
        provisional_identity = mean([1.0 if item else 0.0 for item in target_signature_matches])
        provisional_fidelity = max(0.0, 1.0 - mean(target_scalars))
        selectivity = selectivity_payload(
            recovery_rate=provisional_recovery_rate,
            branch_identity_persistence=provisional_identity,
            closure_fidelity=provisional_fidelity,
            control_mean_distance=mean(control_scalars),
            target_distances=target_scalars,
            epsilons=[float(value) for value in config["selectivity"]["epsilon_sensitivity"]],
            default_epsilon=float(config["selectivity"]["default_epsilon"]),
        )

    vector = recoverability_vector(
        target_distances=target_distances,
        target_signature_matches=target_signature_matches,
        control_distances=control_distances,
        thresholds=config["thresholds"],
        ablation_score=ablation_score,
        selectivity=selectivity,
        scalar_enabled=scalar_enabled and bool(config["scalar_compression"]["enabled"]),
    )
    verdict, reasons, warnings = verdict_from_metrics(
        vector=vector,
        target_distances=target_distances,
        control_distances=control_distances,
        target_seed_count=len(config["seeds"]),
        thresholds=config["thresholds"],
        weight_sensitivity=weight_sensitivity,
        control_family_stats=negative_control_family_stats,
        unavailable_controls=unavailable_controls,
        forced_fail_probe_passed=forced_fail_probe(),
    )
    weight_sensitivity_payload = {**weight_sensitivity, "ablation_details": ablation_details}
    result_payload = {
        "verdict": verdict,
        "reasons": reasons,
        "warnings": warnings,
        "recoverability_vector": vector.to_dict(),
        "thresholds": config["thresholds"],
        "weight_sensitivity": weight_sensitivity_payload,
        "ablation_details": ablation_details,
        "control_family_stats": negative_control_family_stats,
        "unavailable_controls": unavailable_controls,
        "artifact_hashes": source_hashes,
        "config_hash": config_hash,
    }
    result_hash = stable_hash(result_payload, prefix="cst_result_")
    result = CSTBenchmarkResult(
        verdict=verdict,
        reasons=reasons,
        warnings=warnings,
        recoverability_vector=vector,
        thresholds=config["thresholds"],
        weight_sensitivity=weight_sensitivity_payload,
        unavailable_controls=unavailable_controls,
        artifact_hashes=source_hashes,
        config_hash=config_hash,
        result_hash=result_hash,
    )

    observations_by_id = {observation.observation_id: observation for observation in observations}
    all_distances: list[ClosureDistanceComponents] = []
    for left, right in itertools.product(observations, observations):
        all_distances.append(
            compute_closure_distance(
                left,
                right,
                weights=weights,
                enabled_components=enabled_components,
                scalar_enabled=scalar_enabled,
            )
        )

    control_distribution_rows: list[dict[str, Any]] = []
    for label, stats in negative_control_family_stats.items():
        control_distribution_rows.append(
            {
                "control_label": label,
                "availability": stats["availability"],
                "role": stats["role"],
                "count": stats["count"],
                "mean_distance": stats["mean_distance"],
                "std_distance": stats["std_distance"],
                "min_distance": stats["min_distance"],
                "max_distance": stats["max_distance"],
                "reason": "",
            }
        )
    for label in sorted({right.control_group_label for _, right in invariance_pairs}):
        relevant = [
            distance
            for distance, (_, right) in zip(invariance_distances, invariance_pairs)
            if right.control_group_label == label
        ]
        stats = summary_stats(scalar_values(relevant))
        control_distribution_rows.append(
            {
                "control_label": label,
                "availability": "available",
                "role": next(right.control_role for _, right in negative_pairs + invariance_pairs if right.control_group_label == label),
                "count": stats["count"],
                "mean_distance": stats["mean"],
                "std_distance": stats["std"],
                "min_distance": stats["min"],
                "max_distance": stats["max"],
                "reason": "",
            }
        )
    for label, payload in sorted(unavailable_controls.items()):
        control_distribution_rows.append(
            {
                "control_label": label,
                "availability": "unavailable",
                "role": "unavailable",
                "count": 0,
                "mean_distance": None,
                "std_distance": None,
                "min_distance": None,
                "max_distance": None,
                "reason": payload.get("reason", ""),
            }
        )

    output_paths = {
        "cst_runs": repo_rel(output_dir / "cst_runs.json"),
        "closure_signatures": repo_rel(output_dir / "closure_signatures.json"),
        "closure_distance_matrix": repo_rel(output_dir / "closure_distance_matrix.csv"),
        "closure_distance_components": repo_rel(output_dir / "closure_distance_components.csv"),
        "recoverability_vectors": repo_rel(output_dir / "recoverability_vectors.csv"),
        "control_distributions": repo_rel(output_dir / "control_distributions.csv"),
        "seed_manifest": repo_rel(output_dir / "seed_manifest.json"),
        "benchmark_result": repo_rel(output_dir / "benchmark_result.json"),
        "benchmark_report": repo_rel(output_dir / "benchmark_report.md"),
    }
    benchmark_payload = result.to_dict()
    benchmark_payload["ablation_details"] = ablation_details
    benchmark_payload["control_family_stats"] = negative_control_family_stats
    benchmark_payload["result_hash"] = result_hash

    if write_outputs:
        write_json(
            output_dir / "cst_runs.json",
            {
                "benchmark_label": config["benchmark_label"],
                "config_hash": config_hash,
                "observations": [observation.to_dict() for observation in observations],
            },
        )
        write_json(
            output_dir / "closure_signatures.json",
            {
                "benchmark_label": config["benchmark_label"],
                "signatures": [
                    {
                        "observation_id": observation.observation_id,
                        **observation.closure_signature.to_dict(),
                    }
                    for observation in observations
                ],
            },
        )
        write_csv(
            output_dir / "closure_distance_matrix.csv",
            [distance_row(distance, observations_by_id) for distance in all_distances],
            [
                "left_observation_id",
                "right_observation_id",
                "left_kind",
                "right_kind",
                "left_control_group",
                "right_control_group",
                "left_signature_id",
                "right_signature_id",
                "scalar_distance",
                "scalar_enabled",
                "disabled_components",
                "warnings",
            ],
        )
        write_csv(
            output_dir / "closure_distance_components.csv",
            [row for distance in all_distances for row in component_rows(distance)],
            [
                "left_observation_id",
                "right_observation_id",
                "component",
                "value",
                "availability",
                "normalization_method",
                "source_fields",
                "warnings",
            ],
        )
        write_csv(
            output_dir / "recoverability_vectors.csv",
            [
                {
                    "verdict": verdict,
                    **vector.to_dict(),
                    "target_distance_mean": mean(target_scalars),
                    "target_distance_std": population_std(target_scalars),
                    "control_distance_mean": mean(control_scalars),
                    "control_distance_std": population_std(control_scalars),
                }
            ],
            [
                "verdict",
                "recovery_rate",
                "branch_identity_persistence",
                "closure_fidelity",
                "control_margin",
                "variance_bound",
                "ablation_balance",
                "authoritative_mode",
                "optional_scalar",
                "scalar_formula",
                "scalar_warnings",
                "generalized_cst_selectivity",
                "target_distance_mean",
                "target_distance_std",
                "control_distance_mean",
                "control_distance_std",
            ],
        )
        write_csv(
            output_dir / "control_distributions.csv",
            control_distribution_rows,
            [
                "control_label",
                "availability",
                "role",
                "count",
                "mean_distance",
                "std_distance",
                "min_distance",
                "max_distance",
                "reason",
            ],
        )
        write_json(
            output_dir / "seed_manifest.json",
            {
                "benchmark_label": config["benchmark_label"],
                "config_hash": config_hash,
                "seeds": config["seeds"],
                "n_sides": config["n_sides"],
                "target_run_instance_ids": {
                    f"n{n_side}_seed{seed}": observation.provenance.run_instance_id
                    for (n_side, seed), observation in sorted(target_observations.items())
                },
                "source_artifacts": source_artifact_paths(),
                "source_artifact_hashes": source_hashes,
            },
        )
        write_json(output_dir / "benchmark_result.json", benchmark_payload)
        (output_dir / "benchmark_report.md").write_text(
            render_cst_report(
                result=result,
                target_distances=target_distances,
                control_distances=control_distances,
                unavailable_controls=unavailable_controls,
                output_paths=output_paths,
            ),
            encoding="utf-8",
        )

    return {
        "config": config,
        "config_hash": config_hash,
        "observations": observations,
        "target_distances": target_distances,
        "control_distances": control_distances,
        "invariance_distances": invariance_distances,
        "all_distances": all_distances,
        "result": result,
        "benchmark_payload": benchmark_payload,
        "control_distributions": control_distribution_rows,
        "output_paths": output_paths,
        "canonical_payload": canonical_json(benchmark_payload),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the experimental HAOS-IIP CST toy benchmark.")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--no-scalar", action="store_true", help="disable optional scalar compression")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = run_cst_benchmark(
        config_path=args.config,
        output_dir=args.output_dir,
        scalar_enabled_override=False if args.no_scalar else None,
        write_outputs=True,
    )
    benchmark = result["result"]
    print(f"CST verdict: {benchmark.verdict}")
    print(f"result_hash: {benchmark.result_hash}")
    print(f"report: {result['output_paths']['benchmark_report']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
