from __future__ import annotations

import copy
from dataclasses import replace
from typing import Any

from experiments.cst.cst_io import stable_hash
from experiments.cst.cst_metrics import branch_signature_id
from experiments.cst.cst_types import (
    CSTRunProvenance,
    ClosureSignature,
    PerturbationSpec,
    RecoveryObservation,
)


NEGATIVE_CONTROL_LABELS = {
    "periodic_diagonal_augmented_control",
    "shuffled_structure_control",
    "randomized_edge_control",
    "degraded_signal_control",
    "label_permutation_control",
    "generic_graph_operator_control",
    "parameter_matched_null_control",
}


def chain_labels(signature_fields: dict[str, Any]) -> list[str]:
    chain = signature_fields.get("event_chain_summary", {}).get("chain_signature", "")
    return str(chain).split(">") if chain else []


def edges_from_chain(labels: list[str]) -> list[str]:
    return [f"{labels[index]}->{labels[index + 1]}" for index in range(len(labels) - 1)]


def update_chain_fields(fields: dict[str, Any], labels: list[str]) -> None:
    chain = ">".join(labels)
    fields.setdefault("event_chain_summary", {})["chain_signature"] = chain
    fields["event_chain_summary"]["dominant_chain_prefix"] = labels[:2]
    fields.setdefault("influence_edge_summary", {})["edges"] = edges_from_chain(labels)
    fields["influence_edge_summary"]["edge_count"] = max(
        int(fields["influence_edge_summary"].get("edge_count", len(labels) - 1)),
        len(labels) - 1,
    )
    fields["influence_edge_summary"]["stable_edge_count"] = min(
        int(fields["influence_edge_summary"].get("stable_edge_count", max(len(labels) - 2, 0))),
        fields["influence_edge_summary"]["edge_count"],
    )
    fields.setdefault("address_summary", {})["dominant_chain_prefix"] = labels[:2]


def recompute_signature(
    fields: dict[str, Any],
    *,
    signature_version: str,
    float_precision: int,
    source_artifacts: list[str],
) -> ClosureSignature:
    return ClosureSignature(
        branch_signature_id=branch_signature_id(fields, float_precision=float_precision),
        signature_version=signature_version,
        operational_equivalence_note=(
            "Experimental CST branch identity: stable hash of the operational address summary. "
            "It is a reproducibility convention, not an ontological claim."
        ),
        fields=fields,
        unavailable_components={},
        source_artifacts=list(source_artifacts),
        warnings=list(fields.get("warnings", [])),
    )


def control_provenance(
    source: RecoveryObservation,
    *,
    control_label: str,
    config_hash: str,
) -> CSTRunProvenance:
    payload = {
        "source_run_instance_id": source.provenance.run_instance_id,
        "source_observation_id": source.observation_id,
        "control_label": control_label,
        "config_hash": config_hash,
    }
    run_instance_id = stable_hash(payload, prefix="run_")
    return replace(
        source.provenance,
        run_instance_id=run_instance_id,
        hierarchy_label=source.hierarchy_label,
        control_group=control_label,
    )


def clone_control_observation(
    source: RecoveryObservation,
    *,
    control_label: str,
    control_role: str,
    fields: dict[str, Any],
    config_hash: str,
    signature_version: str,
    float_precision: int,
    source_artifacts: list[str],
    notes: list[str],
) -> RecoveryObservation:
    signature = recompute_signature(
        fields,
        signature_version=signature_version,
        float_precision=float_precision,
        source_artifacts=source_artifacts,
    )
    provenance = control_provenance(source, control_label=control_label, config_hash=config_hash)
    observation_id = stable_hash(
        {
            "run_instance_id": provenance.run_instance_id,
            "branch_signature_id": signature.branch_signature_id,
            "control_label": control_label,
        },
        prefix="obs_",
    )
    return RecoveryObservation(
        observation_id=observation_id,
        observation_kind="derived_control",
        control_role=control_role,
        control_group_label=control_label,
        hierarchy_label=source.hierarchy_label,
        candidate_label=source.candidate_label,
        n_side=source.n_side,
        ensemble_size=source.ensemble_size,
        seed=source.seed,
        probe=source.probe,
        perturbation=PerturbationSpec(
            name=control_label,
            family="cst_matched_control",
            parameters={"source_observation_id": source.observation_id},
            source="derived_signature_transform",
        ),
        provenance=provenance,
        closure_signature=signature,
        notes=notes,
    )


def label_permutation(fields: dict[str, Any]) -> dict[str, Any]:
    labels = chain_labels(fields)
    mapping = {label: f"perm_{index:02d}" for index, label in enumerate(sorted(set(labels)))}
    permuted = [mapping[label] for label in labels]
    fields = copy.deepcopy(fields)
    update_chain_fields(fields, permuted)
    fields.setdefault("control_transform", {})["label_permutation"] = mapping
    return fields


def shuffled_structure(fields: dict[str, Any]) -> dict[str, Any]:
    labels = sorted(chain_labels(fields))
    fields = copy.deepcopy(fields)
    update_chain_fields(fields, labels)
    fields.setdefault("control_transform", {})["shuffle_rule"] = "lexicographic event-label sort"
    return fields


def randomized_edge(fields: dict[str, Any]) -> dict[str, Any]:
    labels = chain_labels(fields)
    if labels:
        keyed = sorted(labels, key=lambda label: stable_hash({"label": label, "salt": "cst_edge"}, prefix=""))
    else:
        keyed = labels
    fields = copy.deepcopy(fields)
    update_chain_fields(fields, keyed)
    fields.setdefault("control_transform", {})["randomization_rule"] = "deterministic stable-hash event order"
    return fields


def degraded_signal(fields: dict[str, Any]) -> dict[str, Any]:
    labels = list(reversed(chain_labels(fields)))
    fields = copy.deepcopy(fields)
    update_chain_fields(fields, labels)
    for section_name in ("invariant_summary", "spectral_low_mode_summary"):
        section = fields.setdefault(section_name, {})
        for key in ("terminal_low_k_fraction", "terminal_low_k_power"):
            if isinstance(section.get(key), (int, float)):
                section[key] = max(0.0, float(section[key]) * 0.5)
    causal = fields.setdefault("causal_summary", {})
    if isinstance(causal.get("acyclicity_score"), (int, float)):
        causal["acyclicity_score"] = max(0.0, float(causal["acyclicity_score"]) * 0.5)
    if isinstance(causal.get("mismatch_rate"), (int, float)):
        causal["mismatch_rate"] = min(1.0, float(causal["mismatch_rate"]) * 4.0 + 0.1)
    fields.setdefault("shell_ordering_summary", {})["monotonic_shell_ordering"] = False
    fields.setdefault("address_summary", {})["monotonic_shell_ordering"] = False
    fields.setdefault("control_transform", {})["degradation_rule"] = "reverse chain and attenuate low-mode descriptors"
    return fields


def generic_graph_operator(fields: dict[str, Any]) -> dict[str, Any]:
    labels = ["generic_source", "generic_mid", "generic_sink"]
    fields = copy.deepcopy(fields)
    update_chain_fields(fields, labels)
    fields["influence_edge_summary"]["edge_count"] = 2
    fields["influence_edge_summary"]["stable_edge_count"] = 0
    fields["causal_summary"]["front_depth_profile"] = {"front_near_depth": 1.0, "front_far_depth": 1.0}
    fields["address_summary"]["front_depth_profile"] = {"front_near_depth": 1.0, "front_far_depth": 1.0}
    fields.setdefault("control_transform", {})["generic_rule"] = "three-node generic operator scaffold"
    return fields


def parameter_matched_null(fields: dict[str, Any]) -> dict[str, Any]:
    fields = copy.deepcopy(fields)
    update_chain_fields(fields, ["null_source", "null_response"])
    for section_name in ("invariant_summary", "spectral_low_mode_summary"):
        section = fields.setdefault(section_name, {})
        for key, value in list(section.items()):
            if isinstance(value, (int, float)):
                section[key] = 0.0
    causal = fields.setdefault("causal_summary", {})
    causal["acyclicity_score"] = 0.0
    causal["mean_causal_depth"] = 0.0
    causal["mismatch_rate"] = 1.0
    fields.setdefault("control_transform", {})["null_rule"] = "same n_side/seed/probe metadata with zeroed descriptor payload"
    return fields


def seed_repeat(fields: dict[str, Any]) -> dict[str, Any]:
    fields = copy.deepcopy(fields)
    fields.setdefault("control_transform", {})["seed_repeat_rule"] = "exact signature replay with independent run_instance_id"
    return fields


CONTROL_TRANSFORMS = {
    "shuffled_structure_control": shuffled_structure,
    "randomized_edge_control": randomized_edge,
    "degraded_signal_control": degraded_signal,
    "label_permutation_control": label_permutation,
    "generic_graph_operator_control": generic_graph_operator,
    "parameter_matched_null_control": parameter_matched_null,
    "seed_repeat_control": seed_repeat,
}


def build_derived_control(
    source: RecoveryObservation,
    *,
    control_label: str,
    control_role: str,
    config_hash: str,
    signature_version: str,
    float_precision: int,
) -> RecoveryObservation:
    if control_label not in CONTROL_TRANSFORMS:
        raise KeyError(f"unsupported CST derived control: {control_label}")
    transform = CONTROL_TRANSFORMS[control_label]
    fields = transform(source.closure_signature.fields)
    return clone_control_observation(
        source,
        control_label=control_label,
        control_role=control_role,
        fields=fields,
        config_hash=config_hash,
        signature_version=signature_version,
        float_precision=float_precision,
        source_artifacts=source.closure_signature.source_artifacts,
        notes=[f"derived by {control_label} from {source.observation_id}"],
    )
