from __future__ import annotations

import hashlib
import json
from typing import Any

from experiments.hbp.hbp_types import (
    BRIDGE_LEVELS,
    MAPPING_STATUSES,
    BridgeAssessment,
    BridgeContract,
    BridgeProvenance,
    BaselineSpec,
    ControlSpec,
    DomainStateSpec,
    DynamicsSpec,
    FalsificationSpec,
    HAOSMapping,
    ObservationMap,
)


PREDICTIVE_MAPPING_STATUSES = {"DERIVED", "CALIBRATED"}


def stable_hash(payload: Any, prefix: str) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return f"{prefix}{hashlib.sha256(encoded).hexdigest()[:24]}"


def level_index(level: str) -> int:
    if level not in BRIDGE_LEVELS:
        raise ValueError(f"unknown bridge level: {level}")
    return BRIDGE_LEVELS.index(level)


def min_level(left: str, right: str) -> str:
    return BRIDGE_LEVELS[min(level_index(left), level_index(right))]


def _as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def contract_from_dict(payload: dict[str, Any]) -> BridgeContract:
    state_payload = payload.get("state_schema")
    dynamics_payload = payload.get("dynamics")
    falsification_payload = payload.get("falsification")
    provenance_payload = payload.get("provenance") or {}
    return BridgeContract(
        bridge_id=str(payload["bridge_id"]),
        domain=str(payload.get("domain", "")),
        classification=str(payload.get("classification", "FORMAL_CORRESPONDENCE")),
        state_schema=DomainStateSpec(**state_payload) if isinstance(state_payload, dict) else None,
        units=dict(payload.get("units") or {}),
        dynamics=DynamicsSpec(**dynamics_payload) if isinstance(dynamics_payload, dict) else None,
        boundary_conditions=str(payload.get("boundary_conditions", "")),
        symmetries=[str(item) for item in _as_list(payload.get("symmetries"))],
        haos_mapping=[HAOSMapping(**item) for item in _as_list(payload.get("haos_mapping"))],
        observation_map=[ObservationMap(**item) for item in _as_list(payload.get("observation_map"))],
        prediction_target=str(payload.get("prediction_target", "")),
        calibration_split=str(payload.get("calibration_split", "")),
        holdout_split=str(payload.get("holdout_split", "")),
        controls=[ControlSpec(**item) for item in _as_list(payload.get("controls"))],
        baselines=[BaselineSpec(**item) for item in _as_list(payload.get("baselines"))],
        uncertainty=str(payload.get("uncertainty", "")),
        falsification=FalsificationSpec(**falsification_payload) if isinstance(falsification_payload, dict) else None,
        verdict_logic=str(payload.get("verdict_logic", "")),
        provenance=BridgeProvenance(
            source_artifacts=[str(item) for item in _as_list(provenance_payload.get("source_artifacts"))],
            code_artifacts=[str(item) for item in _as_list(provenance_payload.get("code_artifacts"))],
            precommitment_artifact=provenance_payload.get("precommitment_artifact"),
            external_data_status=str(provenance_payload.get("external_data_status", "unknown")),
            level_history=[str(item) for item in _as_list(provenance_payload.get("level_history"))],
            notes=[str(item) for item in _as_list(provenance_payload.get("notes"))],
        ),
    )


def missing_required_fields(contract: BridgeContract) -> list[str]:
    missing: list[str] = []
    if contract.state_schema is None or not contract.state_schema.variables:
        missing.append("state_schema")
    if not contract.units:
        missing.append("units")
    if contract.dynamics is None or not contract.dynamics.update_law:
        missing.append("dynamics")
    if not contract.observation_map:
        missing.append("observation_map")
    if contract.falsification is None or not contract.falsification.fail_conditions:
        missing.append("falsification")
    return missing


def unit_validation_errors(contract: BridgeContract) -> list[str]:
    errors: list[str] = []
    if contract.state_schema:
        for variable in contract.state_schema.variables:
            if variable not in contract.state_schema.units:
                errors.append(f"missing state unit for {variable}")
            if variable not in contract.state_schema.valid_ranges:
                errors.append(f"missing valid range for {variable}")
    for mapping in contract.haos_mapping:
        if not mapping.units_before or not mapping.units_after:
            errors.append(f"mapping {mapping.domain_target} missing units")
    for observation in contract.observation_map:
        if not observation.units or not observation.valid_range:
            errors.append(f"observation {observation.observable} missing units or valid range")
    return errors


def classification_ceiling(contract: BridgeContract) -> tuple[str, list[str], list[str]]:
    reasons: list[str] = []
    labels: list[str] = []
    ceiling = "PHYSICAL_MECHANISM_CANDIDATE"

    missing = missing_required_fields(contract)
    if missing:
        ceiling = min_level(ceiling, "FORMAL_CORRESPONDENCE")
        labels.append("CONTRACT_INCOMPLETE")
        reasons.append(f"missing required bridge fields: {', '.join(missing)}")

    mapping_statuses = {mapping.domain_target: mapping.mapping_status for mapping in contract.haos_mapping}
    invalid_statuses = sorted({status for status in mapping_statuses.values() if status not in MAPPING_STATUSES})
    if invalid_statuses:
        ceiling = min_level(ceiling, "FORMAL_CORRESPONDENCE")
        labels.append("MAPPING_INVALID")
        reasons.append(f"invalid mapping statuses: {', '.join(invalid_statuses)}")

    if any(status in {"HEURISTIC", "ANALOGICAL"} for status in mapping_statuses.values()):
        ceiling = min_level(ceiling, "OPERATIONAL_MAPPING")
        labels.append("OPERATIONAL_MAPPING_PARTIAL")
        reasons.append("heuristic or analogical mapping cannot support a higher bridge classification")

    if not any(status in PREDICTIVE_MAPPING_STATUSES for status in mapping_statuses.values()):
        ceiling = min_level(ceiling, "OPERATIONAL_MAPPING")
        labels.append("OPERATIONAL_MAPPING_PARTIAL")
        reasons.append("no derived or independently frozen calibrated mapping is available")

    if not contract.holdout_split or contract.holdout_split.lower() in {"none", "not_applicable", "unavailable"}:
        ceiling = min_level(ceiling, "OPERATIONAL_MAPPING")
        reasons.append("no frozen unseen holdout prediction is declared")

    if (
        not contract.provenance.external_data_status.startswith("external_raw_data")
        or not contract.controls
        or not contract.uncertainty
        or "replication" not in contract.verdict_logic.lower()
    ):
        ceiling = min_level(ceiling, "PREDICTIVE_BRIDGE")
        labels.append("PHYSICAL_MECHANISM_NOT_ESTABLISHED")
        reasons.append("external raw data, controls, uncertainty, or replication are insufficient for empirical support")

    if "intervention" not in contract.verdict_logic.lower() or "unique generative mechanism" not in contract.verdict_logic.lower():
        ceiling = min_level(ceiling, "EMPIRICALLY_SUPPORTED_BRIDGE")
        labels.append("PHYSICAL_MECHANISM_NOT_ESTABLISHED")
        reasons.append("no intervention-backed unique generative mechanism is declared")

    if unit_validation_errors(contract):
        ceiling = min_level(ceiling, "FORMAL_CORRESPONDENCE")
        labels.append("DIMENSIONAL_ANALYSIS_FAIL")
        reasons.extend(unit_validation_errors(contract))

    return ceiling, labels, reasons


def level_skip_errors(contract: BridgeContract) -> list[str]:
    history = contract.provenance.level_history
    if not history:
        return []
    errors: list[str] = []
    for item in history:
        if item not in BRIDGE_LEVELS:
            errors.append(f"unknown level in history: {item}")
    if errors:
        return errors
    requested = level_index(contract.classification)
    highest_history = max(level_index(item) for item in history)
    if requested > highest_history + 1:
        errors.append("requested classification skips a bridge level")
    return errors


def assess_contract(contract: BridgeContract) -> BridgeAssessment:
    ceiling, labels, reasons = classification_ceiling(contract)
    missing = missing_required_fields(contract)
    skip_errors = level_skip_errors(contract)
    if skip_errors:
        ceiling = min_level(ceiling, "FORMAL_CORRESPONDENCE")
        labels.append("CONTRACT_INCOMPLETE")
        reasons.extend(skip_errors)

    requested = contract.classification
    effective = min_level(requested, ceiling)
    if effective != requested:
        reasons.append(f"classification demoted from {requested} to {effective}")

    mapping_statuses = {mapping.domain_target: mapping.mapping_status for mapping in contract.haos_mapping}
    verdict_authorized = (
        effective in {"PREDICTIVE_BRIDGE", "EMPIRICALLY_SUPPORTED_BRIDGE", "PHYSICAL_MECHANISM_CANDIDATE"}
        and any(status in PREDICTIVE_MAPPING_STATUSES for status in mapping_statuses.values())
        and not missing
        and not skip_errors
    )

    if effective == "FORMAL_CORRESPONDENCE":
        labels.append("FORMAL_CORRESPONDENCE_ONLY")
    elif effective == "OPERATIONAL_MAPPING":
        labels.append("OPERATIONAL_MAPPING_VALID" if not missing else "OPERATIONAL_MAPPING_PARTIAL")

    labels = sorted(set(labels))
    assessment_payload = {
        "bridge_id": contract.bridge_id,
        "requested_classification": requested,
        "effective_classification": effective,
        "verdict_authorized": verdict_authorized,
        "labels": labels,
        "reasons": reasons,
        "missing_required_fields": missing,
        "mapping_statuses": mapping_statuses,
        "ceiling_applied": ceiling,
    }
    return BridgeAssessment(
        result_hash=stable_hash(assessment_payload, "hbp_assessment_"),
        **assessment_payload,
    )


def assess_contracts(contracts: list[BridgeContract]) -> list[BridgeAssessment]:
    return [assess_contract(contract) for contract in contracts]
