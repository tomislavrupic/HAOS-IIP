from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


JsonDict = dict[str, Any]


BRIDGE_LEVELS = (
    "FORMAL_CORRESPONDENCE",
    "ANALOGY_BRIDGE",
    "OPERATIONAL_MAPPING",
    "PREDICTIVE_BRIDGE",
    "EMPIRICALLY_SUPPORTED_BRIDGE",
    "PHYSICAL_MECHANISM_CANDIDATE",
)

MAPPING_STATUSES = ("DERIVED", "CALIBRATED", "HEURISTIC", "ANALOGICAL", "UNAVAILABLE")


@dataclass(frozen=True)
class DomainStateSpec:
    variables: list[str]
    state_space: str
    units: dict[str, str]
    valid_ranges: dict[str, str]
    missing_data_policy: str

    def to_dict(self) -> JsonDict:
        return asdict(self)


@dataclass(frozen=True)
class DynamicsSpec:
    update_law: str
    parameters: JsonDict
    perturbations: list[str]
    time_axis: str
    boundary_conditions: str

    def to_dict(self) -> JsonDict:
        return asdict(self)


@dataclass(frozen=True)
class HAOSMapping:
    haos_source: str
    domain_target: str
    mapping_function: str
    units_before: str
    units_after: str
    semantic_justification: str
    mapping_status: str
    uncertainty: str
    failure_conditions: list[str]

    def to_dict(self) -> JsonDict:
        return asdict(self)


@dataclass(frozen=True)
class ObservationMap:
    observable: str
    measurement_rule: str
    units: str
    valid_range: str
    target_role: str

    def to_dict(self) -> JsonDict:
        return asdict(self)


@dataclass(frozen=True)
class BaselineSpec:
    name: str
    family: str
    prediction_rule: str
    training_scope: str
    leakage_risk: str

    def to_dict(self) -> JsonDict:
        return asdict(self)


@dataclass(frozen=True)
class ControlSpec:
    name: str
    purpose: str
    preserves: list[str]
    destroys: list[str]
    invalid_if: str

    def to_dict(self) -> JsonDict:
        return asdict(self)


@dataclass(frozen=True)
class FalsificationSpec:
    fail_conditions: list[str]
    invalidation_conditions: list[str]
    underpowered_policy: str
    leakage_policy: str

    def to_dict(self) -> JsonDict:
        return asdict(self)


@dataclass(frozen=True)
class BridgeProvenance:
    source_artifacts: list[str]
    code_artifacts: list[str]
    precommitment_artifact: str | None
    external_data_status: str
    level_history: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> JsonDict:
        return asdict(self)


@dataclass(frozen=True)
class BridgeContract:
    bridge_id: str
    domain: str
    classification: str
    state_schema: DomainStateSpec | None
    units: dict[str, str]
    dynamics: DynamicsSpec | None
    boundary_conditions: str
    symmetries: list[str]
    haos_mapping: list[HAOSMapping]
    observation_map: list[ObservationMap]
    prediction_target: str
    calibration_split: str
    holdout_split: str
    controls: list[ControlSpec]
    baselines: list[BaselineSpec]
    uncertainty: str
    falsification: FalsificationSpec | None
    verdict_logic: str
    provenance: BridgeProvenance

    def to_dict(self) -> JsonDict:
        payload = asdict(self)
        payload["state_schema"] = self.state_schema.to_dict() if self.state_schema else None
        payload["dynamics"] = self.dynamics.to_dict() if self.dynamics else None
        payload["haos_mapping"] = [item.to_dict() for item in self.haos_mapping]
        payload["observation_map"] = [item.to_dict() for item in self.observation_map]
        payload["controls"] = [item.to_dict() for item in self.controls]
        payload["baselines"] = [item.to_dict() for item in self.baselines]
        payload["falsification"] = self.falsification.to_dict() if self.falsification else None
        payload["provenance"] = self.provenance.to_dict()
        return payload


@dataclass(frozen=True)
class BridgeAssessment:
    bridge_id: str
    requested_classification: str
    effective_classification: str
    verdict_authorized: bool
    labels: list[str]
    reasons: list[str]
    missing_required_fields: list[str]
    mapping_statuses: dict[str, str]
    ceiling_applied: str
    result_hash: str

    def to_dict(self) -> JsonDict:
        return asdict(self)
