from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


JsonDict = dict[str, Any]


@dataclass(frozen=True)
class PerturbationSpec:
    name: str
    family: str
    parameters: JsonDict
    source: str
    status: str = "available"
    unavailable_reason: str | None = None

    def to_dict(self) -> JsonDict:
        return asdict(self)


@dataclass(frozen=True)
class CSTRunProvenance:
    run_instance_id: str
    source_artifacts: dict[str, str]
    source_artifact_hashes: dict[str, str]
    config_version: str
    config_hash: str
    code_version: str
    candidate_label: str
    hierarchy_label: str
    n_side: int
    ensemble_size: int
    seed: int | None
    probe: str
    perturbation: str
    control_group: str

    def to_dict(self) -> JsonDict:
        return asdict(self)


@dataclass(frozen=True)
class ClosureSignature:
    branch_signature_id: str
    signature_version: str
    operational_equivalence_note: str
    fields: JsonDict
    unavailable_components: JsonDict = field(default_factory=dict)
    source_artifacts: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> JsonDict:
        return asdict(self)


@dataclass(frozen=True)
class RecoveryObservation:
    observation_id: str
    observation_kind: str
    control_role: str
    control_group_label: str
    hierarchy_label: str
    candidate_label: str
    n_side: int
    ensemble_size: int
    seed: int | None
    probe: str
    perturbation: PerturbationSpec
    provenance: CSTRunProvenance
    closure_signature: ClosureSignature
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> JsonDict:
        payload = asdict(self)
        payload["perturbation"] = self.perturbation.to_dict()
        payload["provenance"] = self.provenance.to_dict()
        payload["closure_signature"] = self.closure_signature.to_dict()
        return payload


@dataclass(frozen=True)
class ComponentDistance:
    name: str
    value: float | None
    availability: str
    normalization_method: str
    source_fields: list[str]
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> JsonDict:
        return asdict(self)


@dataclass(frozen=True)
class ClosureDistanceComponents:
    left_observation_id: str
    right_observation_id: str
    left_signature_id: str
    right_signature_id: str
    components: dict[str, ComponentDistance]
    scalar_distance: float | None
    scalar_enabled: bool
    weights_used: dict[str, float]
    disabled_components: list[str]
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> JsonDict:
        payload = asdict(self)
        payload["components"] = {
            name: component.to_dict()
            for name, component in self.components.items()
        }
        return payload


@dataclass(frozen=True)
class RecoverabilityVector:
    recovery_rate: float
    branch_identity_persistence: float
    closure_fidelity: float
    control_margin: float
    variance_bound: float
    ablation_balance: float
    authoritative_mode: str = "vector"
    optional_scalar: float | None = None
    scalar_formula: str | None = None
    scalar_warnings: list[str] = field(default_factory=list)
    generalized_cst_selectivity: JsonDict | None = None

    def to_dict(self) -> JsonDict:
        return asdict(self)


@dataclass(frozen=True)
class CSTBenchmarkResult:
    verdict: str
    reasons: list[str]
    warnings: list[str]
    recoverability_vector: RecoverabilityVector
    thresholds: JsonDict
    weight_sensitivity: JsonDict
    unavailable_controls: JsonDict
    artifact_hashes: dict[str, str]
    config_hash: str
    result_hash: str

    def to_dict(self) -> JsonDict:
        payload = asdict(self)
        payload["recoverability_vector"] = self.recoverability_vector.to_dict()
        return payload
