"""Auditable generative operators layered above frozen HAOS-IIP telemetry."""

from .spectral_tuning import (
    CandidateAudit,
    ExpansionAudit,
    RecoveryProfile,
    SpectralCandidate,
    TuningConfig,
    evaluate_address,
    expand_spectral_addresses,
    first_sustained_collapse,
    propose_spectral_candidates,
)

__all__ = [
    "CandidateAudit",
    "ExpansionAudit",
    "RecoveryProfile",
    "SpectralCandidate",
    "TuningConfig",
    "evaluate_address",
    "expand_spectral_addresses",
    "first_sustained_collapse",
    "propose_spectral_candidates",
]

from .hostile_validation import (
    HostileCandidateAudit,
    HostileValidationAudit,
    ValidationConfig,
    validate_expansion_hostile,
)

__all__ += [
    "HostileCandidateAudit",
    "HostileValidationAudit",
    "ValidationConfig",
    "validate_expansion_hostile",
]

from .proposal_tournament import (
    Proposal,
    Regime,
    SurvivorLedgerRow,
    TournamentAudit,
    TournamentConfig,
    propose_multi_seed_interference,
    propose_non_local_spectral_mixes,
    propose_recovery_guided_walks,
    run_proposal_tournament,
)

__all__ += [
    "Proposal",
    "Regime",
    "SurvivorLedgerRow",
    "TournamentAudit",
    "TournamentConfig",
    "propose_multi_seed_interference",
    "propose_non_local_spectral_mixes",
    "propose_recovery_guided_walks",
    "run_proposal_tournament",
]

from .structural_intervention import (
    EdgeInterventionConfig,
    EdgeUpdate,
    InterventionLedgerRow,
    StructuralInterventionAudit,
    run_edge_intervention_experiment,
    select_directed_edge_updates,
)

__all__ += [
    "EdgeInterventionConfig",
    "EdgeUpdate",
    "InterventionLedgerRow",
    "StructuralInterventionAudit",
    "run_edge_intervention_experiment",
    "select_directed_edge_updates",
]

from .status_semantics import ResultSemantics, interpret_result_status

__all__ += ["ResultSemantics", "interpret_result_status"]
