"""HAOS-GEN V0.2: sealed held-out validation and hostile null families."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable

import numpy as np

from telemetry.frozen_metrics import SurvivalThresholds, overlap

from .spectral_tuning import (
    ExpansionAudit,
    RecoveryProfile,
    TuningConfig,
    evaluate_address,
    expand_spectral_addresses,
    propose_spectral_candidates,
)


_EPS = 1.0e-12


@dataclass(frozen=True)
class ValidationConfig:
    """Predeclared V0.2 gates; held-out data never change proposal ranking."""

    min_heldout_mean_gain: float = 0.01
    min_persistence_gain: float = 0.0
    min_k_star_displacement: int = 1
    hostile_control_margin: float = 0.01
    random_seed: int = 20260724


@dataclass(frozen=True)
class HostileCandidateAudit:
    candidate_id: str
    tuning_accepted: bool
    heldout_accepted: bool
    heldout_profile: RecoveryProfile
    coordinate_null_profile: RecoveryProfile
    degree_matched_null_profile: RecoveryProfile
    spectrum_preserving_null_profile: RecoveryProfile
    mean_recovery_gain: float
    persistence_gain: float
    k_star_displacement: int
    amplification_pass: bool
    accepted: bool
    failure_modes: tuple[str, ...]


@dataclass(frozen=True)
class RandomBaselineAudit:
    baseline_id: str
    mode_index: int
    tuning_profile: RecoveryProfile
    heldout_profile: RecoveryProfile
    hostile_null_ceiling: float
    accepted: bool
    failure_modes: tuple[str, ...]


@dataclass(frozen=True)
class HostileValidationAudit:
    tuning_audit: ExpansionAudit
    heldout_baseline_profile: RecoveryProfile
    candidates: tuple[HostileCandidateAudit, ...]
    random_baseline: tuple[RandomBaselineAudit, ...]
    yield_summary: dict[str, object]
    retained_addresses: tuple[np.ndarray, ...]
    contract: dict[str, object]

    def to_jsonable(self) -> dict[str, object]:
        return {
            "tuning_audit": self.tuning_audit.to_jsonable(),
            "heldout_baseline_profile": asdict(self.heldout_baseline_profile),
            "candidates": [asdict(candidate) for candidate in self.candidates],
            "random_baseline": [asdict(candidate) for candidate in self.random_baseline],
            "yield_summary": self.yield_summary,
            "retained_count": len(self.retained_addresses),
            "contract": self.contract,
        }


def _normalize(vector: np.ndarray) -> np.ndarray:
    values = np.asarray(vector, dtype=float)
    norm = float(np.linalg.norm(values))
    if norm <= _EPS:
        raise ValueError("address vectors must have non-zero norm")
    return values / norm


def _degree_matched_null(operator: np.ndarray, address: np.ndarray) -> np.ndarray:
    """Permute coefficients only within equal diagonal/degree classes."""

    diagonal = np.diag(np.asarray(operator, dtype=float))
    null = np.asarray(address, dtype=float).copy()
    groups: dict[float, list[int]] = {}
    for index, value in enumerate(diagonal):
        groups.setdefault(round(float(value), 10), []).append(index)
    for indices in groups.values():
        if len(indices) > 1:
            values = null[indices].copy()
            null[indices] = np.roll(values, 1)
    if overlap(address, null) > 1.0 - 1.0e-10:
        null = np.roll(null, 1)
    return _normalize(null)


def _coordinate_null(address: np.ndarray, index: int) -> np.ndarray:
    size = address.size
    return _normalize(np.roll(address, 1 + (2 * index) % max(size - 1, 1)))


def _orthogonal_scrambler(size: int) -> np.ndarray:
    """Fixed orthogonal basis change; deterministic and seed-free."""

    rows = np.arange(1, size + 1, dtype=float)[:, None]
    cols = np.arange(1, size + 1, dtype=float)[None, :]
    source = np.sin(rows * cols) + np.cos(rows * (cols + 0.5))
    orthogonal, _ = np.linalg.qr(source)
    return orthogonal


def _spectrum_preserving_nulls(operators: tuple[np.ndarray, ...]) -> tuple[np.ndarray, ...]:
    scrambler = _orthogonal_scrambler(operators[0].shape[0])
    nulls: list[np.ndarray] = []
    for operator in operators:
        eigenvalues = np.linalg.eigvalsh(operator)
        null = scrambler @ np.diag(eigenvalues) @ scrambler.T
        nulls.append(0.5 * (null + null.T))
    return tuple(nulls)


def _k_star_displacement(
    candidate: RecoveryProfile,
    baseline: RecoveryProfile,
) -> int:
    candidate_index = len(candidate.scores) if candidate.k_star is None else int(candidate.k_star)
    baseline_index = len(baseline.scores) if baseline.k_star is None else int(baseline.k_star)
    return candidate_index - baseline_index


def _amplifies(
    profile: RecoveryProfile,
    baseline: RecoveryProfile,
    config: ValidationConfig,
) -> tuple[bool, float, float, int]:
    mean_gain = profile.mean_recovery - baseline.mean_recovery
    persistence_gain = profile.persistence_time - baseline.persistence_time
    k_displacement = _k_star_displacement(profile, baseline)
    passes = (
        mean_gain >= config.min_heldout_mean_gain
        or persistence_gain > config.min_persistence_gain + _EPS
        or k_displacement >= config.min_k_star_displacement
    )
    return passes, float(mean_gain), float(persistence_gain), int(k_displacement)


def _profile_non_regression(
    profile: RecoveryProfile,
    baseline: RecoveryProfile,
    tuning_config: TuningConfig,
) -> bool:
    return (
        profile.mean_recovery
        >= baseline.mean_recovery - tuning_config.max_mean_regression
        and profile.persistence_time + _EPS >= baseline.persistence_time
        and _k_star_displacement(profile, baseline) >= 0
    )


def _random_eigenmodes(
    operator: np.ndarray,
    budget: int,
    seed: int,
    excluded_indices: frozenset[int] = frozenset(),
) -> tuple[tuple[int, np.ndarray], ...]:
    eigenvalues, eigenvectors = np.linalg.eigh(np.asarray(operator, dtype=float))
    start = 1 if eigenvalues.size > 1 and abs(float(eigenvalues[0])) <= 1.0e-10 else 0
    indices = np.asarray(
        [index for index in range(start, eigenvalues.size) if index not in excluded_indices],
        dtype=int,
    )
    if budget <= 0 or indices.size == 0:
        return ()
    generator = np.random.default_rng(seed)
    selected = generator.choice(indices, size=min(budget, indices.size), replace=False)
    return tuple((int(index), _normalize(eigenvectors[:, int(index)])) for index in selected)


def _yield_curve(audits: list[HostileCandidateAudit]) -> list[dict[str, float | int]]:
    curve: list[dict[str, float | int]] = []
    survivors = 0
    for proposals, audit in enumerate(audits, start=1):
        survivors += int(audit.accepted)
        curve.append(
            {
                "proposals": proposals,
                "survivors": survivors,
                "cumulative_yield": float(survivors / proposals),
            }
        )
    return curve


def validate_expansion_hostile(
    operator: np.ndarray,
    seed_address: np.ndarray,
    tuning_operators: Iterable[np.ndarray],
    tuning_levels: Iterable[float],
    heldout_operators: Iterable[np.ndarray],
    heldout_levels: Iterable[float],
    thresholds: SurvivalThresholds,
    tuning_config: TuningConfig = TuningConfig(),
    validation_config: ValidationConfig = ValidationConfig(),
) -> HostileValidationAudit:
    """Run V0 generation, then require sealed-family amplification and null wins."""

    matrix = np.asarray(operator, dtype=float)
    tuning_ops = tuple(np.asarray(value, dtype=float) for value in tuning_operators)
    tuning_tau = tuple(float(value) for value in tuning_levels)
    heldout_ops = tuple(np.asarray(value, dtype=float) for value in heldout_operators)
    heldout_tau = tuple(float(value) for value in heldout_levels)
    if not heldout_ops:
        raise ValueError("heldout_operators must not be empty")
    if len(heldout_ops) != len(heldout_tau):
        raise ValueError("heldout_operators and heldout_levels must have equal length")

    tuning_audit = expand_spectral_addresses(
        matrix,
        seed_address,
        tuning_ops,
        tuning_tau,
        thresholds,
        tuning_config,
    )
    heldout_audit = expand_spectral_addresses(
        matrix,
        seed_address,
        heldout_ops,
        heldout_tau,
        thresholds,
        tuning_config,
    )
    heldout_by_id = {candidate.candidate_id: candidate for candidate in heldout_audit.candidates}
    tuning_by_id = {candidate.candidate_id: candidate for candidate in tuning_audit.candidates}
    proposals = propose_spectral_candidates(matrix, seed_address, tuning_config)
    spectrum_nulls = _spectrum_preserving_nulls(heldout_ops)

    audits: list[HostileCandidateAudit] = []
    retained: list[np.ndarray] = []
    for proposal in proposals:
        tuning_result = tuning_by_id[proposal.candidate_id]
        heldout_result = heldout_by_id[proposal.candidate_id]
        degree_null = _degree_matched_null(matrix, proposal.address)
        degree_profile = evaluate_address(
            degree_null, heldout_ops, heldout_tau, thresholds, tuning_config
        )
        spectrum_profile = evaluate_address(
            proposal.address, spectrum_nulls, heldout_tau, thresholds, tuning_config
        )
        amplification, mean_gain, persistence_gain, k_displacement = _amplifies(
            heldout_result.candidate_profile,
            heldout_audit.baseline_profile,
            validation_config,
        )
        hostile_ceiling = max(
            heldout_result.null_profile.mean_recovery,
            degree_profile.mean_recovery,
            spectrum_profile.mean_recovery,
        )
        failures: list[str] = []
        if not tuning_result.accepted:
            failures.append("TUNING_GATE_FAILED")
        if not heldout_result.accepted:
            failures.append("HELDOUT_RECOVERABILITY_GATE_FAILED")
        if (
            heldout_result.candidate_profile.mean_recovery
            < hostile_ceiling + validation_config.hostile_control_margin
        ):
            failures.append("HOSTILE_NULL_SEPARATION_FAILED")
        if not amplification:
            failures.append("NO_HELDOUT_AMPLIFICATION")
        accepted = not failures
        if accepted:
            retained.append(proposal.address.copy())
        audits.append(
            HostileCandidateAudit(
                candidate_id=proposal.candidate_id,
                tuning_accepted=tuning_result.accepted,
                heldout_accepted=heldout_result.accepted,
                heldout_profile=heldout_result.candidate_profile,
                coordinate_null_profile=heldout_result.null_profile,
                degree_matched_null_profile=degree_profile,
                spectrum_preserving_null_profile=spectrum_profile,
                mean_recovery_gain=mean_gain,
                persistence_gain=persistence_gain,
                k_star_displacement=k_displacement,
                amplification_pass=amplification,
                accepted=accepted,
                failure_modes=tuple(failures),
            )
        )

    random_audits: list[RandomBaselineAudit] = []
    proposed_mode_indices = frozenset(
        int(proposal.candidate_id.rsplit("_", 1)[1]) for proposal in proposals
    )
    for index, address in _random_eigenmodes(
        matrix,
        budget=len(proposals),
        seed=validation_config.random_seed,
        excluded_indices=proposed_mode_indices,
    ):
        tuning_profile = evaluate_address(
            address, tuning_ops, tuning_tau, thresholds, tuning_config
        )
        heldout_profile = evaluate_address(
            address, heldout_ops, heldout_tau, thresholds, tuning_config
        )
        amplification, _, _, _ = _amplifies(
            heldout_profile, heldout_audit.baseline_profile, validation_config
        )
        coordinate_profile = evaluate_address(
            _coordinate_null(address, index),
            heldout_ops,
            heldout_tau,
            thresholds,
            tuning_config,
        )
        degree_profile = evaluate_address(
            _degree_matched_null(matrix, address),
            heldout_ops,
            heldout_tau,
            thresholds,
            tuning_config,
        )
        spectrum_profile = evaluate_address(
            address, spectrum_nulls, heldout_tau, thresholds, tuning_config
        )
        hostile_ceiling = max(
            coordinate_profile.mean_recovery,
            degree_profile.mean_recovery,
            spectrum_profile.mean_recovery,
        )
        failures: list[str] = []
        if not _profile_non_regression(
            tuning_profile, tuning_audit.baseline_profile, tuning_config
        ):
            failures.append("TUNING_GATE_FAILED")
        if not _profile_non_regression(
            heldout_profile, heldout_audit.baseline_profile, tuning_config
        ):
            failures.append("HELDOUT_RECOVERABILITY_GATE_FAILED")
        if (
            heldout_profile.mean_recovery
            < hostile_ceiling + validation_config.hostile_control_margin
        ):
            failures.append("HOSTILE_NULL_SEPARATION_FAILED")
        if not amplification:
            failures.append("NO_HELDOUT_AMPLIFICATION")
        accepted = not failures
        random_audits.append(
            RandomBaselineAudit(
                baseline_id=f"random_eigenmode_{index:03d}",
                mode_index=index,
                tuning_profile=tuning_profile,
                heldout_profile=heldout_profile,
                hostile_null_ceiling=float(hostile_ceiling),
                accepted=accepted,
                failure_modes=tuple(failures),
            )
        )

    proposed_count = len(audits)
    tuning_survivors = sum(int(item.tuning_accepted) for item in audits)
    heldout_survivors = sum(
        int(item.tuning_accepted and item.heldout_accepted) for item in audits
    )
    final_survivors = sum(int(item.accepted) for item in audits)
    random_survivors = sum(int(item.accepted) for item in random_audits)
    directed_yield = float(final_survivors / max(proposed_count, 1))
    random_yield = float(random_survivors / max(len(random_audits), 1))
    directed_advantage = directed_yield - random_yield
    if final_survivors == 0:
        probe_status = "FAIL_NO_HOSTILE_SURVIVORS"
    elif directed_advantage <= 0.0:
        probe_status = "OPEN_NO_DIRECTED_YIELD_ADVANTAGE"
    else:
        probe_status = "PASS_BOUNDED_DIRECTED_YIELD_ADVANTAGE"
    yield_summary: dict[str, object] = {
        "probe_status": probe_status,
        "proposed_count": proposed_count,
        "tuning_survivors": tuning_survivors,
        "heldout_survivors": heldout_survivors,
        "final_survivors": final_survivors,
        "tuning_yield": float(tuning_survivors / max(proposed_count, 1)),
        "heldout_yield": float(heldout_survivors / max(proposed_count, 1)),
        "final_yield": directed_yield,
        "proposals_per_survivor": (
            None if final_survivors == 0 else float(proposed_count / final_survivors)
        ),
        "random_budget": len(random_audits),
        "random_survivors": random_survivors,
        "random_yield": random_yield,
        "directed_yield_advantage": directed_advantage,
        "cumulative_yield_curve": _yield_curve(audits),
    }
    contract = {
        "validator": "hostile_validation_v0_2",
        "proposal_family": "directed_spectral_tuning_v0",
        "heldout_family_used_for_proposals": False,
        "null_families": [
            "coordinate_permutation",
            "degree_class_permutation",
            "spectrum_preserving_basis_scramble",
        ],
        "random_baseline": "seeded equal-budget random eigenmodes excluding directed proposals",
        "amplification_required": True,
        "frozen_metrics_used": ["recovery_score", "persistence_time"],
        "mutates_frozen_state": False,
    }
    return HostileValidationAudit(
        tuning_audit=tuning_audit,
        heldout_baseline_profile=heldout_audit.baseline_profile,
        candidates=tuple(audits),
        random_baseline=tuple(random_audits),
        yield_summary=yield_summary,
        retained_addresses=tuple(retained),
        contract=contract,
    )
