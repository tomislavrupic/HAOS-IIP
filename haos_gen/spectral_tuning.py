"""Directed spectral-address expansion with frozen recoverability telemetry.

This module does not modify ``haos_core`` or any frozen phase artifact.  It
proposes addresses from the eigensystem of a declared symmetric operator and
audits their reconstruction across caller-supplied perturbations.

"Resonance" has one operational meaning here: normalized overlap between an
existing address and an operator eigenvector.  It carries no physical claim.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable

import numpy as np

from telemetry.frozen_metrics import (
    SurvivalThresholds,
    overlap,
    persistence_time,
    recovery_score,
)


_EPS = 1.0e-12


@dataclass(frozen=True)
class TuningConfig:
    """Predeclared generation and selection gates."""

    max_candidates: int = 6
    support_fraction: float = 0.40
    collapse_threshold: float = 0.55
    sustain_steps: int = 2
    min_resonance: float = 0.10
    min_novelty: float = 0.05
    min_mean_gain: float = 0.01
    max_mean_regression: float = 0.005
    control_margin: float = 0.01


@dataclass(frozen=True)
class SpectralCandidate:
    candidate_id: str
    eigenvalue: float
    resonance: float
    locality: float
    novelty: float
    address: np.ndarray
    null_address: np.ndarray


@dataclass(frozen=True)
class RecoveryProfile:
    scores: tuple[float, ...]
    mean_recovery: float
    minimum_recovery: float
    persistence_time: float
    k_star: int | None


@dataclass(frozen=True)
class CandidateAudit:
    candidate_id: str
    eigenvalue: float
    resonance: float
    locality: float
    novelty: float
    candidate_profile: RecoveryProfile
    null_profile: RecoveryProfile
    accepted: bool
    failure_modes: tuple[str, ...]


@dataclass(frozen=True)
class ExpansionAudit:
    baseline_profile: RecoveryProfile
    candidates: tuple[CandidateAudit, ...]
    retained_addresses: tuple[np.ndarray, ...]
    contract: dict[str, object]

    def to_jsonable(self) -> dict[str, object]:
        """Return metrics and decisions without embedding dense address arrays."""

        return {
            "baseline_profile": asdict(self.baseline_profile),
            "candidates": [asdict(candidate) for candidate in self.candidates],
            "retained_count": len(self.retained_addresses),
            "contract": self.contract,
        }


def _normalize(vector: np.ndarray) -> np.ndarray:
    values = np.asarray(vector, dtype=float)
    norm = float(np.linalg.norm(values))
    if norm <= _EPS:
        raise ValueError("address vectors must have non-zero norm")
    return values / norm


def _validate_operator(operator: np.ndarray) -> np.ndarray:
    values = np.asarray(operator, dtype=float)
    if values.ndim != 2 or values.shape[0] != values.shape[1]:
        raise ValueError("operator must be a square matrix")
    if not np.allclose(values, values.T, atol=1.0e-10):
        raise ValueError("spectral tuning requires a symmetric operator")
    return values


def _support_mask(address: np.ndarray, fraction: float) -> np.ndarray:
    if not 0.0 < float(fraction) <= 1.0:
        raise ValueError("support_fraction must lie in (0, 1]")
    count = max(1, int(np.ceil(address.size * float(fraction))))
    indices = np.argsort(np.abs(address))[-count:]
    mask = np.zeros(address.size, dtype=bool)
    mask[indices] = True
    return mask


def _locality(address: np.ndarray) -> float:
    probabilities = np.abs(_normalize(address)) ** 2
    return float(np.sum(probabilities * probabilities))


def _matched_null(address: np.ndarray, index: int) -> np.ndarray:
    """Deterministic coordinate-null preserving norm and coefficient values."""

    size = address.size
    shift = 1 + (2 * index) % max(size - 1, 1)
    null = np.roll(address, shift)
    if overlap(address, null) > 1.0 - 1.0e-10 and size > 2:
        null = address[np.r_[1:size:2, 0:size:2]]
    return _normalize(null)


def propose_spectral_candidates(
    operator: np.ndarray,
    seed_address: np.ndarray,
    config: TuningConfig = TuningConfig(),
) -> tuple[SpectralCandidate, ...]:
    """Propose directed candidates from modes aligned with the seed address.

    Candidates are ranked by resonance times locality.  Random sampling is not
    used.  The constant/zero mode is skipped when other modes are available.
    """

    matrix = _validate_operator(operator)
    seed = _normalize(seed_address)
    if seed.size != matrix.shape[0]:
        raise ValueError("seed_address length must match operator dimension")

    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    start = 1 if eigenvalues.size > 1 and abs(float(eigenvalues[0])) <= 1.0e-10 else 0
    ranked: list[tuple[float, SpectralCandidate]] = []
    for index in range(start, eigenvalues.size):
        address = _normalize(eigenvectors[:, index])
        signed_overlap = float(np.dot(seed, address))
        if signed_overlap < 0.0:
            address = -address
        resonance = overlap(seed, address)
        novelty = 1.0 - resonance
        locality = _locality(address)
        if resonance + _EPS < config.min_resonance or novelty + _EPS < config.min_novelty:
            continue
        candidate = SpectralCandidate(
            candidate_id=f"spectral_mode_{index:03d}",
            eigenvalue=float(eigenvalues[index]),
            resonance=float(resonance),
            locality=float(locality),
            novelty=float(novelty),
            address=address,
            null_address=_matched_null(address, index),
        )
        ranked.append((resonance * locality, candidate))

    ranked.sort(key=lambda item: (-item[0], item[1].eigenvalue, item[1].candidate_id))
    return tuple(candidate for _, candidate in ranked[: config.max_candidates])


def first_sustained_collapse(
    scores: Iterable[float],
    threshold: float,
    sustain_steps: int,
) -> int | None:
    """Return the first index beginning a sustained below-threshold run."""

    values = tuple(float(value) for value in scores)
    if sustain_steps < 1:
        raise ValueError("sustain_steps must be at least one")
    for index in range(0, len(values) - sustain_steps + 1):
        window = values[index : index + sustain_steps]
        if all(value < float(threshold) for value in window):
            return index
    return None


def _reconstruct(operator: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Reidentify the perturbed eigenvector with maximum reference overlap."""

    _, eigenvectors = np.linalg.eigh(_validate_operator(operator))
    overlaps = np.abs(eigenvectors.T @ _normalize(reference))
    reconstructed = eigenvectors[:, int(np.argmax(overlaps))]
    if float(np.dot(reference, reconstructed)) < 0.0:
        reconstructed = -reconstructed
    return _normalize(reconstructed)


def _profile(
    address: np.ndarray,
    perturbations: tuple[np.ndarray, ...],
    stress_levels: tuple[float, ...],
    thresholds: SurvivalThresholds,
    config: TuningConfig,
) -> RecoveryProfile:
    reference = _normalize(address)
    coords = np.arange(reference.size, dtype=float)[:, None]
    mask = _support_mask(reference, config.support_fraction)
    history: list[dict[str, object]] = []
    scores: list[float] = []
    for perturbed_operator in perturbations:
        state = _reconstruct(perturbed_operator, reference)
        scores.append(recovery_score(reference, state, coords, mask))
        history.append(
            {
                "reference": reference,
                "state": state,
                "coords": coords,
                "mask": mask,
            }
        )
    return RecoveryProfile(
        scores=tuple(float(value) for value in scores),
        mean_recovery=float(np.mean(scores)),
        minimum_recovery=float(np.min(scores)),
        persistence_time=persistence_time(history, list(stress_levels), thresholds),
        k_star=first_sustained_collapse(
            scores,
            threshold=config.collapse_threshold,
            sustain_steps=config.sustain_steps,
        ),
    )


def evaluate_address(
    address: np.ndarray,
    perturbed_operators: Iterable[np.ndarray],
    stress_levels: Iterable[float],
    thresholds: SurvivalThresholds,
    config: TuningConfig = TuningConfig(),
) -> RecoveryProfile:
    """Evaluate one declared address with the same frozen V0 telemetry path."""

    perturbations = tuple(_validate_operator(value) for value in perturbed_operators)
    levels = tuple(float(value) for value in stress_levels)
    if not perturbations:
        raise ValueError("at least one perturbed operator is required")
    if len(perturbations) != len(levels):
        raise ValueError("perturbed_operators and stress_levels must have equal length")
    if any(right < left for left, right in zip(levels, levels[1:])):
        raise ValueError("stress_levels must be monotonically non-decreasing")
    return _profile(address, perturbations, levels, thresholds, config)


def _not_earlier(candidate: int | None, baseline: int | None) -> bool:
    if baseline is None:
        return candidate is None
    return candidate is None or candidate >= baseline


def expand_spectral_addresses(
    operator: np.ndarray,
    seed_address: np.ndarray,
    perturbed_operators: Iterable[np.ndarray],
    stress_levels: Iterable[float],
    thresholds: SurvivalThresholds,
    config: TuningConfig = TuningConfig(),
) -> ExpansionAudit:
    """Generate, stress, and select spectral addresses.

    Input contract:
      symmetric baseline operator, non-zero seed address, ordered symmetric
      perturbed operators, matching monotonically increasing stress levels,
      frozen survival thresholds, and predeclared selection gates.

    Output contract:
      immutable baseline/candidate/null profiles, decisions, failure labels,
      and retained address vectors.  No input or frozen artifact is modified.
    """

    matrix = _validate_operator(operator)
    seed = _normalize(seed_address)
    perturbations = tuple(_validate_operator(value) for value in perturbed_operators)
    levels = tuple(float(value) for value in stress_levels)
    if not perturbations:
        raise ValueError("at least one perturbed operator is required")
    if len(perturbations) != len(levels):
        raise ValueError("perturbed_operators and stress_levels must have equal length")
    if any(right < left for left, right in zip(levels, levels[1:])):
        raise ValueError("stress_levels must be monotonically non-decreasing")
    if any(value.shape != matrix.shape for value in perturbations):
        raise ValueError("all perturbed operators must match baseline shape")

    baseline = _profile(seed, perturbations, levels, thresholds, config)
    audits: list[CandidateAudit] = []
    retained: list[np.ndarray] = []
    for candidate in propose_spectral_candidates(matrix, seed, config):
        profile = _profile(candidate.address, perturbations, levels, thresholds, config)
        null_profile = _profile(candidate.null_address, perturbations, levels, thresholds, config)
        failures: list[str] = []

        if profile.mean_recovery < baseline.mean_recovery - config.max_mean_regression:
            failures.append("RECOVERABILITY_REGRESSION")
        if profile.persistence_time + _EPS < baseline.persistence_time:
            failures.append("PERSISTENCE_REGRESSION")
        if not _not_earlier(profile.k_star, baseline.k_star):
            failures.append("EARLIER_SUSTAINED_COLLAPSE")
        if profile.mean_recovery < null_profile.mean_recovery + config.control_margin:
            failures.append("NO_MATCHED_CONTROL_SEPARATION")

        improves = profile.mean_recovery >= baseline.mean_recovery + config.min_mean_gain
        stabilizes = (
            profile.mean_recovery >= baseline.mean_recovery - config.max_mean_regression
            and profile.persistence_time >= baseline.persistence_time
            and _not_earlier(profile.k_star, baseline.k_star)
        )
        if not (improves or stabilizes):
            failures.append("NO_RECOVERABLE_EXPANSION")

        accepted = not failures
        if accepted:
            retained.append(candidate.address.copy())
        audits.append(
            CandidateAudit(
                candidate_id=candidate.candidate_id,
                eigenvalue=candidate.eigenvalue,
                resonance=candidate.resonance,
                locality=candidate.locality,
                novelty=candidate.novelty,
                candidate_profile=profile,
                null_profile=null_profile,
                accepted=accepted,
                failure_modes=tuple(failures),
            )
        )

    contract = {
        "generator": "directed_spectral_tuning_v0",
        "resonance_definition": "absolute normalized overlap with seed address",
        "candidate_source": "baseline symmetric-operator eigensystem",
        "matched_null": "deterministic coordinate permutation preserving coefficient multiset and norm",
        "selection": "non-regression in recovery, persistence, and k_star plus matched-control separation",
        "frozen_metrics_used": ["recovery_score", "persistence_time"],
        "mutates_frozen_state": False,
    }
    return ExpansionAudit(
        baseline_profile=baseline,
        candidates=tuple(audits),
        retained_addresses=tuple(retained),
        contract=contract,
    )
