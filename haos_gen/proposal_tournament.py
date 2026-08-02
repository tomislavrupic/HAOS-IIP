"""HAOS-GEN V0.3 proposal-operator tournament.

V0.3 changes proposal construction only.  Candidate evaluation delegates to
the frozen V0/V0.2 metric and null helpers and records a parity-comparable
ledger for directed operators and both mandatory baselines.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from itertools import combinations
from typing import Iterable, Mapping

import numpy as np

from telemetry.frozen_metrics import SurvivalThresholds, overlap

from .hostile_validation import (
    ValidationConfig,
    _amplifies,
    _coordinate_null,
    _degree_matched_null,
    _profile_non_regression,
    _random_eigenmodes,
    _spectrum_preserving_nulls,
)
from .spectral_tuning import (
    RecoveryProfile,
    TuningConfig,
    evaluate_address,
    propose_spectral_candidates,
)


_EPS = 1.0e-12
_DIRECTED_OPERATORS = (
    "non_local_spectral_mixing_v0_3",
    "multi_seed_interference_v0_3",
    "recovery_guided_graph_walk_v0_3",
)
_BASELINES = ("v0_local_resonance", "random_eigenmode")


@dataclass(frozen=True)
class TournamentConfig:
    proposal_budget: int = 6
    nonlocal_alphas: tuple[float, ...] = (0.3, 0.5, 0.7)
    min_spectral_separation: float = 0.20
    multi_seed_k: int = 2
    multi_seed_rule: str = "equal_weight_signed"
    candidate_seed_overlap_limit: float = 0.94
    walk_length: int = 4
    walk_bias_recovery_weight: float = 0.75
    walk_extraction_radius: int = 1
    address_distance_threshold: float = 0.08
    subspace_overlap_threshold: float = 0.92
    minimum_positive_seeds: int = 2
    hard_regime_name: str = "hard"
    random_seed: int = 20260724


@dataclass(frozen=True)
class Regime:
    name: str
    tuning_operators: tuple[np.ndarray, ...]
    tuning_levels: tuple[float, ...]
    heldout_operators: tuple[np.ndarray, ...]
    heldout_levels: tuple[float, ...]


@dataclass(frozen=True)
class Proposal:
    operator: str
    candidate_id: str
    address: np.ndarray
    source_ids: tuple[str, ...]
    metadata: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class SurvivorLedgerRow:
    operator: str
    seed_id: str
    regime: str
    candidate_id: str
    tuning_decision: str
    heldout_decision: str
    delta_recovery: float
    delta_persistence: float
    k_star_displacement: int
    coordinate_null_margin: float
    degree_null_margin: float
    spectrum_null_margin: float
    baseline_match: str
    unique_flag: bool
    accepted: bool
    notes: tuple[str, ...]
    address: np.ndarray = field(repr=False, compare=False)

    def to_jsonable(self) -> dict[str, object]:
        payload = asdict(self)
        payload.pop("address")
        return payload


@dataclass(frozen=True)
class TournamentAudit:
    rows: tuple[SurvivorLedgerRow, ...]
    aggregate: dict[str, object]
    contract: dict[str, object]

    def to_jsonable(self) -> dict[str, object]:
        return {
            "rows": [row.to_jsonable() for row in self.rows],
            "aggregate": self.aggregate,
            "contract": self.contract,
        }


def _normalize(vector: np.ndarray) -> np.ndarray:
    values = np.asarray(vector, dtype=float)
    norm = float(np.linalg.norm(values))
    if norm <= _EPS:
        raise ValueError("candidate address must have non-zero norm")
    return values / norm


def _rayleigh(operator: np.ndarray, address: np.ndarray) -> float:
    vector = _normalize(address)
    return float(vector @ np.asarray(operator, dtype=float) @ vector)


def _sign_align(reference: np.ndarray, candidate: np.ndarray) -> np.ndarray:
    values = _normalize(candidate)
    return -values if float(np.dot(reference, values)) < 0.0 else values


def _deduplicate(proposals: Iterable[Proposal]) -> tuple[Proposal, ...]:
    retained: list[Proposal] = []
    for proposal in proposals:
        if any(overlap(proposal.address, prior.address) > 1.0 - 1.0e-9 for prior in retained):
            continue
        retained.append(proposal)
    return tuple(retained)


def _fill_budget(proposals: tuple[Proposal, ...], budget: int) -> tuple[Proposal, ...]:
    """Fail instead of silently assigning unequal proposal budgets."""

    if budget < 1:
        raise ValueError("proposal_budget must be at least one")
    if len(proposals) < budget:
        raise ValueError(
            f"operator produced {len(proposals)} distinct candidates; "
            f"fixed budget requires {budget}"
        )
    return proposals[:budget]


def propose_non_local_spectral_mixes(
    operator: np.ndarray,
    seed_id: str,
    seed_address: np.ndarray,
    tuning_operators: tuple[np.ndarray, ...],
    tuning_levels: tuple[float, ...],
    thresholds: SurvivalThresholds,
    tuning_config: TuningConfig,
    config: TournamentConfig,
) -> tuple[Proposal, ...]:
    """Mix the seed with spectrally distant, tuning-recoverable modes."""

    matrix = np.asarray(operator, dtype=float)
    seed = _normalize(seed_address)
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    span = max(float(np.ptp(eigenvalues)), _EPS)
    seed_lambda = _rayleigh(matrix, seed)
    seed_profile = evaluate_address(
        seed, tuning_operators, tuning_levels, thresholds, tuning_config
    )
    distant: list[tuple[float, float, int, np.ndarray]] = []
    for index in range(1, eigenvalues.size):
        separation = abs(float(eigenvalues[index]) - seed_lambda) / span
        if separation + _EPS < config.min_spectral_separation:
            continue
        mode = _sign_align(seed, eigenvectors[:, index])
        profile = evaluate_address(
            mode, tuning_operators, tuning_levels, thresholds, tuning_config
        )
        distant.append((profile.mean_recovery, separation, index, mode))
    if not distant:
        raise ValueError("no modes satisfy the declared non-local separation")
    eligible = [entry for entry in distant if entry[0] >= seed_profile.mean_recovery]
    pool = eligible if eligible else distant
    pool.sort(key=lambda entry: (-entry[0], -entry[1], entry[2]))

    proposals: list[Proposal] = []
    serial = 0
    for _, separation, mode_index, mode in pool:
        for alpha in config.nonlocal_alphas:
            mixed = _normalize(float(alpha) * seed + (1.0 - float(alpha)) * mode)
            proposals.append(
                Proposal(
                    operator=_DIRECTED_OPERATORS[0],
                    candidate_id=f"{seed_id}_nonlocal_{serial:03d}",
                    address=mixed,
                    source_ids=(seed_id, f"mode_{mode_index:03d}"),
                    metadata={
                        "alpha": float(alpha),
                        "mode_index": mode_index,
                        "null_index": serial,
                        "normalized_spectral_separation": float(separation),
                    },
                )
            )
            serial += 1
    return _fill_budget(_deduplicate(proposals), config.proposal_budget)


def propose_multi_seed_interference(
    seed_id: str,
    seed_addresses: Mapping[str, np.ndarray],
    config: TournamentConfig,
) -> tuple[Proposal, ...]:
    """Construct fixed-k signed equal-weight superpositions of admitted seeds."""

    if config.multi_seed_rule != "equal_weight_signed":
        raise ValueError("unsupported multi_seed_rule")
    if config.multi_seed_k < 2:
        raise ValueError("multi_seed_k must be at least two")
    if seed_id not in seed_addresses:
        raise ValueError(f"unknown seed_id: {seed_id}")
    others = sorted(identifier for identifier in seed_addresses if identifier != seed_id)
    choose_count = config.multi_seed_k - 1
    if len(others) < choose_count:
        raise ValueError("not enough admitted seeds for declared multi_seed_k")

    proposals: list[Proposal] = []
    serial = 0
    for partner_ids in combinations(others, choose_count):
        source_ids = (seed_id, *partner_ids)
        source_vectors = [_normalize(seed_addresses[identifier]) for identifier in source_ids]
        for sign_pattern in (1.0, -1.0):
            components = [source_vectors[0]]
            components.extend(
                sign_pattern * vector if index % 2 == 0 else vector
                for index, vector in enumerate(source_vectors[1:])
            )
            candidate = _normalize(np.sum(components, axis=0))
            maximum_seed_overlap = max(overlap(candidate, vector) for vector in source_vectors)
            if maximum_seed_overlap >= config.candidate_seed_overlap_limit:
                continue
            proposals.append(
                Proposal(
                    operator=_DIRECTED_OPERATORS[1],
                    candidate_id=f"{seed_id}_interference_{serial:03d}",
                    address=candidate,
                    source_ids=source_ids,
                    metadata={
                        "k": config.multi_seed_k,
                        "rule": config.multi_seed_rule,
                        "sign_pattern": sign_pattern,
                        "null_index": serial,
                        "max_input_overlap": float(maximum_seed_overlap),
                    },
                )
            )
            serial += 1
    return _fill_budget(_deduplicate(proposals), config.proposal_budget)


def _local_recovery_field(
    size: int,
    tuning_operators: tuple[np.ndarray, ...],
    tuning_levels: tuple[float, ...],
    thresholds: SurvivalThresholds,
    tuning_config: TuningConfig,
) -> np.ndarray:
    values: list[float] = []
    for node in range(size):
        address = np.zeros(size, dtype=float)
        address[node] = 1.0
        profile = evaluate_address(
            address, tuning_operators, tuning_levels, thresholds, tuning_config
        )
        persistence_scale = max(tuning_levels[-1], _EPS)
        values.append(
            0.75 * profile.mean_recovery
            + 0.25 * min(profile.persistence_time / persistence_scale, 1.0)
        )
    result = np.asarray(values, dtype=float)
    spread = float(np.ptp(result))
    return np.zeros_like(result) if spread <= _EPS else (result - np.min(result)) / spread


def propose_recovery_guided_walks(
    operator: np.ndarray,
    seed_id: str,
    seed_address: np.ndarray,
    tuning_operators: tuple[np.ndarray, ...],
    tuning_levels: tuple[float, ...],
    thresholds: SurvivalThresholds,
    tuning_config: TuningConfig,
    config: TournamentConfig,
) -> tuple[Proposal, ...]:
    """Walk deterministically toward high local recovery and extract support modes."""

    matrix = np.asarray(operator, dtype=float)
    size = matrix.shape[0]
    seed = _normalize(seed_address)
    adjacency = np.maximum(-(matrix - np.diag(np.diag(matrix))), 0.0)
    field = _local_recovery_field(
        size, tuning_operators, tuning_levels, thresholds, tuning_config
    )
    starts = np.argsort(-np.abs(seed))
    proposals: list[Proposal] = []
    serial = 0
    for start in starts:
        path = [int(start)]
        for _ in range(config.walk_length):
            current = path[-1]
            neighbors = np.flatnonzero(adjacency[current] > 0.0)
            if neighbors.size == 0:
                break
            edge_scale = max(float(np.max(adjacency[current, neighbors])), _EPS)
            ranked = []
            for neighbor in neighbors:
                recovery_term = float(field[int(neighbor)])
                edge_term = float(adjacency[current, neighbor] / edge_scale)
                score = (
                    config.walk_bias_recovery_weight * recovery_term
                    + (1.0 - config.walk_bias_recovery_weight) * edge_term
                )
                revisit = int(neighbor) in path
                ranked.append((revisit, -score, int(neighbor)))
            ranked.sort()
            path.append(ranked[0][2])

        for stop in range(1, len(path)):
            selected = set(path[: stop + 1])
            frontier = set(selected)
            for _ in range(config.walk_extraction_radius):
                frontier = {
                    int(neighbor)
                    for node in frontier
                    for neighbor in np.flatnonzero(adjacency[node] > 0.0)
                }
                selected.update(frontier)
            indices = np.asarray(sorted(selected), dtype=int)
            if indices.size < 2:
                continue
            restricted = matrix[np.ix_(indices, indices)]
            _, vectors = np.linalg.eigh(restricted)
            local_index = 1 if vectors.shape[1] > 1 else 0
            address = np.zeros(size, dtype=float)
            address[indices] = vectors[:, local_index]
            address = _sign_align(seed, address)
            proposals.append(
                Proposal(
                    operator=_DIRECTED_OPERATORS[2],
                    candidate_id=f"{seed_id}_walk_{serial:03d}",
                    address=address,
                    source_ids=(seed_id,),
                    metadata={
                        "start_node": int(start),
                        "stop_step": stop,
                        "null_index": serial,
                        "path": tuple(path[: stop + 1]),
                        "support": tuple(int(value) for value in indices),
                    },
                )
            )
            serial += 1
    return _fill_budget(_deduplicate(proposals), config.proposal_budget)


def _v0_decision(
    profile: RecoveryProfile,
    coordinate_null: RecoveryProfile,
    baseline: RecoveryProfile,
    tuning_config: TuningConfig,
) -> tuple[bool, tuple[str, ...]]:
    """Mirror the frozen V0 selection gate for arbitrary candidate addresses."""

    failures: list[str] = []
    if profile.mean_recovery < baseline.mean_recovery - tuning_config.max_mean_regression:
        failures.append("RECOVERABILITY_REGRESSION")
    if profile.persistence_time + _EPS < baseline.persistence_time:
        failures.append("PERSISTENCE_REGRESSION")
    candidate_k = len(profile.scores) if profile.k_star is None else profile.k_star
    baseline_k = len(baseline.scores) if baseline.k_star is None else baseline.k_star
    if candidate_k < baseline_k:
        failures.append("EARLIER_SUSTAINED_COLLAPSE")
    if profile.mean_recovery < coordinate_null.mean_recovery + tuning_config.control_margin:
        failures.append("NO_MATCHED_CONTROL_SEPARATION")
    improves = profile.mean_recovery >= baseline.mean_recovery + tuning_config.min_mean_gain
    stabilizes = _profile_non_regression(profile, baseline, tuning_config)
    if not (improves or stabilizes):
        failures.append("NO_RECOVERABLE_EXPANSION")
    return not failures, tuple(failures)


def _evaluate_proposal(
    proposal: Proposal,
    operator: np.ndarray,
    seed_id: str,
    regime: Regime,
    thresholds: SurvivalThresholds,
    tuning_config: TuningConfig,
    validation_config: ValidationConfig,
    tuning_baseline: RecoveryProfile,
    heldout_baseline: RecoveryProfile,
    spectrum_nulls: tuple[np.ndarray, ...],
) -> SurvivorLedgerRow:
    null_index = int(proposal.metadata.get("null_index", 0))
    tuning_profile = evaluate_address(
        proposal.address,
        regime.tuning_operators,
        regime.tuning_levels,
        thresholds,
        tuning_config,
    )
    tuning_coordinate = evaluate_address(
        _coordinate_null(proposal.address, null_index),
        regime.tuning_operators,
        regime.tuning_levels,
        thresholds,
        tuning_config,
    )
    tuning_pass, tuning_failures = _v0_decision(
        tuning_profile, tuning_coordinate, tuning_baseline, tuning_config
    )
    heldout_profile = evaluate_address(
        proposal.address,
        regime.heldout_operators,
        regime.heldout_levels,
        thresholds,
        tuning_config,
    )
    coordinate_profile = evaluate_address(
        _coordinate_null(proposal.address, null_index),
        regime.heldout_operators,
        regime.heldout_levels,
        thresholds,
        tuning_config,
    )
    heldout_pass, heldout_failures = _v0_decision(
        heldout_profile, coordinate_profile, heldout_baseline, tuning_config
    )
    degree_profile = evaluate_address(
        _degree_matched_null(operator, proposal.address),
        regime.heldout_operators,
        regime.heldout_levels,
        thresholds,
        tuning_config,
    )
    spectrum_profile = evaluate_address(
        proposal.address,
        spectrum_nulls,
        regime.heldout_levels,
        thresholds,
        tuning_config,
    )
    amplification, recovery_gain, persistence_gain, k_displacement = _amplifies(
        heldout_profile, heldout_baseline, validation_config
    )
    margins = (
        heldout_profile.mean_recovery - coordinate_profile.mean_recovery,
        heldout_profile.mean_recovery - degree_profile.mean_recovery,
        heldout_profile.mean_recovery - spectrum_profile.mean_recovery,
    )
    hostile_pass = min(margins) >= validation_config.hostile_control_margin
    notes = list(tuning_failures)
    notes.extend(f"HELDOUT_{failure}" for failure in heldout_failures)
    if not hostile_pass:
        notes.append("HOSTILE_NULL_SEPARATION_FAILED")
    if not amplification:
        notes.append("NO_HELDOUT_AMPLIFICATION")
    accepted = tuning_pass and heldout_pass and hostile_pass and amplification
    return SurvivorLedgerRow(
        operator=proposal.operator,
        seed_id=seed_id,
        regime=regime.name,
        candidate_id=proposal.candidate_id,
        tuning_decision="PASS" if tuning_pass else "FAIL",
        heldout_decision="PASS" if accepted else "FAIL",
        delta_recovery=float(recovery_gain),
        delta_persistence=float(persistence_gain),
        k_star_displacement=int(k_displacement),
        coordinate_null_margin=float(margins[0]),
        degree_null_margin=float(margins[1]),
        spectrum_null_margin=float(margins[2]),
        baseline_match="not_evaluated",
        unique_flag=False,
        accepted=accepted,
        notes=tuple(dict.fromkeys(notes)),
        address=proposal.address.copy(),
    )


def _address_distance(operator: np.ndarray, left: np.ndarray, right: np.ndarray) -> float:
    eigenvalues = np.linalg.eigvalsh(np.asarray(operator, dtype=float))
    span = max(float(np.ptp(eigenvalues)), _EPS)
    return abs(_rayleigh(operator, left) - _rayleigh(operator, right)) / span


def _assign_uniqueness(
    rows: list[SurvivorLedgerRow],
    operator: np.ndarray,
    config: TournamentConfig,
) -> list[SurvivorLedgerRow]:
    output: list[SurvivorLedgerRow] = []
    groups = {(row.seed_id, row.regime) for row in rows}
    baseline_survivors: dict[tuple[str, str], list[SurvivorLedgerRow]] = {}
    for key in groups:
        baseline_survivors[key] = [
            row
            for row in rows
            if (row.seed_id, row.regime) == key
            and row.operator in _BASELINES
            and row.accepted
        ]
    for row in rows:
        if row.operator in _BASELINES:
            output.append(
                SurvivorLedgerRow(
                    **{
                        **row.__dict__,
                        "baseline_match": "self",
                    }
                )
            )
            continue
        if not row.accepted:
            output.append(
                SurvivorLedgerRow(
                    **{
                        **row.__dict__,
                        "baseline_match": "not_applicable",
                    }
                )
            )
            continue
        matches: set[str] = set()
        for baseline in baseline_survivors[(row.seed_id, row.regime)]:
            distance = _address_distance(operator, row.address, baseline.address)
            address_overlap = overlap(row.address, baseline.address)
            if (
                distance <= config.address_distance_threshold
                or address_overlap >= config.subspace_overlap_threshold
            ):
                matches.add(baseline.operator)
        unique = not matches
        match_label = "none" if unique else "+".join(sorted(matches))
        notes = row.notes if unique else (*row.notes, "BASELINE_RECOVERABLE")
        output.append(
            SurvivorLedgerRow(
                **{
                    **row.__dict__,
                    "baseline_match": match_label,
                    "unique_flag": unique,
                    "notes": notes,
                }
            )
        )
    return output


def _local_baseline_proposals(
    operator: np.ndarray,
    seed_id: str,
    seed: np.ndarray,
    tuning_config: TuningConfig,
    budget: int,
) -> tuple[Proposal, ...]:
    local_config = TuningConfig(
        **{**tuning_config.__dict__, "max_candidates": budget}
    )
    candidates = propose_spectral_candidates(operator, seed, local_config)
    proposals = tuple(
        Proposal(
            operator=_BASELINES[0],
            candidate_id=f"{seed_id}_v0_{candidate.candidate_id}",
            address=candidate.address,
            source_ids=(seed_id,),
            metadata={
                "source_candidate_id": candidate.candidate_id,
                "null_index": int(candidate.candidate_id.rsplit("_", 1)[1]),
            },
        )
        for candidate in candidates
    )
    return _fill_budget(_deduplicate(proposals), budget)


def _random_baseline_proposals(
    operator: np.ndarray,
    seed_id: str,
    budget: int,
    random_seed: int,
    excluded_indices: frozenset[int],
) -> tuple[Proposal, ...]:
    modes = _random_eigenmodes(
        operator,
        budget,
        random_seed,
        excluded_indices=excluded_indices,
    )
    proposals = tuple(
        Proposal(
            operator=_BASELINES[1],
            candidate_id=f"{seed_id}_random_mode_{index:03d}",
            address=address,
            source_ids=(seed_id,),
            metadata={"mode_index": index, "null_index": index},
        )
        for index, address in modes
    )
    return _fill_budget(proposals, budget)


def run_proposal_tournament(
    operator: np.ndarray,
    seed_addresses: Mapping[str, np.ndarray],
    regimes: Iterable[Regime],
    thresholds: SurvivalThresholds,
    tuning_config: TuningConfig,
    validation_config: ValidationConfig,
    config: TournamentConfig = TournamentConfig(),
) -> TournamentAudit:
    """Run all proposal operators and baselines under the frozen V0.2 judge."""

    matrix = np.asarray(operator, dtype=float)
    seeds = {identifier: _normalize(address) for identifier, address in seed_addresses.items()}
    regime_list = tuple(regimes)
    if len(seeds) < 3:
        raise ValueError("V0.3 requires at least three admitted seed addresses")
    if config.hard_regime_name not in {regime.name for regime in regime_list}:
        raise ValueError("declared hard regime is missing")
    if len({regime.name for regime in regime_list}) != len(regime_list):
        raise ValueError("regime names must be unique")

    rows: list[SurvivorLedgerRow] = []
    for regime_index, regime in enumerate(regime_list):
        spectrum_nulls = _spectrum_preserving_nulls(regime.heldout_operators)
        for seed_index, (seed_id, seed) in enumerate(sorted(seeds.items())):
            tuning_baseline = evaluate_address(
                seed,
                regime.tuning_operators,
                regime.tuning_levels,
                thresholds,
                tuning_config,
            )
            heldout_baseline = evaluate_address(
                seed,
                regime.heldout_operators,
                regime.heldout_levels,
                thresholds,
                tuning_config,
            )
            local = _local_baseline_proposals(
                matrix, seed_id, seed, tuning_config, config.proposal_budget
            )
            local_indices = frozenset(
                int(proposal.metadata["source_candidate_id"].rsplit("_", 1)[1])
                for proposal in local
            )
            random = _random_baseline_proposals(
                matrix,
                seed_id,
                config.proposal_budget,
                config.random_seed + 1009 * regime_index + seed_index,
                local_indices,
            )
            proposal_sets = (
                local,
                random,
                propose_non_local_spectral_mixes(
                    matrix,
                    seed_id,
                    seed,
                    regime.tuning_operators,
                    regime.tuning_levels,
                    thresholds,
                    tuning_config,
                    config,
                ),
                propose_multi_seed_interference(seed_id, seeds, config),
                propose_recovery_guided_walks(
                    matrix,
                    seed_id,
                    seed,
                    regime.tuning_operators,
                    regime.tuning_levels,
                    thresholds,
                    tuning_config,
                    config,
                ),
            )
            for proposals in proposal_sets:
                if len(proposals) != config.proposal_budget:
                    raise RuntimeError("proposal-budget parity violated")
                rows.extend(
                    _evaluate_proposal(
                        proposal,
                        matrix,
                        seed_id,
                        regime,
                        thresholds,
                        tuning_config,
                        validation_config,
                        tuning_baseline,
                        heldout_baseline,
                        spectrum_nulls,
                    )
                    for proposal in proposals
                )

    rows = _assign_uniqueness(rows, matrix, config)
    aggregate: dict[str, object] = {}
    directed_positive = False
    for operator_name in (*_BASELINES, *_DIRECTED_OPERATORS):
        operator_rows = [row for row in rows if row.operator == operator_name]
        survivors = [row for row in operator_rows if row.accepted]
        unique_survivors = [row for row in survivors if row.unique_flag]
        unique_seed_ids = sorted({row.seed_id for row in unique_survivors})
        hard_unique = [
            row for row in unique_survivors if row.regime == config.hard_regime_name
        ]
        if (
            operator_name in _DIRECTED_OPERATORS
            and len(unique_seed_ids) >= config.minimum_positive_seeds
            and bool(hard_unique)
        ):
            directed_positive = True
        aggregate[operator_name] = {
            "proposal_count": len(operator_rows),
            "survivor_count": len(survivors),
            "yield": float(len(survivors) / max(len(operator_rows), 1)),
            "proposals_per_survivor": (
                None if not survivors else float(len(operator_rows) / len(survivors))
            ),
            "unique_survivor_count": len(unique_survivors),
            "unique_seed_count": len(unique_seed_ids),
            "unique_survivors_by_regime": {
                regime.name: sum(
                    int(row.regime == regime.name) for row in unique_survivors
                )
                for regime in regime_list
            },
        }
    aggregate["status"] = (
        "PASS_BOUNDED_UNIQUE_DIRECTED_ADVANTAGE"
        if directed_positive
        else "OPEN_NO_UNIQUE_DIRECTED_ADVANTAGE"
    )
    aggregate["judge"] = "V0.2_FROZEN"

    contract = {
        "version": "HAOS_GEN_V0_3",
        "judge": "V0.2 byte-for-byte reference implementation",
        "proposal_budget_per_operator_per_seed": config.proposal_budget,
        "seed_count": len(seeds),
        "regimes": [regime.name for regime in regime_list],
        "operators": list(_DIRECTED_OPERATORS),
        "baselines": list(_BASELINES),
        "uniqueness": {
            "normalized_rayleigh_distance_min": config.address_distance_threshold,
            "maximum_absolute_overlap": config.subspace_overlap_threshold,
            "both_conditions_required": True,
        },
        "success_requires_multiple_seeds": config.minimum_positive_seeds,
        "success_requires_hard_regime": config.hard_regime_name,
        "mutates_v0_2_judge": False,
    }
    return TournamentAudit(rows=tuple(rows), aggregate=aggregate, contract=contract)
