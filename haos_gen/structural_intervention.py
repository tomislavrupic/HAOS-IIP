"""HAOS-GEN V0.4 Class 1: bounded local edge-weight interventions."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Iterable, Mapping

import numpy as np

from telemetry.frozen_metrics import SurvivalThresholds, overlap

from .hostile_validation import (
    ValidationConfig,
    _spectrum_preserving_nulls,
    validate_expansion_hostile,
)
from .proposal_tournament import Proposal, Regime, _evaluate_proposal
from .spectral_tuning import (
    RecoveryProfile,
    TuningConfig,
    evaluate_address,
    propose_spectral_candidates,
)


_EPS = 1.0e-12


@dataclass(frozen=True)
class EdgeInterventionConfig:
    edge_budget: int = 2
    max_weight_change: float = 0.08
    retention_tolerance: float = 0.02
    reversal_tolerance: float = 1.0e-10
    subspace_equivalence_threshold: float = 0.92
    random_trials: int = 4
    random_seed: int = 20260724
    minimum_expansion_substrates: int = 2
    hard_regime_name: str = "hard"


@dataclass(frozen=True)
class EdgeUpdate:
    src: int
    dst: int
    delta_weight: float
    tuning_score: float


@dataclass(frozen=True)
class DiscoveredAddress:
    address_id: str
    seed_id: str
    candidate_id: str
    mode_index: int
    address: np.ndarray = field(repr=False, compare=False)
    heldout_profile: RecoveryProfile


@dataclass(frozen=True)
class InterventionLedgerRow:
    substrate_id: str
    regime: str
    intervention_type: str
    edges_touched: int
    total_l1_change: float
    pre_address_count: int
    post_address_count: int
    existing_retention_ok: bool
    new_addresses_count: int
    new_unique_heldout_count: int
    subspace_equivalence_flag: bool
    delta_mean_recovery: float
    delta_persistence: float
    k_star_shift: float
    random_baseline_new_count: int
    random_baseline_delta_recovery: float
    reversal_ok: bool
    outcome_label: str
    notes: tuple[str, ...]
    updates: tuple[EdgeUpdate, ...]

    def to_jsonable(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class StructuralInterventionAudit:
    ledger: tuple[InterventionLedgerRow, ...]
    aggregate: dict[str, object]
    contract: dict[str, object]

    def to_jsonable(self) -> dict[str, object]:
        return {
            "ledger": [row.to_jsonable() for row in self.ledger],
            "aggregate": self.aggregate,
            "contract": self.contract,
        }


def _normalize(vector: np.ndarray) -> np.ndarray:
    values = np.asarray(vector, dtype=float)
    norm = float(np.linalg.norm(values))
    if norm <= _EPS:
        raise ValueError("address must have non-zero norm")
    return values / norm


def _edge_list(operator: np.ndarray) -> tuple[tuple[int, int, float], ...]:
    matrix = np.asarray(operator, dtype=float)
    edges: list[tuple[int, int, float]] = []
    for src in range(matrix.shape[0]):
        for dst in range(src + 1, matrix.shape[0]):
            weight = -float(matrix[src, dst])
            if weight > _EPS:
                edges.append((src, dst, weight))
    return tuple(edges)


def _edge_laplacian(size: int, src: int, dst: int, delta_weight: float) -> np.ndarray:
    incidence = np.zeros(size, dtype=float)
    incidence[src] = 1.0
    incidence[dst] = -1.0
    return float(delta_weight) * np.outer(incidence, incidence)


def _update_matrix(size: int, updates: Iterable[EdgeUpdate]) -> np.ndarray:
    result = np.zeros((size, size), dtype=float)
    for update in updates:
        result += _edge_laplacian(
            size, update.src, update.dst, update.delta_weight
        )
    return result


def _shift_family(
    operators: tuple[np.ndarray, ...],
    update_matrix: np.ndarray,
) -> tuple[np.ndarray, ...]:
    return tuple(np.asarray(operator, dtype=float) + update_matrix for operator in operators)


def _mean_seed_objective(
    seeds: Mapping[str, np.ndarray],
    operators: tuple[np.ndarray, ...],
    levels: tuple[float, ...],
    thresholds: SurvivalThresholds,
    tuning_config: TuningConfig,
) -> float:
    profiles = [
        evaluate_address(address, operators, levels, thresholds, tuning_config)
        for address in seeds.values()
    ]
    level_scale = max(levels[-1], _EPS)
    return float(
        np.mean(
            [
                profile.mean_recovery
                + 0.25 * min(profile.persistence_time / level_scale, 1.0)
                for profile in profiles
            ]
        )
    )


def select_directed_edge_updates(
    operator: np.ndarray,
    admitted_seeds: Mapping[str, np.ndarray],
    tuning_operators: tuple[np.ndarray, ...],
    tuning_levels: tuple[float, ...],
    thresholds: SurvivalThresholds,
    tuning_config: TuningConfig,
    config: EdgeInterventionConfig,
) -> tuple[EdgeUpdate, ...]:
    """Select signed edge changes using tuning-family recovery only."""

    if config.edge_budget < 1:
        raise ValueError("edge_budget must be at least one")
    if config.max_weight_change <= 0.0:
        raise ValueError("max_weight_change must be positive")
    edges = _edge_list(operator)
    if len(edges) < config.edge_budget:
        raise ValueError("edge_budget exceeds available graph edges")
    normalized_seeds = {
        identifier: _normalize(address) for identifier, address in admitted_seeds.items()
    }
    baseline = _mean_seed_objective(
        normalized_seeds,
        tuning_operators,
        tuning_levels,
        thresholds,
        tuning_config,
    )
    scored: list[tuple[float, int, int, float]] = []
    size = np.asarray(operator).shape[0]
    for src, dst, weight in edges:
        signed_trials = (config.max_weight_change, -min(config.max_weight_change, weight))
        trial_results: list[tuple[float, float]] = []
        for delta in signed_trials:
            trial_update = _edge_laplacian(size, src, dst, delta)
            score = _mean_seed_objective(
                normalized_seeds,
                _shift_family(tuning_operators, trial_update),
                tuning_levels,
                thresholds,
                tuning_config,
            )
            trial_results.append((score - baseline, delta))
        improvement, delta = max(trial_results, key=lambda item: (item[0], item[1]))
        scored.append((float(improvement), src, dst, float(delta)))
    scored.sort(key=lambda item: (-item[0], item[1], item[2]))
    return tuple(
        EdgeUpdate(src=src, dst=dst, delta_weight=delta, tuning_score=score)
        for score, src, dst, delta in scored[: config.edge_budget]
    )


def _matched_random_updates(
    operator: np.ndarray,
    directed: tuple[EdgeUpdate, ...],
    trials: int,
    seed: int,
) -> tuple[tuple[EdgeUpdate, ...], ...]:
    edges = _edge_list(operator)
    directed_pairs = {(update.src, update.dst) for update in directed}
    pool = [(src, dst, weight) for src, dst, weight in edges if (src, dst) not in directed_pairs]
    if len(pool) < len(directed):
        pool = list(edges)
    generator = np.random.default_rng(seed)
    magnitudes = tuple(update.delta_weight for update in directed)
    results: list[tuple[EdgeUpdate, ...]] = []
    for trial in range(trials):
        selected_indices = generator.choice(
            len(pool), size=len(directed), replace=False
        )
        permuted = generator.permutation(np.asarray(magnitudes, dtype=float))
        updates: list[EdgeUpdate] = []
        for position, edge_index in enumerate(selected_indices):
            src, dst, weight = pool[int(edge_index)]
            delta = float(permuted[position])
            if delta < 0.0 and abs(delta) > weight + _EPS:
                raise ValueError("matched suppression exceeds random edge weight")
            updates.append(
                EdgeUpdate(
                    src=src,
                    dst=dst,
                    delta_weight=delta,
                    tuning_score=float("nan"),
                )
            )
        results.append(tuple(updates))
    return tuple(results)


def _deduplicate_addresses(
    addresses: Iterable[DiscoveredAddress],
    threshold: float,
) -> tuple[DiscoveredAddress, ...]:
    retained: list[DiscoveredAddress] = []
    for address in addresses:
        if any(overlap(address.address, prior.address) >= threshold for prior in retained):
            continue
        retained.append(address)
    return tuple(retained)


def _discover_addresses(
    operator: np.ndarray,
    seeds: Mapping[str, np.ndarray],
    regime: Regime,
    thresholds: SurvivalThresholds,
    tuning_config: TuningConfig,
    validation_config: ValidationConfig,
    equivalence_threshold: float,
) -> tuple[DiscoveredAddress, ...]:
    discovered: list[DiscoveredAddress] = []
    for seed_id, seed in sorted(seeds.items()):
        audit = validate_expansion_hostile(
            operator,
            seed,
            regime.tuning_operators,
            regime.tuning_levels,
            regime.heldout_operators,
            regime.heldout_levels,
            thresholds,
            tuning_config,
            validation_config,
        )
        proposals = {
            proposal.candidate_id: proposal
            for proposal in propose_spectral_candidates(operator, seed, tuning_config)
        }
        for candidate in audit.candidates:
            if not candidate.accepted:
                continue
            proposal = proposals[candidate.candidate_id]
            mode_index = int(candidate.candidate_id.rsplit("_", 1)[1])
            discovered.append(
                DiscoveredAddress(
                    address_id=f"{seed_id}:{candidate.candidate_id}",
                    seed_id=seed_id,
                    candidate_id=candidate.candidate_id,
                    mode_index=mode_index,
                    address=proposal.address.copy(),
                    heldout_profile=candidate.heldout_profile,
                )
            )
    return _deduplicate_addresses(discovered, equivalence_threshold)


def _profile_summary(addresses: tuple[DiscoveredAddress, ...]) -> tuple[float, float, float]:
    if not addresses:
        return 0.0, 0.0, 0.0
    mean_recovery = float(
        np.mean([address.heldout_profile.mean_recovery for address in addresses])
    )
    mean_persistence = float(
        np.mean([address.heldout_profile.persistence_time for address in addresses])
    )
    mean_k = float(
        np.mean(
            [
                len(address.heldout_profile.scores)
                if address.heldout_profile.k_star is None
                else address.heldout_profile.k_star
                for address in addresses
            ]
        )
    )
    return mean_recovery, mean_persistence, mean_k


def _new_unique_addresses(
    pre: tuple[DiscoveredAddress, ...],
    post: tuple[DiscoveredAddress, ...],
    threshold: float,
) -> tuple[DiscoveredAddress, ...]:
    return tuple(
        candidate
        for candidate in post
        if all(overlap(candidate.address, prior.address) < threshold for prior in pre)
    )


def _retention_check(
    operator: np.ndarray,
    pre_addresses: tuple[DiscoveredAddress, ...],
    seeds: Mapping[str, np.ndarray],
    pre_regime: Regime,
    post_regime: Regime,
    thresholds: SurvivalThresholds,
    tuning_config: TuningConfig,
    validation_config: ValidationConfig,
    tolerance: float,
) -> tuple[bool, tuple[str, ...]]:
    spectrum_nulls = _spectrum_preserving_nulls(post_regime.heldout_operators)
    notes: list[str] = []
    for seed_id, seed in sorted(seeds.items()):
        pre_seed_profile = evaluate_address(
            seed,
            pre_regime.heldout_operators,
            pre_regime.heldout_levels,
            thresholds,
            tuning_config,
        )
        post_seed_profile = evaluate_address(
            seed,
            post_regime.heldout_operators,
            post_regime.heldout_levels,
            thresholds,
            tuning_config,
        )
        if (
            post_seed_profile.mean_recovery
            < pre_seed_profile.mean_recovery - tolerance
            or post_seed_profile.persistence_time
            < pre_seed_profile.persistence_time - tolerance
        ):
            notes.append(f"ADMITTED_SEED_REGRESSION:{seed_id}")
    for address in pre_addresses:
        seed = seeds[address.seed_id]
        tuning_baseline = evaluate_address(
            seed,
            post_regime.tuning_operators,
            post_regime.tuning_levels,
            thresholds,
            tuning_config,
        )
        heldout_baseline = evaluate_address(
            seed,
            post_regime.heldout_operators,
            post_regime.heldout_levels,
            thresholds,
            tuning_config,
        )
        proposal = Proposal(
            operator="v0_4_retention_probe",
            candidate_id=address.address_id,
            address=address.address,
            source_ids=(address.seed_id,),
            metadata={"null_index": address.mode_index},
        )
        row = _evaluate_proposal(
            proposal,
            operator,
            address.seed_id,
            post_regime,
            thresholds,
            tuning_config,
            validation_config,
            tuning_baseline,
            heldout_baseline,
            spectrum_nulls,
        )
        pre_profile = evaluate_address(
            address.address,
            pre_regime.heldout_operators,
            pre_regime.heldout_levels,
            thresholds,
            tuning_config,
        )
        post_profile = evaluate_address(
            address.address,
            post_regime.heldout_operators,
            post_regime.heldout_levels,
            thresholds,
            tuning_config,
        )
        regressed = (
            post_profile.mean_recovery < pre_profile.mean_recovery - tolerance
            or post_profile.persistence_time < pre_profile.persistence_time - tolerance
        )
        if not row.accepted or regressed:
            notes.append(f"RETENTION_FAILED:{address.address_id}")
    return not notes, tuple(notes)


@dataclass(frozen=True)
class _Evaluation:
    pre: tuple[DiscoveredAddress, ...]
    post: tuple[DiscoveredAddress, ...]
    raw_new_count: int
    new_unique: tuple[DiscoveredAddress, ...]
    retention_ok: bool
    delta_recovery: float
    delta_persistence: float
    k_star_shift: float
    notes: tuple[str, ...]


def _evaluate_update(
    operator: np.ndarray,
    seeds: Mapping[str, np.ndarray],
    regime: Regime,
    updates: tuple[EdgeUpdate, ...],
    thresholds: SurvivalThresholds,
    tuning_config: TuningConfig,
    validation_config: ValidationConfig,
    config: EdgeInterventionConfig,
    pre: tuple[DiscoveredAddress, ...] | None = None,
) -> _Evaluation:
    update = _update_matrix(operator.shape[0], updates)
    modified_operator = np.asarray(operator, dtype=float) + update
    modified_regime = Regime(
        name=regime.name,
        tuning_operators=_shift_family(regime.tuning_operators, update),
        tuning_levels=regime.tuning_levels,
        heldout_operators=_shift_family(regime.heldout_operators, update),
        heldout_levels=regime.heldout_levels,
    )
    pre_addresses = pre or _discover_addresses(
        operator,
        seeds,
        regime,
        thresholds,
        tuning_config,
        validation_config,
        config.subspace_equivalence_threshold,
    )
    post_addresses = _discover_addresses(
        modified_operator,
        seeds,
        modified_regime,
        thresholds,
        tuning_config,
        validation_config,
        config.subspace_equivalence_threshold,
    )
    new_unique = _new_unique_addresses(
        pre_addresses, post_addresses, config.subspace_equivalence_threshold
    )
    pre_ids = {address.address_id for address in pre_addresses}
    raw_new_count = sum(
        int(address.address_id not in pre_ids) for address in post_addresses
    )
    retention_ok, notes = _retention_check(
        modified_operator,
        pre_addresses,
        seeds,
        regime,
        modified_regime,
        thresholds,
        tuning_config,
        validation_config,
        config.retention_tolerance,
    )
    pre_recovery, pre_persistence, pre_k = _profile_summary(pre_addresses)
    post_recovery, post_persistence, post_k = _profile_summary(post_addresses)
    return _Evaluation(
        pre=pre_addresses,
        post=post_addresses,
        raw_new_count=raw_new_count,
        new_unique=new_unique,
        retention_ok=retention_ok,
        delta_recovery=float(post_recovery - pre_recovery),
        delta_persistence=float(post_persistence - pre_persistence),
        k_star_shift=float(post_k - pre_k),
        notes=notes,
    )


def _reversal_ok(
    operator: np.ndarray,
    seeds: Mapping[str, np.ndarray],
    regime: Regime,
    updates: tuple[EdgeUpdate, ...],
    pre_addresses: tuple[DiscoveredAddress, ...],
    thresholds: SurvivalThresholds,
    tuning_config: TuningConfig,
    validation_config: ValidationConfig,
    equivalence_threshold: float,
    tolerance: float,
) -> bool:
    update = _update_matrix(operator.shape[0], updates)
    modified_operator = np.asarray(operator, dtype=float) + update
    reversed_operator = modified_operator - update
    if not np.allclose(reversed_operator, operator, atol=tolerance, rtol=0.0):
        return False
    for family in (regime.tuning_operators, regime.heldout_operators):
        for original in family:
            reversed_family_member = (np.asarray(original) + update) - update
            if not np.allclose(
                reversed_family_member, original, atol=tolerance, rtol=0.0
            ):
                return False
    reversed_regime = Regime(
        name=regime.name,
        tuning_operators=tuple(
            (np.asarray(member) + update) - update
            for member in regime.tuning_operators
        ),
        tuning_levels=regime.tuning_levels,
        heldout_operators=tuple(
            (np.asarray(member) + update) - update
            for member in regime.heldout_operators
        ),
        heldout_levels=regime.heldout_levels,
    )
    reversed_addresses = _discover_addresses(
        reversed_operator,
        seeds,
        reversed_regime,
        thresholds,
        tuning_config,
        validation_config,
        equivalence_threshold,
    )
    if len(reversed_addresses) != len(pre_addresses):
        return False
    pre_summary = _profile_summary(pre_addresses)
    reversed_summary = _profile_summary(reversed_addresses)
    return all(
        abs(left - right) <= tolerance
        for left, right in zip(pre_summary, reversed_summary)
    )


def _classify(
    directed: _Evaluation,
    random_evaluations: tuple[_Evaluation, ...],
    reversal_ok: bool,
) -> tuple[str, tuple[str, ...]]:
    notes = list(directed.notes)
    max_random_new = max(
        (
            len(evaluation.new_unique)
            for evaluation in random_evaluations
            if evaluation.retention_ok and len(evaluation.post) >= len(evaluation.pre)
        ),
        default=0,
    )
    max_random_recovery = max(
        (
            evaluation.delta_recovery
            for evaluation in random_evaluations
            if evaluation.retention_ok and len(evaluation.post) >= len(evaluation.pre)
        ),
        default=0.0,
    )
    max_random_persistence = max(
        (
            evaluation.delta_persistence
            for evaluation in random_evaluations
            if evaluation.retention_ok and len(evaluation.post) >= len(evaluation.pre)
        ),
        default=0.0,
    )
    max_random_k_shift = max(
        (
            evaluation.k_star_shift
            for evaluation in random_evaluations
            if evaluation.retention_ok and len(evaluation.post) >= len(evaluation.pre)
        ),
        default=0.0,
    )
    if not directed.retention_ok or len(directed.post) < len(directed.pre):
        return "DEGRADATION", tuple((*notes, "EXISTING_SET_REGRESSED"))
    if not reversal_ok:
        return "RELABELING / NULL", tuple((*notes, "REVERSAL_FAILED"))
    if len(directed.new_unique) > 0:
        if len(directed.new_unique) <= max_random_new:
            return "RELABELING / NULL", tuple(
                (*notes, "MATCHED_OR_EXCEEDED_BY_RANDOM")
            )
        return "STRUCTURAL_EXPANSION", tuple(notes)
    strict_stabilization = (
        directed.delta_recovery > max(0.0, max_random_recovery) + _EPS
        or directed.delta_persistence > max(0.0, max_random_persistence) + _EPS
        or directed.k_star_shift > max(0.0, max_random_k_shift) + _EPS
    )
    if strict_stabilization:
        return "STABILIZATION", tuple(notes)
    return "RELABELING / NULL", tuple((*notes, "NO_NET_STRUCTURAL_GAIN"))


def run_edge_intervention_experiment(
    substrates: Mapping[str, tuple[np.ndarray, Mapping[str, np.ndarray]]],
    regimes_by_substrate: Mapping[str, tuple[Regime, ...]],
    thresholds: SurvivalThresholds,
    tuning_config: TuningConfig,
    validation_config: ValidationConfig,
    config: EdgeInterventionConfig = EdgeInterventionConfig(),
) -> StructuralInterventionAudit:
    """Run the fixed Class-1 intervention matrix across substrates and regimes."""

    if len(substrates) < 3:
        raise ValueError("V0.4 requires at least three distinct substrates")
    ledger: list[InterventionLedgerRow] = []
    valid_random_expansion_trials = 0
    for substrate_index, (substrate_id, payload) in enumerate(sorted(substrates.items())):
        operator, seed_mapping = payload
        seeds = {identifier: _normalize(address) for identifier, address in seed_mapping.items()}
        regimes = regimes_by_substrate.get(substrate_id, ())
        if config.hard_regime_name not in {regime.name for regime in regimes}:
            raise ValueError(f"{substrate_id} is missing the declared hard regime")
        directed_updates = select_directed_edge_updates(
            operator,
            seeds,
            regimes[0].tuning_operators,
            regimes[0].tuning_levels,
            thresholds,
            tuning_config,
            config,
        )
        random_updates = _matched_random_updates(
            operator,
            directed_updates,
            config.random_trials,
            config.random_seed + substrate_index,
        )
        for regime in regimes:
            pre = _discover_addresses(
                operator,
                seeds,
                regime,
                thresholds,
                tuning_config,
                validation_config,
                config.subspace_equivalence_threshold,
            )
            directed = _evaluate_update(
                operator,
                seeds,
                regime,
                directed_updates,
                thresholds,
                tuning_config,
                validation_config,
                config,
                pre=pre,
            )
            random_evaluations = tuple(
                _evaluate_update(
                    operator,
                    seeds,
                    regime,
                    updates,
                    thresholds,
                    tuning_config,
                    validation_config,
                    config,
                    pre=pre,
                )
                for updates in random_updates
            )
            valid_random_expansion_trials += sum(
                int(
                    evaluation.retention_ok
                    and len(evaluation.post) >= len(evaluation.pre)
                    and len(evaluation.new_unique) > 0
                )
                for evaluation in random_evaluations
            )
            reversed_ok = _reversal_ok(
                operator,
                seeds,
                regime,
                directed_updates,
                pre,
                thresholds,
                tuning_config,
                validation_config,
                config.subspace_equivalence_threshold,
                config.reversal_tolerance,
            )
            outcome, notes = _classify(directed, random_evaluations, reversed_ok)
            random_new = max(
                (
                    len(value.new_unique)
                    for value in random_evaluations
                    if value.retention_ok and len(value.post) >= len(value.pre)
                ),
                default=0,
            )
            random_recovery = max(
                (
                    value.delta_recovery
                    for value in random_evaluations
                    if value.retention_ok and len(value.post) >= len(value.pre)
                ),
                default=0.0,
            )
            ledger.append(
                InterventionLedgerRow(
                    substrate_id=substrate_id,
                    regime=regime.name,
                    intervention_type="local_edge_weight_v0_4",
                    edges_touched=len(directed_updates),
                    total_l1_change=float(
                        sum(abs(update.delta_weight) for update in directed_updates)
                    ),
                    pre_address_count=len(directed.pre),
                    post_address_count=len(directed.post),
                    existing_retention_ok=directed.retention_ok,
                    new_addresses_count=directed.raw_new_count,
                    new_unique_heldout_count=len(directed.new_unique),
                    subspace_equivalence_flag=bool(
                        directed.raw_new_count > len(directed.new_unique)
                    ),
                    delta_mean_recovery=directed.delta_recovery,
                    delta_persistence=directed.delta_persistence,
                    k_star_shift=directed.k_star_shift,
                    random_baseline_new_count=random_new,
                    random_baseline_delta_recovery=float(random_recovery),
                    reversal_ok=reversed_ok,
                    outcome_label=outcome,
                    notes=notes,
                    updates=directed_updates,
                )
            )

    expansions = [row for row in ledger if row.outcome_label == "STRUCTURAL_EXPANSION"]
    expansion_substrates = {row.substrate_id for row in expansions}
    hard_expansion_substrates = {
        row.substrate_id
        for row in expansions
        if row.regime == config.hard_regime_name
    }
    confirmed = (
        len(expansion_substrates) >= config.minimum_expansion_substrates
        and len(hard_expansion_substrates) >= config.minimum_expansion_substrates
    )
    random_trial_total = len(ledger) * config.random_trials
    aggregate = {
        "directed_expansion_count": len(expansions),
        "directed_expansion_rate": float(len(expansions) / max(len(ledger), 1)),
        "random_expansion_trial_count": valid_random_expansion_trials,
        "random_trial_total": random_trial_total,
        "random_expansion_rate": float(
            valid_random_expansion_trials / max(random_trial_total, 1)
        ),
        "expansion_substrate_count": len(expansion_substrates),
        "hard_regime_expansion_substrate_count": len(hard_expansion_substrates),
        "status": (
            "STRUCTURAL_EXPANSION_CONFIRMED"
            if confirmed
            else "OPEN_NO_STRUCTURAL_INTERVENTION_ADVANTAGE"
        ),
        "judge": "V0.2_FROZEN",
    }
    contract = {
        "version": "HAOS_GEN_V0_4_CLASS_1",
        "intervention": "local_edge_weight",
        "edge_budget": config.edge_budget,
        "max_weight_change": config.max_weight_change,
        "max_total_l1_change": config.edge_budget * config.max_weight_change,
        "selector_uses_heldout": False,
        "matched_random_trials": config.random_trials,
        "random_uses_identical_signed_magnitude_multiset": True,
        "subspace_equivalence_threshold": config.subspace_equivalence_threshold,
        "retention_tolerance": config.retention_tolerance,
        "reversal_tolerance": config.reversal_tolerance,
        "substrate_count": len(substrates),
        "mutates_frozen_judge": False,
    }
    return StructuralInterventionAudit(
        ledger=tuple(ledger), aggregate=aggregate, contract=contract
    )
