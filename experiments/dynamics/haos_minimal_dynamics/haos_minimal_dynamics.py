"""Stripped-back HAOS relational-address dynamics probe."""

from __future__ import annotations

from dataclasses import dataclass
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parent.parent.parent
EPS = 1.0e-12
RECOVERABILITY_FLOOR = 0.72
CONTROL_MARGIN = 0.05


@dataclass(frozen=True)
class MinimalConfig:
    output_dir: Path = ROOT / "outputs"
    seed: int = 20260425
    steps: int = 120
    dt: float = 0.08
    diffusion: float = 0.15
    address_gain: float = 0.45
    address_mode: str = "spectral"
    spectral_modes: int = 8
    hybrid_spectral_weight: float = 0.7
    invariant_gain: float = 0.025
    perturbation_step: int = 20
    perturbation_scale: float = 0.85
    permutation_trials: int = 64
    identity_bins: int = 4
    null_level: int = 3
    spectral_null_candidates: int = 8
    focus_address: bool = False
    max_nodes: int = 96


@dataclass(frozen=True)
class GraphData:
    adjacency: np.ndarray
    positions: np.ndarray
    source_label: str
    source_kind: str


@dataclass(frozen=True)
class MinimalRun:
    label: str
    states: np.ndarray
    rows: list[dict[str, float | int | str]]
    summary: dict[str, Any]


@dataclass(frozen=True)
class RuleComponents:
    target_address: np.ndarray
    target_shell_variance: np.ndarray
    diffusion_gain: float
    address_gain: float
    invariant_gain: float
    address_mode: str
    spectral_basis: np.ndarray | None = None
    spectral_target: np.ndarray | None = None
    shell_weights: np.ndarray | None = None
    hybrid_spectral_weight: float = 0.7


def run_minimal_probe(config: MinimalConfig) -> dict[str, Any]:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(config.seed)
    graph = load_geometry_graph(config)
    reference = build_reference_field(graph.positions)
    perturbation_nodes = choose_perturbation_nodes(graph.positions, rng)

    observed = run_single("observed", graph.adjacency, reference, perturbation_nodes, config, rng)
    controls = [
        run_single("node_shuffle_control", graph.adjacency, reference, perturbation_nodes, config, rng, initial_override=rng.permutation(reference)),
        run_single("edge_shuffle_control", edge_shuffle(graph.adjacency, rng), reference, perturbation_nodes, config, rng),
        run_single("degree_preserving_rewire_control", degree_preserving_rewire(graph.adjacency, rng), reference, perturbation_nodes, config, rng),
        run_single("topology_randomization_control", topology_randomize(graph.adjacency, rng), reference, perturbation_nodes, config, rng),
    ]
    ablations = run_rule_ablations(graph.adjacency, reference, perturbation_nodes, config)
    status = classify(graph, observed, controls)
    status["address_mode"] = config.address_mode
    status["spectral_modes"] = int(config.spectral_modes)
    status["hybrid_spectral_weight"] = float(config.hybrid_spectral_weight)
    status["null_level"] = int(config.null_level)
    status["spectral_null_candidates"] = int(config.spectral_null_candidates)
    status["ablation_summaries"] = {run.label: run.summary for run in ablations}
    status["ablation_drops"] = summarize_ablation_drops(ablations)
    write_outputs(config.output_dir, status, observed, controls, ablations)
    return status


def run_single(
    label: str,
    adjacency: np.ndarray,
    reference: np.ndarray,
    perturbation_nodes: np.ndarray,
    config: MinimalConfig,
    rng: np.random.Generator,
    initial_override: np.ndarray | None = None,
    rule_variant: str = "full",
) -> MinimalRun:
    transition = build_transition(adjacency)
    laplacian = np.eye(adjacency.shape[0]) - transition
    components = build_rule_components(reference, adjacency, config, rule_variant)
    initial = reference + 0.04 * rng.normal(size=reference.shape) if initial_override is None else initial_override
    states = np.zeros((config.steps + 1, reference.size), dtype=float)
    states[0] = normalize_field(initial, reference)

    for step in range(1, config.steps + 1):
        prior = states[step - 1].copy()
        if step == config.perturbation_step:
            prior[perturbation_nodes] += config.perturbation_scale * rng.normal(size=perturbation_nodes.size)
            prior = normalize_field(prior, reference)
        address_pull = address_restoration_pull(prior, reference, adjacency, components)
        invariant_pull = shell_variance_restoration(prior, adjacency, components.target_shell_variance)
        update = prior - config.dt * components.diffusion_gain * (laplacian @ prior)
        update += config.dt * components.address_gain * address_pull
        update += config.dt * components.invariant_gain * invariant_pull
        states[step] = normalize_field(update, reference)

    rows = compute_rows(label, states, reference, adjacency, config)
    summary = summarize_run(label, rows, states[-1], reference, adjacency, config)
    return MinimalRun(label, states, rows, summary)


def build_rule_components(
    reference: np.ndarray,
    adjacency: np.ndarray,
    config: MinimalConfig,
    rule_variant: str,
) -> RuleComponents:
    target_address = local_address(reference, adjacency)
    target_shell_variance = local_shell_variance(reference, adjacency)
    diffusion_gain = config.diffusion
    address_gain = config.address_gain
    invariant_gain = config.invariant_gain
    address_mode = config.address_mode
    spectral_basis: np.ndarray | None = None
    spectral_target: np.ndarray | None = None
    shell_weights: np.ndarray | None = None
    hybrid_spectral_weight = config.hybrid_spectral_weight

    if rule_variant == "no_address":
        address_gain = 0.0
    elif rule_variant == "no_invariant":
        invariant_gain = 0.0
    elif rule_variant == "diffusion_only":
        address_gain = 0.0
        invariant_gain = 0.0
    elif rule_variant == "address_only":
        diffusion_gain = 0.0
        invariant_gain = 0.0
    elif rule_variant == "address_local":
        address_mode = "local"
    elif rule_variant == "address_spectral":
        address_mode = "spectral"
    elif rule_variant == "address_hybrid":
        address_mode = "hybrid"
    elif rule_variant == "address_phase":
        address_mode = "phase"
    elif rule_variant == "address_multi_scale":
        address_mode = "multi_scale"
        shell_weights = shell_address_weights(reference, adjacency, config.identity_bins)
    elif rule_variant.startswith("address_weight_"):
        try:
            multiplier = float(rule_variant.removeprefix("address_weight_"))
        except ValueError as exc:
            raise ValueError(f"Invalid address weight variant: {rule_variant}") from exc
        address_gain = config.address_gain * multiplier
    elif rule_variant.startswith("address_gain_"):
        try:
            address_gain = float(rule_variant.removeprefix("address_gain_"))
        except ValueError as exc:
            raise ValueError(f"Invalid address gain variant: {rule_variant}") from exc
    elif rule_variant == "randomized_targets":
        rng = np.random.default_rng(config.seed + 50_000)
        target_address = rng.permutation(target_address)
        target_shell_variance = rng.permutation(target_shell_variance)
    elif rule_variant != "full":
        raise ValueError(f"Unknown rule_variant: {rule_variant}")

    if address_mode not in {"local", "spectral", "hybrid", "phase", "multi_scale"}:
        raise ValueError(f"Unknown address_mode: {address_mode}")
    if address_mode in {"spectral", "hybrid"}:
        spectral_basis = spectral_address_basis(adjacency, config.spectral_modes)
        spectral_target = spectral_basis.T @ reference

    return RuleComponents(
        target_address=target_address,
        target_shell_variance=target_shell_variance,
        diffusion_gain=diffusion_gain,
        address_gain=address_gain,
        invariant_gain=invariant_gain,
        address_mode=address_mode,
        spectral_basis=spectral_basis,
        spectral_target=spectral_target,
        shell_weights=shell_weights,
        hybrid_spectral_weight=hybrid_spectral_weight,
    )


def run_rule_ablations(
    adjacency: np.ndarray,
    reference: np.ndarray,
    perturbation_nodes: np.ndarray,
    config: MinimalConfig,
) -> list[MinimalRun]:
    variants = [
        ("ablation_full", "full"),
        ("ablation_no_address", "no_address"),
        ("ablation_no_invariant", "no_invariant"),
        ("ablation_diffusion_only", "diffusion_only"),
        ("ablation_randomized_targets", "randomized_targets"),
        ("ablation_address_only", "address_only"),
        ("ablation_address_local", "address_local"),
        ("ablation_address_spectral", "address_spectral"),
        ("ablation_address_hybrid_spectral_local", "address_hybrid"),
        ("ablation_address_phase", "address_phase"),
        ("ablation_address_multi_scale", "address_multi_scale"),
        ("ablation_address_gain_0_000", "address_gain_0.0"),
        ("ablation_address_gain_0_075", "address_gain_0.075"),
        ("ablation_address_gain_0_150", "address_gain_0.15"),
        ("ablation_address_gain_0_300", "address_gain_0.30"),
        ("ablation_address_gain_0_450", "address_gain_0.45"),
    ]
    if config.focus_address:
        variants = [
            item
            for item in variants
            if item[0] in {"ablation_full", "ablation_no_address", "ablation_randomized_targets"} or item[0].startswith("ablation_address_")
        ]
    runs: list[MinimalRun] = []
    for label, variant in variants:
        # Identical initial noise and perturbation noise across ablations: only
        # the rule component changes.
        rng = np.random.default_rng(config.seed + 60_000)
        runs.append(run_single(label, adjacency, reference, perturbation_nodes, config, rng, rule_variant=variant))
    return runs


def summarize_ablation_drops(ablations: list[MinimalRun]) -> dict[str, dict[str, float]]:
    baseline = next((run for run in ablations if run.label == "ablation_full"), None)
    if baseline is None:
        return {}
    base = baseline.summary
    drops: dict[str, dict[str, float]] = {}
    for run in ablations:
        summary = run.summary
        drops[run.label] = {
            "recoverability_score_drop": float(base["recoverability_score"]) - float(summary["recoverability_score"]),
            "address_z_drop": float(base["address_specificity_z"]) - float(summary["address_specificity_z"]),
            "invariant_z_drop": float(base["invariant_specificity_z"]) - float(summary["invariant_specificity_z"]),
            "branch_identity_z_drop": float(base["branch_identity_z"]) - float(summary["branch_identity_z"]),
        }
    return drops


def local_address(field: np.ndarray, adjacency: np.ndarray) -> np.ndarray:
    """Weighted neighbor-difference signature for each node."""

    degree = np.maximum(adjacency.sum(axis=1), EPS)
    return (adjacency @ field) / degree - field


def address_restoration_pull(
    field: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    components: RuleComponents,
) -> np.ndarray:
    """Compute the active harmonic-address pull for a rule variant."""

    degree = np.maximum(adjacency.sum(axis=1), 1.0)
    local_pull = (components.target_address - local_address(field, adjacency)) / degree
    if components.address_mode == "local":
        return local_pull
    if components.address_mode == "spectral":
        return spectral_address_pull(field, components)
    if components.address_mode == "hybrid":
        spectral_pull = spectral_address_pull(field, components)
        weight = float(np.clip(components.hybrid_spectral_weight, 0.0, 1.0))
        return weight * spectral_pull + (1.0 - weight) * local_pull
    if components.address_mode == "phase":
        target_phase = np.angle(np.exp(1j * reference))
        current_phase = np.angle(np.exp(1j * field))
        return np.sin(target_phase - current_phase)
    if components.address_mode == "multi_scale":
        weights = components.shell_weights if components.shell_weights is not None else np.ones_like(field)
        return weights * local_pull
    raise ValueError(f"Unknown address mode: {components.address_mode}")


def spectral_address_pull(field: np.ndarray, components: RuleComponents) -> np.ndarray:
    if components.spectral_basis is None or components.spectral_target is None:
        raise ValueError("Spectral address mode requires basis and target coefficients.")
    return components.spectral_basis @ (components.spectral_target - components.spectral_basis.T @ field)


def spectral_address_basis(adjacency: np.ndarray, modes: int = 8) -> np.ndarray:
    """Low-frequency harmonic basis of the frozen normalized Laplacian."""

    transition = build_transition(adjacency)
    normalized_laplacian = np.eye(adjacency.shape[0]) - 0.5 * (transition + transition.T)
    eigenvalues, eigenvectors = np.linalg.eigh(normalized_laplacian)
    order = np.argsort(eigenvalues)
    stop = min(max(2, modes + 1), eigenvectors.shape[1])
    return eigenvectors[:, order[1:stop]]


def shell_address_weights(reference: np.ndarray, adjacency: np.ndarray, bins: int) -> np.ndarray:
    """Moderate local-address gain by frozen shell-variance strata."""

    labels = quantile_bins(local_shell_variance(reference, adjacency), max(2, bins))
    scale = labels.astype(float) / max(float(np.max(labels)), 1.0)
    return 0.75 + 0.5 * scale


def local_shell_variance(field: np.ndarray, adjacency: np.ndarray) -> np.ndarray:
    """Frozen branch-local contrast energy around each node.

    This is deliberately small-scope: for each node, measure the weighted
    variance of neighbor offsets relative to that node. It gives the dynamics a
    local invariant target that is not just "return to the reference value".
    """

    degree = np.maximum(adjacency.sum(axis=1), EPS)
    offsets = field[None, :] - field[:, None]
    return np.sum(adjacency * offsets * offsets, axis=1) / degree


def shell_variance_restoration(
    field: np.ndarray,
    adjacency: np.ndarray,
    target_shell_variance: np.ndarray,
) -> np.ndarray:
    """Return a bounded local pull toward the frozen shell-variance invariant."""

    neighbor_mean = field + local_address(field, adjacency)
    current = local_shell_variance(field, adjacency)
    error = target_shell_variance - current
    scale = np.maximum(np.abs(target_shell_variance) + np.abs(current), 1.0)
    direction = field - neighbor_mean
    return np.clip(error / scale, -1.0, 1.0) * direction


def address_retention(field: np.ndarray, reference: np.ndarray, adjacency: np.ndarray) -> float:
    score = address_error_score(field, reference, adjacency)
    return float(np.clip(1.0 + score, 0.0, 1.0))


def address_error_score(field: np.ndarray, reference: np.ndarray, adjacency: np.ndarray) -> float:
    ref = local_address(reference, adjacency)
    cur = local_address(field, adjacency)
    distance = float(np.linalg.norm(cur - ref) / max(np.linalg.norm(ref), EPS))
    return float(-distance)


def invariant_retention(field: np.ndarray, reference: np.ndarray, adjacency: np.ndarray) -> float:
    score = invariant_error_score(field, reference, adjacency)
    return float(np.clip(1.0 + score, 0.0, 1.0))


def invariant_error_score(field: np.ndarray, reference: np.ndarray, adjacency: np.ndarray) -> float:
    ref = local_shell_variance(reference, adjacency)
    cur = local_shell_variance(field, adjacency)
    distance = float(np.linalg.norm(cur - ref) / max(np.linalg.norm(ref), EPS))
    return float(-distance)


def address_specificity_test(
    field: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    trials: int,
    seed: int,
) -> dict[str, float | bool]:
    observed = address_error_score(field, reference, adjacency)
    rng = np.random.default_rng(seed)
    nulls = np.zeros(trials, dtype=float)
    for idx in range(trials):
        nulls[idx] = address_error_score(rng.permutation(field), reference, adjacency)
    mean_null = float(np.mean(nulls))
    std_null = float(np.std(nulls))
    z_score = float((observed - mean_null) / (std_null if std_null > 1.0e-8 else 1.0))
    p_value = float((np.count_nonzero(nulls >= observed) + 1) / (trials + 1))
    return {
        "address_specificity": observed,
        "address_retention": address_retention(field, reference, adjacency),
        "address_specificity_null_mean": mean_null,
        "address_specificity_p": p_value,
        "address_specificity_z": z_score,
        "address_specificity_pass": bool(p_value < 0.05 and z_score > 1.5),
    }


def invariant_specificity_test(
    field: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    trials: int,
    seed: int,
) -> dict[str, float | bool]:
    observed = invariant_error_score(field, reference, adjacency)
    rng = np.random.default_rng(seed)
    nulls = np.zeros(trials, dtype=float)
    for idx in range(trials):
        nulls[idx] = invariant_error_score(rng.permutation(field), reference, adjacency)
    mean_null = float(np.mean(nulls))
    std_null = float(np.std(nulls))
    z_score = float((observed - mean_null) / (std_null if std_null > 1.0e-8 else 1.0))
    p_value = float((np.count_nonzero(nulls >= observed) + 1) / (trials + 1))
    return {
        "invariant_specificity": observed,
        "invariant_retention": invariant_retention(field, reference, adjacency),
        "invariant_specificity_null_mean": mean_null,
        "invariant_specificity_p": p_value,
        "invariant_specificity_z": z_score,
        "invariant_specificity_pass": bool(p_value < 0.05 and z_score > 1.5),
    }


def branch_identity_specificity_test(
    field: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    trials: int,
    seed: int,
    bins: int,
) -> dict[str, float | bool | int]:
    """Specificity test with degree/shell-stat preserving node shuffles.

    The global permutation tests ask whether the final field keeps any relational
    structure. This stricter null asks whether the signal survives after
    preserving two cheap local summaries: weighted degree and frozen shell
    variance. Only values inside matched degree/shell buckets are permuted, so
    an observed win is harder to explain as generic density or contrast effects.
    """

    observed = branch_identity_score(field, reference, adjacency)
    rng = np.random.default_rng(seed)
    nulls = np.zeros(trials, dtype=float)
    for idx in range(trials):
        shuffled = degree_shell_stratified_permutation(field, reference, adjacency, rng, bins)
        nulls[idx] = branch_identity_score(shuffled, reference, adjacency)
    mean_null = float(np.mean(nulls))
    std_null = float(np.std(nulls))
    z_score = float((observed - mean_null) / (std_null if std_null > 1.0e-8 else 1.0))
    p_value = float((np.count_nonzero(nulls >= observed) + 1) / (trials + 1))
    return {
        "branch_identity_specificity": observed,
        "branch_identity_null_mean": mean_null,
        "branch_identity_p": p_value,
        "branch_identity_z": z_score,
        "branch_identity_bins": int(bins),
        "branch_identity_pass": bool(p_value < 0.05 and z_score > 1.5),
    }


def spectral_aware_specificity_test(
    field: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    trials: int,
    seed: int,
    bins: int,
    modes: int,
    candidate_pool: int,
    include_autocorrelation: bool = False,
) -> dict[str, float | bool | int]:
    """Branch specificity against a spectral-aware stratified null.

    Each null draw still shuffles only within degree/shell buckets, but it now
    keeps the candidate whose low-mode spectral-energy signature is closest to
    the observed final field. At null level 3 the selector also matches one-hop
    local autocorrelation. This directly pressures spectral address smoothing:
    the null is allowed to keep the cheap spectral shape while breaking node
    identity.
    """

    observed = branch_identity_score(field, reference, adjacency)
    observed_signature = spectral_field_signature(field, adjacency, modes)
    observed_autocorr = local_autocorrelation(field, adjacency)
    rng = np.random.default_rng(seed)
    nulls = np.zeros(trials, dtype=float)
    pool = max(1, int(candidate_pool))
    for idx in range(trials):
        best_distance = math.inf
        best_candidate = field
        for _ in range(pool):
            candidate = degree_shell_stratified_permutation(field, reference, adjacency, rng, bins)
            distance = float(np.linalg.norm(spectral_field_signature(candidate, adjacency, modes) - observed_signature))
            if include_autocorrelation:
                distance += abs(local_autocorrelation(candidate, adjacency) - observed_autocorr)
            if distance < best_distance:
                best_distance = distance
                best_candidate = candidate
        nulls[idx] = branch_identity_score(best_candidate, reference, adjacency)
    mean_null = float(np.mean(nulls))
    std_null = float(np.std(nulls))
    z_score = float((observed - mean_null) / (std_null if std_null > 1.0e-8 else 1.0))
    p_value = float((np.count_nonzero(nulls >= observed) + 1) / (trials + 1))
    prefix = "autocorr_aware" if include_autocorrelation else "spectral_aware"
    return {
        f"{prefix}_specificity": observed,
        f"{prefix}_null_mean": mean_null,
        f"{prefix}_p": p_value,
        f"{prefix}_z": z_score,
        f"{prefix}_candidate_pool": int(pool),
        f"{prefix}_pass": bool(p_value < 0.05 and z_score > 1.5),
    }


def higher_order_specificity_test(
    field: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    trials: int,
    seed: int,
    bins: int,
    modes: int,
    candidate_pool: int,
) -> dict[str, float | bool | int]:
    """Specificity against a 2-hop/triangle-aware stratified null."""

    observed = branch_identity_score(field, reference, adjacency)
    observed_signature = higher_order_signature(field, adjacency, modes)
    rng = np.random.default_rng(seed)
    nulls = np.zeros(trials, dtype=float)
    pool = max(1, int(candidate_pool))
    for idx in range(trials):
        best_distance = math.inf
        best_candidate = field
        for _ in range(pool):
            candidate = degree_shell_stratified_permutation(field, reference, adjacency, rng, bins)
            distance = float(np.linalg.norm(higher_order_signature(candidate, adjacency, modes) - observed_signature))
            if distance < best_distance:
                best_distance = distance
                best_candidate = candidate
        nulls[idx] = branch_identity_score(best_candidate, reference, adjacency)
    mean_null = float(np.mean(nulls))
    std_null = float(np.std(nulls))
    z_score = float((observed - mean_null) / (std_null if std_null > 1.0e-8 else 1.0))
    p_value = float((np.count_nonzero(nulls >= observed) + 1) / (trials + 1))
    return {
        "higher_order_specificity": observed,
        "higher_order_null_mean": mean_null,
        "higher_order_p": p_value,
        "higher_order_z": z_score,
        "higher_order_candidate_pool": int(pool),
        "higher_order_pass": bool(p_value < 0.05 and z_score > 1.5),
    }


def branch_identity_score(field: np.ndarray, reference: np.ndarray, adjacency: np.ndarray) -> float:
    """Combined address + invariant score used by the stricter identity null."""

    return float(0.5 * (address_error_score(field, reference, adjacency) + invariant_error_score(field, reference, adjacency)))


def spectral_field_signature(field: np.ndarray, adjacency: np.ndarray, modes: int) -> np.ndarray:
    basis = spectral_address_basis(adjacency, modes)
    coeffs = basis.T @ field
    energy = coeffs * coeffs
    total = max(float(np.sum(energy)), EPS)
    return energy / total


def local_autocorrelation(field: np.ndarray, adjacency: np.ndarray) -> float:
    weights = adjacency[adjacency > 0.0]
    if weights.size == 0:
        return 0.0
    centered = field - float(np.mean(field))
    numerator = float(np.sum(adjacency * centered[:, None] * centered[None, :]))
    denominator = float(np.sum(adjacency) * max(np.var(centered), EPS))
    return numerator / max(denominator, EPS)


def higher_order_signature(field: np.ndarray, adjacency: np.ndarray, modes: int) -> np.ndarray:
    two_hop = two_hop_adjacency(adjacency)
    triangle = triangle_supported_adjacency(adjacency)
    return np.concatenate(
        [
            spectral_field_signature(field, adjacency, modes),
            np.array(
                [
                    local_autocorrelation(field, adjacency),
                    local_autocorrelation(field, two_hop),
                    local_autocorrelation(field, triangle),
                ],
                dtype=float,
            ),
        ]
    )


def two_hop_adjacency(adjacency: np.ndarray) -> np.ndarray:
    binary = (adjacency > 0.0).astype(float)
    two_hop = binary @ binary
    two_hop[binary > 0.0] = 0.0
    np.fill_diagonal(two_hop, 0.0)
    return normalize_adjacency(two_hop)


def triangle_supported_adjacency(adjacency: np.ndarray) -> np.ndarray:
    binary = (adjacency > 0.0).astype(float)
    common = binary @ binary
    supported = adjacency * common
    np.fill_diagonal(supported, 0.0)
    return normalize_adjacency(supported)


def degree_shell_stratified_permutation(
    field: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    rng: np.random.Generator,
    bins: int,
) -> np.ndarray:
    """Shuffle field values within matched weighted-degree/shell-variance buckets."""

    labels = degree_shell_strata(reference, adjacency, bins)
    out = np.asarray(field, dtype=float).copy()
    moved = False
    for label in np.unique(labels):
        idx = np.flatnonzero(labels == label)
        if idx.size < 2:
            continue
        permuted_idx = rng.permutation(idx)
        if np.array_equal(permuted_idx, idx):
            permuted_idx = np.roll(idx, 1)
        out[idx] = field[permuted_idx]
        moved = True
    return out if moved else rng.permutation(field)


def degree_shell_strata(reference: np.ndarray, adjacency: np.ndarray, bins: int) -> np.ndarray:
    safe_bins = max(2, int(bins))
    degree = adjacency.sum(axis=1)
    shell = local_shell_variance(reference, adjacency)
    return quantile_bins(degree, safe_bins) * safe_bins + quantile_bins(shell, safe_bins)


def quantile_bins(values: np.ndarray, bins: int) -> np.ndarray:
    if values.size == 0:
        return np.zeros(0, dtype=int)
    edges = np.quantile(values, np.linspace(0.0, 1.0, bins + 1)[1:-1])
    return np.digitize(values, edges, right=True).astype(int)


def compute_rows(
    label: str,
    states: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    config: MinimalConfig,
) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    prev = None
    for step, state in enumerate(states):
        recoverability = state_recoverability(state, reference)
        address_score = address_retention(state, reference, adjacency)
        invariant_score = invariant_retention(state, reference, adjacency)
        delta = 0.0 if prev is None else recoverability - prev
        prev = recoverability
        rows.append(
            {
                "run": label,
                "step": step,
                "time": float(step * config.dt),
                "recoverability": recoverability,
                "delta_persistence": float(delta),
                "address_retention": address_score,
                "invariant_retention": invariant_score,
            }
        )
    return rows


def summarize_run(
    label: str,
    rows: list[dict[str, float | int | str]],
    final_state: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    config: MinimalConfig,
) -> dict[str, Any]:
    recoverability = np.array([float(row["recoverability"]) for row in rows])
    address_values = np.array([float(row["address_retention"]) for row in rows])
    invariant_values = np.array([float(row["invariant_retention"]) for row in rows])
    start = min(config.perturbation_step, len(rows) - 1)
    recovery_window = recoverability[max(start, len(rows) // 2) :]
    specificity = address_specificity_test(final_state, reference, adjacency, config.permutation_trials, specificity_seed(config, label, 0))
    invariant_specificity = invariant_specificity_test(
        final_state,
        reference,
        adjacency,
        config.permutation_trials,
        specificity_seed(config, label, 10_000),
    )
    branch_identity_specificity = branch_identity_specificity_test(
        final_state,
        reference,
        adjacency,
        config.permutation_trials,
        specificity_seed(config, label, 20_000),
        config.identity_bins,
    )
    spectral_aware_specificity = spectral_aware_specificity_test(
        final_state,
        reference,
        adjacency,
        config.permutation_trials,
        specificity_seed(config, label, 30_000),
        config.identity_bins,
        config.spectral_modes,
        config.spectral_null_candidates,
        include_autocorrelation=False,
    )
    autocorr_aware_specificity = spectral_aware_specificity_test(
        final_state,
        reference,
        adjacency,
        config.permutation_trials,
        specificity_seed(config, label, 40_000),
        config.identity_bins,
        config.spectral_modes,
        config.spectral_null_candidates,
        include_autocorrelation=True,
    )
    higher_order_specificity = higher_order_specificity_test(
        final_state,
        reference,
        adjacency,
        config.permutation_trials,
        specificity_seed(config, label, 50_000),
        config.identity_bins,
        config.spectral_modes,
        config.spectral_null_candidates,
    )
    combined_specificity_pass = bool(specificity["address_specificity_pass"] and invariant_specificity["invariant_specificity_pass"])
    null_level = int(np.clip(config.null_level, 1, 4))
    strict_specificity_pass = bool(combined_specificity_pass and branch_identity_specificity["branch_identity_pass"])
    if null_level >= 2:
        strict_specificity_pass = bool(strict_specificity_pass and spectral_aware_specificity["spectral_aware_pass"])
    if null_level >= 3:
        strict_specificity_pass = bool(strict_specificity_pass and autocorr_aware_specificity["autocorr_aware_pass"])
    if null_level >= 4:
        strict_specificity_pass = bool(strict_specificity_pass and higher_order_specificity["higher_order_pass"])
    return {
        "run": label,
        "recoverability_score": float(np.mean(recovery_window)),
        "final_recoverability": float(recoverability[-1]),
        "min_recoverability": float(np.min(recoverability)),
        "final_address_retention": float(address_values[-1]),
        "min_address_retention": float(np.min(address_values)),
        "final_invariant_retention": float(invariant_values[-1]),
        "min_invariant_retention": float(np.min(invariant_values)),
        "recovered": bool(float(recoverability[-1]) >= RECOVERABILITY_FLOOR),
        **specificity,
        **invariant_specificity,
        **branch_identity_specificity,
        **spectral_aware_specificity,
        **autocorr_aware_specificity,
        **higher_order_specificity,
        "combined_specificity_pass": combined_specificity_pass,
        "null_level": null_level,
        "strict_specificity_pass": strict_specificity_pass,
    }


def classify(graph: GraphData, observed: MinimalRun, controls: list[MinimalRun]) -> dict[str, Any]:
    best_control = max(float(run.summary["recoverability_score"]) for run in controls)
    control_contrast = float(observed.summary["recoverability_score"]) >= best_control + CONTROL_MARGIN
    control_specificity_passes = sum(1 for run in controls if bool(run.summary["combined_specificity_pass"]))
    control_strict_specificity_passes = sum(1 for run in controls if bool(run.summary["strict_specificity_pass"]))
    control_spectral_passes = sum(1 for run in controls if bool(run.summary["spectral_aware_pass"]))
    control_autocorr_passes = sum(1 for run in controls if bool(run.summary["autocorr_aware_pass"]))
    control_higher_order_passes = sum(1 for run in controls if bool(run.summary["higher_order_pass"]))
    observed_specific = bool(observed.summary["strict_specificity_pass"])
    observed_recovered = bool(observed.summary["recovered"])
    null_level = int(observed.summary.get("null_level", 1))

    if graph.source_kind == "OPEN_NO_DATA_SYNTHETIC":
        status = "OPEN_NO_DATA_SYNTHETIC"
        failure = "NO_HAOS_GRAPH_ARTIFACT_FOUND"
    elif observed_recovered and observed_specific and control_contrast and control_strict_specificity_passes == 0:
        status = "PASS"
        failure = "RECOVERY_WITH_SPECTRAL_AWARE_SPECIFICITY"
    elif observed_recovered and observed_specific and control_strict_specificity_passes > 0:
        status = "MARGINAL"
        if null_level >= 4:
            failure = "HIGHER_ORDER_CONTROL_MATCH"
        elif null_level >= 3:
            failure = "AUTOCORR_AWARE_CONTROL_MATCH"
        else:
            failure = "SPECTRAL_AWARE_CONTROL_MATCH"
    elif observed_recovered and observed_specific:
        status = "MARGINAL"
        failure = "STRICT_SIGNAL_WITHOUT_CONTROL_CONTRAST"
    elif not observed_recovered:
        status = "FAIL"
        failure = "RECOVERABILITY_COLLAPSE"
    else:
        status = "FAIL"
        failure = "STRICT_BRANCH_IDENTITY_SPECIFICITY_FAILED"

    return {
        "bridge_status": status,
        "failure_mode": failure,
        "graph_source": graph.source_label,
        "graph_source_kind": graph.source_kind,
        "nodes": int(graph.adjacency.shape[0]),
        "edges": int(np.count_nonzero(np.triu(graph.adjacency, k=1))),
        "control_contrast": control_contrast,
        "control_combined_specificity_pass_count": control_specificity_passes,
        "control_strict_specificity_pass_count": control_strict_specificity_passes,
        "control_spectral_aware_pass_count": control_spectral_passes,
        "control_autocorr_aware_pass_count": control_autocorr_passes,
        "control_higher_order_pass_count": control_higher_order_passes,
        "leakage_controls": [run.label for run in controls if bool(run.summary["strict_specificity_pass"])],
        "best_control_recoverability_score": best_control,
        "observed_summary": observed.summary,
        "control_summaries": {run.label: run.summary for run in controls},
    }


def write_outputs(
    output_dir: Path,
    status: dict[str, Any],
    observed: MinimalRun,
    controls: list[MinimalRun],
    ablations: list[MinimalRun],
) -> None:
    write_json(output_dir / "bridge_status.json", status)
    write_csv(output_dir / "minimal_timeseries.csv", [observed, *controls])
    write_report(output_dir / "minimal_dynamics_report.md", status, observed, controls)
    write_ablation_report(output_dir / "ablation_report.md", ablations, status["ablation_drops"])
    write_plots(output_dir, observed, controls, ablations)


def write_report(path: Path, status: dict[str, Any], observed: MinimalRun, controls: list[MinimalRun]) -> None:
    lines = [
        "# HAOS Minimal Dynamics Report",
        "",
        f"- bridge_status: {status['bridge_status']}",
        f"- failure_mode: {status['failure_mode']}",
        f"- graph_source: {status['graph_source']}",
        f"- address_mode: {status.get('address_mode', 'unknown')}",
        f"- spectral_modes: {status.get('spectral_modes', 'unknown')}",
        f"- null_level: {status.get('null_level', 'unknown')}",
        f"- control_contrast: {status['control_contrast']}",
        f"- control_combined_specificity_pass_count: {status['control_combined_specificity_pass_count']}",
        f"- control_strict_specificity_pass_count: {status['control_strict_specificity_pass_count']}",
        f"- control_spectral_aware_pass_count: {status['control_spectral_aware_pass_count']}",
        f"- control_autocorr_aware_pass_count: {status['control_autocorr_aware_pass_count']}",
        f"- control_higher_order_pass_count: {status['control_higher_order_pass_count']}",
        "",
        "## Spectral Smoothing Diagnosis",
        "",
        "The strict gate now includes a spectral-aware null. Null candidates preserve degree/shell buckets and are selected to match low-mode spectral energy before branch identity is scored.",
        "",
        "## Leakage Analysis",
        "",
        f"- leakage_controls: {', '.join(status.get('leakage_controls', [])) or 'none'}",
        "",
        "| run | recoverability_score | final_recoverability | branch_z | spectral_z | autocorr_z | higher_z | strict_pass |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for run in [observed, *controls]:
        summary = run.summary
        lines.append(
            "| {run} | {score:.6f} | {final:.6f} | {branch_z:.6f} | {spectral_z:.6f} | {autocorr_z:.6f} | {higher_z:.6f} | {passed} |".format(
                run=run.label,
                score=float(summary["recoverability_score"]),
                final=float(summary["final_recoverability"]),
                branch_z=float(summary["branch_identity_z"]),
                spectral_z=float(summary["spectral_aware_z"]),
                autocorr_z=float(summary["autocorr_aware_z"]),
                higher_z=float(summary["higher_order_z"]),
                passed=summary["strict_specificity_pass"],
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_ablation_report(path: Path, ablations: list[MinimalRun], drops: dict[str, dict[str, float]]) -> None:
    best_address = best_address_variant(ablations)
    full = next((run for run in ablations if run.label == "ablation_full"), None)
    old_local = next((run for run in ablations if run.label == "ablation_address_local"), None)
    spectral_delta = (
        float(full.summary["branch_identity_z"]) - float(old_local.summary["branch_identity_z"])
        if full is not None and old_local is not None
        else None
    )
    lines = [
        "# HAOS Minimal Dynamics Ablation Report",
        "",
        "Ablations reuse the same initial noise and perturbation noise. Only the dynamics rule component changes.",
        "",
        "## Address Term Diagnosis",
        "",
        f"- best_address_variant: {best_address.label if best_address else 'none'}",
        f"- best_address_branch_z: {float(best_address.summary['branch_identity_z']):.6f}" if best_address else "- best_address_branch_z: n/a",
        f"- full_branch_z: {float(full.summary['branch_identity_z']):.6f}" if full else "- full_branch_z: n/a",
        f"- old_local_branch_z: {float(old_local.summary['branch_identity_z']):.6f}" if old_local else "- old_local_branch_z: n/a",
        f"- spectral_default_delta_vs_local: {spectral_delta:.6f}" if spectral_delta is not None else "- spectral_default_delta_vs_local: n/a",
        f"- address_beats_full: {bool(best_address and full and float(best_address.summary['branch_identity_z']) > float(full.summary['branch_identity_z']))}",
        "",
        "| run | recoverability_score | address_z | invariant_z | branch_z | branch_z_drop | strict_pass |",
        "| --- | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for run in ablations:
        summary = run.summary
        run_drops = drops.get(run.label, {})
        lines.append(
            "| {run} | {score:.6f} | {address_z:.6f} | {invariant_z:.6f} | {branch_z:.6f} | {branch_drop:.6f} | {passed} |".format(
                run=run.label,
                score=float(summary["recoverability_score"]),
                address_z=float(summary["address_specificity_z"]),
                invariant_z=float(summary["invariant_specificity_z"]),
                branch_z=float(summary["branch_identity_z"]),
                branch_drop=float(run_drops.get("branch_identity_z_drop", 0.0)),
                passed=summary["strict_specificity_pass"],
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def best_address_variant(ablations: list[MinimalRun]) -> MinimalRun | None:
    candidates = [run for run in ablations if run.label.startswith("ablation_address_")]
    if not candidates:
        return None
    return max(candidates, key=lambda run: float(run.summary["branch_identity_z"]))


def write_csv(path: Path, runs: list[MinimalRun]) -> None:
    fields = ["run", "step", "time", "recoverability", "delta_persistence", "address_retention", "invariant_retention"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for run in runs:
            writer.writerows(run.rows)


def write_plots(output_dir: Path, observed: MinimalRun, controls: list[MinimalRun], ablations: list[MinimalRun]) -> None:
    import matplotlib.pyplot as plt

    for key, filename, ylabel in (
        ("recoverability", "recoverability.png", "recoverability"),
        ("address_retention", "address_specificity.png", "address retention"),
        ("invariant_retention", "invariant_retention.png", "invariant retention"),
    ):
        plt.figure(figsize=(8, 4))
        for run in [observed, *controls]:
            plt.plot([float(row["time"]) for row in run.rows], [float(row[key]) for row in run.rows], label=run.label)
        plt.xlabel("time")
        plt.ylabel(ylabel)
        plt.title(f"Minimal dynamics {ylabel}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_dir / filename, dpi=160)
        plt.close()

    labels = [run.label.replace("ablation_", "") for run in ablations]
    branch_z = [float(run.summary["branch_identity_z"]) for run in ablations]
    plt.figure(figsize=(8, 4))
    plt.bar(labels, branch_z)
    plt.ylabel("branch identity z")
    plt.title("Dynamics rule ablation")
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    plt.savefig(output_dir / "ablation_branch_identity.png", dpi=160)
    plt.close()

    runs = [observed, *controls]
    labels = [run.label for run in runs]
    x = np.arange(len(runs), dtype=float)
    width = 0.26
    plt.figure(figsize=(9, 4))
    plt.bar(x - width, [float(run.summary["branch_identity_z"]) for run in runs], width, label="degree/shell")
    plt.bar(x, [float(run.summary["spectral_aware_z"]) for run in runs], width, label="spectral-aware")
    plt.bar(x + width, [float(run.summary["higher_order_z"]) for run in runs], width, label="higher-order")
    plt.ylabel("specificity z")
    plt.title("Null erosion")
    plt.xticks(x, labels, rotation=25, ha="right")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "null_erosion_spectral_vs_local.png", dpi=160)
    plt.close()


def load_geometry_graph(config: MinimalConfig) -> GraphData:
    path = REPO_ROOT / "geometry_emergence" / "experiments" / "default_geometry_probe_config.json"
    if not path.exists():
        return synthetic_ring_graph(min(48, config.max_nodes))
    payload = json.loads(path.read_text(encoding="utf-8"))
    n_nodes = min(int(payload.get("n_nodes", 72)), config.max_nodes)
    dim = int(payload.get("embedding_dim", 2))
    clusters = int(payload.get("cluster_count", 3))
    spread = float(payload.get("cluster_spread", 0.075))
    widths = payload.get("kernel_widths", [0.13])
    kernel_width = float(widths[len(widths) // 2])
    radius = float(payload.get("locality_radius_factor", 3.0)) * kernel_width
    rng = np.random.default_rng(int(payload.get("seed", 17)))
    centers = rng.uniform(0.15, 0.85, size=(clusters, dim))
    counts = np.full(clusters, n_nodes // clusters, dtype=int)
    counts[: n_nodes % clusters] += 1
    assignments = np.repeat(np.arange(clusters), counts)
    rng.shuffle(assignments)
    positions = np.clip(centers[assignments] + rng.normal(0.0, spread, size=(n_nodes, dim)), 0.0, 1.0)
    distances = pairwise_distances(positions)
    adjacency = np.exp(-(distances ** 2) / max(2.0 * kernel_width * kernel_width, EPS))
    adjacency[distances > radius] = 0.0
    np.fill_diagonal(adjacency, 0.0)
    return GraphData(normalize_adjacency(adjacency), positions, "geometry_emergence/default_geometry_probe_config.json", "haos_config_reconstruction")


def build_reference_field(positions: np.ndarray) -> np.ndarray:
    x = positions[:, 0]
    y = positions[:, 1] if positions.shape[1] > 1 else np.zeros_like(x)
    field = np.sin(2.0 * math.pi * x) + 0.55 * np.cos(2.0 * math.pi * y)
    return normalize_field(field, field)


def choose_perturbation_nodes(positions: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    center = positions[int(rng.integers(0, positions.shape[0]))]
    distances = np.linalg.norm(positions - center, axis=1)
    return np.argsort(distances)[: max(4, int(round(0.12 * positions.shape[0])))]


def state_recoverability(state: np.ndarray, reference: np.ndarray) -> float:
    distance = float(np.linalg.norm(state - reference) / max(math.sqrt(reference.size), EPS))
    return float(np.clip(1.0 - distance, 0.0, 1.0))


def build_transition(adjacency: np.ndarray) -> np.ndarray:
    row_sums = adjacency.sum(axis=1, keepdims=True)
    transition = np.divide(adjacency, np.maximum(row_sums, EPS), out=np.zeros_like(adjacency), where=row_sums > EPS)
    isolated = np.flatnonzero(row_sums[:, 0] <= EPS)
    transition[isolated, isolated] = 1.0
    return transition


def normalize_field(field: np.ndarray, reference: np.ndarray) -> np.ndarray:
    centered = np.asarray(field, dtype=float) - float(np.mean(field))
    target = np.asarray(reference, dtype=float) - float(np.mean(reference))
    return centered * (max(float(np.linalg.norm(target)), EPS) / max(float(np.linalg.norm(centered)), EPS))


def normalize_adjacency(adjacency: np.ndarray) -> np.ndarray:
    arr = 0.5 * (np.asarray(adjacency, dtype=float) + np.asarray(adjacency, dtype=float).T)
    np.fill_diagonal(arr, 0.0)
    positive = arr[arr > 0.0]
    return arr / max(float(np.median(positive)), EPS) if positive.size else arr


def edge_shuffle(adjacency: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    upper = adjacency[np.triu_indices(adjacency.shape[0], k=1)]
    shuffled = rng.permutation(upper)
    out = np.zeros_like(adjacency)
    out[np.triu_indices(adjacency.shape[0], k=1)] = shuffled
    return normalize_adjacency(out + out.T)


def topology_randomize(adjacency: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    n = adjacency.shape[0]
    edge_count = int(np.count_nonzero(np.triu(adjacency, k=1)))
    weights = adjacency[np.triu_indices(n, k=1)]
    weights = weights[weights > 0.0]
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    selected = rng.choice(len(pairs), size=min(edge_count, len(pairs)), replace=False)
    sampled = rng.choice(weights, size=len(selected), replace=True) if weights.size else np.ones(len(selected))
    out = np.zeros_like(adjacency)
    for pair_idx, weight in zip(selected, sampled):
        i, j = pairs[int(pair_idx)]
        out[i, j] = float(weight)
        out[j, i] = float(weight)
    return normalize_adjacency(out)


def degree_preserving_rewire(adjacency: np.ndarray, rng: np.random.Generator, swap_factor: int = 8) -> np.ndarray:
    """Randomize topology with double-edge swaps while preserving unweighted degrees."""

    binary = np.asarray(adjacency > 0.0, dtype=bool)
    n = binary.shape[0]
    edge_i, edge_j = np.triu_indices(n, k=1)
    edges = [(int(i), int(j)) for i, j in zip(edge_i, edge_j) if binary[i, j]]
    edge_set = set(edges)
    if len(edges) < 2:
        return adjacency.copy()

    attempts = max(len(edges) * swap_factor * 4, 1)
    target_swaps = len(edges) * swap_factor
    swaps = 0
    for _ in range(attempts):
        if swaps >= target_swaps:
            break
        idx_a, idx_b = rng.choice(len(edges), size=2, replace=False)
        a, b = edges[int(idx_a)]
        c, d = edges[int(idx_b)]
        if len({a, b, c, d}) < 4:
            continue
        if rng.random() < 0.5:
            new_a = tuple(sorted((a, d)))
            new_b = tuple(sorted((c, b)))
        else:
            new_a = tuple(sorted((a, c)))
            new_b = tuple(sorted((b, d)))
        if new_a[0] == new_a[1] or new_b[0] == new_b[1] or new_a == new_b:
            continue
        if new_a in edge_set or new_b in edge_set:
            continue
        old_a = edges[int(idx_a)]
        old_b = edges[int(idx_b)]
        edge_set.remove(old_a)
        edge_set.remove(old_b)
        edge_set.add(new_a)
        edge_set.add(new_b)
        edges[int(idx_a)] = new_a
        edges[int(idx_b)] = new_b
        swaps += 1

    weights = adjacency[np.triu_indices(n, k=1)]
    weights = weights[weights > 0.0]
    shuffled_weights = rng.permutation(weights)
    out = np.zeros_like(adjacency, dtype=float)
    for (i, j), weight in zip(edges, shuffled_weights):
        out[i, j] = float(weight)
        out[j, i] = float(weight)
    return normalize_adjacency(out)


def synthetic_ring_graph(node_count: int) -> GraphData:
    angles = np.linspace(0.0, 2.0 * math.pi, node_count, endpoint=False)
    positions = np.column_stack([np.cos(angles), np.sin(angles)])
    adjacency = np.zeros((node_count, node_count), dtype=float)
    for idx in range(node_count):
        adjacency[idx, (idx - 1) % node_count] = 1.0
        adjacency[idx, (idx + 1) % node_count] = 1.0
    return GraphData(adjacency, positions, "synthetic_ring_fallback", "OPEN_NO_DATA_SYNTHETIC")


def pairwise_distances(points: np.ndarray) -> np.ndarray:
    delta = points[:, None, :] - points[None, :, :]
    return np.sqrt(np.sum(delta * delta, axis=2))


def stable_label_offset(label: str) -> int:
    return sum((idx + 1) * ord(char) for idx, char in enumerate(label))


def specificity_seed(config: MinimalConfig, label: str, offset: int) -> int:
    seed_label = "ablation_common" if label.startswith("ablation_") else label
    return int(config.seed + offset + stable_label_offset(seed_label))


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(to_jsonable(data), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def to_jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(item) for item in value]
    return value
