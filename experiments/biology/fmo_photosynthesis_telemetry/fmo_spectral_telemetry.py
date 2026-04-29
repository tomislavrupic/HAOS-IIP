from __future__ import annotations

from dataclasses import dataclass
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parent
EPS = 1.0e-12


@dataclass(frozen=True)
class FMOTelemetryConfig:
    output_dir: Path = ROOT / "outputs"
    seed: int = 20260429
    steps: int = 128
    dt: float = 0.06
    diffusion: float = 0.10
    address_gain: float = 0.48
    local_address_gain: float = 0.18
    sink_gain: float = 0.0
    flux_gain: float = 0.0
    environment_assist_gain: float = 0.0
    invariant_gain: float = 0.035
    spectral_modes: int = 4
    address_mode: str = "spectral"
    perturbation_step: int = 24
    thermal_noise: float = 0.08
    disorder_scale: float = 0.18
    damage_scale: float = 0.42
    permutation_trials: int = 32
    null_candidates: int = 4
    null_level: int = 5


@dataclass(frozen=True)
class FMORun:
    label: str
    states: np.ndarray
    summary: dict[str, Any]


def run_fmo_spectral_telemetry(config: FMOTelemetryConfig) -> dict[str, Any]:
    """Run a bounded HAOS-style spectral telemetry probe on a toy FMO-like network."""
    config.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(config.seed)
    adjacency = normalize_adjacency(fmo_adjacency())
    reference = normalize_field(fmo_reference_field(), fmo_reference_field())
    damage_nodes = np.asarray([0, 1, 5], dtype=int)

    observed = run_single("observed", adjacency, reference, damage_nodes, config, rng)
    shuffled_reference = rng.permutation(reference)
    controls = [
        run_single(
            "site_shuffle_control",
            adjacency,
            reference,
            damage_nodes,
            config,
            rng,
            initial_override=shuffled_reference,
            target_reference=shuffled_reference,
        ),
        run_single("edge_shuffle_control", edge_shuffle(adjacency, rng), reference, damage_nodes, config, rng),
        run_single("topology_randomization_control", topology_randomize(adjacency, rng), reference, damage_nodes, config, rng),
    ]
    status = classify(observed, controls, config, adjacency)
    write_outputs(config, status, observed, controls, adjacency, reference)
    return status


def fmo_hamiltonian_cm() -> np.ndarray:
    """Adolphs-Renger-style 7-site FMO Hamiltonian in cm^-1, used only as topology prior."""
    return np.asarray(
        [
            [240.0, -87.7, 5.5, -5.9, 6.7, -13.7, -9.9],
            [-87.7, 315.0, 30.8, 8.2, 0.7, 11.8, 4.3],
            [5.5, 30.8, 0.0, -53.5, -2.2, -9.6, 6.0],
            [-5.9, 8.2, -53.5, 130.0, -70.7, -17.0, -63.3],
            [6.7, 0.7, -2.2, -70.7, 285.0, 81.1, -1.3],
            [-13.7, 11.8, -9.6, -17.0, 81.1, 435.0, 39.7],
            [-9.9, 4.3, 6.0, -63.3, -1.3, 39.7, 245.0],
        ],
        dtype=float,
    )


def fmo_adjacency() -> np.ndarray:
    coupling = np.abs(fmo_hamiltonian_cm())
    np.fill_diagonal(coupling, 0.0)
    threshold = 4.0
    coupling[coupling < threshold] = 0.0
    return coupling


def fmo_reference_field() -> np.ndarray:
    """Toy transfer-address field: donor side high, reaction-center side structured low."""
    donor = np.asarray([1.00, 0.78, 0.22, -0.18, -0.50, 0.62, -0.72], dtype=float)
    reaction_center_pull = np.asarray([0.0, 0.05, 0.26, 0.45, 0.38, -0.04, 0.52], dtype=float)
    return donor - reaction_center_pull


def run_single(
    label: str,
    adjacency: np.ndarray,
    reference: np.ndarray,
    damage_nodes: np.ndarray,
    config: FMOTelemetryConfig,
    rng: np.random.Generator,
    initial_override: np.ndarray | None = None,
    target_reference: np.ndarray | None = None,
) -> FMORun:
    target = reference if target_reference is None else normalize_field(target_reference, reference)
    basis = spectral_basis(adjacency, config.spectral_modes)
    target_coeffs = basis.T @ target
    target_shell = local_shell_variance(target, adjacency)
    sink_profile = reaction_center_sink_profile(reference)
    target_flux = pathway_flux(target, adjacency)
    transition = build_transition(adjacency)
    laplacian = np.eye(adjacency.shape[0]) - transition
    initial = target + 0.035 * rng.normal(size=reference.shape) if initial_override is None else initial_override
    states = np.zeros((config.steps + 1, reference.size), dtype=float)
    states[0] = normalize_field(initial, reference)

    for step in range(1, config.steps + 1):
        prior = states[step - 1].copy()
        if step == config.perturbation_step:
            prior[damage_nodes] += config.damage_scale * rng.normal(size=damage_nodes.size)
        if config.thermal_noise > 0.0:
            noise = rng.normal(size=prior.shape)
            prior += config.thermal_noise * 0.02 * noise
            if config.environment_assist_gain > 0.0:
                prior += config.thermal_noise * config.environment_assist_gain * sink_profile * np.tanh(target - prior)
        if config.disorder_scale > 0.0:
            prior += config.disorder_scale * 0.002 * np.sign(reference) * rng.normal(size=prior.shape)
        spectral_pull = basis @ (target_coeffs - basis.T @ prior)
        local_pull = local_address_restoration(prior, target, adjacency)
        sink_pull = sink_profile * (target - prior)
        flux_pull = pathway_flux_restoration(prior, target_flux, adjacency)
        invariant_pull = shell_variance_restoration(prior, adjacency, target_shell)
        update = prior - config.dt * config.diffusion * (laplacian @ prior)
        update += config.dt * address_pull(prior, spectral_pull, local_pull, sink_pull, flux_pull, config)
        update += config.dt * config.invariant_gain * invariant_pull
        states[step] = normalize_field(update, reference)

    return FMORun(label, states, summarize_run(label, states, reference, adjacency, config))


def summarize_run(
    label: str,
    states: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    config: FMOTelemetryConfig,
) -> dict[str, Any]:
    final = states[-1]
    recoverability = state_recoverability(final, reference)
    site_identity = site_identity_retention(final, reference)
    pathway_identity = pathway_identity_retention(final, reference, adjacency)
    persistence = float(recoverability - state_recoverability(states[config.perturbation_step], reference))
    nulls = null_ladder(final, states, reference, adjacency, config, config.seed + stable_label_offset(label))
    strict = bool(recoverability >= 0.74 and site_identity >= 0.70 and pathway_identity >= 0.68 and nulls["active_null_pass"])
    return {
        "run": label,
        "recoverability_score": recoverability,
        "site_identity_retention": site_identity,
        "pathway_identity_retention": pathway_identity,
        "delta_persistence": persistence,
        "safety_margin": float(recoverability - run_null_mean(final, reference, adjacency, config)),
        **nulls,
        "strict_pass": strict,
    }


def address_pull(
    prior: np.ndarray,
    spectral_pull: np.ndarray,
    local_pull: np.ndarray,
    sink_pull: np.ndarray,
    flux_pull: np.ndarray,
    config: FMOTelemetryConfig,
) -> np.ndarray:
    if config.address_mode == "spectral":
        return config.address_gain * spectral_pull + config.sink_gain * sink_pull + config.flux_gain * flux_pull
    if config.address_mode == "local":
        return config.local_address_gain * local_pull + config.sink_gain * sink_pull + config.flux_gain * flux_pull
    if config.address_mode == "hybrid":
        return config.address_gain * spectral_pull + config.local_address_gain * local_pull + config.sink_gain * sink_pull + config.flux_gain * flux_pull
    if config.address_mode == "sink":
        return config.address_gain * spectral_pull + config.sink_gain * sink_pull + config.flux_gain * flux_pull
    if config.address_mode == "pathway_flux":
        return config.address_gain * spectral_pull + config.sink_gain * sink_pull + config.flux_gain * flux_pull
    if config.address_mode == "environment_assisted":
        return config.address_gain * spectral_pull + config.local_address_gain * local_pull + config.sink_gain * sink_pull + config.flux_gain * flux_pull
    raise ValueError(f"unknown address_mode: {config.address_mode}")


def local_address_restoration(field: np.ndarray, target: np.ndarray, adjacency: np.ndarray) -> np.ndarray:
    current = local_address(field, adjacency)
    desired = local_address(target, adjacency)
    neighbor_mean = field + current
    return (desired - current) + 0.25 * (target - neighbor_mean)


def reaction_center_sink_profile(reference: np.ndarray) -> np.ndarray:
    profile = np.zeros_like(reference, dtype=float)
    profile[[2, 3, 4, 6]] = np.asarray([0.35, 1.0, 0.60, 0.80], dtype=float)
    return profile / max(float(np.linalg.norm(profile)), EPS)


def pathway_edges() -> np.ndarray:
    return np.asarray([(0, 1), (1, 2), (2, 3), (5, 4), (4, 3), (3, 6)], dtype=int)


def pathway_flux(field: np.ndarray, adjacency: np.ndarray) -> np.ndarray:
    return np.asarray([adjacency[i, j] * (field[i] - field[j]) for i, j in pathway_edges()], dtype=float)


def pathway_flux_restoration(field: np.ndarray, target_flux: np.ndarray, adjacency: np.ndarray) -> np.ndarray:
    pull = np.zeros_like(field, dtype=float)
    current_flux = pathway_flux(field, adjacency)
    for (i, j), delta in zip(pathway_edges(), target_flux - current_flux):
        weight = adjacency[i, j]
        if weight <= EPS:
            continue
        correction = float(delta) / max(2.0 * weight, EPS)
        pull[i] += correction
        pull[j] -= correction
    norm = max(float(np.linalg.norm(pull)), 1.0)
    return np.clip(pull / norm, -1.0, 1.0)


def null_ladder(
    final: np.ndarray,
    states: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    config: FMOTelemetryConfig,
    seed: int,
) -> dict[str, float | bool]:
    rng = np.random.default_rng(seed)
    observed = fmo_identity_score(final, reference, adjacency)
    scores = {"branch": [], "spectral": [], "higher": [], "trajectory": []}
    signatures = {
        "spectral": spectral_field_signature(final, adjacency, config.spectral_modes),
        "higher": higher_order_signature(final, adjacency, config.spectral_modes),
        "trajectory": trajectory_signature(states, adjacency, config.spectral_modes),
    }
    for _ in range(config.permutation_trials):
        # FMO has only seven sites. Degree/shell buckets often become singleton
        # strata, which makes a stratified null unable to break site identity at
        # all. Use full site permutations, then let the stricter null levels
        # choose candidates that best preserve spectral/trajectory signatures.
        candidates = [rng.permutation(final.size) for _ in range(max(1, config.null_candidates))]
        scores["branch"].append(fmo_identity_score(final[candidates[0]], reference, adjacency))
        for key in ("spectral", "higher", "trajectory"):
            best = min(
                candidates,
                key=lambda idx: float(np.linalg.norm(signature_for(key, states[:, idx], final[idx], adjacency, config) - signatures[key])),
            )
            scores[key].append(fmo_identity_score(final[best], reference, adjacency))

    active = "branch"
    if config.null_level >= 2:
        active = "spectral"
    if config.null_level >= 4:
        active = "higher"
    if config.null_level >= 5:
        active = "trajectory"

    out: dict[str, float | bool] = {}
    for key, values in scores.items():
        arr = np.asarray(values, dtype=float)
        mean = float(np.mean(arr))
        std = float(np.std(arr))
        z = float((observed - mean) / (std if std > 1.0e-8 else 1.0))
        p = float((np.count_nonzero(arr >= observed) + 1) / (arr.size + 1))
        out[f"{key}_null_z"] = z
        out[f"{key}_null_p"] = p
        out[f"{key}_null_pass"] = bool(p < 0.05 and z > 1.5)
    out["active_null"] = active
    out["active_null_z"] = float(out[f"{active}_null_z"])
    out["active_null_pass"] = bool(out[f"{active}_null_pass"])
    return out


def signature_for(key: str, states: np.ndarray, final: np.ndarray, adjacency: np.ndarray, config: FMOTelemetryConfig) -> np.ndarray:
    if key == "spectral":
        return spectral_field_signature(final, adjacency, config.spectral_modes)
    if key == "higher":
        return higher_order_signature(final, adjacency, config.spectral_modes)
    if key == "trajectory":
        return trajectory_signature(states, adjacency, config.spectral_modes)
    raise ValueError(key)


def classify(observed: FMORun, controls: list[FMORun], config: FMOTelemetryConfig, adjacency: np.ndarray) -> dict[str, Any]:
    best_control = max(float(run.summary["recoverability_score"]) for run in controls)
    observed_composite = specificity_composite(observed.summary)
    best_control_composite = max(specificity_composite(run.summary) for run in controls)
    control_contrast = observed_composite >= best_control_composite + 0.04
    control_pass_count = sum(1 for run in controls if bool(run.summary["strict_pass"]))
    if bool(observed.summary["strict_pass"]) and control_contrast and control_pass_count == 0:
        status = "PASS"
        failure = "RECOVERABLE_FMO_SPECTRAL_TELEMETRY"
    elif bool(observed.summary["strict_pass"]):
        status = "MARGINAL"
        failure = "FMO_CONTROL_MATCH"
    else:
        status = "FAIL"
        failure = "FMO_TELEMETRY_SPECIFICITY_FAILED"
    return {
        "bridge_status": status,
        "failure_mode": failure,
        "nodes": int(adjacency.shape[0]),
        "edges": int(np.count_nonzero(np.triu(adjacency, k=1))),
        "null_level": int(config.null_level),
        "control_contrast": control_contrast,
        "control_strict_pass_count": control_pass_count,
        "best_control_recoverability_score": best_control,
        "observed_specificity_composite": observed_composite,
        "best_control_specificity_composite": best_control_composite,
        "observed_summary": observed.summary,
        "control_summaries": {run.label: run.summary for run in controls},
    }


def specificity_composite(summary: dict[str, Any]) -> float:
    return float(
        0.34 * float(summary["recoverability_score"])
        + 0.33 * float(summary["site_identity_retention"])
        + 0.33 * float(summary["pathway_identity_retention"])
    )


def site_identity_retention(final: np.ndarray, reference: np.ndarray) -> float:
    corr = np.corrcoef(final, reference)[0, 1] if np.std(final) > EPS and np.std(reference) > EPS else 0.0
    return float(np.clip(0.5 + 0.5 * corr, 0.0, 1.0))


def pathway_identity_retention(final: np.ndarray, reference: np.ndarray, adjacency: np.ndarray) -> float:
    return float(np.clip(1.0 - pathway_distance(final, reference, adjacency), 0.0, 1.0))


def pathway_distance(field: np.ndarray, reference: np.ndarray, adjacency: np.ndarray) -> float:
    field_flux = pathway_flux(field, adjacency)
    ref_flux = pathway_flux(reference, adjacency)
    return float(np.linalg.norm(field_flux - ref_flux) / max(np.linalg.norm(ref_flux), EPS))


def fmo_identity_score(field: np.ndarray, reference: np.ndarray, adjacency: np.ndarray) -> float:
    address_distance = np.linalg.norm(local_address(field, adjacency) - local_address(reference, adjacency))
    address_norm = max(float(np.linalg.norm(local_address(reference, adjacency))), EPS)
    shell_distance = np.linalg.norm(local_shell_variance(field, adjacency) - local_shell_variance(reference, adjacency))
    shell_norm = max(float(np.linalg.norm(local_shell_variance(reference, adjacency))), EPS)
    path_distance = pathway_distance(field, reference, adjacency)
    return float(-0.40 * address_distance / address_norm - 0.35 * shell_distance / shell_norm - 0.25 * path_distance)


def run_null_mean(final: np.ndarray, reference: np.ndarray, adjacency: np.ndarray, config: FMOTelemetryConfig) -> float:
    rng = np.random.default_rng(config.seed + 99)
    buckets = degree_shell_strata(reference, adjacency, 3)
    values = [state_recoverability(final[stratified_indices(buckets, rng)], reference) for _ in range(8)]
    return float(np.mean(values))


def state_recoverability(state: np.ndarray, reference: np.ndarray) -> float:
    distance = float(np.linalg.norm(state - reference) / max(np.sqrt(reference.size), EPS))
    return float(np.clip(1.0 - distance, 0.0, 1.0))


def local_address(field: np.ndarray, adjacency: np.ndarray) -> np.ndarray:
    degree = np.maximum(adjacency.sum(axis=1), EPS)
    return (adjacency @ field) / degree - field


def local_shell_variance(field: np.ndarray, adjacency: np.ndarray) -> np.ndarray:
    degree = np.maximum(adjacency.sum(axis=1), EPS)
    offsets = field[None, :] - field[:, None]
    return np.sum(adjacency * offsets * offsets, axis=1) / degree


def shell_variance_restoration(field: np.ndarray, adjacency: np.ndarray, target: np.ndarray) -> np.ndarray:
    neighbor_mean = field + local_address(field, adjacency)
    current = local_shell_variance(field, adjacency)
    scale = np.maximum(np.abs(target) + np.abs(current), 1.0)
    return np.clip((target - current) / scale, -1.0, 1.0) * (field - neighbor_mean)


def spectral_basis(adjacency: np.ndarray, modes: int) -> np.ndarray:
    transition = build_transition(adjacency)
    laplacian = np.eye(adjacency.shape[0]) - 0.5 * (transition + transition.T)
    eigenvalues, eigenvectors = np.linalg.eigh(laplacian)
    order = np.argsort(eigenvalues)
    stop = min(max(2, modes + 1), eigenvectors.shape[1])
    return eigenvectors[:, order[1:stop]]


def spectral_field_signature(field: np.ndarray, adjacency: np.ndarray, modes: int) -> np.ndarray:
    basis = spectral_basis(adjacency, modes)
    coeffs = basis.T @ field
    energy = coeffs * coeffs
    return energy / max(float(np.sum(energy)), EPS)


def higher_order_signature(field: np.ndarray, adjacency: np.ndarray, modes: int) -> np.ndarray:
    return np.concatenate(
        [
            spectral_field_signature(field, adjacency, modes),
            np.asarray(
                [
                    local_autocorrelation(field, adjacency),
                    local_autocorrelation(field, two_hop_adjacency(adjacency)),
                    local_autocorrelation(field, triangle_supported_adjacency(adjacency)),
                ],
                dtype=float,
            ),
        ]
    )


def trajectory_signature(states: np.ndarray, adjacency: np.ndarray, modes: int) -> np.ndarray:
    delta = np.diff(states, axis=0)
    energy = np.mean(delta * delta, axis=0)
    lag = np.mean(states[1:] * states[:-1], axis=0)
    return np.concatenate(
        [
            higher_order_signature(states[-1], adjacency, modes),
            np.asarray(
                [
                    local_autocorrelation(energy, adjacency),
                    local_autocorrelation(lag, adjacency),
                    float(np.mean(energy)),
                    float(np.std(energy)),
                    float(np.mean(lag)),
                    float(np.std(lag)),
                ],
                dtype=float,
            ),
        ]
    )


def local_autocorrelation(field: np.ndarray, adjacency: np.ndarray) -> float:
    if not np.any(adjacency > 0.0):
        return 0.0
    centered = field - float(np.mean(field))
    numerator = float(np.sum(adjacency * centered[:, None] * centered[None, :]))
    denominator = float(np.sum(adjacency) * max(np.var(centered), EPS))
    return numerator / max(denominator, EPS)


def two_hop_adjacency(adjacency: np.ndarray) -> np.ndarray:
    binary = (adjacency > 0.0).astype(float)
    out = binary @ binary
    out[binary > 0.0] = 0.0
    np.fill_diagonal(out, 0.0)
    return normalize_adjacency(out)


def triangle_supported_adjacency(adjacency: np.ndarray) -> np.ndarray:
    binary = (adjacency > 0.0).astype(float)
    out = adjacency * (binary @ binary)
    np.fill_diagonal(out, 0.0)
    return normalize_adjacency(out)


def degree_shell_strata(reference: np.ndarray, adjacency: np.ndarray, bins: int) -> np.ndarray:
    degree = adjacency.sum(axis=1)
    shell = local_shell_variance(reference, adjacency)
    return quantile_bins(degree, bins) * bins + quantile_bins(shell, bins)


def stratified_indices(labels: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    indices = np.arange(labels.size)
    for label in np.unique(labels):
        idx = np.flatnonzero(labels == label)
        if idx.size > 1:
            indices[idx] = rng.permutation(idx)
    return indices


def quantile_bins(values: np.ndarray, bins: int) -> np.ndarray:
    edges = np.quantile(values, np.linspace(0.0, 1.0, bins + 1)[1:-1])
    return np.digitize(values, edges, right=True).astype(int)


def edge_shuffle(adjacency: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    upper_idx = np.triu_indices(adjacency.shape[0], k=1)
    upper = adjacency[upper_idx]
    shuffled = rng.permutation(upper)
    out = np.zeros_like(adjacency)
    out[upper_idx] = shuffled
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


def write_outputs(
    config: FMOTelemetryConfig,
    status: dict[str, Any],
    observed: FMORun,
    controls: list[FMORun],
    adjacency: np.ndarray,
    reference: np.ndarray,
) -> None:
    write_json(config.output_dir / "bridge_status.json", status)
    write_report(config.output_dir / "fmo_telemetry_report.md", status, observed, controls)
    write_csv(config.output_dir / "fmo_telemetry_summary.csv", [observed, *controls])
    write_plots(config, observed, controls, adjacency, reference)


def write_report(path: Path, status: dict[str, Any], observed: FMORun, controls: list[FMORun]) -> None:
    lines = [
        "# FMO Spectral Telemetry Report",
        "",
        f"- bridge_status: {status['bridge_status']}",
        f"- failure_mode: {status['failure_mode']}",
        f"- nodes: {status['nodes']}",
        f"- edges: {status['edges']}",
        f"- null_level: {status['null_level']}",
        f"- control_contrast: {status['control_contrast']}",
        f"- control_strict_pass_count: {status['control_strict_pass_count']}",
        "",
        "| run | recoverability | site_identity | pathway_identity | delta_persistence | safety_margin | active_null_z | strict_pass |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for run in [observed, *controls]:
        summary = run.summary
        lines.append(
            "| {run} | {rec:.6f} | {site:.6f} | {path:.6f} | {delta:.6f} | {margin:.6f} | {z:.6f} | {passed} |".format(
                run=run.label,
                rec=float(summary["recoverability_score"]),
                site=float(summary["site_identity_retention"]),
                path=float(summary["pathway_identity_retention"]),
                delta=float(summary["delta_persistence"]),
                margin=float(summary["safety_margin"]),
                z=float(summary["active_null_z"]),
                passed=summary["strict_pass"],
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_csv(path: Path, runs: list[FMORun]) -> None:
    fields = [
        "run",
        "recoverability_score",
        "site_identity_retention",
        "pathway_identity_retention",
        "delta_persistence",
        "safety_margin",
        "active_null_z",
        "strict_pass",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for run in runs:
            writer.writerow({field: run.summary[field] for field in fields})


def write_plots(config: FMOTelemetryConfig, observed: FMORun, controls: list[FMORun], adjacency: np.ndarray, reference: np.ndarray) -> None:
    import matplotlib.pyplot as plt

    x = np.arange(reference.size)
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(x, reference, marker="o", label="reference")
    ax.plot(x, observed.states[-1], marker="s", label="observed final")
    ax.set_xlabel("FMO site")
    ax.set_ylabel("normalized address field")
    ax.set_title("FMO spectral telemetry final site field")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(config.output_dir / "final_site_field.png", dpi=170)
    plt.close(fig)

    labels = [run.label for run in [observed, *controls]]
    z_values = [float(run.summary["active_null_z"]) for run in [observed, *controls]]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(labels, z_values)
    ax.set_ylabel("active null z")
    ax.set_title("FMO telemetry null comparison")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(config.output_dir / "null_comparison.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 4))
    im = ax.imshow(adjacency, cmap="magma")
    ax.set_title("FMO-like weighted interaction graph")
    ax.set_xlabel("site")
    ax.set_ylabel("site")
    fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    fig.savefig(config.output_dir / "fmo_weighted_graph.png", dpi=170)
    plt.close(fig)


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


def stable_label_offset(label: str) -> int:
    return sum((idx + 1) * ord(char) for idx, char in enumerate(label))
