from __future__ import annotations

from dataclasses import dataclass
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from microtubule_lattice_model import build_microtubule_lattice, build_weighted_adjacency


ROOT = Path(__file__).resolve().parent
EPS = 1.0e-12


@dataclass(frozen=True)
class SpectralTelemetryConfig:
    output_dir: Path = ROOT / "spectral_outputs"
    seed: int = 20260429
    protofilaments: int = 13
    dimers_per_protofilament: int = 30
    steps: int = 96
    dt: float = 0.08
    diffusion: float = 0.12
    address_gain: float = 0.42
    invariant_gain: float = 0.025
    spectral_modes: int = 10
    perturbation_step: int = 18
    thermal_noise: float = 0.18
    damage_scale: float = 0.72
    permutation_trials: int = 32
    null_candidates: int = 4
    null_level: int = 4


@dataclass(frozen=True)
class TelemetryRun:
    label: str
    states: np.ndarray
    summary: dict[str, Any]


def run_microtubule_spectral_telemetry(config: SpectralTelemetryConfig) -> dict[str, Any]:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(config.seed)
    lattice = build_microtubule_lattice(dimers_per_protofilament=config.dimers_per_protofilament)
    adjacency = normalize_adjacency(build_weighted_adjacency(lattice))
    reference = microtubule_reference_field(config)
    damage_nodes = choose_cap_damage_nodes(config)

    observed = run_single("observed", adjacency, reference, damage_nodes, config, rng)
    controls = [
        run_single("node_shuffle_control", adjacency, reference, damage_nodes, config, rng, initial_override=rng.permutation(reference)),
        run_single("edge_shuffle_control", edge_shuffle(adjacency, rng), reference, damage_nodes, config, rng),
        run_single("topology_randomization_control", topology_randomize(adjacency, rng), reference, damage_nodes, config, rng),
    ]
    status = classify(observed, controls, config, adjacency)
    write_outputs(config, status, observed, controls)
    return status


def run_single(
    label: str,
    adjacency: np.ndarray,
    reference: np.ndarray,
    damage_nodes: np.ndarray,
    config: SpectralTelemetryConfig,
    rng: np.random.Generator,
    initial_override: np.ndarray | None = None,
) -> TelemetryRun:
    basis = spectral_basis(adjacency, config.spectral_modes)
    target_coeffs = basis.T @ reference
    target_shell = local_shell_variance(reference, adjacency)
    transition = build_transition(adjacency)
    laplacian = np.eye(adjacency.shape[0]) - transition
    initial = reference + 0.03 * rng.normal(size=reference.shape) if initial_override is None else initial_override
    states = np.zeros((config.steps + 1, reference.size), dtype=float)
    states[0] = normalize_field(initial, reference)

    for step in range(1, config.steps + 1):
        prior = states[step - 1].copy()
        if step == config.perturbation_step:
            prior[damage_nodes] += config.damage_scale * rng.normal(size=damage_nodes.size)
        if config.thermal_noise > 0.0:
            prior += config.thermal_noise * 0.01 * rng.normal(size=prior.shape)
        spectral_pull = basis @ (target_coeffs - basis.T @ prior)
        invariant_pull = shell_variance_restoration(prior, adjacency, target_shell)
        update = prior - config.dt * config.diffusion * (laplacian @ prior)
        update += config.dt * config.address_gain * spectral_pull
        update += config.dt * config.invariant_gain * invariant_pull
        states[step] = normalize_field(update, reference)

    summary = summarize_run(label, states, reference, adjacency, config)
    return TelemetryRun(label, states, summary)


def microtubule_reference_field(config: SpectralTelemetryConfig) -> np.ndarray:
    values = np.zeros(config.protofilaments * config.dimers_per_protofilament, dtype=float)
    for p in range(config.protofilaments):
        theta = 2.0 * math.pi * p / config.protofilaments
        for z in range(config.dimers_per_protofilament):
            axial = z / max(config.dimers_per_protofilament - 1, 1)
            cap = math.exp(-((1.0 - axial) ** 2) / 0.018)
            helical = math.sin(theta + 2.0 * math.pi * axial)
            lateral = 0.35 * math.cos(2.0 * theta)
            values[p * config.dimers_per_protofilament + z] = helical + lateral + 0.7 * cap
    return normalize_field(values, values)


def choose_cap_damage_nodes(config: SpectralTelemetryConfig) -> np.ndarray:
    nodes: list[int] = []
    start = max(0, config.dimers_per_protofilament - 5)
    for p in range(3, 7):
        for z in range(start, config.dimers_per_protofilament):
            nodes.append((p % config.protofilaments) * config.dimers_per_protofilament + z)
    return np.asarray(nodes, dtype=int)


def summarize_run(
    label: str,
    states: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    config: SpectralTelemetryConfig,
) -> dict[str, Any]:
    final = states[-1]
    recoverability = state_recoverability(final, reference)
    branch_retention = protofilament_identity_retention(final, reference, config)
    persistence = float(recoverability - state_recoverability(states[config.perturbation_step], reference))
    safety_margin = float(recoverability - max(run_null_mean(final, reference, adjacency, config), 0.0))
    nulls = null_ladder(final, states, reference, adjacency, config, config.seed + stable_label_offset(label))
    return {
        "run": label,
        "recoverability_score": recoverability,
        "protofilament_identity_retention": branch_retention,
        "delta_persistence": persistence,
        "safety_margin": safety_margin,
        **nulls,
        "strict_pass": bool(recoverability >= 0.72 and nulls["active_null_pass"]),
    }


def null_ladder(
    final: np.ndarray,
    states: np.ndarray,
    reference: np.ndarray,
    adjacency: np.ndarray,
    config: SpectralTelemetryConfig,
    seed: int,
) -> dict[str, float | bool]:
    rng = np.random.default_rng(seed)
    observed = branch_identity_score(final, reference, adjacency)
    buckets = degree_shell_strata(reference, adjacency, 4)
    scores = {"branch": [], "spectral": [], "higher": [], "trajectory": []}
    signatures = {
        "spectral": spectral_field_signature(final, adjacency, config.spectral_modes),
        "higher": higher_order_signature(final, adjacency, config.spectral_modes),
        "trajectory": trajectory_signature(states, adjacency, config.spectral_modes),
    }
    for _ in range(config.permutation_trials):
        candidates = [stratified_indices(buckets, rng) for _ in range(max(1, config.null_candidates))]
        branch_idx = candidates[0]
        scores["branch"].append(branch_identity_score(final[branch_idx], reference, adjacency))
        for key in ("spectral", "higher", "trajectory"):
            best = min(
                candidates,
                key=lambda idx: float(np.linalg.norm(signature_for(key, states[:, idx], final[idx], adjacency, config) - signatures[key])),
            )
            scores[key].append(branch_identity_score(final[best], reference, adjacency))

    out: dict[str, float | bool] = {}
    active = "branch"
    if config.null_level >= 2:
        active = "spectral"
    if config.null_level >= 4:
        active = "higher"
    if config.null_level >= 5:
        active = "trajectory"
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


def signature_for(
    key: str,
    states: np.ndarray,
    final: np.ndarray,
    adjacency: np.ndarray,
    config: SpectralTelemetryConfig,
) -> np.ndarray:
    if key == "spectral":
        return spectral_field_signature(final, adjacency, config.spectral_modes)
    if key == "higher":
        return higher_order_signature(final, adjacency, config.spectral_modes)
    if key == "trajectory":
        return trajectory_signature(states, adjacency, config.spectral_modes)
    raise ValueError(key)


def classify(observed: TelemetryRun, controls: list[TelemetryRun], config: SpectralTelemetryConfig, adjacency: np.ndarray) -> dict[str, Any]:
    best_control = max(float(run.summary["recoverability_score"]) for run in controls)
    control_contrast = float(observed.summary["recoverability_score"]) >= best_control + 0.05
    control_pass_count = sum(1 for run in controls if bool(run.summary["strict_pass"]))
    if bool(observed.summary["strict_pass"]) and control_contrast and control_pass_count == 0:
        status = "PASS"
        failure = "RECOVERABLE_MICROTUBULE_SPECTRAL_TELEMETRY"
    elif bool(observed.summary["strict_pass"]):
        status = "MARGINAL"
        failure = "BIOLOGY_CONTROL_MATCH"
    else:
        status = "FAIL"
        failure = "BIOLOGY_TELEMETRY_SPECIFICITY_FAILED"
    return {
        "bridge_status": status,
        "failure_mode": failure,
        "nodes": int(adjacency.shape[0]),
        "edges": int(np.count_nonzero(np.triu(adjacency, k=1))),
        "null_level": int(config.null_level),
        "control_contrast": control_contrast,
        "control_strict_pass_count": control_pass_count,
        "best_control_recoverability_score": best_control,
        "observed_summary": observed.summary,
        "control_summaries": {run.label: run.summary for run in controls},
    }


def run_null_mean(final: np.ndarray, reference: np.ndarray, adjacency: np.ndarray, config: SpectralTelemetryConfig) -> float:
    rng = np.random.default_rng(config.seed + 99)
    buckets = degree_shell_strata(reference, adjacency, 4)
    values = [branch_identity_score(final[stratified_indices(buckets, rng)], reference, adjacency) for _ in range(8)]
    return float(np.mean(values))


def protofilament_identity_retention(final: np.ndarray, reference: np.ndarray, config: SpectralTelemetryConfig) -> float:
    final_profile = final.reshape(config.protofilaments, config.dimers_per_protofilament).mean(axis=1)
    ref_profile = reference.reshape(config.protofilaments, config.dimers_per_protofilament).mean(axis=1)
    distance = float(np.linalg.norm(final_profile - ref_profile) / max(np.linalg.norm(ref_profile), EPS))
    return float(np.clip(1.0 - distance, 0.0, 1.0))


def branch_identity_score(field: np.ndarray, reference: np.ndarray, adjacency: np.ndarray) -> float:
    address_distance = np.linalg.norm(local_address(field, adjacency) - local_address(reference, adjacency))
    address_norm = max(float(np.linalg.norm(local_address(reference, adjacency))), EPS)
    shell_distance = np.linalg.norm(local_shell_variance(field, adjacency) - local_shell_variance(reference, adjacency))
    shell_norm = max(float(np.linalg.norm(local_shell_variance(reference, adjacency))), EPS)
    return float(-0.5 * (address_distance / address_norm + shell_distance / shell_norm))


def state_recoverability(state: np.ndarray, reference: np.ndarray) -> float:
    distance = float(np.linalg.norm(state - reference) / max(math.sqrt(reference.size), EPS))
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
            np.array(
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
            np.array(
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


def write_outputs(config: SpectralTelemetryConfig, status: dict[str, Any], observed: TelemetryRun, controls: list[TelemetryRun]) -> None:
    write_json(config.output_dir / "bridge_status.json", status)
    write_report(config.output_dir / "spectral_telemetry_report.md", status, observed, controls)
    write_csv(config.output_dir / "spectral_telemetry_summary.csv", [observed, *controls])
    write_plots(config, observed, controls)


def write_report(path: Path, status: dict[str, Any], observed: TelemetryRun, controls: list[TelemetryRun]) -> None:
    lines = [
        "# Microtubule Spectral Telemetry Report",
        "",
        f"- bridge_status: {status['bridge_status']}",
        f"- failure_mode: {status['failure_mode']}",
        f"- nodes: {status['nodes']}",
        f"- edges: {status['edges']}",
        f"- null_level: {status['null_level']}",
        f"- control_contrast: {status['control_contrast']}",
        f"- control_strict_pass_count: {status['control_strict_pass_count']}",
        "",
        "| run | recoverability | protofilament_identity | delta_persistence | safety_margin | active_null_z | strict_pass |",
        "| --- | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for run in [observed, *controls]:
        summary = run.summary
        lines.append(
            "| {run} | {rec:.6f} | {pid:.6f} | {delta:.6f} | {margin:.6f} | {z:.6f} | {passed} |".format(
                run=run.label,
                rec=float(summary["recoverability_score"]),
                pid=float(summary["protofilament_identity_retention"]),
                delta=float(summary["delta_persistence"]),
                margin=float(summary["safety_margin"]),
                z=float(summary["active_null_z"]),
                passed=summary["strict_pass"],
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_csv(path: Path, runs: list[TelemetryRun]) -> None:
    fields = ["run", "recoverability_score", "protofilament_identity_retention", "delta_persistence", "safety_margin", "active_null_z", "strict_pass"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for run in runs:
            writer.writerow({field: run.summary[field] for field in fields})


def write_plots(config: SpectralTelemetryConfig, observed: TelemetryRun, controls: list[TelemetryRun]) -> None:
    import matplotlib.pyplot as plt

    final = observed.states[-1].reshape(config.protofilaments, config.dimers_per_protofilament)
    fig, ax = plt.subplots(figsize=(8, 4))
    im = ax.imshow(final, aspect="auto", cmap="viridis")
    ax.set_title("Spectral telemetry final field")
    ax.set_xlabel("dimer index")
    ax.set_ylabel("protofilament")
    fig.colorbar(im, ax=ax, fraction=0.03)
    fig.tight_layout()
    fig.savefig(config.output_dir / "final_field_heatmap.png", dpi=160)
    plt.close(fig)

    labels = [run.label for run in [observed, *controls]]
    z_values = [float(run.summary["active_null_z"]) for run in [observed, *controls]]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(labels, z_values)
    ax.set_ylabel("active null z")
    ax.set_title("Microtubule telemetry null comparison")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(config.output_dir / "null_comparison.png", dpi=160)
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
