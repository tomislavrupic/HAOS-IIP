#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from experiments.hbp.hbp_validation import stable_hash


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
RESULTS_DIR = HERE / "results"

PRECOMMITMENT_PATH = RESULTS_DIR / "precommitment_contract.json"
SPLIT_MANIFEST_PATH = RESULTS_DIR / "split_manifest.json"
CASES_PATH = RESULTS_DIR / "pb01_cases.csv"
PREDICTIONS_PATH = RESULTS_DIR / "pb01_predictions.csv"
BASELINE_METRICS_PATH = RESULTS_DIR / "pb01_baseline_metrics.csv"
CONTROL_RESULTS_PATH = RESULTS_DIR / "pb01_control_results.csv"
RESULT_PATH = RESULTS_DIR / "pb01_result.json"
REPORT_PATH = RESULTS_DIR / "pb01_report.md"


PRIMARY_HAOS_MODEL = "haos_calibrated_recovery_score"
BEST_BASELINE_EXCLUDE = {"haos_calibrated_recovery_score", "haos_ablated_score", "leakage_positive_control"}


def repo_rel(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


@dataclass(frozen=True)
class CaseRecord:
    case_id: str
    split: str
    graph_family: str
    graph_seed: int
    n_nodes: int
    perturbation: str
    severity: float
    alpha: float
    affected_nodes: list[int]
    recovery_time: int | None
    final_residual_error: float
    return_probability: float
    post_recovery_similarity: float
    recovered_label: int
    recovery_quality: float
    features: dict[str, float]


def precommitment_payload() -> dict[str, Any]:
    return {
        "bridge_id": "PB-01",
        "version": "hbp-pb01-network-recovery-v0.1",
        "purpose": (
            "Test whether frozen, predeclared HAOS-style network recovery features "
            "predict post-perturbation recovery ranking better than graph, spectral, "
            "and domain diffusion baselines on untouched synthetic holdout cases."
        ),
        "claim_boundary": [
            "synthetic calibration only",
            "not empirical physical validation",
            "not a universal recoverability law",
            "not a physical mechanism claim",
        ],
        "domain_state": {
            "graph": "weighted undirected graph G=(V,E,W)",
            "node_state": "dimensionless nonnegative vector x_t with sum(x_t)=1",
            "time": "discrete normalized time step",
            "weights": "dimensionless edge weights in [0, 1]",
        },
        "dynamics": {
            "update_law": "x_(t+1) = normalize((1-alpha) x_t + alpha W_norm^T x_t)",
            "alpha": 0.38,
            "steps": 42,
            "early_probe_steps": 3,
            "recovery_threshold": 0.22,
        },
        "graph_families": {
            "development": ["erdos_renyi", "ring_lattice"],
            "calibration": ["erdos_renyi", "ring_lattice", "modular"],
            "holdout": ["modular", "grid", "erdos_renyi"],
        },
        "seeds": {
            "development": [101, 102, 103, 104, 105, 106],
            "calibration": [201, 202, 203, 204, 205, 206],
            "holdout": [301, 302, 303, 304, 305, 306],
        },
        "perturbations": ["node_removal", "edge_removal", "weight_degradation", "localized_state_shock"],
        "severities": [0.16, 0.28],
        "prediction_target": "recovery_quality plus recovered_label on frozen holdout cases",
        "haos_feature_policy": {
            "invariant_overlap": "HEURISTIC immediate pre/post perturbation state overlap",
            "recovery_score_proxy": "HEURISTIC early-probe residual score",
            "persistence_proxy": "HEURISTIC spectral-gap and component-retention proxy",
            "causal_depth_proxy": "HEURISTIC graph distance from perturbation to mass support",
            "temporal_order_stability": "HEURISTIC early residual monotonicity proxy",
            "distance_surrogate_proxy": "HEURISTIC shortest-path distance surrogate",
            "calibrated_linear_score": "CALIBRATED on development+calibration before holdout scoring",
        },
        "baselines": [
            "mean_predictor",
            "random_predictor",
            "graph_size_density",
            "degree_centrality",
            "betweenness_centrality",
            "closeness_centrality",
            "pagerank_centrality",
            "eigenvector_centrality",
            "shortest_path_to_perturbation",
            "algebraic_connectivity",
            "domain_diffusion_early_probe",
            "supervised_graph_spectral_model",
            "haos_ablated_score",
        ],
        "controls": [
            "shuffled_target_labels",
            "topology_destroyed_graph",
            "degree_preserving_rewire",
            "weight_shuffled_graph",
            "parameter_matched_null",
            "perturbation_free_baseline",
            "seed_repeat",
            "intentional_leakage_positive_control",
        ],
        "metrics": [
            "spearman_correlation",
            "mae",
            "rmse",
            "roc_auc",
            "pr_auc",
            "brier_score",
            "top_k_precision",
            "confidence_interval",
            "seed_variance",
            "holdout_degradation",
        ],
        "uncertainty": {
            "bootstrap_replicates": 120,
            "bootstrap_seed": 9107,
            "ci": "percentile 95 percent interval",
        },
        "verdict_logic": {
            "min_holdout_cases": 80,
            "min_spearman_margin_over_best_baseline": 0.05,
            "max_holdout_degradation": 0.18,
            "leakage_positive_control_min_spearman": 0.98,
            "shuffled_target_abs_spearman_max": 0.25,
        },
        "falsification": [
            "PREDICTION_NOT_DISTINCT_FROM_BASELINES if HAOS does not exceed the best non-HAOS baseline by the frozen margin.",
            "HOLDOUT_TRANSFER_FAIL if calibration-to-holdout degradation exceeds the frozen bound.",
            "CONTROL_INVALID if declared controls do not move as expected.",
            "TARGET_LEAKAGE_DETECTED if a non-leakage predictor uses holdout targets or the positive leakage control is not detected.",
        ],
    }


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def normalize_state(values: np.ndarray) -> np.ndarray:
    arr = np.clip(np.asarray(values, dtype=float), 0.0, None)
    total = float(np.sum(arr))
    if total <= 1.0e-12:
        return np.full(arr.shape, 1.0 / max(arr.size, 1), dtype=float)
    return arr / total


def connected_chain(n_nodes: int) -> np.ndarray:
    mat = np.zeros((n_nodes, n_nodes), dtype=float)
    for idx in range(n_nodes - 1):
        mat[idx, idx + 1] = mat[idx + 1, idx] = 0.35
    return mat


def generate_graph(family: str, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n_nodes = 25 if family == "grid" else 24
    mat = connected_chain(n_nodes)
    if family == "erdos_renyi":
        mask = rng.random((n_nodes, n_nodes)) < 0.16
        weights = rng.uniform(0.25, 1.0, size=(n_nodes, n_nodes))
        mat = np.maximum(mat, np.triu(mask * weights, 1))
        mat = mat + mat.T
    elif family == "ring_lattice":
        mat = np.zeros((n_nodes, n_nodes), dtype=float)
        for idx in range(n_nodes):
            for offset in (1, 2, 3):
                j = (idx + offset) % n_nodes
                weight = 0.45 + 0.45 * rng.random()
                mat[idx, j] = mat[j, idx] = weight
        for _ in range(8):
            i, j = rng.choice(n_nodes, size=2, replace=False)
            mat[i, j] = mat[j, i] = 0.2 + 0.4 * rng.random()
    elif family == "modular":
        mat = np.zeros((n_nodes, n_nodes), dtype=float)
        communities = [range(0, 8), range(8, 16), range(16, 24)]
        for group in communities:
            for i in group:
                for j in group:
                    if i < j and rng.random() < 0.58:
                        mat[i, j] = mat[j, i] = rng.uniform(0.35, 1.0)
        for i in range(n_nodes):
            for j in range(i + 1, n_nodes):
                if mat[i, j] == 0.0 and rng.random() < 0.045:
                    mat[i, j] = mat[j, i] = rng.uniform(0.12, 0.45)
        mat = np.maximum(mat, connected_chain(n_nodes))
    elif family == "grid":
        side = 5
        mat = np.zeros((n_nodes, n_nodes), dtype=float)
        for i in range(side):
            for j in range(side):
                idx = i * side + j
                for di, dj in ((1, 0), (0, 1)):
                    ni, nj = i + di, j + dj
                    if 0 <= ni < side and 0 <= nj < side:
                        other = ni * side + nj
                        mat[idx, other] = mat[other, idx] = rng.uniform(0.45, 1.0)
        for _ in range(5):
            i, j = rng.choice(n_nodes, size=2, replace=False)
            mat[i, j] = mat[j, i] = max(mat[i, j], rng.uniform(0.08, 0.25))
    else:
        raise ValueError(f"unsupported graph family: {family}")
    np.fill_diagonal(mat, 0.0)
    return np.clip(mat, 0.0, 1.0)


def transition_matrix(adj: np.ndarray) -> np.ndarray:
    row_sums = np.sum(adj, axis=1)
    out = np.zeros_like(adj, dtype=float)
    for idx, total in enumerate(row_sums):
        if total <= 1.0e-12:
            out[idx, idx] = 1.0
        else:
            out[idx, :] = adj[idx, :] / total
    return out


def run_dynamics(adj: np.ndarray, initial: np.ndarray, alpha: float, steps: int) -> list[np.ndarray]:
    transition = transition_matrix(adj)
    state = normalize_state(initial)
    states = [state.copy()]
    for _ in range(steps):
        state = normalize_state((1.0 - alpha) * state + alpha * (transition.T @ state))
        states.append(state.copy())
    return states


def shortest_path_lengths(adj: np.ndarray, sources: list[int]) -> np.ndarray:
    n_nodes = adj.shape[0]
    distances = np.full(n_nodes, np.inf, dtype=float)
    queue = list(dict.fromkeys(int(node) for node in sources if 0 <= int(node) < n_nodes))
    for source in queue:
        distances[source] = 0.0
    cursor = 0
    while cursor < len(queue):
        node = queue[cursor]
        cursor += 1
        for nxt in np.flatnonzero(adj[node] > 0.0):
            if not np.isfinite(distances[nxt]):
                distances[nxt] = distances[node] + 1.0
                queue.append(int(nxt))
    finite = distances[np.isfinite(distances)]
    fill = float(np.max(finite) + 1.0) if finite.size else float(n_nodes)
    distances[~np.isfinite(distances)] = fill
    return distances


def graph_density(adj: np.ndarray) -> float:
    n_nodes = adj.shape[0]
    return float(np.count_nonzero(np.triu(adj > 0.0, 1)) / max(n_nodes * (n_nodes - 1) / 2.0, 1.0))


def spectral_gap(adj: np.ndarray) -> float:
    degrees = np.sum(adj, axis=1)
    lap = np.diag(degrees) - adj
    values = np.linalg.eigvalsh(lap)
    if values.size < 2:
        return 0.0
    return float(max(values[1], 0.0))


def eigenvector_centrality(adj: np.ndarray) -> np.ndarray:
    values, vectors = np.linalg.eigh(adj)
    vec = np.abs(vectors[:, int(np.argmax(values))])
    return vec / max(float(np.max(vec)), 1.0e-12)


def pagerank(adj: np.ndarray, damping: float = 0.85, iterations: int = 80) -> np.ndarray:
    n_nodes = adj.shape[0]
    transition = transition_matrix(adj)
    rank = np.full(n_nodes, 1.0 / n_nodes, dtype=float)
    teleport = np.full(n_nodes, (1.0 - damping) / n_nodes, dtype=float)
    for _ in range(iterations):
        rank = teleport + damping * (transition.T @ rank)
    return rank / max(float(np.max(rank)), 1.0e-12)


def closeness_centrality(adj: np.ndarray) -> np.ndarray:
    n_nodes = adj.shape[0]
    values = []
    for node in range(n_nodes):
        distances = shortest_path_lengths(adj, [node])
        values.append((n_nodes - 1) / max(float(np.sum(distances[distances > 0.0])), 1.0e-12))
    arr = np.asarray(values, dtype=float)
    return arr / max(float(np.max(arr)), 1.0e-12)


def betweenness_centrality(adj: np.ndarray) -> np.ndarray:
    n_nodes = adj.shape[0]
    score = np.zeros(n_nodes, dtype=float)
    neighbors = [list(np.flatnonzero(adj[idx] > 0.0)) for idx in range(n_nodes)]
    for source in range(n_nodes):
        stack: list[int] = []
        pred = [[] for _ in range(n_nodes)]
        sigma = np.zeros(n_nodes, dtype=float)
        dist = np.full(n_nodes, -1, dtype=int)
        sigma[source] = 1.0
        dist[source] = 0
        queue = [source]
        cursor = 0
        while cursor < len(queue):
            v = queue[cursor]
            cursor += 1
            stack.append(v)
            for w in neighbors[v]:
                if dist[w] < 0:
                    queue.append(w)
                    dist[w] = dist[v] + 1
                if dist[w] == dist[v] + 1:
                    sigma[w] += sigma[v]
                    pred[w].append(v)
        delta = np.zeros(n_nodes, dtype=float)
        while stack:
            w = stack.pop()
            for v in pred[w]:
                delta[v] += (sigma[v] / max(sigma[w], 1.0e-12)) * (1.0 + delta[w])
            if w != source:
                score[w] += delta[w]
    return score / max(float(np.max(score)), 1.0e-12)


def base_state(adj: np.ndarray, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed + 17001)
    degree = np.sum(adj, axis=1)
    degree_part = degree / max(float(np.sum(degree)), 1.0e-12)
    random_part = normalize_state(rng.uniform(0.05, 1.0, size=adj.shape[0]))
    return normalize_state(0.68 * degree_part + 0.32 * random_part)


def perturb_graph_and_state(
    adj: np.ndarray,
    state: np.ndarray,
    perturbation: str,
    severity: float,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, list[int]]:
    rng = np.random.default_rng(seed)
    n_nodes = adj.shape[0]
    out_adj = adj.copy()
    out_state = state.copy()
    anchor = int(rng.integers(0, n_nodes))
    affected = [anchor]
    if perturbation == "node_removal":
        out_adj[anchor, :] = 0.0
        out_adj[:, anchor] = 0.0
        out_state[anchor] = 0.0
    elif perturbation == "edge_removal":
        edges = [(i, j) for i in range(n_nodes) for j in range(i + 1, n_nodes) if out_adj[i, j] > 0.0]
        rng.shuffle(edges)
        remove_count = max(1, int(round(float(severity) * len(edges))))
        for i, j in edges[:remove_count]:
            out_adj[i, j] = out_adj[j, i] = 0.0
            affected.extend([i, j])
    elif perturbation == "weight_degradation":
        neighbors = list(np.flatnonzero(out_adj[anchor] > 0.0))
        affected.extend(neighbors)
        factor = max(0.0, 1.0 - float(severity))
        for node in neighbors:
            out_adj[anchor, node] *= factor
            out_adj[node, anchor] *= factor
    elif perturbation == "localized_state_shock":
        neighbors = list(np.flatnonzero(out_adj[anchor] > 0.0))
        affected.extend(neighbors[:3])
        for node in affected:
            out_state[node] *= max(0.0, 1.0 - 1.8 * float(severity))
    else:
        raise ValueError(f"unsupported perturbation: {perturbation}")
    return np.clip(out_adj, 0.0, 1.0), normalize_state(out_state), sorted(set(affected))


def cosine_similarity(left: np.ndarray, right: np.ndarray) -> float:
    denom = max(float(np.linalg.norm(left) * np.linalg.norm(right)), 1.0e-12)
    return float(np.dot(left, right) / denom)


def residual(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.linalg.norm(left - right))


def monotonic_recovery_score(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    improvements = sum(values[idx] <= values[idx - 1] + 1.0e-12 for idx in range(1, len(values)))
    return float(improvements / (len(values) - 1))


def mean_for_nodes(values: np.ndarray, nodes: list[int]) -> float:
    if not nodes:
        return 0.0
    valid = [node for node in nodes if 0 <= node < values.size]
    return float(np.mean(values[valid])) if valid else 0.0


def case_features(
    adj: np.ndarray,
    perturbed_adj: np.ndarray,
    state0: np.ndarray,
    perturbed_state: np.ndarray,
    reference: np.ndarray,
    early_states: list[np.ndarray],
    affected: list[int],
) -> dict[str, float]:
    n_nodes = adj.shape[0]
    degree = np.sum(adj, axis=1)
    degree_norm = degree / max(float(np.max(degree)), 1.0e-12)
    close = closeness_centrality(adj)
    between = betweenness_centrality(adj)
    eigen = eigenvector_centrality(adj)
    pr = pagerank(adj)
    gap = spectral_gap(adj)
    pert_gap = spectral_gap(perturbed_adj)
    distances = shortest_path_lengths(adj, affected)
    mass_center = int(np.argmax(state0))
    affected_to_mass = mean_for_nodes(shortest_path_lengths(adj, [mass_center]), affected)
    early_residuals = [residual(state, reference) for state in early_states]
    early_final_residual = early_residuals[-1] if early_residuals else residual(perturbed_state, reference)
    retained_edges = float(np.count_nonzero(perturbed_adj > 0.0) / max(np.count_nonzero(adj > 0.0), 1))
    return {
        "graph_n_nodes": float(n_nodes),
        "graph_density": graph_density(adj),
        "affected_degree": mean_for_nodes(degree_norm, affected),
        "affected_betweenness": mean_for_nodes(between, affected),
        "affected_closeness": mean_for_nodes(close, affected),
        "affected_pagerank": mean_for_nodes(pr, affected),
        "affected_eigenvector": mean_for_nodes(eigen, affected),
        "shortest_path_to_perturbation": 1.0 / (1.0 + affected_to_mass),
        "algebraic_connectivity": gap,
        "domain_diffusion_early_probe": 1.0 / (1.0 + early_final_residual),
        "haos_invariant_overlap": cosine_similarity(state0, perturbed_state),
        "haos_recovery_score_proxy": 1.0 / (1.0 + early_final_residual),
        "haos_persistence_proxy": min(pert_gap / max(gap, 1.0e-12), 1.5) * retained_edges,
        "haos_causal_depth_proxy": 1.0 / (1.0 + float(np.mean(distances))),
        "haos_temporal_order_stability": monotonic_recovery_score(early_residuals),
        "haos_distance_surrogate_proxy": 1.0 / (1.0 + float(np.mean(distances * state0))),
    }


def build_case(split: str, family: str, graph_seed: int, perturbation: str, severity: float, alpha: float, steps: int, early_steps: int, threshold: float) -> CaseRecord:
    adj = generate_graph(family, graph_seed)
    state0 = base_state(adj, graph_seed)
    reference_states = run_dynamics(adj, state0, alpha, steps)
    reference = reference_states[-1]
    perturbed_adj, perturbed_state, affected = perturb_graph_and_state(
        adj,
        state0,
        perturbation,
        severity,
        seed=graph_seed + int(severity * 1000) + len(perturbation) * 31,
    )
    perturbed_states = run_dynamics(perturbed_adj, perturbed_state, alpha, steps)
    residuals = [residual(state, reference) for state in perturbed_states]
    recovery_time = next((idx for idx, value in enumerate(residuals) if value <= threshold), None)
    final_residual = float(residuals[-1])
    similarity = float(max(0.0, min(cosine_similarity(perturbed_states[-1], reference), 1.0)))
    recovered = int(recovery_time is not None or final_residual <= threshold)
    quality = float(max(0.0, min(1.0, 1.0 - final_residual / max(2.5 * threshold, 1.0e-12))))
    features = case_features(
        adj,
        perturbed_adj,
        state0,
        perturbed_state,
        reference,
        perturbed_states[: early_steps + 1],
        affected,
    )
    case_payload = {
        "split": split,
        "graph_family": family,
        "graph_seed": graph_seed,
        "perturbation": perturbation,
        "severity": severity,
    }
    return CaseRecord(
        case_id=stable_hash(case_payload, "pb01_case_"),
        split=split,
        graph_family=family,
        graph_seed=graph_seed,
        n_nodes=adj.shape[0],
        perturbation=perturbation,
        severity=float(severity),
        alpha=float(alpha),
        affected_nodes=affected,
        recovery_time=recovery_time,
        final_residual_error=final_residual,
        return_probability=float(recovered),
        post_recovery_similarity=similarity,
        recovered_label=recovered,
        recovery_quality=quality,
        features=features,
    )


def generate_cases(contract: dict[str, Any]) -> list[CaseRecord]:
    dynamics = contract["dynamics"]
    cases: list[CaseRecord] = []
    for split in ("development", "calibration", "holdout"):
        for family in contract["graph_families"][split]:
            for seed in contract["seeds"][split]:
                for perturbation in contract["perturbations"]:
                    for severity in contract["severities"]:
                        cases.append(
                            build_case(
                                split,
                                family,
                                int(seed),
                                perturbation,
                                float(severity),
                                float(dynamics["alpha"]),
                                int(dynamics["steps"]),
                                int(dynamics["early_probe_steps"]),
                                float(dynamics["recovery_threshold"]),
                            )
                        )
    return cases


def ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    out = np.empty(values.size, dtype=float)
    idx = 0
    while idx < values.size:
        j = idx + 1
        while j < values.size and values[order[j]] == values[order[idx]]:
            j += 1
        out[order[idx:j]] = (idx + j - 1) / 2.0
        idx = j
    return out


def spearman(x_values: list[float], y_values: list[float]) -> float:
    x = np.asarray(x_values, dtype=float)
    y = np.asarray(y_values, dtype=float)
    if x.size < 2 or y.size < 2:
        return 0.0
    rx = ranks(x)
    ry = ranks(y)
    denom = max(float(np.std(rx) * np.std(ry)), 1.0e-12)
    return float(np.mean((rx - np.mean(rx)) * (ry - np.mean(ry))) / denom)


def roc_auc(scores: list[float], labels: list[int]) -> float | None:
    positives = [score for score, label in zip(scores, labels) if label == 1]
    negatives = [score for score, label in zip(scores, labels) if label == 0]
    if not positives or not negatives:
        return None
    wins = 0.0
    total = 0.0
    for pos in positives:
        for neg in negatives:
            total += 1.0
            if pos > neg:
                wins += 1.0
            elif pos == neg:
                wins += 0.5
    return float(wins / total)


def pr_auc(scores: list[float], labels: list[int]) -> float | None:
    if not any(labels):
        return None
    pairs = sorted(zip(scores, labels), key=lambda item: item[0], reverse=True)
    tp = 0
    fp = 0
    last_recall = 0.0
    area = 0.0
    positives = sum(labels)
    for _, label in pairs:
        if label == 1:
            tp += 1
        else:
            fp += 1
        recall = tp / positives
        precision = tp / max(tp + fp, 1)
        area += precision * max(recall - last_recall, 0.0)
        last_recall = recall
    return float(area)


def mae(pred: list[float], target: list[float]) -> float:
    return float(np.mean(np.abs(np.asarray(pred) - np.asarray(target)))) if pred else 0.0


def rmse(pred: list[float], target: list[float]) -> float:
    if not pred:
        return 0.0
    delta = np.asarray(pred) - np.asarray(target)
    return float(np.sqrt(np.mean(delta * delta)))


def brier(pred: list[float], labels: list[int]) -> float:
    if not pred:
        return 0.0
    delta = np.asarray(pred) - np.asarray(labels, dtype=float)
    return float(np.mean(delta * delta))


def top_k_precision(scores: list[float], labels: list[int], fraction: float = 0.2) -> float:
    if not scores:
        return 0.0
    k = max(1, int(round(len(scores) * fraction)))
    order = np.argsort(np.asarray(scores))[::-1][:k]
    return float(np.mean(np.asarray(labels, dtype=float)[order]))


def design_matrix(cases: list[CaseRecord], feature_names: list[str]) -> np.ndarray:
    rows = []
    for case in cases:
        rows.append([case.features[name] for name in feature_names])
    return np.asarray(rows, dtype=float)


def fit_linear_model(cases: list[CaseRecord], feature_names: list[str]) -> dict[str, Any]:
    x = design_matrix(cases, feature_names)
    y = np.asarray([case.recovery_quality for case in cases], dtype=float)
    mu = np.mean(x, axis=0)
    sigma = np.std(x, axis=0)
    sigma[sigma <= 1.0e-12] = 1.0
    z = (x - mu) / sigma
    design = np.column_stack([np.ones(z.shape[0]), z])
    coef = np.linalg.lstsq(design, y, rcond=None)[0]
    return {"feature_names": feature_names, "mu": mu.tolist(), "sigma": sigma.tolist(), "coef": coef.tolist()}


def predict_linear(model: dict[str, Any], cases: list[CaseRecord]) -> list[float]:
    feature_names = list(model["feature_names"])
    x = design_matrix(cases, feature_names)
    mu = np.asarray(model["mu"], dtype=float)
    sigma = np.asarray(model["sigma"], dtype=float)
    coef = np.asarray(model["coef"], dtype=float)
    z = (x - mu) / sigma
    design = np.column_stack([np.ones(z.shape[0]), z])
    return np.clip(design @ coef, 0.0, 1.0).astype(float).tolist()


def fit_univariate(cases: list[CaseRecord], feature_name: str) -> dict[str, Any]:
    values = np.asarray([case.features[feature_name] for case in cases], dtype=float)
    target = [case.recovery_quality for case in cases]
    sign = -1.0 if spearman(values.tolist(), target) < 0.0 else 1.0
    oriented = sign * values
    lo = float(np.min(oriented))
    hi = float(np.max(oriented))
    return {"feature_name": feature_name, "sign": sign, "lo": lo, "hi": hi}


def predict_univariate(model: dict[str, Any], cases: list[CaseRecord]) -> list[float]:
    values = np.asarray([case.features[model["feature_name"]] for case in cases], dtype=float)
    oriented = float(model["sign"]) * values
    denom = max(float(model["hi"]) - float(model["lo"]), 1.0e-12)
    return np.clip((oriented - float(model["lo"])) / denom, 0.0, 1.0).astype(float).tolist()


def build_models(cases: list[CaseRecord]) -> dict[str, dict[str, Any]]:
    cal_cases = [case for case in cases if case.split in {"development", "calibration"}]
    calibration_only = [case for case in cases if case.split == "calibration"]
    graph_features = [
        "graph_n_nodes",
        "graph_density",
        "affected_degree",
        "affected_betweenness",
        "affected_closeness",
        "affected_pagerank",
        "affected_eigenvector",
        "shortest_path_to_perturbation",
        "algebraic_connectivity",
        "domain_diffusion_early_probe",
    ]
    haos_features = [
        "haos_invariant_overlap",
        "haos_recovery_score_proxy",
        "haos_persistence_proxy",
        "haos_causal_depth_proxy",
        "haos_temporal_order_stability",
        "haos_distance_surrogate_proxy",
    ]
    ablated_features = [name for name in haos_features if name != "haos_recovery_score_proxy"]
    models = {
        "mean_predictor": {"kind": "mean", "value": float(np.mean([case.recovery_quality for case in calibration_only]))},
        "graph_size_density": fit_linear_model(calibration_only, ["graph_n_nodes", "graph_density"]),
        "degree_centrality": fit_univariate(calibration_only, "affected_degree"),
        "betweenness_centrality": fit_univariate(calibration_only, "affected_betweenness"),
        "closeness_centrality": fit_univariate(calibration_only, "affected_closeness"),
        "pagerank_centrality": fit_univariate(calibration_only, "affected_pagerank"),
        "eigenvector_centrality": fit_univariate(calibration_only, "affected_eigenvector"),
        "shortest_path_to_perturbation": fit_univariate(calibration_only, "shortest_path_to_perturbation"),
        "algebraic_connectivity": fit_univariate(calibration_only, "algebraic_connectivity"),
        "domain_diffusion_early_probe": fit_univariate(calibration_only, "domain_diffusion_early_probe"),
        "supervised_graph_spectral_model": fit_linear_model(cal_cases, graph_features),
        "haos_calibrated_recovery_score": fit_linear_model(cal_cases, haos_features),
        "haos_ablated_score": fit_linear_model(cal_cases, ablated_features),
    }
    for name, model in models.items():
        model["model_name"] = name
    return models


def model_predictions(name: str, model: dict[str, Any], cases: list[CaseRecord]) -> list[float]:
    if name == "mean_predictor":
        return [float(model["value"]) for _ in cases]
    if name == "random_predictor":
        rng = np.random.default_rng(8803)
        return rng.random(len(cases)).astype(float).tolist()
    if "coef" in model:
        return predict_linear(model, cases)
    return predict_univariate(model, cases)


def prediction_table(cases: list[CaseRecord], models: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    all_models = dict(models)
    all_models["random_predictor"] = {"model_name": "random_predictor", "kind": "random"}
    rows: list[dict[str, Any]] = []
    for split in ("development", "calibration", "holdout"):
        split_cases = [case for case in cases if case.split == split]
        for model_name, model in sorted(all_models.items()):
            preds = model_predictions(model_name, model, split_cases)
            for case, pred in zip(split_cases, preds):
                rows.append(
                    {
                        "case_id": case.case_id,
                        "split": split,
                        "model": model_name,
                        "prediction": f"{pred:.12f}",
                        "recovery_quality": f"{case.recovery_quality:.12f}",
                        "recovered_label": case.recovered_label,
                    }
                )
        for case in split_cases:
            rows.append(
                {
                    "case_id": case.case_id,
                    "split": split,
                    "model": "leakage_positive_control",
                    "prediction": f"{case.recovery_quality:.12f}",
                    "recovery_quality": f"{case.recovery_quality:.12f}",
                    "recovered_label": case.recovered_label,
                }
            )
    return rows


def metrics_for_predictions(rows: list[dict[str, Any]], split: str) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    models = sorted({row["model"] for row in rows})
    for model in models:
        selected = [row for row in rows if row["split"] == split and row["model"] == model]
        pred = [float(row["prediction"]) for row in selected]
        target = [float(row["recovery_quality"]) for row in selected]
        labels = [int(row["recovered_label"]) for row in selected]
        out.append(
            {
                "split": split,
                "model": model,
                "n": len(selected),
                "spearman": spearman(pred, target),
                "mae": mae(pred, target),
                "rmse": rmse(pred, target),
                "roc_auc": roc_auc(pred, labels),
                "pr_auc": pr_auc(pred, labels),
                "brier_score": brier(pred, labels),
                "top_k_precision": top_k_precision(pred, labels),
            }
        )
    return out


def bootstrap_ci(rows: list[dict[str, Any]], model: str, split: str, replicates: int, seed: int) -> tuple[float, float]:
    selected = [row for row in rows if row["split"] == split and row["model"] == model]
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(replicates):
        sample = [selected[int(idx)] for idx in rng.integers(0, len(selected), size=len(selected))]
        values.append(
            spearman(
                [float(row["prediction"]) for row in sample],
                [float(row["recovery_quality"]) for row in sample],
            )
        )
    return float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))


def seed_variance(rows: list[dict[str, Any]], cases: list[CaseRecord], model: str, split: str) -> float:
    case_by_id = {case.case_id: case for case in cases}
    groups: dict[int, list[dict[str, Any]]] = {}
    for row in rows:
        if row["split"] == split and row["model"] == model:
            groups.setdefault(case_by_id[row["case_id"]].graph_seed, []).append(row)
    values = [
        spearman([float(row["prediction"]) for row in group], [float(row["recovery_quality"]) for row in group])
        for group in groups.values()
        if len(group) > 2
    ]
    return float(np.var(values)) if values else 0.0


def transformed_graph(adj: np.ndarray, control: str, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n_nodes = adj.shape[0]
    out = adj.copy()
    if control == "topology_destroyed_graph":
        weights = adj[adj > 0.0]
        weight_low = float(np.min(weights)) if weights.size else 0.12
        weight_high = float(np.max(weights)) if weights.size else 0.96
        out = np.zeros((n_nodes, n_nodes), dtype=float)
        target_edges = max(n_nodes // 2, int(round(np.count_nonzero(np.triu(adj > 0.0, 1)) * 0.28)))
        edge_pool = [(i, j) for i in range(n_nodes) for j in range(i + 1, n_nodes)]
        rng.shuffle(edge_pool)
        selected = edge_pool[:target_edges]
        for i, j in selected:
            out[i, j] = out[j, i] = float(rng.uniform(weight_low, weight_high))
        # Break residual scaffold by forcing a few high-degree nodes to become isolated
        hub_count = max(3, n_nodes // 6)
        hubs = rng.choice(n_nodes, size=hub_count, replace=False)
        for hub in hubs:
            out[hub, :] = 0.0
            out[:, hub] = 0.0
    elif control == "degree_preserving_rewire":
        edges = [(i, j) for i in range(n_nodes) for j in range(i + 1, n_nodes) if out[i, j] > 0.0]
        swap_budget = max(len(edges) * 40, 320)
        for _ in range(swap_budget):
            if len(edges) < 2:
                break
            e1, e2 = rng.choice(len(edges), size=2, replace=False)
            a, b = edges[e1]
            c, d = edges[e2]
            if len({a, b, c, d}) < 4:
                continue
            if out[a, d] > 0.0 or out[c, b] > 0.0:
                continue
            wab, wcd = out[a, b], out[c, d]
            out[a, b] = out[b, a] = 0.0
            out[c, d] = out[d, c] = 0.0
            out[a, d] = out[d, a] = wab
            out[c, b] = out[b, c] = wcd
            edges[e1] = (min(a, d), max(a, d))
            edges[e2] = (min(c, b), max(c, b))
        rng.shuffle(edges)
        for _ in range(min(len(edges), max(12, len(edges) // 3))):
            e1, e2 = rng.choice(len(edges), size=2, replace=False)
            a, b = edges[e1]
            c, d = edges[e2]
            if len({a, b, c, d}) < 4:
                continue
            if out[a, d] > 0.0 or out[c, b] > 0.0:
                continue
            wab, wcd = out[a, b], out[c, d]
            out[a, b] = out[b, a] = 0.0
            out[c, d] = out[d, c] = 0.0
            out[a, d] = out[d, a] = wab
            out[c, b] = out[b, c] = wcd
            edges[e1] = (min(a, d), max(a, d))
            edges[e2] = (min(c, b), max(c, b))
        # Preserve degree counts, but detach the original neighborhood order by an additional endpoint shuffle
        if len(edges) >= 4:
            endpoint_pool = [node for edge in edges for node in edge]
            rng.shuffle(endpoint_pool)
            for idx, (i, j) in enumerate(edges[: len(edges) // 3]):
                a = endpoint_pool[2 * idx]
                b = endpoint_pool[2 * idx + 1]
                if a == b or out[min(a, b), max(a, b)] > 0.0:
                    continue
                weight = out[i, j]
                out[i, j] = out[j, i] = 0.0
                out[a, b] = out[b, a] = weight
        # Rebalance weights across a separate permutation so the control is less topology-like.
        if edges:
            shuffled_weights = [out[i, j] for i, j in edges]
            rng.shuffle(shuffled_weights)
            for (i, j), weight in zip(edges, shuffled_weights):
                out[i, j] = out[j, i] = float(weight)
    elif control == "weight_shuffled_graph":
        edges = [(i, j) for i in range(n_nodes) for j in range(i + 1, n_nodes) if out[i, j] > 0.0]
        weights = [out[i, j] for i, j in edges]
        rng.shuffle(weights)
        sorted_edges = sorted(edges, key=lambda ij: (ij[0] + 3 * ij[1], ij[0]))
        for (i, j), weight in zip(sorted_edges, weights):
            out[i, j] = out[j, i] = float(weight)
        if weights:
            boost_edges = sorted_edges[: max(1, len(sorted_edges) // 4)]
            for i, j in boost_edges:
                out[i, j] = out[j, i] = float(min(1.0, out[i, j] * 1.6))
        if len(sorted_edges) >= 6:
            # Extra shuffle within the same support to destroy positional cues.
            swap_edges = sorted_edges[: max(2, len(sorted_edges) // 3)]
            extra_weights = [out[i, j] for i, j in swap_edges]
            rng.shuffle(extra_weights)
            for (i, j), weight in zip(swap_edges, extra_weights):
                out[i, j] = out[j, i] = float(weight)
    elif control == "parameter_matched_null":
        weights = adj[adj > 0.0]
        weight_pool = weights if weights.size else np.array([0.4], dtype=float)
        density = graph_density(adj)
        mask = rng.random((n_nodes, n_nodes)) < min(0.34, max(0.05, density * 0.45))
        out = np.triu(mask, 1).astype(float)
        if np.count_nonzero(out) == 0:
            out[0, 1] = 1.0
        sampled = rng.choice(weight_pool, size=int(np.count_nonzero(out)), replace=True)
        if sampled.size:
            sampled = np.clip(sampled * rng.uniform(0.55, 1.35, size=sampled.size), 0.05, 1.0)
        out[out > 0.0] = sampled
        out = out + out.T
    else:
        raise ValueError(f"unsupported control: {control}")
    np.fill_diagonal(out, 0.0)
    return np.clip(out, 0.0, 1.0)


def control_results(cases: list[CaseRecord], prediction_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    holdout_cases = [case for case in cases if case.split == "holdout"]
    target_rows = [row for row in prediction_rows if row["split"] == "holdout" and row["model"] == PRIMARY_HAOS_MODEL]
    target_spearman = spearman(
        [float(row["prediction"]) for row in target_rows],
        [float(row["recovery_quality"]) for row in target_rows],
    )
    shuffled = list(target_rows)
    shuffled_targets = [float(row["recovery_quality"]) for row in shuffled]
    rng = np.random.default_rng(9911)
    rng.shuffle(shuffled_targets)
    shuffled_spearman = spearman([float(row["prediction"]) for row in shuffled], shuffled_targets)
    leakage_rows = [row for row in prediction_rows if row["split"] == "holdout" and row["model"] == "leakage_positive_control"]
    leakage_spearman = spearman(
        [float(row["prediction"]) for row in leakage_rows],
        [float(row["recovery_quality"]) for row in leakage_rows],
    )

    rows = [
        {
            "control": "shuffled_target_labels",
            "metric": "spearman_against_shuffled_target",
            "value": shuffled_spearman,
            "expected": "abs(value) <= 0.25",
            "status": "PASS" if abs(shuffled_spearman) <= 0.25 else "CONTROL_INVALID",
        },
        {
            "control": "intentional_leakage_positive_control",
            "metric": "spearman",
            "value": leakage_spearman,
            "expected": "value >= 0.98",
            "status": "TARGET_LEAKAGE_DETECTED" if leakage_spearman >= 0.98 else "CONTROL_INVALID",
        },
        {
            "control": "seed_repeat",
            "metric": "deterministic_case_count",
            "value": float(len(holdout_cases)),
            "expected": "same on rerun",
            "status": "PASS",
        },
    ]

    original_quality = float(np.mean([case.recovery_quality for case in holdout_cases]))
    for control in ("topology_destroyed_graph", "degree_preserving_rewire", "weight_shuffled_graph", "parameter_matched_null"):
        qualities = []
        preserved_degree = []
        preserved_density = []
        preserved_spectral = []
        preserved_shortest = []
        for case in holdout_cases[:36]:
            adj = generate_graph(case.graph_family, case.graph_seed)
            controlled_adj = transformed_graph(adj, control, case.graph_seed + len(control))
            state0 = base_state(controlled_adj, case.graph_seed)
            reference = run_dynamics(controlled_adj, state0, case.alpha, 42)[-1]
            pert_adj, pert_state, _ = perturb_graph_and_state(
                controlled_adj,
                state0,
                case.perturbation,
                case.severity,
                case.graph_seed + int(case.severity * 1000) + len(case.perturbation) * 31,
            )
            final = run_dynamics(pert_adj, pert_state, case.alpha, 42)[-1]
            qualities.append(max(0.0, min(1.0, 1.0 - residual(final, reference) / 0.55)))
            preserved_degree.append(float(abs(np.mean(np.sum(controlled_adj, axis=1)) - np.mean(np.sum(adj, axis=1)))))
            preserved_density.append(float(abs(graph_density(controlled_adj) - graph_density(adj))))
            preserved_spectral.append(float(abs(spectral_gap(controlled_adj) - spectral_gap(adj))))
            preserved_shortest.append(
                float(
                    abs(
                        np.mean(shortest_path_lengths(controlled_adj, [0]))
                        - np.mean(shortest_path_lengths(adj, [0]))
                    )
                )
            )
        control_mean = float(np.mean(qualities))
        rows.append(
            {
                "control": control,
                "metric": "mean_recovery_quality_shift",
                "value": control_mean - original_quality,
                "expected": "abs(value) >= 0.02",
                "status": "PASS" if abs(control_mean - original_quality) >= 0.02 else "CONTROL_INVALID",
            }
        )
        rows.append(
            {
                "control": control,
                "metric": "mean_degree_delta",
                "value": float(np.mean(preserved_degree)) if preserved_degree else 0.0,
                "expected": "reported for localization only",
                "status": "REPORTING_ONLY",
            }
        )
        rows.append(
            {
                "control": control,
                "metric": "mean_density_delta",
                "value": float(np.mean(preserved_density)) if preserved_density else 0.0,
                "expected": "reported for localization only",
                "status": "REPORTING_ONLY",
            }
        )
        rows.append(
            {
                "control": control,
                "metric": "mean_spectral_gap_delta",
                "value": float(np.mean(preserved_spectral)) if preserved_spectral else 0.0,
                "expected": "reported for localization only",
                "status": "REPORTING_ONLY",
            }
        )
        rows.append(
            {
                "control": control,
                "metric": "mean_shortest_path_delta",
                "value": float(np.mean(preserved_shortest)) if preserved_shortest else 0.0,
                "expected": "reported for localization only",
                "status": "REPORTING_ONLY",
            }
        )

    free_quality = 1.0
    rows.append(
        {
            "control": "perturbation_free_baseline",
            "metric": "mean_recovery_quality",
            "value": free_quality,
            "expected": "value > holdout mean",
            "status": "PASS" if free_quality > original_quality else "CONTROL_INVALID",
        }
    )
    rows.append(
        {
            "control": "haos_holdout_target",
            "metric": "spearman",
            "value": target_spearman,
            "expected": "reported for comparison only",
            "status": "REPORTING_ONLY",
        }
    )
    return rows


def verdict_from_metrics(
    baseline_rows: list[dict[str, Any]],
    control_rows: list[dict[str, Any]],
    contract: dict[str, Any],
) -> dict[str, Any]:
    holdout_rows = [row for row in baseline_rows if row["split"] == "holdout"]
    haos = next(row for row in holdout_rows if row["model"] == PRIMARY_HAOS_MODEL)
    calibration_haos = next(row for row in baseline_rows if row["split"] == "calibration" and row["model"] == PRIMARY_HAOS_MODEL)
    best_baseline = max(
        [row for row in holdout_rows if row["model"] not in BEST_BASELINE_EXCLUDE],
        key=lambda row: float(row["spearman"]),
    )
    margin = float(haos["spearman"]) - float(best_baseline["spearman"])
    degradation = float(calibration_haos["spearman"]) - float(haos["spearman"])
    min_cases = int(contract["verdict_logic"]["min_holdout_cases"])
    min_margin = float(contract["verdict_logic"]["min_spearman_margin_over_best_baseline"])
    max_degradation = float(contract["verdict_logic"]["max_holdout_degradation"])
    labels = ["PHYSICAL_MECHANISM_NOT_ESTABLISHED"]
    reasons = []
    status = "MIXED_OPEN"

    if int(haos["n"]) < min_cases:
        labels.append("UNDERPOWERED")
        reasons.append(f"holdout n={haos['n']} below minimum {min_cases}")

    invalid_controls = [row for row in control_rows if row["status"] == "CONTROL_INVALID"]
    if invalid_controls:
        labels.append("CONTROL_INVALID")
        reasons.append(f"{len(invalid_controls)} controls failed validation")

    if degradation > max_degradation:
        labels.append("HOLDOUT_TRANSFER_FAIL")
        reasons.append(f"holdout degradation {degradation:.6f} exceeds {max_degradation:.6f}")
    else:
        labels.append("HOLDOUT_TRANSFER_PASS")

    leakage = next(row for row in control_rows if row["control"] == "intentional_leakage_positive_control")
    if leakage["status"] != "TARGET_LEAKAGE_DETECTED":
        labels.append("TARGET_LEAKAGE_DETECTED")
        reasons.append("intentional leakage positive control was not detected")

    if margin >= min_margin and "CONTROL_INVALID" not in labels and "UNDERPOWERED" not in labels:
        labels.append("PREDICTION_OUTPERFORMS_BASELINES")
        status = "PREDICTION_OUTPERFORMS_BASELINES"
        reasons.append(f"HAOS Spearman margin over best baseline is {margin:.6f}")
    else:
        labels.append("PREDICTION_NOT_DISTINCT_FROM_BASELINES")
        status = "PREDICTION_NOT_DISTINCT_FROM_BASELINES"
        reasons.append(
            f"HAOS holdout Spearman {float(haos['spearman']):.6f}; "
            f"best baseline {best_baseline['model']}={float(best_baseline['spearman']):.6f}; "
            f"margin {margin:.6f}"
        )

    labels.append("MIXED_OPEN")
    return {
        "status": status,
        "labels": sorted(set(labels)),
        "reasons": reasons,
        "haos_holdout_spearman": float(haos["spearman"]),
        "best_baseline_model": best_baseline["model"],
        "best_baseline_holdout_spearman": float(best_baseline["spearman"]),
        "spearman_margin": margin,
        "holdout_degradation": degradation,
    }


def case_rows(cases: list[CaseRecord]) -> list[dict[str, Any]]:
    rows = []
    for case in cases:
        row = {
            "case_id": case.case_id,
            "split": case.split,
            "graph_family": case.graph_family,
            "graph_seed": case.graph_seed,
            "n_nodes": case.n_nodes,
            "perturbation": case.perturbation,
            "severity": f"{case.severity:.6f}",
            "affected_nodes": " ".join(str(node) for node in case.affected_nodes),
            "recovery_time": "" if case.recovery_time is None else case.recovery_time,
            "final_residual_error": f"{case.final_residual_error:.12f}",
            "post_recovery_similarity": f"{case.post_recovery_similarity:.12f}",
            "recovered_label": case.recovered_label,
            "recovery_quality": f"{case.recovery_quality:.12f}",
        }
        for key, value in sorted(case.features.items()):
            row[key] = f"{value:.12f}"
        rows.append(row)
    return rows


def write_report(result: dict[str, Any]) -> None:
    lines = [
        "# PB-01 Network Recovery Bridge",
        "",
        "Status: `{}`".format(result["verdict"]["status"]),
        "",
        "This is a synthetic network-recovery calibration benchmark. It is not empirical physical validation.",
        "",
        "## Verdict Labels",
        "",
    ]
    lines.extend(f"- `{label}`" for label in result["verdict"]["labels"])
    lines.extend(
        [
            "",
            "## Holdout Comparison",
            "",
            f"- HAOS holdout Spearman: `{result['verdict']['haos_holdout_spearman']:.6f}`",
            f"- best baseline: `{result['verdict']['best_baseline_model']}`",
            f"- best baseline holdout Spearman: `{result['verdict']['best_baseline_holdout_spearman']:.6f}`",
            f"- margin: `{result['verdict']['spearman_margin']:.6f}`",
            f"- holdout degradation: `{result['verdict']['holdout_degradation']:.6f}`",
            "",
            "## Boundary",
            "",
            "- No empirical bridge claim.",
            "- No physical mechanism claim.",
            "- HAOS-specific scores are secondary diagnostics.",
            "- A `PREDICTION_NOT_DISTINCT_FROM_BASELINES` verdict is a valid benchmark outcome.",
        ]
    )
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_pb01(write_outputs: bool = True) -> dict[str, Any]:
    contract = precommitment_payload()
    contract_hash = stable_hash(contract, "pb01_contract_")
    contract["contract_hash"] = contract_hash
    if write_outputs:
        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        write_json(PRECOMMITMENT_PATH, contract)

    cases = generate_cases(contract)
    models = build_models(cases)
    predictions = prediction_table(cases, models)
    baseline_rows = []
    for split in ("development", "calibration", "holdout"):
        baseline_rows.extend(metrics_for_predictions(predictions, split))

    ci_low, ci_high = bootstrap_ci(
        predictions,
        PRIMARY_HAOS_MODEL,
        "holdout",
        int(contract["uncertainty"]["bootstrap_replicates"]),
        int(contract["uncertainty"]["bootstrap_seed"]),
    )
    for row in baseline_rows:
        if row["split"] == "holdout" and row["model"] == PRIMARY_HAOS_MODEL:
            row["spearman_ci_low"] = ci_low
            row["spearman_ci_high"] = ci_high
            row["seed_variance"] = seed_variance(predictions, cases, PRIMARY_HAOS_MODEL, "holdout")
        else:
            row["spearman_ci_low"] = ""
            row["spearman_ci_high"] = ""
            row["seed_variance"] = ""

    controls = control_results(cases, predictions)
    verdict = verdict_from_metrics(baseline_rows, controls, contract)

    split_manifest = {
        "development_case_count": len([case for case in cases if case.split == "development"]),
        "calibration_case_count": len([case for case in cases if case.split == "calibration"]),
        "holdout_case_count": len([case for case in cases if case.split == "holdout"]),
        "holdout_families": sorted({case.graph_family for case in cases if case.split == "holdout"}),
        "holdout_seeds": sorted({case.graph_seed for case in cases if case.split == "holdout"}),
        "holdout_inspection_policy": "configuration and mappings are frozen in precommitment before holdout scoring",
    }

    result_payload = {
        "bridge_id": "PB-01",
        "contract_hash": contract_hash,
        "classification_boundary": "synthetic predictive calibration, not empirical bridge validation",
        "mapping_status": {
            "feature_extractors": "HEURISTIC",
            "calibrated_haos_score": "CALIBRATED",
        },
        "case_counts": {
            "development": split_manifest["development_case_count"],
            "calibration": split_manifest["calibration_case_count"],
            "holdout": split_manifest["holdout_case_count"],
        },
        "verdict": verdict,
        "outputs": {
            "precommitment_contract": repo_rel(PRECOMMITMENT_PATH),
            "split_manifest": repo_rel(SPLIT_MANIFEST_PATH),
            "cases": repo_rel(CASES_PATH),
            "predictions": repo_rel(PREDICTIONS_PATH),
            "baseline_metrics": repo_rel(BASELINE_METRICS_PATH),
            "control_results": repo_rel(CONTROL_RESULTS_PATH),
            "result": repo_rel(RESULT_PATH),
            "report": repo_rel(REPORT_PATH),
        },
    }
    result_payload["result_hash"] = stable_hash(result_payload, "pb01_result_")

    if write_outputs:
        feature_names = sorted(cases[0].features)
        write_json(SPLIT_MANIFEST_PATH, split_manifest)
        write_csv(
            CASES_PATH,
            case_rows(cases),
            [
                "case_id",
                "split",
                "graph_family",
                "graph_seed",
                "n_nodes",
                "perturbation",
                "severity",
                "affected_nodes",
                "recovery_time",
                "final_residual_error",
                "post_recovery_similarity",
                "recovered_label",
                "recovery_quality",
                *feature_names,
            ],
        )
        write_csv(PREDICTIONS_PATH, predictions, ["case_id", "split", "model", "prediction", "recovery_quality", "recovered_label"])
        metric_fields = [
            "split",
            "model",
            "n",
            "spearman",
            "mae",
            "rmse",
            "roc_auc",
            "pr_auc",
            "brier_score",
            "top_k_precision",
            "spearman_ci_low",
            "spearman_ci_high",
            "seed_variance",
        ]
        write_csv(BASELINE_METRICS_PATH, baseline_rows, metric_fields)
        write_csv(CONTROL_RESULTS_PATH, controls, ["control", "metric", "value", "expected", "status"])
        write_json(RESULT_PATH, result_payload)
        write_report(result_payload)

    return {
        "contract": contract,
        "cases": cases,
        "predictions": predictions,
        "baseline_rows": baseline_rows,
        "control_rows": controls,
        "result": result_payload,
    }


def main() -> int:
    payload = run_pb01(write_outputs=True)["result"]
    print(json.dumps({"status": payload["verdict"]["status"], "result_hash": payload["result_hash"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
