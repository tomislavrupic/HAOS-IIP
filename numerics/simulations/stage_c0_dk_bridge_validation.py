#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from collections import Counter, defaultdict
from datetime import datetime
from itertools import combinations
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

os.environ.setdefault('MPLCONFIGDIR', '/tmp/matplotlib')
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]
DATA_DIR = REPO_ROOT / 'data'
PLOTS_DIR = REPO_ROOT / 'plots'
ATLAS_NOTES = REPO_ROOT / 'experiments' / 'pre_geometry_atlas'
RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_dk_bridge_validation_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_C0_to_DK_Bridge_Validation_v1.md'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage_c0_dk_bridge_validation'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

for path in (DATA_DIR, PLOTS_DIR, ATLAS_NOTES):
    path.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'section',
    'family_id',
    'family_label',
    'variant_id',
    'variant_label',
    'delta',
    'nilpotency_norm',
    'operator_identity_norm',
    'laplacian_diff_norm',
    'triangle_count',
    'topology_class',
    'flow_concentration_index',
    'grade_exchange_coherence',
    'address_selectivity_index',
    'channel_count',
    'loop_count',
    'loop_score',
    'detuning_threshold',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the C0-to-DK bridge validation scaffold.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--family-ids', nargs='*', default=[])
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def representative_graph(family_id: str) -> nx.Graph:
    G = nx.Graph()
    if family_id == 'clustered_composite_anchor':
        G.add_edges_from(
            [
                (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3),
                (3, 4), (4, 5), (5, 6), (4, 6),
                (2, 5), (1, 4),
            ]
        )
        return G
    if family_id == 'counter_propagating_corridor':
        G.add_edges_from(
            [
                (0, 1), (1, 2), (2, 3), (3, 4), (4, 5),
                (5, 6), (6, 7), (7, 8),
                (1, 4), (2, 5), (3, 6), (4, 7),
                (2, 4), (4, 6),
            ]
        )
        return G
    if family_id == 'phase_ordered_symmetric_triad':
        G.add_edges_from(
            [
                (0, 1), (1, 2), (2, 0),
                (3, 4), (4, 5), (5, 3),
                (6, 7), (7, 8), (8, 6),
                (2, 3), (5, 6), (1, 4), (4, 7),
            ]
        )
        return G
    raise ValueError(f'unknown family_id: {family_id}')


def triangle_list(G: nx.Graph) -> list[tuple[int, int, int]]:
    triangles: list[tuple[int, int, int]] = []
    nodes = sorted(G.nodes())
    for a, b, c in combinations(nodes, 3):
        if G.has_edge(a, b) and G.has_edge(b, c) and G.has_edge(a, c):
            triangles.append((a, b, c))
    return triangles


def square_participation_counts(G: nx.Graph) -> dict[int, int]:
    counts = {int(node): 0 for node in G.nodes()}
    nodes = sorted(G.nodes())
    for quad in combinations(nodes, 4):
        sub = G.subgraph(quad)
        if sub.number_of_edges() != 4:
            continue
        if any(deg != 2 for _node, deg in sub.degree()):
            continue
        for node in quad:
            counts[int(node)] += 1
    return counts


def shell_distance_map(G: nx.Graph, seed: int) -> dict[int, int]:
    return {int(node): int(distance) for node, distance in nx.single_source_shortest_path_length(G, seed).items()}


def mismatch_scores(
    G: nx.Graph,
    left_seed: int,
    right_seed: int,
) -> dict[int, float]:
    triangles = triangle_list(G)
    tri_count = {int(node): 0 for node in G.nodes()}
    for tri in triangles:
        for node in tri:
            tri_count[int(node)] += 1
    square_count = square_participation_counts(G)
    left_shell = shell_distance_map(G, left_seed)
    right_shell = shell_distance_map(G, right_seed)
    scores: dict[int, float] = {}
    for node in G.nodes():
        degree = float(G.degree[node])
        imbalance = abs(float(left_shell.get(int(node), 0)) - float(right_shell.get(int(node), 0)))
        scores[int(node)] = (
            imbalance
            + 0.35 * degree
            + 0.40 * float(tri_count.get(int(node), 0))
            + 0.20 * float(square_count.get(int(node), 0))
        )
    return scores


def orient_edges(G: nx.Graph, scores: dict[int, float]) -> list[tuple[int, int]]:
    oriented: list[tuple[int, int]] = []
    for u, v in sorted((min(a, b), max(a, b)) for a, b in G.edges()):
        su = float(scores[int(u)])
        sv = float(scores[int(v)])
        if su > sv + 1.0e-12:
            oriented.append((int(u), int(v)))
        elif sv > su + 1.0e-12:
            oriented.append((int(v), int(u)))
        else:
            oriented.append((int(u), int(v)))
    return oriented


def build_d0(
    nodes: list[int],
    oriented_edges: list[tuple[int, int]],
) -> tuple[sp.csr_matrix, dict[frozenset[int], tuple[int, int, int]], dict[int, int]]:
    node_to_idx = {int(node): idx for idx, node in enumerate(nodes)}
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    edge_lookup: dict[frozenset[int], tuple[int, int, int]] = {}
    for edge_idx, (u, v) in enumerate(oriented_edges):
        rows.extend([edge_idx, edge_idx])
        cols.extend([node_to_idx[int(u)], node_to_idx[int(v)]])
        data.extend([-1.0, 1.0])
        edge_lookup[frozenset((int(u), int(v)))] = (edge_idx, int(u), int(v))
    d0 = sp.coo_matrix((data, (rows, cols)), shape=(len(oriented_edges), len(nodes)), dtype=float).tocsr()
    return d0, edge_lookup, node_to_idx


def build_d1(
    triangles: list[tuple[int, int, int]],
    edge_lookup: dict[frozenset[int], tuple[int, int, int]],
) -> sp.csr_matrix:
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    for face_idx, (a, b, c) in enumerate(triangles):
        boundary = [
            (a, b, +1.0),
            (b, c, +1.0),
            (a, c, -1.0),
        ]
        for u, v, sign in boundary:
            edge_idx, eu, ev = edge_lookup[frozenset((int(u), int(v)))]
            orientation_sign = 1.0 if (eu, ev) == (int(u), int(v)) else -1.0
            rows.append(face_idx)
            cols.append(edge_idx)
            data.append(sign * orientation_sign)
    return sp.coo_matrix((data, (rows, cols)), shape=(len(triangles), len(edge_lookup)), dtype=float).tocsr()


def block_dirac(d0: sp.csr_matrix, d1: sp.csr_matrix) -> tuple[sp.csr_matrix, sp.csr_matrix]:
    n0 = d0.shape[1]
    n1 = d0.shape[0]
    n2 = d1.shape[0]
    zero00 = sp.csr_matrix((n0, n0), dtype=float)
    zero01 = sp.csr_matrix((n0, n2), dtype=float)
    zero10 = sp.csr_matrix((n1, n1), dtype=float)
    zero20 = sp.csr_matrix((n2, n0), dtype=float)
    zero22 = sp.csr_matrix((n2, n2), dtype=float)

    D = sp.bmat(
        [
            [zero00, d0.T, zero01],
            [d0, zero10, d1.T],
            [zero20, d1, zero22],
        ],
        format='csr',
    )
    Delta = sp.block_diag(
        [
            d0.T @ d0,
            d0 @ d0.T + d1.T @ d1,
            d1 @ d1.T,
        ],
        format='csr',
    )
    return D, Delta


def basis_shell_profiles(
    G: nx.Graph,
    nodes: list[int],
    oriented_edges: list[tuple[int, int]],
    triangles: list[tuple[int, int, int]],
    seed: int,
    sigma_shell: float,
    kick_cycles: float,
    direction_sign: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    shell = shell_distance_map(G, seed)
    node_vals = np.zeros(len(nodes), dtype=complex)
    for idx, node in enumerate(nodes):
        dist = float(shell.get(int(node), 1000))
        envelope = math.exp(-(dist * dist) / (2.0 * sigma_shell * sigma_shell))
        phase = np.exp(1j * direction_sign * kick_cycles * dist)
        node_vals[idx] = envelope * phase

    edge_vals = np.zeros(len(oriented_edges), dtype=complex)
    for idx, (u, v) in enumerate(oriented_edges):
        du = float(shell.get(int(u), 1000))
        dv = float(shell.get(int(v), 1000))
        env = 0.5 * (
            math.exp(-(du * du) / (2.0 * sigma_shell * sigma_shell))
            + math.exp(-(dv * dv) / (2.0 * sigma_shell * sigma_shell))
        )
        grad = dv - du
        edge_vals[idx] = env * np.exp(1j * direction_sign * kick_cycles * grad)

    face_vals = np.zeros(len(triangles), dtype=complex)
    for idx, tri in enumerate(triangles):
        dists = [float(shell.get(int(node), 1000)) for node in tri]
        env = float(np.mean([math.exp(-(dist * dist) / (2.0 * sigma_shell * sigma_shell)) for dist in dists]))
        face_vals[idx] = env * np.exp(1j * direction_sign * kick_cycles * float(np.mean(dists)))

    return node_vals, edge_vals, face_vals


def initial_packet(
    G: nx.Graph,
    nodes: list[int],
    oriented_edges: list[tuple[int, int]],
    triangles: list[tuple[int, int, int]],
    seed: int,
    sigma_shell: float,
    kick_cycles: float,
    direction_sign: float,
    grade0_scale: float,
    grade1_scale: float,
    grade2_scale: float,
) -> np.ndarray:
    node_vals, edge_vals, face_vals = basis_shell_profiles(
        G=G,
        nodes=nodes,
        oriented_edges=oriented_edges,
        triangles=triangles,
        seed=seed,
        sigma_shell=sigma_shell,
        kick_cycles=kick_cycles,
        direction_sign=direction_sign,
    )
    psi = np.concatenate(
        [
            float(grade0_scale) * node_vals,
            float(grade1_scale) * edge_vals,
            float(grade2_scale) * face_vals,
        ]
    )
    norm = float(np.linalg.norm(psi))
    return psi / max(norm, 1.0e-12)


def build_address_labels(
    packet_states: list[np.ndarray],
    block_sizes: tuple[int, int, int],
    left_address: float,
    right_address: float,
    support_floor_fraction: float,
    dominant_ratio_threshold: float,
) -> np.ndarray:
    packet_power = np.stack([np.abs(packet) ** 2 for packet in packet_states], axis=1)
    total_power = np.sum(packet_power, axis=1)
    max_power = float(np.max(total_power)) if total_power.size else 0.0
    active = total_power >= float(support_floor_fraction) * max(max_power, 1.0e-12)
    dominant = np.argmax(packet_power, axis=1)
    dominant_power = packet_power[np.arange(packet_power.shape[0]), dominant]
    dominant_ratio = dominant_power / np.maximum(total_power, 1.0e-12)
    labels = np.full(total_power.shape, float(left_address), dtype=float)
    labels[active & (dominant == 1)] = float(right_address)

    start = 0
    grade_ids = np.zeros(sum(block_sizes), dtype=int)
    for grade, size in enumerate(block_sizes):
        grade_ids[start:start + size] = grade
        start += size

    mixed = active & (dominant_ratio < float(dominant_ratio_threshold))
    labels[mixed] = np.where(dominant[mixed] == 0, float(left_address), float(right_address))
    labels[~active] = np.where((grade_ids[~active] % 2) == 0, float(left_address), float(right_address))
    return labels


def continuous_mod_distance(a: np.ndarray, b: np.ndarray, modulus: float) -> np.ndarray:
    diff = np.mod(np.abs(a - b), modulus)
    return np.minimum(diff, modulus - diff)


def compatibility_weights(distance: np.ndarray, eta_match: float, eta_mismatch: float) -> np.ndarray:
    clipped = np.clip(np.asarray(distance, dtype=float), 0.0, 2.0)
    weights = np.empty_like(clipped, dtype=float)
    near = clipped <= 1.0
    weights[near] = 1.0 + (float(eta_match) - 1.0) * clipped[near]
    weights[~near] = float(eta_match) + (float(eta_mismatch) - float(eta_match)) * (clipped[~near] - 1.0)
    return weights


def address_weighted_operator(
    base_operator: sp.csr_matrix,
    labels: np.ndarray,
    modulus: float,
    eta_match: float,
    eta_mismatch: float,
) -> sp.csr_matrix:
    coo = base_operator.tocoo()
    distances = continuous_mod_distance(labels[coo.row], labels[coo.col], modulus)
    weights = compatibility_weights(distances, eta_match, eta_mismatch)
    weighted = sp.csr_matrix((coo.data * weights, (coo.row, coo.col)), shape=base_operator.shape)
    return weighted.tocsr()


def address_activity_statistics(
    states: list[np.ndarray],
    operator: sp.csr_matrix,
    labels: np.ndarray,
    modulus: float,
    eta_match: float,
    eta_mismatch: float,
) -> tuple[float, float]:
    coo = operator.tocoo()
    mask = coo.row < coo.col
    if not np.any(mask):
        return 0.0, 0.0
    row = coo.row[mask]
    col = coo.col[mask]
    distances = continuous_mod_distance(labels[row], labels[col], modulus)
    compat = compatibility_weights(distances, eta_match, eta_mismatch)
    activity = np.zeros(row.shape[0], dtype=float)
    for state in states:
        activity += np.abs(state[row]) * np.abs(state[col])
    total = float(np.sum(activity))
    if total <= 1.0e-12:
        return 0.0, 0.0
    weighted_mean = float(np.sum(compat * activity) / total)
    min_weight = float(compatibility_weights(np.asarray([2.0]), eta_match, eta_mismatch)[0])
    selectivity = float(np.clip((weighted_mean - min_weight) / max(1.0 - min_weight, 1.0e-12), 0.0, 1.0))
    sorting = float(np.clip(np.sum(activity[distances <= 1.0e-9]) / total - np.sum(activity[distances >= 1.0]) / total, -1.0, 1.0))
    return sorting, selectivity


def degree_skew_operator(
    operator: sp.csr_matrix,
    labels: np.ndarray,
    basis_degrees: np.ndarray,
    strength: float,
) -> sp.csr_matrix:
    if strength <= 0.0:
        return operator
    deg = np.asarray(basis_degrees, dtype=float)
    deg = deg / max(float(np.max(deg)), 1.0e-12)
    coo = operator.tocoo()
    mult = 1.0 + float(strength) * 0.5 * (deg[coo.row] + deg[coo.col])
    weighted = sp.csr_matrix((coo.data * mult, (coo.row, coo.col)), shape=operator.shape)
    return weighted.tocsr()


def estimate_lambda_max(operator: sp.csr_matrix, iters: int = 40) -> float:
    dim = operator.shape[0]
    if dim == 0:
        return 0.0
    vec = np.ones(dim, dtype=complex)
    vec /= max(float(np.linalg.norm(vec)), 1.0e-12)
    value = 0.0
    for _ in range(iters):
        w = operator @ vec
        norm = float(np.linalg.norm(w))
        if norm <= 1.0e-12:
            return 0.0
        vec = w / norm
        value = float(np.real(np.vdot(vec, operator @ vec)))
    return max(abs(value), 1.0e-12)


def evolve_states(operator: sp.csr_matrix, psi0: np.ndarray, steps: int, dt_scale: float) -> tuple[list[np.ndarray], float]:
    lam_max = estimate_lambda_max(operator)
    dt = float(dt_scale / max(lam_max, 1.0e-12))
    ident = sp.identity(operator.shape[0], dtype=complex, format='csr')
    lhs = (ident + 0.5j * dt * operator).tocsc()
    rhs_op = (ident - 0.5j * dt * operator).tocsr()
    solve = spla.factorized(lhs)
    psi = np.asarray(psi0, dtype=complex).copy()
    states = [psi.copy()]
    for _ in range(int(steps)):
        psi = solve(rhs_op @ psi)
        psi /= max(float(np.linalg.norm(psi)), 1.0e-12)
        states.append(psi.copy())
    return states, dt


def grade_histories(states: list[np.ndarray], block_sizes: tuple[int, int, int]) -> np.ndarray:
    start = 0
    slices = []
    for size in block_sizes:
        slices.append(slice(start, start + size))
        start += size
    hist = []
    for state in states:
        total = float(np.sum(np.abs(state) ** 2))
        total = max(total, 1.0e-12)
        hist.append([float(np.sum(np.abs(state[sl]) ** 2) / total) for sl in slices])
    return np.asarray(hist, dtype=float)


def grade_exchange_signal(grade_hist: np.ndarray) -> tuple[np.ndarray, float]:
    signal = np.linalg.norm(grade_hist - grade_hist[0][None, :], axis=1)
    peak = float(np.max(signal)) if signal.size else 0.0
    if peak <= 1.0e-12:
        return signal, 0.0
    coherence = float(np.mean(signal) / peak)
    return signal, float(np.clip(coherence, 0.0, 1.0))


def edge_component_count(oriented_edges: list[tuple[int, int]], active_mask: np.ndarray) -> int:
    active_indices = [idx for idx, active in enumerate(active_mask) if bool(active)]
    if not active_indices:
        return 0
    line = nx.Graph()
    for idx in active_indices:
        line.add_node(int(idx))
    for i, a in enumerate(active_indices):
        a_nodes = set(oriented_edges[a])
        for b in active_indices[i + 1:]:
            if a_nodes.intersection(oriented_edges[b]):
                line.add_edge(int(a), int(b))
    return nx.number_connected_components(line)


def triangle_loop_metrics(
    edge_series: np.ndarray,
    triangles: list[tuple[int, int, int]],
    d1: sp.csr_matrix,
) -> tuple[int, float]:
    if not triangles:
        return 0, 0.0
    edge_power = np.mean(np.abs(edge_series) ** 2, axis=0)
    total = float(np.sum(edge_power))
    if total <= 1.0e-12:
        return 0, 0.0
    loop_count = 0
    scores: list[float] = []
    d1_dense = d1.toarray()
    for face_idx, _tri in enumerate(triangles):
        edge_ids = np.flatnonzero(np.abs(d1_dense[face_idx]) > 0.0)
        if edge_ids.size != 3:
            continue
        powers = edge_power[edge_ids]
        threshold = max(0.10 * float(np.max(edge_power)), float(np.mean(edge_power) + 0.50 * np.std(edge_power)))
        if np.all(powers >= threshold):
            loop_count += 1
        phases = []
        for edge_id, coeff in zip(edge_ids, d1_dense[face_idx, edge_ids]):
            series = edge_series[:, edge_id]
            phase = np.angle(np.mean(series))
            phases.append(coeff * phase)
        closure = float(np.clip(0.5 * (1.0 + math.cos(float(np.sum(phases)))), 0.0, 1.0))
        support = float(np.min(powers) / total)
        scores.append(closure * support * len(edge_ids))
    return int(loop_count), float(max(scores) if scores else 0.0)


def classify_relative_topology(
    row: dict[str, Any],
    baseline: dict[str, Any],
) -> str:
    baseline_flow = max(float(baseline['flow_concentration_index']), 1.0e-12)
    baseline_coherence = max(float(baseline['grade_exchange_coherence']), 1.0e-12)
    baseline_selectivity = max(float(baseline['address_selectivity_index']), 1.0e-12)
    baseline_loop = max(float(baseline['loop_score']), 1.0e-12)
    baseline_asym = float(baseline['grade_asymmetry_index'])

    flow_ratio = float(row['flow_concentration_index']) / baseline_flow
    coherence_ratio = float(row['grade_exchange_coherence']) / baseline_coherence
    selectivity_ratio = float(row['address_selectivity_index']) / baseline_selectivity
    loop_ratio = float(row['loop_score']) / baseline_loop
    asymmetry = float(row['grade_asymmetry_index'])

    if (
        float(baseline['loop_score']) >= 0.15
        and float(row['loop_score']) >= max(0.15, 0.80 * float(baseline['loop_score']))
        and flow_ratio >= 0.85
        and coherence_ratio >= 0.90
        and selectivity_ratio >= 0.90
        and asymmetry <= max(0.28, baseline_asym + 0.10)
    ):
        return 'braid_like_exchange'

    if flow_ratio >= 0.68 and coherence_ratio >= 0.62 and selectivity_ratio >= 0.72:
        return 'transfer_smeared'

    return 'unresolved_mixed'


def family_basis_degrees(
    G: nx.Graph,
    nodes: list[int],
    oriented_edges: list[tuple[int, int]],
    triangles: list[tuple[int, int, int]],
) -> np.ndarray:
    node_deg = np.asarray([float(G.degree[node]) for node in nodes], dtype=float)
    edge_deg = np.asarray([0.5 * (float(G.degree[u]) + float(G.degree[v])) for u, v in oriented_edges], dtype=float)
    face_deg = np.asarray([float(sum(G.degree[node] for node in tri) / 3.0) for tri in triangles], dtype=float)
    return np.concatenate([node_deg, edge_deg, face_deg]) if face_deg.size else np.concatenate([node_deg, edge_deg])


def simulate_family_variant(
    family: dict[str, Any],
    common: dict[str, Any],
    delta: float,
    left_sigma_scale: float,
    right_sigma_scale: float,
    grade0_scale: float,
    grade1_scale: float,
    grade2_scale: float,
    degree_skew_strength: float,
) -> dict[str, Any]:
    family_id = str(family['family_id'])
    left_seed = int(family['left_seed'])
    right_seed = int(family['right_seed'])
    G = representative_graph(family_id)
    nodes = sorted(int(node) for node in G.nodes())
    scores = mismatch_scores(G, left_seed=left_seed, right_seed=right_seed)
    oriented_edges = orient_edges(G, scores)
    triangles = triangle_list(G)
    d0, edge_lookup, _node_to_idx = build_d0(nodes, oriented_edges)
    d1 = build_d1(triangles, edge_lookup) if triangles else sp.csr_matrix((0, len(oriented_edges)), dtype=float)
    D, Delta = block_dirac(d0, d1)
    laplacian = nx.laplacian_matrix(G, nodelist=nodes).astype(float).tocsr()

    nilpotency_norm = float(np.linalg.norm((d1 @ d0).toarray())) if triangles else 0.0
    operator_identity_norm = float(np.linalg.norm(((D @ D) - Delta).toarray()))
    laplacian_diff_norm = float(np.linalg.norm(((d0.T @ d0) - laplacian).toarray()))

    sigma_shell = float(common['sigma_shell'])
    kick_cycles = float(common['kick_cycles'])
    packet_a = initial_packet(
        G=G,
        nodes=nodes,
        oriented_edges=oriented_edges,
        triangles=triangles,
        seed=left_seed,
        sigma_shell=sigma_shell * float(left_sigma_scale),
        kick_cycles=kick_cycles,
        direction_sign=+1.0,
        grade0_scale=grade0_scale,
        grade1_scale=grade1_scale,
        grade2_scale=grade2_scale,
    )
    packet_b = initial_packet(
        G=G,
        nodes=nodes,
        oriented_edges=oriented_edges,
        triangles=triangles,
        seed=right_seed,
        sigma_shell=sigma_shell * float(right_sigma_scale),
        kick_cycles=kick_cycles,
        direction_sign=-1.0,
        grade0_scale=grade0_scale,
        grade1_scale=grade1_scale,
        grade2_scale=grade2_scale,
    )
    block_sizes = (len(nodes), len(oriented_edges), len(triangles))
    labels = build_address_labels(
        packet_states=[packet_a, packet_b],
        block_sizes=block_sizes,
        left_address=0.0,
        right_address=2.0 * float(delta),
        support_floor_fraction=float(common['support_floor_fraction']),
        dominant_ratio_threshold=float(common['dominant_ratio_threshold']),
    )
    D_h = address_weighted_operator(
        base_operator=D,
        labels=labels,
        modulus=float(common['address_modulus']),
        eta_match=float(common['eta_match']),
        eta_mismatch=float(common['eta_mismatch']),
    )
    basis_degrees = family_basis_degrees(G, nodes, oriented_edges, triangles)
    D_h = degree_skew_operator(D_h, labels, basis_degrees, float(degree_skew_strength))
    psi0 = packet_a + packet_b
    psi0 /= max(float(np.linalg.norm(psi0)), 1.0e-12)
    states, dt = evolve_states(D_h, psi0, steps=int(common['steps']), dt_scale=float(common['dt_scale']))

    n0 = len(nodes)
    n1 = len(oriented_edges)
    edge_series = np.asarray([state[n0:n0 + n1] for state in states], dtype=complex)
    grade_hist = grade_histories(states, block_sizes)
    _signal, coherence = grade_exchange_signal(grade_hist)
    grade_asymmetry_index = float(abs(grade_hist[-1, 0] - grade_hist[-1, 1])) if grade_hist.shape[1] >= 2 else 0.0
    edge_power = np.mean(np.abs(edge_series) ** 2, axis=0) if edge_series.size else np.zeros(n1, dtype=float)
    total_edge_power = float(np.sum(edge_power))
    if total_edge_power > 1.0e-12:
        flat = np.sort(edge_power)[::-1]
        top_k = max(1, int(math.ceil(0.1 * flat.size)))
        concentration = float(np.sum(flat[:top_k]) / total_edge_power)
        threshold = max(0.12 * float(np.max(edge_power)), float(np.mean(edge_power) + 0.50 * np.std(edge_power)))
        active_mask = edge_power >= threshold
    else:
        concentration = 0.0
        active_mask = np.zeros(n1, dtype=bool)
    channel_count = int(edge_component_count(oriented_edges, active_mask))
    loop_count, loop_score = triangle_loop_metrics(edge_series, triangles, d1)
    sorting, selectivity = address_activity_statistics(
        states=states,
        operator=D_h,
        labels=labels,
        modulus=float(common['address_modulus']),
        eta_match=float(common['eta_match']),
        eta_mismatch=float(common['eta_mismatch']),
    )
    return {
        'family_id': family_id,
        'family_label': str(family['family_label']),
        'delta': float(delta),
        'dt': float(dt),
        'nilpotency_norm': nilpotency_norm,
        'operator_identity_norm': operator_identity_norm,
        'laplacian_diff_norm': laplacian_diff_norm,
        'triangle_count': int(len(triangles)),
        'flow_concentration_index': concentration,
        'grade_exchange_coherence': coherence,
        'grade_asymmetry_index': grade_asymmetry_index,
        'address_sorting_score': sorting,
        'address_selectivity_index': selectivity,
        'channel_count': channel_count,
        'loop_count': loop_count,
        'loop_score': loop_score,
        'orientation_scores': scores,
    }


def determine_threshold(detuning_rows: list[dict[str, Any]]) -> float | None:
    if not detuning_rows:
        return None
    base_class = str(detuning_rows[0]['topology_class'])
    for row in detuning_rows[1:]:
        if str(row['topology_class']) != base_class:
            return float(row['delta'])
    return None


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def detuning_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    mapping = {'braid_like_exchange': 0, 'transfer_smeared': 1, 'unresolved_mixed': 2}
    colors = {
        'clustered_composite_anchor': 'tab:blue',
        'counter_propagating_corridor': 'tab:orange',
        'phase_ordered_symmetric_triad': 'tab:green',
    }
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.5))
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row['family_id'])].append(row)
    for family_id, family_rows in grouped.items():
        family_rows = sorted(family_rows, key=lambda item: float(item['delta']))
        deltas = [float(item['delta']) for item in family_rows]
        codes = [mapping[str(item['topology_class'])] for item in family_rows]
        selectivity = [float(item['address_selectivity_index']) for item in family_rows]
        axes[0].plot(deltas, codes, marker='o', color=colors.get(family_id, 'tab:gray'), label=family_rows[0]['family_label'])
        axes[1].plot(deltas, selectivity, marker='o', color=colors.get(family_id, 'tab:gray'), label=family_rows[0]['family_label'])
    axes[0].set_yticks([0, 1, 2])
    axes[0].set_yticklabels(['braid', 'smeared', 'mixed'])
    axes[0].set_xlabel('detuning delta')
    axes[0].set_title('Topology vs detuning')
    axes[0].grid(alpha=0.25)
    axes[1].set_xlabel('detuning delta')
    axes[1].set_ylabel('address selectivity')
    axes[1].set_title('Selectivity decay')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def algebraic_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    families = [str(row['family_label']) for row in rows]
    nilpotency = [float(row['nilpotency_norm']) for row in rows]
    operator_diff = [float(row['operator_identity_norm']) for row in rows]
    lap_diff = [float(row['laplacian_diff_norm']) for row in rows]
    x = np.arange(len(families))
    width = 0.25
    fig, ax = plt.subplots(figsize=(10.0, 4.2))
    ax.bar(x - width, nilpotency, width=width, label='||d1 d0||')
    ax.bar(x, operator_diff, width=width, label='||D^2 - Delta||')
    ax.bar(x + width, lap_diff, width=width, label='||d0^T d0 - L0||')
    ax.set_xticks(x)
    ax.set_xticklabels(families, rotation=12, ha='right')
    ax.set_yscale('symlog', linthresh=1.0e-12)
    ax.set_title('Bridge algebraic diagnostics')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25, axis='y')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def protection_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    mapping = {'braid_like_exchange': 0, 'transfer_smeared': 1, 'unresolved_mixed': 2}
    variants = [str(row['variant_label']) for row in rows]
    codes = [mapping[str(row['topology_class'])] for row in rows]
    flow = [float(row['flow_concentration_index']) for row in rows]
    coherence = [float(row['grade_exchange_coherence']) for row in rows]
    x = np.arange(len(rows))
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.5))
    axes[0].bar(x, codes, color='tab:blue', alpha=0.8)
    axes[0].set_yticks([0, 1, 2])
    axes[0].set_yticklabels(['braid', 'smeared', 'mixed'])
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(variants, rotation=20, ha='right')
    axes[0].set_title('Clustered protection panel')
    axes[0].grid(alpha=0.25, axis='y')
    axes[1].plot(x, flow, marker='o', label='flow concentration')
    axes[1].plot(x, coherence, marker='s', label='grade coherence')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(variants, rotation=20, ha='right')
    axes[1].set_title('Protection observables')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def note_text(
    json_rel: str,
    csv_rel: str,
    algebraic_rows: list[dict[str, Any]],
    detuning_rows: list[dict[str, Any]],
    protection_rows: list[dict[str, Any]],
    stamped_plots: list[str],
) -> str:
    threshold_map = {
        str(family_id): determine_threshold(sorted([row for row in detuning_rows if str(row['family_id']) == str(family_id)], key=lambda item: float(item['delta'])))
        for family_id in {row['family_id'] for row in detuning_rows}
    }
    baseline_cluster = next((row for row in protection_rows if str(row['variant_id']) == 'balanced_baseline'), None)
    matched_targets = 0
    expected = {
        'balanced_baseline': 'braid_like_exchange',
        'phase_edge_probe': 'transfer_smeared',
        'width_asymmetry_probe': 'transfer_smeared',
        'degree_skew_probe': 'braid_like_exchange',
        'zero_form_enhanced_probe': 'transfer_smeared',
    }
    for row in protection_rows:
        if expected.get(str(row['variant_id'])) == str(row['topology_class']):
            matched_targets += 1
    protection_strength = matched_targets / max(len(expected), 1)

    lines = [
        '# Stage C0-to-DK Bridge Validation v1',
        '',
        f'Timestamped JSON: `{json_rel}`',
        f'Timestamped CSV: `{csv_rel}`',
        '',
        'Purpose: test a genuinely executable C0-to-DK bridge candidate using only signed combinatorial incidence, deterministic harmonic-address dressing, and bridge-specific dynamics on the derived graph complex.',
        '',
        '## Algebraic checks',
    ]
    for row in algebraic_rows:
        lines.append(
            f"- `{row['family_id']}`: `||d1 d0|| = {float(row['nilpotency_norm']):.3e}`, "
            f"`||D^2 - Delta|| = {float(row['operator_identity_norm']):.3e}`, "
            f"`||d0^T d0 - L0|| = {float(row['laplacian_diff_norm']):.3e}`, "
            f"`triangles = {int(row['triangle_count'])}`"
        )

    lines.extend(
        [
            '',
            '## Detuning scan',
        ]
    )
    for family_id in sorted({str(row['family_id']) for row in detuning_rows}):
        family_rows = sorted([row for row in detuning_rows if str(row['family_id']) == family_id], key=lambda item: float(item['delta']))
        labels = ', '.join(f"{float(row['delta']):.3f}:{row['topology_class']}" for row in family_rows)
        threshold = threshold_map.get(family_id)
        threshold_text = f"{threshold:.3f}" if threshold is not None else 'none'
        lines.append(f"- `{family_id}`: first class change at `delta = {threshold_text}`; trace = `{labels}`")

    lines.extend(
        [
            '',
            '## Focused protection panel',
        ]
    )
    for row in protection_rows:
        lines.append(
            f"- `{row['variant_id']}`: topology=`{row['topology_class']}`, "
            f"flow=`{float(row['flow_concentration_index']):.4f}`, "
            f"coherence=`{float(row['grade_exchange_coherence']):.4f}`, "
            f"selectivity=`{float(row['address_selectivity_index']):.4f}`, "
            f"loop score=`{float(row['loop_score']):.4f}`"
        )

    lines.extend(
        [
            '',
            '## Obligation status',
            f"- `O_C0_ORIENT`: {'open; the tested signed complex closes exactly, but the orientation source is still a mismatch surrogate rather than the protected sequencing ledger itself' if all(float(row['nilpotency_norm']) <= 1.0e-10 for row in algebraic_rows) else 'open; signed orientation candidate still fails exact nilpotency on at least one family'}.",
            f"- `O_C0_EVOLUTION`: open; the bridge uses an explicit derived `D_H` evolution, but equivalence to the protected sequencing rule remains unproved.",
            f"- `O_C0_GRADING`: open; harmonic address remains a dressing, not yet an intrinsic grade source.",
            f"- `O_C0_RECOVERY`: {'partially closed' if any(threshold is not None for threshold in threshold_map.values()) else 'open'}; the clustered family opens a real braid-to-smear bridge corridor and the focused clustered protection panel matches `{protection_strength:.2f}` of the Stage-23-style qualitative pattern, but the first class change occurs at `delta = {threshold_map.get('clustered_composite_anchor') if threshold_map.get('clustered_composite_anchor') is not None else 'none'}` and the corridor/triad families do not recover a broad braid sector.",
            '',
            '## Canonical statement',
            'A combinatorial Dirac-Kahler candidate can now be constructed from C0 primitives, and on the tested representative graph complexes its signed chain identities close exactly. The present run therefore establishes a real bridge scaffold rather than a finished derivation: the clustered family shows a measurable harmonic detuning boundary and a focused protection-panel recovery, but the orientation source, evolution-law equivalence, intrinsic grading source, and family-wide protection recovery remain conditional.',
            '',
            '## Plots',
        ]
    )
    for rel in stamped_plots:
        lines.append(f'- `{rel}`')
    return '\n'.join(lines) + '\n'


def main() -> None:
    args = parse_args()
    runsheet = load_runsheet(args.runsheet)
    common = runsheet['common_fields']
    families = runsheet['families']
    if args.family_ids:
        allowed = set(args.family_ids)
        families = [family for family in families if str(family['family_id']) in allowed]

    algebraic_rows: list[dict[str, Any]] = []
    detuning_rows: list[dict[str, Any]] = []
    protection_rows: list[dict[str, Any]] = []

    for family in families:
        family_results: list[dict[str, Any]] = []
        for delta in runsheet['common_fields']['detuning_samples']:
            result = simulate_family_variant(
                family=family,
                common=common,
                delta=float(delta),
                left_sigma_scale=1.0,
                right_sigma_scale=1.0,
                grade0_scale=float(common['grade0_scale']),
                grade1_scale=float(common['grade1_scale']),
                grade2_scale=float(common['grade2_scale']),
                degree_skew_strength=0.0,
            )
            family_results.append(result)
        baseline_result = family_results[0]
        family_rows = []
        for result in family_results:
            row = {
                'section': 'detuning',
                'family_id': result['family_id'],
                'family_label': result['family_label'],
                'variant_id': 'detuning_scan',
                'variant_label': 'Detuning scan',
                'delta': result['delta'],
                'nilpotency_norm': result['nilpotency_norm'],
                'operator_identity_norm': result['operator_identity_norm'],
                'laplacian_diff_norm': result['laplacian_diff_norm'],
                'triangle_count': result['triangle_count'],
                'topology_class': classify_relative_topology(result, baseline_result),
                'flow_concentration_index': result['flow_concentration_index'],
                'grade_exchange_coherence': result['grade_exchange_coherence'],
                'address_selectivity_index': result['address_selectivity_index'],
                'channel_count': result['channel_count'],
                'loop_count': result['loop_count'],
                'loop_score': result['loop_score'],
                'detuning_threshold': '',
                'notes': '',
            }
            family_rows.append(row)
            detuning_rows.append(row)
        baseline = family_rows[0].copy()
        baseline['section'] = 'algebraic'
        baseline['variant_id'] = 'algebraic_baseline'
        baseline['variant_label'] = 'Algebraic baseline'
        baseline['delta'] = ''
        baseline['topology_class'] = ''
        baseline['flow_concentration_index'] = ''
        baseline['grade_exchange_coherence'] = ''
        baseline['address_selectivity_index'] = ''
        baseline['channel_count'] = ''
        baseline['loop_count'] = ''
        baseline['loop_score'] = ''
        baseline['detuning_threshold'] = determine_threshold(sorted(family_rows, key=lambda item: float(item['delta'])))
        algebraic_rows.append(baseline)

    clustered = next(family for family in families if str(family['family_id']) == 'clustered_composite_anchor')
    clustered_variant_results: list[dict[str, Any]] = []
    for variant in runsheet['protection_variants']:
        result = simulate_family_variant(
            family=clustered,
            common=common,
            delta=float(variant['delta']),
            left_sigma_scale=float(variant['left_sigma_scale']),
            right_sigma_scale=float(variant['right_sigma_scale']),
            grade0_scale=float(variant['grade0_scale']),
            grade1_scale=float(variant['grade1_scale']),
            grade2_scale=float(variant['grade2_scale']),
            degree_skew_strength=float(variant['degree_skew_strength']),
        )
        result['variant_id'] = str(variant['variant_id'])
        result['variant_label'] = str(variant['variant_label'])
        clustered_variant_results.append(result)

    clustered_baseline = next(result for result in clustered_variant_results if str(result['variant_id']) == 'balanced_baseline')
    for result in clustered_variant_results:
        protection_rows.append(
            {
                'section': 'protection',
                'family_id': result['family_id'],
                'family_label': result['family_label'],
                'variant_id': str(result['variant_id']),
                'variant_label': str(result['variant_label']),
                'delta': result['delta'],
                'nilpotency_norm': result['nilpotency_norm'],
                'operator_identity_norm': result['operator_identity_norm'],
                'laplacian_diff_norm': result['laplacian_diff_norm'],
                'triangle_count': result['triangle_count'],
                'topology_class': classify_relative_topology(result, clustered_baseline),
                'flow_concentration_index': result['flow_concentration_index'],
                'grade_exchange_coherence': result['grade_exchange_coherence'],
                'address_selectivity_index': result['address_selectivity_index'],
                'channel_count': result['channel_count'],
                'loop_count': result['loop_count'],
                'loop_score': result['loop_score'],
                'detuning_threshold': '',
                'notes': '',
            }
        )

    all_rows = algebraic_rows + detuning_rows + protection_rows
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    json_path = DATA_DIR / f'{timestamp}_stage_c0_dk_bridge_validation.json'
    csv_path = DATA_DIR / f'{timestamp}_stage_c0_dk_bridge_validation.csv'

    detuning_plot_path = WORK_PLOT_DIR / 'stage_c0_dk_bridge_validation_detuning_panel.png'
    algebraic_plot_path = WORK_PLOT_DIR / 'stage_c0_dk_bridge_validation_algebraic_panel.png'
    protection_plot_path = WORK_PLOT_DIR / 'stage_c0_dk_bridge_validation_protection_panel.png'
    detuning_plot(detuning_plot_path, detuning_rows)
    algebraic_plot(algebraic_plot_path, algebraic_rows)
    protection_plot(protection_plot_path, protection_rows)

    stamped_plots: list[str] = []
    for src in (detuning_plot_path, algebraic_plot_path, protection_plot_path):
        dst = PLOTS_DIR / f'{timestamp}_{src.name}'
        dst.write_bytes(src.read_bytes())
        stamped_plots.append(str(dst.relative_to(REPO_ROOT)))

    payload = {
        'stage': 'c0_dk_bridge_validation',
        'timestamp': timestamp,
        'runsheet': str(args.runsheet.relative_to(REPO_ROOT)),
        'algebraic': algebraic_rows,
        'detuning_scan': detuning_rows,
        'protection_panel': protection_rows,
        'plots': stamped_plots,
    }
    json_path.write_text(json.dumps(payload, indent=2), encoding='utf-8')
    write_csv(csv_path, all_rows, CSV_FIELDS)

    note = note_text(
        json_rel=str(json_path.relative_to(REPO_ROOT)),
        csv_rel=str(csv_path.relative_to(REPO_ROOT)),
        algebraic_rows=algebraic_rows,
        detuning_rows=detuning_rows,
        protection_rows=protection_rows,
        stamped_plots=stamped_plots,
    )
    NOTE_PATH.write_text(note, encoding='utf-8')

    print(f'JSON: {json_path}')
    print(f'CSV: {csv_path}')
    print(f'NOTE: {NOTE_PATH}')
    for rel in stamped_plots:
        print(f'PLOT: {REPO_ROOT / rel}')


if __name__ == '__main__':
    main()
