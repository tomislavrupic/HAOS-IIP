#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import scipy.sparse as sp

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt, save_atlas_payload
from stage_c0_dk_bridge_validation import (
    address_activity_statistics,
    address_weighted_operator,
    block_dirac,
    build_d0,
    build_d1,
    edge_component_count,
    grade_exchange_signal,
    grade_histories,
    triangle_loop_metrics,
)


RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_11_integrated_austerity_resolution_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_11_Integrated_Austerity_Resolution_Probe_v1.md'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage_11_integrated_austerity'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'run_id',
    'family_id',
    'family_label',
    'detuning',
    'protection_strength',
    'ledger_alignment',
    'delta_equiv_coarse_mean',
    'delta_equiv_fine_mean',
    'delta_equiv_refine_ratio',
    'sector_entropy_coarse',
    'sector_entropy_fine',
    'base_survival_time',
    'degree_skew_survival_time',
    'motif_survival_time',
    'radius_minus_survival_time',
    'radius_center_survival_time',
    'radius_plus_survival_time',
    'ordering_rank',
    'resolution_mode',
    'final_topology_class',
    'trace_plot',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run Stage 11 integrated austerity resolution probe.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def matched_family_graph(family_id: str) -> nx.Graph:
    G = nx.cycle_graph(12)
    if family_id == 'clustered_braid_seed':
        matching = [(0, 2), (1, 11), (3, 5), (4, 6), (7, 9), (8, 10)]
    elif family_id == 'corridor_transport_seed':
        matching = [(0, 6), (1, 7), (2, 8), (3, 9), (4, 10), (5, 11)]
    elif family_id == 'triad_competition_seed':
        matching = [(0, 2), (1, 4), (3, 5), (6, 8), (7, 10), (9, 11)]
    else:
        raise ValueError(f'unknown family_id: {family_id}')
    G.add_edges_from(matching)
    return G


def triangle_list(G: nx.Graph) -> list[tuple[int, int, int]]:
    triangles: list[tuple[int, int, int]] = []
    nodes = sorted(G.nodes())
    for a in nodes:
        nbrs = sorted(node for node in G.neighbors(a) if node > a)
        for i, b in enumerate(nbrs):
            for c in nbrs[i + 1:]:
                if G.has_edge(b, c):
                    tri = tuple(sorted((int(a), int(b), int(c))))
                    if tri not in triangles:
                        triangles.append(tri)
    return triangles


def node_triangle_counts(G: nx.Graph, triangles: list[tuple[int, int, int]]) -> dict[int, int]:
    counts = {int(node): 0 for node in G.nodes()}
    for tri in triangles:
        for node in tri:
            counts[int(node)] += 1
    return counts


def motif_nodes_for_family(G: nx.Graph, triangles: list[tuple[int, int, int]], left_seed: int, right_seed: int) -> set[int]:
    motif_nodes: set[int] = set()
    for tri in triangles:
        tri_set = set(int(node) for node in tri)
        if left_seed in tri_set or right_seed in tri_set:
            motif_nodes.update(tri_set)
    if not motif_nodes:
        motif_nodes.update(nx.shortest_path(G, source=left_seed, target=right_seed))
    return motif_nodes


def shell_distance_map(G: nx.Graph, seed: int) -> dict[int, int]:
    return {int(node): int(distance) for node, distance in nx.single_source_shortest_path_length(G, seed).items()}


def intrinsic_node_addresses(
    G: nx.Graph,
    left_seed: int,
    right_seed: int,
    modulus: float,
) -> tuple[np.ndarray, np.ndarray]:
    triangles = triangle_list(G)
    tri_counts = node_triangle_counts(G, triangles)
    left_shell = shell_distance_map(G, left_seed)
    right_shell = shell_distance_map(G, right_seed)
    nodes = sorted(int(node) for node in G.nodes())
    addresses = np.zeros(len(nodes), dtype=float)
    sectors = np.zeros(len(nodes), dtype=float)
    for idx, node in enumerate(nodes):
        ls = int(left_shell.get(node, 0))
        rs = int(right_shell.get(node, 0))
        deg = int(G.degree[node])
        tri = int(tri_counts.get(node, 0))
        addresses[idx] = float((ls + 2 * rs + deg + tri) % int(modulus))
        sectors[idx] = 1.0 if rs < ls else 0.0
    return addresses, sectors


def orient_edges_from_mismatch(
    G: nx.Graph,
    mismatch: np.ndarray,
    node_order: list[int],
    tri_counts: dict[int, int],
) -> list[tuple[int, int]]:
    mismatch_map = {int(node): float(mismatch[idx]) for idx, node in enumerate(node_order)}
    oriented: list[tuple[int, int]] = []
    for u, v in sorted((min(a, b), max(a, b)) for a, b in G.edges()):
        mu = mismatch_map[int(u)]
        mv = mismatch_map[int(v)]
        if mu > mv + 1.0e-12:
            oriented.append((int(u), int(v)))
        elif mv > mu + 1.0e-12:
            oriented.append((int(v), int(u)))
        else:
            pu = tri_counts.get(int(u), 0) % 2
            pv = tri_counts.get(int(v), 0) % 2
            if pu < pv:
                oriented.append((int(u), int(v)))
            elif pv < pu:
                oriented.append((int(v), int(u)))
            else:
                oriented.append((int(u), int(v)))
    return oriented


def basis_supports(
    nodes: list[int],
    oriented_edges: list[tuple[int, int]],
    triangles: list[tuple[int, int, int]],
) -> list[tuple[int, ...]]:
    supports: list[tuple[int, ...]] = []
    supports.extend([(int(node),) for node in nodes])
    supports.extend([tuple(sorted((int(u), int(v)))) for u, v in oriented_edges])
    supports.extend([tuple(sorted(int(node) for node in tri)) for tri in triangles])
    return supports


def basis_labels_from_intrinsic(
    nodes: list[int],
    oriented_edges: list[tuple[int, int]],
    triangles: list[tuple[int, int, int]],
    node_addresses: np.ndarray,
    node_sectors: np.ndarray,
    detuning: float,
    modulus: float,
) -> np.ndarray:
    node_index = {int(node): idx for idx, node in enumerate(nodes)}
    labels: list[float] = []
    for node in nodes:
        idx = node_index[int(node)]
        labels.append(float((node_addresses[idx] + 2.0 * detuning * node_sectors[idx]) % modulus))
    for u, v in oriented_edges:
        ui = node_index[int(u)]
        vi = node_index[int(v)]
        sector = 0.5 * (node_sectors[ui] + node_sectors[vi])
        base = 0.5 * (node_addresses[ui] + node_addresses[vi])
        labels.append(float((base + 2.0 * detuning * sector) % modulus))
    for tri in triangles:
        idxs = [node_index[int(node)] for node in tri]
        sector = float(np.mean(node_sectors[idxs]))
        base = float(np.mean(node_addresses[idxs]))
        labels.append(float((base + 2.0 * detuning * sector) % modulus))
    return np.asarray(labels, dtype=float)


def initial_packet(
    G: nx.Graph,
    nodes: list[int],
    oriented_edges: list[tuple[int, int]],
    triangles: list[tuple[int, int, int]],
    seed: int,
    kernel_width: float,
    direction_sign: float,
    grade0_scale: float,
    grade1_scale: float,
    grade2_scale: float,
) -> np.ndarray:
    shell = shell_distance_map(G, seed)
    node_vals = np.zeros(len(nodes), dtype=complex)
    for idx, node in enumerate(nodes):
        dist = float(shell.get(int(node), 999))
        envelope = math.exp(-(dist * dist) / (2.0 * kernel_width * kernel_width))
        phase = np.exp(1j * direction_sign * 0.35 * dist)
        node_vals[idx] = envelope * phase

    edge_vals = np.zeros(len(oriented_edges), dtype=complex)
    for idx, (u, v) in enumerate(oriented_edges):
        du = float(shell.get(int(u), 999))
        dv = float(shell.get(int(v), 999))
        env = 0.5 * (
            math.exp(-(du * du) / (2.0 * kernel_width * kernel_width))
            + math.exp(-(dv * dv) / (2.0 * kernel_width * kernel_width))
        )
        edge_vals[idx] = env * np.exp(1j * direction_sign * 0.35 * (dv - du))

    face_vals = np.zeros(len(triangles), dtype=complex)
    for idx, tri in enumerate(triangles):
        dists = [float(shell.get(int(node), 999)) for node in tri]
        env = float(np.mean([math.exp(-(dist * dist) / (2.0 * kernel_width * kernel_width)) for dist in dists]))
        face_vals[idx] = env * np.exp(1j * direction_sign * 0.35 * float(np.mean(dists)))

    psi = np.concatenate(
        [
            float(grade0_scale) * node_vals,
            float(grade1_scale) * edge_vals,
            float(grade2_scale) * face_vals,
        ]
    )
    norm = max(float(np.linalg.norm(psi)), 1.0e-12)
    return psi / norm


def node_kernel_matrix(G: nx.Graph, radius: int) -> np.ndarray:
    nodes = sorted(int(node) for node in G.nodes())
    index = {node: idx for idx, node in enumerate(nodes)}
    mat = np.zeros((len(nodes), len(nodes)), dtype=float)
    for source in nodes:
        lengths = nx.single_source_shortest_path_length(G, source=source, cutoff=radius)
        for target, distance in lengths.items():
            mat[index[source], index[int(target)]] = 1.0 / (float(distance) + 1.0)
        row_sum = float(np.sum(mat[index[source]]))
        if row_sum > 0.0:
            mat[index[source]] /= row_sum
    return mat


def support_distance_vector(
    G: nx.Graph,
    supports: list[tuple[int, ...]],
    row_ids: np.ndarray,
    col_ids: np.ndarray,
) -> np.ndarray:
    all_lengths = dict(nx.all_pairs_shortest_path_length(G))
    distances = np.zeros(len(row_ids), dtype=float)
    for idx, (r, c) in enumerate(zip(row_ids, col_ids)):
        sup_r = supports[int(r)]
        sup_c = supports[int(c)]
        best = min(float(all_lengths[int(u)][int(v)]) for u in sup_r for v in sup_c)
        distances[idx] = best
    return distances


def radius_operator(base_operator: sp.csr_matrix, support_distances: np.ndarray, radius: int) -> sp.csr_matrix:
    coo = base_operator.tocoo()
    mask = support_distances <= float(radius)
    weights = np.where(mask, 1.0 / (support_distances + 1.0), 0.0)
    return sp.csr_matrix((coo.data * weights, (coo.row, coo.col)), shape=base_operator.shape).tocsr()


def apply_degree_skew(operator: sp.csr_matrix, basis_degrees: np.ndarray, strength: float) -> sp.csr_matrix:
    if strength <= 0.0:
        return operator
    deg = np.asarray(basis_degrees, dtype=float)
    deg = deg / max(float(np.max(deg)), 1.0e-12)
    coo = operator.tocoo()
    mult = 1.0 + float(strength) * 0.5 * (deg[coo.row] + deg[coo.col])
    return sp.csr_matrix((coo.data * mult, (coo.row, coo.col)), shape=operator.shape).tocsr()


def apply_motif_boost(operator: sp.csr_matrix, supports: list[tuple[int, ...]], motif_nodes: set[int], strength: float) -> sp.csr_matrix:
    if strength <= 0.0:
        return operator
    motif_flags = np.asarray([1.0 if motif_nodes.intersection(support) else 0.0 for support in supports], dtype=float)
    coo = operator.tocoo()
    mult = 1.0 + float(strength) * 0.5 * (motif_flags[coo.row] + motif_flags[coo.col])
    return sp.csr_matrix((coo.data * mult, (coo.row, coo.col)), shape=operator.shape).tocsr()


def basis_degrees(nodes: list[int], oriented_edges: list[tuple[int, int]], triangles: list[tuple[int, int, int]], G: nx.Graph) -> np.ndarray:
    node_deg = np.asarray([float(G.degree[node]) for node in nodes], dtype=float)
    edge_deg = np.asarray([0.5 * (float(G.degree[u]) + float(G.degree[v])) for u, v in oriented_edges], dtype=float)
    face_deg = np.asarray([float(np.mean([G.degree[node] for node in tri])) for tri in triangles], dtype=float) if triangles else np.zeros(0, dtype=float)
    return np.concatenate([node_deg, edge_deg, face_deg]) if face_deg.size else np.concatenate([node_deg, edge_deg])


def field_update(operator: sp.csr_matrix, psi: np.ndarray, epsilon: float) -> np.ndarray:
    next_state = psi + 1j * float(epsilon) * (operator @ psi)
    next_state = np.asarray(next_state, dtype=complex)
    norm = max(float(np.linalg.norm(next_state)), 1.0e-12)
    return next_state / norm


def node_activity_from_state(
    psi: np.ndarray,
    nodes: list[int],
    oriented_edges: list[tuple[int, int]],
    triangles: list[tuple[int, int, int]],
) -> np.ndarray:
    n0 = len(nodes)
    n1 = len(oriented_edges)
    node_power = np.abs(psi[:n0]) ** 2
    edge_power = np.abs(psi[n0:n0 + n1]) ** 2
    face_power = np.abs(psi[n0 + n1:]) ** 2
    activity = node_power.copy()
    node_index = {int(node): idx for idx, node in enumerate(nodes)}
    for edge_idx, (u, v) in enumerate(oriented_edges):
        activity[node_index[int(u)]] += 0.5 * edge_power[edge_idx]
        activity[node_index[int(v)]] += 0.5 * edge_power[edge_idx]
    for face_idx, tri in enumerate(triangles):
        for node in tri:
            activity[node_index[int(node)]] += face_power[face_idx] / 3.0
    return activity


def classify_topology(
    flow_concentration: float,
    coherence: float,
    selectivity: float,
    loop_score: float,
    channel_count: int,
    capture_localization: float,
) -> str:
    if loop_score >= 0.18 and flow_concentration >= 0.24 and coherence >= 0.42 and selectivity >= 0.90 and channel_count <= 2:
        return 'braid_like_exchange'
    if capture_localization >= 0.34 and channel_count <= 1 and selectivity >= 0.60:
        return 'localized_capture'
    if flow_concentration >= 0.20 and coherence >= 0.28 and selectivity >= 0.72:
        return 'transfer_smeared'
    return 'dispersive_pass'


def metrics_for_state(
    states_window: list[np.ndarray],
    operator: sp.csr_matrix,
    labels: np.ndarray,
    nodes: list[int],
    oriented_edges: list[tuple[int, int]],
    triangles: list[tuple[int, int, int]],
    d1: sp.csr_matrix,
    block_sizes: tuple[int, int, int],
    modulus: float,
    eta_match: float,
    eta_mismatch: float,
) -> dict[str, float | int | str]:
    n0, n1, _n2 = block_sizes
    edge_series = np.asarray([state[n0:n0 + n1] for state in states_window], dtype=complex)
    grade_hist = grade_histories(states_window, block_sizes)
    _signal, coherence = grade_exchange_signal(grade_hist)
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
        states=states_window,
        operator=operator,
        labels=labels,
        modulus=modulus,
        eta_match=eta_match,
        eta_mismatch=eta_mismatch,
    )
    node_power = np.abs(states_window[-1][:n0]) ** 2
    capture_localization = float(np.max(node_power) / max(float(np.sum(node_power)), 1.0e-12))
    topology = classify_topology(
        flow_concentration=concentration,
        coherence=float(coherence),
        selectivity=float(selectivity),
        loop_score=float(loop_score),
        channel_count=channel_count,
        capture_localization=capture_localization,
    )
    return {
        'topology_class': topology,
        'flow_concentration_index': concentration,
        'grade_exchange_coherence': float(coherence),
        'address_selectivity_index': float(selectivity),
        'address_sorting_score': float(sorting),
        'loop_count': int(loop_count),
        'loop_score': float(loop_score),
        'channel_count': int(channel_count),
        'capture_localization_index': capture_localization,
    }


def conditional_entropy(labels: list[str], bands: list[str]) -> float:
    if not labels or not bands or len(labels) != len(bands):
        return 0.0
    total = len(labels)
    by_band: dict[str, list[str]] = defaultdict(list)
    for label, band in zip(labels, bands):
        by_band[str(band)].append(str(label))
    entropy = 0.0
    for band, values in by_band.items():
        weight = len(values) / total
        counts = Counter(values)
        probs = np.asarray([count / len(values) for count in counts.values()], dtype=float)
        band_entropy = -float(np.sum(probs * np.log(np.maximum(probs, 1.0e-12)))) / max(math.log(4.0), 1.0)
        entropy += weight * band_entropy
    return float(entropy)


def correlation_or_zero(a: list[float], b: list[float]) -> float:
    if len(a) < 3 or len(b) < 3:
        return 0.0
    arr_a = np.asarray(a, dtype=float)
    arr_b = np.asarray(b, dtype=float)
    if float(np.std(arr_a)) <= 1.0e-12 or float(np.std(arr_b)) <= 1.0e-12:
        return 0.0
    return float(np.corrcoef(arr_a, arr_b)[0, 1])


def trace_plot(path: Path, run_id: str, topology_codes: list[int], mismatch_gradient: list[float], delta_equiv: list[float]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.2))
    axes[0].plot(range(len(topology_codes)), topology_codes, marker='o')
    axes[0].set_yticks([0, 1, 2, 3])
    axes[0].set_yticklabels(['braid', 'smear', 'capture', 'disp'])
    axes[0].set_title('Topology trace')
    axes[0].grid(alpha=0.25)
    axes[1].plot(range(len(mismatch_gradient)), mismatch_gradient, marker='o', color='tab:orange')
    axes[1].set_title('Mismatch gradient')
    axes[1].grid(alpha=0.25)
    axes[2].plot(range(len(delta_equiv)), delta_equiv, marker='o', color='tab:red')
    axes[2].set_title('Evolution mismatch')
    axes[2].grid(alpha=0.25)
    fig.suptitle(run_id)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def summary_plot(path: Path, rows: list[dict[str, Any]], field: str, title: str) -> None:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row['family_label'])].append(row)
    fig, ax = plt.subplots(figsize=(10.5, 4.5))
    for family, family_rows in grouped.items():
        family_rows = sorted(family_rows, key=lambda item: (float(item['detuning']), float(item['protection_strength'])))
        x = np.arange(len(family_rows))
        y = [float(item[field]) for item in family_rows]
        ax.plot(x, y, marker='o', label=family)
    ax.set_title(title)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def determine_resolution_mode(
    orient_success: bool,
    evolution_success: bool,
    grading_success: bool,
    recovery_success: bool,
) -> str:
    if orient_success and evolution_success and grading_success and recovery_success:
        return 'strong_resolution'
    if evolution_success and grading_success:
        return 'partial_structuralization'
    return 'negative_resolution'


def note_text(json_rel: str, csv_rel: str, rows: list[dict[str, Any]], plots: list[str], summary: dict[str, Any]) -> str:
    lines = [
        '# Stage 11 Integrated Austerity Resolution Probe v1',
        '',
        f'Timestamped JSON: `{json_rel}`',
        f'Timestamped CSV: `{csv_rel}`',
        '',
        f"Global read: `{summary['resolution_mode']}`",
        '',
        '## Obligation diagnostics',
        f"- `O_C0_ORIENT`: median positive-family count = `{summary['positive_alignment_families']}`; success = `{summary['orient_success']}`",
        f"- `O_C0_EVOLUTION`: median refine ratio = `{summary['median_refine_ratio']:.4f}`; success = `{summary['evolution_success']}`",
        f"- `O_C0_GRADING`: entropy coarse/fine = `{summary['mean_entropy_coarse']:.4f}` / `{summary['mean_entropy_fine']:.4f}`; success = `{summary['grading_success']}`",
        f"- `O_C0_RECOVERY`: ordered settings fraction = `{summary['ordered_fraction']:.4f}`; success = `{summary['recovery_success']}`",
        '',
        '## Family medians',
    ]
    for family_id, values in summary['family_medians'].items():
        lines.append(
            f"- `{family_id}`: ledger_alignment=`{values['ledger_alignment']:.4f}`, "
            f"delta_equiv_ratio=`{values['delta_equiv_refine_ratio']:.4f}`, "
            f"sector_entropy_fine=`{values['sector_entropy_fine']:.4f}`, "
            f"base_survival=`{values['base_survival_time']:.2f}`"
        )
    lines.extend(['', '## Per-run summary'])
    for row in rows:
        lines.append(
            f"- `{row['run_id']}`: family=`{row['family_label']}`, detuning=`{float(row['detuning']):.2f}`, "
            f"protect=`{float(row['protection_strength']):.2f}`, align=`{float(row['ledger_alignment']):.4f}`, "
            f"equiv_ratio=`{float(row['delta_equiv_refine_ratio']):.4f}`, entropy_fine=`{float(row['sector_entropy_fine']):.4f}`, "
            f"survival=`{float(row['base_survival_time']):.2f}`, final=`{row['final_topology_class']}`"
        )
    lines.extend(['', '## Plots'])
    for rel in plots:
        lines.append(f'- `{rel}`')
    return '\n'.join(lines) + '\n'


def main() -> None:
    args = parse_args()
    runsheet = load_runsheet(args.runsheet)
    common = runsheet['common_fields']
    selected: list[tuple[str, dict[str, Any], float, float]] = []
    for family in runsheet['families']:
        for detuning in runsheet['detuning_samples']:
            for protection in runsheet['protection_strengths']:
                run_id = f"S11_{family['family_id']}_d{int(round(detuning * 100)):03d}_p{int(round(protection * 100)):03d}"
                selected.append((run_id, family, float(detuning), float(protection)))
    if args.run_ids:
        wanted = set(args.run_ids)
        selected = [item for item in selected if item[0] in wanted]

    rows: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    topology_map = {'braid_like_exchange': 0, 'transfer_smeared': 1, 'localized_capture': 2, 'dispersive_pass': 3}

    for run_id, family, detuning, protection_strength in selected:
        G = matched_family_graph(str(family['family_id']))
        nodes = sorted(int(node) for node in G.nodes())
        triangles = triangle_list(G)
        tri_counts = node_triangle_counts(G, triangles)
        motif_nodes = motif_nodes_for_family(G, triangles, int(family['left_seed']), int(family['right_seed']))
        node_addresses, node_sectors = intrinsic_node_addresses(
            G=G,
            left_seed=int(family['left_seed']),
            right_seed=int(family['right_seed']),
            modulus=float(common['address_modulus']),
        )
        mismatch = np.zeros(len(nodes), dtype=float)
        left_packet = initial_packet(
            G=G,
            nodes=nodes,
            oriented_edges=orient_edges_from_mismatch(G, mismatch, nodes, tri_counts),
            triangles=triangles,
            seed=int(family['left_seed']),
            kernel_width=float(common['kernel_width']),
            direction_sign=+1.0,
            grade0_scale=float(common['grade0_scale']),
            grade1_scale=float(common['grade1_scale']),
            grade2_scale=float(common['grade2_scale']),
        )
        right_packet = initial_packet(
            G=G,
            nodes=nodes,
            oriented_edges=orient_edges_from_mismatch(G, mismatch, nodes, tri_counts),
            triangles=triangles,
            seed=int(family['right_seed']),
            kernel_width=float(common['kernel_width']),
            direction_sign=-1.0,
            grade0_scale=float(common['grade0_scale']),
            grade1_scale=float(common['grade1_scale']),
            grade2_scale=float(common['grade2_scale']),
        )
        psi = left_packet + right_packet
        psi /= max(float(np.linalg.norm(psi)), 1.0e-12)

        step_topology: list[str] = []
        step_topology_fine: list[str] = []
        dominant_bands: list[str] = []
        mismatch_gradients: list[float] = []
        protected_sequence_index: list[float] = []
        delta_equiv_coarse: list[float] = []
        delta_equiv_fine: list[float] = []
        topology_codes: list[int] = []
        coarse_states: list[np.ndarray] = [psi.copy()]
        fine_states: list[np.ndarray] = [psi.copy()]
        current_protected = 0.0
        path_nodes = nx.shortest_path(G, source=int(family['left_seed']), target=int(family['right_seed']))
        path_mask = np.asarray([1.0 if node in path_nodes else 0.0 for node in nodes], dtype=float)
        kernel_center = int(common['kernel_radius_center'])
        radius_values = [max(1, kernel_center + int(offset)) for offset in common['radius_offsets']]
        survival_counts = {'base': 0.0, 'skew': 0.0, 'motif': 0.0, 'rminus': 0.0, 'rcenter': 0.0, 'rplus': 0.0}

        for _step in range(int(common['steps'])):
            node_activity = node_activity_from_state(
                psi=psi,
                nodes=nodes,
                oriented_edges=orient_edges_from_mismatch(G, mismatch, nodes, tri_counts),
                triangles=triangles,
            )
            node_kernel = node_kernel_matrix(G, radius=kernel_center)
            motif_mask = np.asarray([1.0 if node in motif_nodes else 0.0 for node in nodes], dtype=float)
            forcing = node_kernel @ node_activity
            forcing += float(protection_strength) * (0.35 + 0.65 * current_protected) * (0.6 * path_mask + 0.4 * motif_mask)
            delta_mismatch = float(common['epsilon_mismatch']) * forcing
            mismatch = mismatch + delta_mismatch

            oriented_edges = orient_edges_from_mismatch(G, mismatch, nodes, tri_counts)
            d0, edge_lookup, _node_to_idx = build_d0(nodes, oriented_edges)
            d1 = build_d1(triangles, edge_lookup)
            D_signed, _Delta = block_dirac(d0, d1)
            supports = basis_supports(nodes, oriented_edges, triangles)
            coo = D_signed.tocoo()
            support_distances = support_distance_vector(G, supports, coo.row, coo.col)
            operator_rminus = radius_operator(D_signed, support_distances, radius_values[0])
            operator_rcenter = radius_operator(D_signed, support_distances, radius_values[1])
            operator_rplus = radius_operator(D_signed, support_distances, radius_values[2])
            labels = basis_labels_from_intrinsic(
                nodes=nodes,
                oriented_edges=oriented_edges,
                triangles=triangles,
                node_addresses=node_addresses,
                node_sectors=node_sectors,
                detuning=detuning,
                modulus=float(common['address_modulus']),
            )
            operator_base = address_weighted_operator(
                base_operator=operator_rcenter,
                labels=labels,
                modulus=float(common['address_modulus']),
                eta_match=float(common['eta_match']),
                eta_mismatch=float(common['eta_mismatch']),
            )
            basis_deg = basis_degrees(nodes, oriented_edges, triangles, G)
            operator_skew = apply_degree_skew(operator_base, basis_deg, float(common['degree_skew_strength']))
            operator_motif = apply_motif_boost(operator_base, supports, motif_nodes, float(common['motif_boost_strength']))
            operator_rminus = address_weighted_operator(
                base_operator=operator_rminus,
                labels=labels,
                modulus=float(common['address_modulus']),
                eta_match=float(common['eta_match']),
                eta_mismatch=float(common['eta_mismatch']),
            )
            operator_rplus = address_weighted_operator(
                base_operator=operator_rplus,
                labels=labels,
                modulus=float(common['address_modulus']),
                eta_match=float(common['eta_match']),
                eta_mismatch=float(common['eta_mismatch']),
            )

            psi_next = field_update(operator_base, psi, float(common['epsilon_field']))
            fine_mid = field_update(operator_base, psi, 0.5 * float(common['epsilon_field']))
            psi_next_fine = field_update(operator_base, fine_mid, 0.5 * float(common['epsilon_field']))

            n0 = len(nodes)
            n1 = len(oriented_edges)
            mismatch_projection = np.asarray(d0 @ delta_mismatch, dtype=float)
            coarse_flux = np.abs((psi_next[n0:n0 + n1] - psi[n0:n0 + n1]) / max(float(common['epsilon_field']), 1.0e-12))
            fine_flux = np.abs((psi_next_fine[n0:n0 + n1] - psi[n0:n0 + n1]) / max(float(common['epsilon_field']), 1.0e-12))
            delta_equiv_coarse.append(float(np.linalg.norm(np.abs(mismatch_projection) - coarse_flux)))
            delta_equiv_fine.append(float(np.linalg.norm(np.abs(mismatch_projection) - fine_flux)))

            coarse_states.append(psi_next.copy())
            fine_states.append(psi_next_fine.copy())
            states_window = coarse_states[max(0, len(coarse_states) - 12):]
            fine_window = fine_states[max(0, len(fine_states) - 12):]
            base_metrics = metrics_for_state(
                states_window=states_window,
                operator=operator_base,
                labels=labels,
                nodes=nodes,
                oriented_edges=oriented_edges,
                triangles=triangles,
                d1=d1,
                block_sizes=(len(nodes), len(oriented_edges), len(triangles)),
                modulus=float(common['address_modulus']),
                eta_match=float(common['eta_match']),
                eta_mismatch=float(common['eta_mismatch']),
            )
            fine_metrics = metrics_for_state(
                states_window=fine_window,
                operator=operator_base,
                labels=labels,
                nodes=nodes,
                oriented_edges=oriented_edges,
                triangles=triangles,
                d1=d1,
                block_sizes=(len(nodes), len(oriented_edges), len(triangles)),
                modulus=float(common['address_modulus']),
                eta_match=float(common['eta_match']),
                eta_mismatch=float(common['eta_mismatch']),
            )
            skew_metrics = metrics_for_state(
                states_window=states_window,
                operator=operator_skew,
                labels=labels,
                nodes=nodes,
                oriented_edges=oriented_edges,
                triangles=triangles,
                d1=d1,
                block_sizes=(len(nodes), len(oriented_edges), len(triangles)),
                modulus=float(common['address_modulus']),
                eta_match=float(common['eta_match']),
                eta_mismatch=float(common['eta_mismatch']),
            )
            motif_metrics = metrics_for_state(
                states_window=states_window,
                operator=operator_motif,
                labels=labels,
                nodes=nodes,
                oriented_edges=oriented_edges,
                triangles=triangles,
                d1=d1,
                block_sizes=(len(nodes), len(oriented_edges), len(triangles)),
                modulus=float(common['address_modulus']),
                eta_match=float(common['eta_match']),
                eta_mismatch=float(common['eta_mismatch']),
            )
            rminus_metrics = metrics_for_state(
                states_window=states_window,
                operator=operator_rminus,
                labels=labels,
                nodes=nodes,
                oriented_edges=oriented_edges,
                triangles=triangles,
                d1=d1,
                block_sizes=(len(nodes), len(oriented_edges), len(triangles)),
                modulus=float(common['address_modulus']),
                eta_match=float(common['eta_match']),
                eta_mismatch=float(common['eta_mismatch']),
            )
            rplus_metrics = metrics_for_state(
                states_window=states_window,
                operator=operator_rplus,
                labels=labels,
                nodes=nodes,
                oriented_edges=oriented_edges,
                triangles=triangles,
                d1=d1,
                block_sizes=(len(nodes), len(oriented_edges), len(triangles)),
                modulus=float(common['address_modulus']),
                eta_match=float(common['eta_match']),
                eta_mismatch=float(common['eta_mismatch']),
            )

            base_topology = str(base_metrics['topology_class'])
            fine_topology = str(fine_metrics['topology_class'])
            step_topology.append(base_topology)
            step_topology_fine.append(fine_topology)
            topology_codes.append(topology_map[base_topology])
            current_protected = 1.0 if base_topology in {'braid_like_exchange', 'transfer_smeared'} else 0.0
            protected_sequence_index.append((protected_sequence_index[-1] if protected_sequence_index else 0.0) + current_protected)
            edge_grad = np.asarray([abs(mismatch[nodes.index(u)] - mismatch[nodes.index(v)]) for u, v in oriented_edges], dtype=float)
            mismatch_gradients.append(float(np.mean(edge_grad)) if edge_grad.size else 0.0)

            if base_metrics['address_selectivity_index'] >= 0.92:
                dominant_bands.append('match')
            elif base_metrics['address_selectivity_index'] >= 0.75:
                dominant_bands.append('near')
            else:
                dominant_bands.append('mismatch')

            if base_topology in {'braid_like_exchange', 'transfer_smeared'}:
                survival_counts['base'] += 1.0
            if str(skew_metrics['topology_class']) in {'braid_like_exchange', 'transfer_smeared'}:
                survival_counts['skew'] += 1.0
            if str(motif_metrics['topology_class']) in {'braid_like_exchange', 'transfer_smeared'}:
                survival_counts['motif'] += 1.0
            if str(rminus_metrics['topology_class']) in {'braid_like_exchange', 'transfer_smeared'}:
                survival_counts['rminus'] += 1.0
            if str(base_metrics['topology_class']) in {'braid_like_exchange', 'transfer_smeared'}:
                survival_counts['rcenter'] += 1.0
            if str(rplus_metrics['topology_class']) in {'braid_like_exchange', 'transfer_smeared'}:
                survival_counts['rplus'] += 1.0

            psi = psi_next

        coarse_entropy = conditional_entropy(step_topology, dominant_bands)
        fine_entropy = conditional_entropy(step_topology_fine, dominant_bands)
        refine_ratio = float(np.mean(delta_equiv_fine) / max(float(np.mean(delta_equiv_coarse)), 1.0e-12))

        row = {
            'run_id': run_id,
            'family_id': str(family['family_id']),
            'family_label': str(family['family_label']),
            'detuning': float(detuning),
            'protection_strength': float(protection_strength),
            'ledger_alignment': correlation_or_zero(protected_sequence_index, mismatch_gradients),
            'delta_equiv_coarse_mean': float(np.mean(delta_equiv_coarse)),
            'delta_equiv_fine_mean': float(np.mean(delta_equiv_fine)),
            'delta_equiv_refine_ratio': refine_ratio,
            'sector_entropy_coarse': coarse_entropy,
            'sector_entropy_fine': fine_entropy,
            'base_survival_time': float(survival_counts['base']),
            'degree_skew_survival_time': float(survival_counts['skew']),
            'motif_survival_time': float(survival_counts['motif']),
            'radius_minus_survival_time': float(survival_counts['rminus']),
            'radius_center_survival_time': float(survival_counts['rcenter']),
            'radius_plus_survival_time': float(survival_counts['rplus']),
            'ordering_rank': 0.0,
            'resolution_mode': '',
            'final_topology_class': step_topology[-1] if step_topology else 'dispersive_pass',
            'notes': '',
        }
        trace_path = WORK_PLOT_DIR / f'{run_id}_trace.png'
        trace_plot(trace_path, run_id, topology_codes, mismatch_gradients, delta_equiv_coarse)
        row['trace_plot'] = str(trace_path)
        rows.append(row)
        plot_paths.append(trace_path)

    grouped_family: dict[str, list[dict[str, Any]]] = defaultdict(list)
    grouped_setting: dict[tuple[float, float], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped_family[str(row['family_id'])].append(row)
        grouped_setting[(float(row['detuning']), float(row['protection_strength']))].append(row)

    family_medians: dict[str, dict[str, float]] = {}
    for family_id, family_rows in grouped_family.items():
        family_medians[family_id] = {
            'ledger_alignment': float(np.median([float(item['ledger_alignment']) for item in family_rows])),
            'delta_equiv_refine_ratio': float(np.median([float(item['delta_equiv_refine_ratio']) for item in family_rows])),
            'sector_entropy_fine': float(np.median([float(item['sector_entropy_fine']) for item in family_rows])),
            'base_survival_time': float(np.median([float(item['base_survival_time']) for item in family_rows])),
        }

    positive_alignment_families = sum(1 for values in family_medians.values() if values['ledger_alignment'] > 0.15)
    orient_success = positive_alignment_families >= 2

    median_refine_ratio = float(np.median([float(row['delta_equiv_refine_ratio']) for row in rows])) if rows else 0.0
    evolution_success = 0.5 <= median_refine_ratio <= 1.5

    mean_entropy_coarse = float(np.mean([float(row['sector_entropy_coarse']) for row in rows])) if rows else 0.0
    mean_entropy_fine = float(np.mean([float(row['sector_entropy_fine']) for row in rows])) if rows else 0.0
    entropy_by_detuning: dict[float, float] = {}
    for detuning in sorted({float(row['detuning']) for row in rows}):
        entropy_by_detuning[detuning] = float(np.mean([float(row['sector_entropy_fine']) for row in rows if float(row['detuning']) == detuning]))
    grading_success = mean_entropy_fine <= mean_entropy_coarse and entropy_by_detuning.get(1.0, mean_entropy_fine) <= entropy_by_detuning.get(0.0, mean_entropy_fine)

    ordered_settings = 0
    for setting, setting_rows in grouped_setting.items():
        lookup = {str(item['family_id']): float(item['base_survival_time']) for item in setting_rows}
        if (
            lookup.get('clustered_braid_seed', -1.0) >= lookup.get('corridor_transport_seed', -1.0)
            and lookup.get('corridor_transport_seed', -1.0) >= lookup.get('triad_competition_seed', -1.0)
        ):
            ordered_settings += 1
        for item in setting_rows:
            item['ordering_rank'] = 1.0 if (
                lookup.get('clustered_braid_seed', -1.0) >= lookup.get('corridor_transport_seed', -1.0)
                and lookup.get('corridor_transport_seed', -1.0) >= lookup.get('triad_competition_seed', -1.0)
            ) else 0.0
    ordered_fraction = float(ordered_settings / max(len(grouped_setting), 1))
    recovery_success = ordered_fraction >= 0.70

    resolution_mode = determine_resolution_mode(orient_success, evolution_success, grading_success, recovery_success)
    for row in rows:
        row['resolution_mode'] = resolution_mode

    alignment_plot = WORK_PLOT_DIR / 'stage_11_alignment_panel.png'
    entropy_plot = WORK_PLOT_DIR / 'stage_11_entropy_panel.png'
    survival_plot = WORK_PLOT_DIR / 'stage_11_survival_panel.png'
    summary_plot(alignment_plot, rows, 'ledger_alignment', 'Stage 11 ledger alignment by run')
    summary_plot(entropy_plot, rows, 'sector_entropy_fine', 'Stage 11 fine sector entropy by run')
    summary_plot(survival_plot, rows, 'base_survival_time', 'Stage 11 protected survival by run')
    plot_paths.extend([alignment_plot, entropy_plot, survival_plot])

    summary = {
        'positive_alignment_families': positive_alignment_families,
        'orient_success': orient_success,
        'median_refine_ratio': median_refine_ratio,
        'evolution_success': evolution_success,
        'mean_entropy_coarse': mean_entropy_coarse,
        'mean_entropy_fine': mean_entropy_fine,
        'grading_success': grading_success,
        'ordered_fraction': ordered_fraction,
        'recovery_success': recovery_success,
        'family_medians': family_medians,
        'resolution_mode': resolution_mode,
    }

    result_payload = {
        'stage': 'stage_11_integrated_austerity_resolution_probe',
        'runsheet': str(Path(args.runsheet).relative_to(REPO_ROOT)),
        'rows': rows,
        'summary': summary,
    }
    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug='stage_11_integrated_austerity_resolution_probe',
        result=result_payload,
        csv_rows=rows,
        csv_fieldnames=CSV_FIELDS,
        plot_paths=plot_paths,
    )
    NOTE_PATH.write_text(
        note_text(
            str(json_path.relative_to(REPO_ROOT)),
            str(csv_path.relative_to(REPO_ROOT)),
            rows,
            stamped_plots,
            summary,
        ),
        encoding='utf-8',
    )


if __name__ == '__main__':
    main()
