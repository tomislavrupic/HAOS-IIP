#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict, deque
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import ATLAS_NOTES, REPO_ROOT, load_stage10_defaults, plt, save_atlas_payload
from stage11_collective_wave_interaction import build_component_packet
from stage_c0_1_combinatorial_kernel_probe import build_line_graph_cache, get_cached_setup

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_5_emergent_compression_runs.json'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage_c0_5_emergent_compression_preorder'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'run_id',
    'representative_id',
    'compression_mode_id',
    'compression_mode_label',
    'reference_kernel_class_id',
    'reference_interaction_label',
    'reference_coarse_label',
    'resolution',
    'initial_active_support',
    'initial_seed_count',
    'steps_completed',
    'final_survivor_count',
    'final_survivor_fraction',
    'final_equivalence_class_count',
    'final_largest_class_size',
    'tau_stabilization_depth',
    'delta_quiet_step',
    'z2_arrow_step',
    'persistent_step',
    'obstruction_rank_final',
    'nonfunctoriality_initial',
    'nonfunctoriality_final',
    'max_defect_density',
    'final_defect_density',
    'pass_blocked_merge_total',
    'defect_ledger_zero_before_pass',
    'compression_class',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.5 emergent compression preorder scan.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def selected_runs(runs: list[dict[str, Any]], run_ids: list[str]) -> list[dict[str, Any]]:
    if not run_ids:
        return runs
    wanted = set(run_ids)
    return [run for run in runs if run['run_id'] in wanted]


def lookup_by_id(items: list[dict[str, Any]], key: str) -> dict[str, dict[str, Any]]:
    return {str(item[key]): item for item in items}


def compression_mode_order(mode_id: str) -> int:
    return {
        'adjacency_preorder': 0,
        'graph_shell_preorder': 1,
        'alternating_pass_preorder': 2,
    }.get(mode_id, 99)


def load_c01_reference_rows() -> dict[tuple[str, str], dict[str, str]]:
    pattern = REPO_ROOT / 'data'
    candidates = sorted(pattern.glob('*_stage_c0_1_combinatorial_kernel_probe.csv'))
    if not candidates:
        return {}
    rows: dict[tuple[str, str], dict[str, str]] = {}
    with candidates[-1].open('r', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            key = (str(row['representative_id']), str(row['kernel_class_id']))
            rows[key] = {
                'interaction_label': str(row['interaction_label']),
                'coarse_label': str(row['coarse_label']),
            }
    return rows


def initial_edge_weights(
    rep: dict[str, Any],
    resolution: int,
    boundary_type: str,
    defaults: dict[str, Any],
    base: dict[str, Any],
) -> tuple[Any, np.ndarray]:
    data, projector = get_cached_setup(resolution, defaults, boundary_type)
    packet_qs: list[np.ndarray] = []
    for center, amp, phase, kick_sign in zip(
        rep['packet_centers'],
        rep['packet_amplitudes'],
        rep['phase_offsets_rad'],
        rep['kick_signs'],
    ):
        q0, _v0 = build_component_packet(
            center=np.asarray(center, dtype=float),
            amplitude=float(amp),
            phase_offset=float(phase),
            kick_sign=float(kick_sign),
            bandwidth=float(base['bandwidth']),
            central_k=float(base['central_k']),
            phase_pattern=str(base['phase_pattern']),
            kick_axis=int(defaults['kick_axis']),
            boundary_type=boundary_type,
            midpoints=data.midpoints,
            edge_axes=data.edge_axes,
            projector=projector,
        )
        packet_qs.append(q0)
    total_q = np.sum(packet_qs, axis=0)
    return data, np.abs(total_q) ** 2


def select_active_support(weights: np.ndarray, sigma_scale: float, min_active_count: int) -> list[int]:
    mu = float(np.mean(weights))
    sigma = float(np.std(weights))
    threshold = mu + sigma_scale * sigma
    active = np.flatnonzero(weights >= threshold).tolist()
    if len(active) < min_active_count:
        order = np.argsort(weights)[::-1]
        active = [int(idx) for idx in order[:min_active_count]]
    return sorted(set(int(idx) for idx in active))


def top_seed_edges(weights: np.ndarray, support: list[int], seed_count: int) -> list[int]:
    order = sorted(support, key=lambda idx: float(weights[idx]), reverse=True)
    return [int(idx) for idx in order[: min(seed_count, len(order))]]


def adjacency_sets(adjacency_lists: list[list[int]]) -> list[set[int]]:
    return [set(items) for items in adjacency_lists]


def multi_source_shell_distance(adjacency_lists: list[list[int]], seed_nodes: list[int], cap: int) -> np.ndarray:
    distance = np.full(len(adjacency_lists), cap + 1, dtype=int)
    if not seed_nodes:
        return distance
    queue: deque[int] = deque()
    for node in seed_nodes:
        distance[int(node)] = 0
        queue.append(int(node))
    while queue:
        cur = queue.popleft()
        if distance[cur] >= cap:
            continue
        for nbr in adjacency_lists[cur]:
            if distance[nbr] <= distance[cur] + 1:
                continue
            distance[nbr] = distance[cur] + 1
            queue.append(nbr)
    return distance


def canonical_partition(classes: list[list[int]] | list[tuple[int, ...]]) -> tuple[tuple[int, ...], ...]:
    cleaned = [tuple(sorted(int(value) for value in cls)) for cls in classes if cls]
    return tuple(sorted(cleaned, key=lambda cls: (cls[0], len(cls))))


def partition_map(partition: tuple[tuple[int, ...], ...]) -> dict[int, tuple[int, ...]]:
    mapping: dict[int, tuple[int, ...]] = {}
    for cls in partition:
        for node in cls:
            mapping[int(node)] = cls
    return mapping


def survivor_partition(partition: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(cls for cls in partition if len(cls) > 1)


def survivor_support(partition: tuple[tuple[int, ...], ...]) -> list[int]:
    nodes: list[int] = []
    for cls in partition:
        if len(cls) > 1:
            nodes.extend(int(node) for node in cls)
    return sorted(nodes)


def focused_partition(
    partition: tuple[tuple[int, ...], ...],
    focus_nodes: list[int],
) -> tuple[tuple[int, ...], ...]:
    focus_set = set(int(node) for node in focus_nodes)
    if focus_set:
        classes = [cls for cls in partition if len(cls) > 1 and any(int(node) in focus_set for node in cls)]
        if classes:
            return canonical_partition(list(classes))
    return survivor_partition(partition)


def restricted_neighbors(support: list[int], full_adj_sets: list[set[int]]) -> dict[int, set[int]]:
    support_set = set(int(node) for node in support)
    return {int(node): set(full_adj_sets[int(node)]).intersection(support_set) for node in support}


def local_triangle_counts(restricted: dict[int, set[int]]) -> dict[int, int]:
    counts: dict[int, int] = {}
    for node, nbrs in restricted.items():
        nbr_list = sorted(nbrs)
        tri = 0
        for idx, u in enumerate(nbr_list[:-1]):
            nbr_u = restricted[u]
            for v in nbr_list[idx + 1 :]:
                if v in nbr_u:
                    tri += 1
        counts[int(node)] = int(tri)
    return counts


def local_square_counts(restricted: dict[int, set[int]]) -> dict[int, int]:
    counts: dict[int, int] = {}
    for node, nbrs in restricted.items():
        nbr_list = sorted(nbrs)
        square = 0
        for idx, u in enumerate(nbr_list[:-1]):
            nbr_u = restricted[u]
            for v in nbr_list[idx + 1 :]:
                commons = nbr_u.intersection(restricted[v]) - {int(node)}
                square += len(commons)
        counts[int(node)] = int(square // 2)
    return counts


def quantize_metric(value: float, scale: float, levels: int = 12) -> int:
    if scale <= 1.0e-12:
        return 0
    clipped = max(0.0, min(1.0, float(value) / scale))
    return int(round(levels * clipped))


def pairwise_partition_mismatch(
    support: list[int],
    adjacency_map: dict[int, tuple[int, ...]],
    shell_map: dict[int, tuple[int, ...]],
) -> tuple[int, dict[int, int]]:
    nodes = [int(node) for node in support]
    contributions = {int(node): 0 for node in nodes}
    mismatch = 0
    for idx, u in enumerate(nodes[:-1]):
        for v in nodes[idx + 1 :]:
            same_adj = adjacency_map.get(u) == adjacency_map.get(v)
            same_shell = shell_map.get(u) == shell_map.get(v)
            if same_adj != same_shell:
                mismatch += 1
                contributions[u] += 1
                contributions[v] += 1
    return mismatch, contributions


def line_graph_cycle_rank(restricted: dict[int, set[int]]) -> int:
    nodes = list(restricted.keys())
    if not nodes:
        return 0
    visited: set[int] = set()
    components = 0
    for start in nodes:
        if start in visited:
            continue
        components += 1
        queue: deque[int] = deque([start])
        visited.add(start)
        while queue:
            cur = queue.popleft()
            for nbr in restricted[cur]:
                if nbr not in visited:
                    visited.add(nbr)
                    queue.append(nbr)
    edge_count = sum(len(nbrs) for nbrs in restricted.values()) // 2
    return max(edge_count - len(nodes) + components, 0)


def analyze_support(
    support: list[int],
    weights: np.ndarray,
    full_adj_sets: list[set[int]],
    full_shells: list[dict[int, int]],
    seed_shell_distance: np.ndarray,
    params: dict[str, Any],
) -> dict[str, Any]:
    support = sorted(set(int(node) for node in support))
    if not support:
        empty_partition: tuple[tuple[int, ...], ...] = tuple()
        return {
            'support': [],
            'support_size': 0,
            'restricted': {},
            'motif_signatures': {},
            'adjacency_partition': empty_partition,
            'graph_shell_partition': empty_partition,
            'adjacency_map': {},
            'graph_shell_map': {},
            'adjacency_survivor_count': 0,
            'graph_shell_survivor_count': 0,
            'adjacency_equivalence_count': 0,
            'graph_shell_equivalence_count': 0,
            'nonfunctoriality_count': 0,
            'nonfunctoriality_ratio': 0.0,
            'degree_mean': 0.0,
            'cycle_rank': 0,
        }

    max_shell = int(params.get('max_shell', 3))
    graph_shell_weights = {
        int(key): float(value)
        for key, value in dict(params.get('graph_shell_weights', {'1': 1.0, '2': 0.5, '3': 0.25})).items()
    }
    restricted = restricted_neighbors(support, full_adj_sets)
    triangle_counts = local_triangle_counts(restricted)
    square_counts = local_square_counts(restricted)
    degrees = {int(node): len(restricted[int(node)]) for node in support}
    degree_values = np.asarray([degrees[int(node)] for node in support], dtype=float)

    adjacency_mass: dict[int, float] = {}
    shell_mass: dict[int, float] = {}
    shell_counts: dict[int, tuple[int, int, int]] = {}
    motif_signatures: dict[int, tuple[int, ...]] = {}
    adjacency_signatures: dict[int, tuple[int, ...]] = {}
    shell_signatures: dict[int, tuple[int, ...]] = {}

    support_set = set(support)
    max_weight = float(np.max(weights[support]))
    for node in support:
        adj_mass = float(np.sum([weights[nbr] for nbr in restricted[node]])) if restricted[node] else 0.0
        shell_count_values = [0, 0, 0]
        weighted_shell_mass = 0.0
        for nbr, dist in full_shells[int(node)].items():
            if nbr not in support_set or dist < 1 or dist > max_shell:
                continue
            shell_count_values[dist - 1] += 1
            weighted_shell_mass += graph_shell_weights.get(dist, 0.0) * float(weights[nbr])
        adjacency_mass[int(node)] = adj_mass
        shell_mass[int(node)] = weighted_shell_mass
        shell_counts[int(node)] = tuple(int(value) for value in shell_count_values)

    max_adj_mass = max(adjacency_mass.values()) if adjacency_mass else 0.0
    max_shell_mass = max(shell_mass.values()) if shell_mass else 0.0
    for node in support:
        activity_bin = quantize_metric(float(weights[node]), max_weight)
        adj_mass_bin = quantize_metric(float(adjacency_mass[node]), max_adj_mass)
        shell_mass_bin = quantize_metric(float(shell_mass[node]), max_shell_mass)
        seed_dist = int(min(int(seed_shell_distance[int(node)]), max_shell + 1))
        tri = int(triangle_counts[int(node)])
        square = int(square_counts[int(node)])
        deg = int(degrees[int(node)])
        shell1, shell2, shell3 = shell_counts[int(node)]
        motif_sig = (seed_dist, deg, tri, square, shell1, shell2, shell3)
        motif_signatures[int(node)] = motif_sig
        adjacency_signatures[int(node)] = motif_sig + (activity_bin, adj_mass_bin)
        shell_signatures[int(node)] = (seed_dist, shell1, shell2, shell3, deg, tri, square, activity_bin, shell_mass_bin)

    adjacency_groups: dict[tuple[int, ...], list[int]] = defaultdict(list)
    shell_groups: dict[tuple[int, ...], list[int]] = defaultdict(list)
    for node in support:
        adjacency_groups[adjacency_signatures[int(node)]].append(int(node))
        shell_groups[shell_signatures[int(node)]].append(int(node))

    adjacency_partition = canonical_partition(list(adjacency_groups.values()))
    graph_shell_partition = canonical_partition(list(shell_groups.values()))
    adjacency_map = partition_map(adjacency_partition)
    graph_shell_map = partition_map(graph_shell_partition)
    nonfun_count, _contrib = pairwise_partition_mismatch(support, adjacency_map, graph_shell_map)
    pair_total = max(len(support) * max(len(support) - 1, 0) // 2, 1)

    return {
        'support': support,
        'support_size': len(support),
        'restricted': restricted,
        'motif_signatures': motif_signatures,
        'adjacency_partition': adjacency_partition,
        'graph_shell_partition': graph_shell_partition,
        'adjacency_map': adjacency_map,
        'graph_shell_map': graph_shell_map,
        'adjacency_survivor_count': len(survivor_support(adjacency_partition)),
        'graph_shell_survivor_count': len(survivor_support(graph_shell_partition)),
        'adjacency_equivalence_count': sum(1 for cls in adjacency_partition if len(cls) > 1),
        'graph_shell_equivalence_count': sum(1 for cls in graph_shell_partition if len(cls) > 1),
        'nonfunctoriality_count': int(nonfun_count),
        'nonfunctoriality_ratio': float(nonfun_count / pair_total),
        'degree_mean': float(np.mean(degree_values)) if degree_values.size else 0.0,
        'cycle_rank': int(line_graph_cycle_rank(restricted)),
    }


def apply_pass_separation(
    raw_partition: tuple[tuple[int, ...], ...],
    inherited_partition: tuple[tuple[int, ...], ...],
) -> tuple[tuple[tuple[int, ...], ...], int]:
    if not inherited_partition:
        return raw_partition, 0
    inherited_map = partition_map(inherited_partition)
    merged_classes = 0
    next_classes: list[list[int]] = []
    for raw_cls in raw_partition:
        by_inherited: dict[tuple[int, ...], list[int]] = defaultdict(list)
        for node in raw_cls:
            key = inherited_map.get(int(node), (int(node),))
            by_inherited[key].append(int(node))
        if len(by_inherited) > 1:
            merged_classes += 1
        next_classes.extend(list(by_inherited.values()))
    return canonical_partition(next_classes), merged_classes


def next_partition_for_kernel(state: dict[str, Any], kernel_name: str) -> tuple[tuple[int, ...], ...]:
    if kernel_name == 'adjacency':
        return state['adjacency_partition']
    if kernel_name == 'graph_shell':
        return state['graph_shell_partition']
    raise ValueError(f'unsupported kernel name: {kernel_name}')


def class_map_for_kernel(state: dict[str, Any], kernel_name: str) -> dict[int, tuple[int, ...]]:
    if kernel_name == 'adjacency':
        return state['adjacency_map']
    if kernel_name == 'graph_shell':
        return state['graph_shell_map']
    raise ValueError(f'unsupported kernel name: {kernel_name}')


def evaluate_candidate(
    state: dict[str, Any],
    current_partition: tuple[tuple[int, ...], ...],
    focus_nodes: list[int],
    kernel_name: str,
    pass_active: bool,
    weights: np.ndarray,
    full_adj_sets: list[set[int]],
    full_shells: list[dict[int, int]],
    seed_shell_distance: np.ndarray,
    params: dict[str, Any],
) -> dict[str, Any]:
    raw_partition = next_partition_for_kernel(state, kernel_name)
    blocked_merges = 0
    final_partition = raw_partition
    if pass_active:
        final_partition, blocked_merges = apply_pass_separation(raw_partition, current_partition)
    final_partition = focused_partition(final_partition, focus_nodes)
    next_support = survivor_support(final_partition)
    next_state = analyze_support(next_support, weights, full_adj_sets, full_shells, seed_shell_distance, params)
    return {
        'kernel_name': kernel_name,
        'raw_partition': raw_partition,
        'final_partition': survivor_partition(final_partition),
        'next_support': next_support,
        'next_state': next_state,
        'blocked_merges': int(blocked_merges),
        'defect_reduction': int(state['nonfunctoriality_count']) - int(next_state['nonfunctoriality_count']),
    }


def defect_ledger(
    prev_state: dict[str, Any],
    prev_partition: tuple[tuple[int, ...], ...],
    next_state: dict[str, Any],
    kernel_name: str,
) -> list[int]:
    next_support = set(int(node) for node in next_state['support'])
    if not next_support or not prev_partition:
        return []
    restricted_prev = canonical_partition(
        [
            tuple(sorted(node for node in cls if int(node) in next_support))
            for cls in prev_partition
            if any(int(node) in next_support for node in cls)
        ]
    )
    prev_map = partition_map(restricted_prev)
    next_map = class_map_for_kernel(next_state, kernel_name)
    ledger: list[int] = []
    for node in sorted(next_support):
        prev_cls = prev_map.get(int(node))
        next_cls = next_map.get(int(node))
        if prev_cls is None or next_cls is None:
            continue
        motif_prev = prev_state['motif_signatures'].get(int(node))
        motif_next = next_state['motif_signatures'].get(int(node))
        if motif_prev != motif_next and prev_cls == next_cls:
            ledger.append(int(node))
    return ledger


def classify_compression(summary: dict[str, Any]) -> str:
    if int(summary['defect_ledger_zero_before_pass']) == 1 and int(summary['z2_arrow_step']) < 0:
        return 'static_regime'
    if int(summary['z2_arrow_step']) >= 0 and int(summary['persistent_step']) >= 0:
        return 'pass_ordered_persistent'
    if int(summary['z2_arrow_step']) >= 0:
        return 'pass_ordered_transient'
    if int(summary['persistent_step']) >= 0:
        return 'compression_stable_without_arrow'
    return 'unresolved_preorder'


def survivor_overlay_plot(
    path: Path,
    run_id: str,
    midpoints: np.ndarray,
    weights: np.ndarray,
    initial_support: list[int],
    final_support: list[int],
    seed_edges: list[int],
    boundary_type: str,
) -> None:
    coords = np.asarray(midpoints, dtype=float)
    slice_value = 0.5 if boundary_type == 'periodic' else 0.5
    distances = np.abs(coords[:, 2] - slice_value)
    threshold = np.min(distances) + 1.0e-9
    mask = distances <= threshold

    initial_mask = np.zeros(len(coords), dtype=bool)
    initial_mask[np.asarray(initial_support, dtype=int)] = True
    final_mask = np.zeros(len(coords), dtype=bool)
    if final_support:
        final_mask[np.asarray(final_support, dtype=int)] = True
    seed_mask = np.zeros(len(coords), dtype=bool)
    if seed_edges:
        seed_mask[np.asarray(seed_edges, dtype=int)] = True

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    sc = ax.scatter(coords[mask, 0], coords[mask, 1], c=weights[mask], s=14, cmap='Greys', alpha=0.25)
    if np.any(mask & initial_mask):
        ax.scatter(
            coords[mask & initial_mask, 0],
            coords[mask & initial_mask, 1],
            c=weights[mask & initial_mask],
            s=24,
            cmap='viridis',
            alpha=0.65,
            label='initial active support',
        )
    if np.any(mask & final_mask):
        ax.scatter(
            coords[mask & final_mask, 0],
            coords[mask & final_mask, 1],
            s=52,
            facecolors='none',
            edgecolors='tab:red',
            linewidths=1.3,
            label='final survivors',
        )
    if np.any(mask & seed_mask):
        ax.scatter(
            coords[mask & seed_mask, 0],
            coords[mask & seed_mask, 1],
            marker='x',
            s=70,
            c='tab:orange',
            linewidths=1.8,
            label='seed edges',
        )
    ax.set_title(f'Survivor support: {run_id}')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.legend(fontsize=8)
    plt.colorbar(sc, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def preorder_trace_plot(path: Path, run_id: str, trace: list[dict[str, Any]]) -> None:
    steps = [int(item['step_index']) for item in trace]
    w_values = [int(item['survivor_count']) for item in trace]
    largest = [int(item['largest_class_size']) for item in trace]
    defect = [float(item['defect_density']) for item in trace]
    nonfun = [int(item['nonfunctoriality_after']) for item in trace]
    rho = [int(item['obstruction_rank_after']) for item in trace]
    blocked = [int(item['pass_blocked_merges']) for item in trace]

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.2))
    axes[0].plot(steps, w_values, marker='o', color='tab:blue', label='W survivor count')
    axes[0].plot(steps, largest, marker='s', color='tab:orange', label='largest class')
    axes[0].set_title('Behavioral-congruence survivors')
    axes[0].set_xlabel('compression step')
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=8)

    axes[1].plot(steps, defect, marker='o', color='tab:red', label='defect density')
    axes[1].plot(steps, nonfun, marker='s', color='tab:purple', label='non-functoriality')
    axes[1].set_title('Defect ledger and mismatch count')
    axes[1].set_xlabel('compression step')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)

    axes[2].plot(steps, rho, marker='o', color='tab:green', label='cycle rank')
    axes[2].bar(steps, blocked, alpha=0.35, color='tab:gray', label='PASS blocks')
    axes[2].set_title('Obstruction rank and PASS blocks')
    axes[2].set_xlabel('compression step')
    axes[2].grid(alpha=0.25)
    axes[2].legend(fontsize=8)

    fig.suptitle(run_id)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def simulate_run(
    run: dict[str, Any],
    rep: dict[str, Any],
    compression_mode: dict[str, Any],
    runsheet: dict[str, Any],
    defaults: dict[str, Any],
    base: dict[str, Any],
    c01_reference: dict[tuple[str, str], dict[str, str]],
) -> tuple[dict[str, Any], list[Path]]:
    params = dict(runsheet.get('parameter_policy', {}))
    resolution = int(run['resolution'])
    boundary_type = str(run['boundary_type'])
    mode_id = str(run['compression_mode_id'])

    data, weights = initial_edge_weights(rep, resolution, boundary_type, defaults, base)
    line_graph = build_line_graph_cache(data, max_shell=int(params.get('max_shell', 3)))
    full_adj_lists = line_graph['adjacency_lists']
    full_adj_sets = adjacency_sets(full_adj_lists)
    full_shells = line_graph['shells']

    initial_support = select_active_support(
        weights=weights,
        sigma_scale=float(params.get('active_threshold_sigma', 1.0)),
        min_active_count=int(params.get('min_active_count', 48)),
    )
    seed_edges = top_seed_edges(weights, initial_support, int(params.get('seed_edge_count', 12)))
    seed_shell_distance = multi_source_shell_distance(
        full_adj_lists,
        seed_edges,
        cap=int(params.get('seed_shell_cap', 4)),
    )

    initial_state = analyze_support(initial_support, weights, full_adj_sets, full_shells, seed_shell_distance, params)
    start_kernel = 'adjacency' if mode_id == 'adjacency_preorder' else 'graph_shell'
    start_partition = survivor_partition(next_partition_for_kernel(initial_state, start_kernel))
    current_support = survivor_support(start_partition)
    current_partition = start_partition
    current_state = analyze_support(current_support, weights, full_adj_sets, full_shells, seed_shell_distance, params)
    current_focus_nodes = [int(node) for node in seed_edges if int(node) in current_support] or list(current_support)

    trace: list[dict[str, Any]] = []
    arrow_step: int | None = None
    tau_step: int | None = None
    persistent_step: int | None = None
    delta_quiet_step: int | None = None
    defect_zero_before_pass = False
    pass_blocked_total = 0
    last_cycle_rank: int | None = None
    stable_cycle_count = 0
    defect_threshold = float(params.get('defect_density_threshold', 0.004))
    rho_window = int(params.get('rho_stability_window', 2))
    max_steps = int(params.get('max_steps', 8))
    current_kernel = start_kernel

    initial_ledger = defect_ledger(initial_state, start_partition, current_state, start_kernel)
    initial_delta = float(len(initial_ledger) / max(len(full_adj_lists), 1))
    trace.append(
        {
            'step_index': 0,
            'kernel_name': start_kernel,
            'support_before': int(initial_state['support_size']),
            'survivor_count': int(current_state['support_size']),
            'equivalence_class_count': int(sum(1 for cls in current_partition if len(cls) > 1)),
            'largest_class_size': int(max((len(cls) for cls in current_partition), default=0)),
            'nonfunctoriality_before': int(initial_state['nonfunctoriality_count']),
            'nonfunctoriality_after': int(current_state['nonfunctoriality_count']),
            'defect_ledger_size': int(len(initial_ledger)),
            'defect_density': initial_delta,
            'pass_blocked_merges': 0,
            'obstruction_rank_after': int(current_state['cycle_rank']),
        }
    )
    if initial_delta <= defect_threshold:
        delta_quiet_step = 0
    if len(initial_ledger) == 0:
        defect_zero_before_pass = True
    last_cycle_rank = int(current_state['cycle_rank'])
    stable_cycle_count = 1
    if rho_window <= 1:
        persistent_step = 0

    for step_idx in range(1, max_steps + 1):
        if not current_support:
            break
        if mode_id == 'adjacency_preorder':
            candidate = evaluate_candidate(
                current_state,
                current_partition,
                current_focus_nodes,
                'adjacency',
                False,
                weights,
                full_adj_sets,
                full_shells,
                seed_shell_distance,
                params,
            )
        elif mode_id == 'graph_shell_preorder':
            candidate = evaluate_candidate(
                current_state,
                current_partition,
                current_focus_nodes,
                'graph_shell',
                False,
                weights,
                full_adj_sets,
                full_shells,
                seed_shell_distance,
                params,
            )
        else:
            adjacency_candidate = evaluate_candidate(
                current_state,
                current_partition,
                current_focus_nodes,
                'adjacency',
                True,
                weights,
                full_adj_sets,
                full_shells,
                seed_shell_distance,
                params,
            )
            shell_candidate = evaluate_candidate(
                current_state,
                current_partition,
                current_focus_nodes,
                'graph_shell',
                True,
                weights,
                full_adj_sets,
                full_shells,
                seed_shell_distance,
                params,
            )
            candidate = sorted(
                [adjacency_candidate, shell_candidate],
                key=lambda item: (
                    -int(item['defect_reduction']),
                    int(item['next_state']['nonfunctoriality_count']),
                    int(item['blocked_merges']),
                    int(item['next_state']['support_size']),
                    0 if item['kernel_name'] == current_kernel else 1,
                ),
            )[0]

        next_support = list(candidate['next_support'])
        next_partition = tuple(tuple(int(node) for node in cls) for cls in candidate['final_partition'])
        next_state = candidate['next_state']
        ledger = defect_ledger(current_state, current_partition, next_state, str(candidate['kernel_name']))
        defect_density = float(len(ledger) / max(len(full_adj_lists), 1))
        pass_blocked_total += int(candidate['blocked_merges'])
        if int(candidate['blocked_merges']) > 0 and arrow_step is None:
            arrow_step = int(step_idx)
            defect_zero_before_pass = False
        if delta_quiet_step is None and defect_density <= defect_threshold:
            delta_quiet_step = int(step_idx)
        if arrow_step is None and len(ledger) == 0:
            defect_zero_before_pass = True
        if last_cycle_rank == int(next_state['cycle_rank']):
            stable_cycle_count += 1
        else:
            stable_cycle_count = 1
        last_cycle_rank = int(next_state['cycle_rank'])
        if persistent_step is None and stable_cycle_count >= rho_window:
            persistent_step = int(step_idx)
        if tau_step is None and next_support == current_support and next_partition == current_partition:
            tau_step = int(step_idx)

        trace.append(
            {
                'step_index': int(step_idx),
                'kernel_name': str(candidate['kernel_name']),
                'support_before': int(current_state['support_size']),
                'survivor_count': int(len(next_support)),
                'equivalence_class_count': int(sum(1 for cls in next_partition if len(cls) > 1)),
                'largest_class_size': int(max((len(cls) for cls in next_partition), default=0)),
                'nonfunctoriality_before': int(current_state['nonfunctoriality_count']),
                'nonfunctoriality_after': int(next_state['nonfunctoriality_count']),
                'defect_ledger_size': int(len(ledger)),
                'defect_density': defect_density,
                'pass_blocked_merges': int(candidate['blocked_merges']),
                'obstruction_rank_after': int(next_state['cycle_rank']),
            }
        )

        stabilized = next_support == current_support and next_partition == current_partition
        current_support = next_support
        current_partition = next_partition
        current_state = next_state
        current_kernel = str(candidate['kernel_name'])
        current_focus_nodes = [int(node) for node in seed_edges if int(node) in current_support] or list(current_support)
        if stabilized:
            break

    if tau_step is None:
        tau_step = max((int(item['step_index']) for item in trace), default=0)

    if not current_support:
        final_equivalence_count = 0
        final_largest_class = 0
    else:
        final_equivalence_count = int(sum(1 for cls in current_partition if len(cls) > 1))
        final_largest_class = int(max((len(cls) for cls in current_partition), default=0))

    reference_kernel = 'pure_adjacency' if mode_id == 'adjacency_preorder' else 'graph_shell'
    reference = c01_reference.get((str(rep['representative_id']), reference_kernel), {})
    summary = {
        'run_id': str(run['run_id']),
        'representative_id': str(run['representative_id']),
        'compression_mode_id': mode_id,
        'compression_mode_label': str(compression_mode['label']),
        'reference_kernel_class_id': reference_kernel,
        'reference_interaction_label': str(reference.get('interaction_label', 'unknown')),
        'reference_coarse_label': str(reference.get('coarse_label', 'unknown')),
        'resolution': resolution,
        'initial_active_support': int(len(initial_support)),
        'initial_seed_count': int(len(seed_edges)),
        'steps_completed': int(trace[-1]['step_index']) if trace else 0,
        'final_survivor_count': int(len(current_support)),
        'final_survivor_fraction': float(len(current_support) / max(len(initial_support), 1)),
        'final_equivalence_class_count': final_equivalence_count,
        'final_largest_class_size': final_largest_class,
        'tau_stabilization_depth': int(tau_step),
        'delta_quiet_step': int(delta_quiet_step) if delta_quiet_step is not None else -1,
        'z2_arrow_step': int(arrow_step) if arrow_step is not None else -1,
        'persistent_step': int(persistent_step) if persistent_step is not None else -1,
        'obstruction_rank_final': int(current_state['cycle_rank']),
        'nonfunctoriality_initial': int(trace[0]['nonfunctoriality_before']) if trace else 0,
        'nonfunctoriality_final': int(trace[-1]['nonfunctoriality_after']) if trace else 0,
        'max_defect_density': float(max((float(item['defect_density']) for item in trace), default=0.0)),
        'final_defect_density': float(trace[-1]['defect_density']) if trace else 0.0,
        'pass_blocked_merge_total': int(pass_blocked_total),
        'defect_ledger_zero_before_pass': int(defect_zero_before_pass),
        'compression_class': '',
        'notes': str(run['notes']),
    }
    summary['compression_class'] = classify_compression(summary)

    plot_paths: list[Path] = []
    trace_path = WORK_PLOT_DIR / f"{run['run_id']}_compression_trace.png"
    preorder_trace_plot(trace_path, str(run['run_id']), trace)
    plot_paths.append(trace_path)

    overlay_path = WORK_PLOT_DIR / f"{run['run_id']}_survivor_overlay.png"
    survivor_overlay_plot(
        overlay_path,
        str(run['run_id']),
        np.asarray(data.midpoints, dtype=float),
        np.asarray(weights, dtype=float),
        initial_support,
        current_support,
        seed_edges,
        boundary_type,
    )
    plot_paths.append(overlay_path)

    result = {
        'run': run,
        'compression_mode': compression_mode,
        'summary': summary,
        'trace': trace,
        'seed_edges': [int(node) for node in seed_edges],
        'initial_active_support': [int(node) for node in initial_support],
        'final_support': [int(node) for node in current_support],
    }
    return result, plot_paths


def create_summary_plots(rows: list[dict[str, Any]]) -> list[Path]:
    plot_paths: list[Path] = []
    reps = sorted({str(row['representative_id']) for row in rows})
    modes = sorted({str(row['compression_mode_id']) for row in rows}, key=compression_mode_order)
    row_map = {(str(row['representative_id']), str(row['compression_mode_id'])): row for row in rows}

    class_scores = {
        'static_regime': 0,
        'compression_stable_without_arrow': 1,
        'pass_ordered_transient': 2,
        'pass_ordered_persistent': 3,
        'unresolved_preorder': 4,
    }
    matrix_path = WORK_PLOT_DIR / 'stage_c0_5_compression_class_matrix.png'
    matrix = np.zeros((len(reps), len(modes)), dtype=float)
    for i, rep in enumerate(reps):
        for j, mode in enumerate(modes):
            matrix[i, j] = float(class_scores.get(str(row_map[(rep, mode)]['compression_class']), 0))
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    im = ax.imshow(matrix, cmap='viridis', vmin=0, vmax=max(class_scores.values()))
    ax.set_xticks(range(len(modes)))
    ax.set_xticklabels(modes, rotation=20, ha='right')
    ax.set_yticks(range(len(reps)))
    ax.set_yticklabels(reps)
    ax.set_title('Compression class matrix')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(matrix_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(matrix_path)

    depth_path = WORK_PLOT_DIR / 'stage_c0_5_preorder_depth_panel.png'
    x = np.arange(len(reps))
    width = 0.22
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.6))
    for idx, mode in enumerate(modes):
        tau_values = [int(row_map[(rep, mode)]['tau_stabilization_depth']) for rep in reps]
        persistent_values = [max(int(row_map[(rep, mode)]['persistent_step']), 0) for rep in reps]
        axes[0].bar(x + (idx - 1) * width, tau_values, width=width, label=mode)
        axes[1].bar(x + (idx - 1) * width, persistent_values, width=width, label=mode)
    axes[0].set_title('Tau stabilization depth')
    axes[1].set_title('Persistent-step depth')
    for axis in axes:
        axis.set_xticks(x)
        axis.set_xticklabels(reps, rotation=20, ha='right')
        axis.grid(alpha=0.25, axis='y')
    axes[0].legend(fontsize=8)
    fig.savefig(depth_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(depth_path)

    defect_path = WORK_PLOT_DIR / 'stage_c0_5_defect_density_panel.png'
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.6))
    for idx, mode in enumerate(modes):
        max_delta = [float(row_map[(rep, mode)]['max_defect_density']) for rep in reps]
        final_delta = [float(row_map[(rep, mode)]['final_defect_density']) for rep in reps]
        arrow_values = [max(int(row_map[(rep, mode)]['z2_arrow_step']), 0) for rep in reps]
        axes[0].bar(x + (idx - 1) * width, max_delta, width=width, label=f'{mode} max')
        axes[1].bar(x + (idx - 1) * width, arrow_values, width=width, label=f'{mode} arrow')
        axes[0].plot(x + (idx - 1) * width, final_delta, marker='o', linewidth=1.0)
    axes[0].set_title('Defect density by mode')
    axes[1].set_title('Z2 arrow step by mode')
    for axis in axes:
        axis.set_xticks(x)
        axis.set_xticklabels(reps, rotation=20, ha='right')
        axis.grid(alpha=0.25, axis='y')
    axes[0].legend(fontsize=8)
    fig.savefig(defect_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(defect_path)

    return plot_paths


def write_note(path: Path, json_rel: str, csv_rel: str, rows: list[dict[str, Any]], runsheet: dict[str, Any]) -> None:
    lines = [
        '# Stage C0.5 Emergent Compression Preorder v1',
        '',
        f'Timestamped summary: `{json_rel}`',
        f'Timestamped run table: `{csv_rel}`',
        '',
        f'Architecture notice: {runsheet["architecture_notice"]}',
        f'Baseline notice: {runsheet["baseline_notice"]}',
        '',
        'Per-run summary:',
        '',
    ]
    for row in rows:
        lines.extend(
            [
                f"- `{row['compression_mode_label']}` / `{row['representative_id']}`",
                f"  - reference interaction label: `{row['reference_interaction_label']}`",
                f"  - compression class: `{row['compression_class']}`",
                f"  - tau stabilization depth: `{row['tau_stabilization_depth']}`",
                f"  - Z2 arrow step: `{row['z2_arrow_step']}`",
                f"  - persistent step: `{row['persistent_step']}`",
                f"  - final survivor count: `{row['final_survivor_count']}`",
                f"  - final defect density: `{row['final_defect_density']:.4f}`",
                f"  - final obstruction rank: `{row['obstruction_rank_final']}`",
            ]
        )
    lines.extend(
        [
            '',
            'Interpretation boundary:',
            '- this scan replaces the packet-evolution step with a combinatorial survivor-refinement preorder',
            '- the step index is a defect-ledger compression count, not an external physical time parameter',
            '- it does not claim geometry emergence, gravity, particles, or universal irreversibility',
        ]
    )
    path.write_text('\n'.join(lines), encoding='utf-8')


def main() -> None:
    args = parse_args()
    runsheet = load_runsheet(args.runsheet)
    defaults = load_stage10_defaults()
    base = dict(runsheet['base_seed_reference'])
    representatives = lookup_by_id(runsheet['representatives'], 'representative_id')
    compression_modes = lookup_by_id(runsheet['compression_modes'], 'compression_mode_id')
    c01_reference = load_c01_reference_rows()

    results: list[dict[str, Any]] = []
    csv_rows: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    for run in selected_runs(runsheet['runs'], args.run_ids):
        rep = representatives[str(run['representative_id'])]
        compression_mode = compression_modes[str(run['compression_mode_id'])]
        result, run_plots = simulate_run(run, rep, compression_mode, runsheet, defaults, base, c01_reference)
        results.append(result)
        csv_rows.append(dict(result['summary']))
        plot_paths.extend(run_plots)

    csv_rows.sort(key=lambda row: (str(row['representative_id']), compression_mode_order(str(row['compression_mode_id']))))
    summary_plots = create_summary_plots(csv_rows)
    plot_paths.extend(summary_plots)

    payload = {
        'stage': str(runsheet['stage']),
        'description': str(runsheet['description']),
        'architecture_notice': str(runsheet['architecture_notice']),
        'baseline_notice': str(runsheet['baseline_notice']),
        'results': results,
    }
    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug=str(runsheet.get('experiment_slug', 'stage_c0_5_emergent_compression_preorder')),
        result=payload,
        csv_rows=csv_rows,
        csv_fieldnames=CSV_FIELDS,
        plot_paths=plot_paths,
    )
    note_path = ATLAS_NOTES / str(runsheet['note_name'])
    write_note(
        note_path,
        str(json_path.relative_to(REPO_ROOT)),
        str(csv_path.relative_to(REPO_ROOT)),
        csv_rows,
        runsheet,
    )
    payload['plots'] = stamped_plots
    json_path.write_text(json.dumps(payload, indent=2), encoding='utf-8')
    print(json.dumps({'json': str(json_path), 'csv': str(csv_path), 'note': str(note_path)}, indent=2))


if __name__ == '__main__':
    main()
