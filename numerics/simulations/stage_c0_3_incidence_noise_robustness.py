#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, deque
from pathlib import Path
from typing import Any

import numpy as np
import scipy.sparse as sp

from stage10_common import (
    ATLAS_NOTES,
    REPO_ROOT,
    anisotropy_ratio,
    coherence_score,
    create_field_snapshot,
    displacement,
    load_stage10_defaults,
    packet_width,
    plt,
    save_atlas_payload,
    state_energy,
    suggested_dt,
    weighted_center,
)
from stage11_collective_wave_interaction import build_component_packet, classify_collective, make_leakage_fn
from stage12_coarse_field_structure import (
    ALPHA,
    SIGMAS,
    classify_coarse,
    component_mask,
    connected_components_periodic,
    edge_field_to_grid,
    gaussian_smooth_periodic,
    jaccard,
)
from stage9_common import estimate_lambda_max
from stage_c0_1_combinatorial_kernel_probe import (
    build_line_graph_cache,
    get_cached_setup,
    pairwise_distance_stats,
    pairwise_overlap_mean,
    pairwise_phase_alignment,
    path_length,
    render_grid_slice,
)

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_3_incidence_noise_runs.json'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage_c0_3_incidence_noise'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

PERTURBATION_CACHE: dict[tuple[int, str, str, int, int, int], dict[str, Any]] = {}
ANALYSIS_STRIDE = 4

CSV_FIELDS = [
    'run_id',
    'representative_id',
    'perturbation_id',
    'perturbation_label',
    'topology_class',
    'interaction_label',
    'coarse_label',
    'braid_like_exchange',
    'smeared_transfer',
    'unresolved',
    'flow_concentration_index',
    'exchange_coherence',
    'transport_span',
    'channel_count',
    'loop_count',
    'recurrence_indicator',
    'degree_spectrum_shift',
    'clustering_coefficient_shift',
    'mean_degree',
    'degree_std',
    'max_mean_overlap',
    'phase_lock_indicator',
    'constraint_max',
    'sector_leakage',
    'mutation_completed',
    'mutation_attempts',
    'touched_node_fraction',
    'edge_change_fraction',
    'delta_vs_baseline_concentration',
    'delta_vs_baseline_coherence',
    'delta_vs_baseline_transport_span',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.3 incidence-noise robustness scan.')
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
    return {item[key]: item for item in items}


def perturbation_order(perturbation_id: str) -> int:
    return {
        'baseline_no_noise': 0,
        'ultra_weak_rewiring': 1,
        'weak_rewiring': 2,
        'moderate_rewiring': 3,
        'moderate_deletion_insertion': 4,
        'clustered_rewiring_patch': 5,
        'random_localized_defect_pair': 6,
        'degree_biased_perturbation': 7,
        'connectivity_edge_stress_test': 8,
    }.get(perturbation_id, 99)


def copy_adjacency(adjacency_lists: list[list[int]]) -> list[set[int]]:
    return [set(items) for items in adjacency_lists]


def adjacency_components(adjacency_lists: list[list[int]]) -> list[list[int]]:
    visited = np.zeros(len(adjacency_lists), dtype=bool)
    comps: list[list[int]] = []
    for start in range(len(adjacency_lists)):
        if visited[start]:
            continue
        queue: deque[int] = deque([start])
        visited[start] = True
        comp: list[int] = []
        while queue:
            cur = queue.popleft()
            comp.append(cur)
            for nbr in adjacency_lists[cur]:
                if not visited[nbr]:
                    visited[nbr] = True
                    queue.append(nbr)
        comps.append(comp)
    return comps


def adjacency_edge_set(adjacency_lists: list[list[int]]) -> set[tuple[int, int]]:
    return {
        (idx, nbr)
        for idx, nbrs in enumerate(adjacency_lists)
        for nbr in nbrs
        if nbr > idx
    }


def compute_shells_from_adjacency(adjacency_lists: list[list[int]], max_shell: int) -> list[dict[int, int]]:
    shells: list[dict[int, int]] = []
    for start in range(len(adjacency_lists)):
        visited = {start: 0}
        queue: deque[int] = deque([start])
        local: dict[int, int] = {}
        while queue:
            cur = queue.popleft()
            depth = visited[cur]
            if depth >= max_shell:
                continue
            for nbr in adjacency_lists[cur]:
                if nbr in visited:
                    continue
                visited[nbr] = depth + 1
                local[nbr] = depth + 1
                queue.append(nbr)
        shells.append(local)
    return shells


def degree_stats(adjacency_lists: list[list[int]]) -> dict[str, Any]:
    degrees = np.asarray([len(items) for items in adjacency_lists], dtype=float)
    return {
        'degrees': degrees,
        'mean': float(np.mean(degrees)),
        'std': float(np.std(degrees)),
        'min': float(np.min(degrees)),
        'max': float(np.max(degrees)),
    }


def local_clustering_values(adjacency_lists: list[list[int]]) -> np.ndarray:
    adj_sets = [set(items) for items in adjacency_lists]
    values = np.zeros(len(adjacency_lists), dtype=float)
    for idx, nbrs in enumerate(adjacency_lists):
        k = len(nbrs)
        if k < 2:
            continue
        links = 0
        for offset, u in enumerate(nbrs[:-1]):
            u_neighbors = adj_sets[u]
            for v in nbrs[offset + 1:]:
                if v in u_neighbors:
                    links += 1
        values[idx] = 2.0 * links / float(k * (k - 1))
    return values


def bfs_patch(adjacency_lists: list[list[int]], seed: int, radius: int) -> list[int]:
    visited = {seed: 0}
    queue: deque[int] = deque([seed])
    nodes = [seed]
    while queue:
        cur = queue.popleft()
        depth = visited[cur]
        if depth >= radius:
            continue
        for nbr in adjacency_lists[cur]:
            if nbr in visited:
                continue
            visited[nbr] = depth + 1
            nodes.append(nbr)
            queue.append(nbr)
    return sorted(nodes)


def remove_edge(adj: list[set[int]], u: int, v: int) -> None:
    adj[u].remove(v)
    adj[v].remove(u)


def add_edge(adj: list[set[int]], u: int, v: int) -> None:
    adj[u].add(v)
    adj[v].add(u)


def endpoint_rewire_step(
    adj: list[set[int]],
    rng: np.random.Generator,
    min_degree_floor: int,
    source_pool: list[int] | None = None,
    target_pool: list[int] | None = None,
) -> tuple[int, int, int] | None:
    nodes = source_pool if source_pool else list(range(len(adj)))
    targets = target_pool if target_pool else list(range(len(adj)))
    for _ in range(400):
        u = int(rng.choice(nodes))
        if len(adj[u]) <= min_degree_floor:
            continue
        removable = [nbr for nbr in adj[u] if len(adj[nbr]) > min_degree_floor]
        if not removable:
            continue
        v = int(rng.choice(removable))
        candidates = [w for w in targets if w != u and w != v and w not in adj[u]]
        if not candidates:
            continue
        w = int(rng.choice(candidates))
        remove_edge(adj, u, v)
        add_edge(adj, u, w)
        return u, v, w
    return None


def delete_insert_step(
    adj: list[set[int]],
    rng: np.random.Generator,
    min_degree_floor: int,
    source_pool: list[int] | None = None,
    target_pool: list[int] | None = None,
) -> tuple[int, int, int, int] | None:
    nodes = source_pool if source_pool else list(range(len(adj)))
    targets = target_pool if target_pool else list(range(len(adj)))
    for _ in range(500):
        u = int(rng.choice(nodes))
        if len(adj[u]) <= min_degree_floor:
            continue
        removable = [nbr for nbr in adj[u] if len(adj[nbr]) > min_degree_floor]
        if not removable:
            continue
        v = int(rng.choice(removable))

        a = int(rng.choice(targets))
        candidate_bs = [b for b in targets if b != a and b not in adj[a]]
        if not candidate_bs:
            continue
        b = int(rng.choice(candidate_bs))
        if {a, b} == {u, v}:
            continue

        remove_edge(adj, u, v)
        add_edge(adj, a, b)
        return u, v, a, b
    return None


def low_support_edge(adjacency_lists: list[list[int]], rng: np.random.Generator) -> tuple[int, int] | None:
    best: list[tuple[int, int]] = []
    best_score: int | None = None
    adj_sets = [set(items) for items in adjacency_lists]
    for u, nbrs in enumerate(adjacency_lists):
        for v in nbrs:
            if v <= u:
                continue
            score = len(adj_sets[u].intersection(adj_sets[v]))
            if best_score is None or score < best_score:
                best_score = score
                best = [(u, v)]
            elif score == best_score:
                best.append((u, v))
    if not best:
        return None
    return best[int(rng.integers(len(best)))]


def touch_nodes(touch_counts: np.ndarray, *nodes: int) -> None:
    for node in nodes:
        touch_counts[int(node)] += 1


def finalize_perturbation(
    adjacency_sets: list[set[int]],
    baseline_edges: set[tuple[int, int]],
    touch_counts: np.ndarray,
    mutations: list[dict[str, Any]],
    attempts: int,
    seed_nodes: list[int],
) -> dict[str, Any]:
    adjacency_lists = [sorted(items) for items in adjacency_sets]
    changed_edges = adjacency_edge_set(adjacency_lists).symmetric_difference(baseline_edges)
    return {
        'adjacency_lists': adjacency_lists,
        'touch_counts': touch_counts,
        'mutations': mutations,
        'mutation_completed': len(mutations),
        'mutation_attempts': attempts,
        'touched_node_fraction': float(np.mean(touch_counts > 0)),
        'edge_change_fraction': float(len(changed_edges) / max(len(baseline_edges), 1)),
        'seed_nodes': [int(value) for value in seed_nodes],
        'connectivity_ok': int(len(adjacency_components(adjacency_lists)) == 1),
    }


def build_incidence_perturbation(
    data: Any,
    perturbation_id: str,
    params: dict[str, Any],
) -> dict[str, Any]:
    max_shell = int(params.get('max_shell', 3))
    random_seed = int(params.get('random_seed', 3303))
    min_degree_floor = int(params.get('min_degree_floor', 3))
    patch_radius = int(params.get('patch_radius', 2))
    localized_radius = int(params.get('localized_radius', 1))
    cache_key = (
        int(data.n_side),
        str(data.boundary_type),
        str(perturbation_id),
        max_shell,
        random_seed,
        min_degree_floor,
    )
    if cache_key in PERTURBATION_CACHE:
        return PERTURBATION_CACHE[cache_key]

    baseline = build_line_graph_cache(data, max_shell=max_shell)
    baseline_adj = baseline['adjacency_lists']
    baseline_edges = adjacency_edge_set(baseline_adj)
    node_count = len(baseline_adj)

    if perturbation_id == 'baseline_no_noise':
        payload = finalize_perturbation(
            adjacency_sets=copy_adjacency(baseline_adj),
            baseline_edges=baseline_edges,
            touch_counts=np.zeros(node_count, dtype=int),
            mutations=[],
            attempts=0,
            seed_nodes=[],
        )
        PERTURBATION_CACHE[cache_key] = payload
        return payload

    for restart in range(10):
        rng = np.random.default_rng(random_seed + 1000 * perturbation_order(perturbation_id) + restart)
        adj = copy_adjacency(baseline_adj)
        touch_counts = np.zeros(node_count, dtype=int)
        mutations: list[dict[str, Any]] = []
        attempts = 0
        seed_nodes: list[int] = []

        if perturbation_id == 'ultra_weak_rewiring':
            op_count = int(params.get('ultra_weak_rewire_ops', 16))
            for _ in range(op_count):
                attempts += 1
                moved = endpoint_rewire_step(adj, rng, min_degree_floor=min_degree_floor)
                if moved is None:
                    continue
                u, v, w = moved
                touch_nodes(touch_counts, u, v, w)
                mutations.append({'kind': 'endpoint_rewire', 'u': u, 'v': v, 'w': w})
        elif perturbation_id == 'weak_rewiring':
            op_count = int(params.get('weak_rewire_ops', 48))
            for _ in range(op_count):
                attempts += 1
                moved = endpoint_rewire_step(adj, rng, min_degree_floor=min_degree_floor)
                if moved is None:
                    continue
                u, v, w = moved
                touch_nodes(touch_counts, u, v, w)
                mutations.append({'kind': 'endpoint_rewire', 'u': u, 'v': v, 'w': w})
        elif perturbation_id == 'moderate_rewiring':
            op_count = int(params.get('moderate_rewire_ops', 144))
            for _ in range(op_count):
                attempts += 1
                moved = endpoint_rewire_step(adj, rng, min_degree_floor=min_degree_floor)
                if moved is None:
                    continue
                u, v, w = moved
                touch_nodes(touch_counts, u, v, w)
                mutations.append({'kind': 'endpoint_rewire', 'u': u, 'v': v, 'w': w})
        elif perturbation_id == 'moderate_deletion_insertion':
            op_count = int(params.get('moderate_delete_insert_ops', 96))
            for _ in range(op_count):
                attempts += 1
                moved = delete_insert_step(adj, rng, min_degree_floor=min_degree_floor)
                if moved is None:
                    continue
                u, v, a, b = moved
                touch_nodes(touch_counts, u, v, a, b)
                mutations.append({'kind': 'delete_insert', 'u': u, 'v': v, 'a': a, 'b': b})
        elif perturbation_id == 'clustered_rewiring_patch':
            seed = int(rng.integers(node_count))
            patch_nodes = bfs_patch(baseline_adj, seed, patch_radius)
            seed_nodes = [seed]
            op_count = int(params.get('cluster_patch_ops', 96))
            for _ in range(op_count):
                attempts += 1
                moved = endpoint_rewire_step(
                    adj,
                    rng,
                    min_degree_floor=min_degree_floor,
                    source_pool=patch_nodes,
                    target_pool=patch_nodes,
                )
                if moved is None:
                    continue
                u, v, w = moved
                touch_nodes(touch_counts, u, v, w)
                mutations.append({'kind': 'patch_rewire', 'u': u, 'v': v, 'w': w})
        elif perturbation_id == 'random_localized_defect_pair':
            seed_a = int(rng.integers(node_count))
            far_candidates = [idx for idx in range(node_count) if idx != seed_a and idx not in baseline_adj[seed_a]]
            seed_b = int(rng.choice(far_candidates)) if far_candidates else int((seed_a + node_count // 2) % node_count)
            patch_a = bfs_patch(baseline_adj, seed_a, localized_radius)
            patch_b = bfs_patch(baseline_adj, seed_b, localized_radius)
            local_nodes = sorted(set(patch_a + patch_b))
            seed_nodes = [seed_a, seed_b]
            op_count = int(params.get('localized_defect_ops', 72))
            for _ in range(op_count):
                attempts += 1
                moved = delete_insert_step(
                    adj,
                    rng,
                    min_degree_floor=min_degree_floor,
                    source_pool=local_nodes,
                    target_pool=local_nodes,
                )
                if moved is None:
                    continue
                u, v, a, b = moved
                touch_nodes(touch_counts, u, v, a, b)
                mutations.append({'kind': 'localized_defect', 'u': u, 'v': v, 'a': a, 'b': b})
        elif perturbation_id == 'degree_biased_perturbation':
            op_count = int(params.get('degree_biased_ops', 120))
            for _ in range(op_count):
                attempts += 1
                degrees = np.asarray([len(items) for items in adj], dtype=float)
                source_pool = [idx for idx, value in enumerate(degrees) if value <= np.percentile(degrees, 40.0)]
                target_pool = [idx for idx, value in enumerate(degrees) if value >= np.percentile(degrees, 75.0)]
                moved = endpoint_rewire_step(
                    adj,
                    rng,
                    min_degree_floor=min_degree_floor,
                    source_pool=source_pool if source_pool else None,
                    target_pool=target_pool if target_pool else None,
                )
                if moved is None:
                    continue
                u, v, w = moved
                touch_nodes(touch_counts, u, v, w)
                mutations.append({'kind': 'degree_biased', 'u': u, 'v': v, 'w': w})
        elif perturbation_id == 'connectivity_edge_stress_test':
            op_count = int(params.get('connectivity_stress_ops', 96))
            for _ in range(op_count):
                attempts += 1
                temp_lists = [sorted(items) for items in adj]
                picked = low_support_edge(temp_lists, rng)
                if picked is None:
                    continue
                u, v = picked
                if len(adj[u]) <= min_degree_floor or len(adj[v]) <= min_degree_floor:
                    continue
                candidates = [w for w in range(node_count) if w not in adj[u] and w != u and w != v]
                if not candidates:
                    continue
                w = int(rng.choice(candidates))
                remove_edge(adj, u, v)
                add_edge(adj, u, w)
                touch_nodes(touch_counts, u, v, w)
                mutations.append({'kind': 'connectivity_stress', 'u': u, 'v': v, 'w': w})
        else:
            raise ValueError(f'unsupported perturbation id: {perturbation_id}')

        payload = finalize_perturbation(adj, baseline_edges, touch_counts, mutations, attempts, seed_nodes)
        if payload['connectivity_ok']:
            PERTURBATION_CACHE[cache_key] = payload
            return payload

    raise RuntimeError(f'failed to build connected perturbation for {perturbation_id}')


def build_operator_from_adjacency(
    data: Any,
    adjacency_lists: list[list[int]],
    params: dict[str, Any],
    target_lambda: float,
) -> dict[str, Any]:
    max_shell = int(params.get('max_shell', 3))
    graph_shell_weights = {
        int(key): float(value)
        for key, value in dict(params.get('graph_shell_weights', {'1': 1.0, '2': 0.5, '3': 0.25})).items()
    }
    shells = compute_shells_from_adjacency(adjacency_lists, max_shell=max_shell)

    rows: list[int] = []
    cols: list[int] = []
    values: list[float] = []
    for edge_idx, shell_map in enumerate(shells):
        for nbr_idx, shell in shell_map.items():
            if nbr_idx <= edge_idx:
                continue
            weight = float(graph_shell_weights.get(shell, 0.0))
            if weight <= 0.0:
                continue
            rows.extend([edge_idx, nbr_idx])
            cols.extend([nbr_idx, edge_idx])
            values.extend([weight, weight])

    size = len(adjacency_lists)
    weight_matrix = sp.coo_matrix((values, (rows, cols)), shape=(size, size), dtype=float).tocsr()
    weighted_degree_vec = np.asarray(weight_matrix.sum(axis=1)).ravel()
    operator = (sp.diags(weighted_degree_vec) - weight_matrix).tocsr()
    raw_lambda = estimate_lambda_max(operator)
    if target_lambda > 0.0 and raw_lambda > 0.0:
        scale_factor = float(target_lambda / raw_lambda)
        operator = (scale_factor * operator).tocsr()
        lambda_max = estimate_lambda_max(operator)
    else:
        scale_factor = 1.0
        lambda_max = raw_lambda
    stats = degree_stats(adjacency_lists)
    clustering_values = local_clustering_values(adjacency_lists)
    return {
        'adjacency_lists': adjacency_lists,
        'shells': shells,
        'operator': operator,
        'weight_matrix': weight_matrix,
        'graph_degrees': stats['degrees'],
        'degree_mean': float(stats['mean']),
        'degree_std': float(stats['std']),
        'degree_min': float(stats['min']),
        'degree_max': float(stats['max']),
        'clustering_values': clustering_values,
        'mean_clustering': float(np.mean(clustering_values)),
        'raw_lambda_max': float(raw_lambda),
        'lambda_max': float(lambda_max),
        'scale_factor': float(scale_factor),
        'support_entries': int(weight_matrix.nnz),
        'max_shell': max_shell,
    }


def recurrence_indicator(series: np.ndarray, close_threshold: float) -> float:
    if len(series) < 6:
        return 0.0
    tail = series[len(series) // 3 :]
    if tail.size < 4:
        return 0.0
    deriv = np.diff(tail)
    signs = np.sign(deriv)
    signs = signs[signs != 0.0]
    if signs.size < 2:
        return 0.0
    turns = int(np.sum(signs[1:] * signs[:-1] < 0.0))
    closeness = float(np.mean(tail <= close_threshold))
    return float((turns / max(1, signs.size - 1)) * closeness)


def active_graph_metrics(avg_weights: np.ndarray, adjacency_lists: list[list[int]], sigma_scale: float) -> tuple[int, int]:
    if avg_weights.size == 0:
        return 0, 0
    threshold = float(np.mean(avg_weights) + sigma_scale * np.std(avg_weights))
    threshold = max(threshold, 0.2 * float(np.max(avg_weights)))
    active = avg_weights >= threshold
    active_nodes = [idx for idx, flag in enumerate(active) if flag]
    if not active_nodes:
        return 0, 0
    active_set = set(active_nodes)
    visited: set[int] = set()
    components = 0
    edge_count = 0
    for node in active_nodes:
        edge_count += sum(1 for nbr in adjacency_lists[node] if nbr in active_set)
        if node in visited:
            continue
        components += 1
        queue: deque[int] = deque([node])
        visited.add(node)
        while queue:
            cur = queue.popleft()
            for nbr in adjacency_lists[cur]:
                if nbr in active_set and nbr not in visited:
                    visited.add(nbr)
                    queue.append(nbr)
    edge_count //= 2
    loop_count = max(edge_count - len(active_nodes) + components, 0)
    return components, loop_count


def exchange_coherence_metric(overlap_series: list[float], phase_alignment_series: list[float]) -> float:
    if not overlap_series and not phase_alignment_series:
        return 0.0
    overlap_term = float(np.max(overlap_series)) if overlap_series else 0.0
    phase_term = float(np.mean(phase_alignment_series)) if phase_alignment_series else 0.0
    return float(0.5 * overlap_term + 0.5 * phase_term)


def perturbation_overlay_plot(
    path: Path,
    run_id: str,
    midpoints: np.ndarray,
    touch_counts: np.ndarray,
    avg_weights: np.ndarray,
    boundary_type: str,
) -> None:
    coords = np.asarray(midpoints, dtype=float)
    slice_value = 0.5 if boundary_type == 'periodic' else 0.5
    distances = np.abs(coords[:, 2] - slice_value)
    threshold = np.min(distances) + 1.0e-9
    mask = distances <= threshold
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.6))
    sc0 = axes[0].scatter(coords[mask, 0], coords[mask, 1], c=touch_counts[mask], s=18, cmap='magma')
    axes[0].set_title('Incidence perturbation footprint')
    axes[0].set_xlabel('x')
    axes[0].set_ylabel('y')
    plt.colorbar(sc0, ax=axes[0], fraction=0.045, pad=0.04)
    sc1 = axes[1].scatter(coords[mask, 0], coords[mask, 1], c=avg_weights[mask], s=18, cmap='viridis')
    axes[1].set_title('Average flow intensity')
    axes[1].set_xlabel('x')
    axes[1].set_ylabel('y')
    plt.colorbar(sc1, ax=axes[1], fraction=0.045, pad=0.04)
    fig.suptitle(run_id)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def simulate_run(
    run: dict[str, Any],
    rep: dict[str, Any],
    perturbation: dict[str, Any],
    runsheet: dict[str, Any],
    defaults: dict[str, Any],
    base: dict[str, Any],
) -> tuple[dict[str, Any], list[Path]]:
    resolution = int(run['resolution'])
    boundary_type = str(base['boundary_type'])
    data, projector = get_cached_setup(resolution, defaults, boundary_type)
    target_lambda = estimate_lambda_max(data.L1) if bool(runsheet.get('parameter_policy', {}).get('match_lambda_to_frozen_L1', True)) else 0.0

    perturbation_info = build_incidence_perturbation(data, str(run['perturbation_id']), runsheet.get('parameter_policy', {}))
    operator_info = build_operator_from_adjacency(
        data=data,
        adjacency_lists=perturbation_info['adjacency_lists'],
        params=runsheet.get('parameter_policy', {}),
        target_lambda=target_lambda,
    )
    baseline_operator_info = build_operator_from_adjacency(
        data=data,
        adjacency_lists=build_line_graph_cache(data, max_shell=int(runsheet.get('parameter_policy', {}).get('max_shell', 3)))['adjacency_lists'],
        params=runsheet.get('parameter_policy', {}),
        target_lambda=target_lambda,
    )

    operator = operator_info['operator']
    dt = suggested_dt(operator, float(defaults['dt_safety']))
    steps = max(2, int(np.ceil(float(defaults['t_final']) / dt)))
    leakage_fn = make_leakage_fn(projector)
    kick_axis = int(defaults['kick_axis'])

    packet_qs: list[np.ndarray] = []
    packet_vs: list[np.ndarray] = []
    for center, amp, phase, kick_sign in zip(rep['packet_centers'], rep['packet_amplitudes'], rep['phase_offsets_rad'], rep['kick_signs']):
        q0, v0 = build_component_packet(
            center=np.asarray(center, dtype=float),
            amplitude=float(amp),
            phase_offset=float(phase),
            kick_sign=float(kick_sign),
            bandwidth=float(base['bandwidth']),
            central_k=float(base['central_k']),
            phase_pattern=str(base['phase_pattern']),
            kick_axis=kick_axis,
            boundary_type=boundary_type,
            midpoints=data.midpoints,
            edge_axes=data.edge_axes,
            projector=projector,
        )
        packet_qs.append(q0)
        packet_vs.append(v0)

    times: list[float] = []
    total_centers: list[list[float]] = []
    widths: list[float] = []
    energies: list[float] = []
    anisotropies: list[float] = []
    sector_leakages: list[float] = []
    constraint_norms: list[float] = []
    coherences: list[float] = []
    pair_distance_series: list[float] = []
    overlap_series: list[float] = []
    phase_alignment_series: list[float] = []
    occupancy_accumulator = np.zeros(len(data.edges), dtype=float)
    occupancy_samples = 0
    coarse: dict[float, dict[str, Any]] = {
        sigma: {
            'dominant_area_fractions': [],
            'basin_counts': [],
            'stable_flags': [],
            'dominant_masks': [],
            'envelope_variances': [],
            'last_env': None,
        }
        for sigma in SIGMAS
    }

    def record(t: float) -> None:
        nonlocal occupancy_samples
        total_q = np.sum(packet_qs, axis=0)
        total_v = np.sum(packet_vs, axis=0)
        total_weights = np.abs(total_q) ** 2
        occupancy_accumulator[:] += total_weights
        occupancy_samples += 1
        total_center = weighted_center(data.midpoints, total_weights, boundary_type)
        times.append(float(t))
        total_centers.append(total_center.tolist())
        widths.append(packet_width(data.midpoints, total_weights, total_center, boundary_type))
        energies.append(state_energy(operator, total_q, total_v))
        anisotropies.append(anisotropy_ratio(data.midpoints, total_weights, total_center, boundary_type))
        sector_leakages.append(float(leakage_fn(total_q)))
        constraint_norms.append(float(np.linalg.norm(data.d0.T @ total_q)))
        coherences.append(coherence_score(total_weights))

        component_centers: list[np.ndarray] = []
        for comp_q in packet_qs:
            comp_weights = np.abs(comp_q) ** 2
            component_centers.append(weighted_center(data.midpoints, comp_weights, boundary_type))
        mean_pair, _min_pair = pairwise_distance_stats(component_centers, boundary_type)
        pair_distance_series.append(mean_pair)
        overlap_series.append(pairwise_overlap_mean(packet_qs))
        phase_alignment_series.append(pairwise_phase_alignment(packet_qs))

        point_grid = edge_field_to_grid(data.midpoints, np.abs(total_q), resolution)
        for sigma in SIGMAS:
            env = np.maximum(gaussian_smooth_periodic(point_grid, sigma), 0.0)
            coarse[sigma]['last_env'] = env
            max_env = float(np.max(env))
            norm_env = env / max(max_env, 1.0e-12)
            coarse[sigma]['envelope_variances'].append(float(np.var(norm_env)))
            mask = env > (ALPHA * max_env) if max_env > 0.0 else np.zeros_like(env, dtype=bool)
            comps = connected_components_periodic(mask)
            coarse[sigma]['basin_counts'].append(float(len(comps)))
            if comps:
                dominant = max(comps, key=len)
                dominant_mask = component_mask(dominant, env.shape[0])
                dom_frac = float(len(dominant) / env.size)
            else:
                dominant_mask = None
                dom_frac = 0.0
            prev_mask = coarse[sigma]['dominant_masks'][-1] if coarse[sigma]['dominant_masks'] else None
            if dominant_mask is None:
                stable = 0.0
            elif prev_mask is None:
                stable = 1.0
            else:
                stable = 1.0 if jaccard(prev_mask, dominant_mask) >= 0.5 else 0.0
            coarse[sigma]['dominant_area_fractions'].append(dom_frac)
            coarse[sigma]['stable_flags'].append(stable)
            coarse[sigma]['dominant_masks'].append(dominant_mask)

    record(0.0)
    for step_idx in range(steps):
        next_qs: list[np.ndarray] = []
        next_vs: list[np.ndarray] = []
        for q, v in zip(packet_qs, packet_vs):
            acc = np.asarray(projector(-(operator @ q)), dtype=float)
            v_half = v + 0.5 * dt * acc
            q_new = np.asarray(projector(q + dt * v_half), dtype=float)
            acc_new = np.asarray(projector(-(operator @ q_new)), dtype=float)
            v_new = np.asarray(projector(v_half + 0.5 * dt * acc_new), dtype=float)
            next_qs.append(q_new)
            next_vs.append(v_new)
        packet_qs = next_qs
        packet_vs = next_vs
        if ((step_idx + 1) % ANALYSIS_STRIDE == 0) or (step_idx == steps - 1):
            record((step_idx + 1) * dt)

    sample_dt = float(np.median(np.diff(times))) if len(times) >= 2 else dt
    centers_arr = np.asarray(total_centers, dtype=float)
    avg_weights = occupancy_accumulator / max(occupancy_samples, 1)
    avg_weights = avg_weights / max(float(np.sum(avg_weights)), 1.0e-12)
    flow_concentration = float(len(avg_weights) * coherence_score(avg_weights))
    exchange_coherence = exchange_coherence_metric(overlap_series, phase_alignment_series)
    transport_span = path_length(centers_arr, boundary_type)
    recurrence = recurrence_indicator(np.asarray(pair_distance_series, dtype=float), close_threshold=0.18)
    active_sigma = float(runsheet.get('parameter_policy', {}).get('active_threshold_sigma', 1.0))
    channel_count, loop_count = active_graph_metrics(avg_weights, operator_info['adjacency_lists'], active_sigma)

    interaction_summary = {
        'packet_count': int(rep['packet_count']),
        'center_shift': float(np.linalg.norm(displacement(centers_arr[-1][None, :], centers_arr[0], boundary_type)[0])),
        'width_ratio': float(widths[-1] / max(widths[0], 1.0e-12)),
        'max_mean_overlap': float(np.max(overlap_series)) if overlap_series else 0.0,
        'phase_lock_indicator': float(np.mean(phase_alignment_series)) if phase_alignment_series else 0.0,
        'composite_lifetime': float(np.sum(np.asarray(overlap_series) >= 0.35) * sample_dt),
        'oscillation_count': 0,
        'initial_mean_pair_distance': float(pair_distance_series[0]) if pair_distance_series else 0.0,
        'min_mean_pair_distance': float(np.min(pair_distance_series)) if pair_distance_series else 0.0,
        'exchange_or_merger_flag': 'none',
        't_final': float(times[-1]),
    }
    interaction_label = classify_collective(interaction_summary)

    coarse_summary: dict[str, float] = {'t_final': float(times[-1])}
    for sigma in SIGMAS:
        coarse_summary[f'mean_basin_count_sigma{int(sigma)}'] = float(np.mean(coarse[sigma]['basin_counts'])) if coarse[sigma]['basin_counts'] else 0.0
        coarse_summary[f'mean_dominant_area_fraction_sigma{int(sigma)}'] = float(np.mean(coarse[sigma]['dominant_area_fractions'])) if coarse[sigma]['dominant_area_fractions'] else 0.0
        coarse_summary[f'coarse_persistence_sigma{int(sigma)}'] = float(np.mean(coarse[sigma]['stable_flags'])) if coarse[sigma]['stable_flags'] else 0.0
        coarse_summary[f'mean_envelope_variance_sigma{int(sigma)}'] = float(np.mean(coarse[sigma]['envelope_variances'])) if coarse[sigma]['envelope_variances'] else 0.0
        coarse_summary[f'basin_lifetime_sigma{int(sigma)}'] = float(np.sum(np.asarray(coarse[sigma]['dominant_area_fractions']) >= 0.05) * sample_dt)
    coarse_label = classify_coarse(coarse_summary)

    baseline_degrees = np.asarray(baseline_operator_info['graph_degrees'], dtype=float)
    run_degrees = np.asarray(operator_info['graph_degrees'], dtype=float)
    degree_spectrum_shift = float(np.linalg.norm(np.sort(run_degrees) - np.sort(baseline_degrees)) / math.sqrt(max(len(run_degrees), 1)))
    clustering_shift = float(operator_info['mean_clustering'] - baseline_operator_info['mean_clustering'])

    row = {
        'run_id': run['run_id'],
        'representative_id': run['representative_id'],
        'perturbation_id': run['perturbation_id'],
        'perturbation_label': perturbation['label'],
        'topology_class': 'unclassified',
        'interaction_label': interaction_label,
        'coarse_label': coarse_label,
        'braid_like_exchange': 0,
        'smeared_transfer': 0,
        'unresolved': 0,
        'flow_concentration_index': flow_concentration,
        'exchange_coherence': exchange_coherence,
        'transport_span': transport_span,
        'channel_count': int(channel_count),
        'loop_count': int(loop_count),
        'recurrence_indicator': recurrence,
        'degree_spectrum_shift': degree_spectrum_shift,
        'clustering_coefficient_shift': clustering_shift,
        'mean_degree': float(operator_info['degree_mean']),
        'degree_std': float(operator_info['degree_std']),
        'max_mean_overlap': float(np.max(overlap_series)) if overlap_series else 0.0,
        'phase_lock_indicator': float(np.mean(phase_alignment_series)) if phase_alignment_series else 0.0,
        'constraint_max': float(np.max(constraint_norms)) if constraint_norms else 0.0,
        'sector_leakage': float(np.max(sector_leakages)) if sector_leakages else 0.0,
        'mutation_completed': int(perturbation_info['mutation_completed']),
        'mutation_attempts': int(perturbation_info['mutation_attempts']),
        'touched_node_fraction': float(perturbation_info['touched_node_fraction']),
        'edge_change_fraction': float(perturbation_info['edge_change_fraction']),
        'delta_vs_baseline_concentration': 0.0,
        'delta_vs_baseline_coherence': 0.0,
        'delta_vs_baseline_transport_span': 0.0,
        'notes': run['notes'],
    }

    plot_paths: list[Path] = []

    trace_path = WORK_PLOT_DIR / f"{run['run_id']}_coherence_flow_trace.png"
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))
    axes[0].plot(times, overlap_series, color='tab:red', label='overlap')
    axes[0].plot(times, phase_alignment_series, color='tab:purple', label='phase')
    axes[0].set_title('Exchange coherence trace')
    axes[0].set_xlabel('time')
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=8)
    center_shift_trace = [
        float(np.linalg.norm(displacement(np.asarray([center], dtype=float), centers_arr[0], boundary_type)[0]))
        for center in centers_arr
    ]
    axes[1].plot(times, center_shift_trace, color='tab:orange', label='center shift')
    axes[1].set_title('Transport trace')
    axes[1].set_xlabel('time')
    axes[1].grid(alpha=0.25)
    fig.suptitle(run['run_id'])
    fig.savefig(trace_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(trace_path)

    overlay_path = WORK_PLOT_DIR / f"{run['run_id']}_perturbation_flow_overlay.png"
    perturbation_overlay_plot(
        path=overlay_path,
        run_id=run['run_id'],
        midpoints=data.midpoints,
        touch_counts=np.asarray(perturbation_info['touch_counts'], dtype=float),
        avg_weights=avg_weights,
        boundary_type=boundary_type,
    )
    plot_paths.append(overlay_path)

    field_path = WORK_PLOT_DIR / f"{run['run_id']}_field_snapshot.png"
    create_field_snapshot(field_path, run['run_id'], data.midpoints, np.sum(packet_qs, axis=0), boundary_type, int(defaults['slice_axis']))
    plot_paths.append(field_path)

    env = coarse[4.0]['last_env']
    if env is not None:
        env_path = WORK_PLOT_DIR / f"{run['run_id']}_coarse_envelope.png"
        fig, ax = plt.subplots(figsize=(5.4, 4.6))
        render_grid_slice(ax, env, f"Envelope sigma=4: {run['run_id']}")
        fig.savefig(env_path, dpi=180, bbox_inches='tight')
        plt.close(fig)
        plot_paths.append(env_path)

    result = {
        'run': run,
        'perturbation': perturbation,
        'perturbation_summary': {
            'mutation_completed': int(perturbation_info['mutation_completed']),
            'mutation_attempts': int(perturbation_info['mutation_attempts']),
            'touched_node_fraction': float(perturbation_info['touched_node_fraction']),
            'edge_change_fraction': float(perturbation_info['edge_change_fraction']),
            'seed_nodes': [int(value) for value in perturbation_info['seed_nodes']],
            'connectivity_ok': int(perturbation_info['connectivity_ok']),
            'mutation_log_excerpt': perturbation_info['mutations'][:20],
        },
        'operator_summary': {
            'lambda_max': float(operator_info['lambda_max']),
            'raw_lambda_max': float(operator_info['raw_lambda_max']),
            'scale_factor': float(operator_info['scale_factor']),
            'support_entries': int(operator_info['support_entries']),
            'degree_mean': float(operator_info['degree_mean']),
            'degree_std': float(operator_info['degree_std']),
            'mean_clustering': float(operator_info['mean_clustering']),
        },
        'summary': row,
        'times': times,
        'pair_distance_series': pair_distance_series,
        'overlap_series': overlap_series,
        'phase_alignment_series': phase_alignment_series,
        'coherences': coherences,
    }
    return result, plot_paths


def classify_against_baseline(rows: list[dict[str, Any]], runsheet: dict[str, Any]) -> list[dict[str, Any]]:
    baseline = next((row for row in rows if str(row['perturbation_id']) == 'baseline_no_noise'), None)
    if baseline is None:
        return rows
    flow_tol = float(runsheet.get('parameter_policy', {}).get('baseline_flow_tolerance', 5.0))
    coherence_tol = float(runsheet.get('parameter_policy', {}).get('baseline_coherence_tolerance', 0.04))
    transport_floor = float(runsheet.get('parameter_policy', {}).get('braid_transport_floor', 0.08))
    baseline_flow = float(baseline['flow_concentration_index'])
    baseline_coherence = float(baseline['exchange_coherence'])
    baseline_span = float(baseline['transport_span'])

    for row in rows:
        row['delta_vs_baseline_concentration'] = float(row['flow_concentration_index']) - baseline_flow
        row['delta_vs_baseline_coherence'] = float(row['exchange_coherence']) - baseline_coherence
        row['delta_vs_baseline_transport_span'] = float(row['transport_span']) - baseline_span

        braid_condition = (
            float(row['flow_concentration_index']) >= baseline_flow - flow_tol
            and float(row['exchange_coherence']) >= max(baseline_coherence - coherence_tol, 0.12)
            and float(row['transport_span']) >= max(baseline_span + transport_floor, transport_floor)
            and int(row['channel_count']) <= 2
            and int(row['loop_count']) == 0
            and float(row['recurrence_indicator']) <= 0.1
        )
        smeared_condition = (
            float(row['flow_concentration_index']) >= max(0.6 * baseline_flow, 20.0)
            and float(row['exchange_coherence']) >= max(0.5 * baseline_coherence, 0.08)
            and float(row['recurrence_indicator']) <= 0.35
        )

        if braid_condition:
            topology = 'braid_like_exchange'
        elif smeared_condition:
            topology = 'smeared_transfer'
        else:
            topology = 'unresolved'

        row['topology_class'] = topology
        row['braid_like_exchange'] = int(topology == 'braid_like_exchange')
        row['smeared_transfer'] = int(topology == 'smeared_transfer')
        row['unresolved'] = int(topology == 'unresolved')
    return rows


def create_summary_plots(rows: list[dict[str, Any]]) -> list[Path]:
    plot_paths: list[Path] = []
    ordered = sorted(rows, key=lambda item: perturbation_order(str(item['perturbation_id'])))
    labels = [str(row['perturbation_label']) for row in ordered]
    x = np.arange(len(ordered))
    topology_scores = {
        'unresolved': 0,
        'smeared_transfer': 1,
        'braid_like_exchange': 2,
    }

    topo_path = WORK_PLOT_DIR / 'stage_c0_3_topology_classification_summary.png'
    data = np.asarray([[topology_scores.get(str(row['topology_class']), 0) for row in ordered]], dtype=float)
    fig, ax = plt.subplots(figsize=(11.2, 2.8))
    im = ax.imshow(data, cmap='viridis', vmin=0, vmax=max(topology_scores.values()))
    ax.set_yticks([0])
    ax.set_yticklabels(['corridor'])
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=25, ha='right')
    ax.set_title('Topology classification summary')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(topo_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(topo_path)

    flow_path = WORK_PLOT_DIR / 'stage_c0_3_flow_concentration_panel.png'
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.6))
    axes[0].bar(x, [float(row['flow_concentration_index']) for row in ordered], color='tab:blue')
    axes[0].set_title('Flow concentration')
    axes[1].bar(x, [float(row['exchange_coherence']) for row in ordered], color='tab:orange')
    axes[1].set_title('Exchange coherence')
    for axis in axes:
        axis.set_xticks(x)
        axis.set_xticklabels(labels, rotation=25, ha='right')
        axis.grid(alpha=0.25, axis='y')
    fig.savefig(flow_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(flow_path)

    drift_path = WORK_PLOT_DIR / 'stage_c0_3_degree_spectrum_drift_panel.png'
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.6))
    axes[0].bar(x, [float(row['degree_spectrum_shift']) for row in ordered], color='tab:green')
    axes[0].set_title('Degree spectrum shift')
    axes[1].bar(x, [float(row['clustering_coefficient_shift']) for row in ordered], color='tab:red')
    axes[1].set_title('Clustering coefficient shift')
    for axis in axes:
        axis.set_xticks(x)
        axis.set_xticklabels(labels, rotation=25, ha='right')
        axis.grid(alpha=0.25, axis='y')
    fig.savefig(drift_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(drift_path)

    return plot_paths


def write_note(
    path: Path,
    json_rel: str,
    csv_rel: str,
    rows: list[dict[str, Any]],
    runsheet: dict[str, Any],
) -> None:
    lines = [
        '# Stage C0.3 Topology-Exchange Robustness Under Incidence Noise v1',
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
    for row in sorted(rows, key=lambda item: perturbation_order(str(item['perturbation_id']))):
        lines.extend([
            f"- `{row['perturbation_label']}`",
            f"  - topology class: `{row['topology_class']}`",
            f"  - flow concentration: `{row['flow_concentration_index']:.4f}`",
            f"  - exchange coherence: `{row['exchange_coherence']:.4f}`",
            f"  - channel count: `{row['channel_count']}`",
            f"  - loop count: `{row['loop_count']}`",
            f"  - degree spectrum shift: `{row['degree_spectrum_shift']:.4f}`",
            f"  - clustering coefficient shift: `{row['clustering_coefficient_shift']:.4f}`",
        ])
    lines.extend([
        '',
        'Interpretation boundary:',
        '- this scan asks whether no-distance combinatorial exchange-like transport is robust to controlled incidence perturbations',
        '- it does not claim geometry emergence, particles, or bound-state physics',
        '',
    ])
    path.write_text('\n'.join(lines), encoding='utf-8')


def summarize_outcome(rows: list[dict[str, Any]]) -> tuple[str, str]:
    ordered = sorted(rows, key=lambda item: perturbation_order(str(item['perturbation_id'])))
    baseline = ordered[0]
    braid_count = sum(int(row['braid_like_exchange']) for row in ordered)
    smeared_count = sum(int(row['smeared_transfer']) for row in ordered)
    unresolved_count = sum(int(row['unresolved']) for row in ordered)

    if str(baseline['topology_class']) == 'braid_like_exchange':
        weak_ids = {'ultra_weak_rewiring', 'weak_rewiring'}
        weak_preserved = all(
            next(row for row in ordered if row['perturbation_id'] == perturb_id)['topology_class'] == 'braid_like_exchange'
            for perturb_id in weak_ids
        )
        if braid_count >= 7:
            observation = 'the braid-like class survives across most incidence perturbations in the 9-run scan'
            conclusion = 'Stage C0.3 reads as strong combinatorial protection on the current branch'
        elif weak_preserved and braid_count >= 3:
            observation = 'the braid-like class survives weak incidence noise but collapses under stronger perturbations'
            conclusion = 'Stage C0.3 reads as threshold protection tied to mesoscopic incidence regularity'
        else:
            observation = 'the braid-like class collapses rapidly once incidence noise is introduced'
            conclusion = 'Stage C0.3 reads as fragile under small combinatorial perturbation'
    else:
        if str(baseline['topology_class']) == 'smeared_transfer' and smeared_count >= 6 and unresolved_count <= 3:
            observation = 'the no-noise baseline remains a smeared-transfer surrogate and most perturbed runs stay inside that surrogate class'
            conclusion = 'this branch does not yet exhibit a braid-like exchange baseline; C0.3 currently shows robustness of the smeared exchange surrogate rather than protected braid topology'
        elif str(baseline['topology_class']) == 'smeared_transfer' and unresolved_count >= 4:
            observation = 'the no-noise baseline starts in the smeared-transfer surrogate but several perturbations push the run into the unresolved class'
            conclusion = 'the current combinatorial exchange surrogate is only weakly robust to incidence noise and no protected braid family is present'
        else:
            observation = 'the scan remains outside the braid-like class and does not produce a protected exchange topology under the current incidence-noise criteria'
            conclusion = 'C0.3 is negative for braid protection on the current no-distance branch'
    return observation, conclusion


def main() -> None:
    args = parse_args()
    defaults = load_stage10_defaults()
    runsheet = load_runsheet(args.runsheet)
    base = runsheet['base_seed_reference']
    runs = selected_runs(runsheet['runs'], args.run_ids)
    if not runs:
        raise SystemExit('no Stage C0.3 runs selected')

    rep_lookup = lookup_by_id(runsheet['representatives'], 'representative_id')
    perturbation_lookup = lookup_by_id(runsheet['perturbations'], 'perturbation_id')
    note_name = str(runsheet.get('note_name', 'Stage_C0_3_Topology_Exchange_Robustness_Under_Incidence_Noise_v1.md'))
    note_path = ATLAS_NOTES / note_name

    rows: list[dict[str, Any]] = []
    payload_runs: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    for run in sorted(runs, key=lambda item: perturbation_order(str(item['perturbation_id']))):
        rep = rep_lookup[run['representative_id']]
        perturbation = perturbation_lookup[run['perturbation_id']]
        payload, run_plots = simulate_run(run, rep, perturbation, runsheet, defaults, base)
        rows.append(payload['summary'])
        payload_runs.append(payload)
        plot_paths.extend(run_plots)

    rows = classify_against_baseline(rows, runsheet)
    row_lookup = {str(row['run_id']): row for row in rows}
    for payload in payload_runs:
        payload['summary'] = row_lookup[str(payload['summary']['run_id'])]

    plot_paths.extend(create_summary_plots(rows))
    observation, conclusion = summarize_outcome(rows)

    result = {
        'stage': runsheet['stage'],
        'description': runsheet['description'],
        'architecture_notice': runsheet['architecture_notice'],
        'baseline_notice': runsheet['baseline_notice'],
        'kernel_mode_family': runsheet['kernel_mode_family'],
        'base_seed_reference': base,
        'parameter_policy': runsheet['parameter_policy'],
        'runs': payload_runs,
        'topology_counts': dict(Counter(row['topology_class'] for row in rows)),
        'interaction_label_counts': dict(Counter(row['interaction_label'] for row in rows)),
        'coarse_label_counts': dict(Counter(row['coarse_label'] for row in rows)),
        'observation': observation,
        'conclusion': conclusion,
    }

    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug=str(runsheet.get('experiment_slug', 'stage_c0_3_incidence_noise_robustness')),
        result=result,
        csv_rows=rows,
        csv_fieldnames=CSV_FIELDS,
        plot_paths=plot_paths,
    )

    write_note(
        path=note_path,
        json_rel=str(json_path.relative_to(REPO_ROOT)),
        csv_rel=str(csv_path.relative_to(REPO_ROOT)),
        rows=rows,
        runsheet=runsheet,
    )

    print(f'wrote {json_path}')
    print(f'wrote {csv_path}')
    for item in stamped_plots:
        print(item)
    print(observation)
    print(conclusion)


if __name__ == '__main__':
    main()
