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

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_2_degree_spectrum_runs.json'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage_c0_2_degree_spectrum'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

FAMILY_CACHE: dict[tuple[int, str, str, int, int, int, int, int], dict[str, Any]] = {}
ANALYSIS_STRIDE = 4

CSV_FIELDS = [
    'run_id',
    'representative_id',
    'graph_family_id',
    'graph_family_label',
    'kernel_class_id',
    'topology_label',
    'interaction_label',
    'coarse_label',
    'braid_like_exchange',
    'transfer_smeared',
    'trapped',
    'dispersive',
    'defect_pinned',
    'channel_split',
    'flow_concentration_index',
    'transport_span',
    'coherence_metric',
    'overlap_persistence',
    'degree_localization_correlation',
    'defect_pinning_score',
    'hub_attraction_score',
    'hub_mass_fraction',
    'low_degree_mass_fraction',
    'degree_mean',
    'degree_std',
    'degree_min',
    'degree_max',
    'initial_mean_pair_distance',
    'min_mean_pair_distance',
    'final_mean_pair_distance',
    'max_mean_overlap',
    'phase_lock_indicator',
    'constraint_max',
    'sector_leakage',
    'operator_lambda_max',
    'operator_scale_factor',
    'delta_vs_uniform_transport_span',
    'delta_vs_uniform_flow_concentration',
    'delta_vs_uniform_degree_localization',
    'topology_changed_vs_uniform',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.2 degree-spectrum protection test.')
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


def family_order(graph_family_id: str) -> int:
    return {
        'degree_uniform_baseline': 0,
        'mild_degree_spread': 1,
        'strong_degree_heterogeneity': 2,
    }.get(graph_family_id, 99)


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


def build_mild_degree_spread(
    baseline_adj: list[list[int]],
    base_shells: list[dict[int, int]],
    operations: int,
    min_degree_floor: int,
    seed: int,
) -> tuple[list[list[int]], dict[str, Any]]:
    shell_candidates = [
        [idx for idx, shell in shell_map.items() if shell in (2, 3)]
        for shell_map in base_shells
    ]

    for attempt_id in range(8):
        rng = np.random.default_rng(seed + attempt_id)
        adj = copy_adjacency(baseline_adj)
        completed = 0
        max_attempts = max(operations * 40, 2000)
        attempts = 0
        while completed < operations and attempts < max_attempts:
            attempts += 1
            u = int(rng.integers(len(adj)))
            if len(adj[u]) <= min_degree_floor + 1:
                continue
            removable = [nbr for nbr in adj[u] if len(adj[nbr]) > min_degree_floor]
            if not removable:
                continue
            v = int(rng.choice(removable))
            candidates = [w for w in shell_candidates[u] if w not in adj[u] and w != u]
            if not candidates:
                continue
            w = int(rng.choice(candidates))
            adj[u].remove(v)
            adj[v].remove(u)
            adj[u].add(w)
            adj[w].add(u)
            completed += 1

        adjacency_lists = [sorted(items) for items in adj]
        if len(adjacency_components(adjacency_lists)) == 1:
            stats = degree_stats(adjacency_lists)
            return adjacency_lists, {
                'mutation_completed': completed,
                'attempts': attempts,
                'hub_nodes': [],
                'degrees': stats['degrees'],
                'degree_mean': stats['mean'],
                'degree_std': stats['std'],
                'degree_min': stats['min'],
                'degree_max': stats['max'],
            }
    raise RuntimeError('failed to build connected mild degree-spread graph family')


def build_strong_degree_heterogeneity(
    baseline_adj: list[list[int]],
    operations: int,
    hub_count: int,
    min_degree_floor: int,
    seed: int,
) -> tuple[list[list[int]], dict[str, Any]]:
    node_count = len(baseline_adj)
    for attempt_id in range(8):
        rng = np.random.default_rng(seed + attempt_id)
        adj = copy_adjacency(baseline_adj)
        hubs = [int(value) for value in rng.choice(node_count, size=hub_count, replace=False)]
        completed = 0
        max_attempts = max(operations * 60, 3000)
        attempts = 0
        while completed < operations and attempts < max_attempts:
            attempts += 1
            hub = int(rng.choice(hubs))
            donors = [idx for idx in range(node_count) if idx not in hubs and len(adj[idx]) > min_degree_floor]
            if not donors:
                break
            u = int(rng.choice(donors))
            movable = [
                nbr for nbr in adj[u]
                if nbr not in hubs and len(adj[nbr]) > min_degree_floor and hub not in adj[nbr] and nbr != hub
            ]
            if not movable:
                continue
            v = int(rng.choice(movable))
            adj[u].remove(v)
            adj[v].remove(u)
            adj[hub].add(v)
            adj[v].add(hub)
            completed += 1

        adjacency_lists = [sorted(items) for items in adj]
        if len(adjacency_components(adjacency_lists)) == 1:
            stats = degree_stats(adjacency_lists)
            return adjacency_lists, {
                'mutation_completed': completed,
                'attempts': attempts,
                'hub_nodes': hubs,
                'degrees': stats['degrees'],
                'degree_mean': stats['mean'],
                'degree_std': stats['std'],
                'degree_min': stats['min'],
                'degree_max': stats['max'],
            }
    raise RuntimeError('failed to build connected strong degree-heterogeneity graph family')


def build_graph_family_operator(
    data: Any,
    graph_family_id: str,
    params: dict[str, Any],
    target_lambda: float,
) -> dict[str, Any]:
    max_shell = int(params.get('max_shell', 3))
    graph_shell_weights = {
        int(key): float(value)
        for key, value in dict(params.get('graph_shell_weights', {'1': 1.0, '2': 0.5, '3': 0.25})).items()
    }
    random_seed = int(params.get('random_seed', 2202))
    mild_ops = int(params.get('mild_rewire_operations', 480))
    strong_ops = int(params.get('strong_rewire_operations', 2400))
    hub_count = int(params.get('strong_hub_count', 4))
    min_degree_floor = int(params.get('min_degree_floor', 3))

    cache_key = (
        int(data.n_side),
        str(data.boundary_type),
        str(graph_family_id),
        max_shell,
        mild_ops,
        strong_ops,
        hub_count,
        random_seed,
    )
    if cache_key in FAMILY_CACHE:
        return FAMILY_CACHE[cache_key]

    baseline = build_line_graph_cache(data, max_shell=max_shell)
    baseline_adj = baseline['adjacency_lists']
    base_shells = baseline['shells']

    if graph_family_id == 'degree_uniform_baseline':
        adjacency_lists = [list(items) for items in baseline_adj]
        stats = degree_stats(adjacency_lists)
        family_meta = {
            'mutation_completed': 0,
            'attempts': 0,
            'hub_nodes': [],
            'degrees': stats['degrees'],
            'degree_mean': stats['mean'],
            'degree_std': stats['std'],
            'degree_min': stats['min'],
            'degree_max': stats['max'],
        }
    elif graph_family_id == 'mild_degree_spread':
        adjacency_lists, family_meta = build_mild_degree_spread(
            baseline_adj=baseline_adj,
            base_shells=base_shells,
            operations=mild_ops,
            min_degree_floor=min_degree_floor,
            seed=random_seed + 100,
        )
    elif graph_family_id == 'strong_degree_heterogeneity':
        adjacency_lists, family_meta = build_strong_degree_heterogeneity(
            baseline_adj=baseline_adj,
            operations=strong_ops,
            hub_count=hub_count,
            min_degree_floor=min_degree_floor,
            seed=random_seed + 200,
        )
    else:
        raise ValueError(f'unsupported graph family: {graph_family_id}')

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

    family_payload = {
        'graph_family_id': graph_family_id,
        'adjacency_lists': adjacency_lists,
        'shells': shells,
        'operator': operator,
        'weight_matrix': weight_matrix,
        'graph_degrees': np.asarray(family_meta['degrees'], dtype=float),
        'weighted_degree_vec': weighted_degree_vec,
        'hub_nodes': [int(value) for value in family_meta.get('hub_nodes', [])],
        'degree_mean': float(family_meta['degree_mean']),
        'degree_std': float(family_meta['degree_std']),
        'degree_min': float(family_meta['degree_min']),
        'degree_max': float(family_meta['degree_max']),
        'mutation_completed': int(family_meta['mutation_completed']),
        'mutation_attempts': int(family_meta['attempts']),
        'raw_lambda_max': float(raw_lambda),
        'lambda_max': float(lambda_max),
        'scale_factor': float(scale_factor),
        'support_entries': int(weight_matrix.nnz),
        'max_shell': max_shell,
    }
    FAMILY_CACHE[cache_key] = family_payload
    return family_payload


def classify_topology(
    interaction_label: str,
    transport_span: float,
    binding_persistence: float,
    max_overlap: float,
    mean_basin_count_sigma2: float,
    defect_pinning_score: float,
    hub_attraction_score: float,
) -> str:
    if defect_pinning_score >= 0.22 and hub_attraction_score >= 0.6 and transport_span <= 0.08:
        return 'defect_pinned'
    if interaction_label == 'metastable composite regime' or (binding_persistence >= 0.75 and transport_span <= 0.08):
        return 'trapped'
    if mean_basin_count_sigma2 >= 1.2:
        return 'channel_split'
    if transport_span >= 0.12 and binding_persistence < 0.2 and max_overlap < 0.35:
        return 'dispersive'
    return 'transfer_smeared'


def overlay_plot(
    path: Path,
    run_id: str,
    midpoints: np.ndarray,
    degree_values: np.ndarray,
    avg_weights: np.ndarray,
    boundary_type: str,
) -> None:
    coords = np.asarray(midpoints, dtype=float)
    slice_value = 0.5 if boundary_type == 'periodic' else 0.5
    distances = np.abs(coords[:, 2] - slice_value)
    threshold = np.min(distances) + 1.0e-9
    mask = distances <= threshold
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.6))
    sc0 = axes[0].scatter(coords[mask, 0], coords[mask, 1], c=degree_values[mask], s=18, cmap='plasma')
    axes[0].set_title('Node degree overlay')
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
    graph_family: dict[str, Any],
    runsheet: dict[str, Any],
    defaults: dict[str, Any],
    base: dict[str, Any],
) -> tuple[dict[str, Any], list[Path]]:
    resolution = int(run['resolution'])
    boundary_type = str(base['boundary_type'])
    data, projector = get_cached_setup(resolution, defaults, boundary_type)
    target_lambda = estimate_lambda_max(data.L1) if bool(runsheet.get('parameter_policy', {}).get('match_lambda_to_frozen_L1', True)) else 0.0
    family_info = build_graph_family_operator(
        data=data,
        graph_family_id=str(run['graph_family_id']),
        params=runsheet.get('parameter_policy', {}),
        target_lambda=target_lambda,
    )
    operator = family_info['operator']
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

    graph_degrees = np.asarray(family_info['graph_degrees'], dtype=float)
    degree_std = float(np.std(graph_degrees))
    defect_quantile = float(runsheet.get('parameter_policy', {}).get('defect_quantile', 0.1))
    if degree_std > 1.0e-12:
        count = max(1, int(round(defect_quantile * len(graph_degrees))))
        order = np.argsort(graph_degrees, kind='mergesort')
        high_mask = np.zeros_like(graph_degrees, dtype=bool)
        low_mask = np.zeros_like(graph_degrees, dtype=bool)
        low_mask[order[:count]] = True
        high_mask[order[-count:]] = True
    else:
        high_mask = np.zeros_like(graph_degrees, dtype=bool)
        low_mask = np.zeros_like(graph_degrees, dtype=bool)

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
    mean_degree_trace: list[float] = []
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

        norm_weights = total_weights / max(float(np.sum(total_weights)), 1.0e-12)
        mean_degree_trace.append(float(np.sum(norm_weights * graph_degrees)))

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
    coherence_metric = float(np.mean(coherences)) if coherences else 0.0
    binding_mask = (np.asarray(pair_distance_series) <= 0.25) & (np.asarray(overlap_series) >= 0.35)
    overlap_persistence = float(np.mean(binding_mask.astype(float))) if binding_mask.size else 0.0
    if degree_std > 1.0e-12 and np.std(avg_weights) > 1.0e-18:
        degree_localization_correlation = float(np.corrcoef(graph_degrees, avg_weights)[0, 1])
        hub_attraction = float((np.sum(avg_weights * graph_degrees) - np.mean(graph_degrees)) / degree_std)
    else:
        degree_localization_correlation = 0.0
        hub_attraction = 0.0
    hub_mass_fraction = float(np.sum(avg_weights[high_mask])) if np.any(high_mask) else 0.0
    low_degree_mass_fraction = float(np.sum(avg_weights[low_mask])) if np.any(low_mask) else 0.0
    defect_pinning_score = max(hub_mass_fraction, low_degree_mass_fraction)

    transport_span = path_length(centers_arr, boundary_type)
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

    topology_label = classify_topology(
        interaction_label=interaction_label,
        transport_span=transport_span,
        binding_persistence=overlap_persistence,
        max_overlap=float(np.max(overlap_series)) if overlap_series else 0.0,
        mean_basin_count_sigma2=float(coarse_summary['mean_basin_count_sigma2']),
        defect_pinning_score=defect_pinning_score,
        hub_attraction_score=hub_attraction,
    )

    row = {
        'run_id': run['run_id'],
        'representative_id': run['representative_id'],
        'graph_family_id': run['graph_family_id'],
        'graph_family_label': graph_family['label'],
        'kernel_class_id': str(runsheet['kernel_mode_family']['kernel_class_id']),
        'topology_label': topology_label,
        'interaction_label': interaction_label,
        'coarse_label': coarse_label,
        'braid_like_exchange': int(topology_label == 'braid_like_exchange'),
        'transfer_smeared': int(topology_label == 'transfer_smeared'),
        'trapped': int(topology_label == 'trapped'),
        'dispersive': int(topology_label == 'dispersive'),
        'defect_pinned': int(topology_label == 'defect_pinned'),
        'channel_split': int(topology_label == 'channel_split'),
        'flow_concentration_index': flow_concentration,
        'transport_span': transport_span,
        'coherence_metric': coherence_metric,
        'overlap_persistence': overlap_persistence,
        'degree_localization_correlation': degree_localization_correlation,
        'defect_pinning_score': defect_pinning_score,
        'hub_attraction_score': hub_attraction,
        'hub_mass_fraction': hub_mass_fraction,
        'low_degree_mass_fraction': low_degree_mass_fraction,
        'degree_mean': float(family_info['degree_mean']),
        'degree_std': float(family_info['degree_std']),
        'degree_min': float(family_info['degree_min']),
        'degree_max': float(family_info['degree_max']),
        'initial_mean_pair_distance': float(pair_distance_series[0]) if pair_distance_series else 0.0,
        'min_mean_pair_distance': float(np.min(pair_distance_series)) if pair_distance_series else 0.0,
        'final_mean_pair_distance': float(pair_distance_series[-1]) if pair_distance_series else 0.0,
        'max_mean_overlap': float(np.max(overlap_series)) if overlap_series else 0.0,
        'phase_lock_indicator': float(np.mean(phase_alignment_series)) if phase_alignment_series else 0.0,
        'constraint_max': float(np.max(constraint_norms)) if constraint_norms else 0.0,
        'sector_leakage': float(np.max(sector_leakages)) if sector_leakages else 0.0,
        'operator_lambda_max': float(family_info['lambda_max']),
        'operator_scale_factor': float(family_info['scale_factor']),
        'delta_vs_uniform_transport_span': 0.0,
        'delta_vs_uniform_flow_concentration': 0.0,
        'delta_vs_uniform_degree_localization': 0.0,
        'topology_changed_vs_uniform': 0,
        'notes': run['notes'],
    }

    plot_paths: list[Path] = []

    trace_path = WORK_PLOT_DIR / f"{run['run_id']}_coherence_span_trace.png"
    center_shift_trace = [
        float(np.linalg.norm(displacement(np.asarray([center], dtype=float), centers_arr[0], boundary_type)[0]))
        for center in centers_arr
    ]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))
    axes[0].plot(times, coherences, color='tab:blue', label='coherence')
    axes[0].set_title('Coherence trace')
    axes[0].set_xlabel('time')
    axes[0].grid(alpha=0.25)
    axes[1].plot(times, center_shift_trace, color='tab:orange', label='center shift')
    axes[1].set_title('Transport span trace')
    axes[1].set_xlabel('time')
    axes[1].grid(alpha=0.25)
    fig.suptitle(run['run_id'])
    fig.savefig(trace_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(trace_path)

    pair_path = WORK_PLOT_DIR / f"{run['run_id']}_pair_overlap_trace.png"
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))
    axes[0].plot(times, pair_distance_series, color='tab:green')
    axes[0].set_title('Mean pair distance')
    axes[0].set_xlabel('time')
    axes[0].grid(alpha=0.25)
    axes[1].plot(times, overlap_series, color='tab:red', label='overlap')
    axes[1].plot(times, phase_alignment_series, color='tab:purple', label='phase')
    axes[1].set_title('Overlap and phase')
    axes[1].set_xlabel('time')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(pair_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(pair_path)

    overlay_path = WORK_PLOT_DIR / f"{run['run_id']}_degree_flow_overlay.png"
    overlay_plot(
        path=overlay_path,
        run_id=run['run_id'],
        midpoints=data.midpoints,
        degree_values=graph_degrees,
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
        'graph_family': graph_family,
        'graph_family_summary': {
            'degree_mean': float(family_info['degree_mean']),
            'degree_std': float(family_info['degree_std']),
            'degree_min': float(family_info['degree_min']),
            'degree_max': float(family_info['degree_max']),
            'mutation_completed': int(family_info['mutation_completed']),
            'mutation_attempts': int(family_info['mutation_attempts']),
            'hub_nodes': [int(value) for value in family_info['hub_nodes']],
        },
        'operator_summary': {
            'lambda_max': float(family_info['lambda_max']),
            'raw_lambda_max': float(family_info['raw_lambda_max']),
            'scale_factor': float(family_info['scale_factor']),
            'support_entries': int(family_info['support_entries']),
            'max_shell': int(family_info['max_shell']),
        },
        'summary': row,
        'times': times,
        'pair_distance_series': pair_distance_series,
        'overlap_series': overlap_series,
        'phase_alignment_series': phase_alignment_series,
        'coherences': coherences,
        'mean_degree_trace': mean_degree_trace,
    }
    return result, plot_paths


def compare_to_uniform(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    uniform_map = {
        str(row['representative_id']): row
        for row in rows
        if row['graph_family_id'] == 'degree_uniform_baseline'
    }
    for row in rows:
        baseline = uniform_map.get(str(row['representative_id']))
        if baseline is None:
            continue
        row['delta_vs_uniform_transport_span'] = float(row['transport_span']) - float(baseline['transport_span'])
        row['delta_vs_uniform_flow_concentration'] = float(row['flow_concentration_index']) - float(baseline['flow_concentration_index'])
        row['delta_vs_uniform_degree_localization'] = float(row['degree_localization_correlation']) - float(baseline['degree_localization_correlation'])
        row['topology_changed_vs_uniform'] = int(str(row['topology_label']) != str(baseline['topology_label']))
    return rows


def create_summary_plots(rows: list[dict[str, Any]]) -> list[Path]:
    plot_paths: list[Path] = []
    reps = sorted({str(row['representative_id']) for row in rows})
    families = sorted({str(row['graph_family_id']) for row in rows}, key=family_order)
    family_labels = {
        'degree_uniform_baseline': 'uniform',
        'mild_degree_spread': 'mild',
        'strong_degree_heterogeneity': 'strong',
    }
    row_map = {(str(row['representative_id']), str(row['graph_family_id'])): row for row in rows}
    x = np.arange(len(reps))
    width = 0.22
    offsets = np.arange(len(families), dtype=float) - 0.5 * (len(families) - 1)

    topology_scores = {
        'transfer_smeared': 0,
        'trapped': 1,
        'dispersive': 2,
        'defect_pinned': 3,
        'channel_split': 4,
        'braid_like_exchange': 5,
    }

    matrix_path = WORK_PLOT_DIR / 'stage_c0_2_topology_degree_matrix.png'
    data = np.zeros((len(reps), len(families)), dtype=float)
    for i, rep in enumerate(reps):
        for j, family in enumerate(families):
            row = row_map[(rep, family)]
            data[i, j] = float(topology_scores.get(str(row['topology_label']), 0))
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    im = ax.imshow(data, cmap='viridis', vmin=0, vmax=max(topology_scores.values()))
    ax.set_xticks(range(len(families)))
    ax.set_xticklabels([family_labels.get(f, f) for f in families])
    ax.set_yticks(range(len(reps)))
    ax.set_yticklabels(reps)
    ax.set_title('Topology vs degree-spectrum matrix')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(matrix_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(matrix_path)

    concentration_path = WORK_PLOT_DIR / 'stage_c0_2_flow_concentration_by_family.png'
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.8))
    for idx, family in enumerate(families):
        concentration = []
        transport = []
        for rep in reps:
            row = row_map[(rep, family)]
            concentration.append(float(row['flow_concentration_index']))
            transport.append(float(row['transport_span']))
        axes[0].bar(x + offsets[idx] * width, concentration, width=width, label=family_labels.get(family, family))
        axes[1].bar(x + offsets[idx] * width, transport, width=width, label=family_labels.get(family, family))
    axes[0].set_title('Flow concentration by graph family')
    axes[1].set_title('Transport span by graph family')
    for axis in axes:
        axis.set_xticks(x)
        axis.set_xticklabels(reps, rotation=20, ha='right')
        axis.grid(alpha=0.25, axis='y')
        axis.legend(fontsize=8)
    fig.savefig(concentration_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(concentration_path)

    correlation_path = WORK_PLOT_DIR / 'stage_c0_2_degree_localization_correlation.png'
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.8))
    for idx, family in enumerate(families):
        corr_values = []
        hub_values = []
        for rep in reps:
            row = row_map[(rep, family)]
            corr_values.append(float(row['degree_localization_correlation']))
            hub_values.append(float(row['hub_attraction_score']))
        axes[0].bar(x + offsets[idx] * width, corr_values, width=width, label=family_labels.get(family, family))
        axes[1].bar(x + offsets[idx] * width, hub_values, width=width, label=family_labels.get(family, family))
    axes[0].axhline(0.0, color='k', linewidth=1.0)
    axes[1].axhline(0.0, color='k', linewidth=1.0)
    axes[0].set_title('Degree-localization correlation')
    axes[1].set_title('Hub attraction / avoidance score')
    for axis in axes:
        axis.set_xticks(x)
        axis.set_xticklabels(reps, rotation=20, ha='right')
        axis.grid(alpha=0.25, axis='y')
        axis.legend(fontsize=8)
    fig.savefig(correlation_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(correlation_path)

    pinning_path = WORK_PLOT_DIR / 'stage_c0_2_hub_defect_pinning_summary.png'
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.8))
    for idx, family in enumerate(families):
        defect_values = []
        hub_mass_values = []
        for rep in reps:
            row = row_map[(rep, family)]
            defect_values.append(float(row['defect_pinning_score']))
            hub_mass_values.append(float(row['hub_mass_fraction']))
        axes[0].bar(x + offsets[idx] * width, defect_values, width=width, label=family_labels.get(family, family))
        axes[1].bar(x + offsets[idx] * width, hub_mass_values, width=width, label=family_labels.get(family, family))
    axes[0].set_title('Defect pinning score')
    axes[1].set_title('Hub mass fraction')
    for axis in axes:
        axis.set_xticks(x)
        axis.set_xticklabels(reps, rotation=20, ha='right')
        axis.grid(alpha=0.25, axis='y')
        axis.legend(fontsize=8)
    fig.savefig(pinning_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(pinning_path)

    return plot_paths


def write_note(
    path: Path,
    json_rel: str,
    csv_rel: str,
    rows: list[dict[str, Any]],
    runsheet: dict[str, Any],
) -> None:
    lines = [
        '# Stage C0.2 Degree Spectrum Protection Test v1',
        '',
        f'Timestamped summary: `{json_rel}`',
        f'Timestamped run table: `{csv_rel}`',
        '',
        f'Architecture notice: {runsheet["architecture_notice"]}',
        '',
        'Per-run summary:',
        '',
    ]
    for row in sorted(rows, key=lambda item: (item['representative_id'], family_order(str(item['graph_family_id'])))):
        lines.extend([
            f"- `{row['representative_id']}` / `{row['graph_family_label']}`",
            f"  - topology label: `{row['topology_label']}`",
            f"  - interaction label: `{row['interaction_label']}`",
            f"  - flow concentration: `{row['flow_concentration_index']:.4f}`",
            f"  - transport span: `{row['transport_span']:.4f}`",
            f"  - degree-localization correlation: `{row['degree_localization_correlation']:.4f}`",
            f"  - defect pinning score: `{row['defect_pinning_score']:.4f}`",
            f"  - hub attraction score: `{row['hub_attraction_score']:.4f}`",
        ])
    lines.extend([
        '',
        'Interpretation boundary:',
        '- this scan asks whether combinatorial transport topology is protected by graph degree-spectrum structure',
        '- it does not claim geometry emergence, gravity, or particle trapping',
        '',
    ])
    path.write_text('\n'.join(lines), encoding='utf-8')


def summarize_outcome(rows: list[dict[str, Any]]) -> tuple[str, str]:
    reps = sorted({str(row['representative_id']) for row in rows})
    row_map = {(str(row['representative_id']), str(row['graph_family_id'])): row for row in rows}
    required_families = {'degree_uniform_baseline', 'mild_degree_spread', 'strong_degree_heterogeneity'}
    present_families = {str(row['graph_family_id']) for row in rows}
    if not required_families.issubset(present_families):
        observation = 'partial degree-spectrum subset executed successfully'
        conclusion = 'full protection interpretation requires the complete uniform/mild/strong matrix'
        return observation, conclusion
    mild_preserved = all(
        str(row_map[(rep, 'mild_degree_spread')]['topology_label']) == str(row_map[(rep, 'degree_uniform_baseline')]['topology_label'])
        for rep in reps
    )
    strong_changed = any(
        str(row_map[(rep, 'strong_degree_heterogeneity')]['topology_label']) != str(row_map[(rep, 'degree_uniform_baseline')]['topology_label'])
        for rep in reps
    )
    pinning_detected = any(float(row['defect_pinning_score']) >= 0.2 for row in rows if row['graph_family_id'] != 'degree_uniform_baseline')
    degree_signal = any(abs(float(row['degree_localization_correlation'])) >= 0.1 or abs(float(row['hub_attraction_score'])) >= 0.5 for row in rows if row['graph_family_id'] != 'degree_uniform_baseline')
    transport_sensitive = any(
        abs(float(row['delta_vs_uniform_transport_span'])) >= 0.05
        or abs(float(row['delta_vs_uniform_flow_concentration'])) >= 5.0
        for row in rows
        if row['graph_family_id'] != 'degree_uniform_baseline'
    )

    if mild_preserved and strong_changed:
        observation = 'topology classes survive mild degree spread but at least one representative changes class under strong degree heterogeneity'
        conclusion = 'combinatorial transport topology is protected by near-regular degree structure and becomes degree-sensitive once heterogeneity is strong enough'
    elif mild_preserved and transport_sensitive:
        observation = 'topology labels remain stable across the scan, but transport span or flow concentration shifts materially once the degree spectrum broadens'
        conclusion = 'degree spectrum is not overturning the coarse topology classes here, but it is reshaping transport texture and localization inside those classes'
    elif pinning_detected or degree_signal:
        observation = 'degree-sensitive localization and pinning diagnostics appear on the heterogeneous graph families even when the coarse class does not fully collapse'
        conclusion = 'the combinatorial branch shows degree-spectrum sensitivity through flow localization and hub/defect response'
    else:
        observation = 'all graph families remain in the same topology envelope with weak correlation to degree structure'
        conclusion = 'degree spectrum is not the active protection axis in this first combinatorial scan'
    return observation, conclusion


def main() -> None:
    args = parse_args()
    defaults = load_stage10_defaults()
    runsheet = load_runsheet(args.runsheet)
    base = runsheet['base_seed_reference']
    runs = selected_runs(runsheet['runs'], args.run_ids)
    if not runs:
        raise SystemExit('no Stage C0.2 runs selected')

    rep_lookup = lookup_by_id(runsheet['representatives'], 'representative_id')
    family_lookup = lookup_by_id(runsheet['graph_families'], 'graph_family_id')
    note_name = str(runsheet.get('note_name', 'Stage_C0_2_Degree_Spectrum_Protection_Test_v1.md'))
    note_path = ATLAS_NOTES / note_name

    rows: list[dict[str, Any]] = []
    payload_runs: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    for run in sorted(runs, key=lambda item: (rep_lookup[item['representative_id']]['representative_id'], family_order(str(item['graph_family_id'])))):
        rep = rep_lookup[run['representative_id']]
        graph_family = family_lookup[run['graph_family_id']]
        payload, run_plots = simulate_run(run, rep, graph_family, runsheet, defaults, base)
        rows.append(payload['summary'])
        payload_runs.append(payload)
        plot_paths.extend(run_plots)

    rows = compare_to_uniform(rows)
    plot_paths.extend(create_summary_plots(rows))
    observation, conclusion = summarize_outcome(rows)

    result = {
        'stage': runsheet['stage'],
        'description': runsheet['description'],
        'architecture_notice': runsheet['architecture_notice'],
        'kernel_mode_family': runsheet['kernel_mode_family'],
        'base_seed_reference': base,
        'parameter_policy': runsheet['parameter_policy'],
        'runs': payload_runs,
        'topology_counts': dict(Counter(row['topology_label'] for row in rows)),
        'interaction_label_counts': dict(Counter(row['interaction_label'] for row in rows)),
        'coarse_label_counts': dict(Counter(row['coarse_label'] for row in rows)),
        'observation': observation,
        'conclusion': conclusion,
    }

    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug=str(runsheet.get('experiment_slug', 'stage_c0_2_degree_spectrum_protection')),
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
