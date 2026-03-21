#!/usr/bin/env python3

from __future__ import annotations

import argparse
import itertools
import json
import math
from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import Any

import numpy as np

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
from stage_c0_3_incidence_noise_robustness import (
    active_graph_metrics,
    adjacency_edge_set,
    build_operator_from_adjacency,
    recurrence_indicator,
)

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_4_motif_injection_runs.json'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage_c0_4_controlled_motif_injection'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

ANALYSIS_STRIDE = 4
BASELINE_CACHE: dict[tuple[str, int, str], dict[str, Any]] = {}
MOTIF_CACHE: dict[tuple[str, str, int, str], dict[str, Any]] = {}

CSV_FIELDS = [
    'run_id',
    'representative_id',
    'motif_id',
    'motif_label',
    'baseline_topology',
    'topology_label',
    'effect_class',
    'interaction_label',
    'coarse_label',
    'braid_like_exchange',
    'smeared_transfer',
    'split_channel_exchange',
    'defect_pinned',
    'trapped_local',
    'unresolved',
    'flow_concentration_index',
    'exchange_coherence',
    'transport_span',
    'channel_count',
    'loop_count',
    'motif_occupancy_correlation',
    'motif_mass_fraction',
    'pinning_score',
    'routing_asymmetry',
    'recurrence_indicator',
    'max_mean_overlap',
    'degree_spectrum_shift',
    'clustering_coefficient_shift',
    'delta_vs_baseline_concentration',
    'delta_vs_baseline_coherence',
    'delta_vs_baseline_channel_count',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.4 controlled motif injection scan.')
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


def motif_order(motif_id: str) -> int:
    return {
        'triangle_cluster': 0,
        'diamond': 1,
        'hub_micro_star': 2,
    }.get(motif_id, 99)


def copy_adjacency(adjacency_lists: list[list[int]]) -> list[set[int]]:
    return [set(items) for items in adjacency_lists]


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


def adjacency_connected_sets(adjacency_sets: list[set[int]]) -> bool:
    if not adjacency_sets:
        return True
    visited = {0}
    queue: deque[int] = deque([0])
    while queue:
        cur = queue.popleft()
        for nbr in adjacency_sets[cur]:
            if nbr not in visited:
                visited.add(nbr)
                queue.append(nbr)
    return len(visited) == len(adjacency_sets)


def remove_edge(adj: list[set[int]], u: int, v: int) -> None:
    adj[u].remove(v)
    adj[v].remove(u)


def add_edge(adj: list[set[int]], u: int, v: int) -> None:
    adj[u].add(v)
    adj[v].add(u)


def edge_key(u: int, v: int) -> tuple[int, int]:
    return (u, v) if u < v else (v, u)


def edge_support(adjacency_sets: list[set[int]], u: int, v: int) -> int:
    return len(adjacency_sets[u].intersection(adjacency_sets[v]))


def choose_seed_patch(
    adjacency_lists: list[list[int]],
    avg_weights: np.ndarray,
    min_nodes: int,
    params: dict[str, Any],
) -> tuple[int, list[int], list[int]]:
    patch_radius = int(params.get('motif_patch_radius', 2))
    max_radius = int(params.get('motif_max_radius', 4))
    candidate_count = int(params.get('motif_candidate_count', 8))
    order = np.argsort(avg_weights)[::-1]
    for seed in order[: min(len(order), 48)]:
        for radius in range(patch_radius, max_radius + 1):
            patch = bfs_patch(adjacency_lists, int(seed), radius)
            if len(patch) < min_nodes:
                continue
            ranked = sorted(patch, key=lambda idx: float(avg_weights[idx]), reverse=True)
            return int(seed), patch, ranked[: max(candidate_count, min_nodes)]
    seed = int(order[0])
    patch = bfs_patch(adjacency_lists, seed, max_radius)
    ranked = sorted(patch, key=lambda idx: float(avg_weights[idx]), reverse=True)
    return seed, patch, ranked[: max(candidate_count, min_nodes)]


def triangle_specs(adjacency_lists: list[list[int]], avg_weights: np.ndarray, params: dict[str, Any]) -> list[dict[str, Any]]:
    seed, patch, candidates = choose_seed_patch(adjacency_lists, avg_weights, min_nodes=3, params=params)
    specs: list[tuple[float, dict[str, Any]]] = []
    for triple in itertools.combinations(candidates, 3):
        desired = {edge_key(triple[0], triple[1]), edge_key(triple[0], triple[2]), edge_key(triple[1], triple[2])}
        missing = [edge for edge in desired if edge[1] not in adjacency_lists[edge[0]]]
        if not missing:
            continue
        score = sum(float(avg_weights[idx]) for idx in triple) - 0.15 * max(len(missing) - 1, 0)
        specs.append(
            (
                score,
                {
                    'motif_nodes': sorted(int(idx) for idx in triple),
                    'desired_edges': sorted(desired),
                    'seed_nodes': [seed],
                    'patch_nodes': patch,
                    'missing_edges': sorted(missing),
                },
            )
        )
    specs.sort(key=lambda item: item[0], reverse=True)
    return [item[1] for item in specs[:20]]


def diamond_specs(adjacency_lists: list[list[int]], avg_weights: np.ndarray, params: dict[str, Any]) -> list[dict[str, Any]]:
    seed, patch, candidates = choose_seed_patch(adjacency_lists, avg_weights, min_nodes=4, params=params)
    subset = candidates[: min(len(candidates), 7)]
    specs: list[tuple[float, dict[str, Any]]] = []
    seen: set[tuple[tuple[int, int], ...]] = set()
    for quad in itertools.combinations(subset, 4):
        for perm in itertools.permutations(quad, 4):
            a, b, c, d = [int(value) for value in perm]
            desired = {
                edge_key(a, b),
                edge_key(a, c),
                edge_key(b, c),
                edge_key(b, d),
                edge_key(c, d),
            }
            if tuple(sorted(desired)) in seen:
                continue
            seen.add(tuple(sorted(desired)))
            missing = [edge for edge in desired if edge[1] not in adjacency_lists[edge[0]]]
            if not missing:
                continue
            if len(missing) > 3:
                continue
            score = sum(float(avg_weights[idx]) for idx in quad) - 0.1 * len(missing)
            specs.append(
                (
                    score,
                    {
                        'motif_nodes': sorted(int(idx) for idx in quad),
                        'desired_edges': sorted(desired),
                        'seed_nodes': [seed, b, c],
                        'patch_nodes': patch,
                        'missing_edges': sorted(missing),
                    },
                )
            )
    specs.sort(key=lambda item: item[0], reverse=True)
    return [item[1] for item in specs[:20]]


def hub_specs(adjacency_lists: list[list[int]], avg_weights: np.ndarray, params: dict[str, Any]) -> list[dict[str, Any]]:
    patch_radius = int(params.get('motif_patch_radius', 2))
    max_radius = int(params.get('motif_max_radius', 4))
    leaf_count = int(params.get('hub_leaf_count', 4))
    order = np.argsort(avg_weights)[::-1]
    specs: list[tuple[float, dict[str, Any]]] = []
    for center in order[: min(len(order), 32)]:
        center = int(center)
        for radius in range(patch_radius, max_radius + 1):
            patch = bfs_patch(adjacency_lists, center, radius)
            leaves = [
                int(idx)
                for idx in sorted(patch, key=lambda item: float(avg_weights[item]), reverse=True)
                if idx != center and idx not in adjacency_lists[center]
            ]
            if len(leaves) < 2:
                continue
            leaves = leaves[: min(len(leaves), leaf_count)]
            desired = {edge_key(center, leaf) for leaf in leaves}
            missing = [edge for edge in desired if edge[1] not in adjacency_lists[edge[0]]]
            if not missing:
                continue
            score = float(avg_weights[center]) + sum(float(avg_weights[idx]) for idx in leaves)
            specs.append(
                (
                    score,
                    {
                        'motif_nodes': sorted([center] + leaves),
                        'desired_edges': sorted(desired),
                        'seed_nodes': [center],
                        'patch_nodes': patch,
                        'missing_edges': sorted(missing),
                    },
                )
            )
            break
    specs.sort(key=lambda item: item[0], reverse=True)
    return [item[1] for item in specs[:20]]


def ranked_removal_candidates(
    adjacency_sets: list[set[int]],
    protected_nodes: set[int],
    protected_edges: set[tuple[int, int]],
    avoid_nodes: set[int],
    min_degree_floor: int,
) -> list[tuple[int, int, tuple[int, int]]]:
    candidates: list[tuple[int, int, tuple[int, int]]] = []
    for u in range(len(adjacency_sets)):
        for v in adjacency_sets[u]:
            if v <= u:
                continue
            edge = (u, v)
            if edge in protected_edges:
                continue
            if u in protected_nodes and v in protected_nodes:
                continue
            if len(adjacency_sets[u]) <= min_degree_floor or len(adjacency_sets[v]) <= min_degree_floor:
                continue
            outside_bonus = int(u not in avoid_nodes and v not in avoid_nodes)
            support = edge_support(adjacency_sets, u, v)
            candidates.append((outside_bonus, support, edge))
    candidates.sort(key=lambda item: (item[0], item[1]), reverse=True)
    return candidates


def inject_balanced_motif(
    baseline_adj: list[list[int]],
    spec: dict[str, Any],
    min_degree_floor: int,
) -> dict[str, Any] | None:
    adjacency_sets = copy_adjacency(baseline_adj)
    protected_nodes = set(int(value) for value in spec['motif_nodes'])
    patch_nodes = set(int(value) for value in spec['patch_nodes'])
    protected_edges = set(edge_key(*edge) for edge in spec['desired_edges'])
    added_edges: list[tuple[int, int]] = []
    removed_edges: list[tuple[int, int]] = []
    for u, v in spec['desired_edges']:
        if v in adjacency_sets[u]:
            continue
        add_edge(adjacency_sets, u, v)
        chosen_remove: tuple[int, int] | None = None
        for _outside_bonus, _support, edge in ranked_removal_candidates(
            adjacency_sets,
            protected_nodes=protected_nodes,
            protected_edges=protected_edges.union({edge_key(u, v)}),
            avoid_nodes=patch_nodes,
            min_degree_floor=min_degree_floor,
        ):
            x, y = edge
            remove_edge(adjacency_sets, x, y)
            if adjacency_connected_sets(adjacency_sets):
                chosen_remove = edge
                break
            add_edge(adjacency_sets, x, y)
        if chosen_remove is None:
            remove_edge(adjacency_sets, u, v)
            for x, y in removed_edges:
                add_edge(adjacency_sets, x, y)
            for x, y in added_edges:
                if y in adjacency_sets[x]:
                    remove_edge(adjacency_sets, x, y)
            return None
        added_edges.append(edge_key(u, v))
        removed_edges.append(chosen_remove)
    adjacency_lists = [sorted(items) for items in adjacency_sets]
    return {
        'adjacency_lists': adjacency_lists,
        'motif_nodes': sorted(protected_nodes),
        'seed_nodes': [int(value) for value in spec['seed_nodes']],
        'patch_nodes': sorted(patch_nodes),
        'desired_edges': [list(edge) for edge in sorted(protected_edges)],
        'added_edges': [list(edge) for edge in added_edges],
        'removed_edges': [list(edge) for edge in removed_edges],
    }


def build_motif_injection(
    data: Any,
    representative_id: str,
    motif_id: str,
    baseline_avg_weights: np.ndarray,
    params: dict[str, Any],
) -> dict[str, Any]:
    cache_key = (representative_id, motif_id, int(data.n_side), str(data.boundary_type))
    if cache_key in MOTIF_CACHE:
        return MOTIF_CACHE[cache_key]

    baseline_adj = build_line_graph_cache(data, max_shell=int(params.get('max_shell', 3)))['adjacency_lists']
    if motif_id == 'triangle_cluster':
        specs = triangle_specs(baseline_adj, baseline_avg_weights, params)
    elif motif_id == 'diamond':
        specs = diamond_specs(baseline_adj, baseline_avg_weights, params)
    elif motif_id == 'hub_micro_star':
        specs = hub_specs(baseline_adj, baseline_avg_weights, params)
    else:
        raise ValueError(f'unsupported motif_id: {motif_id}')

    min_degree_floor = int(params.get('min_degree_floor', 3))
    for spec in specs:
        payload = inject_balanced_motif(baseline_adj, spec, min_degree_floor=min_degree_floor)
        if payload is not None:
            MOTIF_CACHE[cache_key] = payload
            return payload
    raise RuntimeError(f'failed to inject motif {motif_id} for {representative_id}')


def exchange_coherence_metric(overlap_series: list[float], phase_alignment_series: list[float]) -> float:
    if not overlap_series and not phase_alignment_series:
        return 0.0
    overlap_term = float(np.max(overlap_series)) if overlap_series else 0.0
    phase_term = float(np.mean(phase_alignment_series)) if phase_alignment_series else 0.0
    return float(0.5 * overlap_term + 0.5 * phase_term)


def measure_configuration(
    rep: dict[str, Any],
    data: Any,
    projector: Any,
    operator_info: dict[str, Any],
    defaults: dict[str, Any],
    base: dict[str, Any],
    motif_nodes: list[int] | None = None,
) -> dict[str, Any]:
    operator = operator_info['operator']
    boundary_type = str(base['boundary_type'])
    dt = suggested_dt(operator, float(defaults['dt_safety']))
    steps = max(2, int(np.ceil(float(defaults['t_final']) / dt)))
    leakage_fn = make_leakage_fn(projector)
    kick_axis = int(defaults['kick_axis'])
    motif_idx = np.asarray(sorted(motif_nodes or []), dtype=int)

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
    sector_leakages: list[float] = []
    constraint_norms: list[float] = []
    coherences: list[float] = []
    pair_distance_series: list[float] = []
    overlap_series: list[float] = []
    phase_alignment_series: list[float] = []
    occupancy_accumulator = np.zeros(len(data.edges), dtype=float)
    component_occupancy_accumulator = [np.zeros(len(data.edges), dtype=float) for _ in packet_qs]
    occupancy_samples = 0
    motif_mass_trace: list[float] = []
    component_motif_mass_trace: list[list[float]] = [[] for _ in packet_qs]
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
        total_norm = max(float(np.sum(total_weights)), 1.0e-12)
        occupancy_accumulator[:] += total_weights
        occupancy_samples += 1
        total_center = weighted_center(data.midpoints, total_weights, boundary_type)
        times.append(float(t))
        total_centers.append(total_center.tolist())
        widths.append(packet_width(data.midpoints, total_weights, total_center, boundary_type))
        _energy = state_energy(operator, total_q, total_v)
        _anisotropy = anisotropy_ratio(data.midpoints, total_weights, total_center, boundary_type)
        sector_leakages.append(float(leakage_fn(total_q)))
        constraint_norms.append(float(np.linalg.norm(data.d0.T @ total_q)))
        coherences.append(coherence_score(total_weights))

        component_centers: list[np.ndarray] = []
        for idx, comp_q in enumerate(packet_qs):
            comp_weights = np.abs(comp_q) ** 2
            component_occupancy_accumulator[idx][:] += comp_weights
            comp_norm = max(float(np.sum(comp_weights)), 1.0e-12)
            component_centers.append(weighted_center(data.midpoints, comp_weights, boundary_type))
            if motif_idx.size:
                component_motif_mass_trace[idx].append(float(np.sum(comp_weights[motif_idx]) / comp_norm))
        if motif_idx.size:
            motif_mass_trace.append(float(np.sum(total_weights[motif_idx]) / total_norm))

        mean_pair, _min_pair = pairwise_distance_stats(component_centers, boundary_type)
        pair_distance_series.append(mean_pair)
        overlap_series.append(pairwise_overlap_mean(packet_qs))
        phase_alignment_series.append(pairwise_phase_alignment(packet_qs))

        point_grid = edge_field_to_grid(data.midpoints, np.abs(total_q), int(data.n_side))
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
    component_avg_weights: list[np.ndarray] = []
    for accum in component_occupancy_accumulator:
        normed = accum / max(occupancy_samples, 1)
        normed = normed / max(float(np.sum(normed)), 1.0e-12)
        component_avg_weights.append(normed)

    flow_concentration = float(len(avg_weights) * coherence_score(avg_weights))
    exchange_coherence = exchange_coherence_metric(overlap_series, phase_alignment_series)
    transport_span = path_length(centers_arr, boundary_type)
    channel_count, loop_count = active_graph_metrics(
        avg_weights,
        operator_info['adjacency_lists'],
        float(1.0),
    )
    recurrence = recurrence_indicator(np.asarray(pair_distance_series, dtype=float), close_threshold=0.18)
    binding_mask = (np.asarray(pair_distance_series) <= 0.25) & (np.asarray(overlap_series) >= 0.35)
    overlap_persistence = float(np.mean(binding_mask.astype(float))) if binding_mask.size else 0.0

    coarse_summary: dict[str, float] = {'t_final': float(times[-1])}
    for sigma in SIGMAS:
        coarse_summary[f'mean_basin_count_sigma{int(sigma)}'] = float(np.mean(coarse[sigma]['basin_counts'])) if coarse[sigma]['basin_counts'] else 0.0
        coarse_summary[f'mean_dominant_area_fraction_sigma{int(sigma)}'] = float(np.mean(coarse[sigma]['dominant_area_fractions'])) if coarse[sigma]['dominant_area_fractions'] else 0.0
        coarse_summary[f'coarse_persistence_sigma{int(sigma)}'] = float(np.mean(coarse[sigma]['stable_flags'])) if coarse[sigma]['stable_flags'] else 0.0
        coarse_summary[f'mean_envelope_variance_sigma{int(sigma)}'] = float(np.mean(coarse[sigma]['envelope_variances'])) if coarse[sigma]['envelope_variances'] else 0.0
        coarse_summary[f'basin_lifetime_sigma{int(sigma)}'] = float(np.sum(np.asarray(coarse[sigma]['dominant_area_fractions']) >= 0.05) * sample_dt)
    coarse_label = classify_coarse(coarse_summary)

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

    return {
        'times': times,
        'total_centers': total_centers,
        'avg_weights': avg_weights,
        'component_avg_weights': component_avg_weights,
        'pair_distance_series': pair_distance_series,
        'overlap_series': overlap_series,
        'phase_alignment_series': phase_alignment_series,
        'coherences': coherences,
        'motif_mass_trace': motif_mass_trace,
        'component_motif_mass_trace': component_motif_mass_trace,
        'interaction_label': interaction_label,
        'coarse_label': coarse_label,
        'flow_concentration_index': flow_concentration,
        'exchange_coherence': exchange_coherence,
        'transport_span': transport_span,
        'channel_count': int(channel_count),
        'loop_count': int(loop_count),
        'recurrence_indicator': recurrence,
        'overlap_persistence': overlap_persistence,
        'max_mean_overlap': float(np.max(overlap_series)) if overlap_series else 0.0,
        'phase_lock_indicator': float(np.mean(phase_alignment_series)) if phase_alignment_series else 0.0,
        'constraint_max': float(np.max(constraint_norms)) if constraint_norms else 0.0,
        'sector_leakage': float(np.max(sector_leakages)) if sector_leakages else 0.0,
        'final_state': np.sum(packet_qs, axis=0),
        'last_env': coarse[4.0]['last_env'],
    }


def augment_with_motif_metrics(measurement: dict[str, Any], motif_nodes: list[int]) -> dict[str, float]:
    motif_idx = np.asarray(sorted(motif_nodes), dtype=int)
    avg_weights = np.asarray(measurement['avg_weights'], dtype=float)
    motif_indicator = np.zeros_like(avg_weights)
    motif_indicator[motif_idx] = 1.0
    motif_node_fraction = float(len(motif_nodes) / max(len(avg_weights), 1))
    motif_mass_fraction = float(np.sum(avg_weights[motif_idx])) if motif_idx.size else 0.0
    if motif_idx.size and np.std(avg_weights) > 1.0e-18 and np.std(motif_indicator) > 1.0e-18:
        motif_corr = float(np.corrcoef(motif_indicator, avg_weights)[0, 1])
    else:
        motif_corr = 0.0
    pinning_score = float(motif_mass_fraction) if motif_idx.size else 0.0
    component_masses = [float(np.sum(np.asarray(comp, dtype=float)[motif_idx])) for comp in measurement['component_avg_weights']] if motif_idx.size else []
    total_component_mass = float(np.sum(component_masses)) if component_masses else 0.0
    if total_component_mass > 1.0e-12:
        routing_asymmetry = float((max(component_masses) - min(component_masses)) / total_component_mass)
    else:
        routing_asymmetry = 0.0
    return {
        'motif_occupancy_correlation': motif_corr,
        'motif_mass_fraction': motif_mass_fraction,
        'pinning_score': pinning_score,
        'routing_asymmetry': routing_asymmetry,
    }


def classify_baseline_topology(metrics: dict[str, float]) -> str:
    if (
        float(metrics['transport_span']) >= 0.12
        and float(metrics['exchange_coherence']) >= 0.28
        and float(metrics['flow_concentration_index']) >= 40.0
        and int(metrics['loop_count']) <= 1
        and float(metrics['recurrence_indicator']) <= 0.1
    ):
        return 'braid_like_exchange'
    if int(metrics['channel_count']) >= 12 or (int(metrics['loop_count']) >= 1 and float(metrics['transport_span']) >= 0.04):
        return 'split_channel_exchange'
    if float(metrics['max_mean_overlap']) >= 0.55 and float(metrics['transport_span']) <= 0.08:
        return 'trapped_local'
    if (
        float(metrics['pinning_score']) >= 0.12
        and float(metrics['motif_occupancy_correlation']) >= 0.08
        and float(metrics['transport_span']) <= 0.06
    ):
        return 'defect_pinned'
    if (
        float(metrics['exchange_coherence']) >= 0.15
        and float(metrics['flow_concentration_index']) >= 25.0
        and float(metrics['recurrence_indicator']) <= 0.35
    ):
        return 'smeared_transfer'
    return 'unresolved'


def classify_run_topology(run_metrics: dict[str, float], baseline_metrics: dict[str, float], params: dict[str, Any]) -> str:
    flow_tol = float(params.get('topology_flow_tolerance', 5.0))
    coh_tol = float(params.get('topology_coherence_tolerance', 0.03))
    if (
        float(run_metrics['transport_span']) >= max(float(baseline_metrics['transport_span']) + 0.08, 0.08)
        and float(run_metrics['exchange_coherence']) >= max(float(baseline_metrics['exchange_coherence']) + 0.04, 0.26)
        and float(run_metrics['flow_concentration_index']) >= max(float(baseline_metrics['flow_concentration_index']) - flow_tol, 35.0)
        and int(run_metrics['loop_count']) <= 1
        and float(run_metrics['recurrence_indicator']) <= 0.1
    ):
        return 'braid_like_exchange'
    if (
        int(run_metrics['channel_count']) >= max(int(baseline_metrics['channel_count']) + 2, 12)
        or (int(run_metrics['loop_count']) >= 1 and float(run_metrics['transport_span']) >= max(float(baseline_metrics['transport_span']) + 0.02, 0.04))
    ):
        return 'split_channel_exchange'
    if (
        float(run_metrics['pinning_score']) >= max(float(baseline_metrics['pinning_score']) + 0.03, 0.12)
        and float(run_metrics['motif_occupancy_correlation']) >= 0.08
        and float(run_metrics['transport_span']) <= max(float(baseline_metrics['transport_span']) + 0.03, 0.08)
    ):
        return 'defect_pinned'
    if (
        float(run_metrics['max_mean_overlap']) >= max(float(baseline_metrics['max_mean_overlap']), 0.55)
        and float(run_metrics['transport_span']) <= 0.08
    ):
        return 'trapped_local'
    if (
        float(run_metrics['exchange_coherence']) >= max(float(baseline_metrics['exchange_coherence']) - coh_tol, 0.15)
        and float(run_metrics['flow_concentration_index']) >= max(0.75 * float(baseline_metrics['flow_concentration_index']), 25.0)
        and float(run_metrics['recurrence_indicator']) <= 0.35
    ):
        return 'smeared_transfer'
    return 'unresolved'


def classify_effect(run_metrics: dict[str, float], baseline_metrics: dict[str, float], params: dict[str, Any]) -> str:
    flow_tol = float(params.get('effect_flow_tolerance', 5.0))
    coh_tol = float(params.get('effect_coherence_tolerance', 0.03))
    route_tol = float(params.get('effect_routing_tolerance', 0.2))
    baseline_topology = str(baseline_metrics['topology_label'])
    run_topology = str(run_metrics['topology_label'])
    if run_topology == 'braid_like_exchange' and baseline_topology != 'braid_like_exchange':
        return 'braid-enhancing'
    if baseline_topology == 'braid_like_exchange' and run_topology != 'braid_like_exchange':
        return 'braid-breaking'
    if run_topology == 'split_channel_exchange' and int(run_metrics['channel_count']) > int(baseline_metrics['channel_count']):
        return 'channel-splitting'
    if run_topology in {'defect_pinned', 'trapped_local'} and (
        float(run_metrics['pinning_score']) > float(baseline_metrics['pinning_score']) + 0.02
        or float(run_metrics['motif_mass_fraction']) > float(baseline_metrics['motif_mass_fraction']) + 0.02
    ):
        return 'local-pinning'
    if run_topology == baseline_topology:
        if (
            abs(float(run_metrics['delta_vs_baseline_concentration'])) <= flow_tol
            and abs(float(run_metrics['delta_vs_baseline_coherence'])) <= coh_tol
            and abs(float(run_metrics['routing_asymmetry']) - float(baseline_metrics['routing_asymmetry'])) <= route_tol
        ):
            return 'neutral / no significant effect'
        return 'topology-preserving'
    return 'neutral / no significant effect'


def motif_overlay_plot(
    path: Path,
    run_id: str,
    midpoints: np.ndarray,
    avg_weights: np.ndarray,
    motif_nodes: list[int],
    seed_nodes: list[int],
    boundary_type: str,
) -> None:
    coords = np.asarray(midpoints, dtype=float)
    slice_value = 0.5 if boundary_type == 'periodic' else 0.5
    distances = np.abs(coords[:, 2] - slice_value)
    threshold = np.min(distances) + 1.0e-9
    mask = distances <= threshold
    motif_mask = np.zeros(len(coords), dtype=bool)
    motif_mask[np.asarray(motif_nodes, dtype=int)] = True
    seed_mask = np.zeros(len(coords), dtype=bool)
    if seed_nodes:
        seed_mask[np.asarray(seed_nodes, dtype=int)] = True

    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    sc = ax.scatter(coords[mask, 0], coords[mask, 1], c=avg_weights[mask], s=16, cmap='viridis', alpha=0.65)
    if np.any(mask & motif_mask):
        ax.scatter(coords[mask & motif_mask, 0], coords[mask & motif_mask, 1], s=44, facecolors='none', edgecolors='tab:red', linewidths=1.4, label='motif nodes')
    if np.any(mask & seed_mask):
        ax.scatter(coords[mask & seed_mask, 0], coords[mask & seed_mask, 1], s=72, marker='x', c='tab:orange', linewidths=1.8, label='seed')
    ax.set_title(f'Motif highlight: {run_id}')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if np.any(mask & motif_mask) or np.any(mask & seed_mask):
        ax.legend(fontsize=8)
    plt.colorbar(sc, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def baseline_key(run: dict[str, Any]) -> tuple[str, int, str]:
    return (str(run['representative_id']), int(run['resolution']), str(run['boundary_type']))


def get_baseline_measurement(
    run: dict[str, Any],
    rep: dict[str, Any],
    runsheet: dict[str, Any],
    defaults: dict[str, Any],
    base: dict[str, Any],
) -> dict[str, Any]:
    key = baseline_key(run)
    if key in BASELINE_CACHE:
        return BASELINE_CACHE[key]

    resolution = int(run['resolution'])
    data, projector = get_cached_setup(resolution, defaults, str(base['boundary_type']))
    target_lambda = estimate_lambda_max(data.L1) if bool(runsheet.get('parameter_policy', {}).get('match_lambda_to_frozen_L1', True)) else 0.0
    baseline_adj = build_line_graph_cache(data, max_shell=int(runsheet.get('parameter_policy', {}).get('max_shell', 3)))['adjacency_lists']
    operator_info = build_operator_from_adjacency(
        data=data,
        adjacency_lists=baseline_adj,
        params=runsheet.get('parameter_policy', {}),
        target_lambda=target_lambda,
    )
    measurement = measure_configuration(rep, data, projector, operator_info, defaults, base, motif_nodes=None)
    payload = {
        'data': data,
        'projector': projector,
        'operator_info': operator_info,
        'measurement': measurement,
    }
    BASELINE_CACHE[key] = payload
    return payload


def simulate_run(
    run: dict[str, Any],
    rep: dict[str, Any],
    motif: dict[str, Any],
    runsheet: dict[str, Any],
    defaults: dict[str, Any],
    base: dict[str, Any],
) -> tuple[dict[str, Any], list[Path]]:
    baseline = get_baseline_measurement(run, rep, runsheet, defaults, base)
    data = baseline['data']
    projector = baseline['projector']
    baseline_operator = baseline['operator_info']
    baseline_measurement = baseline['measurement']
    motif_injection = build_motif_injection(
        data=data,
        representative_id=str(rep['representative_id']),
        motif_id=str(motif['motif_id']),
        baseline_avg_weights=np.asarray(baseline_measurement['avg_weights'], dtype=float),
        params=runsheet.get('parameter_policy', {}),
    )
    target_lambda = estimate_lambda_max(data.L1) if bool(runsheet.get('parameter_policy', {}).get('match_lambda_to_frozen_L1', True)) else 0.0
    motif_operator = build_operator_from_adjacency(
        data=data,
        adjacency_lists=motif_injection['adjacency_lists'],
        params=runsheet.get('parameter_policy', {}),
        target_lambda=target_lambda,
    )

    baseline_motif_metrics = augment_with_motif_metrics(baseline_measurement, motif_injection['motif_nodes'])
    baseline_summary = {
        'flow_concentration_index': float(baseline_measurement['flow_concentration_index']),
        'exchange_coherence': float(baseline_measurement['exchange_coherence']),
        'transport_span': float(baseline_measurement['transport_span']),
        'channel_count': int(baseline_measurement['channel_count']),
        'loop_count': int(baseline_measurement['loop_count']),
        'recurrence_indicator': float(baseline_measurement['recurrence_indicator']),
        'max_mean_overlap': float(baseline_measurement['max_mean_overlap']),
        **baseline_motif_metrics,
    }
    baseline_summary['topology_label'] = classify_baseline_topology(baseline_summary)

    measurement = measure_configuration(
        rep=rep,
        data=data,
        projector=projector,
        operator_info=motif_operator,
        defaults=defaults,
        base=base,
        motif_nodes=motif_injection['motif_nodes'],
    )
    run_motif_metrics = augment_with_motif_metrics(measurement, motif_injection['motif_nodes'])
    summary = {
        'run_id': run['run_id'],
        'representative_id': run['representative_id'],
        'motif_id': motif['motif_id'],
        'motif_label': motif['label'],
        'baseline_topology': baseline_summary['topology_label'],
        'topology_label': '',
        'effect_class': '',
        'braid_like_exchange': 0,
        'smeared_transfer': 0,
        'split_channel_exchange': 0,
        'defect_pinned': 0,
        'trapped_local': 0,
        'unresolved': 0,
        'flow_concentration_index': float(measurement['flow_concentration_index']),
        'exchange_coherence': float(measurement['exchange_coherence']),
        'transport_span': float(measurement['transport_span']),
        'channel_count': int(measurement['channel_count']),
        'loop_count': int(measurement['loop_count']),
        'motif_occupancy_correlation': float(run_motif_metrics['motif_occupancy_correlation']),
        'motif_mass_fraction': float(run_motif_metrics['motif_mass_fraction']),
        'pinning_score': float(run_motif_metrics['pinning_score']),
        'routing_asymmetry': float(run_motif_metrics['routing_asymmetry']),
        'recurrence_indicator': float(measurement['recurrence_indicator']),
        'degree_spectrum_shift': float(np.linalg.norm(np.sort(motif_operator['graph_degrees']) - np.sort(baseline_operator['graph_degrees'])) / math.sqrt(max(len(motif_operator['graph_degrees']), 1))),
        'clustering_coefficient_shift': float(motif_operator['mean_clustering'] - baseline_operator['mean_clustering']),
        'delta_vs_baseline_concentration': float(measurement['flow_concentration_index']) - float(baseline_summary['flow_concentration_index']),
        'delta_vs_baseline_coherence': float(measurement['exchange_coherence']) - float(baseline_summary['exchange_coherence']),
        'delta_vs_baseline_channel_count': int(measurement['channel_count']) - int(baseline_summary['channel_count']),
        'notes': run['notes'],
        'interaction_label': str(measurement['interaction_label']),
        'coarse_label': str(measurement['coarse_label']),
        'max_mean_overlap': float(measurement['max_mean_overlap']),
    }
    topology = classify_run_topology(summary, baseline_summary, runsheet.get('parameter_policy', {}))
    summary['topology_label'] = topology
    summary['braid_like_exchange'] = int(topology == 'braid_like_exchange')
    summary['smeared_transfer'] = int(topology == 'smeared_transfer')
    summary['split_channel_exchange'] = int(topology == 'split_channel_exchange')
    summary['defect_pinned'] = int(topology == 'defect_pinned')
    summary['trapped_local'] = int(topology == 'trapped_local')
    summary['unresolved'] = int(topology == 'unresolved')
    summary['effect_class'] = classify_effect(summary, baseline_summary, runsheet.get('parameter_policy', {}))

    plot_paths: list[Path] = []
    trace_path = WORK_PLOT_DIR / f"{run['run_id']}_coherence_span_trace.png"
    centers_arr = np.asarray(measurement['total_centers'], dtype=float)
    center_shift_trace = [
        float(np.linalg.norm(displacement(np.asarray([center], dtype=float), centers_arr[0], str(base['boundary_type']))[0]))
        for center in centers_arr
    ]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))
    axes[0].plot(measurement['times'], measurement['overlap_series'], color='tab:red', label='overlap')
    axes[0].plot(measurement['times'], measurement['phase_alignment_series'], color='tab:purple', label='phase')
    axes[0].set_title('Exchange coherence trace')
    axes[0].set_xlabel('time')
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=8)
    axes[1].plot(measurement['times'], center_shift_trace, color='tab:orange')
    axes[1].set_title('Transport span trace')
    axes[1].set_xlabel('time')
    axes[1].grid(alpha=0.25)
    fig.suptitle(run['run_id'])
    fig.savefig(trace_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(trace_path)

    motif_trace_path = WORK_PLOT_DIR / f"{run['run_id']}_motif_occupancy_trace.png"
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(measurement['times'], measurement['motif_mass_trace'], color='tab:blue', label='total motif occupancy')
    for idx, series in enumerate(measurement['component_motif_mass_trace']):
        ax.plot(measurement['times'], series, linewidth=1.0, alpha=0.7, label=f'component {idx + 1}')
    ax.set_title('Motif occupancy trace')
    ax.set_xlabel('time')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(motif_trace_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(motif_trace_path)

    overlay_path = WORK_PLOT_DIR / f"{run['run_id']}_motif_graph_overlay.png"
    motif_overlay_plot(
        path=overlay_path,
        run_id=run['run_id'],
        midpoints=data.midpoints,
        avg_weights=np.asarray(measurement['avg_weights'], dtype=float),
        motif_nodes=motif_injection['motif_nodes'],
        seed_nodes=motif_injection['seed_nodes'],
        boundary_type=str(base['boundary_type']),
    )
    plot_paths.append(overlay_path)

    field_path = WORK_PLOT_DIR / f"{run['run_id']}_field_snapshot.png"
    create_field_snapshot(field_path, run['run_id'], data.midpoints, measurement['final_state'], str(base['boundary_type']), int(defaults['slice_axis']))
    plot_paths.append(field_path)

    env = measurement['last_env']
    if env is not None:
        env_path = WORK_PLOT_DIR / f"{run['run_id']}_coarse_envelope.png"
        fig, ax = plt.subplots(figsize=(5.4, 4.6))
        render_grid_slice(ax, env, f"Envelope sigma=4: {run['run_id']}")
        fig.savefig(env_path, dpi=180, bbox_inches='tight')
        plt.close(fig)
        plot_paths.append(env_path)

    result = {
        'run': run,
        'motif': motif,
        'baseline_summary': baseline_summary,
        'motif_injection': {
            'motif_nodes': motif_injection['motif_nodes'],
            'seed_nodes': motif_injection['seed_nodes'],
            'patch_nodes': motif_injection['patch_nodes'],
            'desired_edges': motif_injection['desired_edges'],
            'added_edges': motif_injection['added_edges'],
            'removed_edges': motif_injection['removed_edges'],
        },
        'operator_summary': {
            'lambda_max': float(motif_operator['lambda_max']),
            'raw_lambda_max': float(motif_operator['raw_lambda_max']),
            'scale_factor': float(motif_operator['scale_factor']),
            'support_entries': int(motif_operator['support_entries']),
            'degree_mean': float(motif_operator['degree_mean']),
            'degree_std': float(motif_operator['degree_std']),
            'mean_clustering': float(motif_operator['mean_clustering']),
        },
        'summary': summary,
        'times': measurement['times'],
        'motif_mass_trace': measurement['motif_mass_trace'],
        'pair_distance_series': measurement['pair_distance_series'],
        'overlap_series': measurement['overlap_series'],
        'phase_alignment_series': measurement['phase_alignment_series'],
    }
    return result, plot_paths


def create_summary_plots(rows: list[dict[str, Any]]) -> list[Path]:
    plot_paths: list[Path] = []
    motifs = sorted({str(row['motif_id']) for row in rows}, key=motif_order)
    reps = sorted({str(row['representative_id']) for row in rows})
    row_map = {(str(row['representative_id']), str(row['motif_id'])): row for row in rows}
    topology_scores = {
        'unresolved': 0,
        'smeared_transfer': 1,
        'trapped_local': 2,
        'defect_pinned': 3,
        'split_channel_exchange': 4,
        'braid_like_exchange': 5,
    }

    matrix_path = WORK_PLOT_DIR / 'stage_c0_4_motif_encounter_topology_matrix.png'
    data = np.zeros((len(reps), len(motifs)), dtype=float)
    for i, rep in enumerate(reps):
        for j, motif_id in enumerate(motifs):
            data[i, j] = float(topology_scores.get(str(row_map[(rep, motif_id)]['topology_label']), 0))
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    im = ax.imshow(data, cmap='viridis', vmin=0, vmax=max(topology_scores.values()))
    ax.set_xticks(range(len(motifs)))
    ax.set_xticklabels(motifs, rotation=20, ha='right')
    ax.set_yticks(range(len(reps)))
    ax.set_yticklabels(reps)
    ax.set_title('Motif x encounter topology matrix')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(matrix_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(matrix_path)

    effect_path = WORK_PLOT_DIR / 'stage_c0_4_braid_effect_panel.png'
    effect_classes = [
        'braid-enhancing',
        'braid-breaking',
        'channel-splitting',
        'local-pinning',
        'topology-preserving',
        'neutral / no significant effect',
    ]
    x = np.arange(len(motifs))
    width = 0.12
    fig, ax = plt.subplots(figsize=(11.2, 4.8))
    for idx, effect in enumerate(effect_classes):
        values = [sum(1 for row in rows if row['motif_id'] == motif_id and row['effect_class'] == effect) for motif_id in motifs]
        ax.bar(x + (idx - 0.5 * (len(effect_classes) - 1)) * width, values, width=width, label=effect)
    ax.set_xticks(x)
    ax.set_xticklabels(motifs, rotation=20, ha='right')
    ax.set_title('Braid-enhancement / motif-effect panel')
    ax.grid(alpha=0.25, axis='y')
    ax.legend(fontsize=8)
    fig.savefig(effect_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(effect_path)

    occupancy_path = WORK_PLOT_DIR / 'stage_c0_4_motif_occupancy_vs_topology.png'
    color_map = {
        'unresolved': 'tab:gray',
        'smeared_transfer': 'tab:blue',
        'trapped_local': 'tab:green',
        'defect_pinned': 'tab:red',
        'split_channel_exchange': 'tab:orange',
        'braid_like_exchange': 'tab:purple',
    }
    fig, ax = plt.subplots(figsize=(6.8, 5.2))
    for row in rows:
        ax.scatter(
            float(row['motif_occupancy_correlation']),
            float(row['pinning_score']),
            color=color_map.get(str(row['topology_label']), 'black'),
            s=48,
        )
        ax.text(
            float(row['motif_occupancy_correlation']) + 0.005,
            float(row['pinning_score']) + 0.02,
            f"{row['motif_id'].split('_')[0]}:{row['representative_id'].split('_')[0]}",
            fontsize=7,
        )
    ax.set_title('Motif occupancy vs topology summary')
    ax.set_xlabel('motif occupancy correlation')
    ax.set_ylabel('pinning score')
    ax.grid(alpha=0.25)
    fig.savefig(occupancy_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(occupancy_path)
    return plot_paths


def write_note(path: Path, json_rel: str, csv_rel: str, rows: list[dict[str, Any]], runsheet: dict[str, Any]) -> None:
    lines = [
        '# Stage C0.4 Controlled Motif Injection v1',
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
    for row in sorted(rows, key=lambda item: (motif_order(str(item['motif_id'])), str(item['representative_id']))):
        lines.extend([
            f"- `{row['motif_label']}` / `{row['representative_id']}`",
            f"  - matched baseline topology: `{row['baseline_topology']}`",
            f"  - motif run topology: `{row['topology_label']}`",
            f"  - effect class: `{row['effect_class']}`",
            f"  - flow concentration: `{row['flow_concentration_index']:.4f}`",
            f"  - exchange coherence: `{row['exchange_coherence']:.4f}`",
            f"  - motif occupancy correlation: `{row['motif_occupancy_correlation']:.4f}`",
            f"  - pinning score: `{row['pinning_score']:.4f}`",
        ])
    lines.extend([
        '',
        'Interpretation boundary:',
        '- this scan asks whether local combinatorial motifs seed reproducible exchange-topology changes in the no-distance branch',
        '- it does not claim emergent particles, topological matter, confinement, or true bound states',
        '',
    ])
    path.write_text('\n'.join(lines), encoding='utf-8')


def summarize_outcome(rows: list[dict[str, Any]]) -> tuple[str, str]:
    effect_counts = Counter(str(row['effect_class']) for row in rows)
    motif_effects: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        motif_effects[str(row['motif_id'])][str(row['effect_class'])] += 1

    reproducible_positive = any(
        counts['braid-enhancing'] >= 1 or counts['channel-splitting'] >= 2 or counts['local-pinning'] >= 2
        for counts in motif_effects.values()
    )
    all_neutral = all(
        str(row['effect_class']) in {'neutral / no significant effect', 'topology-preserving'}
        for row in rows
    )

    if reproducible_positive:
        observation = 'at least one motif class produces a repeated non-neutral topology shift across the encounter family scan'
        conclusion = 'C0.4 is positive: specific local combinatorial motifs act as topology-seeding or topology-redirecting structures in the current no-distance branch'
    elif all_neutral:
        observation = 'motif injection does not systematically change topology labels and the measured shifts stay inside the matched-baseline envelope'
        conclusion = 'C0.4 is negative on this first pass: local combinatorial motifs are not the active topology-seeding mechanism here'
    else:
        observation = 'motif injection causes isolated quantitative shifts, but no motif class produces a repeated new topology across the 3 encounter families'
        conclusion = 'C0.4 is weakly exploratory rather than positive: motif effects appear local and case-specific, not yet reproducible as a seeded topology class'
    return observation, conclusion


def main() -> None:
    args = parse_args()
    defaults = load_stage10_defaults()
    runsheet = load_runsheet(args.runsheet)
    base = runsheet['base_seed_reference']
    runs = selected_runs(runsheet['runs'], args.run_ids)
    if not runs:
        raise SystemExit('no Stage C0.4 runs selected')

    rep_lookup = lookup_by_id(runsheet['representatives'], 'representative_id')
    motif_lookup = lookup_by_id(runsheet['motifs'], 'motif_id')
    note_name = str(runsheet.get('note_name', 'Stage_C0_4_Controlled_Motif_Injection_v1.md'))
    note_path = ATLAS_NOTES / note_name

    rows: list[dict[str, Any]] = []
    payload_runs: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    for run in sorted(runs, key=lambda item: (motif_order(str(item['motif_id'])), str(item['representative_id']))):
        rep = rep_lookup[run['representative_id']]
        motif = motif_lookup[run['motif_id']]
        payload, run_plots = simulate_run(run, rep, motif, runsheet, defaults, base)
        rows.append(payload['summary'])
        payload_runs.append(payload)
        plot_paths.extend(run_plots)

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
        'topology_counts': dict(Counter(row['topology_label'] for row in rows)),
        'effect_class_counts': dict(Counter(row['effect_class'] for row in rows)),
        'observation': observation,
        'conclusion': conclusion,
    }

    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug=str(runsheet.get('experiment_slug', 'stage_c0_4_controlled_motif_injection')),
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
