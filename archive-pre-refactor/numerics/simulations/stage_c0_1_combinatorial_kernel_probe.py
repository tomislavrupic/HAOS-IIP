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
    build_transverse_setup,
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

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_1_combinatorial_kernel_runs.json'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage_c0_1_combinatorial_kernel'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

SETUP_CACHE: dict[tuple[int, str, float, float], tuple[Any, Any]] = {}
LINE_GRAPH_CACHE: dict[tuple[int, str, int], dict[str, Any]] = {}
OPERATOR_CACHE: dict[tuple[int, str, str, float, float, float, bool], dict[str, Any]] = {}

ANALYSIS_STRIDE = 4

CSV_FIELDS = [
    'run_id',
    'representative_id',
    'kernel_class_id',
    'kernel_label',
    'resolution',
    'interaction_label',
    'coarse_label',
    'composite_lifetime',
    'binding_persistence_score',
    'coarse_persistence_sigma4',
    'transport_span',
    'corridor_dwell_time',
    'initial_mean_pair_distance',
    'min_mean_pair_distance',
    'final_mean_pair_distance',
    'max_mean_overlap',
    'phase_lock_indicator',
    'constraint_max',
    'sector_leakage',
    'operator_lambda_max',
    'operator_scale_factor',
    'delta_vs_gaussian_composite_lifetime',
    'delta_vs_gaussian_binding_persistence',
    'delta_vs_gaussian_coarse_persistence',
    'interaction_changed_vs_gaussian',
    'coarse_changed_vs_gaussian',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.1 combinatorial kernel austerity probe.')
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


def kernel_order(kernel_class_id: str) -> int:
    return {
        'gaussian_midpoint': 0,
        'pure_adjacency': 1,
        'graph_shell': 2,
    }.get(kernel_class_id, 99)


def get_cached_setup(n_side: int, defaults: dict[str, Any], boundary_type: str) -> tuple[Any, Any]:
    key = (int(n_side), str(boundary_type), float(defaults['epsilon']), float(defaults['harmonic_tol']))
    if key not in SETUP_CACHE:
        SETUP_CACHE[key] = build_transverse_setup(
            n_side=int(n_side),
            epsilon=float(defaults['epsilon']),
            harmonic_tol=float(defaults['harmonic_tol']),
            boundary_type=str(boundary_type),
        )
    return SETUP_CACHE[key]


def build_line_graph_cache(data: Any, max_shell: int) -> dict[str, Any]:
    key = (int(data.n_side), str(data.boundary_type), int(max_shell))
    if key in LINE_GRAPH_CACHE:
        return LINE_GRAPH_CACHE[key]

    node_count = int(data.points.shape[0])
    incident: list[list[int]] = [[] for _ in range(node_count)]
    for edge_idx, (u, v) in enumerate(data.edges):
        incident[int(u)].append(edge_idx)
        incident[int(v)].append(edge_idx)

    adjacency_lists: list[list[int]] = []
    degrees = np.zeros(len(data.edges), dtype=float)
    for edge_idx, (u, v) in enumerate(data.edges):
        nbrs = set(incident[int(u)])
        nbrs.update(incident[int(v)])
        nbrs.discard(edge_idx)
        ordered = sorted(nbrs)
        adjacency_lists.append(ordered)
        degrees[edge_idx] = float(len(ordered))

    shells: list[dict[int, int]] = []
    for start in range(len(data.edges)):
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

    payload = {
        'adjacency_lists': adjacency_lists,
        'degrees': degrees,
        'shells': shells,
    }
    LINE_GRAPH_CACHE[key] = payload
    return payload


def build_edge_graph_operator(
    data: Any,
    kernel_class_id: str,
    params: dict[str, Any],
    target_lambda: float,
) -> dict[str, Any]:
    sigma_multiple = float(params.get('gaussian_sigma_multiple_h', 1.5))
    max_shell = int(params.get('max_shell', 3))
    graph_shell_weights = {
        int(key): float(value)
        for key, value in dict(params.get('graph_shell_weights', {'1': 1.0, '2': 0.5, '3': 0.25})).items()
    }
    degree_normalized = bool(params.get('adjacency_degree_normalized', True))
    key = (
        int(data.n_side),
        str(data.boundary_type),
        str(kernel_class_id),
        sigma_multiple,
        float(graph_shell_weights.get(1, 1.0)),
        float(graph_shell_weights.get(2, 0.5)),
        degree_normalized,
    )
    if key in OPERATOR_CACHE:
        return OPERATOR_CACHE[key]

    line_graph = build_line_graph_cache(data, max_shell=max_shell)
    shells: list[dict[int, int]] = line_graph['shells']
    degrees = np.asarray(line_graph['degrees'], dtype=float)
    h = 1.0 / float(data.n_side)
    sigma = sigma_multiple * h

    rows: list[int] = []
    cols: list[int] = []
    values: list[float] = []

    for edge_idx, shell_map in enumerate(shells):
        for nbr_idx, shell in shell_map.items():
            if nbr_idx <= edge_idx:
                continue
            weight = 0.0
            if kernel_class_id == 'gaussian_midpoint':
                delta = displacement(np.asarray([data.midpoints[nbr_idx]], dtype=float), data.midpoints[edge_idx], str(data.boundary_type))[0]
                r2 = float(np.dot(delta, delta))
                weight = math.exp(-r2 / max(2.0 * sigma * sigma, 1.0e-12))
            elif kernel_class_id == 'pure_adjacency':
                if shell == 1:
                    if degree_normalized:
                        weight = 1.0 / math.sqrt(max(degrees[edge_idx] * degrees[nbr_idx], 1.0e-12))
                    else:
                        weight = 1.0
            elif kernel_class_id == 'graph_shell':
                weight = float(graph_shell_weights.get(shell, 0.0))
            else:
                raise ValueError(f'unsupported kernel class: {kernel_class_id}')

            if weight <= 0.0:
                continue
            rows.extend([edge_idx, nbr_idx])
            cols.extend([nbr_idx, edge_idx])
            values.extend([weight, weight])

    size = len(data.edges)
    weight_matrix = sp.coo_matrix((values, (rows, cols)), shape=(size, size), dtype=float).tocsr()
    degree_vec = np.asarray(weight_matrix.sum(axis=1)).ravel()
    operator = (sp.diags(degree_vec) - weight_matrix).tocsr()

    raw_lambda = estimate_lambda_max(operator)
    if target_lambda > 0.0 and raw_lambda > 0.0:
        scale_factor = float(target_lambda / raw_lambda)
        operator = (scale_factor * operator).tocsr()
        lambda_max = estimate_lambda_max(operator)
    else:
        scale_factor = 1.0
        lambda_max = raw_lambda

    payload = {
        'operator': operator,
        'weight_matrix': weight_matrix,
        'degree_vec': degree_vec,
        'raw_lambda_max': raw_lambda,
        'lambda_max': lambda_max,
        'scale_factor': scale_factor,
        'support_entries': int(weight_matrix.nnz),
        'sigma': sigma,
        'max_shell': max_shell,
    }
    OPERATOR_CACHE[key] = payload
    return payload


def pairwise_distance_stats(centers: list[np.ndarray], boundary_type: str) -> tuple[float, float]:
    values: list[float] = []
    for i in range(len(centers)):
        for j in range(i + 1, len(centers)):
            delta = displacement(np.asarray([centers[j]], dtype=float), np.asarray(centers[i], dtype=float), boundary_type)[0]
            values.append(float(np.linalg.norm(delta)))
    if not values:
        return 0.0, 0.0
    arr = np.asarray(values, dtype=float)
    return float(np.mean(arr)), float(np.min(arr))


def pairwise_overlap_mean(component_qs: list[np.ndarray]) -> float:
    values: list[float] = []
    weights = [np.abs(q) ** 2 for q in component_qs]
    for i in range(len(weights)):
        for j in range(i + 1, len(weights)):
            wi = weights[i]
            wj = weights[j]
            denom = math.sqrt(max(float(np.sum(wi)) * float(np.sum(wj)), 1.0e-12))
            values.append(float(np.sum(np.sqrt(wi * wj)) / denom))
    return float(np.mean(values)) if values else 0.0


def pairwise_phase_alignment(component_qs: list[np.ndarray]) -> float:
    values: list[float] = []
    for i in range(len(component_qs)):
        for j in range(i + 1, len(component_qs)):
            qi = component_qs[i]
            qj = component_qs[j]
            denom = max(float(np.linalg.norm(qi)) * float(np.linalg.norm(qj)), 1.0e-12)
            values.append(abs(float(np.vdot(qi, qj).real / denom)))
    return float(np.mean(values)) if values else 0.0


def path_length(centers: np.ndarray, boundary_type: str) -> float:
    if len(centers) < 2:
        return 0.0
    total = 0.0
    for prev, cur in zip(centers[:-1], centers[1:]):
        total += float(np.linalg.norm(displacement(cur[None, :], prev, boundary_type)[0]))
    return total


def render_grid_slice(ax, grid: np.ndarray, title: str, cmap: str = 'viridis') -> None:
    slice_idx = grid.shape[2] // 2
    img = grid[:, :, slice_idx].T
    im = ax.imshow(img, origin='lower', cmap=cmap)
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)


def simulate_run(
    run: dict[str, Any],
    rep: dict[str, Any],
    kernel_class: dict[str, Any],
    runsheet: dict[str, Any],
    defaults: dict[str, Any],
    base: dict[str, Any],
) -> tuple[dict[str, Any], list[Path]]:
    resolution = int(run['resolution'])
    boundary_type = str(base['boundary_type'])
    data, projector = get_cached_setup(resolution, defaults, boundary_type)
    target_lambda = estimate_lambda_max(data.L1) if bool(runsheet.get('parameter_policy', {}).get('match_lambda_to_frozen_L1', True)) else 0.0
    operator_info = build_edge_graph_operator(
        data=data,
        kernel_class_id=str(run['kernel_class_id']),
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
        total_q = np.sum(packet_qs, axis=0)
        total_v = np.sum(packet_vs, axis=0)
        total_weights = np.abs(total_q) ** 2
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

    binding_mask = (np.asarray(pair_distance_series) <= 0.25) & (np.asarray(overlap_series) >= 0.35)
    binding_persistence = float(np.mean(binding_mask.astype(float))) if binding_mask.size else 0.0
    initial_distance = float(pair_distance_series[0]) if pair_distance_series else 0.0
    corridor_mask = np.asarray(pair_distance_series) <= 1.05 * max(initial_distance, 1.0e-12)
    corridor_dwell_time = float(np.sum(corridor_mask.astype(float)) * sample_dt) if corridor_mask.size else 0.0
    transport_span = path_length(centers_arr, boundary_type)

    row = {
        'run_id': run['run_id'],
        'representative_id': run['representative_id'],
        'kernel_class_id': run['kernel_class_id'],
        'kernel_label': kernel_class['label'],
        'resolution': resolution,
        'interaction_label': interaction_label,
        'coarse_label': coarse_label,
        'composite_lifetime': float(interaction_summary['composite_lifetime']),
        'binding_persistence_score': binding_persistence,
        'coarse_persistence_sigma4': float(coarse_summary['coarse_persistence_sigma4']),
        'transport_span': transport_span,
        'corridor_dwell_time': corridor_dwell_time,
        'initial_mean_pair_distance': initial_distance,
        'min_mean_pair_distance': float(np.min(pair_distance_series)) if pair_distance_series else 0.0,
        'final_mean_pair_distance': float(pair_distance_series[-1]) if pair_distance_series else 0.0,
        'max_mean_overlap': float(np.max(overlap_series)) if overlap_series else 0.0,
        'phase_lock_indicator': float(np.mean(phase_alignment_series)) if phase_alignment_series else 0.0,
        'constraint_max': float(np.max(constraint_norms)) if constraint_norms else 0.0,
        'sector_leakage': float(np.max(sector_leakages)) if sector_leakages else 0.0,
        'operator_lambda_max': float(operator_info['lambda_max']),
        'operator_scale_factor': float(operator_info['scale_factor']),
        'delta_vs_gaussian_composite_lifetime': 0.0,
        'delta_vs_gaussian_binding_persistence': 0.0,
        'delta_vs_gaussian_coarse_persistence': 0.0,
        'interaction_changed_vs_gaussian': 0,
        'coarse_changed_vs_gaussian': 0,
        'notes': run['notes'],
    }

    plot_paths: list[Path] = []
    trace_path = WORK_PLOT_DIR / f"{run['run_id']}_pair_distance_overlap.png"
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))
    axes[0].plot(times, pair_distance_series, color='tab:blue')
    axes[0].set_title('Mean pair distance')
    axes[0].set_xlabel('time')
    axes[0].grid(alpha=0.25)
    axes[1].plot(times, overlap_series, color='tab:orange', label='overlap')
    axes[1].plot(times, phase_alignment_series, color='tab:green', label='phase')
    axes[1].set_title('Overlap and phase')
    axes[1].set_xlabel('time')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.suptitle(run['run_id'])
    fig.savefig(trace_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(trace_path)

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
        'kernel_class': kernel_class,
        'operator_summary': {
            'lambda_max': float(operator_info['lambda_max']),
            'raw_lambda_max': float(operator_info['raw_lambda_max']),
            'scale_factor': float(operator_info['scale_factor']),
            'support_entries': int(operator_info['support_entries']),
            'sigma': float(operator_info['sigma']),
            'max_shell': int(operator_info['max_shell']),
        },
        'summary': row,
        'times': times,
        'pair_distance_series': pair_distance_series,
        'overlap_series': overlap_series,
        'phase_alignment_series': phase_alignment_series,
        'coarse_summary': coarse_summary,
    }
    return result, plot_paths


def compare_to_gaussian(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    gaussian_map = {
        row['representative_id']: row
        for row in rows
        if row['kernel_class_id'] == 'gaussian_midpoint'
    }
    for row in rows:
        baseline = gaussian_map.get(row['representative_id'])
        if baseline is None:
            continue
        row['delta_vs_gaussian_composite_lifetime'] = float(row['composite_lifetime']) - float(baseline['composite_lifetime'])
        row['delta_vs_gaussian_binding_persistence'] = float(row['binding_persistence_score']) - float(baseline['binding_persistence_score'])
        row['delta_vs_gaussian_coarse_persistence'] = float(row['coarse_persistence_sigma4']) - float(baseline['coarse_persistence_sigma4'])
        row['interaction_changed_vs_gaussian'] = int(str(row['interaction_label']) != str(baseline['interaction_label']))
        row['coarse_changed_vs_gaussian'] = int(str(row['coarse_label']) != str(baseline['coarse_label']))
    return rows


def create_summary_plots(rows: list[dict[str, Any]]) -> list[Path]:
    plot_paths: list[Path] = []
    reps = sorted({str(row['representative_id']) for row in rows})
    kernels = sorted({str(row['kernel_class_id']) for row in rows}, key=kernel_order)
    kernel_labels = {
        'gaussian_midpoint': 'Gaussian',
        'pure_adjacency': 'Adjacency',
        'graph_shell': 'Shell',
    }
    row_map = {(str(row['representative_id']), str(row['kernel_class_id'])): row for row in rows}
    if not reps or not kernels:
        return plot_paths
    x = np.arange(len(reps))
    width = 0.22

    metric_plot = WORK_PLOT_DIR / 'stage_c0_1_metric_comparison.png'
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.8))
    metric_keys = [
        ('composite_lifetime', 'Composite lifetime'),
        ('binding_persistence_score', 'Binding persistence'),
        ('coarse_persistence_sigma4', 'Coarse persistence (sigma=4)'),
    ]
    for axis, (metric_key, title) in zip(axes, metric_keys):
        offsets = np.arange(len(kernels), dtype=float) - 0.5 * (len(kernels) - 1)
        for idx, kernel in enumerate(kernels):
            values = []
            for rep in reps:
                row = row_map.get((rep, kernel))
                values.append(float(row[metric_key]) if row is not None else 0.0)
            axis.bar(x + offsets[idx] * width, values, width=width, label=kernel_labels.get(kernel, kernel))
        axis.set_title(title)
        axis.set_xticks(x)
        axis.set_xticklabels(reps, rotation=20, ha='right')
        axis.grid(alpha=0.25, axis='y')
        axis.legend(fontsize=8)
    fig.suptitle('Stage C0.1 kernel-class comparison')
    fig.savefig(metric_plot, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(metric_plot)

    delta_plot = WORK_PLOT_DIR / 'stage_c0_1_delta_vs_gaussian.png'
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8))
    non_gaussian = [kernel for kernel in kernels if kernel != 'gaussian_midpoint']
    offsets = np.arange(len(non_gaussian), dtype=float) - 0.5 * (len(non_gaussian) - 1) if non_gaussian else np.asarray([], dtype=float)
    for idx, kernel in enumerate(non_gaussian):
        comp_values = []
        bind_values = []
        for rep in reps:
            row = row_map.get((rep, kernel))
            comp_values.append(float(row['delta_vs_gaussian_composite_lifetime']) if row is not None else 0.0)
            bind_values.append(float(row['delta_vs_gaussian_binding_persistence']) if row is not None else 0.0)
        axes[0].bar(x + offsets[idx] * width, comp_values, width=width, label=kernel_labels.get(kernel, kernel))
        axes[1].bar(x + offsets[idx] * width, bind_values, width=width, label=kernel_labels.get(kernel, kernel))
    axes[0].axhline(0.0, color='k', linewidth=1.0)
    axes[1].axhline(0.0, color='k', linewidth=1.0)
    axes[0].set_title('Composite lifetime delta vs Gaussian')
    axes[1].set_title('Binding persistence delta vs Gaussian')
    for axis in axes:
        axis.set_xticks(x)
        axis.set_xticklabels(reps, rotation=20, ha='right')
        axis.grid(alpha=0.25, axis='y')
        if non_gaussian:
            axis.legend(fontsize=8)
    fig.savefig(delta_plot, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(delta_plot)

    response_plot = WORK_PLOT_DIR / 'stage_c0_1_response_matrix.png'
    score = {
        'dispersive independent regime': 0,
        'transient binding regime': 1,
        'metastable composite regime': 2,
        'oscillatory exchange regime': 3,
        'large-scale drift field regime': 4,
    }
    data = np.zeros((len(reps), len(kernels)), dtype=float)
    for i, rep in enumerate(reps):
        for j, kernel in enumerate(kernels):
            row = row_map.get((rep, kernel))
            data[i, j] = float(score.get(str(row['interaction_label']), 0)) if row is not None else 0.0
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    im = ax.imshow(data, cmap='viridis', vmin=0, vmax=4)
    ax.set_xticks(range(len(kernels)))
    ax.set_xticklabels([kernel_labels.get(k, k) for k in kernels])
    ax.set_yticks(range(len(reps)))
    ax.set_yticklabels(reps)
    ax.set_title('Interaction-regime response matrix')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(response_plot, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(response_plot)

    return plot_paths


def write_note(
    path: Path,
    json_rel: str,
    csv_rel: str,
    rows: list[dict[str, Any]],
    runsheet: dict[str, Any],
) -> None:
    lines = [
        '# Stage C0.1 Combinatorial Kernel Austerity Probe v1',
        '',
        f'Timestamped summary: `{json_rel}`',
        f'Timestamped run table: `{csv_rel}`',
        '',
        f'Architecture notice: {runsheet["architecture_notice"]}',
        '',
        'Kernel comparison table:',
        '',
    ]
    for row in sorted(rows, key=lambda item: (item['representative_id'], kernel_order(str(item['kernel_class_id'])))):
        lines.extend([
            f"- `{row['representative_id']}` / `{row['kernel_label']}`",
            f"  - interaction label: `{row['interaction_label']}`",
            f"  - coarse label: `{row['coarse_label']}`",
            f"  - composite lifetime: `{row['composite_lifetime']:.4f}`",
            f"  - binding persistence: `{row['binding_persistence_score']:.4f}`",
            f"  - coarse persistence sigma=4: `{row['coarse_persistence_sigma4']:.4f}`",
            f"  - delta vs Gaussian composite lifetime: `{row['delta_vs_gaussian_composite_lifetime']:.4f}`",
            f"  - delta vs Gaussian binding persistence: `{row['delta_vs_gaussian_binding_persistence']:.4f}`",
        ])
    lines.extend([
        '',
        'Interpretation boundary:',
        '- this is an explicit edge-graph control branch, not a replacement of the frozen Hodge operator stack',
        '- on the frozen regular-lattice branch, nearest-neighbor edge weights are already uniform, so adjacency-only is not a meaningful falsification unless a separate edge-graph operator is instantiated',
        '- positive outcomes here mean graph-structural survival inside this control branch, not a blanket claim about every existing stage',
        '',
    ])
    path.write_text('\n'.join(lines), encoding='utf-8')


def summarize_outcome(rows: list[dict[str, Any]]) -> tuple[str, str]:
    gaussian_map = {
        str(row['representative_id']): row
        for row in rows
        if row['kernel_class_id'] == 'gaussian_midpoint'
    }
    non_gaussian = [row for row in rows if row['kernel_class_id'] != 'gaussian_midpoint']
    all_survive = bool(non_gaussian) and all(
        gaussian_map.get(str(row['representative_id'])) is not None
        and int(row['interaction_changed_vs_gaussian']) == 0
        and int(row['coarse_changed_vs_gaussian']) == 0
        and abs(float(row['delta_vs_gaussian_composite_lifetime'])) <= 1.0e-9
        and abs(float(row['delta_vs_gaussian_binding_persistence'])) <= 1.0e-9
        and abs(float(row['delta_vs_gaussian_coarse_persistence'])) <= 1.0e-9
        for row in non_gaussian
    )
    survives = any(
        gaussian_map.get(str(row['representative_id'])) is not None
        and float(row['composite_lifetime']) >= 0.9 * float(gaussian_map[str(row['representative_id'])]['composite_lifetime'])
        and int(row['interaction_changed_vs_gaussian']) == 0
        for row in non_gaussian
    )
    shell_only = any(
        row['kernel_class_id'] == 'graph_shell'
        and int(row['interaction_changed_vs_gaussian']) == 0
        for row in non_gaussian
    )
    if all_survive:
        observation = 'all non-Gaussian kernels preserve the Gaussian interaction and coarse-field classes on the projected-transverse edge-graph control branch'
        conclusion = 'the primary persistence observables in this explicit control branch survive removal of coordinate-distance weighting, although secondary transport texture can still shift'
    elif survives:
        observation = 'at least one non-Gaussian kernel preserves the Gaussian interaction class on the projected-transverse edge-graph control branch'
        conclusion = 'some of the observed packet morphology survives removal of coordinate-distance weighting inside the explicit edge-graph control branch'
    elif shell_only:
        observation = 'pure adjacency changes the response class, but graph-shell support retains at least one Gaussian-class morphology'
        conclusion = 'nearest-neighbor austerity is too severe here, while discrete shell structure still carries part of the morphology without coordinate weighting'
    else:
        observation = 'all non-Gaussian kernels change the Gaussian response class on this control branch'
        conclusion = 'within the explicit edge-graph control branch, the Gaussian metric weighting is doing material structural work'
    return observation, conclusion


def main() -> None:
    args = parse_args()
    defaults = load_stage10_defaults()
    runsheet = load_runsheet(args.runsheet)
    base = runsheet['base_seed_reference']
    runs = selected_runs(runsheet['runs'], args.run_ids)
    if not runs:
        raise SystemExit('no Stage C0.1 runs selected')

    rep_lookup = lookup_by_id(runsheet['representatives'], 'representative_id')
    kernel_lookup = lookup_by_id(runsheet['kernel_classes'], 'kernel_class_id')
    note_name = str(runsheet.get('note_name', 'Stage_C0_1_Combinatorial_Kernel_Austerity_Probe_v1.md'))
    note_path = ATLAS_NOTES / note_name

    rows: list[dict[str, Any]] = []
    payload_runs: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    for run in sorted(runs, key=lambda item: (rep_lookup[item['representative_id']]['representative_id'], kernel_order(str(item['kernel_class_id'])))):
        rep = rep_lookup[run['representative_id']]
        kernel_class = kernel_lookup[run['kernel_class_id']]
        payload, run_plots = simulate_run(run, rep, kernel_class, runsheet, defaults, base)
        payload_runs.append(payload)
        rows.append(payload['summary'])
        plot_paths.extend(run_plots)

    rows = compare_to_gaussian(rows)
    plot_paths.extend(create_summary_plots(rows))
    observation, conclusion = summarize_outcome(rows)

    result = {
        'stage': runsheet['stage'],
        'description': runsheet['description'],
        'architecture_notice': runsheet['architecture_notice'],
        'base_seed_reference': base,
        'parameter_policy': runsheet['parameter_policy'],
        'runs': payload_runs,
        'interaction_label_counts': dict(Counter(row['interaction_label'] for row in rows)),
        'coarse_label_counts': dict(Counter(row['coarse_label'] for row in rows)),
        'observation': observation,
        'conclusion': conclusion,
    }

    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug=str(runsheet.get('experiment_slug', 'stage_c0_1_combinatorial_kernel_probe')),
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
