#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import (
    ATLAS_NOTES,
    REPO_ROOT,
    append_log,
    anisotropy_ratio,
    build_transverse_setup,
    coherence_score,
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
from stage15_relational_geometry_probe import (
    ANALYSIS_STRIDE,
    first_arrival_time,
    matrix_mean_absolute_asymmetry,
    ordering_persistence,
    overlap_matrix,
    phase_lock_matrix,
    sample_component_envelopes,
    transport_order_breakdown,
    triangle_stats,
)

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage16_relational_feedback_runs.json'
SETUP_CACHE: dict[tuple[int, str, float, float], tuple[Any, Any]] = {}
STABILIZATION_MIN_LIFETIME_DELTA = 0.05
STABILIZATION_MIN_BINDING_DELTA = 0.03
STABILIZATION_MIN_COARSE_DELTA = 0.02

CSV_FIELDS = [
    'run_id',
    'phase',
    'case_id',
    'condition_id',
    'condition_label',
    'relational_driver',
    'driver_strength_class',
    'feedback_eps',
    'n_side',
    'interaction_label',
    'coarse_label',
    'composite_lifetime',
    'binding_persistence_score',
    'center_shift',
    'morphology_drift_reduction',
    'initial_mean_pair_distance',
    'min_mean_pair_distance',
    'final_mean_pair_distance',
    'max_mean_overlap',
    'phase_lock_indicator',
    'dominant_basin_area_fraction_sigma4',
    'coarse_persistence_sigma4',
    'mean_basin_count_sigma4',
    'metric_violation_rate',
    'mean_absolute_asymmetry',
    'ordering_persistence_under_refinement',
    'transport_ordering_breakdown',
    'constraint_max',
    'sector_leakage',
    'relational_stabilization_score',
    'compare_to_null_lifetime_delta',
    'compare_to_null_binding_delta',
    'compare_to_null_coarse_persistence_delta',
    'feedback_induced_regime_change',
    'output_label',
    'promoted_followup',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage 16 relational feedback dynamics pilot.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--case-ids', nargs='*', default=[])
    parser.add_argument('--condition-ids', nargs='*', default=[])
    parser.add_argument('--skip-followup', action='store_true')
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def selected_items(items: list[dict[str, Any]], ids: list[str], key: str) -> list[dict[str, Any]]:
    if not ids:
        return items
    wanted = set(ids)
    return [item for item in items if item[key] in wanted]


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


def interaction_rank(label: str) -> int:
    order = {
        'dispersive independent regime': 0,
        'transient binding regime': 1,
        'metastable composite regime': 2,
        'oscillatory exchange regime': 3,
        'large-scale drift field regime': 3,
    }
    return order.get(label, 0)


def pairwise_distance_stats(centers: list[np.ndarray], boundary_type: str) -> tuple[float, float, float]:
    values: list[float] = []
    for i in range(len(centers)):
        for j in range(i + 1, len(centers)):
            delta = displacement(np.asarray([centers[j]], dtype=float), np.asarray(centers[i], dtype=float), boundary_type)[0]
            values.append(float(np.linalg.norm(delta)))
    if not values:
        return 0.0, 0.0, 0.0
    arr = np.asarray(values, dtype=float)
    return float(np.mean(arr)), float(np.min(arr)), float(np.max(arr))


def pairwise_phase_metric(component_qs: list[np.ndarray]) -> float:
    mat = phase_lock_matrix(component_qs)
    if mat.shape[0] < 2:
        return 0.0
    values = [float(mat[i, j]) for i in range(mat.shape[0]) for j in range(i + 1, mat.shape[0])]
    return float(np.mean(values)) if values else 0.0


def row_normalize(weights: np.ndarray) -> np.ndarray:
    mat = np.asarray(weights, dtype=float).copy()
    np.fill_diagonal(mat, 0.0)
    mat = np.clip(mat, 0.0, None)
    row_sums = np.sum(mat, axis=1, keepdims=True)
    valid = row_sums[:, 0] > 1.0e-12
    mat[valid] /= row_sums[valid]
    return mat


def driver_weight_matrix(driver_id: str, component_qs: list[np.ndarray], env_samples: np.ndarray) -> np.ndarray:
    if driver_id == 'persistence_distance':
        weights = overlap_matrix(component_qs)
        return row_normalize(weights)
    if driver_id == 'phase_correlation_distance':
        weights = phase_lock_matrix(component_qs)
        return row_normalize(weights)
    if driver_id == 'transport_delay_distance':
        return row_normalize(env_samples)
    raise ValueError(f'unsupported relational driver: {driver_id}')


def directional_feedback_term(
    q: np.ndarray,
    center: np.ndarray,
    centers: list[np.ndarray],
    weight_row: np.ndarray,
    bandwidth: float,
    boundary_type: str,
    midpoints: np.ndarray,
    projector,
) -> np.ndarray:
    force_vec = np.zeros(3, dtype=float)
    for j, weight in enumerate(weight_row):
        if weight <= 0.0:
            continue
        delta = displacement(np.asarray([centers[j]], dtype=float), np.asarray(center, dtype=float), boundary_type)[0]
        force_vec += float(weight) * delta
    if float(np.linalg.norm(force_vec)) <= 1.0e-12:
        return np.zeros_like(q)
    local_delta = displacement(midpoints, center, boundary_type)
    directional = np.sum(local_delta * force_vec[None, :], axis=1) / max(float(bandwidth * bandwidth), 1.0e-12)
    directional = np.tanh(directional)
    return np.asarray(projector(directional * q), dtype=float)


def build_distance_payload(
    packet_count: int,
    overlap_series: dict[tuple[int, int], list[float]],
    phase_series: dict[tuple[int, int], list[float]],
    transport_series: dict[tuple[int, int], list[float]],
    times: list[float],
) -> dict[str, dict[str, Any]]:
    persistence_sym = np.eye(packet_count, dtype=float)
    phase_sym = np.eye(packet_count, dtype=float)
    transport_raw = np.zeros((packet_count, packet_count), dtype=float)

    for i in range(packet_count):
        for j in range(i + 1, packet_count):
            mean_overlap = float(np.mean(overlap_series[(i, j)])) if overlap_series[(i, j)] else 0.0
            persistence_sym[i, j] = 1.0 / (0.001 + mean_overlap)
            persistence_sym[j, i] = persistence_sym[i, j]
            mean_phase = float(np.mean(np.abs(phase_series[(i, j)]))) if phase_series[(i, j)] else 0.0
            phase_sym[i, j] = 1.0 - mean_phase
            phase_sym[j, i] = phase_sym[i, j]

    for i in range(packet_count):
        for j in range(packet_count):
            if i == j:
                continue
            transport_raw[i, j] = first_arrival_time(transport_series[(i, j)], times, 0.15)
    transport_sym = 0.5 * (transport_raw + transport_raw.T)

    payloads: dict[str, dict[str, Any]] = {}
    for distance_id, raw_mat, sym_mat in [
        ('persistence_distance', persistence_sym, persistence_sym),
        ('phase_correlation_distance', phase_sym, phase_sym),
        ('transport_delay_distance', transport_raw, transport_sym),
    ]:
        mean_asym, max_asym = matrix_mean_absolute_asymmetry(raw_mat)
        tri_rate, tri_mean, _tri_max = triangle_stats(sym_mat)
        payloads[distance_id] = {
            'raw_matrix': raw_mat.tolist(),
            'sym_matrix': sym_mat.tolist(),
            'mean_absolute_asymmetry': mean_asym,
            'max_absolute_asymmetry': max_asym,
            'triangle_violation_rate': tri_rate,
            'mean_violation_magnitude': tri_mean,
        }
    return payloads


def render_grid_slice(ax, grid: np.ndarray, title: str, cmap: str = 'viridis') -> None:
    slice_idx = grid.shape[2] // 2
    img = grid[:, :, slice_idx].T
    im = ax.imshow(img, origin='lower', cmap=cmap)
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)


def simulate_run(case: dict[str, Any], condition: dict[str, Any], resolution: int, defaults: dict[str, Any], base: dict[str, Any], work_plot_dir: Path) -> tuple[dict[str, Any], list[Path]]:
    data, projector = get_cached_setup(resolution, defaults, str(base['boundary_type']))
    kick_axis = int(defaults['kick_axis'])
    bandwidth = float(base['bandwidth'])
    feedback_eps = float(condition['feedback_eps'])
    driver_id = str(case['relational_driver'])

    packet_qs: list[np.ndarray] = []
    packet_vs: list[np.ndarray] = []
    for center, amp, phase, kick_sign in zip(case['packet_centers'], case['packet_amplitudes'], case['phase_offsets_rad'], case['kick_signs']):
        q0, v0 = build_component_packet(
            center=np.asarray(center, dtype=float),
            amplitude=float(amp),
            phase_offset=float(phase),
            kick_sign=float(kick_sign),
            bandwidth=bandwidth,
            central_k=float(base['central_k']),
            phase_pattern=str(base['phase_pattern']),
            kick_axis=kick_axis,
            boundary_type=str(base['boundary_type']),
            midpoints=data.midpoints,
            edge_axes=data.edge_axes,
            projector=projector,
        )
        packet_qs.append(q0)
        packet_vs.append(v0)

    operator = data.L1
    dt = suggested_dt(operator, float(defaults['dt_safety']))
    steps = max(2, int(np.ceil(float(defaults['t_final']) / dt)))
    leakage_fn = make_leakage_fn(projector)

    times: list[float] = []
    total_centers: list[list[float]] = []
    widths: list[float] = []
    norms: list[float] = []
    energies: list[float] = []
    anisotropies: list[float] = []
    sector_leakages: list[float] = []
    constraint_norms: list[float] = []
    coherences: list[float] = []
    pair_distance_series: list[float] = []
    overlap_series_scalar: list[float] = []
    phase_alignment_series: list[float] = []
    overlap_series: dict[tuple[int, int], list[float]] = {
        (i, j): [] for i in range(int(case['packet_count'])) for j in range(i + 1, int(case['packet_count']))
    }
    phase_series: dict[tuple[int, int], list[float]] = {
        (i, j): [] for i in range(int(case['packet_count'])) for j in range(i + 1, int(case['packet_count']))
    }
    transport_series: dict[tuple[int, int], list[float]] = {
        (i, j): [] for i in range(int(case['packet_count'])) for j in range(int(case['packet_count'])) if i != j
    }
    coarse: dict[float, dict[str, Any]] = {
        sigma: {
            'basin_counts': [],
            'dominant_area_fractions': [],
            'envelope_variances': [],
            'stable_flags': [],
            'dominant_masks': [],
            'last_env': None,
        }
        for sigma in SIGMAS
    }

    def record(t: float) -> None:
        total_q = np.sum(packet_qs, axis=0)
        total_v = np.sum(packet_vs, axis=0)
        total_weights = np.abs(total_q) ** 2
        total_center = weighted_center(data.midpoints, total_weights, str(base['boundary_type']))
        times.append(float(t))
        total_centers.append(total_center.tolist())
        widths.append(packet_width(data.midpoints, total_weights, total_center, str(base['boundary_type'])))
        norms.append(float(np.linalg.norm(total_q)))
        energies.append(state_energy(operator, total_q, total_v))
        anisotropies.append(anisotropy_ratio(data.midpoints, total_weights, total_center, str(base['boundary_type'])))
        sector_leakages.append(float(leakage_fn(total_q)))
        constraint_norms.append(float(np.linalg.norm(data.d0.T @ total_q)))
        coherences.append(coherence_score(total_weights))

        component_centers: list[np.ndarray] = []
        for comp_q in packet_qs:
            comp_weights = np.abs(comp_q) ** 2
            component_centers.append(weighted_center(data.midpoints, comp_weights, str(base['boundary_type'])))
        mean_pair_dist, _min_pair_dist, _max_pair_dist = pairwise_distance_stats(component_centers, str(base['boundary_type']))
        pair_distance_series.append(mean_pair_dist)

        ov = overlap_matrix(packet_qs)
        ph = phase_lock_matrix(packet_qs)
        env_samples = sample_component_envelopes(data.midpoints, packet_qs, resolution, case['packet_centers'])
        overlap_values: list[float] = []
        phase_values: list[float] = []
        for i in range(int(case['packet_count'])):
            for j in range(i + 1, int(case['packet_count'])):
                overlap_values.append(float(ov[i, j]))
                phase_values.append(float(ph[i, j]))
                overlap_series[(i, j)].append(float(ov[i, j]))
                phase_series[(i, j)].append(float(ph[i, j]))
        overlap_series_scalar.append(float(np.mean(overlap_values)) if overlap_values else 0.0)
        phase_alignment_series.append(float(np.mean(phase_values)) if phase_values else 0.0)
        for i in range(int(case['packet_count'])):
            for j in range(int(case['packet_count'])):
                if i == j:
                    continue
                transport_series[(i, j)].append(float(env_samples[i, j]))

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
        component_centers = []
        for comp_q in packet_qs:
            comp_weights = np.abs(comp_q) ** 2
            component_centers.append(weighted_center(data.midpoints, comp_weights, str(base['boundary_type'])))
        env_samples = sample_component_envelopes(data.midpoints, packet_qs, resolution, [center.tolist() for center in component_centers])
        weight_mat = driver_weight_matrix(driver_id, packet_qs, env_samples)
        next_qs: list[np.ndarray] = []
        next_vs: list[np.ndarray] = []
        for idx, (q, v, center) in enumerate(zip(packet_qs, packet_vs, component_centers)):
            rel_term = directional_feedback_term(
                q=q,
                center=center,
                centers=component_centers,
                weight_row=weight_mat[idx],
                bandwidth=bandwidth,
                boundary_type=str(base['boundary_type']),
                midpoints=data.midpoints,
                projector=projector,
            )
            acc = np.asarray(projector(-(operator @ q) + feedback_eps * rel_term), dtype=float)
            v_half = v + 0.5 * dt * acc
            q_new = np.asarray(projector(q + dt * v_half), dtype=float)
            rel_term_new = directional_feedback_term(
                q=q_new,
                center=center,
                centers=component_centers,
                weight_row=weight_mat[idx],
                bandwidth=bandwidth,
                boundary_type=str(base['boundary_type']),
                midpoints=data.midpoints,
                projector=projector,
            )
            acc_new = np.asarray(projector(-(operator @ q_new) + feedback_eps * rel_term_new), dtype=float)
            v_new = np.asarray(projector(v_half + 0.5 * dt * acc_new), dtype=float)
            next_qs.append(q_new)
            next_vs.append(v_new)
        packet_qs = next_qs
        packet_vs = next_vs
        if ((step_idx + 1) % ANALYSIS_STRIDE == 0) or (step_idx == steps - 1):
            record((step_idx + 1) * dt)

    sample_dt = float(np.median(np.diff(times))) if len(times) >= 2 else dt
    centers_arr = np.asarray(total_centers, dtype=float)
    center_shift = float(np.linalg.norm(displacement(centers_arr[-1][None, :], centers_arr[0], str(base['boundary_type']))[0]))
    mean_pair_dist, min_pair_dist, _max_pair_dist = pairwise_distance_stats(
        [np.asarray(center, dtype=float) for center in total_centers[-min(int(case['packet_count']), len(total_centers)):]],
        str(base['boundary_type']),
    )
    interaction_summary = {
        'packet_count': int(case['packet_count']),
        'center_shift': center_shift,
        'width_ratio': float(widths[-1] / max(widths[0], 1.0e-12)),
        'max_mean_overlap': float(np.max(overlap_series_scalar)) if overlap_series_scalar else 0.0,
        'phase_lock_indicator': float(np.mean(phase_alignment_series)) if phase_alignment_series else 0.0,
        'composite_lifetime': float(np.sum(np.asarray(overlap_series_scalar) >= 0.35) * sample_dt),
        'oscillation_count': 0,
        'initial_mean_pair_distance': float(pair_distance_series[0]) if pair_distance_series else 0.0,
        'min_mean_pair_distance': float(np.min(pair_distance_series)) if pair_distance_series else min_pair_dist,
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

    binding_mask = (np.asarray(pair_distance_series) <= 0.25) & (np.asarray(overlap_series_scalar) >= 0.35)
    binding_persistence = float(np.mean(binding_mask.astype(float))) if binding_mask.size else 0.0

    distance_payloads = build_distance_payload(int(case['packet_count']), overlap_series, phase_series, transport_series, times)
    driver_payload = distance_payloads[driver_id]
    row = {
        'run_id': f"{case['case_id']}__{condition['condition_id']}__n{resolution}",
        'phase': 'base',
        'case_id': case['case_id'],
        'condition_id': condition['condition_id'],
        'condition_label': condition['label'],
        'relational_driver': driver_id,
        'driver_strength_class': case['driver_strength_class'],
        'feedback_eps': feedback_eps,
        'n_side': int(resolution),
        'interaction_label': interaction_label,
        'coarse_label': coarse_label,
        'composite_lifetime': float(interaction_summary['composite_lifetime']),
        'binding_persistence_score': binding_persistence,
        'center_shift': center_shift,
        'morphology_drift_reduction': 0.0,
        'initial_mean_pair_distance': float(pair_distance_series[0]) if pair_distance_series else 0.0,
        'min_mean_pair_distance': float(np.min(pair_distance_series)) if pair_distance_series else 0.0,
        'final_mean_pair_distance': float(pair_distance_series[-1]) if pair_distance_series else 0.0,
        'max_mean_overlap': float(np.max(overlap_series_scalar)) if overlap_series_scalar else 0.0,
        'phase_lock_indicator': float(np.mean(phase_alignment_series)) if phase_alignment_series else 0.0,
        'dominant_basin_area_fraction_sigma4': coarse_summary['mean_dominant_area_fraction_sigma4'],
        'coarse_persistence_sigma4': coarse_summary['coarse_persistence_sigma4'],
        'mean_basin_count_sigma4': coarse_summary['mean_basin_count_sigma4'],
        'metric_violation_rate': float(driver_payload['triangle_violation_rate']) if driver_payload['triangle_violation_rate'] is not None else 0.0,
        'mean_absolute_asymmetry': float(driver_payload['mean_absolute_asymmetry']),
        'ordering_persistence_under_refinement': None,
        'transport_ordering_breakdown': None,
        'constraint_max': float(np.max(constraint_norms)) if constraint_norms else 0.0,
        'sector_leakage': float(np.max(sector_leakages)) if sector_leakages else 0.0,
        'relational_stabilization_score': 0.0,
        'compare_to_null_lifetime_delta': 0.0,
        'compare_to_null_binding_delta': 0.0,
        'compare_to_null_coarse_persistence_delta': 0.0,
        'feedback_induced_regime_change': 'none',
        'output_label': 'no relational dynamic effect',
        'promoted_followup': 0,
        'notes': case['selection_note'],
    }

    payload = {
        'run_id': row['run_id'],
        'case': case,
        'condition': condition,
        'n_side': int(resolution),
        'interaction_label': interaction_label,
        'coarse_label': coarse_label,
        'summary': row,
        'times': times,
        'pair_distance_series': pair_distance_series,
        'overlap_series': overlap_series_scalar,
        'phase_alignment_series': phase_alignment_series,
        'driver_payload': driver_payload,
        'driver_distance_payloads': distance_payloads,
        'coarse_sigma4': {
            'dominant_area_fractions': coarse[4.0]['dominant_area_fractions'],
            'basin_counts': coarse[4.0]['basin_counts'],
            'stable_flags': coarse[4.0]['stable_flags'],
            'last_env': coarse[4.0]['last_env'],
        },
    }

    plot_paths: list[Path] = []
    trace_path = work_plot_dir / f"stage16_{row['run_id']}_pair_distance_overlap.png"
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))
    axes[0].plot(times, pair_distance_series, color='tab:blue')
    axes[0].set_title('Mean pair distance')
    axes[0].set_xlabel('time')
    axes[0].grid(alpha=0.25)
    axes[1].plot(times, overlap_series_scalar, color='tab:orange')
    axes[1].set_title('Mean overlap')
    axes[1].set_xlabel('time')
    axes[1].grid(alpha=0.25)
    fig.suptitle(f"Stage 16: {row['run_id']}")
    fig.savefig(trace_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(trace_path)

    basin_path = work_plot_dir / f"stage16_{row['run_id']}_basin_trace.png"
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))
    axes[0].plot(times, coarse[4.0]['dominant_area_fractions'], color='tab:green')
    axes[0].set_title('Dominant basin area (sigma=4)')
    axes[0].set_xlabel('time')
    axes[0].grid(alpha=0.25)
    axes[1].plot(times, coarse[4.0]['basin_counts'], color='tab:red')
    axes[1].set_title('Basin count (sigma=4)')
    axes[1].set_xlabel('time')
    axes[1].grid(alpha=0.25)
    fig.suptitle(f"Basin trace: {row['run_id']}")
    fig.savefig(basin_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(basin_path)

    env = coarse[4.0]['last_env']
    if env is not None:
        env_path = work_plot_dir / f"stage16_{row['run_id']}_envelope_snapshot.png"
        fig, ax = plt.subplots(figsize=(5.4, 4.6))
        render_grid_slice(ax, env, f"Envelope sigma=4: {row['run_id']}")
        fig.savefig(env_path, dpi=180, bbox_inches='tight')
        plt.close(fig)
        plot_paths.append(env_path)

    return payload, plot_paths


def compare_to_null(base_payloads: dict[tuple[str, str], dict[str, Any]], rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, Any] | None]:
    promoted: dict[str, Any] | None = None
    best_score = None
    for row in rows:
        if row['condition_id'] == 'null_control' or row['phase'] != 'base':
            continue
        null_row = base_payloads[(row['case_id'], 'null_control')]['summary']
        lifetime_delta = float(row['composite_lifetime']) - float(null_row['composite_lifetime'])
        binding_delta = float(row['binding_persistence_score']) - float(null_row['binding_persistence_score'])
        coarse_delta = float(row['coarse_persistence_sigma4']) - float(null_row['coarse_persistence_sigma4'])
        drift_reduction = float(null_row['center_shift']) - float(row['center_shift'])
        row['compare_to_null_lifetime_delta'] = lifetime_delta
        row['compare_to_null_binding_delta'] = binding_delta
        row['compare_to_null_coarse_persistence_delta'] = coarse_delta
        row['morphology_drift_reduction'] = drift_reduction
        row['relational_stabilization_score'] = (
            0.5 * lifetime_delta
            + 0.3 * binding_delta
            + 0.15 * coarse_delta
            + 0.05 * max(drift_reduction, 0.0)
        )
        changes = []
        if str(row['interaction_label']) != str(null_row['interaction_label']):
            changes.append(f"interaction:{null_row['interaction_label']}->{row['interaction_label']}")
        if str(row['coarse_label']) != str(null_row['coarse_label']):
            changes.append(f"coarse:{null_row['coarse_label']}->{row['coarse_label']}")
        row['feedback_induced_regime_change'] = '; '.join(changes) if changes else 'none'

    strong_nonstable = any(
        row['condition_id'] != 'null_control'
        and row['driver_strength_class'] != 'scale-stable metric-like ordering'
        and (
            float(row['compare_to_null_lifetime_delta']) >= STABILIZATION_MIN_LIFETIME_DELTA
            or float(row['compare_to_null_binding_delta']) >= STABILIZATION_MIN_BINDING_DELTA
        )
        for row in rows
    )

    for row in rows:
        if row['condition_id'] == 'null_control' or row['phase'] != 'base':
            continue
        lifetime_delta = float(row['compare_to_null_lifetime_delta'])
        binding_delta = float(row['compare_to_null_binding_delta'])
        coarse_delta = float(row['compare_to_null_coarse_persistence_delta'])
        driver_stable = str(row['driver_strength_class']) == 'scale-stable metric-like ordering'
        if abs(lifetime_delta) < 0.01 and abs(binding_delta) < 0.01 and abs(coarse_delta) < 0.01 and row['feedback_induced_regime_change'] == 'none':
            row['output_label'] = 'no relational dynamic effect'
        elif (not driver_stable) and strong_nonstable and interaction_rank(str(row['interaction_label'])) >= 2:
            row['output_label'] = 'global artificial binding'
        elif (
            driver_stable
            and lifetime_delta >= STABILIZATION_MIN_LIFETIME_DELTA
            and binding_delta >= STABILIZATION_MIN_BINDING_DELTA
            and coarse_delta >= STABILIZATION_MIN_COARSE_DELTA
            and interaction_rank(str(row['interaction_label'])) >= interaction_rank(str(base_payloads[(row['case_id'], 'null_control')]['summary']['interaction_label']))
        ):
            row['output_label'] = 'selective regime stabilization'
        else:
            row['output_label'] = 'weak stabilization'

        if row['output_label'] == 'selective regime stabilization':
            score = float(row['relational_stabilization_score'])
            if best_score is None or score > best_score:
                best_score = score
                promoted = {
                    'case_id': row['case_id'],
                    'condition_id': row['condition_id'],
                }
    if promoted is not None:
        for row in rows:
            if row['case_id'] == promoted['case_id'] and row['condition_id'] == promoted['condition_id']:
                row['promoted_followup'] = 1
    return rows, promoted


def add_refinement_metrics(
    base_payload: dict[str, Any],
    refined_payload: dict[str, Any],
    row: dict[str, Any],
) -> None:
    driver_id = str(row['relational_driver'])
    base_driver = base_payload['driver_distance_payloads'][driver_id]
    refined_driver = refined_payload['driver_distance_payloads'][driver_id]
    ordering_flip, _rank_corr = ordering_persistence(
        np.asarray(base_driver['sym_matrix'], dtype=float),
        np.asarray(refined_driver['sym_matrix'], dtype=float),
    )
    row['ordering_persistence_under_refinement'] = 1.0 - float(ordering_flip)
    if driver_id == 'transport_delay_distance':
        row['transport_ordering_breakdown'] = transport_order_breakdown(
            np.asarray(base_driver['raw_matrix'], dtype=float),
            np.asarray(refined_driver['raw_matrix'], dtype=float),
        )
    if (
        row['output_label'] in {'selective regime stabilization', 'weak stabilization'}
        and row['ordering_persistence_under_refinement'] is not None
        and float(row['ordering_persistence_under_refinement']) < 0.75
    ):
        row['output_label'] = 'weak stabilization'


def create_response_matrix(path: Path, rows: list[dict[str, Any]]) -> None:
    base_rows = [row for row in rows if row['phase'] == 'base']
    cases = sorted({row['case_id'] for row in base_rows})
    conditions = sorted({row['condition_id'] for row in base_rows}, key=lambda item: ['null_control', 'rel_feedback_eps001', 'rel_feedback_eps002'].index(item) if item in ['null_control', 'rel_feedback_eps001', 'rel_feedback_eps002'] else 99)
    label_map = {
        'no relational dynamic effect': 0,
        'weak stabilization': 1,
        'selective regime stabilization': 2,
        'global artificial binding': 3,
    }
    data = np.zeros((len(cases), len(conditions)), dtype=float)
    for i, case_id in enumerate(cases):
        for j, cond in enumerate(conditions):
            row = next(item for item in base_rows if item['case_id'] == case_id and item['condition_id'] == cond)
            data[i, j] = label_map.get(str(row['output_label']), 0)
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    im = ax.imshow(data, cmap='viridis', vmin=0, vmax=3)
    ax.set_xticks(range(len(conditions)))
    xticklabels = {
        'null_control': 'null',
        'rel_feedback_eps001': '0.01',
        'rel_feedback_eps002': '0.02',
    }
    ax.set_xticklabels([xticklabels.get(cond, cond) for cond in conditions])
    ax.set_yticks(range(len(cases)))
    ax.set_yticklabels(cases)
    ax.set_title('Stage 16 relational-feedback response matrix')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_comparison_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    base_rows = [row for row in rows if row['phase'] == 'base']
    cases = sorted({row['case_id'] for row in base_rows})
    conditions = sorted({row['condition_id'] for row in base_rows}, key=lambda item: ['null_control', 'rel_feedback_eps001', 'rel_feedback_eps002'].index(item) if item in ['null_control', 'rel_feedback_eps001', 'rel_feedback_eps002'] else 99)
    labels = {'null_control': 'null', 'rel_feedback_eps001': '0.01', 'rel_feedback_eps002': '0.02'}
    x = np.arange(len(cases))
    width = 0.8 / max(len(conditions), 1)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    for idx, cond in enumerate(conditions):
        lifetimes = []
        coarse_pers = []
        for case_id in cases:
            row = next(item for item in base_rows if item['case_id'] == case_id and item['condition_id'] == cond)
            lifetimes.append(float(row['composite_lifetime']))
            coarse_pers.append(float(row['coarse_persistence_sigma4']))
        offset = idx - 0.5 * (len(conditions) - 1)
        axes[0].bar(x + offset * width, lifetimes, width=width, label=labels.get(cond, cond))
        axes[1].bar(x + offset * width, coarse_pers, width=width, label=labels.get(cond, cond))
    axes[0].set_title('Composite lifetime')
    axes[1].set_title('Coarse persistence (sigma=4)')
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(cases, rotation=20, ha='right')
        ax.grid(alpha=0.25, axis='y')
        ax.legend(fontsize=8)
    fig.suptitle('Stage 16 null vs relational feedback comparison')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_score_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    feedback_rows = [row for row in rows if row['phase'] == 'base' and row['condition_id'] != 'null_control']
    if not feedback_rows:
        fig, ax = plt.subplots(figsize=(6.4, 4.0))
        ax.text(0.5, 0.5, 'No feedback rows selected', ha='center', va='center')
        ax.set_axis_off()
        fig.savefig(path, dpi=180, bbox_inches='tight')
        plt.close(fig)
        return
    labels = [f"{row['case_id']}\n{row['condition_id']}" for row in feedback_rows]
    scores = [float(row['relational_stabilization_score']) for row in feedback_rows]
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    ax.bar(range(len(labels)), scores, color='tab:purple')
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=20, ha='right')
    ax.set_ylabel('stabilization score')
    ax.set_title('Stage 16 relational stabilization scores')
    ax.grid(alpha=0.25, axis='y')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_refinement_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.4))
    labels = [f"{row['condition_id']}\nn={row['n_side']}" for row in rows]
    lifetimes = [float(row['composite_lifetime']) for row in rows]
    orderings = [float(row['ordering_persistence_under_refinement'] or 0.0) for row in rows]
    axes[0].bar(range(len(labels)), lifetimes, color='tab:blue')
    axes[0].set_xticks(range(len(labels)))
    axes[0].set_xticklabels(labels)
    axes[0].set_title('Composite lifetime across refinement')
    axes[0].grid(alpha=0.25, axis='y')
    axes[1].bar(range(len(labels)), orderings, color='tab:orange')
    axes[1].set_xticks(range(len(labels)))
    axes[1].set_xticklabels(labels)
    axes[1].set_title('Ordering persistence under refinement')
    axes[1].grid(alpha=0.25, axis='y')
    fig.suptitle('Stage 16 refinement follow-up')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(
    path: Path,
    json_rel: str,
    csv_rel: str,
    rows: list[dict[str, Any]],
    label_counts: dict[str, int],
    promoted: dict[str, Any] | None,
) -> None:
    lines = [
        '# Stage 16 Relational Feedback Dynamics v1',
        '',
        f'Timestamped summary: `{json_rel}`',
        f'Timestamped run table: `{csv_rel}`',
        '',
        f'Output-label counts: {label_counts}',
        '',
        'This pilot tests whether Stage 15 relational ordering can act back on dynamics as a weak organizing influence without changing the frozen base operator architecture.',
        '',
        'Base comparisons:',
        '',
    ]
    for row in [item for item in rows if item['phase'] == 'base']:
        lines.extend([
            f"- `{row['case_id']}` / `{row['condition_id']}`",
            f"  - driver: `{row['relational_driver']}` ({row['driver_strength_class']})",
            f"  - interaction: `{row['interaction_label']}`",
            f"  - coarse: `{row['coarse_label']}`",
            f"  - lifetime delta vs null: `{row['compare_to_null_lifetime_delta']:.4f}`",
            f"  - binding delta vs null: `{row['compare_to_null_binding_delta']:.4f}`",
            f"  - coarse persistence delta vs null: `{row['compare_to_null_coarse_persistence_delta']:.4f}`",
            f"  - output: `{row['output_label']}`",
        ])
    lines.extend(['', f'Promoted follow-up: `{promoted}`' if promoted else 'Promoted follow-up: `none`', ''])
    refined_rows = [row for row in rows if row['phase'] == 'refinement']
    if refined_rows:
        lines.append('Refinement follow-up:')
        lines.append('')
        for row in refined_rows:
            lines.extend([
                f"- `{row['case_id']}` / `{row['condition_id']}` / `n={row['n_side']}`",
                f"  - ordering persistence under refinement: `{row['ordering_persistence_under_refinement']}`",
                f"  - output: `{row['output_label']}`",
            ])
        lines.append('')
    lines.extend([
        'Interpretation boundary:',
        '- no spacetime or physical-force claim is made here',
        '- any positive effect is read only as weak relational organization acting back on morphology',
        '',
    ])
    path.write_text('\n'.join(lines), encoding='utf-8')


def main() -> None:
    args = parse_args()
    defaults = load_stage10_defaults()
    runsheet = load_runsheet(args.runsheet)
    base = runsheet['base_seed_reference']
    cases = selected_items(runsheet['cases'], args.case_ids, 'case_id')
    conditions = selected_items(runsheet['conditions'], args.condition_ids, 'condition_id')
    work_plot_dir = Path('/tmp') / 'haos_iip_stage16_feedback_plots'
    work_plot_dir.mkdir(parents=True, exist_ok=True)

    base_payloads: dict[tuple[str, str], dict[str, Any]] = {}
    base_rows: list[dict[str, Any]] = []
    plot_paths: list[Path] = []

    for case in cases:
        for condition in conditions:
            payload, run_plots = simulate_run(
                case=case,
                condition=condition,
                resolution=int(runsheet['resolution_policy']['base_resolution']),
                defaults=defaults,
                base=base,
                work_plot_dir=work_plot_dir,
            )
            base_payloads[(case['case_id'], condition['condition_id'])] = payload
            base_rows.append(payload['summary'])
            plot_paths.extend(run_plots)

    base_rows, promoted = compare_to_null(base_payloads, base_rows)
    refined_rows: list[dict[str, Any]] = []
    if promoted is not None and not args.skip_followup:
        case = next(item for item in cases if item['case_id'] == promoted['case_id'])
        follow_conditions = [
            next(item for item in conditions if item['condition_id'] == 'null_control'),
            next(item for item in conditions if item['condition_id'] == promoted['condition_id']),
        ]
        refined_payloads: dict[str, dict[str, Any]] = {}
        for condition in follow_conditions:
            payload, run_plots = simulate_run(
                case=case,
                condition=condition,
                resolution=int(runsheet['resolution_policy']['refined_resolution']),
                defaults=defaults,
                base=base,
                work_plot_dir=work_plot_dir,
            )
            payload['summary']['phase'] = 'refinement'
            null_base = base_payloads[(case['case_id'], 'null_control')]['summary']
            payload['summary']['compare_to_null_lifetime_delta'] = float(payload['summary']['composite_lifetime']) - float(null_base['composite_lifetime'])
            payload['summary']['compare_to_null_binding_delta'] = float(payload['summary']['binding_persistence_score']) - float(null_base['binding_persistence_score'])
            payload['summary']['compare_to_null_coarse_persistence_delta'] = float(payload['summary']['coarse_persistence_sigma4']) - float(null_base['coarse_persistence_sigma4'])
            payload['summary']['morphology_drift_reduction'] = float(null_base['center_shift']) - float(payload['summary']['center_shift'])
            payload['summary']['relational_stabilization_score'] = (
                0.5 * float(payload['summary']['compare_to_null_lifetime_delta'])
                + 0.3 * float(payload['summary']['compare_to_null_binding_delta'])
                + 0.15 * float(payload['summary']['compare_to_null_coarse_persistence_delta'])
                + 0.05 * max(float(payload['summary']['morphology_drift_reduction']), 0.0)
            )
            refined_payloads[condition['condition_id']] = payload
            refined_rows.append(payload['summary'])
            plot_paths.extend(run_plots)

        if promoted['condition_id'] in refined_payloads:
            add_refinement_metrics(
                base_payload=base_payloads[(case['case_id'], promoted['condition_id'])],
                refined_payload=refined_payloads[promoted['condition_id']],
                row=refined_payloads[promoted['condition_id']]['summary'],
            )
        if 'null_control' in refined_payloads:
            refined_payloads['null_control']['summary']['ordering_persistence_under_refinement'] = 1.0

    all_rows = base_rows + refined_rows
    label_counts = dict(Counter(row['output_label'] for row in all_rows))

    comparison_plot = work_plot_dir / 'stage16_relational_feedback_comparison.png'
    create_comparison_plot(comparison_plot, base_rows)
    plot_paths.append(comparison_plot)

    response_plot = work_plot_dir / 'stage16_relational_feedback_response_matrix.png'
    create_response_matrix(response_plot, base_rows)
    plot_paths.append(response_plot)

    score_plot = work_plot_dir / 'stage16_relational_feedback_scores.png'
    create_score_plot(score_plot, base_rows)
    plot_paths.append(score_plot)

    if refined_rows:
        refinement_plot = work_plot_dir / 'stage16_relational_feedback_refinement.png'
        create_refinement_plot(refinement_plot, refined_rows)
        plot_paths.append(refinement_plot)

    result = {
        'stage': 'Stage 16',
        'description': runsheet['description'],
        'base_seed_reference': base,
        'resolution_policy': runsheet['resolution_policy'],
        'cases': cases,
        'base_runs': [
            {
                'case_id': row['case_id'],
                'condition_id': row['condition_id'],
                'summary': row,
            }
            for row in base_rows
        ],
        'refinement_runs': [
            {
                'case_id': row['case_id'],
                'condition_id': row['condition_id'],
                'summary': row,
            }
            for row in refined_rows
        ],
        'output_label_counts': label_counts,
        'promoted_followup': promoted,
        'conclusion': 'Stage 16 tests whether Stage 15 relational ordering can act back on dynamics as a weak organizing influence while preserving the frozen operator architecture.',
    }

    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug='stage16_relational_feedback',
        result=result,
        csv_rows=all_rows,
        csv_fieldnames=CSV_FIELDS,
        plot_paths=plot_paths,
    )
    note_path = ATLAS_NOTES / runsheet['note_name']
    write_note(note_path, str(json_path.relative_to(REPO_ROOT)), str(csv_path.relative_to(REPO_ROOT)), all_rows, label_counts, promoted)
    append_log(
        title=f"Stage 16 Relational Feedback Dynamics ({json_path.stem})",
        config_summary=f"cases={[case['case_id'] for case in cases]}, conditions={[cond['condition_id'] for cond in conditions]}",
        result_path=json_path,
        stamped_plots=stamped_plots,
        observation=f"output_label_counts={label_counts}",
        conclusion='the Stage 16 pilot tests whether metric-like relational ordering can act as a weak organizing influence on collective morphology',
    )

    print(f'wrote {json_path}')
    print(f'wrote {csv_path}')
    for item in stamped_plots:
        print(item)
    print(f'output_label_counts={label_counts}')
    print(f'promoted_followup={promoted}')


if __name__ == '__main__':
    main()
