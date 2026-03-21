#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
from scipy.signal import hilbert

from stage10_common import (
    ATLAS_NOTES,
    REPO_ROOT,
    append_log,
    anisotropy_ratio,
    build_transverse_setup,
    coherence_score,
    displacement,
    gaussian_packet,
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

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage17_phase_sync_runs.json'
SETUP_CACHE: dict[tuple[int, str, float, float], tuple[Any, Any]] = {}
R_THRESHOLD = 0.7
CLUSTER_PHASE_THRESHOLD = 0.6
SUSTAINED_GAIN_DELTA = 0.03
FRACTION_GAIN_DELTA = 0.08
VARIANCE_DROP_DELTA = -0.03
FREQ_SPREAD_DROP_DELTA = -0.01
CLUSTER_DWELL_GAIN_DELTA = 0.05
MORPH_DELTA_MIN = 0.01

CSV_FIELDS = [
    'run_id',
    'phase',
    'representative_id',
    'condition_id',
    'phase_proxy',
    'control_phase_proxy',
    'phase_coupling_eps',
    'resolution',
    'boundary_type',
    'operator_sector',
    'packet_count',
    'interaction_label',
    'coarse_label',
    'mean_R_canonical',
    'median_R_canonical',
    'R_fraction_canonical',
    'mean_R_control',
    'median_R_control',
    'R_fraction_control',
    'phase_variance_mean_canonical',
    'phase_variance_mean_control',
    'frequency_spread_mean_canonical',
    'frequency_spread_mean_control',
    'cluster_count_mean_canonical',
    'cluster_count_mean_control',
    'cluster_dwell_fraction_canonical',
    'cluster_dwell_fraction_control',
    'mean_R_delta',
    'median_R_delta',
    'R_fraction_delta',
    'phase_variance_delta',
    'frequency_spread_delta',
    'cluster_dwell_delta',
    'composite_lifetime',
    'binding_persistence',
    'coarse_persistence',
    'composite_lifetime_delta',
    'binding_persistence_delta',
    'coarse_persistence_delta',
    'morphology_correlate_flag',
    'sync_label',
    'promoted_followup',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage 17 phase synchronization pilot.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    parser.add_argument('--skip-followup', action='store_true')
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def selected_runs(runs: list[dict[str, Any]], run_ids: list[str]) -> list[dict[str, Any]]:
    if not run_ids:
        return runs
    wanted = set(run_ids)
    return [run for run in runs if run['run_id'] in wanted]


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


def representative_lookup(runsheet: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {item['representative_id']: item for item in runsheet['representatives']}


def condition_lookup(runsheet: dict[str, Any]) -> dict[float, dict[str, Any]]:
    return {float(item['phase_coupling_eps']): item for item in runsheet['conditions']}


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


def local_packet_signal(midpoints: np.ndarray, q: np.ndarray, center: np.ndarray, bandwidth: float, boundary_type: str) -> float:
    envelope, _delta = gaussian_packet(
        points=midpoints,
        center=center,
        sigma=bandwidth,
        boundary_type=boundary_type,
        central_k=0.0,
        phase_pattern='flat',
        kick_axis=0,
    )
    weight_sum = max(float(np.sum(envelope)), 1.0e-12)
    return float(np.sum(envelope * q) / weight_sum)


def dominant_mode_phase(
    center: np.ndarray,
    q: np.ndarray,
    bandwidth: float,
    central_k: float,
    phase_pattern: str,
    kick_axis: int,
    boundary_type: str,
    midpoints: np.ndarray,
    edge_axes: list[str],
    projector,
) -> tuple[float, complex]:
    cos_basis, _ = build_component_packet(
        center=center,
        amplitude=1.0,
        phase_offset=0.0,
        kick_sign=1.0,
        bandwidth=bandwidth,
        central_k=central_k,
        phase_pattern=phase_pattern,
        kick_axis=kick_axis,
        boundary_type=boundary_type,
        midpoints=midpoints,
        edge_axes=edge_axes,
        projector=projector,
    )
    sin_basis, _ = build_component_packet(
        center=center,
        amplitude=1.0,
        phase_offset=0.5 * math.pi,
        kick_sign=1.0,
        bandwidth=bandwidth,
        central_k=central_k,
        phase_pattern=phase_pattern,
        kick_axis=kick_axis,
        boundary_type=boundary_type,
        midpoints=midpoints,
        edge_axes=edge_axes,
        projector=projector,
    )
    coeff = complex(float(np.vdot(cos_basis, q).real), float(np.vdot(sin_basis, q).real))
    return float(np.angle(coeff)), coeff


def local_weight_matrix(centers: list[np.ndarray], boundary_type: str, bandwidth: float) -> np.ndarray:
    n = len(centers)
    mat = np.zeros((n, n), dtype=float)
    scale = max(float(bandwidth), 1.0e-12)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            delta = displacement(np.asarray([centers[j]], dtype=float), np.asarray(centers[i], dtype=float), boundary_type)[0]
            dist = float(np.linalg.norm(delta))
            mat[i, j] = math.exp(-(dist * dist) / (2.0 * scale * scale))
    row_sums = np.sum(mat, axis=1, keepdims=True)
    valid = row_sums[:, 0] > 1.0e-12
    mat[valid] /= row_sums[valid]
    return mat


def quadrature_feedback_basis(
    center: np.ndarray,
    bandwidth: float,
    central_k: float,
    phase_pattern: str,
    kick_axis: int,
    boundary_type: str,
    midpoints: np.ndarray,
    edge_axes: list[str],
    projector,
) -> np.ndarray:
    basis, _ = build_component_packet(
        center=center,
        amplitude=1.0,
        phase_offset=0.5 * math.pi,
        kick_sign=1.0,
        bandwidth=bandwidth,
        central_k=central_k,
        phase_pattern=phase_pattern,
        kick_axis=kick_axis,
        boundary_type=boundary_type,
        midpoints=midpoints,
        edge_axes=edge_axes,
        projector=projector,
    )
    return basis


def circular_distance(a: float, b: float) -> float:
    return abs(math.atan2(math.sin(a - b), math.cos(a - b)))


def phase_cluster_signature(phases: np.ndarray, threshold: float) -> tuple[tuple[int, ...], ...]:
    n = len(phases)
    if n == 0:
        return tuple()
    visited = [False] * n
    comps: list[tuple[int, ...]] = []
    for i in range(n):
        if visited[i]:
            continue
        stack = [i]
        visited[i] = True
        comp = []
        while stack:
            idx = stack.pop()
            comp.append(idx)
            for j in range(n):
                if visited[j] or j == idx:
                    continue
                if circular_distance(float(phases[idx]), float(phases[j])) <= threshold:
                    visited[j] = True
                    stack.append(j)
        comps.append(tuple(sorted(comp)))
    comps.sort()
    return tuple(comps)


def phase_metrics(unwrapped_phases: np.ndarray, dt: float) -> dict[str, Any]:
    wrapped = np.angle(np.exp(1j * unwrapped_phases))
    R_t = np.abs(np.mean(np.exp(1j * wrapped), axis=1))
    phase_var = 1.0 - R_t
    if unwrapped_phases.shape[0] >= 2:
        omega = np.diff(unwrapped_phases, axis=0) / max(dt, 1.0e-12)
        freq_spread = np.std(omega, axis=1)
    else:
        freq_spread = np.zeros(1, dtype=float)

    signatures = [phase_cluster_signature(row, CLUSTER_PHASE_THRESHOLD) for row in wrapped]
    counts = np.asarray([len(sig) for sig in signatures], dtype=float)
    dominant = Counter(signatures).most_common(1)[0][0] if signatures else tuple()
    dwell = float(sum(1 for sig in signatures if sig == dominant) / max(len(signatures), 1))

    return {
        'R_t': R_t.tolist(),
        'mean_R': float(np.mean(R_t)),
        'median_R': float(np.median(R_t)),
        'R_fraction': float(np.mean(R_t > R_THRESHOLD)),
        'phase_variance_mean': float(np.mean(phase_var)),
        'phase_variance_median': float(np.median(phase_var)),
        'frequency_spread_mean': float(np.mean(freq_spread)),
        'frequency_spread_median': float(np.median(freq_spread)),
        'cluster_count_mean': float(np.mean(counts)) if counts.size else 0.0,
        'cluster_count_median': float(np.median(counts)) if counts.size else 0.0,
        'cluster_dwell_fraction': dwell,
    }


def interaction_rank(label: str) -> int:
    order = {
        'dispersive independent regime': 0,
        'transient binding regime': 1,
        'metastable composite regime': 2,
        'oscillatory exchange regime': 3,
        'large-scale drift field regime': 3,
    }
    return order.get(label, 0)


def simulate_run(run: dict[str, Any], rep: dict[str, Any], defaults: dict[str, Any], base: dict[str, Any], work_plot_dir: Path) -> tuple[dict[str, Any], list[Path]]:
    resolution = int(run['resolution'])
    data, projector = get_cached_setup(resolution, defaults, str(base['boundary_type']))
    kick_axis = int(defaults['kick_axis'])
    bandwidth = float(base['bandwidth'])
    central_k = float(base['central_k'])
    phase_pattern = str(base['phase_pattern'])
    boundary_type = str(base['boundary_type'])
    feedback_eps = float(run['phase_coupling_eps'])

    packet_qs: list[np.ndarray] = []
    packet_vs: list[np.ndarray] = []
    for center, amp, phase, kick_sign in zip(rep['packet_centers'], rep['packet_amplitudes'], rep['phase_offsets_rad'], rep['kick_signs']):
        q0, v0 = build_component_packet(
            center=np.asarray(center, dtype=float),
            amplitude=float(amp),
            phase_offset=float(phase),
            kick_sign=float(kick_sign),
            bandwidth=bandwidth,
            central_k=central_k,
            phase_pattern=phase_pattern,
            kick_axis=kick_axis,
            boundary_type=boundary_type,
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
    energies: list[float] = []
    anisotropies: list[float] = []
    sector_leakages: list[float] = []
    constraint_norms: list[float] = []
    coherences: list[float] = []
    pair_distance_series: list[float] = []
    overlap_series: list[float] = []
    local_signals: list[list[float]] = [[] for _ in range(int(rep['packet_count']))]
    control_coeffs: list[list[complex]] = [[] for _ in range(int(rep['packet_count']))]
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

        comp_centers: list[np.ndarray] = []
        for idx, comp_q in enumerate(packet_qs):
            comp_weights = np.abs(comp_q) ** 2
            center = weighted_center(data.midpoints, comp_weights, boundary_type)
            comp_centers.append(center)
            local_signals[idx].append(local_packet_signal(data.midpoints, comp_q, center, bandwidth, boundary_type))
            _phase, coeff = dominant_mode_phase(
                center=center,
                q=comp_q,
                bandwidth=bandwidth,
                central_k=central_k,
                phase_pattern=phase_pattern,
                kick_axis=kick_axis,
                boundary_type=boundary_type,
                midpoints=data.midpoints,
                edge_axes=data.edge_axes,
                projector=projector,
            )
            control_coeffs[idx].append(coeff)

        mean_pair_distance, min_pair_distance = pairwise_distance_stats(comp_centers, boundary_type)
        pair_distance_series.append(mean_pair_distance)
        ov = []
        for i in range(len(packet_qs)):
            for j in range(i + 1, len(packet_qs)):
                qi = packet_qs[i]
                qj = packet_qs[j]
                wi = np.abs(qi) ** 2
                wj = np.abs(qj) ** 2
                denom = math.sqrt(max(float(np.sum(wi)) * float(np.sum(wj)), 1.0e-12))
                ov.append(float(np.sum(np.sqrt(wi * wj)) / denom))
        overlap_series.append(float(np.mean(ov)) if ov else 0.0)

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
        centers = []
        phases = []
        quad_bases = []
        for comp_q in packet_qs:
            comp_weights = np.abs(comp_q) ** 2
            center = weighted_center(data.midpoints, comp_weights, boundary_type)
            centers.append(center)
            phase, _coeff = dominant_mode_phase(
                center=center,
                q=comp_q,
                bandwidth=bandwidth,
                central_k=central_k,
                phase_pattern=phase_pattern,
                kick_axis=kick_axis,
                boundary_type=boundary_type,
                midpoints=data.midpoints,
                edge_axes=data.edge_axes,
                projector=projector,
            )
            phases.append(phase)
            quad_bases.append(
                quadrature_feedback_basis(
                    center=center,
                    bandwidth=bandwidth,
                    central_k=central_k,
                    phase_pattern=phase_pattern,
                    kick_axis=kick_axis,
                    boundary_type=boundary_type,
                    midpoints=data.midpoints,
                    edge_axes=data.edge_axes,
                    projector=projector,
                )
            )
        weight_mat = local_weight_matrix(centers, boundary_type, bandwidth)

        next_qs: list[np.ndarray] = []
        next_vs: list[np.ndarray] = []
        for idx, (q, v) in enumerate(zip(packet_qs, packet_vs)):
            coupling_scalar = 0.0
            for j in range(len(packet_qs)):
                if j == idx:
                    continue
                coupling_scalar += float(weight_mat[idx, j]) * math.sin(phases[j] - phases[idx])
            feedback_term = coupling_scalar * max(float(np.linalg.norm(q)), 1.0e-12) * quad_bases[idx]
            acc = np.asarray(projector(-(operator @ q) + feedback_eps * feedback_term), dtype=float)
            v_half = v + 0.5 * dt * acc
            q_new = np.asarray(projector(q + dt * v_half), dtype=float)

            new_center = centers[idx]
            new_phase = phases[idx]
            new_quad = quad_bases[idx]
            if feedback_eps != 0.0:
                new_center = weighted_center(data.midpoints, np.abs(q_new) ** 2, boundary_type)
                new_phase, _ = dominant_mode_phase(
                    center=new_center,
                    q=q_new,
                    bandwidth=bandwidth,
                    central_k=central_k,
                    phase_pattern=phase_pattern,
                    kick_axis=kick_axis,
                    boundary_type=boundary_type,
                    midpoints=data.midpoints,
                    edge_axes=data.edge_axes,
                    projector=projector,
                )
                new_quad = quadrature_feedback_basis(
                    center=new_center,
                    bandwidth=bandwidth,
                    central_k=central_k,
                    phase_pattern=phase_pattern,
                    kick_axis=kick_axis,
                    boundary_type=boundary_type,
                    midpoints=data.midpoints,
                    edge_axes=data.edge_axes,
                    projector=projector,
                )
            new_coupling_scalar = 0.0
            for j in range(len(packet_qs)):
                if j == idx:
                    continue
                new_coupling_scalar += float(weight_mat[idx, j]) * math.sin(phases[j] - new_phase)
            acc_new = np.asarray(projector(-(operator @ q_new) + feedback_eps * (new_coupling_scalar * max(float(np.linalg.norm(q_new)), 1.0e-12) * new_quad)), dtype=float)
            v_new = np.asarray(projector(v_half + 0.5 * dt * acc_new), dtype=float)
            next_qs.append(q_new)
            next_vs.append(v_new)
        packet_qs = next_qs
        packet_vs = next_vs
        if ((step_idx + 1) % 4 == 0) or (step_idx == steps - 1):
            record((step_idx + 1) * dt)

    sample_dt = float(np.median(np.diff(times))) if len(times) >= 2 else dt
    canonical_phase = np.unwrap(np.angle(np.column_stack([hilbert(np.asarray(sig, dtype=float)) for sig in local_signals])), axis=0)
    control_phase = np.unwrap(np.angle(np.column_stack([np.asarray(coeffs, dtype=complex) for coeffs in control_coeffs])), axis=0)
    canonical_metrics = phase_metrics(canonical_phase, sample_dt)
    control_metrics = phase_metrics(control_phase, sample_dt)

    centers_arr = np.asarray(total_centers, dtype=float)
    center_shift = float(np.linalg.norm(displacement(centers_arr[-1][None, :], centers_arr[0], boundary_type)[0]))
    interaction_summary = {
        'packet_count': int(rep['packet_count']),
        'center_shift': center_shift,
        'width_ratio': float(widths[-1] / max(widths[0], 1.0e-12)),
        'max_mean_overlap': float(np.max(overlap_series)) if overlap_series else 0.0,
        'phase_lock_indicator': float(canonical_metrics['mean_R']),
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

    row = {
        'run_id': run['run_id'],
        'phase': 'base',
        'representative_id': run['representative_id'],
        'condition_id': next(item['condition_id'] for item in load_runsheet(RUNSHEET_PATH)['conditions'] if float(item['phase_coupling_eps']) == feedback_eps),
        'phase_proxy': run['phase_proxy'],
        'control_phase_proxy': run['control_phase_proxy'],
        'phase_coupling_eps': feedback_eps,
        'resolution': resolution,
        'boundary_type': run['boundary_type'],
        'operator_sector': run['operator_sector'],
        'packet_count': int(run['packet_count']),
        'interaction_label': interaction_label,
        'coarse_label': coarse_label,
        'mean_R_canonical': canonical_metrics['mean_R'],
        'median_R_canonical': canonical_metrics['median_R'],
        'R_fraction_canonical': canonical_metrics['R_fraction'],
        'mean_R_control': control_metrics['mean_R'],
        'median_R_control': control_metrics['median_R'],
        'R_fraction_control': control_metrics['R_fraction'],
        'phase_variance_mean_canonical': canonical_metrics['phase_variance_mean'],
        'phase_variance_mean_control': control_metrics['phase_variance_mean'],
        'frequency_spread_mean_canonical': canonical_metrics['frequency_spread_mean'],
        'frequency_spread_mean_control': control_metrics['frequency_spread_mean'],
        'cluster_count_mean_canonical': canonical_metrics['cluster_count_mean'],
        'cluster_count_mean_control': control_metrics['cluster_count_mean'],
        'cluster_dwell_fraction_canonical': canonical_metrics['cluster_dwell_fraction'],
        'cluster_dwell_fraction_control': control_metrics['cluster_dwell_fraction'],
        'mean_R_delta': 0.0,
        'median_R_delta': 0.0,
        'R_fraction_delta': 0.0,
        'phase_variance_delta': 0.0,
        'frequency_spread_delta': 0.0,
        'cluster_dwell_delta': 0.0,
        'composite_lifetime': interaction_summary['composite_lifetime'],
        'binding_persistence': binding_persistence,
        'coarse_persistence': coarse_summary['coarse_persistence_sigma4'],
        'composite_lifetime_delta': 0.0,
        'binding_persistence_delta': 0.0,
        'coarse_persistence_delta': 0.0,
        'morphology_correlate_flag': 0,
        'sync_label': 'no synchronization effect',
        'promoted_followup': 0,
        'notes': run['notes'],
    }

    payload = {
        'run': run,
        'summary': row,
        'times': times,
        'canonical_metrics': canonical_metrics,
        'control_metrics': control_metrics,
        'pair_distance_series': pair_distance_series,
        'overlap_series': overlap_series,
        'coarse_sigma4': {
            'dominant_area_fractions': coarse[4.0]['dominant_area_fractions'],
            'basin_counts': coarse[4.0]['basin_counts'],
            'stable_flags': coarse[4.0]['stable_flags'],
            'last_env': coarse[4.0]['last_env'],
        },
    }

    plot_paths: list[Path] = []
    sync_trace = work_plot_dir / f"stage17_{run['run_id']}_sync_trace.png"
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))
    axes[0].plot(times, canonical_metrics['R_t'], color='tab:blue', label='canonical R')
    axes[0].axhline(R_THRESHOLD, color='black', linestyle='--', linewidth=1.0)
    axes[0].set_title('Canonical R(t)')
    axes[0].set_xlabel('time')
    axes[0].grid(alpha=0.25)
    axes[1].plot(times, [1.0 - value for value in canonical_metrics['R_t']], color='tab:orange', label='canonical variance')
    axes[1].set_title('Canonical phase variance')
    axes[1].set_xlabel('time')
    axes[1].grid(alpha=0.25)
    fig.suptitle(f"Stage 17: {run['run_id']}")
    fig.savefig(sync_trace, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(sync_trace)

    basin_trace = work_plot_dir / f"stage17_{run['run_id']}_basin_trace.png"
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))
    axes[0].plot(times, coarse[4.0]['dominant_area_fractions'], color='tab:green')
    axes[0].set_title('Dominant basin area (sigma=4)')
    axes[0].set_xlabel('time')
    axes[0].grid(alpha=0.25)
    axes[1].plot(times, coarse[4.0]['basin_counts'], color='tab:red')
    axes[1].set_title('Basin count (sigma=4)')
    axes[1].set_xlabel('time')
    axes[1].grid(alpha=0.25)
    fig.savefig(basin_trace, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(basin_trace)

    env = coarse[4.0]['last_env']
    if env is not None:
        env_path = work_plot_dir / f"stage17_{run['run_id']}_envelope_snapshot.png"
        fig, ax = plt.subplots(figsize=(5.4, 4.6))
        slice_idx = env.shape[2] // 2
        img = env[:, :, slice_idx].T
        im = ax.imshow(img, origin='lower', cmap='viridis')
        ax.set_title(f"Envelope sigma=4: {run['run_id']}")
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
        fig.savefig(env_path, dpi=180, bbox_inches='tight')
        plt.close(fig)
        plot_paths.append(env_path)

    return payload, plot_paths


def classify_base_rows(rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, Any] | None]:
    promoted: dict[str, Any] | None = None
    best_score = None
    by_rep = {(row['representative_id'], row['condition_id']): row for row in rows}

    for row in rows:
        if row['condition_id'] == 'null_control':
            continue
        null_row = by_rep[(row['representative_id'], 'null_control')]
        row['mean_R_delta'] = float(row['mean_R_canonical']) - float(null_row['mean_R_canonical'])
        row['median_R_delta'] = float(row['median_R_canonical']) - float(null_row['median_R_canonical'])
        row['R_fraction_delta'] = float(row['R_fraction_canonical']) - float(null_row['R_fraction_canonical'])
        row['phase_variance_delta'] = float(row['phase_variance_mean_canonical']) - float(null_row['phase_variance_mean_canonical'])
        row['frequency_spread_delta'] = float(row['frequency_spread_mean_canonical']) - float(null_row['frequency_spread_mean_canonical'])
        row['cluster_dwell_delta'] = float(row['cluster_dwell_fraction_canonical']) - float(null_row['cluster_dwell_fraction_canonical'])
        row['composite_lifetime_delta'] = float(row['composite_lifetime']) - float(null_row['composite_lifetime'])
        row['binding_persistence_delta'] = float(row['binding_persistence']) - float(null_row['binding_persistence'])
        row['coarse_persistence_delta'] = float(row['coarse_persistence']) - float(null_row['coarse_persistence'])
        row['morphology_correlate_flag'] = int(
            row['composite_lifetime_delta'] > MORPH_DELTA_MIN
            or row['binding_persistence_delta'] > MORPH_DELTA_MIN
            or row['coarse_persistence_delta'] > MORPH_DELTA_MIN
        )

        sustained_gain_count = sum(
            [
                row['mean_R_delta'] > SUSTAINED_GAIN_DELTA,
                row['median_R_delta'] > SUSTAINED_GAIN_DELTA,
                row['R_fraction_delta'] > FRACTION_GAIN_DELTA,
            ]
        )
        canonical_gain = sustained_gain_count >= 2
        variance_ok = row['phase_variance_delta'] <= VARIANCE_DROP_DELTA
        cluster_ok = row['cluster_dwell_delta'] >= CLUSTER_DWELL_GAIN_DELTA
        frequency_ok = row['frequency_spread_delta'] <= FREQ_SPREAD_DROP_DELTA
        control_contradiction = (
            (float(row['mean_R_control']) + 0.03 < float(null_row['mean_R_control']))
            or (float(row['R_fraction_control']) + 0.08 < float(null_row['R_fraction_control']))
        )
        symmetry_case = row['representative_id'] in {'counter_propagating_corridor', 'phase_ordered_symmetric_triad'}

        if not canonical_gain or control_contradiction or not variance_ok or not cluster_ok:
            row['sync_label'] = 'no synchronization effect'
        else:
            if symmetry_case and not (frequency_ok and row['morphology_correlate_flag']):
                row['sync_label'] = 'transient synchronization'
            elif not frequency_ok or not row['morphology_correlate_flag']:
                row['sync_label'] = 'transient synchronization'
            else:
                row['sync_label'] = 'scale-fragile phase locking'
                score = (
                    float(row['mean_R_delta'])
                    - float(row['phase_variance_delta'])
                    - float(row['frequency_spread_delta'])
                    + float(row['cluster_dwell_delta'])
                )
                if best_score is None or score > best_score:
                    best_score = score
                    promoted = {'representative_id': row['representative_id'], 'condition_id': row['condition_id']}

    if promoted is not None:
        for row in rows:
            if row['representative_id'] == promoted['representative_id'] and row['condition_id'] == promoted['condition_id']:
                row['promoted_followup'] = 1
    return rows, promoted


def add_refinement_label(base_row: dict[str, Any], refined_row: dict[str, Any]) -> None:
    refined_row['mean_R_delta'] = float(refined_row['mean_R_canonical']) - float(base_row['mean_R_canonical'])
    refined_row['phase_variance_delta'] = float(refined_row['phase_variance_mean_canonical']) - float(base_row['phase_variance_mean_canonical'])
    refined_row['frequency_spread_delta'] = float(refined_row['frequency_spread_mean_canonical']) - float(base_row['frequency_spread_mean_canonical'])
    refined_row['cluster_dwell_delta'] = float(refined_row['cluster_dwell_fraction_canonical']) - float(base_row['cluster_dwell_fraction_canonical'])
    if (
        refined_row['mean_R_delta'] > -0.02
        and refined_row['phase_variance_delta'] <= 0.02
        and refined_row['frequency_spread_delta'] <= 0.02
        and refined_row['cluster_dwell_delta'] >= -0.02
    ):
        refined_row['sync_label'] = 'scale-stable phase locking'
    else:
        refined_row['sync_label'] = 'scale-fragile phase locking'


def create_summary_plots(work_plot_dir: Path, rows: list[dict[str, Any]], refined_rows: list[dict[str, Any]]) -> list[Path]:
    plot_paths: list[Path] = []
    base_rows = [row for row in rows if row['phase'] == 'base']
    reps = sorted({row['representative_id'] for row in base_rows})
    conds = sorted(
        {row['condition_id'] for row in base_rows},
        key=lambda item: ['null_control', 'very_weak_phase_coupling', 'weak_phase_coupling'].index(item)
        if item in ['null_control', 'very_weak_phase_coupling', 'weak_phase_coupling']
        else 99,
    )
    cond_labels = {
        'null_control': 'null',
        'very_weak_phase_coupling': '0.01',
        'weak_phase_coupling': '0.02',
    }

    comparison = work_plot_dir / 'stage17_phase_sync_comparison.png'
    x = np.arange(len(reps))
    width = 0.8 / max(len(conds), 1)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    for idx, cond in enumerate(conds):
        mean_rs = []
        cluster_dwell = []
        for rep in reps:
            row = next(item for item in base_rows if item['representative_id'] == rep and item['condition_id'] == cond)
            mean_rs.append(float(row['mean_R_canonical']))
            cluster_dwell.append(float(row['cluster_dwell_fraction_canonical']))
        offset = idx - 0.5 * (len(conds) - 1)
        axes[0].bar(x + offset * width, mean_rs, width=width, label=cond_labels.get(cond, cond))
        axes[1].bar(x + offset * width, cluster_dwell, width=width, label=cond_labels.get(cond, cond))
    axes[0].set_title('Mean canonical R')
    axes[1].set_title('Cluster dwell fraction')
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(reps, rotation=20, ha='right')
        ax.grid(alpha=0.25, axis='y')
        ax.legend(fontsize=8)
    fig.suptitle('Stage 17 null vs phase coupling')
    fig.savefig(comparison, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(comparison)

    response = work_plot_dir / 'stage17_phase_sync_response_matrix.png'
    label_map = {
        'no synchronization effect': 0,
        'transient synchronization': 1,
        'scale-fragile phase locking': 2,
        'scale-stable phase locking': 3,
    }
    data = np.zeros((len(reps), len(conds)), dtype=float)
    for i, rep in enumerate(reps):
        for j, cond in enumerate(conds):
            row = next(item for item in base_rows if item['representative_id'] == rep and item['condition_id'] == cond)
            data[i, j] = label_map.get(str(row['sync_label']), 0)
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    im = ax.imshow(data, cmap='viridis', vmin=0, vmax=3)
    ax.set_xticks(range(len(conds)))
    ax.set_xticklabels([cond_labels.get(cond, cond) for cond in conds])
    ax.set_yticks(range(len(reps)))
    ax.set_yticklabels(reps)
    ax.set_title('Stage 17 synchronization response matrix')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(response, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(response)

    if refined_rows:
        ref_plot = work_plot_dir / 'stage17_phase_sync_refinement.png'
        labels = [f"{row['condition_id']}\nn={row['resolution']}" for row in refined_rows]
        values = [float(row['mean_R_canonical']) for row in refined_rows]
        fig, ax = plt.subplots(figsize=(7.2, 4.4))
        ax.bar(range(len(labels)), values, color='tab:purple')
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels)
        ax.set_title('Refined mean canonical R')
        ax.grid(alpha=0.25, axis='y')
        fig.savefig(ref_plot, dpi=180, bbox_inches='tight')
        plt.close(fig)
        plot_paths.append(ref_plot)

    return plot_paths


def write_note(path: Path, json_rel: str, csv_rel: str, base_rows: list[dict[str, Any]], refined_rows: list[dict[str, Any]], label_counts: dict[str, int], promoted: dict[str, Any] | None) -> None:
    lines = [
        '# Stage 17 Emergent Phase Synchronization Dynamics v1',
        '',
        f'Timestamped summary: `{json_rel}`',
        f'Timestamped run table: `{csv_rel}`',
        '',
        f'Synchronization label counts: {label_counts}',
        '',
        'This pilot tests whether weak phase-coupled feedback produces self-sustained collective timing structure without modifying the frozen kernel-induced cochain operator architecture.',
        '',
        'Base comparisons:',
        '',
    ]
    for row in base_rows:
        lines.extend([
            f"- `{row['representative_id']}` / `{row['condition_id']}`",
            f"  - mean_R delta: `{row['mean_R_delta']:.4f}`",
            f"  - phase variance delta: `{row['phase_variance_delta']:.4f}`",
            f"  - frequency spread delta: `{row['frequency_spread_delta']:.4f}`",
            f"  - cluster dwell delta: `{row['cluster_dwell_delta']:.4f}`",
            f"  - morphology correlate: `{row['morphology_correlate_flag']}`",
            f"  - label: `{row['sync_label']}`",
        ])
    lines.extend(['', f'Promoted follow-up: `{promoted}`' if promoted else 'Promoted follow-up: `none`', ''])
    if refined_rows:
        lines.append('Refinement follow-up:')
        lines.append('')
        for row in refined_rows:
            lines.extend([
                f"- `{row['representative_id']}` / `{row['condition_id']}` / `n={row['resolution']}`",
                f"  - label: `{row['sync_label']}`",
            ])
        lines.append('')
    lines.extend([
        'Interpretation boundary:',
        '- synchronization is read only as a collective ordering precursor',
        '- no binding, spacetime, or force claim is made here',
        '',
    ])
    path.write_text('\n'.join(lines), encoding='utf-8')


def main() -> None:
    args = parse_args()
    defaults = load_stage10_defaults()
    runsheet = load_runsheet(args.runsheet)
    rep_lookup = representative_lookup(runsheet)
    base = runsheet['base_seed_reference']
    runs = selected_runs(runsheet['runs'], args.run_ids)
    work_plot_dir = Path('/tmp') / 'haos_iip_stage17_phase_sync'
    work_plot_dir.mkdir(parents=True, exist_ok=True)

    base_payloads: dict[tuple[str, str], dict[str, Any]] = {}
    base_rows: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    for run in runs:
        if int(run['resolution']) != int(runsheet['resolution_policy']['base_resolution']):
            continue
        rep = rep_lookup[run['representative_id']]
        payload, run_plots = simulate_run(run, rep, defaults, base, work_plot_dir)
        base_payloads[(run['representative_id'], payload['summary']['condition_id'])] = payload
        base_rows.append(payload['summary'])
        plot_paths.extend(run_plots)

    base_rows, promoted = classify_base_rows(base_rows)
    refined_rows: list[dict[str, Any]] = []
    if promoted is not None and not args.skip_followup:
        chosen_base = next(row for row in base_rows if row['representative_id'] == promoted['representative_id'] and row['condition_id'] == promoted['condition_id'])
        chosen_null = next(row for row in base_rows if row['representative_id'] == promoted['representative_id'] and row['condition_id'] == 'null_control')
        rep = rep_lookup[promoted['representative_id']]
        followups = [
            dict(next(run for run in runs if run['representative_id'] == promoted['representative_id'] and run['condition_id'] == 'null_control'), resolution=int(runsheet['resolution_policy']['refined_resolution']), run_id=f"S17_{promoted['representative_id']}__null_control__n24"),
            dict(next(run for run in runs if run['representative_id'] == promoted['representative_id'] and run['condition_id'] == promoted['condition_id']), resolution=int(runsheet['resolution_policy']['refined_resolution']), run_id=f"S17_{promoted['representative_id']}__{promoted['condition_id']}__n24"),
        ]
        for follow in followups:
            payload, run_plots = simulate_run(follow, rep, defaults, base, work_plot_dir)
            payload['summary']['phase'] = 'refinement'
            payload['summary']['mean_R_delta'] = float(payload['summary']['mean_R_canonical']) - float(chosen_null['mean_R_canonical'])
            payload['summary']['median_R_delta'] = float(payload['summary']['median_R_canonical']) - float(chosen_null['median_R_canonical'])
            payload['summary']['R_fraction_delta'] = float(payload['summary']['R_fraction_canonical']) - float(chosen_null['R_fraction_canonical'])
            payload['summary']['phase_variance_delta'] = float(payload['summary']['phase_variance_mean_canonical']) - float(chosen_null['phase_variance_mean_canonical'])
            payload['summary']['frequency_spread_delta'] = float(payload['summary']['frequency_spread_mean_canonical']) - float(chosen_null['frequency_spread_mean_canonical'])
            payload['summary']['cluster_dwell_delta'] = float(payload['summary']['cluster_dwell_fraction_canonical']) - float(chosen_null['cluster_dwell_fraction_canonical'])
            payload['summary']['composite_lifetime_delta'] = float(payload['summary']['composite_lifetime']) - float(chosen_null['composite_lifetime'])
            payload['summary']['binding_persistence_delta'] = float(payload['summary']['binding_persistence']) - float(chosen_null['binding_persistence'])
            payload['summary']['coarse_persistence_delta'] = float(payload['summary']['coarse_persistence']) - float(chosen_null['coarse_persistence'])
            payload['summary']['morphology_correlate_flag'] = int(
                payload['summary']['composite_lifetime_delta'] > MORPH_DELTA_MIN
                or payload['summary']['binding_persistence_delta'] > MORPH_DELTA_MIN
                or payload['summary']['coarse_persistence_delta'] > MORPH_DELTA_MIN
            )
            if payload['summary']['condition_id'] == promoted['condition_id']:
                add_refinement_label(chosen_base, payload['summary'])
            else:
                payload['summary']['sync_label'] = 'no synchronization effect'
            refined_rows.append(payload['summary'])
            plot_paths.extend(run_plots)

    all_rows = base_rows + refined_rows
    label_counts = dict(Counter(row['sync_label'] for row in all_rows))
    plot_paths.extend(create_summary_plots(work_plot_dir, all_rows, refined_rows))

    result = {
        'stage': 'Stage 17',
        'description': runsheet['description'],
        'base_seed_reference': base,
        'runs': [
            {'run_id': row['run_id'], 'summary': row}
            for row in all_rows
        ],
        'synchronization_label_counts': label_counts,
        'promoted_followup': promoted,
        'system_level_interpretation': (
            'no phase channel'
            if label_counts.get('no synchronization effect', 0) == len(all_rows)
            else 'localized locking only'
        ),
        'conclusion': 'Stage 17 tests whether weak internally derived phase coupling can generate self-sustained synchronization before any binding or geometric claim is made.',
    }

    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug='stage17_phase_sync',
        result=result,
        csv_rows=all_rows,
        csv_fieldnames=CSV_FIELDS,
        plot_paths=plot_paths,
    )
    note_path = ATLAS_NOTES / runsheet['note_name']
    write_note(note_path, str(json_path.relative_to(REPO_ROOT)), str(csv_path.relative_to(REPO_ROOT)), base_rows, refined_rows, label_counts, promoted)
    append_log(
        title=f"Stage 17 Emergent Phase Synchronization ({json_path.stem})",
        config_summary=f"runs={len(base_rows)}, promoted_followup={promoted}",
        result_path=json_path,
        stamped_plots=stamped_plots,
        observation=f"sync_label_counts={label_counts}",
        conclusion='the Stage 17 pilot tests whether weak phase coupling can produce self-sustained collective timing structure on the frozen architecture',
    )

    print(f'wrote {json_path}')
    print(f'wrote {csv_path}')
    for item in stamped_plots:
        print(item)
    print(f'sync_label_counts={label_counts}')
    print(f'promoted_followup={promoted}')


if __name__ == '__main__':
    main()
