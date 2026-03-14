#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from itertools import combinations, permutations
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import (
    ATLAS_NOTES,
    REPO_ROOT,
    append_log,
    build_transverse_setup,
    displacement,
    load_stage10_defaults,
    plt,
    save_atlas_payload,
    suggested_dt,
    weighted_center,
)
from stage11_collective_wave_interaction import build_component_packet, make_leakage_fn
from stage12_coarse_field_structure import edge_field_to_grid, gaussian_smooth_periodic

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage15_relational_geometry_runs.json'
ANALYSIS_STRIDE = 4
ARRIVAL_SMOOTH_SIGMA = 2.0
SETUP_CACHE: dict[tuple[int, str, float, float], tuple[Any, Any]] = {}
CSV_FIELDS = [
    'case_id',
    'distance_id',
    'n_side',
    'mean_absolute_asymmetry',
    'max_absolute_asymmetry',
    'triangle_violation_rate',
    'mean_violation_magnitude',
    'ordering_flip_fraction',
    'pair_rank_correlation',
    'effective_dimension',
    'dimension_drift',
    'shell_distortion_metric',
    'ordering_breakdown_frequency',
    'loop_inconsistency_index',
    'geodesic_deviation_analog',
    'relational_ordering_label',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage 15 emergent relational geometry pilot.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--case-ids', nargs='*', default=[])
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def selected_cases(cases: list[dict[str, Any]], case_ids: list[str]) -> list[dict[str, Any]]:
    if not case_ids:
        return cases
    wanted = set(case_ids)
    return [case for case in cases if case['case_id'] in wanted]


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


def pair_ids(packet_count: int) -> list[tuple[int, int]]:
    return list(combinations(range(packet_count), 2))


def pair_distance_matrix(centers: list[np.ndarray], boundary_type: str) -> np.ndarray:
    n = len(centers)
    mat = np.zeros((n, n), dtype=float)
    for i, j in combinations(range(n), 2):
        delta = displacement(np.asarray([centers[j]], dtype=float), np.asarray(centers[i], dtype=float), boundary_type)[0]
        value = float(np.linalg.norm(delta))
        mat[i, j] = value
        mat[j, i] = value
    return mat


def overlap_matrix(component_qs: list[np.ndarray]) -> np.ndarray:
    n = len(component_qs)
    mat = np.eye(n, dtype=float)
    weights = [np.abs(q) ** 2 for q in component_qs]
    for i, j in combinations(range(n), 2):
        wi = weights[i]
        wj = weights[j]
        denom = math.sqrt(max(float(np.sum(wi)) * float(np.sum(wj)), 1.0e-12))
        value = float(np.sum(np.sqrt(wi * wj)) / denom)
        mat[i, j] = value
        mat[j, i] = value
    return mat


def phase_lock_matrix(component_qs: list[np.ndarray]) -> np.ndarray:
    n = len(component_qs)
    mat = np.eye(n, dtype=float)
    for i, j in combinations(range(n), 2):
        qi = component_qs[i]
        qj = component_qs[j]
        denom = max(float(np.linalg.norm(qi)) * float(np.linalg.norm(qj)), 1.0e-12)
        value = abs(float(np.vdot(qi, qj).real / denom))
        mat[i, j] = value
        mat[j, i] = value
    return mat


def sample_component_envelopes(midpoints: np.ndarray, component_qs: list[np.ndarray], n_side: int, packet_centers: list[list[float]]) -> np.ndarray:
    coords = np.floor(np.asarray(packet_centers, dtype=float) * n_side).astype(int) % n_side
    n = len(component_qs)
    samples = np.zeros((n, n), dtype=float)
    for i, q in enumerate(component_qs):
        grid = edge_field_to_grid(midpoints, np.abs(q), n_side)
        env = np.maximum(gaussian_smooth_periodic(grid, ARRIVAL_SMOOTH_SIGMA), 0.0)
        for j, coord in enumerate(coords):
            samples[i, j] = float(env[coord[0], coord[1], coord[2]])
    return samples


def first_arrival_time(series: list[float], times: list[float], threshold: float) -> float:
    baseline = float(series[0]) if series else 0.0
    peak = float(max(series)) if series else baseline
    if peak <= baseline + 1.0e-12:
        return float(times[-1]) if times else 0.0
    target = baseline + threshold * (peak - baseline)
    for t, value in zip(times, series):
        if float(value) >= target:
            return float(t)
    return float(times[-1]) if times else 0.0


def matrix_mean_absolute_asymmetry(mat: np.ndarray) -> tuple[float, float]:
    diffs = []
    for i, j in combinations(range(mat.shape[0]), 2):
        diffs.append(abs(float(mat[i, j] - mat[j, i])))
    if not diffs:
        return 0.0, 0.0
    return float(np.mean(diffs)), float(np.max(diffs))


def triangle_stats(sym_mat: np.ndarray) -> tuple[float | None, float | None, float | None]:
    n = sym_mat.shape[0]
    if n < 3:
        return None, None, None
    violations: list[float] = []
    for i, j, k in combinations(range(n), 3):
        tri = [
            max(0.0, float(sym_mat[i, k] - (sym_mat[i, j] + sym_mat[j, k]))),
            max(0.0, float(sym_mat[i, j] - (sym_mat[i, k] + sym_mat[k, j]))),
            max(0.0, float(sym_mat[j, k] - (sym_mat[j, i] + sym_mat[i, k]))),
        ]
        violations.extend(tri)
    if not violations:
        return None, None, None
    vals = np.asarray(violations, dtype=float)
    return float(np.mean(vals > 1.0e-9)), float(np.mean(vals)), float(np.max(vals))


def vectorize_upper(sym_mat: np.ndarray) -> tuple[list[tuple[int, int]], np.ndarray]:
    ids = pair_ids(sym_mat.shape[0])
    values = np.asarray([float(sym_mat[i, j]) for i, j in ids], dtype=float)
    return ids, values


def average_ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind='mergesort')
    ranks = np.empty(len(values), dtype=float)
    sorted_vals = values[order]
    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values) and abs(float(sorted_vals[j] - sorted_vals[i])) <= 1.0e-12:
            j += 1
        avg_rank = 0.5 * (i + j - 1)
        ranks[order[i:j]] = avg_rank
        i = j
    return ranks


def ordering_persistence(base_sym: np.ndarray, refined_sym: np.ndarray) -> tuple[float, float]:
    ids, base_vals = vectorize_upper(base_sym)
    _ids, refined_vals = vectorize_upper(refined_sym)
    if len(ids) <= 1:
        return 0.0, 1.0
    flips = 0
    total = 0
    for a, b in combinations(range(len(ids)), 2):
        sign_base = np.sign(base_vals[a] - base_vals[b])
        sign_ref = np.sign(refined_vals[a] - refined_vals[b])
        if sign_base != sign_ref:
            flips += 1
        total += 1
    base_ranks = average_ranks(base_vals)
    refined_ranks = average_ranks(refined_vals)
    if np.std(base_ranks) < 1.0e-12 or np.std(refined_ranks) < 1.0e-12:
        corr = 1.0
    else:
        corr = float(np.corrcoef(base_ranks, refined_ranks)[0, 1])
    return float(flips / max(total, 1)), corr


def ball_growth(distance_sym: np.ndarray) -> tuple[list[float] | None, list[float] | None, float | None]:
    n = distance_sym.shape[0]
    if n < 3:
        return None, None, None
    radii_union = sorted({float(distance_sym[i, j]) for i, j in combinations(range(n), 2) if distance_sym[i, j] > 1.0e-12})
    if len(radii_union) < 2:
        return None, None, None
    mean_volumes: list[float] = []
    for r in radii_union:
        counts = []
        for seed in range(n):
            count = int(np.sum(distance_sym[seed] <= r + 1.0e-12))
            counts.append(float(count))
        mean_volumes.append(float(np.mean(counts)))
    x = np.log(np.asarray(radii_union, dtype=float) + 1.0e-12)
    y = np.log(np.asarray(mean_volumes, dtype=float) + 1.0e-12)
    if len(x) < 2:
        return radii_union, mean_volumes, None
    slope = float(np.polyfit(x, y, 1)[0])
    return radii_union, mean_volumes, slope


def transport_shell_metrics(raw_delay: np.ndarray) -> dict[str, float]:
    n = raw_delay.shape[0]
    distortions: list[float] = []
    for seed in range(n):
        arrivals = [float(raw_delay[seed, j]) for j in range(n) if j != seed]
        if not arrivals:
            continue
        arrivals = sorted(arrivals)
        mean_val = float(np.mean(arrivals))
        distortions.append(float(np.std(arrivals) / max(mean_val, 1.0e-12)))
    if not distortions:
        return {'shell_distortion_metric': 0.0}
    return {'shell_distortion_metric': float(np.mean(distortions))}


def transport_order_breakdown(base_raw: np.ndarray, refined_raw: np.ndarray) -> float:
    n = base_raw.shape[0]
    if n <= 2:
        return 0.0
    changes = 0
    total = 0
    for seed in range(n):
        base_order = tuple(j for j in sorted((j for j in range(n) if j != seed), key=lambda j: base_raw[seed, j]))
        refined_order = tuple(j for j in sorted((j for j in range(n) if j != seed), key=lambda j: refined_raw[seed, j]))
        if base_order != refined_order:
            changes += 1
        total += 1
    return float(changes / max(total, 1))


def loop_defect_stats(sym_mat: np.ndarray) -> tuple[float | None, float | None]:
    n = sym_mat.shape[0]
    if n < 3:
        return None, None
    deviations: list[float] = []
    for i, j, k in permutations(range(n), 3):
        direct = float(sym_mat[i, k])
        indirect = float(sym_mat[i, j] + sym_mat[j, k])
        scale = max(max(direct, indirect), 1.0e-12)
        deviations.append(abs(indirect - direct) / scale)
    if not deviations:
        return None, None
    vals = np.asarray(deviations, dtype=float)
    return float(np.mean(vals)), float(np.max(vals))


def classify_ordering(summary: dict[str, Any], rules: dict[str, Any]) -> str:
    ordering_flip = float(summary['ordering_flip_fraction'])
    symmetry = float(summary['mean_absolute_asymmetry'])
    triangle = summary['triangle_violation_rate']
    shell_breakdown = summary['ordering_breakdown_frequency']
    dimension_drift = summary['dimension_drift']

    if (
        ordering_flip >= float(rules['no_coherent_relational_ordering']['ordering_flip_fraction_min'])
        or symmetry >= float(rules['no_coherent_relational_ordering']['symmetry_mean_min'])
        or (triangle is not None and triangle >= float(rules['no_coherent_relational_ordering']['triangle_violation_rate_min']))
    ):
        return 'no coherent relational ordering'

    stable_ok = (
        ordering_flip <= float(rules['scale_stable_metric_like_ordering']['ordering_flip_fraction_max'])
        and symmetry <= float(rules['scale_stable_metric_like_ordering']['symmetry_mean_max'])
        and ((triangle is None) or (triangle <= float(rules['scale_stable_metric_like_ordering']['triangle_violation_rate_max'])))
        and ((shell_breakdown is None) or (shell_breakdown <= float(rules['scale_stable_metric_like_ordering']['shell_breakdown_frequency_max'])))
        and ((dimension_drift is None) or (dimension_drift <= float(rules['scale_stable_metric_like_ordering']['dimension_drift_max'])))
    )
    if stable_ok:
        return 'scale-stable metric-like ordering'

    fragile_ok = (
        ordering_flip <= float(rules['scale_fragile_metric_like_ordering']['ordering_flip_fraction_max'])
        and symmetry <= float(rules['scale_fragile_metric_like_ordering']['symmetry_mean_max'])
        and ((triangle is None) or (triangle <= float(rules['scale_fragile_metric_like_ordering']['triangle_violation_rate_max'])))
        and ((shell_breakdown is None) or (shell_breakdown <= float(rules['scale_fragile_metric_like_ordering']['shell_breakdown_frequency_max'])))
        and ((dimension_drift is None) or (dimension_drift <= float(rules['scale_fragile_metric_like_ordering']['dimension_drift_max'])))
    )
    if fragile_ok:
        return 'scale-fragile metric-like ordering'

    return 'weak ordering with high metric violation'


def render_matrix(ax, mat: np.ndarray, title: str) -> None:
    im = ax.imshow(mat, cmap='viridis')
    ax.set_title(title, fontsize=9)
    ax.set_xticks(range(mat.shape[0]))
    ax.set_yticks(range(mat.shape[0]))
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


def simulate_case(case: dict[str, Any], n_side: int, defaults: dict[str, Any], base: dict[str, Any]) -> dict[str, Any]:
    data, projector = get_cached_setup(n_side=n_side, defaults=defaults, boundary_type=str(base['boundary_type']))
    kick_axis = int(defaults['kick_axis'])
    packet_qs: list[np.ndarray] = []
    packet_vs: list[np.ndarray] = []
    for center, amp, phase, kick_sign in zip(case['packet_centers'], case['packet_amplitudes'], case['phase_offsets_rad'], case['kick_signs']):
        q0, v0 = build_component_packet(
            center=np.asarray(center, dtype=float),
            amplitude=float(amp),
            phase_offset=float(phase),
            kick_sign=float(kick_sign),
            bandwidth=float(base['bandwidth']),
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
    overlap_series: dict[tuple[int, int], list[float]] = {pid: [] for pid in pair_ids(int(case['packet_count']))}
    phase_series: dict[tuple[int, int], list[float]] = {pid: [] for pid in pair_ids(int(case['packet_count']))}
    transport_series: dict[tuple[int, int], list[float]] = {
        (i, j): [] for i in range(int(case['packet_count'])) for j in range(int(case['packet_count'])) if i != j
    }
    constraint_norms: list[float] = []
    sector_leakages: list[float] = []

    def record(t: float) -> None:
        times.append(float(t))
        ov = overlap_matrix(packet_qs)
        ph = phase_lock_matrix(packet_qs)
        env_samples = sample_component_envelopes(data.midpoints, packet_qs, n_side, case['packet_centers'])
        for i, j in pair_ids(int(case['packet_count'])):
            overlap_series[(i, j)].append(float(ov[i, j]))
            phase_series[(i, j)].append(float(ph[i, j]))
        for i in range(int(case['packet_count'])):
            for j in range(int(case['packet_count'])):
                if i == j:
                    continue
                transport_series[(i, j)].append(float(env_samples[i, j]))
        total_q = np.sum(packet_qs, axis=0)
        constraint_norms.append(float(np.linalg.norm(data.d0.T @ total_q)))
        sector_leakages.append(float(leakage_fn(total_q)))

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

    packet_count = int(case['packet_count'])
    persistence_sym = np.eye(packet_count, dtype=float)
    phase_sym = np.eye(packet_count, dtype=float)
    transport_raw = np.zeros((packet_count, packet_count), dtype=float)

    for i, j in pair_ids(packet_count):
        mean_overlap = float(np.mean(overlap_series[(i, j)]))
        persistence_sym[i, j] = 1.0 / (0.001 + mean_overlap)
        persistence_sym[j, i] = persistence_sym[i, j]
        mean_phase = float(np.mean(np.abs(phase_series[(i, j)])))
        phase_sym[i, j] = 1.0 - mean_phase
        phase_sym[j, i] = phase_sym[i, j]
    for i in range(packet_count):
        for j in range(packet_count):
            if i == j:
                continue
            transport_raw[i, j] = first_arrival_time(transport_series[(i, j)], times, 0.15)
    transport_sym = 0.5 * (transport_raw + transport_raw.T)

    distance_payloads: dict[str, dict[str, Any]] = {}
    for distance_id, raw_mat, sym_mat in [
        ('persistence_distance', persistence_sym, persistence_sym),
        ('phase_correlation_distance', phase_sym, phase_sym),
        ('transport_delay_distance', transport_raw, transport_sym),
    ]:
        mean_asym, max_asym = matrix_mean_absolute_asymmetry(raw_mat)
        tri_rate, tri_mean, tri_max = triangle_stats(sym_mat)
        radii, volumes, eff_dim = (None, None, None)
        if distance_id == 'persistence_distance' and bool(case['supports_ball_growth']):
            radii, volumes, eff_dim = ball_growth(sym_mat)
        loop_idx, geodesic_dev = (None, None)
        if bool(case['supports_triangle_tests']):
            loop_idx, geodesic_dev = loop_defect_stats(sym_mat)
        shell_metrics = {'shell_distortion_metric': None}
        if distance_id == 'transport_delay_distance':
            shell_metrics = transport_shell_metrics(raw_mat)
        distance_payloads[distance_id] = {
            'raw_matrix': raw_mat.tolist(),
            'sym_matrix': sym_mat.tolist(),
            'mean_absolute_asymmetry': mean_asym,
            'max_absolute_asymmetry': max_asym,
            'triangle_violation_rate': tri_rate,
            'mean_violation_magnitude': tri_mean,
            'max_violation_magnitude': tri_max,
            'effective_dimension': eff_dim,
            'ball_growth_radii': radii,
            'ball_growth_volumes': volumes,
            'shell_distortion_metric': shell_metrics['shell_distortion_metric'],
            'loop_inconsistency_index': loop_idx,
            'geodesic_deviation_analog': geodesic_dev,
        }

    return {
        'case_id': case['case_id'],
        'n_side': int(n_side),
        'times': times,
        'distance_payloads': distance_payloads,
        'constraint_max': float(np.max(constraint_norms)) if constraint_norms else 0.0,
        'sector_leakage': float(np.max(sector_leakages)) if sector_leakages else 0.0,
    }


def build_rows(case: dict[str, Any], base_result: dict[str, Any], refined_result: dict[str, Any], runsheet: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    summaries: dict[str, Any] = {}
    order_breakdown_transport = transport_order_breakdown(
        np.asarray(base_result['distance_payloads']['transport_delay_distance']['raw_matrix'], dtype=float),
        np.asarray(refined_result['distance_payloads']['transport_delay_distance']['raw_matrix'], dtype=float),
    )
    for distance in runsheet['distance_candidates']:
        distance_id = distance['distance_id']
        base_payload = base_result['distance_payloads'][distance_id]
        refined_payload = refined_result['distance_payloads'][distance_id]
        ordering_flip, rank_corr = ordering_persistence(
            np.asarray(base_payload['sym_matrix'], dtype=float),
            np.asarray(refined_payload['sym_matrix'], dtype=float),
        )
        dimension_drift = None
        if base_payload['effective_dimension'] is not None and refined_payload['effective_dimension'] is not None:
            dimension_drift = abs(float(refined_payload['effective_dimension']) - float(base_payload['effective_dimension']))
        summary = {
            'ordering_flip_fraction': ordering_flip,
            'pair_rank_correlation': rank_corr,
            'mean_absolute_asymmetry': max(float(base_payload['mean_absolute_asymmetry']), float(refined_payload['mean_absolute_asymmetry'])),
            'max_absolute_asymmetry': max(float(base_payload['max_absolute_asymmetry']), float(refined_payload['max_absolute_asymmetry'])),
            'triangle_violation_rate': max(
                [v for v in [base_payload['triangle_violation_rate'], refined_payload['triangle_violation_rate']] if v is not None],
                default=None,
            ),
            'mean_violation_magnitude': max(
                [v for v in [base_payload['mean_violation_magnitude'], refined_payload['mean_violation_magnitude']] if v is not None],
                default=None,
            ),
            'effective_dimension_base': base_payload['effective_dimension'],
            'effective_dimension_refined': refined_payload['effective_dimension'],
            'dimension_drift': dimension_drift,
            'shell_distortion_metric': max(
                [v for v in [base_payload['shell_distortion_metric'], refined_payload['shell_distortion_metric']] if v is not None],
                default=None,
            ),
            'ordering_breakdown_frequency': order_breakdown_transport if distance_id == 'transport_delay_distance' else None,
            'loop_inconsistency_index': max(
                [v for v in [base_payload['loop_inconsistency_index'], refined_payload['loop_inconsistency_index']] if v is not None],
                default=None,
            ),
            'geodesic_deviation_analog': max(
                [v for v in [base_payload['geodesic_deviation_analog'], refined_payload['geodesic_deviation_analog']] if v is not None],
                default=None,
            ),
        }
        summary['relational_ordering_label'] = classify_ordering(summary, runsheet['classification_rules'])
        summaries[distance_id] = summary
        for result, payload in [(base_result, base_payload), (refined_result, refined_payload)]:
            rows.append({
                'case_id': case['case_id'],
                'distance_id': distance_id,
                'n_side': int(result['n_side']),
                'mean_absolute_asymmetry': payload['mean_absolute_asymmetry'],
                'max_absolute_asymmetry': payload['max_absolute_asymmetry'],
                'triangle_violation_rate': payload['triangle_violation_rate'],
                'mean_violation_magnitude': payload['mean_violation_magnitude'],
                'ordering_flip_fraction': summary['ordering_flip_fraction'],
                'pair_rank_correlation': summary['pair_rank_correlation'],
                'effective_dimension': payload['effective_dimension'],
                'dimension_drift': summary['dimension_drift'],
                'shell_distortion_metric': payload['shell_distortion_metric'],
                'ordering_breakdown_frequency': summary['ordering_breakdown_frequency'],
                'loop_inconsistency_index': payload['loop_inconsistency_index'],
                'geodesic_deviation_analog': payload['geodesic_deviation_analog'],
                'relational_ordering_label': summary['relational_ordering_label'],
            })
    return rows, summaries


def create_distance_matrix_grid(path: Path, case_results: dict[str, dict[str, dict[str, Any]]], runsheet: dict[str, Any]) -> None:
    cases = list(case_results.keys())
    distances = [item['distance_id'] for item in runsheet['distance_candidates']]
    fig, axes = plt.subplots(len(cases), len(distances) * 2, figsize=(3.2 * len(distances) * 2, 3.0 * len(cases)))
    if len(cases) == 1:
        axes = np.asarray([axes])
    for row_idx, case_id in enumerate(cases):
        for dist_idx, distance_id in enumerate(distances):
            for col_offset, res_key in enumerate(['base', 'refined']):
                ax = axes[row_idx, dist_idx * 2 + col_offset]
                mat = np.asarray(case_results[case_id][res_key]['distance_payloads'][distance_id]['sym_matrix'], dtype=float)
                title = f"{case_id}\n{distance_id}\nn={case_results[case_id][res_key]['n_side']}"
                render_matrix(ax, mat, title)
    fig.suptitle('Stage 15 relational distance matrices')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_metric_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [f"{row['case_id']}\n{row['distance_id']}\nn={row['n_side']}" for row in rows]
    x = np.arange(len(labels))
    asym = [float(row['mean_absolute_asymmetry']) for row in rows]
    tri = [0.0 if row['triangle_violation_rate'] is None else float(row['triangle_violation_rate']) for row in rows]
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.2))
    axes[0].bar(x, asym, color='tab:blue')
    axes[0].set_title('Mean symmetry deviation')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=70, ha='right', fontsize=7)
    axes[0].grid(alpha=0.25, axis='y')
    axes[1].bar(x, tri, color='tab:red')
    axes[1].set_title('Triangle violation rate')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=70, ha='right', fontsize=7)
    axes[1].grid(alpha=0.25, axis='y')
    fig.suptitle('Stage 15 metric-violation summary')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_ordering_flip_heatmap(path: Path, summary_rows: list[dict[str, Any]], runsheet: dict[str, Any]) -> None:
    cases = [case['case_id'] for case in runsheet['cases']]
    distances = [item['distance_id'] for item in runsheet['distance_candidates']]
    data = np.zeros((len(cases), len(distances)), dtype=float)
    for i, case_id in enumerate(cases):
        for j, distance_id in enumerate(distances):
            row = next(item for item in summary_rows if item['case_id'] == case_id and item['distance_id'] == distance_id)
            data[i, j] = float(row['ordering_flip_fraction'])
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    im = ax.imshow(data, cmap='magma', vmin=0.0, vmax=1.0)
    ax.set_xticks(range(len(distances)))
    ax.set_xticklabels(distances, rotation=20, ha='right')
    ax.set_yticks(range(len(cases)))
    ax.set_yticklabels(cases)
    ax.set_title('Stage 15 ordering flip fraction')
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_dimension_plot(path: Path, case_results: dict[str, dict[str, dict[str, Any]]]) -> None:
    triad_cases = [case_id for case_id in case_results if case_results[case_id]['base']['distance_payloads']['persistence_distance']['ball_growth_radii']]
    fig, axes = plt.subplots(1, max(len(triad_cases), 1), figsize=(6.0 * max(len(triad_cases), 1), 4.5))
    if not isinstance(axes, np.ndarray):
        axes = np.asarray([axes])
    if not triad_cases:
        axes[0].text(0.5, 0.5, 'No triadic ball-growth data', ha='center', va='center')
        axes[0].axis('off')
    else:
        for ax, case_id in zip(axes, triad_cases):
            for label, color, res_key in [('base', 'tab:blue', 'base'), ('refined', 'tab:orange', 'refined')]:
                payload = case_results[case_id][res_key]['distance_payloads']['persistence_distance']
                radii = payload['ball_growth_radii']
                volumes = payload['ball_growth_volumes']
                if radii and volumes:
                    ax.plot(radii, volumes, marker='o', color=color, label=f"{label} (D={payload['effective_dimension']:.2f})" if payload['effective_dimension'] is not None else label)
            ax.set_title(case_id)
            ax.set_xlabel('relational radius')
            ax.set_ylabel('mean ball volume')
            ax.grid(alpha=0.25)
            ax.legend(fontsize=8)
    fig.suptitle('Stage 15 effective dimensionality growth')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_transport_shells(path: Path, case_results: dict[str, dict[str, dict[str, Any]]]) -> None:
    cases = list(case_results.keys())
    fig, axes = plt.subplots(1, len(cases), figsize=(5.4 * len(cases), 4.5))
    if not isinstance(axes, np.ndarray):
        axes = np.asarray([axes])
    for ax, case_id in zip(axes, cases):
        for label, color, res_key in [('base', 'tab:green', 'base'), ('refined', 'tab:purple', 'refined')]:
            raw = np.asarray(case_results[case_id][res_key]['distance_payloads']['transport_delay_distance']['raw_matrix'], dtype=float)
            arrivals = sorted(float(raw[0, j]) for j in range(raw.shape[0]) if j != 0)
            ax.plot(range(1, len(arrivals) + 1), arrivals, marker='o', color=color, label=label)
        ax.set_title(case_id)
        ax.set_xlabel('shell index from packet 0')
        ax.set_ylabel('arrival time')
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    fig.suptitle('Stage 15 transport-delay shells')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_loop_plot(path: Path, summary_rows: list[dict[str, Any]]) -> None:
    rows = [row for row in summary_rows if row['loop_inconsistency_index'] is not None]
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8))
    if not rows:
        for ax in axes:
            ax.text(0.5, 0.5, 'No loop-defect data', ha='center', va='center')
            ax.axis('off')
    else:
        labels = [f"{row['case_id']}\n{row['distance_id']}" for row in rows]
        x = np.arange(len(labels))
        loop_vals = [float(row['loop_inconsistency_index']) for row in rows]
        geo_vals = [float(row['geodesic_deviation_analog']) for row in rows]
        axes[0].bar(x, loop_vals, color='tab:brown')
        axes[0].set_title('Loop inconsistency index')
        axes[0].set_xticks(x)
        axes[0].set_xticklabels(labels, rotation=60, ha='right', fontsize=8)
        axes[0].grid(alpha=0.25, axis='y')
        axes[1].bar(x, geo_vals, color='tab:gray')
        axes[1].set_title('Geodesic deviation analog')
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(labels, rotation=60, ha='right', fontsize=8)
        axes[1].grid(alpha=0.25, axis='y')
    fig.suptitle('Stage 15 loop-defect summary')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(path: Path, json_rel: str, csv_rel: str, summary_rows: list[dict[str, Any]], label_counts: dict[str, int]) -> None:
    lines = [
        '# Stage 15 Emergent Relational Geometry v1',
        '',
        'This note records the first executable Stage 15 pilot run.',
        '',
        f'Timestamped summary: `{json_rel}`',
        f'Timestamped run table: `{csv_rel}`',
        '',
        f'Relational ordering counts: {label_counts}',
        '',
        'Pilot scope:',
        '- 3 cases',
        '- 2 resolutions (`12 -> 24`)',
        '- 3 induced distance candidates',
        '',
        'Per case / distance summary:',
        '',
    ]
    for row in summary_rows:
        lines.extend([
            f"- `{row['case_id']}` / `{row['distance_id']}`",
            f"  - label: `{row['relational_ordering_label']}`",
            f"  - ordering flip fraction: `{row['ordering_flip_fraction']:.4f}`",
            f"  - pair-rank correlation: `{row['pair_rank_correlation']:.4f}`",
            f"  - symmetry deviation: `{row['mean_absolute_asymmetry']:.4f}`",
            f"  - triangle violation rate: `{0.0 if row['triangle_violation_rate'] is None else row['triangle_violation_rate']:.4f}`",
        ])
        if row['effective_dimension'] is not None:
            lines.append(f"  - refined effective dimension: `{row['effective_dimension']:.4f}`")
        if row['dimension_drift'] is not None:
            lines.append(f"  - dimension drift: `{row['dimension_drift']:.4f}`")
        if row['ordering_breakdown_frequency'] is not None:
            lines.append(f"  - transport shell breakdown: `{row['ordering_breakdown_frequency']:.4f}`")
    lines.extend([
        '',
        'Interpretation boundary:',
        '- relational ordering is treated only as diagnostic bookkeeping',
        '- no spacetime, metric, or gravitational claims are made here',
    ])
    path.write_text('\n'.join(lines), encoding='utf-8')


def main() -> None:
    args = parse_args()
    defaults = load_stage10_defaults()
    runsheet = load_runsheet(args.runsheet)
    cases = selected_cases(runsheet['cases'], args.case_ids)
    base = runsheet['base_seed_reference']
    base_n = int(runsheet['resolution_policy']['base_resolution'])
    refined_n = int(runsheet['resolution_policy']['refined_resolution'])

    case_results: dict[str, dict[str, dict[str, Any]]] = {}
    csv_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []

    for case in cases:
        base_result = simulate_case(case, base_n, defaults, base)
        refined_result = simulate_case(case, refined_n, defaults, base)
        case_results[case['case_id']] = {'base': base_result, 'refined': refined_result}
        rows, summaries = build_rows(case, base_result, refined_result, runsheet)
        csv_rows.extend(rows)
        for distance_id, summary in summaries.items():
            summary_rows.append({
                'case_id': case['case_id'],
                'distance_id': distance_id,
                **summary,
                'effective_dimension': refined_result['distance_payloads'][distance_id]['effective_dimension'],
            })

    label_counts = dict(Counter(row['relational_ordering_label'] for row in summary_rows))

    work_plot_dir = Path('/tmp') / 'haos_iip_stage15_relational_geometry'
    work_plot_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = work_plot_dir / runsheet['plot_layout']['matrix_plot']
    metric_path = work_plot_dir / runsheet['plot_layout']['metric_summary_plot']
    ordering_path = work_plot_dir / runsheet['plot_layout']['ordering_stability_plot']
    dim_path = work_plot_dir / runsheet['plot_layout']['dimension_plot']
    cone_path = work_plot_dir / runsheet['plot_layout']['cone_plot']
    loop_path = work_plot_dir / runsheet['plot_layout']['loop_plot']

    create_distance_matrix_grid(matrix_path, case_results, runsheet)
    create_metric_summary(metric_path, csv_rows)
    create_ordering_flip_heatmap(ordering_path, summary_rows, runsheet)
    create_dimension_plot(dim_path, case_results)
    create_transport_shells(cone_path, case_results)
    create_loop_plot(loop_path, summary_rows)

    result = {
        'stage': 'Stage 15',
        'description': runsheet['description'],
        'base_seed_reference': base,
        'resolution_policy': runsheet['resolution_policy'],
        'cases': case_results,
        'summary_rows': summary_rows,
        'relational_ordering_counts': label_counts,
        'conclusion': 'Stage 15 tests whether frozen collective interaction data induce stable metric-like relational ordering under refinement.',
    }

    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug='stage15_relational_geometry',
        result=result,
        csv_rows=csv_rows,
        csv_fieldnames=CSV_FIELDS,
        plot_paths=[matrix_path, metric_path, ordering_path, dim_path, cone_path, loop_path],
    )
    note_path = ATLAS_NOTES / runsheet['note_name']
    write_note(note_path, str(json_path.relative_to(REPO_ROOT)), str(csv_path.relative_to(REPO_ROOT)), summary_rows, label_counts)
    append_log(
        title=f"Stage 15 Emergent Relational Geometry ({json_path.stem})",
        config_summary=f"cases={[case['case_id'] for case in cases]}, resolutions={[base_n, refined_n]}, distances={[item['distance_id'] for item in runsheet['distance_candidates']]}",
        result_path=json_path,
        stamped_plots=stamped_plots,
        observation=f"relational_ordering_counts={label_counts}",
        conclusion='the Stage 15 pilot tests whether frozen interaction persistence induces metric-like relational ordering under refinement',
    )

    print(f'wrote {json_path}')
    print(f'wrote {csv_path}')
    for item in stamped_plots:
        print(item)
    print(f'relational_ordering_counts={label_counts}')


if __name__ == '__main__':
    main()
