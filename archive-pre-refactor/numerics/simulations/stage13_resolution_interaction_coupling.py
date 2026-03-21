#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, deque
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
    create_field_snapshot,
    displacement,
    load_stage10_defaults,
    packet_width,
    plt,
    save_atlas_payload,
    spectral_moments,
    state_energy,
    suggested_dt,
    weighted_center,
)
from stage11_collective_wave_interaction import build_component_packet, classify_collective, make_leakage_fn
from stage12_coarse_field_structure import classify_coarse, edge_field_to_grid, gaussian_smooth_periodic, connected_components_periodic, component_mask, jaccard

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage13_resolution_runs.json'
CSV_FIELDS = [
    'case_id',
    'axis_set',
    'low_resolution',
    'high_resolution',
    'interaction_label_low',
    'interaction_label_high',
    'coarse_label_low',
    'coarse_label_high',
    'transition_type',
    'output_label',
    'recoverability_flag',
    'composite_lifetime_scaling',
    'overlap_scaling',
    'centroid_drift_scaling',
    'dominant_basin_area_scaling',
    'coarse_persistence_scaling',
    'constraint_norm_scaling',
    'notes',
]
SUMMARY_LABELS = [
    'scale-robust composite',
    'scale-fragile composite',
    'scale-induced diffusion',
    'resolution-dependent morphology',
]
SIGMAS = (2.0, 4.0)
ALPHA = 0.6
_SETUP_CACHE: dict[tuple[int, str, float, float], tuple[Any, Any]] = {}
ANALYSIS_STRIDE = 4


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage 13 resolution-interaction coupling pilot.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--case-ids', nargs='*', default=[])
    parser.add_argument('--max-cases', type=int, default=0)
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def selected_cases(all_cases: list[dict[str, Any]], case_ids: list[str], max_cases: int) -> list[dict[str, Any]]:
    cases = all_cases
    if case_ids:
        wanted = set(case_ids)
        cases = [case for case in cases if case['case_id'] in wanted]
    if max_cases > 0:
        cases = cases[:max_cases]
    return cases


def scale_centers(centers: list[list[float]], factor: float) -> list[list[float]]:
    arr = np.asarray(centers, dtype=float)
    centroid = np.mean(arr, axis=0)
    arr = centroid + factor * (arr - centroid)
    return arr.tolist()


def interaction_rank(label: str) -> int:
    order = {
        'dispersive independent regime': 0,
        'transient binding regime': 1,
        'metastable composite regime': 2,
        'oscillatory exchange regime': 3,
        'large-scale drift field regime': 3,
    }
    return order.get(label, 0)


def coarse_rank(label: str) -> int:
    order = {
        'diffuse coarse field regime': 0,
        'multi-basin fluctuating regime': 1,
        'basin-dominated regime': 2,
    }
    return order.get(label, 0)


def transform_case(case: dict[str, Any], axis_set: str, resolution: int) -> tuple[dict[str, Any], float]:
    run = json.loads(json.dumps(case))
    primary_sigma = 2.0
    if axis_set == 'A':
        if resolution > int(case['base_resolution']):
            primary_sigma = 4.0
    elif axis_set == 'B':
        if resolution > int(case['base_resolution']):
            run['bandwidth'] = float(case['bandwidth']) / 2.0
    elif axis_set == 'C':
        if resolution > int(case['base_resolution']):
            run['packet_centers'] = scale_centers(run['packet_centers'], 0.5)
    run['n_side'] = int(resolution)
    return run, primary_sigma


def get_cached_setup(n_side: int, defaults: dict[str, Any], boundary_type: str) -> tuple[Any, Any]:
    key = (int(n_side), str(boundary_type), float(defaults['epsilon']), float(defaults['harmonic_tol']))
    if key not in _SETUP_CACHE:
        _SETUP_CACHE[key] = build_transverse_setup(
            n_side=int(n_side),
            epsilon=float(defaults['epsilon']),
            harmonic_tol=float(defaults['harmonic_tol']),
            boundary_type=str(boundary_type),
        )
    return _SETUP_CACHE[key]


def simulate_resolution(run: dict[str, Any], defaults: dict[str, Any], base: dict[str, Any]) -> dict[str, Any]:
    n_side = int(run['n_side'])
    bandwidth = float(run.get('bandwidth', base['bandwidth']))
    data, projector = get_cached_setup(n_side=n_side, defaults=defaults, boundary_type=str(base['boundary_type']))
    kick_axis = int(defaults['kick_axis'])
    packet_qs: list[np.ndarray] = []
    packet_vs: list[np.ndarray] = []
    kick_signs = run.get('kick_signs', [1.0] * int(run['packet_count']))
    for center, amp, phase, kick_sign in zip(run['packet_centers'], run['packet_amplitudes'], run['phase_offsets_rad'], kick_signs):
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
    spectral_centroids: list[float] = []
    spectral_spreads: list[float] = []
    sector_leakages: list[float] = []
    constraint_norms: list[float] = []
    coherences: list[float] = []
    pair_distance_series: list[float] = []
    overlap_series: list[float] = []
    phase_alignment_series: list[float] = []
    field_grids: list[np.ndarray] = []
    coarse: dict[float, dict[str, Any]] = {
        sigma: {
            'envelopes': [],
            'basin_counts': [],
            'dominant_area_fractions': [],
            'envelope_variances': [],
            'stable_flags': [],
            'dominant_masks': [],
            'basin_lifetime_flags': [],
        }
        for sigma in SIGMAS
    }

    def pairwise_distance(centers: list[np.ndarray]) -> float:
        if len(centers) < 2:
            return 0.0
        vals = []
        for i in range(len(centers)):
            for j in range(i + 1, len(centers)):
                delta = displacement(np.asarray([centers[j]], dtype=float), np.asarray(centers[i], dtype=float), str(base['boundary_type']))[0]
                vals.append(float(np.linalg.norm(delta)))
        return float(np.mean(vals)) if vals else 0.0

    def pairwise_overlap(component_qs_local: list[np.ndarray]) -> float:
        if len(component_qs_local) < 2:
            return 0.0
        vals = []
        weights = [np.abs(q) ** 2 for q in component_qs_local]
        for i in range(len(weights)):
            for j in range(i + 1, len(weights)):
                wi = weights[i]
                wj = weights[j]
                denom = math.sqrt(max(float(np.sum(wi)) * float(np.sum(wj)), 1.0e-12))
                vals.append(float(np.sum(np.sqrt(wi * wj)) / denom))
        return float(np.mean(vals)) if vals else 0.0

    def pairwise_phase(component_qs_local: list[np.ndarray]) -> float:
        if len(component_qs_local) < 2:
            return 0.0
        vals = []
        for i in range(len(component_qs_local)):
            for j in range(i + 1, len(component_qs_local)):
                qi = component_qs_local[i]
                qj = component_qs_local[j]
                denom = max(float(np.linalg.norm(qi)) * float(np.linalg.norm(qj)), 1.0e-12)
                vals.append(abs(float(np.vdot(qi, qj).real / denom)))
        return float(np.mean(vals)) if vals else 0.0

    def record(t: float) -> None:
        total_q = np.sum(packet_qs, axis=0)
        total_v = np.sum(packet_vs, axis=0)
        total_weights = np.abs(total_q) ** 2
        total_center = weighted_center(data.midpoints, total_weights, str(base['boundary_type']))
        centroid, spread = spectral_moments(operator, total_q)
        times.append(float(t))
        total_centers.append(total_center.tolist())
        widths.append(packet_width(data.midpoints, total_weights, total_center, str(base['boundary_type'])))
        norms.append(float(np.linalg.norm(total_q)))
        energies.append(state_energy(operator, total_q, total_v))
        anisotropies.append(anisotropy_ratio(data.midpoints, total_weights, total_center, str(base['boundary_type'])))
        spectral_centroids.append(centroid)
        spectral_spreads.append(spread)
        sector_leakages.append(float(leakage_fn(total_q)))
        constraint_norms.append(float(np.linalg.norm(data.d0.T @ total_q)))
        coherences.append(coherence_score(total_weights))

        component_centers = []
        for comp_q in packet_qs:
            comp_weights = np.abs(comp_q) ** 2
            comp_center = weighted_center(data.midpoints, comp_weights, str(base['boundary_type']))
            component_centers.append(comp_center)
        pair_distance_series.append(pairwise_distance(component_centers))
        overlap_series.append(pairwise_overlap(packet_qs))
        phase_alignment_series.append(pairwise_phase(packet_qs))

        point_grid = edge_field_to_grid(data.midpoints, np.abs(total_q), n_side)
        field_grids.append(point_grid)
        for sigma in SIGMAS:
            env = gaussian_smooth_periodic(point_grid, sigma)
            env = np.maximum(env, 0.0)
            coarse[sigma]['envelopes'].append(env)
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
            coarse[sigma]['dominant_area_fractions'].append(dom_frac)
            prev_mask = coarse[sigma]['dominant_masks'][-1] if coarse[sigma]['dominant_masks'] else None
            if dominant_mask is None:
                stable = 0.0
            elif prev_mask is None:
                stable = 1.0
            else:
                stable = 1.0 if jaccard(prev_mask, dominant_mask) >= 0.5 else 0.0
            coarse[sigma]['stable_flags'].append(stable)
            coarse[sigma]['basin_lifetime_flags'].append(1.0 if dom_frac >= 0.05 else 0.0)
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

    centers_arr = np.asarray(total_centers, dtype=float)
    shift = float(np.linalg.norm(displacement(centers_arr[-1][None, :], centers_arr[0], str(base['boundary_type']))[0]))
    sample_dt = float(np.median(np.diff(times))) if len(times) >= 2 else dt
    oscillation_count = 0
    if len(pair_distance_series) >= 3:
        deriv = np.diff(np.asarray(pair_distance_series, dtype=float))
        deriv[np.abs(deriv) < 1.0e-4] = 0.0
        signs = [s for s in np.sign(deriv) if s != 0.0]
        oscillation_count = int(sum(1 for a, b in zip(signs[:-1], signs[1:]) if a != b)) if len(signs) >= 2 else 0
    interaction_summary = {
        'packet_count': int(run['packet_count']),
        'center_shift': shift,
        'width_ratio': float(widths[-1] / max(widths[0], 1.0e-12)),
        'max_mean_overlap': float(np.max(overlap_series)) if overlap_series else 0.0,
        'phase_lock_indicator': float(np.mean(phase_alignment_series)) if phase_alignment_series else 0.0,
        'composite_lifetime': float(np.sum(np.asarray(overlap_series) >= 0.35) * sample_dt),
        'oscillation_count': oscillation_count,
        'initial_mean_pair_distance': float(pair_distance_series[0]) if pair_distance_series else 0.0,
        'min_mean_pair_distance': float(np.min(pair_distance_series)) if pair_distance_series else 0.0,
        'exchange_or_merger_flag': 'none',
        't_final': float(times[-1]),
    }
    if interaction_summary['max_mean_overlap'] > 0.65 and interaction_summary['composite_lifetime'] > 0.18:
        interaction_summary['exchange_or_merger_flag'] = 'merger'
    elif oscillation_count >= 2 and interaction_summary['max_mean_overlap'] > 0.18:
        interaction_summary['exchange_or_merger_flag'] = 'exchange'
    interaction_label = classify_collective(interaction_summary)

    coarse_summary: dict[str, float] = {'t_final': float(times[-1])}
    for sigma in SIGMAS:
        stable_flags = np.asarray(coarse[sigma]['stable_flags'], dtype=float)
        lifetime_flags = np.asarray(coarse[sigma]['basin_lifetime_flags'], dtype=float)
        coarse_summary[f'mean_basin_count_sigma{int(sigma)}'] = float(np.mean(coarse[sigma]['basin_counts']))
        coarse_summary[f'mean_dominant_area_fraction_sigma{int(sigma)}'] = float(np.mean(coarse[sigma]['dominant_area_fractions']))
        coarse_summary[f'coarse_persistence_sigma{int(sigma)}'] = float(np.mean(stable_flags)) if stable_flags.size else 0.0
        coarse_summary[f'basin_lifetime_sigma{int(sigma)}'] = float(np.sum(lifetime_flags) * sample_dt)
        coarse_summary[f'mean_envelope_variance_sigma{int(sigma)}'] = float(np.mean(coarse[sigma]['envelope_variances']))
    coarse_label = classify_coarse(coarse_summary)

    summary = {
        'interaction_label': interaction_label,
        'coarse_label': coarse_label,
        'center_shift': shift,
        'width_ratio': float(widths[-1] / max(widths[0], 1.0e-12)),
        'coherence_ratio': float(coherences[-1] / max(coherences[0], 1.0e-12)),
        'max_anisotropy': float(np.max(anisotropies)),
        'max_constraint_norm': float(np.max(constraint_norms)),
        'max_sector_leakage': float(np.max(sector_leakages)),
        'max_mean_overlap': interaction_summary['max_mean_overlap'],
        'composite_lifetime': interaction_summary['composite_lifetime'],
        'phase_lock_indicator': interaction_summary['phase_lock_indicator'],
        'initial_mean_pair_distance': interaction_summary['initial_mean_pair_distance'],
        'min_mean_pair_distance': interaction_summary['min_mean_pair_distance'],
        'mean_basin_count_sigma2': coarse_summary['mean_basin_count_sigma2'],
        'mean_basin_count_sigma4': coarse_summary['mean_basin_count_sigma4'],
        'mean_dominant_area_fraction_sigma2': coarse_summary['mean_dominant_area_fraction_sigma2'],
        'mean_dominant_area_fraction_sigma4': coarse_summary['mean_dominant_area_fraction_sigma4'],
        'coarse_persistence_sigma2': coarse_summary['coarse_persistence_sigma2'],
        'coarse_persistence_sigma4': coarse_summary['coarse_persistence_sigma4'],
        'mean_envelope_variance_sigma2': coarse_summary['mean_envelope_variance_sigma2'],
        'mean_envelope_variance_sigma4': coarse_summary['mean_envelope_variance_sigma4'],
        't_final': float(times[-1]),
    }

    return {
        'run': run,
        'times': times,
        'pair_distance_series': pair_distance_series,
        'overlap_series': overlap_series,
        'phase_alignment_series': phase_alignment_series,
        'field_grids': field_grids,
        'coarse': coarse,
        'summary': summary,
    }


def infer_transition_type(low: dict[str, Any], high: dict[str, Any]) -> str:
    int_low = low['summary']['interaction_label']
    int_high = high['summary']['interaction_label']
    coarse_low = low['summary']['coarse_label']
    coarse_high = high['summary']['coarse_label']
    int_same = int_low == int_high
    coarse_same = coarse_low == coarse_high
    if int_same and coarse_same:
        return 'stable under refinement'
    if int_same and not coarse_same:
        return 'coarse-only drift'
    if not int_same and coarse_same:
        return 'interaction-only drift'
    if interaction_rank(int_high) < interaction_rank(int_low) or coarse_rank(coarse_high) < coarse_rank(coarse_low):
        return 'softening under refinement'
    return 'class flip under refinement'


def infer_output_label(low: dict[str, Any], high: dict[str, Any]) -> str:
    int_low = low['summary']['interaction_label']
    int_high = high['summary']['interaction_label']
    coarse_low = low['summary']['coarse_label']
    coarse_high = high['summary']['coarse_label']
    if int_low == 'metastable composite regime' and int_high == 'metastable composite regime' and coarse_low == 'basin-dominated regime' and coarse_high == 'basin-dominated regime':
        return 'scale-robust composite'
    if int_low == 'metastable composite regime' and (int_high != int_low or coarse_high != coarse_low):
        return 'scale-fragile composite'
    if coarse_high == 'diffuse coarse field regime' and coarse_low != 'diffuse coarse field regime':
        return 'scale-induced diffusion'
    return 'resolution-dependent morphology'


def infer_recoverability(low: dict[str, Any], high: dict[str, Any]) -> str:
    score_low = interaction_rank(low['summary']['interaction_label']) + coarse_rank(low['summary']['coarse_label'])
    score_high = interaction_rank(high['summary']['interaction_label']) + coarse_rank(high['summary']['coarse_label'])
    if score_high > score_low:
        return 'more stable with refinement'
    if score_high < score_low:
        return 'less stable with refinement'
    return 'unchanged with refinement'


def make_comparison_plots(case_id: str, axis_set: str, low: dict[str, Any], high: dict[str, Any], low_sigma: float, high_sigma: float, work_plot_dir: Path) -> list[Path]:
    plot_paths: list[Path] = []
    stem = f'stage13_{case_id}_{axis_set}'

    field_path = work_plot_dir / f'{stem}_field_snapshot.png'
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6))
    for ax, payload, title in [(axes[0], low, f'n={low["run"]["n_side"]}'), (axes[1], high, f'n={high["run"]["n_side"]}')]:
        grid = payload['field_grids'][-1]
        img = grid[:, :, grid.shape[2] // 2].T
        im = ax.imshow(img, origin='lower', cmap='viridis')
        ax.set_title(title)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.suptitle(f'Field snapshot comparison: {case_id} / {axis_set}')
    fig.savefig(field_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(field_path)

    trace_path = work_plot_dir / f'{stem}_overlap_pair_distance.png'
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4))
    axes[0].plot(low['times'], low['overlap_series'], label=f'n={low["run"]["n_side"]}')
    axes[0].plot(high['times'], high['overlap_series'], label=f'n={high["run"]["n_side"]}')
    axes[0].set_title('Mean overlap')
    axes[0].set_xlabel('time')
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=8)
    axes[1].plot(low['times'], low['pair_distance_series'], label=f'n={low["run"]["n_side"]}')
    axes[1].plot(high['times'], high['pair_distance_series'], label=f'n={high["run"]["n_side"]}')
    axes[1].set_title('Mean pair distance')
    axes[1].set_xlabel('time')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.suptitle(f'Interaction traces: {case_id} / {axis_set}')
    fig.savefig(trace_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(trace_path)

    env_path = work_plot_dir / f'{stem}_envelope_snapshot.png'
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6))
    for ax, payload, sigma, title in [(axes[0], low, low_sigma, f'n={low["run"]["n_side"]}, sigma={int(low_sigma)}'), (axes[1], high, high_sigma, f'n={high["run"]["n_side"]}, sigma={int(high_sigma)}')]:
        env = payload['coarse'][sigma]['envelopes'][-1]
        img = env[:, :, env.shape[2] // 2].T
        im = ax.imshow(img, origin='lower', cmap='viridis')
        ax.set_title(title)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.suptitle(f'Envelope comparison: {case_id} / {axis_set}')
    fig.savefig(env_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(env_path)

    basin_path = work_plot_dir / f'{stem}_basin_segmentation.png'
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6))
    for ax, payload, sigma, title in [(axes[0], low, low_sigma, f'n={low["run"]["n_side"]}, sigma={int(low_sigma)}'), (axes[1], high, high_sigma, f'n={high["run"]["n_side"]}, sigma={int(high_sigma)}')]:
        env = payload['coarse'][sigma]['envelopes'][-1]
        max_env = float(np.max(env))
        mask = env > (ALPHA * max_env) if max_env > 0.0 else np.zeros_like(env, dtype=bool)
        img = mask[:, :, mask.shape[2] // 2].T.astype(float)
        im = ax.imshow(img, origin='lower', cmap='magma')
        ax.set_title(title)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.suptitle(f'Basin segmentation: {case_id} / {axis_set}')
    fig.savefig(basin_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(basin_path)

    return plot_paths


def create_transition_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    inter_pairs = Counter((row['interaction_label_low'], row['interaction_label_high']) for row in rows)
    coarse_pairs = Counter((row['coarse_label_low'], row['coarse_label_high']) for row in rows)
    axes[0].bar(range(len(inter_pairs)), list(inter_pairs.values()))
    axes[0].set_xticks(range(len(inter_pairs)))
    axes[0].set_xticklabels([f'{a}\n→\n{b}' for a, b in inter_pairs.keys()], rotation=45, ha='right', fontsize=7)
    axes[0].set_title('Interaction label transitions')
    axes[0].grid(alpha=0.25, axis='y')
    axes[1].bar(range(len(coarse_pairs)), list(coarse_pairs.values()), color='tab:orange')
    axes[1].set_xticks(range(len(coarse_pairs)))
    axes[1].set_xticklabels([f'{a}\n→\n{b}' for a, b in coarse_pairs.keys()], rotation=45, ha='right', fontsize=7)
    axes[1].set_title('Coarse label transitions')
    axes[1].grid(alpha=0.25, axis='y')
    fig.suptitle('Stage 13 resolution-transition matrix')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_persistence_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    counts = Counter(row['transition_type'] for row in rows)
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    ax.bar(range(len(counts)), list(counts.values()), color='tab:green')
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(list(counts.keys()), rotation=30, ha='right')
    ax.set_title('Stage 13 transition stability counts')
    ax.grid(alpha=0.25, axis='y')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_metric_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    colors = {
        'scale-robust composite': 'tab:green',
        'scale-fragile composite': 'tab:red',
        'scale-induced diffusion': 'tab:blue',
        'resolution-dependent morphology': 'tab:orange',
    }
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for row in rows:
        ax.scatter(
            float(row['composite_lifetime_scaling']),
            float(row['dominant_basin_area_scaling']),
            color=colors.get(str(row['output_label']), 'black'),
            s=80,
            alpha=0.85,
            label=str(row['output_label']),
        )
        ax.text(float(row['composite_lifetime_scaling']) + 0.01, float(row['dominant_basin_area_scaling']), row['case_id'], fontsize=7)
    handles, labels = ax.get_legend_handles_labels()
    unique = {}
    for h, l in zip(handles, labels):
        unique.setdefault(l, h)
    ax.legend(unique.values(), unique.keys(), fontsize=8)
    ax.set_xlabel('composite lifetime scaling')
    ax.set_ylabel('dominant basin area scaling')
    ax.set_title('Stage 13 metric drift summary')
    ax.grid(alpha=0.25)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(path: Path, json_rel: str, csv_rel: str, rows: list[dict[str, Any]], output_counts: dict[str, int]) -> None:
    lines = [
        '# Stage 13 Resolution-Interaction Coupling Atlas v1',
        '',
        f'Timestamped summary: `{json_rel}`',
        '',
        f'Timestamped run table: `{csv_rel}`',
        '',
        f'Output-label counts: {output_counts}',
        '',
        'This pilot remains classification-first and interpretation-minimal.',
        'Stage 13 measures how frozen Stage 11 interaction labels and Stage 12 coarse labels behave under coupled scale variation.',
        '',
        'Per-case summary:',
        '',
    ]
    for row in rows:
        lines.extend([
            f"- `{row['case_id']}` / axis `{row['axis_set']}`",
            f"  - interaction: `{row['interaction_label_low']}` -> `{row['interaction_label_high']}`",
            f"  - coarse: `{row['coarse_label_low']}` -> `{row['coarse_label_high']}`",
            f"  - transition: `{row['transition_type']}`",
            f"  - output label: `{row['output_label']}`",
            f"  - recoverability: `{row['recoverability_flag']}`",
        ])
    lines.extend([
        '',
        'Output labels remain descriptive only:',
        '- scale-robust composite',
        '- scale-fragile composite',
        '- scale-induced diffusion',
        '- resolution-dependent morphology',
        '',
        'Do not read these runs as continuum, force, or nonlinear interaction claims. They are coupled scale-morphology diagnostics only.',
        '',
    ])
    path.write_text('\n'.join(lines), encoding='utf-8')


def run_signature(run: dict[str, Any]) -> str:
    return json.dumps(run, sort_keys=True)


def main() -> None:
    args = parse_args()
    defaults = load_stage10_defaults()
    runsheet = load_runsheet(args.runsheet)
    cases = selected_cases(runsheet['cases'], args.case_ids, args.max_cases)
    base = runsheet['base_seed_reference']
    stem = args.runsheet.stem.replace('_runs', '')
    note_path = ATLAS_NOTES / runsheet.get('note_name', 'Stage_13_Resolution_Interaction_Coupling_v1.md')
    work_plot_dir = Path('/tmp') / f'haos_iip_{stem}_plots'
    work_plot_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    payload_cases: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    simulation_cache: dict[str, dict[str, Any]] = {}
    for case in cases:
        for axis_set in runsheet['axis_sets']:
            print(f'[stage13] case={case["case_id"]} axis={axis_set} starting', flush=True)
            low_run, low_sigma = transform_case(case, axis_set, int(case['base_resolution']))
            high_run, high_sigma = transform_case(case, axis_set, int(case['refined_resolution']))
            low_key = run_signature(low_run)
            high_key = run_signature(high_run)
            if low_key not in simulation_cache:
                print(f'[stage13] building low resolution n={low_run["n_side"]} for {case["case_id"]}', flush=True)
                simulation_cache[low_key] = simulate_resolution(low_run, defaults, base)
            if high_key not in simulation_cache:
                print(f'[stage13] building high resolution n={high_run["n_side"]} for {case["case_id"]} axis={axis_set}', flush=True)
                simulation_cache[high_key] = simulate_resolution(high_run, defaults, base)
            low = simulation_cache[low_key]
            high = simulation_cache[high_key]
            transition_type = infer_transition_type(low, high)
            output_label = infer_output_label(low, high)
            recoverability = infer_recoverability(low, high)
            def scale(num: float, den: float) -> float:
                return float(num / max(abs(den), 1.0e-12))
            row = {
                'case_id': case['case_id'],
                'axis_set': axis_set,
                'low_resolution': int(case['base_resolution']),
                'high_resolution': int(case['refined_resolution']),
                'interaction_label_low': low['summary']['interaction_label'],
                'interaction_label_high': high['summary']['interaction_label'],
                'coarse_label_low': low['summary']['coarse_label'],
                'coarse_label_high': high['summary']['coarse_label'],
                'transition_type': transition_type,
                'output_label': output_label,
                'recoverability_flag': recoverability,
                'composite_lifetime_scaling': scale(high['summary']['composite_lifetime'], low['summary']['composite_lifetime']),
                'overlap_scaling': scale(high['summary']['max_mean_overlap'], low['summary']['max_mean_overlap']),
                'centroid_drift_scaling': scale(high['summary']['center_shift'], low['summary']['center_shift']),
                'dominant_basin_area_scaling': scale(high['summary']['mean_dominant_area_fraction_sigma4'], low['summary']['mean_dominant_area_fraction_sigma4']),
                'coarse_persistence_scaling': scale(high['summary']['coarse_persistence_sigma4'], low['summary']['coarse_persistence_sigma4']),
                'constraint_norm_scaling': scale(high['summary']['max_constraint_norm'], low['summary']['max_constraint_norm']),
                'notes': case['selection_note'],
            }
            rows.append(row)
            payload_cases.append({
                'case': case,
                'axis_set': axis_set,
                'low_run': low['run'],
                'high_run': high['run'],
                'low_summary': low['summary'],
                'high_summary': high['summary'],
                'transition_type': transition_type,
                'output_label': output_label,
                'recoverability_flag': recoverability,
            })
            plot_paths.extend(make_comparison_plots(case['case_id'], axis_set, low, high, low_sigma, high_sigma, work_plot_dir))
            print(f'[stage13] case={case["case_id"]} axis={axis_set} done -> {output_label}', flush=True)

    output_counts = dict(Counter(row['output_label'] for row in rows))
    transition_plot = work_plot_dir / f'{stem}_resolution_transition_matrix.png'
    persistence_plot = work_plot_dir / f'{stem}_interaction_coarse_persistence.png'
    metric_plot = work_plot_dir / f'{stem}_metric_drift_summary.png'
    create_transition_plot(transition_plot, rows)
    create_persistence_plot(persistence_plot, rows)
    create_metric_plot(metric_plot, rows)
    plot_paths.extend([transition_plot, persistence_plot, metric_plot])

    result = {
        'stage': 'Stage 13',
        'description': runsheet['description'],
        'base_seed_reference': base,
        'axis_sets': runsheet['axis_sets'],
        'summary_labels': SUMMARY_LABELS,
        'output_label_counts': output_counts,
        'cases': payload_cases,
        'conclusion': 'Stage 13 measures how frozen interaction and coarse-envelope labels behave under coupled scale variation without making continuum or nonlinear claims.',
    }

    json_path, csv_path, stamped_plots, timestamp = save_atlas_payload(
        experiment_slug=stem,
        result=result,
        csv_rows=rows,
        csv_fieldnames=CSV_FIELDS,
        plot_paths=plot_paths,
    )
    write_note(note_path, str(json_path.relative_to(REPO_ROOT)), str(csv_path.relative_to(REPO_ROOT)), rows, output_counts)
    append_log(
        title=f'Stage 13 Resolution-Interaction Coupling ({json_path.stem})',
        config_summary=f'cases={len(cases)}, axis_sets={runsheet["axis_sets"]}, resolutions={[cases[0]["base_resolution"], cases[0]["refined_resolution"]] if cases else []}',
        result_path=json_path,
        stamped_plots=stamped_plots,
        observation=f'output_label_counts={output_counts}',
        conclusion='the Stage 13 pilot couples resolution changes to interaction and coarse-envelope morphology without altering the frozen Stage 11 or Stage 12 classifiers',
    )

    print(f'wrote {json_path}')
    print(f'wrote {csv_path}')
    for item in stamped_plots:
        print(item)
    print(f'output_label_counts={output_counts}')


if __name__ == '__main__':
    main()
