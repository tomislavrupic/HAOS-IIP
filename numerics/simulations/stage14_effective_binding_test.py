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
    spectral_moments,
    state_energy,
    suggested_dt,
    weighted_center,
)
from stage11_collective_wave_interaction import build_component_packet, classify_collective, make_leakage_fn
from stage12_coarse_field_structure import classify_coarse, connected_components_periodic, component_mask, edge_field_to_grid, gaussian_smooth_periodic, jaccard

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage14_binding_runs.json'
SIGMAS = (2.0, 4.0)
ALPHA = 0.6
ANALYSIS_STRIDE = 4
PROMOTION_MIN_LIFETIME_DELTA = 0.05
PROMOTION_MIN_BINDING_DELTA = 0.05
PROMOTION_MIN_COARSE_DELTA = 0.02
SETUP_CACHE: dict[tuple[int, str, float, float], tuple[Any, Any]] = {}

CSV_FIELDS = [
    'run_id',
    'phase',
    'representative_id',
    'condition_id',
    'condition_label',
    'feedback_eps',
    'n_side',
    'interaction_label',
    'coarse_label',
    'composite_lifetime',
    'overlap_persistence_time',
    'binding_persistence_score',
    'initial_mean_pair_distance',
    'min_mean_pair_distance',
    'final_mean_pair_distance',
    'max_mean_overlap',
    'phase_lock_indicator',
    'mean_dominant_area_fraction_sigma4',
    'coarse_persistence_sigma4',
    'mean_basin_count_sigma4',
    'constraint_max',
    'sector_leakage',
    'compare_to_null_lifetime_delta',
    'compare_to_null_lifetime_ratio',
    'compare_to_null_binding_delta',
    'compare_to_null_coarse_persistence_delta',
    'promoted_followup',
    'resolution_binding_ratio',
    'basin_invariance_flag',
    'composite_robustness_label',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage 14A Branch A effective binding pilot.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--representative-ids', nargs='*', default=[])
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


def edge_feedback_profile(midpoints: np.ndarray, total_q: np.ndarray, n_side: int, sigma: float) -> tuple[np.ndarray, np.ndarray]:
    grid = edge_field_to_grid(midpoints, np.abs(total_q), n_side)
    env = np.maximum(gaussian_smooth_periodic(grid, sigma), 0.0)
    max_env = float(np.max(env))
    if max_env > 0.0:
        env = env / max_env
    coords = np.floor(midpoints * n_side).astype(int) % n_side
    profile = env[coords[:, 0], coords[:, 1], coords[:, 2]]
    return profile.astype(float), env


def pairwise_distance(centers: list[np.ndarray], boundary_type: str) -> float:
    if len(centers) < 2:
        return 0.0
    values = []
    for i in range(len(centers)):
        for j in range(i + 1, len(centers)):
            delta = displacement(np.asarray([centers[j]], dtype=float), np.asarray(centers[i], dtype=float), boundary_type)[0]
            values.append(float(np.linalg.norm(delta)))
    return float(np.mean(values)) if values else 0.0


def pairwise_overlap(component_qs: list[np.ndarray]) -> float:
    if len(component_qs) < 2:
        return 0.0
    weights = [np.abs(q) ** 2 for q in component_qs]
    values = []
    for i in range(len(weights)):
        for j in range(i + 1, len(weights)):
            wi = weights[i]
            wj = weights[j]
            denom = math.sqrt(max(float(np.sum(wi)) * float(np.sum(wj)), 1.0e-12))
            values.append(float(np.sum(np.sqrt(wi * wj)) / denom))
    return float(np.mean(values)) if values else 0.0


def pairwise_phase(component_qs: list[np.ndarray]) -> float:
    if len(component_qs) < 2:
        return 0.0
    values = []
    for i in range(len(component_qs)):
        for j in range(i + 1, len(component_qs)):
            qi = component_qs[i]
            qj = component_qs[j]
            denom = max(float(np.linalg.norm(qi)) * float(np.linalg.norm(qj)), 1.0e-12)
            values.append(abs(float(np.vdot(qi, qj).real / denom)))
    return float(np.mean(values)) if values else 0.0


def render_grid_slice(ax, grid: np.ndarray, title: str, cmap: str = 'viridis') -> None:
    slice_idx = grid.shape[2] // 2
    img = grid[:, :, slice_idx].T
    im = ax.imshow(img, origin='lower', cmap=cmap)
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)


def simulate_run(rep: dict[str, Any], condition: dict[str, Any], resolution: int, defaults: dict[str, Any], base: dict[str, Any], feedback_sigma: float, overlap_threshold: float, work_plot_dir: Path) -> tuple[dict[str, Any], list[Path]]:
    data, projector = get_cached_setup(n_side=resolution, defaults=defaults, boundary_type=str(base['boundary_type']))
    kick_axis = int(defaults['kick_axis'])
    feedback_eps = float(condition['feedback_eps'])
    bandwidth = float(rep['bandwidth'])

    packet_qs: list[np.ndarray] = []
    packet_vs: list[np.ndarray] = []
    for center, amp, phase, kick_sign in zip(rep['packet_centers'], rep['packet_amplitudes'], rep['phase_offsets_rad'], rep['kick_signs']):
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
    pair_distance_series: list[float] = []
    overlap_series: list[float] = []
    phase_alignment_series: list[float] = []
    total_centers: list[list[float]] = []
    widths: list[float] = []
    norms: list[float] = []
    energies: list[float] = []
    anisotropies: list[float] = []
    sector_leakages: list[float] = []
    constraint_norms: list[float] = []
    coherences: list[float] = []
    coarse: dict[float, dict[str, Any]] = {
        sigma: {
            'basin_counts': [],
            'dominant_area_fractions': [],
            'envelope_variances': [],
            'stable_flags': [],
            'dominant_masks': [],
            'basin_lifetime_flags': [],
            'last_env': None,
        }
        for sigma in SIGMAS
    }

    def record(t: float, total_q: np.ndarray, total_v: np.ndarray) -> None:
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

        component_centers = []
        for comp_q in packet_qs:
            comp_weights = np.abs(comp_q) ** 2
            component_centers.append(weighted_center(data.midpoints, comp_weights, str(base['boundary_type'])))
        pair_distance_series.append(pairwise_distance(component_centers, str(base['boundary_type'])))
        overlap_series.append(pairwise_overlap(packet_qs))
        phase_alignment_series.append(pairwise_phase(packet_qs))

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
            coarse[sigma]['basin_lifetime_flags'].append(1.0 if dom_frac >= 0.05 else 0.0)
            coarse[sigma]['dominant_masks'].append(dominant_mask)

    total_q = np.sum(packet_qs, axis=0)
    total_v = np.sum(packet_vs, axis=0)
    record(0.0, total_q, total_v)
    for step_idx in range(steps):
        total_q = np.sum(packet_qs, axis=0)
        feedback_profile, _ = edge_feedback_profile(data.midpoints, total_q, resolution, feedback_sigma)
        next_qs: list[np.ndarray] = []
        next_vs: list[np.ndarray] = []
        for q, v in zip(packet_qs, packet_vs):
            feedback_term = feedback_eps * feedback_profile * q
            acc = np.asarray(projector(-(operator @ q) + feedback_term), dtype=float)
            v_half = v + 0.5 * dt * acc
            q_new = np.asarray(projector(q + dt * v_half), dtype=float)
            acc_new = np.asarray(projector(-(operator @ q_new) + feedback_eps * feedback_profile * q_new), dtype=float)
            v_new = np.asarray(projector(v_half + 0.5 * dt * acc_new), dtype=float)
            next_qs.append(q_new)
            next_vs.append(v_new)
        packet_qs = next_qs
        packet_vs = next_vs
        if ((step_idx + 1) % ANALYSIS_STRIDE == 0) or (step_idx == steps - 1):
            total_q = np.sum(packet_qs, axis=0)
            total_v = np.sum(packet_vs, axis=0)
            record((step_idx + 1) * dt, total_q, total_v)

    sample_dt = float(np.median(np.diff(times))) if len(times) >= 2 else dt
    centers_arr = np.asarray(total_centers, dtype=float)
    center_shift = float(np.linalg.norm(displacement(centers_arr[-1][None, :], centers_arr[0], str(base['boundary_type']))[0]))
    interaction_summary = {
        'packet_count': int(rep['packet_count']),
        'center_shift': center_shift,
        'width_ratio': float(widths[-1] / max(widths[0], 1.0e-12)),
        'max_mean_overlap': float(np.max(overlap_series)) if overlap_series else 0.0,
        'phase_lock_indicator': float(np.mean(phase_alignment_series)) if phase_alignment_series else 0.0,
        'composite_lifetime': float(np.sum(np.asarray(overlap_series) >= overlap_threshold) * sample_dt),
        'oscillation_count': 0,
        'initial_mean_pair_distance': float(pair_distance_series[0]) if pair_distance_series else 0.0,
        'min_mean_pair_distance': float(np.min(pair_distance_series)) if pair_distance_series else 0.0,
        'exchange_or_merger_flag': 'none',
        't_final': float(times[-1]),
    }
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

    binding_mask = (np.asarray(pair_distance_series) <= float(rep['binding_distance_threshold'])) & (np.asarray(overlap_series) >= overlap_threshold)
    binding_persistence_score = float(np.mean(binding_mask.astype(float))) if binding_mask.size else 0.0
    overlap_persistence_time = float(np.sum(np.asarray(overlap_series) >= overlap_threshold) * sample_dt)

    row = {
        'run_id': f"{rep['representative_id']}__{condition['condition_id']}__n{resolution}",
        'phase': 'base',
        'representative_id': rep['representative_id'],
        'condition_id': condition['condition_id'],
        'condition_label': condition['label'],
        'feedback_eps': feedback_eps,
        'n_side': resolution,
        'interaction_label': interaction_label,
        'coarse_label': coarse_label,
        'composite_lifetime': interaction_summary['composite_lifetime'],
        'overlap_persistence_time': overlap_persistence_time,
        'binding_persistence_score': binding_persistence_score,
        'initial_mean_pair_distance': interaction_summary['initial_mean_pair_distance'],
        'min_mean_pair_distance': interaction_summary['min_mean_pair_distance'],
        'final_mean_pair_distance': float(pair_distance_series[-1]) if pair_distance_series else 0.0,
        'max_mean_overlap': interaction_summary['max_mean_overlap'],
        'phase_lock_indicator': interaction_summary['phase_lock_indicator'],
        'mean_dominant_area_fraction_sigma4': coarse_summary['mean_dominant_area_fraction_sigma4'],
        'coarse_persistence_sigma4': coarse_summary['coarse_persistence_sigma4'],
        'mean_basin_count_sigma4': coarse_summary['mean_basin_count_sigma4'],
        'constraint_max': float(np.max(constraint_norms)) if constraint_norms else 0.0,
        'sector_leakage': float(np.max(sector_leakages)) if sector_leakages else 0.0,
        'compare_to_null_lifetime_delta': 0.0,
        'compare_to_null_lifetime_ratio': 1.0,
        'compare_to_null_binding_delta': 0.0,
        'compare_to_null_coarse_persistence_delta': 0.0,
        'promoted_followup': 0,
        'resolution_binding_ratio': 0.0,
        'basin_invariance_flag': 0,
        'composite_robustness_label': 'no composite gain',
        'notes': rep['selection_note'],
    }

    payload = {
        'run_id': row['run_id'],
        'representative': rep,
        'condition': condition,
        'resolution': resolution,
        'interaction_label': interaction_label,
        'coarse_label': coarse_label,
        'summary': dict(row),
        'times': times,
        'pair_distance_series': pair_distance_series,
        'overlap_series': overlap_series,
        'phase_alignment_series': phase_alignment_series,
        'coarse_sigma4': {
            'dominant_area_fractions': coarse[4.0]['dominant_area_fractions'],
            'basin_counts': coarse[4.0]['basin_counts'],
            'stable_flags': coarse[4.0]['stable_flags'],
        },
    }

    plot_paths: list[Path] = []
    trace_path = work_plot_dir / f"stage14_{row['run_id']}_pair_distance_overlap.png"
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))
    axes[0].plot(times, pair_distance_series, color='tab:blue')
    axes[0].axhline(float(rep['binding_distance_threshold']), color='black', linestyle='--', linewidth=1.0)
    axes[0].set_title('Mean pair distance')
    axes[0].set_xlabel('time')
    axes[0].grid(alpha=0.25)
    axes[1].plot(times, overlap_series, color='tab:orange')
    axes[1].axhline(overlap_threshold, color='black', linestyle='--', linewidth=1.0)
    axes[1].set_title('Mean overlap')
    axes[1].set_xlabel('time')
    axes[1].grid(alpha=0.25)
    fig.suptitle(f"Stage 14A: {row['run_id']}")
    fig.savefig(trace_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(trace_path)

    basin_path = work_plot_dir / f"stage14_{row['run_id']}_basin_trace.png"
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
        env_path = work_plot_dir / f"stage14_{row['run_id']}_envelope_snapshot.png"
        fig, ax = plt.subplots(figsize=(5.4, 4.6))
        render_grid_slice(ax, env, f"Envelope sigma=4: {row['run_id']}")
        fig.savefig(env_path, dpi=180, bbox_inches='tight')
        plt.close(fig)
        plot_paths.append(env_path)

    return payload, plot_paths


def enrich_vs_null(base_payloads: dict[tuple[str, str], dict[str, Any]], rows: list[dict[str, Any]]) -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
    promoted: dict[str, Any] | None = None
    best_score = None
    for row in rows:
        if row['condition_id'] == 'null_control' or row['phase'] != 'base':
            continue
        null_row = base_payloads[(row['representative_id'], 'null_control')]['summary']
        lifetime_delta = float(row['composite_lifetime']) - float(null_row['composite_lifetime'])
        lifetime_ratio = float(row['composite_lifetime']) / max(float(null_row['composite_lifetime']), 1.0e-12)
        binding_delta = float(row['binding_persistence_score']) - float(null_row['binding_persistence_score'])
        coarse_delta = float(row['coarse_persistence_sigma4']) - float(null_row['coarse_persistence_sigma4'])
        row['compare_to_null_lifetime_delta'] = lifetime_delta
        row['compare_to_null_lifetime_ratio'] = lifetime_ratio
        row['compare_to_null_binding_delta'] = binding_delta
        row['compare_to_null_coarse_persistence_delta'] = coarse_delta
        if lifetime_delta > 0.01 or binding_delta > 0.01:
            row['composite_robustness_label'] = 'weak composite stabilization'
        if (
            lifetime_delta >= PROMOTION_MIN_LIFETIME_DELTA
            and binding_delta >= PROMOTION_MIN_BINDING_DELTA
            and coarse_delta >= PROMOTION_MIN_COARSE_DELTA
            and interaction_rank(str(row['interaction_label'])) >= interaction_rank(str(null_row['interaction_label']))
        ):
            score = (lifetime_delta, binding_delta, coarse_delta)
            if best_score is None or score > best_score:
                best_score = score
                promoted = {
                    'representative_id': row['representative_id'],
                    'condition_id': row['condition_id'],
                }
    if promoted is not None:
        for row in rows:
            if row['representative_id'] == promoted['representative_id'] and row['condition_id'] == promoted['condition_id']:
                row['promoted_followup'] = 1
    return promoted, rows


def add_refinement_labels(base_rows: dict[tuple[str, str, int], dict[str, Any]], refined_rows: list[dict[str, Any]]) -> None:
    for row in refined_rows:
        base_row = base_rows[(row['representative_id'], row['condition_id'], int(row['n_side']) // 2)]
        row['resolution_binding_ratio'] = float(row['composite_lifetime']) / max(float(base_row['composite_lifetime']), 1.0e-12)
        row['basin_invariance_flag'] = 1 if str(row['coarse_label']) == str(base_row['coarse_label']) else 0
        if row['condition_id'] == 'null_control':
            row['composite_robustness_label'] = 'no composite gain'
            continue
        if (
            row['resolution_binding_ratio'] >= 1.1
            and row['basin_invariance_flag'] == 1
            and str(row['interaction_label']) == 'metastable composite regime'
            and str(base_row['interaction_label']) == 'metastable composite regime'
        ):
            row['composite_robustness_label'] = 'candidate scale-robust composite'
        elif row['resolution_binding_ratio'] >= 1.0 and str(row['interaction_label']) in {'metastable composite regime', 'transient binding regime'}:
            row['composite_robustness_label'] = 'partial scale stabilization'
        elif float(row['compare_to_null_lifetime_delta']) > 0.01 or float(row['compare_to_null_binding_delta']) > 0.01:
            row['composite_robustness_label'] = 'weak composite stabilization'
        else:
            row['composite_robustness_label'] = 'no composite gain'


def create_comparison_plot(path: Path, base_rows: list[dict[str, Any]]) -> None:
    reps = sorted({row['representative_id'] for row in base_rows})
    cond_order = ['null_control', 'very_weak_feedback', 'weak_feedback']
    cond_labels = {
        'null_control': 'null',
        'very_weak_feedback': 'very weak',
        'weak_feedback': 'weak',
    }
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    x = np.arange(len(reps))
    width = 0.22
    for idx, cond in enumerate(cond_order):
        vals_life = []
        vals_bind = []
        for rep in reps:
            match = next(row for row in base_rows if row['representative_id'] == rep and row['condition_id'] == cond)
            vals_life.append(float(match['composite_lifetime']))
            vals_bind.append(float(match['binding_persistence_score']))
        axes[0].bar(x + (idx - 1) * width, vals_life, width=width, label=cond_labels[cond])
        axes[1].bar(x + (idx - 1) * width, vals_bind, width=width, label=cond_labels[cond])
    for ax, title in zip(axes, ['Composite lifetime', 'Binding persistence score']):
        ax.set_xticks(x)
        ax.set_xticklabels(reps, rotation=20, ha='right')
        ax.set_title(title)
        ax.grid(alpha=0.25, axis='y')
        ax.legend(fontsize=8)
    fig.suptitle('Stage 14A null vs feedback comparison')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_regime_map(path: Path, rows: list[dict[str, Any]]) -> None:
    reps = sorted({row['representative_id'] for row in rows})
    conds = ['null_control', 'very_weak_feedback', 'weak_feedback']
    coarse_map = {
        'diffuse coarse field regime': 0,
        'multi-basin fluctuating regime': 1,
        'basin-dominated regime': 2,
    }
    data = np.zeros((len(reps), len(conds)), dtype=float)
    for i, rep in enumerate(reps):
        for j, cond in enumerate(conds):
            match = next(row for row in rows if row['representative_id'] == rep and row['condition_id'] == cond and row['phase'] == 'base')
            data[i, j] = coarse_map.get(str(match['coarse_label']), 0)
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    im = ax.imshow(data, cmap='viridis', vmin=0, vmax=2)
    ax.set_xticks(range(len(conds)))
    ax.set_xticklabels(['null', 'very weak', 'weak'])
    ax.set_yticks(range(len(reps)))
    ax.set_yticklabels(reps)
    ax.set_title('Stage 14A coarse morphology regime map')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def create_resolution_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.4))
    labels = []
    lifetimes = []
    bindings = []
    for row in rows:
        labels.append(f"{row['condition_id']}\n n={row['n_side']}")
        lifetimes.append(float(row['composite_lifetime']))
        bindings.append(float(row['binding_persistence_score']))
    axes[0].bar(range(len(labels)), lifetimes, color='tab:purple')
    axes[0].set_xticks(range(len(labels)))
    axes[0].set_xticklabels(labels)
    axes[0].set_title('Composite lifetime across resolution')
    axes[0].grid(alpha=0.25, axis='y')
    axes[1].bar(range(len(labels)), bindings, color='tab:olive')
    axes[1].set_xticks(range(len(labels)))
    axes[1].set_xticklabels(labels)
    axes[1].set_title('Binding persistence across resolution')
    axes[1].grid(alpha=0.25, axis='y')
    fig.suptitle('Stage 14A refinement follow-up')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(path: Path, json_rel: str, csv_rel: str, base_rows: list[dict[str, Any]], promoted: dict[str, Any] | None, refined_rows: list[dict[str, Any]], robustness_counts: dict[str, int]) -> None:
    lines = [
        '# Stage 14 Effective Binding Test v1',
        '',
        f'Timestamped summary: `{json_rel}`',
        f'Timestamped run table: `{csv_rel}`',
        '',
        f'Composite robustness counts: {robustness_counts}',
        '',
        'This pilot tests only Stage 14A Branch A: weak local envelope feedback against frozen Stage 13 null controls.',
        '',
        'Base comparisons:',
        '',
    ]
    for row in base_rows:
        lines.extend([
            f"- `{row['representative_id']}` / `{row['condition_id']}`",
            f"  - interaction: `{row['interaction_label']}`",
            f"  - coarse: `{row['coarse_label']}`",
            f"  - composite lifetime: `{row['composite_lifetime']:.4f}`",
            f"  - binding persistence: `{row['binding_persistence_score']:.4f}`",
            f"  - lifetime delta vs null: `{row['compare_to_null_lifetime_delta']:.4f}`",
            f"  - robustness: `{row['composite_robustness_label']}`",
        ])
    lines.extend(['', f'Promoted follow-up: `{promoted}`' if promoted else 'Promoted follow-up: `none`', ''])
    if refined_rows:
        lines.append('Refinement follow-up:')
        lines.append('')
        for row in refined_rows:
            lines.extend([
                f"- `{row['representative_id']}` / `{row['condition_id']}` / `n={row['n_side']}`",
                f"  - interaction: `{row['interaction_label']}`",
                f"  - coarse: `{row['coarse_label']}`",
                f"  - resolution binding ratio: `{row['resolution_binding_ratio']:.4f}`",
                f"  - basin invariance: `{row['basin_invariance_flag']}`",
                f"  - robustness: `{row['composite_robustness_label']}`",
            ])
        lines.append('')
    lines.extend([
        'Interpretation boundary:',
        '- results are read as modification of superposition morphology only',
        '- any persistence gain is an effective interaction persistence diagnostic, not a force claim',
        '',
    ])
    path.write_text('\n'.join(lines), encoding='utf-8')


def main() -> None:
    args = parse_args()
    defaults = load_stage10_defaults()
    runsheet = load_runsheet(args.runsheet)
    base = runsheet['base_seed_reference']
    representatives = selected_items(runsheet['representatives'], args.representative_ids, 'representative_id')
    conditions = selected_items(runsheet['conditions'], args.condition_ids, 'condition_id')
    work_plot_dir = Path('/tmp') / 'haos_iip_stage14_binding_plots'
    work_plot_dir.mkdir(parents=True, exist_ok=True)

    base_payloads: dict[tuple[str, str], dict[str, Any]] = {}
    base_rows: list[dict[str, Any]] = []
    payload_runs: list[dict[str, Any]] = []
    plot_paths: list[Path] = []

    for rep in representatives:
        resolution = int(rep['base_resolution'])
        for condition in conditions:
            payload, run_plots = simulate_run(
                rep=rep,
                condition=condition,
                resolution=resolution,
                defaults=defaults,
                base=base,
                feedback_sigma=float(runsheet['feedback_sigma']),
                overlap_threshold=float(runsheet['overlap_threshold']),
                work_plot_dir=work_plot_dir,
            )
            summary = payload['summary']
            base_rows.append(summary)
            base_payloads[(rep['representative_id'], condition['condition_id'])] = payload
            payload_runs.append({
                'phase': 'base',
                'representative_id': rep['representative_id'],
                'condition_id': condition['condition_id'],
                'payload': {
                    'interaction_label': payload['interaction_label'],
                    'coarse_label': payload['coarse_label'],
                    'summary': summary,
                },
            })
            plot_paths.extend(run_plots)

    promoted, base_rows = enrich_vs_null(base_payloads, base_rows)

    refined_rows: list[dict[str, Any]] = []
    if promoted is not None and not args.skip_followup:
        rep = next(item for item in representatives if item['representative_id'] == promoted['representative_id'])
        follow_conditions = [
            next(item for item in conditions if item['condition_id'] == 'null_control'),
            next(item for item in conditions if item['condition_id'] == promoted['condition_id']),
        ]
        for condition in follow_conditions:
            payload, run_plots = simulate_run(
                rep=rep,
                condition=condition,
                resolution=int(rep['refined_resolution']),
                defaults=defaults,
                base=base,
                feedback_sigma=float(runsheet['feedback_sigma']),
                overlap_threshold=float(runsheet['overlap_threshold']),
                work_plot_dir=work_plot_dir,
            )
            summary = payload['summary']
            summary['phase'] = 'refinement'
            null_base = base_payloads[(rep['representative_id'], 'null_control')]['summary']
            summary['compare_to_null_lifetime_delta'] = float(summary['composite_lifetime']) - float(null_base['composite_lifetime'])
            summary['compare_to_null_lifetime_ratio'] = float(summary['composite_lifetime']) / max(float(null_base['composite_lifetime']), 1.0e-12)
            summary['compare_to_null_binding_delta'] = float(summary['binding_persistence_score']) - float(null_base['binding_persistence_score'])
            summary['compare_to_null_coarse_persistence_delta'] = float(summary['coarse_persistence_sigma4']) - float(null_base['coarse_persistence_sigma4'])
            refined_rows.append(summary)
            payload_runs.append({
                'phase': 'refinement',
                'representative_id': rep['representative_id'],
                'condition_id': condition['condition_id'],
                'payload': {
                    'interaction_label': payload['interaction_label'],
                    'coarse_label': payload['coarse_label'],
                    'summary': summary,
                },
            })
            plot_paths.extend(run_plots)

        base_row_map = {(row['representative_id'], row['condition_id'], int(row['n_side'])): row for row in base_rows}
        add_refinement_labels(base_row_map, refined_rows)

    all_rows = base_rows + refined_rows
    robustness_counts = dict(Counter(row['composite_robustness_label'] for row in all_rows))

    comparison_plot = work_plot_dir / 'stage14A_null_vs_feedback_comparison.png'
    create_comparison_plot(comparison_plot, base_rows)
    plot_paths.append(comparison_plot)

    regime_map = work_plot_dir / 'stage14A_coarse_regime_map.png'
    create_regime_map(regime_map, base_rows)
    plot_paths.append(regime_map)

    if refined_rows:
        resolution_plot = work_plot_dir / 'stage14A_resolution_binding_comparison.png'
        create_resolution_plot(resolution_plot, refined_rows)
        plot_paths.append(resolution_plot)

    result = {
        'stage': 'Stage 14A',
        'branch': 'Branch A — Local envelope feedback',
        'description': runsheet['description'],
        'base_seed_reference': base,
        'feedback_sigma': runsheet['feedback_sigma'],
        'overlap_threshold': runsheet['overlap_threshold'],
        'base_runs': [
            {
                'representative_id': row['representative_id'],
                'condition_id': row['condition_id'],
                'summary': row,
            }
            for row in base_rows
        ],
        'refinement_runs': [
            {
                'representative_id': row['representative_id'],
                'condition_id': row['condition_id'],
                'summary': row,
            }
            for row in refined_rows
        ],
        'promoted_followup': promoted,
        'composite_robustness_counts': robustness_counts,
        'conclusion': 'Stage 14A compares null control against very weak and weak local envelope feedback while preserving Stage 13 comparability.',
    }

    stem = args.runsheet.stem.replace('_runs', '')
    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug=stem,
        result=result,
        csv_rows=all_rows,
        csv_fieldnames=CSV_FIELDS,
        plot_paths=plot_paths,
    )
    note_path = ATLAS_NOTES / runsheet['note_name']
    write_note(note_path, str(json_path.relative_to(REPO_ROOT)), str(csv_path.relative_to(REPO_ROOT)), base_rows, promoted, refined_rows, robustness_counts)
    append_log(
        title=f"Stage 14A Effective Binding Test ({json_path.stem})",
        config_summary=f"representatives={len(representatives)}, conditions={[c['condition_id'] for c in conditions]}, promoted={promoted}",
        result_path=json_path,
        stamped_plots=stamped_plots,
        observation=f"composite_robustness_counts={robustness_counts}",
        conclusion='the Stage 14A Branch A pilot tests whether weak local envelope feedback can improve composite persistence over frozen Stage 13 null controls',
    )

    print(f'wrote {json_path}')
    print(f'wrote {csv_path}')
    for item in stamped_plots:
        print(item)
    print(f'promoted_followup={promoted}')
    print(f'composite_robustness_counts={robustness_counts}')


if __name__ == '__main__':
    main()
