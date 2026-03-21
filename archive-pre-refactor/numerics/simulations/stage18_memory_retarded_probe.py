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

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage18_memory_retarded_runs.json'
SETUP_CACHE: dict[tuple[int, str, float, float], tuple[Any, Any]] = {}
FEEDBACK_SIGMA = 4.0
ANALYSIS_STRIDE = 4
COMPOSITE_GAIN_MIN = 0.05
BINDING_GAIN_MIN = 0.03
COARSE_GAIN_MIN = 0.02
DELAY_ADVANTAGE_MIN = 0.05

CSV_FIELDS = [
    'run_id',
    'phase',
    'memory_class_id',
    'representative_id',
    'condition_id',
    'memory_strength',
    'tau_over_w',
    'delay_steps',
    'resolution',
    'interaction_label',
    'coarse_label',
    'composite_lifetime',
    'binding_persistence_score',
    'coarse_persistence_sigma4',
    'delayed_response_correlation',
    'delayed_response_advantage',
    'new_ordering_class',
    'composite_lifetime_delta',
    'binding_persistence_delta',
    'coarse_persistence_delta',
    'gate_met',
    'promoted_followup',
    'output_label',
    'constraint_max',
    'sector_leakage',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage 18 Branch-A-first memory retarded probe.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--memory-class-id', default='A_single_lag_retarded_term')
    parser.add_argument('--run-ids', nargs='*', default=[])
    parser.add_argument('--skip-followup', action='store_true')
    parser.add_argument('--append-log', action='store_true')
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def selected_runs(runs: list[dict[str, Any]], run_ids: list[str], memory_class_id: str) -> list[dict[str, Any]]:
    rows = [run for run in runs if run['memory_class_id'] == memory_class_id]
    if not run_ids:
        return rows
    wanted = set(run_ids)
    return [run for run in rows if run['run_id'] in wanted]


def lookup_by_id(items: list[dict[str, Any]], key: str) -> dict[str, dict[str, Any]]:
    return {item[key]: item for item in items}


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


def edge_feedback_profile(midpoints: np.ndarray, q: np.ndarray, n_side: int, sigma: float) -> np.ndarray:
    grid = edge_field_to_grid(midpoints, np.abs(q), n_side)
    env = np.maximum(gaussian_smooth_periodic(grid, sigma), 0.0)
    max_env = float(np.max(env))
    if max_env > 0.0:
        env = env / max_env
    coords = np.floor(midpoints * n_side).astype(int) % n_side
    return env[coords[:, 0], coords[:, 1], coords[:, 2]].astype(float)


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


def safe_corr(x: list[float], y: list[float]) -> float:
    if len(x) < 3 or len(y) < 3:
        return 0.0
    xa = np.asarray(x, dtype=float)
    ya = np.asarray(y, dtype=float)
    if float(np.std(xa)) < 1.0e-12 or float(np.std(ya)) < 1.0e-12:
        return 0.0
    return float(np.corrcoef(xa, ya)[0, 1])


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
    condition: dict[str, Any],
    runsheet: dict[str, Any],
    defaults: dict[str, Any],
    base: dict[str, Any],
    work_plot_dir: Path,
) -> tuple[dict[str, Any], list[Path]]:
    resolution = int(run['resolution'])
    boundary_type = str(base['boundary_type'])
    data, projector = get_cached_setup(resolution, defaults, boundary_type)
    kick_axis = int(defaults['kick_axis'])
    bandwidth = float(base['bandwidth'])
    memory_strength = float(run['memory_strength'])
    tau_over_w = float(runsheet['branch_a_first_variant']['locked_tau_over_w'])

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
    delay_time = tau_over_w * bandwidth
    delay_steps = max(1, int(round(delay_time / max(dt, 1.0e-12))))
    leakage_fn = make_leakage_fn(projector)

    profile_history: list[deque[np.ndarray]] = [deque(maxlen=delay_steps + 1) for _ in range(int(rep['packet_count']))]
    for idx, q in enumerate(packet_qs):
        profile_history[idx].append(edge_feedback_profile(data.midpoints, q, resolution, FEEDBACK_SIGMA))

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
    delayed_profile_norm_series: list[float] = []
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

    def current_delayed_profiles() -> list[np.ndarray]:
        out: list[np.ndarray] = []
        for history in profile_history:
            if memory_strength == 0.0:
                out.append(history[0])
            elif len(history) > delay_steps:
                out.append(history[0])
            else:
                out.append(history[0])
        return out

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
        for comp_q in packet_qs:
            comp_weights = np.abs(comp_q) ** 2
            comp_centers.append(weighted_center(data.midpoints, comp_weights, boundary_type))
        mean_pair, _min_pair = pairwise_distance_stats(comp_centers, boundary_type)
        pair_distance_series.append(mean_pair)
        overlap_series.append(pairwise_overlap_mean(packet_qs))
        delayed_profiles = current_delayed_profiles()
        delayed_profile_norm_series.append(float(np.mean([np.linalg.norm(profile) for profile in delayed_profiles])))

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
        delayed_profiles = current_delayed_profiles()
        next_qs: list[np.ndarray] = []
        next_vs: list[np.ndarray] = []
        for idx, (q, v) in enumerate(zip(packet_qs, packet_vs)):
            delayed_term = delayed_profiles[idx] * q
            acc = np.asarray(projector(-(operator @ q) + memory_strength * delayed_term), dtype=float)
            v_half = v + 0.5 * dt * acc
            q_new = np.asarray(projector(q + dt * v_half), dtype=float)
            acc_new = np.asarray(projector(-(operator @ q_new) + memory_strength * (delayed_profiles[idx] * q_new)), dtype=float)
            v_new = np.asarray(projector(v_half + 0.5 * dt * acc_new), dtype=float)
            next_qs.append(q_new)
            next_vs.append(v_new)
        packet_qs = next_qs
        packet_vs = next_vs
        for idx, q in enumerate(packet_qs):
            profile_history[idx].append(edge_feedback_profile(data.midpoints, q, resolution, FEEDBACK_SIGMA))
        if ((step_idx + 1) % ANALYSIS_STRIDE == 0) or (step_idx == steps - 1):
            record((step_idx + 1) * dt)

    sample_dt = float(np.median(np.diff(times))) if len(times) >= 2 else dt
    centers_arr = np.asarray(total_centers, dtype=float)
    center_shift = float(np.linalg.norm(displacement(centers_arr[-1][None, :], centers_arr[0], boundary_type)[0]))
    interaction_summary = {
        'packet_count': int(rep['packet_count']),
        'center_shift': center_shift,
        'width_ratio': float(widths[-1] / max(widths[0], 1.0e-12)),
        'max_mean_overlap': float(np.max(overlap_series)) if overlap_series else 0.0,
        'phase_lock_indicator': 0.0,
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
    delayed_corr = safe_corr(delayed_profile_norm_series, overlap_series)

    row = {
        'run_id': run['run_id'],
        'phase': 'base',
        'memory_class_id': run['memory_class_id'],
        'representative_id': run['representative_id'],
        'condition_id': run['condition_id'],
        'memory_strength': memory_strength,
        'tau_over_w': 0.0 if memory_strength == 0.0 else tau_over_w,
        'delay_steps': 0 if memory_strength == 0.0 else delay_steps,
        'resolution': resolution,
        'interaction_label': interaction_label,
        'coarse_label': coarse_label,
        'composite_lifetime': float(interaction_summary['composite_lifetime']),
        'binding_persistence_score': binding_persistence,
        'coarse_persistence_sigma4': coarse_summary['coarse_persistence_sigma4'],
        'delayed_response_correlation': delayed_corr,
        'delayed_response_advantage': 0.0,
        'new_ordering_class': 0,
        'composite_lifetime_delta': 0.0,
        'binding_persistence_delta': 0.0,
        'coarse_persistence_delta': 0.0,
        'gate_met': 0,
        'promoted_followup': 0,
        'output_label': 'no memory dynamic effect',
        'constraint_max': float(np.max(constraint_norms)) if constraint_norms else 0.0,
        'sector_leakage': float(np.max(sector_leakages)) if sector_leakages else 0.0,
        'notes': rep['notes'],
    }

    payload = {
        'run': run,
        'summary': row,
        'times': times,
        'pair_distance_series': pair_distance_series,
        'overlap_series': overlap_series,
        'delayed_profile_norm_series': delayed_profile_norm_series,
        'coarse_sigma4': {
            'dominant_area_fractions': coarse[4.0]['dominant_area_fractions'],
            'basin_counts': coarse[4.0]['basin_counts'],
            'stable_flags': coarse[4.0]['stable_flags'],
            'last_env': coarse[4.0]['last_env'],
        },
    }

    plot_paths: list[Path] = []
    trace_path = work_plot_dir / f"stage18_{run['run_id']}_pair_distance_overlap.png"
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))
    axes[0].plot(times, pair_distance_series, color='tab:blue')
    axes[0].set_title('Mean pair distance')
    axes[0].set_xlabel('time')
    axes[0].grid(alpha=0.25)
    axes[1].plot(times, overlap_series, color='tab:orange')
    axes[1].set_title('Mean overlap')
    axes[1].set_xlabel('time')
    axes[1].grid(alpha=0.25)
    fig.suptitle(f"Stage 18: {run['run_id']}")
    fig.savefig(trace_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(trace_path)

    basin_path = work_plot_dir / f"stage18_{run['run_id']}_basin_trace.png"
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))
    axes[0].plot(times, coarse[4.0]['dominant_area_fractions'], color='tab:green')
    axes[0].set_title('Dominant basin area (sigma=4)')
    axes[0].set_xlabel('time')
    axes[0].grid(alpha=0.25)
    axes[1].plot(times, delayed_profile_norm_series, color='tab:purple')
    axes[1].set_title('Delayed profile norm')
    axes[1].set_xlabel('time')
    axes[1].grid(alpha=0.25)
    fig.savefig(basin_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(basin_path)

    env = coarse[4.0]['last_env']
    if env is not None:
        env_path = work_plot_dir / f"stage18_{run['run_id']}_envelope_snapshot.png"
        fig, ax = plt.subplots(figsize=(5.4, 4.6))
        render_grid_slice(ax, env, f"Envelope sigma=4: {run['run_id']}")
        fig.savefig(env_path, dpi=180, bbox_inches='tight')
        plt.close(fig)
        plot_paths.append(env_path)

    return payload, plot_paths


def classify_base_rows(rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, Any] | None]:
    promoted: dict[str, Any] | None = None
    best_score = None
    base_map = {(row['representative_id'], row['condition_id']): row for row in rows}
    for row in rows:
        if row['condition_id'] == 'null_control':
            continue
        null_row = base_map[(row['representative_id'], 'null_control')]
        row['composite_lifetime_delta'] = float(row['composite_lifetime']) - float(null_row['composite_lifetime'])
        row['binding_persistence_delta'] = float(row['binding_persistence_score']) - float(null_row['binding_persistence_score'])
        row['coarse_persistence_delta'] = float(row['coarse_persistence_sigma4']) - float(null_row['coarse_persistence_sigma4'])
        row['delayed_response_advantage'] = float(row['delayed_response_correlation']) - float(null_row['delayed_response_correlation'])
        row['new_ordering_class'] = int(
            interaction_rank(str(row['interaction_label'])) > interaction_rank(str(null_row['interaction_label']))
            or coarse_rank(str(row['coarse_label'])) > coarse_rank(str(null_row['coarse_label']))
        )
        gate_hit = (
            row['composite_lifetime_delta'] >= COMPOSITE_GAIN_MIN
            or row['binding_persistence_delta'] >= BINDING_GAIN_MIN
            or row['coarse_persistence_delta'] >= COARSE_GAIN_MIN
            or row['new_ordering_class'] == 1
        )
        row['gate_met'] = int(gate_hit and row['delayed_response_advantage'] >= DELAY_ADVANTAGE_MIN)
        if row['gate_met']:
            row['output_label'] = 'transient history response'
            score = (
                float(row['composite_lifetime_delta'])
                + float(row['binding_persistence_delta'])
                + float(row['coarse_persistence_delta'])
                + float(row['delayed_response_advantage'])
            )
            if best_score is None or score > best_score:
                best_score = score
                promoted = {
                    'representative_id': row['representative_id'],
                    'condition_id': row['condition_id'],
                }
        elif row['delayed_response_advantage'] >= DELAY_ADVANTAGE_MIN:
            row['output_label'] = 'transient history response'
    if promoted is not None:
        for row in rows:
            if row['representative_id'] == promoted['representative_id'] and row['condition_id'] == promoted['condition_id']:
                row['promoted_followup'] = 1
    return rows, promoted


def add_refinement_label(base_row: dict[str, Any], refined_row: dict[str, Any]) -> None:
    if (
        float(refined_row['composite_lifetime']) >= float(base_row['composite_lifetime'])
        and float(refined_row['binding_persistence_score']) >= float(base_row['binding_persistence_score'])
        and float(refined_row['coarse_persistence_sigma4']) >= float(base_row['coarse_persistence_sigma4'])
    ):
        refined_row['output_label'] = 'scale-stable memory-assisted ordering'
    else:
        refined_row['output_label'] = 'scale-fragile memory-assisted ordering'


def create_summary_plots(work_plot_dir: Path, rows: list[dict[str, Any]], refined_rows: list[dict[str, Any]]) -> list[Path]:
    plot_paths: list[Path] = []
    base_rows = [row for row in rows if row['phase'] == 'base']
    reps = sorted({row['representative_id'] for row in base_rows})
    conds = [cond for cond in ['null_control', 'very_weak_memory_coupling', 'weak_memory_coupling'] if any(row['condition_id'] == cond for row in base_rows)]
    labels = {
        'null_control': 'null',
        'very_weak_memory_coupling': '0.01',
        'weak_memory_coupling': '0.02',
    }

    if not reps or not conds:
        return plot_paths

    comparison = work_plot_dir / 'stage18_memory_comparison.png'
    x = np.arange(len(reps))
    width = 0.22
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    for idx, cond in enumerate(conds):
        lifetimes = []
        coarse_pers = []
        for rep in reps:
            row = next((item for item in base_rows if item['representative_id'] == rep and item['condition_id'] == cond), None)
            lifetimes.append(float(row['composite_lifetime']) if row is not None else 0.0)
            coarse_pers.append(float(row['coarse_persistence_sigma4']) if row is not None else 0.0)
        axes[0].bar(x + (idx - 1) * width, lifetimes, width=width, label=labels[cond])
        axes[1].bar(x + (idx - 1) * width, coarse_pers, width=width, label=labels[cond])
    axes[0].set_title('Composite lifetime')
    axes[1].set_title('Coarse persistence (sigma=4)')
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(reps, rotation=20, ha='right')
        ax.grid(alpha=0.25, axis='y')
        ax.legend(fontsize=8)
    fig.suptitle('Stage 18 Branch-A null vs delayed coupling')
    fig.savefig(comparison, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(comparison)

    response = work_plot_dir / 'stage18_memory_response_matrix.png'
    label_map = {
        'no memory dynamic effect': 0,
        'transient history response': 1,
        'scale-fragile memory-assisted ordering': 2,
        'scale-stable memory-assisted ordering': 3,
    }
    data = np.zeros((len(reps), len(conds)), dtype=float)
    for i, rep in enumerate(reps):
        for j, cond in enumerate(conds):
            row = next((item for item in base_rows if item['representative_id'] == rep and item['condition_id'] == cond), None)
            if row is not None:
                data[i, j] = label_map.get(str(row['output_label']), 0)
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    im = ax.imshow(data, cmap='viridis', vmin=0, vmax=3)
    ax.set_xticks(range(len(conds)))
    ax.set_xticklabels([labels[cond] for cond in conds])
    ax.set_yticks(range(len(reps)))
    ax.set_yticklabels(reps)
    ax.set_title('Stage 18 Branch-A response matrix')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(response, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(response)

    if refined_rows:
        refinement = work_plot_dir / 'stage18_memory_refinement.png'
        labels_x = [f"{row['condition_id']}\nn={row['resolution']}" for row in refined_rows]
        values = [float(row['composite_lifetime']) for row in refined_rows]
        fig, ax = plt.subplots(figsize=(7.2, 4.4))
        ax.bar(range(len(labels_x)), values, color='tab:purple')
        ax.set_xticks(range(len(labels_x)))
        ax.set_xticklabels(labels_x)
        ax.set_title('Composite lifetime across refinement')
        ax.grid(alpha=0.25, axis='y')
        fig.savefig(refinement, dpi=180, bbox_inches='tight')
        plt.close(fig)
        plot_paths.append(refinement)

    return plot_paths


def write_note(path: Path, json_rel: str, csv_rel: str, base_rows: list[dict[str, Any]], refined_rows: list[dict[str, Any]], label_counts: dict[str, int], promoted: dict[str, Any] | None) -> None:
    lines = [
        '# Stage 18 Emergent Memory Retarded Interaction Dynamics v1',
        '',
        f'Timestamped summary: `{json_rel}`',
        f'Timestamped run table: `{csv_rel}`',
        '',
        f'Output-label counts: {label_counts}',
        '',
        'This Branch-A-first pilot tests whether a single-lag retarded feedback term can create collective organization that instantaneous weak couplings could not.',
        '',
        'Base comparisons:',
        '',
    ]
    for row in base_rows:
        lines.extend([
            f"- `{row['representative_id']}` / `{row['condition_id']}`",
            f"  - composite lifetime delta: `{row['composite_lifetime_delta']:.4f}`",
            f"  - binding persistence delta: `{row['binding_persistence_delta']:.4f}`",
            f"  - coarse persistence delta: `{row['coarse_persistence_delta']:.4f}`",
            f"  - delayed response advantage: `{row['delayed_response_advantage']:.4f}`",
            f"  - gate met: `{row['gate_met']}`",
            f"  - label: `{row['output_label']}`",
        ])
    lines.extend(['', f'Promoted follow-up: `{promoted}`' if promoted else 'Promoted follow-up: `none`', ''])
    if refined_rows:
        lines.append('Refinement follow-up:')
        lines.append('')
        for row in refined_rows:
            lines.extend([
                f"- `{row['representative_id']}` / `{row['condition_id']}` / `n={row['resolution']}`",
                f"  - label: `{row['output_label']}`",
            ])
        lines.append('')
    lines.extend([
        'Interpretation boundary:',
        '- any positive effect is read only as history-sensitive ordering',
        '- no force law or geometry claim is made here',
        '',
    ])
    path.write_text('\n'.join(lines), encoding='utf-8')


def main() -> None:
    args = parse_args()
    defaults = load_stage10_defaults()
    runsheet = load_runsheet(args.runsheet)
    base = runsheet['base_seed_reference']
    rep_lookup = lookup_by_id(runsheet['representatives'], 'representative_id')
    cond_lookup = lookup_by_id(runsheet['conditions'], 'condition_id')
    runs = selected_runs(runsheet['runs'], args.run_ids, args.memory_class_id)
    work_plot_dir = Path('/tmp') / 'haos_iip_stage18_memory'
    work_plot_dir.mkdir(parents=True, exist_ok=True)

    base_rows: list[dict[str, Any]] = []
    base_payloads: dict[tuple[str, str], dict[str, Any]] = {}
    plot_paths: list[Path] = []

    for run in runs:
        rep = rep_lookup[run['representative_id']]
        condition = cond_lookup[run['condition_id']]
        payload, run_plots = simulate_run(run, rep, condition, runsheet, defaults, base, work_plot_dir)
        base_rows.append(payload['summary'])
        base_payloads[(run['representative_id'], run['condition_id'])] = payload
        plot_paths.extend(run_plots)

    base_rows, promoted = classify_base_rows(base_rows)
    refined_rows: list[dict[str, Any]] = []
    if promoted is not None and not args.skip_followup:
        base_row = next(row for row in base_rows if row['representative_id'] == promoted['representative_id'] and row['condition_id'] == promoted['condition_id'])
        null_row = next(row for row in base_rows if row['representative_id'] == promoted['representative_id'] and row['condition_id'] == 'null_control')
        for cond_id in ['null_control', promoted['condition_id']]:
            run = next(item for item in runs if item['representative_id'] == promoted['representative_id'] and item['condition_id'] == cond_id)
            follow = dict(run)
            follow['resolution'] = int(runsheet['resolution_policy']['refined_resolution'])
            follow['run_id'] = f"{run['run_id'].replace('__n12', '')}__n24"
            rep = rep_lookup[follow['representative_id']]
            condition = cond_lookup[follow['condition_id']]
            payload, run_plots = simulate_run(follow, rep, condition, runsheet, defaults, base, work_plot_dir)
            row = payload['summary']
            row['phase'] = 'refinement'
            row['composite_lifetime_delta'] = float(row['composite_lifetime']) - float(null_row['composite_lifetime'])
            row['binding_persistence_delta'] = float(row['binding_persistence_score']) - float(null_row['binding_persistence_score'])
            row['coarse_persistence_delta'] = float(row['coarse_persistence_sigma4']) - float(null_row['coarse_persistence_sigma4'])
            row['delayed_response_advantage'] = float(row['delayed_response_correlation']) - float(null_row['delayed_response_correlation'])
            row['new_ordering_class'] = int(
                interaction_rank(str(row['interaction_label'])) > interaction_rank(str(null_row['interaction_label']))
                or coarse_rank(str(row['coarse_label'])) > coarse_rank(str(null_row['coarse_label']))
            )
            if cond_id == promoted['condition_id']:
                add_refinement_label(base_row, row)
            else:
                row['output_label'] = 'no memory dynamic effect'
            refined_rows.append(row)
            plot_paths.extend(run_plots)

    all_rows = base_rows + refined_rows
    label_counts = dict(Counter(row['output_label'] for row in all_rows))
    plot_paths.extend(create_summary_plots(work_plot_dir, all_rows, refined_rows))

    result = {
        'stage': 'Stage 18 Branch A',
        'description': runsheet['description'],
        'memory_class_id': args.memory_class_id,
        'base_seed_reference': base,
        'base_runs': [{'run_id': row['run_id'], 'summary': row} for row in base_rows],
        'refinement_runs': [{'run_id': row['run_id'], 'summary': row} for row in refined_rows],
        'output_label_counts': label_counts,
        'promoted_followup': promoted,
        'conclusion': 'Stage 18 Branch A tests whether a single-lag retarded self-feedback term can generate collective organization on the frozen architecture.',
    }

    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug='stage18_memory_retarded',
        result=result,
        csv_rows=all_rows,
        csv_fieldnames=CSV_FIELDS,
        plot_paths=plot_paths,
    )
    note_path = ATLAS_NOTES / runsheet['note_name']
    write_note(note_path, str(json_path.relative_to(REPO_ROOT)), str(csv_path.relative_to(REPO_ROOT)), base_rows, refined_rows, label_counts, promoted)
    if args.append_log:
        append_log(
            title=f"Stage 18 Memory Retarded Probe ({json_path.stem})",
            config_summary=f"runs={len(base_rows)}, promoted_followup={promoted}",
            result_path=json_path,
            stamped_plots=stamped_plots,
            observation=f"output_label_counts={label_counts}",
            conclusion='the Stage 18 Branch-A-first pilot tests whether single-lag retarded feedback creates memory-assisted collective organization',
        )

    print(f'wrote {json_path}')
    print(f'wrote {csv_path}')
    for item in stamped_plots:
        print(item)
    print(f'output_label_counts={label_counts}')
    print(f'promoted_followup={promoted}')


if __name__ == '__main__':
    main()
