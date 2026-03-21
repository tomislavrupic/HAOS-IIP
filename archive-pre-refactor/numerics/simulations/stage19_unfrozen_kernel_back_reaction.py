#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
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

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage19_unfrozen_kernel_runs.json'
SETUP_CACHE: dict[tuple[int, str, float, float], tuple[Any, Any]] = {}
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage19_kernel'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

ANALYSIS_STRIDE = 4
SOURCE_SIGMA = 4.0
COMPOSITE_GAIN_MIN = 0.05
BINDING_GAIN_MIN = 0.03
COARSE_GAIN_MIN = 0.02

CSV_FIELDS = [
    'run_id',
    'representative_id',
    'condition_id',
    'eps',
    'tau_relax',
    'resolution',
    'interaction_label',
    'coarse_label',
    'composite_lifetime',
    'binding_persistence_score',
    'coarse_persistence_sigma4',
    'kernel_deviation_norm',
    'max_local_kernel_deformation',
    'substrate_relaxation_memory_score',
    'local_stress_persistence',
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
    parser = argparse.ArgumentParser(description='Run the Stage 19A unfrozen-kernel back-reaction pilot.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    parser.add_argument('--append-log', action='store_true')
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


def recover_unweighted_incidence(data: Any) -> tuple[sp.csr_matrix, sp.csr_matrix, sp.csr_matrix, np.ndarray]:
    edge_scale = np.sqrt(np.maximum(np.asarray(data.edge_weights, dtype=float), 1.0e-12))
    face_scale = np.sqrt(np.maximum(np.asarray(data.face_weights, dtype=float), 1.0e-12))
    b0 = (sp.diags(1.0 / edge_scale) @ data.d0).tocsr()
    c = (sp.diags(1.0 / face_scale) @ data.d1).tocsr()
    abs_c = abs(c).tocsr()
    face_counts = np.maximum(np.asarray(abs_c.sum(axis=1)).ravel(), 1.0)
    return b0, c, abs_c, face_counts


def dynamic_operator(b0: sp.csr_matrix, c: sp.csr_matrix, abs_c: sp.csr_matrix, face_counts: np.ndarray, edge_weights: np.ndarray) -> tuple[sp.csr_matrix, np.ndarray]:
    edge_weights = np.asarray(edge_weights, dtype=float)
    face_weights = np.asarray(abs_c @ edge_weights, dtype=float).ravel() / face_counts
    d0 = (sp.diags(np.sqrt(np.maximum(edge_weights, 1.0e-12))) @ b0).tocsr()
    d1 = (sp.diags(np.sqrt(np.maximum(face_weights, 1.0e-12))) @ c).tocsr()
    operator = (d0 @ d0.T + d1.T @ d1).tocsr()
    return operator, face_weights


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


def render_grid_slice(ax, grid: np.ndarray, title: str, cmap: str = 'viridis') -> None:
    slice_idx = grid.shape[2] // 2
    img = grid[:, :, slice_idx].T
    im = ax.imshow(img, origin='lower', cmap=cmap)
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)


def edge_energy_profile(midpoints: np.ndarray, q: np.ndarray, n_side: int, sigma: float) -> np.ndarray:
    grid = edge_field_to_grid(midpoints, np.abs(q) ** 2, n_side)
    env = np.maximum(gaussian_smooth_periodic(grid, sigma), 0.0)
    coords = np.floor(midpoints * n_side).astype(int) % n_side
    return env[coords[:, 0], coords[:, 1], coords[:, 2]].astype(float)


def estimate_crossing_time(base: dict[str, Any], dt: float) -> float:
    width_time = 2.0 * float(base['bandwidth']) / max(float(base['central_k']), 1.0e-12)
    return max(width_time, 4.0 * dt)


def simulate_run(
    run: dict[str, Any],
    rep: dict[str, Any],
    condition: dict[str, Any],
    runsheet: dict[str, Any],
    defaults: dict[str, Any],
    base: dict[str, Any],
) -> tuple[dict[str, Any], list[Path]]:
    resolution = int(run['resolution'])
    boundary_type = str(base['boundary_type'])
    data, projector = get_cached_setup(resolution, defaults, boundary_type)
    b0, c, abs_c, face_counts = recover_unweighted_incidence(data)
    kick_axis = int(defaults['kick_axis'])
    bandwidth = float(base['bandwidth'])
    eps = float(run['eps'])
    kappa_max = float(runsheet['parameter_policy']['kappa_max'])

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

    baseline_edge_weights = np.asarray(data.edge_weights, dtype=float)
    baseline_operator = data.L1
    dt = suggested_dt(baseline_operator, float(defaults['dt_safety']))
    steps = max(2, int(np.ceil(float(defaults['t_final']) / dt)))
    tau_relax = float(condition['tau_relax_multiple']) * estimate_crossing_time(base, dt)
    tau_relax = max(tau_relax, 2.0 * dt)
    decay = max(0.0, 1.0 - dt / tau_relax)
    leakage_fn = make_leakage_fn(projector)

    stress = np.zeros_like(baseline_edge_weights)
    current_edge_weights = baseline_edge_weights.copy()

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
    kernel_deviation_series: list[float] = []
    max_deformation_series: list[float] = []
    stress_norm_series: list[float] = []
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

    def record(t: float, operator: sp.csr_matrix) -> None:
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

        relative_deformation = np.abs(current_edge_weights / np.maximum(baseline_edge_weights, 1.0e-12) - 1.0)
        kernel_deviation_series.append(float(np.linalg.norm(current_edge_weights - baseline_edge_weights) / max(np.linalg.norm(baseline_edge_weights), 1.0e-12)))
        max_deformation_series.append(float(np.max(relative_deformation)))
        stress_norm_series.append(float(np.linalg.norm(stress) / math.sqrt(max(len(stress), 1))))

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

    record(0.0, baseline_operator)
    for step_idx in range(steps):
        total_q = np.sum(packet_qs, axis=0)
        if eps > 0.0:
            rho = edge_energy_profile(data.midpoints, total_q, resolution, SOURCE_SIGMA)
            rho_centered = rho - float(np.mean(rho))
            stress = decay * stress + rho_centered
            stress_cap = kappa_max / max(eps, 1.0e-12)
            stress = np.clip(stress, -stress_cap, stress_cap)
            scale = np.clip(1.0 + eps * stress, 1.0 - kappa_max, 1.0 + kappa_max)
            current_edge_weights = baseline_edge_weights * scale
        else:
            stress[:] = 0.0
            current_edge_weights = baseline_edge_weights.copy()

        operator, _face_weights = dynamic_operator(b0, c, abs_c, face_counts, current_edge_weights)
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
            record((step_idx + 1) * dt, operator)

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
    peak_stress = max(max(stress_norm_series), 1.0e-12)
    substrate_memory_score = float(stress_norm_series[-1] / peak_stress) if stress_norm_series else 0.0
    local_stress_persistence = float(np.mean((np.asarray(stress_norm_series) >= 0.5 * peak_stress).astype(float))) if stress_norm_series else 0.0

    row = {
        'run_id': run['run_id'],
        'representative_id': run['representative_id'],
        'condition_id': run['condition_id'],
        'eps': eps,
        'tau_relax': tau_relax,
        'resolution': resolution,
        'interaction_label': interaction_label,
        'coarse_label': coarse_label,
        'composite_lifetime': float(interaction_summary['composite_lifetime']),
        'binding_persistence_score': binding_persistence,
        'coarse_persistence_sigma4': coarse_summary['coarse_persistence_sigma4'],
        'kernel_deviation_norm': float(np.max(kernel_deviation_series)) if kernel_deviation_series else 0.0,
        'max_local_kernel_deformation': float(np.max(max_deformation_series)) if max_deformation_series else 0.0,
        'substrate_relaxation_memory_score': substrate_memory_score,
        'local_stress_persistence': local_stress_persistence,
        'new_ordering_class': 0,
        'composite_lifetime_delta': 0.0,
        'binding_persistence_delta': 0.0,
        'coarse_persistence_delta': 0.0,
        'gate_met': 0,
        'promoted_followup': 0,
        'output_label': 'no kernel back-reaction effect',
        'constraint_max': float(np.max(constraint_norms)) if constraint_norms else 0.0,
        'sector_leakage': float(np.max(sector_leakages)) if sector_leakages else 0.0,
        'notes': rep['notes'],
    }

    plot_paths: list[Path] = []
    trace_path = WORK_PLOT_DIR / f"stage19_{run['run_id']}_pair_distance_overlap.png"
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))
    axes[0].plot(times, pair_distance_series, color='tab:blue')
    axes[0].set_title('Mean pair distance')
    axes[0].set_xlabel('time')
    axes[0].grid(alpha=0.25)
    axes[1].plot(times, overlap_series, color='tab:orange')
    axes[1].set_title('Mean overlap')
    axes[1].set_xlabel('time')
    axes[1].grid(alpha=0.25)
    fig.suptitle(f"Stage 19A: {run['run_id']}")
    fig.savefig(trace_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(trace_path)

    basin_path = WORK_PLOT_DIR / f"stage19_{run['run_id']}_basin_kernel.png"
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3))
    axes[0].plot(times, coarse[4.0]['dominant_area_fractions'], color='tab:green')
    axes[0].set_title('Dominant basin area (sigma=4)')
    axes[0].set_xlabel('time')
    axes[0].grid(alpha=0.25)
    axes[1].plot(times, max_deformation_series, color='tab:red')
    axes[1].set_title('Max local kernel deformation')
    axes[1].set_xlabel('time')
    axes[1].grid(alpha=0.25)
    fig.savefig(basin_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(basin_path)

    field_path = WORK_PLOT_DIR / f"stage19_{run['run_id']}_field_snapshot.png"
    create_field_snapshot(field_path, run['run_id'], data.midpoints, np.sum(packet_qs, axis=0), boundary_type, int(defaults['slice_axis']))
    plot_paths.append(field_path)

    env = coarse[4.0]['last_env']
    if env is not None:
        env_path = WORK_PLOT_DIR / f"stage19_{run['run_id']}_envelope_snapshot.png"
        fig, ax = plt.subplots(figsize=(5.4, 4.6))
        render_grid_slice(ax, env, f"Envelope sigma=4: {run['run_id']}")
        fig.savefig(env_path, dpi=180, bbox_inches='tight')
        plt.close(fig)
        plot_paths.append(env_path)

    result = {
        'run': run,
        'summary': row,
        'times': times,
        'pair_distance_series': pair_distance_series,
        'overlap_series': overlap_series,
        'kernel_deviation_series': kernel_deviation_series,
        'max_deformation_series': max_deformation_series,
        'stress_norm_series': stress_norm_series,
    }
    return result, plot_paths


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
        row['gate_met'] = int(gate_hit)
        if gate_hit:
            row['output_label'] = 'transient substrate response'
            score = (
                float(row['composite_lifetime_delta'])
                + float(row['binding_persistence_delta'])
                + float(row['coarse_persistence_delta'])
            )
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
    return rows, promoted


def create_summary_plots(rows: list[dict[str, Any]]) -> list[Path]:
    plot_paths: list[Path] = []
    reps = sorted({row['representative_id'] for row in rows})
    conds = ['null_control', 'very_weak_backreaction', 'weak_backreaction']
    labels = {
        'null_control': 'null',
        'very_weak_backreaction': '0.01',
        'weak_backreaction': '0.02',
    }

    comparison = WORK_PLOT_DIR / 'stage19_kernel_comparison.png'
    x = np.arange(len(reps))
    width = 0.22
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    for idx, cond in enumerate(conds):
        lifetimes = []
        coarse_pers = []
        for rep in reps:
            row = next(item for item in rows if item['representative_id'] == rep and item['condition_id'] == cond)
            lifetimes.append(float(row['composite_lifetime']))
            coarse_pers.append(float(row['coarse_persistence_sigma4']))
        axes[0].bar(x + (idx - 1) * width, lifetimes, width=width, label=labels[cond])
        axes[1].bar(x + (idx - 1) * width, coarse_pers, width=width, label=labels[cond])
    axes[0].set_title('Composite lifetime')
    axes[1].set_title('Coarse persistence (sigma=4)')
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(reps, rotation=20, ha='right')
        ax.grid(alpha=0.25, axis='y')
        ax.legend(fontsize=8)
    fig.suptitle('Stage 19A null vs kernel back-reaction')
    fig.savefig(comparison, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(comparison)

    kernel_plot = WORK_PLOT_DIR / 'stage19_kernel_deformation_summary.png'
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    for idx, cond in enumerate(conds):
        deviations = []
        stress_pers = []
        for rep in reps:
            row = next(item for item in rows if item['representative_id'] == rep and item['condition_id'] == cond)
            deviations.append(float(row['kernel_deviation_norm']))
            stress_pers.append(float(row['local_stress_persistence']))
        axes[0].bar(x + (idx - 1) * width, deviations, width=width, label=labels[cond])
        axes[1].bar(x + (idx - 1) * width, stress_pers, width=width, label=labels[cond])
    axes[0].set_title('Kernel deviation norm')
    axes[1].set_title('Local stress persistence')
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(reps, rotation=20, ha='right')
        ax.grid(alpha=0.25, axis='y')
        ax.legend(fontsize=8)
    fig.savefig(kernel_plot, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(kernel_plot)

    response = WORK_PLOT_DIR / 'stage19_kernel_response_matrix.png'
    label_map = {
        'no kernel back-reaction effect': 0,
        'transient substrate response': 1,
        'scale-fragile substrate-assisted persistence': 2,
        'scale-stable substrate-assisted persistence': 3,
    }
    data = np.zeros((len(reps), len(conds)), dtype=float)
    for i, rep in enumerate(reps):
        for j, cond in enumerate(conds):
            row = next(item for item in rows if item['representative_id'] == rep and item['condition_id'] == cond)
            data[i, j] = label_map.get(str(row['output_label']), 0)
    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    im = ax.imshow(data, cmap='viridis', vmin=0, vmax=3)
    ax.set_xticks(range(len(conds)))
    ax.set_xticklabels([labels[cond] for cond in conds])
    ax.set_yticks(range(len(reps)))
    ax.set_yticklabels(reps)
    ax.set_title('Stage 19A response matrix')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(response, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(response)
    return plot_paths


def write_note(path: Path, json_rel: str, csv_rel: str, rows: list[dict[str, Any]], label_counts: dict[str, int], promoted: dict[str, Any] | None) -> None:
    lines = [
        '# Stage 19 Unfrozen Kernel Back-Reaction v1',
        '',
        f'Timestamped summary: `{json_rel}`',
        f'Timestamped run table: `{csv_rel}`',
        '',
        f'Output-label counts: {label_counts}',
        '',
        'This Stage 19A pilot tests whether bounded local energy-density back-reaction on the kernel edge weights can create the first substrate-assisted persistence effect.',
        '',
        'Base comparisons:',
        '',
    ]
    for row in rows:
        lines.extend([
            f"- `{row['representative_id']}` / `{row['condition_id']}`",
            f"  - composite lifetime delta: `{row['composite_lifetime_delta']:.4f}`",
            f"  - binding persistence delta: `{row['binding_persistence_delta']:.4f}`",
            f"  - coarse persistence delta: `{row['coarse_persistence_delta']:.4f}`",
            f"  - kernel deviation norm: `{row['kernel_deviation_norm']:.4f}`",
            f"  - max local kernel deformation: `{row['max_local_kernel_deformation']:.4f}`",
            f"  - gate met: `{row['gate_met']}`",
            f"  - label: `{row['output_label']}`",
        ])
    lines.extend([
        '',
        f'Promoted follow-up: `{promoted}`' if promoted else 'Promoted follow-up: `none`',
        '',
        'Interpretation boundary:',
        '- any positive effect is read only as bounded substrate-assisted persistence',
        '- no geometric or force-law claim is made here',
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
    runs = selected_runs(runsheet['runs'], args.run_ids)
    if not runs:
        raise SystemExit('no Stage 19A runs selected')

    rows: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    for run in runs:
        rep = rep_lookup[run['representative_id']]
        condition = cond_lookup[run['condition_id']]
        payload, run_plots = simulate_run(run, rep, condition, runsheet, defaults, base)
        rows.append(payload['summary'])
        plot_paths.extend(run_plots)

    rows, promoted = classify_base_rows(rows)
    label_counts = dict(Counter(row['output_label'] for row in rows))
    plot_paths.extend(create_summary_plots(rows))

    result = {
        'stage': 'Stage 19A',
        'description': runsheet['description'],
        'kernel_backreaction': runsheet['kernel_backreaction'],
        'parameter_policy': runsheet['parameter_policy'],
        'base_seed_reference': base,
        'base_runs': [{'run_id': row['run_id'], 'summary': row} for row in rows],
        'output_label_counts': label_counts,
        'promoted_followup': promoted,
        'conclusion': 'Stage 19A tests whether minimal bounded kernel back-reaction from local energy density can create substrate-assisted persistence on the projected transverse architecture.',
    }

    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug='stage19_unfrozen_kernel',
        result=result,
        csv_rows=rows,
        csv_fieldnames=CSV_FIELDS,
        plot_paths=plot_paths,
    )
    note_path = ATLAS_NOTES / runsheet['note_name']
    write_note(note_path, str(json_path.relative_to(REPO_ROOT)), str(csv_path.relative_to(REPO_ROOT)), rows, label_counts, promoted)

    print(f'wrote {json_path}')
    print(f'wrote {csv_path}')
    for item in stamped_plots:
        print(item)
    print(f'output_label_counts={label_counts}')
    print(f'promoted_followup={promoted}')


if __name__ == '__main__':
    main()
