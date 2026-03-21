#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np

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

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage21_mutual_registration_runs.json'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage21_mutual_registration'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)
SETUP_CACHE: dict[tuple[int, str, float, float], tuple[Any, Any]] = {}

ANALYSIS_STRIDE = 4
COMPOSITE_GAIN_MIN = 0.05
BINDING_GAIN_MIN = 0.03
COARSE_GAIN_MIN = 0.02
CORRIDOR_GAIN_MIN = 0.05
ACTIVATION_THRESHOLD = 0.25

CSV_FIELDS = [
    'run_id',
    'representative_id',
    'registration_class',
    'condition_id',
    'eps',
    'resolution',
    'interaction_label',
    'coarse_label',
    'composite_lifetime',
    'binding_persistence_score',
    'coarse_persistence_sigma4',
    'mean_registration_strength',
    'max_registration_strength',
    'registration_duty_cycle',
    'effective_interaction_energy_fraction',
    'encounter_dwell_mean',
    'motif_survival_time',
    'corridor_dwell_time',
    'transport_channel_locking',
    'local_transport_anisotropy_index',
    'trajectory_curvature_proxy',
    'energy_dispersion_shift',
    'new_ordering_class',
    'composite_lifetime_delta',
    'binding_persistence_delta',
    'coarse_persistence_delta',
    'corridor_dwell_delta',
    'gate_met',
    'promoted_followup',
    'output_label',
    'constraint_max',
    'sector_leakage',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage 21 mutual registration pilot.')
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


def is_null_condition(condition_id: str) -> bool:
    return condition_id == 'null'


def condition_order(condition_id: str) -> int:
    return {'null': 0, 'very_weak': 1, 'weak': 2}.get(condition_id, 99)


def condition_label(condition_id: str) -> str:
    return {'null': 'null', 'very_weak': '0.01', 'weak': '0.02'}.get(condition_id, condition_id)


def class_order(registration_class: str) -> int:
    return {
        'phase_coherence': 0,
        'spatial_configuration': 1,
        'composite_topology': 2,
    }.get(registration_class, 99)


def class_label(registration_class: str) -> str:
    return {
        'phase_coherence': 'phase',
        'spatial_configuration': 'spatial',
        'composite_topology': 'topology',
    }.get(registration_class, registration_class)


def cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    denom = max(float(np.linalg.norm(a)) * float(np.linalg.norm(b)), 1.0e-12)
    return float(np.clip(np.vdot(a, b).real / denom, -1.0, 1.0))


def pair_overlap(qi: np.ndarray, qj: np.ndarray) -> float:
    wi = np.abs(qi) ** 2
    wj = np.abs(qj) ** 2
    denom = math.sqrt(max(float(np.sum(wi)) * float(np.sum(wj)), 1.0e-12))
    return float(np.sum(np.sqrt(wi * wj)) / denom)


def pairwise_overlap_mean(component_qs: list[np.ndarray]) -> float:
    values = [pair_overlap(component_qs[i], component_qs[j]) for i, j in combinations(range(len(component_qs)), 2)]
    return float(np.mean(values)) if values else 0.0


def pairwise_distance_stats(centers: list[np.ndarray], boundary_type: str) -> tuple[float, float, list[float]]:
    values: list[float] = []
    for i in range(len(centers)):
        for j in range(i + 1, len(centers)):
            delta = displacement(np.asarray([centers[j]], dtype=float), np.asarray(centers[i], dtype=float), boundary_type)[0]
            values.append(float(np.linalg.norm(delta)))
    if not values:
        return 0.0, 0.0, []
    arr = np.asarray(values, dtype=float)
    return float(np.mean(arr)), float(np.min(arr)), values


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


def trajectory_curvature_proxy(centers: np.ndarray, boundary_type: str) -> float:
    if len(centers) < 3:
        return 0.0
    turns: list[float] = []
    for idx in range(1, len(centers) - 1):
        prev_vec = displacement(centers[idx][None, :], centers[idx - 1], boundary_type)[0]
        next_vec = displacement(centers[idx + 1][None, :], centers[idx], boundary_type)[0]
        prev_norm = float(np.linalg.norm(prev_vec))
        next_norm = float(np.linalg.norm(next_vec))
        if prev_norm <= 1.0e-12 or next_norm <= 1.0e-12:
            continue
        cosine = float(np.clip(np.dot(prev_vec, next_vec) / (prev_norm * next_norm), -1.0, 1.0))
        turns.append(float(np.arccos(cosine)))
    return float(np.mean(turns)) if turns else 0.0


def render_grid_slice(ax, grid: np.ndarray, title: str, cmap: str = 'viridis') -> None:
    slice_idx = grid.shape[2] // 2
    img = grid[:, :, slice_idx].T
    im = ax.imshow(img, origin='lower', cmap=cmap)
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)


def mean_upper_triangle(mat: np.ndarray) -> float:
    if mat.shape[0] < 2:
        return 0.0
    vals = mat[np.triu_indices(mat.shape[0], k=1)]
    return float(np.mean(vals)) if vals.size else 0.0


def max_upper_triangle(mat: np.ndarray) -> float:
    if mat.shape[0] < 2:
        return 0.0
    vals = mat[np.triu_indices(mat.shape[0], k=1)]
    return float(np.max(vals)) if vals.size else 0.0


def mean_dwell_duration(active: list[float], dt: float) -> float:
    if not active:
        return 0.0
    durations: list[int] = []
    current = 0
    for flag in active:
        if flag >= 0.5:
            current += 1
        elif current > 0:
            durations.append(current)
            current = 0
    if current > 0:
        durations.append(current)
    if not durations:
        return 0.0
    return float(np.mean(durations) * dt)


def build_centers(data: Any, packet_qs: list[np.ndarray], boundary_type: str) -> list[np.ndarray]:
    centers: list[np.ndarray] = []
    for comp_q in packet_qs:
        comp_weights = np.abs(comp_q) ** 2
        centers.append(weighted_center(data.midpoints, comp_weights, boundary_type))
    return centers


def compute_registration(packet_qs: list[np.ndarray], packet_vs: list[np.ndarray], centers: list[np.ndarray], boundary_type: str, registration_class: str, params: dict[str, Any]) -> tuple[np.ndarray, float]:
    n = len(packet_qs)
    reg = np.zeros((n, n), dtype=float)
    pair_distances: list[float] = []
    pair_alignments: list[float] = []
    pair_overlaps: list[float] = []
    pair_indices: list[tuple[int, int]] = []

    sigma_phi = float(params.get('phase_sigma', 0.35))
    sigma_vel = float(params.get('velocity_sigma', 0.4))
    sigma_d = float(params.get('distance_sigma', 0.25))

    for i in range(n):
        for j in range(i + 1, n):
            align = 0.5 * (1.0 + cosine_similarity(packet_qs[i], packet_qs[j]))
            vel_align = 0.5 * (1.0 + cosine_similarity(packet_vs[i], packet_vs[j]))
            overlap = pair_overlap(packet_qs[i], packet_qs[j])
            delta = displacement(np.asarray([centers[j]], dtype=float), np.asarray(centers[i], dtype=float), boundary_type)[0]
            dist = float(np.linalg.norm(delta))
            pair_indices.append((i, j))
            pair_distances.append(dist)
            pair_alignments.append(align)
            pair_overlaps.append(overlap)
            if registration_class == 'phase_coherence':
                diff = 1.0 - align
                vel_mismatch = 1.0 - vel_align
                value = math.exp(-((diff / max(sigma_phi, 1.0e-12)) ** 2 + (vel_mismatch / max(sigma_vel, 1.0e-12)) ** 2))
            elif registration_class == 'spatial_configuration':
                value = math.exp(-((dist / max(sigma_d, 1.0e-12)) ** 2)) * overlap
            else:
                value = 0.0
            reg[i, j] = value
            reg[j, i] = value

    motif_strength = 0.0
    if registration_class == 'composite_topology' and pair_indices:
        mean_align = float(np.mean(pair_alignments))
        mean_overlap = float(np.mean(pair_overlaps))
        if len(pair_distances) >= 2:
            symmetry = math.exp(-(float(np.std(pair_distances)) / max(sigma_d, 1.0e-12)) ** 2)
        else:
            symmetry = math.exp(-((pair_distances[0] - sigma_d) / max(sigma_d, 1.0e-12)) ** 2)
        motif_strength = float(np.clip(symmetry * (0.5 * mean_align + 0.5 * mean_overlap), 0.0, 1.0))
        for (i, j), align, overlap in zip(pair_indices, pair_alignments, pair_overlaps):
            reg[i, j] = reg[j, i] = float(np.clip((0.5 * align + 0.5 * overlap) * motif_strength, 0.0, 1.0))
    else:
        motif_strength = mean_upper_triangle(reg)

    reg = np.clip(reg, 0.0, 1.0)
    return reg, float(motif_strength)


def compute_coupling_terms(packet_qs: list[np.ndarray], packet_vs: list[np.ndarray], data: Any, projector: Any, boundary_type: str, registration_class: str, params: dict[str, Any]) -> tuple[np.ndarray, list[np.ndarray], float]:
    centers = build_centers(data, packet_qs, boundary_type)
    reg_matrix, motif_strength = compute_registration(packet_qs, packet_vs, centers, boundary_type, registration_class, params)
    terms: list[np.ndarray] = []
    for i in range(len(packet_qs)):
        term = np.zeros_like(packet_qs[i])
        for j in range(len(packet_qs)):
            if i == j:
                continue
            term += float(reg_matrix[i, j]) * np.asarray(projector(packet_qs[j] - packet_qs[i]), dtype=float)
        terms.append(term)
    return reg_matrix, terms, motif_strength


def simulate_run(run: dict[str, Any], rep: dict[str, Any], runsheet: dict[str, Any], defaults: dict[str, Any], base: dict[str, Any]) -> tuple[dict[str, Any], list[Path]]:
    resolution = int(run['resolution'])
    boundary_type = str(base['boundary_type'])
    registration_class = str(run['registration_class'])
    params = runsheet.get('parameter_policy', {})
    coupling_law = runsheet.get('coupling_law', {})
    clamp_ratio = float(coupling_law.get('strict_clamp', 0.25))
    eps = float(run['eps'])
    data, projector = get_cached_setup(resolution, defaults, boundary_type)
    operator = data.L1
    kick_axis = int(defaults['kick_axis'])
    bandwidth = float(base['bandwidth'])
    leakage_fn = make_leakage_fn(projector)

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

    dt = suggested_dt(operator, float(defaults['dt_safety']))
    steps = max(2, int(np.ceil(float(defaults['t_final']) / dt)))

    times: list[float] = []
    total_centers: list[list[float]] = []
    widths: list[float] = []
    energies: list[float] = []
    sector_leakages: list[float] = []
    constraint_norms: list[float] = []
    pair_distance_series: list[float] = []
    overlap_series: list[float] = []
    registration_mean_series: list[float] = []
    registration_max_series: list[float] = []
    registration_active_series: list[float] = []
    effective_fraction_series: list[float] = []
    motif_strength_series: list[float] = []
    anisotropy_series: list[float] = []
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

    def record(t: float, reg_matrix: np.ndarray, motif_strength: float, interaction_fraction: float) -> None:
        total_q = np.sum(packet_qs, axis=0)
        total_v = np.sum(packet_vs, axis=0)
        total_weights = np.abs(total_q) ** 2
        total_center = weighted_center(data.midpoints, total_weights, boundary_type)
        times.append(float(t))
        total_centers.append(total_center.tolist())
        widths.append(packet_width(data.midpoints, total_weights, total_center, boundary_type))
        energies.append(state_energy(operator, total_q, total_v))
        sector_leakages.append(float(leakage_fn(total_q)))
        constraint_norms.append(float(np.linalg.norm(data.d0.T @ total_q)))
        anisotropy_series.append(anisotropy_ratio(data.midpoints, total_weights, total_center, boundary_type))

        comp_centers = build_centers(data, packet_qs, boundary_type)
        mean_pair, _min_pair, _pair_values = pairwise_distance_stats(comp_centers, boundary_type)
        pair_distance_series.append(mean_pair)
        overlap_series.append(pairwise_overlap_mean(packet_qs))
        reg_mean = mean_upper_triangle(reg_matrix)
        reg_max = max_upper_triangle(reg_matrix)
        registration_mean_series.append(reg_mean)
        registration_max_series.append(reg_max)
        registration_active_series.append(float(reg_mean >= ACTIVATION_THRESHOLD))
        effective_fraction_series.append(float(interaction_fraction))
        motif_strength_series.append(float(motif_strength))

        point_grid = edge_field_to_grid(data.midpoints, np.abs(total_q), resolution)
        for sigma in SIGMAS:
            env = np.maximum(gaussian_smooth_periodic(point_grid, sigma), 0.0)
            coarse[sigma]['last_env'] = env
            max_env = float(np.max(env))
            norm_env = env / max(max_env, 1.0e-12)
            coarse[sigma]['envelope_variances'].append(float(np.var(norm_env)))
            mask = env > (ALPHA * max_env) if max_env > 0.0 else np.zeros_like(env, dtype=bool)
            comps = connected_components_periodic(mask)
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
            coarse[sigma]['basin_counts'].append(float(len(comps)))
            coarse[sigma]['dominant_area_fractions'].append(dom_frac)
            coarse[sigma]['stable_flags'].append(stable)
            coarse[sigma]['dominant_masks'].append(dominant_mask)

    reg0, _terms0, motif0 = compute_coupling_terms(packet_qs, packet_vs, data, projector, boundary_type, registration_class, params)
    record(0.0, reg0, motif0, 0.0)
    baseline_energy = energies[0] if energies else 0.0

    for step_idx in range(steps):
        reg_matrix, raw_terms, motif_strength = compute_coupling_terms(packet_qs, packet_vs, data, projector, boundary_type, registration_class, params)
        base_accs: list[np.ndarray] = []
        coupling_accs: list[np.ndarray] = []
        for q, term in zip(packet_qs, raw_terms):
            base_acc = np.asarray(projector(-(operator @ q)), dtype=float)
            coupling = eps * np.asarray(projector(term), dtype=float)
            base_norm = max(float(np.linalg.norm(base_acc)), 1.0e-12)
            cap = clamp_ratio * base_norm
            coupling_norm = float(np.linalg.norm(coupling))
            if coupling_norm > cap:
                coupling *= cap / max(coupling_norm, 1.0e-12)
            base_accs.append(base_acc)
            coupling_accs.append(coupling)

        interaction_fraction = float(np.mean([
            float(np.linalg.norm(c)) / max(float(np.linalg.norm(b)), 1.0e-12)
            for b, c in zip(base_accs, coupling_accs)
        ])) if base_accs else 0.0

        v_halves: list[np.ndarray] = []
        q_preds: list[np.ndarray] = []
        for q, v, base_acc, coupling in zip(packet_qs, packet_vs, base_accs, coupling_accs):
            acc = base_acc + coupling
            v_half = np.asarray(projector(v + 0.5 * dt * acc), dtype=float)
            q_new = np.asarray(projector(q + dt * v_half), dtype=float)
            v_halves.append(v_half)
            q_preds.append(q_new)

        reg_new, raw_terms_new, _motif_new = compute_coupling_terms(q_preds, v_halves, data, projector, boundary_type, registration_class, params)
        next_qs: list[np.ndarray] = []
        next_vs: list[np.ndarray] = []
        for q_new, v_half, raw_term in zip(q_preds, v_halves, raw_terms_new):
            base_acc_new = np.asarray(projector(-(operator @ q_new)), dtype=float)
            coupling_new = eps * np.asarray(projector(raw_term), dtype=float)
            base_norm_new = max(float(np.linalg.norm(base_acc_new)), 1.0e-12)
            cap_new = clamp_ratio * base_norm_new
            coupling_norm_new = float(np.linalg.norm(coupling_new))
            if coupling_norm_new > cap_new:
                coupling_new *= cap_new / max(coupling_norm_new, 1.0e-12)
            v_new = np.asarray(projector(v_half + 0.5 * dt * (base_acc_new + coupling_new)), dtype=float)
            next_qs.append(q_new)
            next_vs.append(v_new)
        packet_qs = next_qs
        packet_vs = next_vs

        if ((step_idx + 1) % ANALYSIS_STRIDE == 0) or (step_idx == steps - 1):
            record((step_idx + 1) * dt, reg_new, motif_strength, interaction_fraction)

    sample_dt = float(np.median(np.diff(times))) if len(times) >= 2 else dt
    centers_arr = np.asarray(total_centers, dtype=float)
    center_shift = float(np.linalg.norm(displacement(centers_arr[-1][None, :], centers_arr[0], boundary_type)[0]))
    interaction_summary = {
        'packet_count': int(rep['packet_count']),
        'center_shift': center_shift,
        'width_ratio': float(widths[-1] / max(widths[0], 1.0e-12)),
        'max_mean_overlap': float(np.max(overlap_series)) if overlap_series else 0.0,
        'phase_lock_indicator': float(np.max(registration_mean_series)) if registration_mean_series else 0.0,
        'composite_lifetime': float(np.sum(np.asarray(overlap_series) >= 0.35) * sample_dt),
        'oscillation_count': 0,
        'initial_mean_pair_distance': float(pair_distance_series[0]) if pair_distance_series else 0.0,
        'min_mean_pair_distance': float(np.min(pair_distance_series)) if pair_distance_series else 0.0,
        'exchange_or_merger_flag': 'none',
        't_final': float(times[-1]) if times else 0.0,
    }
    interaction_label = classify_collective(interaction_summary)

    coarse_summary: dict[str, float] = {'t_final': float(times[-1]) if times else 0.0}
    for sigma in SIGMAS:
        coarse_summary[f'mean_basin_count_sigma{int(sigma)}'] = float(np.mean(coarse[sigma]['basin_counts'])) if coarse[sigma]['basin_counts'] else 0.0
        coarse_summary[f'mean_dominant_area_fraction_sigma{int(sigma)}'] = float(np.mean(coarse[sigma]['dominant_area_fractions'])) if coarse[sigma]['dominant_area_fractions'] else 0.0
        coarse_summary[f'coarse_persistence_sigma{int(sigma)}'] = float(np.mean(coarse[sigma]['stable_flags'])) if coarse[sigma]['stable_flags'] else 0.0
        coarse_summary[f'mean_envelope_variance_sigma{int(sigma)}'] = float(np.mean(coarse[sigma]['envelope_variances'])) if coarse[sigma]['envelope_variances'] else 0.0
        coarse_summary[f'basin_lifetime_sigma{int(sigma)}'] = float(np.sum(np.asarray(coarse[sigma]['dominant_area_fractions']) >= 0.05) * sample_dt)
    coarse_label = classify_coarse(coarse_summary)

    binding_mask = (np.asarray(pair_distance_series) <= 0.25) & (np.asarray(overlap_series) >= 0.35)
    binding_persistence = float(np.mean(binding_mask.astype(float))) if binding_mask.size else 0.0
    registration_duty_cycle = float(np.mean(registration_active_series)) if registration_active_series else 0.0
    encounter_dwell_mean = mean_dwell_duration(registration_active_series, sample_dt)
    motif_survival_time = float(np.sum(np.asarray(motif_strength_series) >= 0.5) * sample_dt) if motif_strength_series else 0.0
    local_transport_anisotropy_index = float(np.mean(anisotropy_series)) if anisotropy_series else 0.0
    corridor_dwell_time = float(np.mean(((np.asarray(registration_active_series) >= 0.5) & (np.asarray(anisotropy_series) >= 0.1)).astype(float))) if registration_active_series else 0.0
    transport_locking = int(
        rep['representative_id'] == 'counter_propagating_corridor'
        and len(pair_distance_series) >= 3
        and float(np.std(pair_distance_series)) <= 0.02
        and corridor_dwell_time >= 0.5
        and local_transport_anisotropy_index >= 0.1
    )
    curvature_proxy = trajectory_curvature_proxy(centers_arr, boundary_type)
    energy_dispersion_shift = float((energies[-1] - baseline_energy) / max(abs(baseline_energy), 1.0e-12)) if energies else 0.0

    row = {
        'run_id': run['run_id'],
        'representative_id': run['representative_id'],
        'registration_class': registration_class,
        'condition_id': run['condition_id'],
        'eps': eps,
        'resolution': resolution,
        'interaction_label': interaction_label,
        'coarse_label': coarse_label,
        'composite_lifetime': float(interaction_summary['composite_lifetime']),
        'binding_persistence_score': binding_persistence,
        'coarse_persistence_sigma4': coarse_summary['coarse_persistence_sigma4'],
        'mean_registration_strength': float(np.mean(registration_mean_series)) if registration_mean_series else 0.0,
        'max_registration_strength': float(np.max(registration_max_series)) if registration_max_series else 0.0,
        'registration_duty_cycle': registration_duty_cycle,
        'effective_interaction_energy_fraction': float(np.mean(effective_fraction_series)) if effective_fraction_series else 0.0,
        'encounter_dwell_mean': encounter_dwell_mean,
        'motif_survival_time': motif_survival_time,
        'corridor_dwell_time': corridor_dwell_time,
        'transport_channel_locking': transport_locking,
        'local_transport_anisotropy_index': local_transport_anisotropy_index,
        'trajectory_curvature_proxy': curvature_proxy,
        'energy_dispersion_shift': energy_dispersion_shift,
        'new_ordering_class': 0,
        'composite_lifetime_delta': 0.0,
        'binding_persistence_delta': 0.0,
        'coarse_persistence_delta': 0.0,
        'corridor_dwell_delta': 0.0,
        'gate_met': 0,
        'promoted_followup': 0,
        'output_label': 'no registration effect',
        'constraint_max': float(np.max(constraint_norms)) if constraint_norms else 0.0,
        'sector_leakage': float(np.max(sector_leakages)) if sector_leakages else 0.0,
        'notes': rep['notes'],
    }

    plot_paths: list[Path] = []
    trace_path = WORK_PLOT_DIR / f"stage21_{run['run_id']}_pair_distance_registration.png"
    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.2))
    axes[0].plot(times, pair_distance_series, color='tab:blue')
    axes[0].set_title('Mean pair distance')
    axes[0].set_xlabel('time')
    axes[0].grid(alpha=0.25)
    axes[1].plot(times, overlap_series, color='tab:orange')
    axes[1].set_title('Mean overlap')
    axes[1].set_xlabel('time')
    axes[1].grid(alpha=0.25)
    axes[2].plot(times, registration_mean_series, color='tab:green')
    axes[2].set_title('Mean registration strength')
    axes[2].set_xlabel('time')
    axes[2].grid(alpha=0.25)
    fig.suptitle(f"Stage 21: {run['run_id']}")
    fig.savefig(trace_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(trace_path)

    field_path = WORK_PLOT_DIR / f"stage21_{run['run_id']}_field_snapshot.png"
    create_field_snapshot(field_path, run['run_id'], data.midpoints, np.sum(packet_qs, axis=0), boundary_type, int(defaults['slice_axis']))
    plot_paths.append(field_path)

    env = coarse[4.0]['last_env']
    if env is not None:
        env_path = WORK_PLOT_DIR / f"stage21_{run['run_id']}_envelope_snapshot.png"
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
        'registration_mean_series': registration_mean_series,
        'effective_fraction_series': effective_fraction_series,
        'motif_strength_series': motif_strength_series,
    }
    return result, plot_paths


def classify_base_rows(rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, Any] | None]:
    promoted: dict[str, Any] | None = None
    best_score = None
    base_map = {(row['representative_id'], row['registration_class'], row['condition_id']): row for row in rows}
    for row in rows:
        if is_null_condition(str(row['condition_id'])):
            continue
        null_row = base_map.get((row['representative_id'], row['registration_class'], 'null'))
        if null_row is None:
            continue
        row['composite_lifetime_delta'] = float(row['composite_lifetime']) - float(null_row['composite_lifetime'])
        row['binding_persistence_delta'] = float(row['binding_persistence_score']) - float(null_row['binding_persistence_score'])
        row['coarse_persistence_delta'] = float(row['coarse_persistence_sigma4']) - float(null_row['coarse_persistence_sigma4'])
        row['corridor_dwell_delta'] = float(row['corridor_dwell_time']) - float(null_row['corridor_dwell_time'])
        row['new_ordering_class'] = int(
            interaction_rank(str(row['interaction_label'])) > interaction_rank(str(null_row['interaction_label']))
            or coarse_rank(str(row['coarse_label'])) > coarse_rank(str(null_row['coarse_label']))
        )
        persistence_gain = (
            row['composite_lifetime_delta'] >= COMPOSITE_GAIN_MIN
            or row['binding_persistence_delta'] >= BINDING_GAIN_MIN
            or row['coarse_persistence_delta'] >= COARSE_GAIN_MIN
            or row['corridor_dwell_delta'] >= CORRIDOR_GAIN_MIN
        )
        gate_hit = (
            persistence_gain
            and float(row['registration_duty_cycle']) > 0.0
            and float(row['effective_interaction_energy_fraction']) > 0.0
        )
        row['gate_met'] = int(gate_hit)
        if gate_hit:
            row['output_label'] = 'transient registration response'
            score = (
                float(row['composite_lifetime_delta'])
                + float(row['binding_persistence_delta'])
                + float(row['coarse_persistence_delta'])
                + float(row['corridor_dwell_delta'])
                + float(row['registration_duty_cycle'])
            )
            if best_score is None or score > best_score:
                best_score = score
                promoted = {
                    'representative_id': row['representative_id'],
                    'registration_class': row['registration_class'],
                    'condition_id': row['condition_id'],
                }
    if promoted is not None:
        for row in rows:
            if row['representative_id'] == promoted['representative_id'] and row['registration_class'] == promoted['registration_class'] and row['condition_id'] == promoted['condition_id']:
                row['promoted_followup'] = 1
    return rows, promoted


def create_summary_plots(rows: list[dict[str, Any]]) -> list[Path]:
    plot_paths: list[Path] = []
    classes = sorted({row['registration_class'] for row in rows}, key=class_order)
    reps = sorted({row['representative_id'] for row in rows})
    conds = sorted({row['condition_id'] for row in rows}, key=condition_order)
    row_map = {(row['registration_class'], row['representative_id'], row['condition_id']): row for row in rows}

    comparison = WORK_PLOT_DIR / 'stage21_registration_persistence_comparison.png'
    fig, axes = plt.subplots(len(classes), 2, figsize=(12.5, 4.0 * len(classes)))
    if len(classes) == 1:
        axes = np.asarray([axes])
    x = np.arange(len(reps))
    width = 0.22
    for class_idx, reg_class in enumerate(classes):
        for cond_idx, cond in enumerate(conds):
            lifetimes = [float(row_map[(reg_class, rep, cond)]['composite_lifetime']) for rep in reps]
            coarse_pers = [float(row_map[(reg_class, rep, cond)]['coarse_persistence_sigma4']) for rep in reps]
            axes[class_idx, 0].bar(x + (cond_idx - 1) * width, lifetimes, width=width, label=condition_label(cond))
            axes[class_idx, 1].bar(x + (cond_idx - 1) * width, coarse_pers, width=width, label=condition_label(cond))
        axes[class_idx, 0].set_title(f'{class_label(reg_class)}: composite lifetime')
        axes[class_idx, 1].set_title(f'{class_label(reg_class)}: coarse persistence')
        for ax in axes[class_idx]:
            ax.set_xticks(x)
            ax.set_xticklabels(reps, rotation=20, ha='right')
            ax.grid(alpha=0.25, axis='y')
            ax.legend(fontsize=8)
    fig.savefig(comparison, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(comparison)

    activity = WORK_PLOT_DIR / 'stage21_registration_activity_summary.png'
    fig, axes = plt.subplots(len(classes), 2, figsize=(12.5, 4.0 * len(classes)))
    if len(classes) == 1:
        axes = np.asarray([axes])
    for class_idx, reg_class in enumerate(classes):
        for cond_idx, cond in enumerate(conds):
            duty = [float(row_map[(reg_class, rep, cond)]['registration_duty_cycle']) for rep in reps]
            eff = [float(row_map[(reg_class, rep, cond)]['effective_interaction_energy_fraction']) for rep in reps]
            axes[class_idx, 0].bar(x + (cond_idx - 1) * width, duty, width=width, label=condition_label(cond))
            axes[class_idx, 1].bar(x + (cond_idx - 1) * width, eff, width=width, label=condition_label(cond))
        axes[class_idx, 0].set_title(f'{class_label(reg_class)}: registration duty cycle')
        axes[class_idx, 1].set_title(f'{class_label(reg_class)}: interaction energy fraction')
        for ax in axes[class_idx]:
            ax.set_xticks(x)
            ax.set_xticklabels(reps, rotation=20, ha='right')
            ax.grid(alpha=0.25, axis='y')
            ax.legend(fontsize=8)
    fig.savefig(activity, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(activity)

    response = WORK_PLOT_DIR / 'stage21_registration_response_matrix.png'
    combos = [(reg_class, rep) for reg_class in classes for rep in reps]
    label_map = {
        'no registration effect': 0,
        'transient registration response': 1,
        'scale-fragile registration-assisted persistence': 2,
        'scale-stable registration-assisted persistence': 3,
    }
    data = np.zeros((len(combos), len(conds)), dtype=float)
    ylabels: list[str] = []
    for i, (reg_class, rep) in enumerate(combos):
        ylabels.append(f'{class_label(reg_class)} / {rep}')
        for j, cond in enumerate(conds):
            data[i, j] = label_map.get(str(row_map[(reg_class, rep, cond)]['output_label']), 0)
    fig, ax = plt.subplots(figsize=(8.2, 7.0))
    im = ax.imshow(data, cmap='viridis', vmin=0, vmax=3)
    ax.set_xticks(range(len(conds)))
    ax.set_xticklabels([condition_label(cond) for cond in conds])
    ax.set_yticks(range(len(combos)))
    ax.set_yticklabels(ylabels)
    ax.set_title('Stage 21 registration response matrix')
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    fig.savefig(response, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(response)
    return plot_paths


def write_note(path: Path, json_rel: str, csv_rel: str, rows: list[dict[str, Any]], label_counts: dict[str, int], promoted: dict[str, Any] | None) -> None:
    lines = [
        '# Stage 21 Mutual Registration / Constraint-Coupled Interaction v1',
        '',
        f'Timestamped summary: `{json_rel}`',
        f'Timestamped run table: `{csv_rel}`',
        '',
        f'Output-label counts: {label_counts}',
        '',
        'This Stage 21 pilot tests whether bounded registration-gated pairwise interaction can create persistence without introducing substrate memory.',
        '',
    ]
    classes = sorted({row['registration_class'] for row in rows}, key=class_order)
    for reg_class in classes:
        lines.extend([
            f'## {reg_class}',
            '',
        ])
        for row in sorted([r for r in rows if r['registration_class'] == reg_class], key=lambda item: (item['representative_id'], condition_order(str(item['condition_id'])))):
            lines.extend([
                f"- `{row['representative_id']}` / `{row['condition_id']}`",
                f"  - composite lifetime delta: `{row['composite_lifetime_delta']:.4f}`",
                f"  - binding persistence delta: `{row['binding_persistence_delta']:.4f}`",
                f"  - coarse persistence delta: `{row['coarse_persistence_delta']:.4f}`",
                f"  - corridor dwell delta: `{row['corridor_dwell_delta']:.4f}`",
                f"  - mean registration strength: `{row['mean_registration_strength']:.4f}`",
                f"  - max registration strength: `{row['max_registration_strength']:.4f}`",
                f"  - registration duty cycle: `{row['registration_duty_cycle']:.4f}`",
                f"  - effective interaction energy fraction: `{row['effective_interaction_energy_fraction']:.4f}`",
                f"  - encounter dwell mean: `{row['encounter_dwell_mean']:.4f}`",
                f"  - motif survival time: `{row['motif_survival_time']:.4f}`",
                f"  - corridor dwell time: `{row['corridor_dwell_time']:.4f}`",
                f"  - transport locking: `{row['transport_channel_locking']}`",
                f"  - gate met: `{row['gate_met']}`",
                f"  - label: `{row['output_label']}`",
            ])
            lines.append('')
    lines.extend([
        f'Promoted follow-up: `{promoted}`' if promoted else 'Promoted follow-up: `none`',
        '',
        'Interpretation boundary:',
        '- any positive effect is read only as bounded registration-assisted persistence',
        '- no substrate memory, geometry, or force-law claim is made here',
        '',
    ])
    path.write_text('\n'.join(lines), encoding='utf-8')


def main() -> None:
    args = parse_args()
    defaults = load_stage10_defaults()
    runsheet = load_runsheet(args.runsheet)
    base = runsheet['base_seed_reference']
    rep_lookup = lookup_by_id(runsheet['representatives'], 'representative_id')
    runs = selected_runs(runsheet['runs'], args.run_ids)
    if not runs:
        raise SystemExit('no Stage 21 runs selected')

    rows: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    for run in runs:
        rep = rep_lookup[run['representative_id']]
        payload, run_plots = simulate_run(run, rep, runsheet, defaults, base)
        rows.append(payload['summary'])
        plot_paths.extend(run_plots)

    rows, promoted = classify_base_rows(rows)
    label_counts = dict(Counter(row['output_label'] for row in rows))
    plot_paths.extend(create_summary_plots(rows))

    result = {
        'stage': runsheet.get('stage', 'Stage 21'),
        'description': runsheet['description'],
        'coupling_law': runsheet.get('coupling_law', {}),
        'parameter_policy': runsheet.get('parameter_policy', {}),
        'base_seed_reference': base,
        'base_runs': [{'run_id': row['run_id'], 'summary': row} for row in rows],
        'output_label_counts': label_counts,
        'promoted_followup': promoted,
        'conclusion': 'Stage 21 tests whether selective interaction eligibility can create persistence without substrate memory on the projected transverse architecture.',
    }

    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug=str(runsheet.get('experiment_slug', 'stage21_mutual_registration')),
        result=result,
        csv_rows=rows,
        csv_fieldnames=CSV_FIELDS,
        plot_paths=plot_paths,
    )
    note_path = ATLAS_NOTES / str(runsheet.get('note_name', 'Stage_21_Mutual_Registration_Constraint_Coupled_Interaction_v1.md'))
    write_note(note_path, str(json_path.relative_to(REPO_ROOT)), str(csv_path.relative_to(REPO_ROOT)), rows, label_counts, promoted)

    print(f'wrote {json_path}')
    print(f'wrote {csv_path}')
    for item in stamped_plots:
        print(item)
    print(f'output_label_counts={label_counts}')
    print(f'promoted_followup={promoted}')


if __name__ == '__main__':
    main()
