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
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt, save_atlas_payload
from stage9b_common import (
    append_log,
    build_dk2d_complex,
    first_order_dt,
    gaussian_profile,
    load_stage9b_config,
    pack_positions,
    periodic_displacement,
)

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage23_1_dk_collision_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_23_1_DK_Collision_Geometry_Scan_v1.md'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage23_1_dk_collision'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'run_id',
    'geometry_id',
    'phase_id',
    'resolution',
    'graph_type',
    'packet_count',
    'initial_peak_count',
    'collision_label',
    'persistence_label',
    'composite_lifetime',
    'binding_persistence',
    'corridor_dwell',
    'coarse_basin_persistence',
    'minimum_separation',
    'final_mean_separation',
    'post_collision_separation_trend',
    'encounter_dwell_time',
    'deflection_angle_proxy',
    'reflection_fraction',
    'grade_transfer_amplitude',
    'omega0_weight_initial',
    'omega1_weight_initial',
    'omega2_weight_initial',
    'omega0_weight_final',
    'omega1_weight_final',
    'omega2_weight_final',
    'localized_circulation_proxy',
    'new_collision_class',
    'gate_met',
    'promoted_followup',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage 23.1 Dirac-Kaehler collision geometry scan.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    parser.add_argument('--append-log', action='store_true')
    parser.add_argument('--note-path', type=Path, default=None)
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


def packet_state(
    positions: np.ndarray,
    block_sizes: tuple[int, ...],
    center: np.ndarray,
    sigma: float,
    amplitude: float,
    phase_offset: float,
    kick_vector: np.ndarray,
    kick_cycles: float,
    active_grades: tuple[int, ...] = (0, 1),
) -> np.ndarray:
    parts: list[np.ndarray] = []
    start = 0
    for grade_idx, size in enumerate(block_sizes):
        coords = positions[start:start + size]
        profile, delta = gaussian_profile(coords, center, sigma)
        phase_arg = 2.0 * math.pi * kick_cycles * (delta @ kick_vector)
        if grade_idx in active_grades:
            part = profile * np.exp(1j * (phase_arg + phase_offset))
        else:
            part = np.zeros(size, dtype=complex)
        parts.append(part.astype(complex))
        start += size
    psi = np.concatenate(parts)
    psi /= max(np.linalg.norm(psi), 1.0e-12)
    return float(amplitude) * psi


def crank_nicolson_states(operator: sp.csr_matrix, psi0: np.ndarray, dt: float, steps: int) -> list[np.ndarray]:
    psi = np.asarray(psi0, dtype=complex).copy()
    ident = sp.identity(operator.shape[0], dtype=complex, format='csr')
    lhs = (ident - 0.5j * dt * operator).tocsc()
    rhs_op = (ident + 0.5j * dt * operator).tocsr()
    solve = spla.factorized(lhs)
    states = [psi.copy()]
    for _ in range(steps):
        psi = solve(rhs_op @ psi)
        states.append(np.asarray(psi, dtype=complex))
    return states


def grade_weights(vec: np.ndarray, block_sizes: tuple[int, ...]) -> list[float]:
    weights: list[float] = []
    start = 0
    total = float(np.sum(np.abs(vec) ** 2)) or 1.0
    for size in block_sizes:
        part = vec[start:start + size]
        weights.append(float(np.sum(np.abs(part) ** 2) / total))
        start += size
    return weights


def weighted_center(positions: np.ndarray, weights: np.ndarray) -> np.ndarray:
    total = float(np.sum(weights))
    if total <= 0.0:
        return np.zeros(2, dtype=float)
    center = np.zeros(2, dtype=float)
    for axis in range(2):
        angles = 2.0 * math.pi * positions[:, axis]
        s = float(np.sum(weights * np.sin(angles)))
        c = float(np.sum(weights * np.cos(angles)))
        center[axis] = ((math.atan2(s, c) / (2.0 * math.pi)) + 1.0) % 1.0
    return center


def field_grid_2d(positions: np.ndarray, vec: np.ndarray, n_side: int) -> np.ndarray:
    grid = np.zeros((n_side, n_side), dtype=float)
    coords = np.floor(np.asarray(positions) * n_side).astype(int) % n_side
    weights = np.abs(vec) ** 2
    for idx, coord in enumerate(coords):
        grid[coord[0], coord[1]] += float(weights[idx])
    return grid


def local_peak_candidates(grid: np.ndarray) -> list[tuple[float, tuple[int, int]]]:
    n_side = grid.shape[0]
    peaks: list[tuple[float, tuple[int, int]]] = []
    for i in range(n_side):
        for j in range(n_side):
            value = float(grid[i, j])
            if value <= 0.0:
                continue
            is_peak = True
            for di in (-1, 0, 1):
                for dj in (-1, 0, 1):
                    if di == 0 and dj == 0:
                        continue
                    ni = (i + di) % n_side
                    nj = (j + dj) % n_side
                    if grid[ni, nj] > value:
                        is_peak = False
                        break
                if not is_peak:
                    break
            if is_peak:
                peaks.append((value, (i, j)))
    if not peaks:
        flat_idx = int(np.argmax(grid))
        peaks.append((float(grid.flat[flat_idx]), np.unravel_index(flat_idx, grid.shape)))
    peaks.sort(key=lambda item: item[0], reverse=True)
    return peaks


def top_peak_positions(grid: np.ndarray, count: int) -> tuple[list[np.ndarray], int]:
    peaks = local_peak_candidates(grid)
    raw_count = len(peaks)
    taken = set()
    positions: list[np.ndarray] = []
    n_side = grid.shape[0]
    for _, (i, j) in peaks:
        if len(positions) >= count:
            break
        taken.add((i, j))
        positions.append(np.array([(i + 0.5) / n_side, (j + 0.5) / n_side], dtype=float))
    if len(positions) < count:
        flat_order = np.argsort(grid.ravel())[::-1]
        for flat_idx in flat_order:
            ij = np.unravel_index(int(flat_idx), grid.shape)
            if ij in taken:
                continue
            positions.append(np.array([(ij[0] + 0.5) / n_side, (ij[1] + 0.5) / n_side], dtype=float))
            taken.add(ij)
            if len(positions) >= count:
                break
    while len(positions) < count:
        positions.append(positions[-1].copy())
    return positions[:count], raw_count


def pairing_cost(points_a: list[np.ndarray], points_b: list[np.ndarray]) -> float:
    total = 0.0
    for a, b in zip(points_a, points_b):
        total += float(np.linalg.norm(periodic_displacement(np.asarray([b]), np.asarray(a))[0]))
    return total


def ordered_peaks(grid: np.ndarray, count: int, anchors: list[np.ndarray] | None, previous: list[np.ndarray] | None) -> tuple[list[np.ndarray], int]:
    peaks, raw_count = top_peak_positions(grid, count)
    if anchors is not None:
        best = min(permutations(peaks, len(peaks)), key=lambda perm: pairing_cost(list(anchors), list(perm)))
        return [np.asarray(p, dtype=float) for p in best], raw_count
    if previous is not None:
        best = min(permutations(peaks, len(peaks)), key=lambda perm: pairing_cost(list(previous), list(perm)))
        return [np.asarray(p, dtype=float) for p in best], raw_count
    peaks.sort(key=lambda p: (p[0], p[1]))
    return peaks, raw_count


def mean_pair_distance(centers: list[np.ndarray]) -> float:
    values = []
    for i, j in combinations(range(len(centers)), 2):
        delta = periodic_displacement(np.asarray([centers[j]], dtype=float), np.asarray(centers[i], dtype=float))[0]
        values.append(float(np.linalg.norm(delta)))
    return float(np.mean(values)) if values else 0.0


def min_pair_distance(centers: list[np.ndarray]) -> float:
    values = []
    for i, j in combinations(range(len(centers)), 2):
        delta = periodic_displacement(np.asarray([centers[j]], dtype=float), np.asarray(centers[i], dtype=float))[0]
        values.append(float(np.linalg.norm(delta)))
    return float(np.min(values)) if values else 0.0


def mean_pair_separation_series(center_histories: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    mean_vals: list[float] = []
    min_vals: list[float] = []
    steps = len(center_histories[0]) if center_histories else 0
    for idx in range(steps):
        centers = [hist[idx] for hist in center_histories]
        mean_vals.append(mean_pair_distance(centers))
        min_vals.append(min_pair_distance(centers))
    return np.asarray(mean_vals, dtype=float), np.asarray(min_vals, dtype=float)


def longest_constant_run(values: list[int]) -> float:
    if not values:
        return 0.0
    best = 1
    current = 1
    for prev, nxt in zip(values, values[1:]):
        if prev == nxt:
            current += 1
        else:
            best = max(best, current)
            current = 1
    best = max(best, current)
    return float(best / len(values))


def basin_id(center: np.ndarray) -> int:
    x = 0 if center[0] < 0.5 else 1
    y = 0 if center[1] < 0.5 else 1
    return 2 * y + x


def angle_between(a: np.ndarray, b: np.ndarray) -> float:
    a_norm = float(np.linalg.norm(a))
    b_norm = float(np.linalg.norm(b))
    if a_norm <= 1.0e-12 or b_norm <= 1.0e-12:
        return 0.0
    cosine = float(np.clip(np.dot(a, b) / (a_norm * b_norm), -1.0, 1.0))
    return float(math.acos(cosine))


def classify_collision(
    geometry_id: str,
    initial_peak_count: int,
    reflection_fraction: float,
    deflection_angle: float,
    grade_transfer_amplitude: float,
    bound_tail_fraction: float,
    encounter_dwell_time: float,
    t_final: float,
) -> str:
    if bound_tail_fraction >= 0.6 and grade_transfer_amplitude >= 0.18:
        if initial_peak_count >= 2:
            return 'oscillatory grade-trapped'
        return 'unresolved / mixed'
    if bound_tail_fraction >= 0.6:
        if initial_peak_count >= 2:
            return 'metastable composite'
        return 'unresolved / mixed'
    if reflection_fraction >= 0.75:
        return 'reflective / exclusion-like'
    if deflection_angle >= 0.35:
        return 'deflective / glancing'
    if grade_transfer_amplitude >= 0.12:
        return 'pass-through with grade exchange'
    if geometry_id == 'symmetric_triad_collision' and encounter_dwell_time >= 0.2 * t_final:
        return 'unresolved / mixed'
    return 'pass-through dispersive'


def classify_persistence(collision_label: str, composite_lifetime: float, encounter_dwell_time: float, t_final: float) -> str:
    if collision_label in {'metastable composite', 'oscillatory grade-trapped'} and composite_lifetime >= 0.2 * t_final:
        return 'candidate bounded collision regime'
    if encounter_dwell_time >= 0.15 * t_final:
        return 'weak persistence gain'
    return 'no persistence gain'


def plot_trajectories(path: Path, run_id: str, center_histories: list[np.ndarray], geometry_id: str) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 5.0))
    colors = ['tab:blue', 'tab:orange', 'tab:green']
    for idx, centers in enumerate(center_histories):
        ax.plot(centers[:, 0], centers[:, 1], color=colors[idx % len(colors)], label=f'peak {idx + 1}')
        ax.scatter([centers[0, 0]], [centers[0, 1]], color=colors[idx % len(colors)], marker='o', s=20)
        ax.scatter([centers[-1, 0]], [centers[-1, 1]], color=colors[idx % len(colors)], marker='x', s=30)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'{geometry_id}: {run_id}')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_trace(path: Path, run_id: str, times: list[float], values: np.ndarray, ylabel: str, title: str) -> None:
    fig, ax = plt.subplots(figsize=(5.8, 4.4))
    ax.plot(times, values, color='tab:blue')
    ax.set_xlabel('time')
    ax.set_ylabel(ylabel)
    ax.set_title(f'{title}: {run_id}')
    ax.grid(alpha=0.25)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_grades(path: Path, run_id: str, times: list[float], grade_hist: np.ndarray) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    ax.stackplot(times, grade_hist.T, labels=[r'$\Omega^0$', r'$\Omega^1$', r'$\Omega^2$'], alpha=0.85)
    ax.set_xlabel('time')
    ax.set_ylabel('grade weight fraction')
    ax.set_title(f'Grade weights: {run_id}')
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8, loc='upper right')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_collision_matrix(path: Path, rows: list[dict[str, Any]]) -> None:
    geometry_order = [
        'head_on_symmetric_pair',
        'counter_propagating_corridor_pair',
        'offset_glancing_collision',
        'tight_clustered_pair',
        'symmetric_triad_collision',
    ]
    phase_order = ['in_phase', 'out_of_phase']
    label_to_value = {
        'pass-through dispersive': 0,
        'pass-through with grade exchange': 1,
        'reflective / exclusion-like': 2,
        'deflective / glancing': 3,
        'metastable composite': 4,
        'oscillatory grade-trapped': 5,
        'unresolved / mixed': 6,
    }
    grid = np.full((len(geometry_order), len(phase_order)), -1, dtype=int)
    labels = {}
    for row in rows:
        i = geometry_order.index(row['geometry_id'])
        j = phase_order.index(row['phase_id'])
        grid[i, j] = label_to_value[row['collision_label']]
        labels[(i, j)] = row['collision_label']

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    im = ax.imshow(grid, cmap='viridis', aspect='auto', vmin=0, vmax=max(label_to_value.values()))
    ax.set_xticks(range(len(phase_order)), phase_order)
    ax.set_yticks(range(len(geometry_order)), geometry_order)
    ax.set_title('Collision class matrix')
    for (i, j), label in labels.items():
        ax.text(j, i, label.split(' / ')[0], ha='center', va='center', fontsize=8, color='white')
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_persistence_panel(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [row['run_id'].replace('S23_1_', '') for row in rows]
    composite = [float(row['composite_lifetime']) for row in rows]
    binding = [float(row['binding_persistence']) for row in rows]
    fig, axes = plt.subplots(2, 1, figsize=(10.0, 6.8), sharex=True)
    axes[0].bar(range(len(rows)), composite, color='tab:blue')
    axes[0].set_ylabel('composite lifetime')
    axes[0].set_title('Persistence comparison panel')
    axes[0].grid(axis='y', alpha=0.25)
    axes[1].bar(range(len(rows)), binding, color='tab:orange')
    axes[1].set_ylabel('binding persistence')
    axes[1].set_xticks(range(len(rows)), labels, rotation=45, ha='right')
    axes[1].grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_grade_transfer_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [row['run_id'].replace('S23_1_', '') for row in rows]
    transfer = [float(row['grade_transfer_amplitude']) for row in rows]
    fig, ax = plt.subplots(figsize=(10.0, 4.6))
    ax.bar(range(len(rows)), transfer, color='tab:green')
    ax.set_ylabel('grade transfer amplitude')
    ax.set_xticks(range(len(rows)), labels, rotation=45, ha='right')
    ax.set_title('Grade-transfer summary')
    ax.grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def simulate_run(run: dict[str, Any], geometry: dict[str, Any], fixed: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], list[Path]]:
    resolution = int(run['resolution'])
    epsilon = float(fixed['epsilon'])
    bandwidth = float(fixed['bandwidth'])
    amplitude = float(fixed['amplitude'])
    t_final = float(fixed['t_final'])
    kick_cycles = float(fixed['kick_cycles'])
    dt_scale = float(fixed['dt_scale'])

    complex_data = build_dk2d_complex(n_side=resolution, epsilon=epsilon)
    positions = pack_positions(complex_data.points, complex_data.edge_midpoints, complex_data.face_centers)
    block_sizes = complex_data.block_sizes
    phase_offset = float(run['phase_offset_rad'])

    packet_states: list[np.ndarray] = []
    kick_vectors = [np.asarray(vec, dtype=float) for vec in geometry['kick_vectors']]
    for idx, center in enumerate(geometry['packet_centers']):
        packet_states.append(
            packet_state(
                positions=positions,
                block_sizes=block_sizes,
                center=np.asarray(center, dtype=float),
                sigma=bandwidth,
                amplitude=amplitude,
                phase_offset=0.0 if idx == 0 else phase_offset,
                kick_vector=kick_vectors[idx],
                kick_cycles=kick_cycles,
            )
        )

    total_psi0 = np.sum(packet_states, axis=0)
    dt = first_order_dt(complex_data.dirac_kahler, dt_scale)
    steps = max(2, int(math.ceil(t_final / dt)))
    states = crank_nicolson_states(complex_data.dirac_kahler, total_psi0, dt, steps)
    times = [idx * dt for idx in range(len(states))]

    grade_hist = np.asarray([grade_weights(state, block_sizes) for state in states], dtype=float)
    total_centers = np.asarray([weighted_center(positions, np.abs(state) ** 2) for state in states], dtype=float)
    basin_ids = [basin_id(center) for center in total_centers]
    coarse_basin_persistence = longest_constant_run(basin_ids)

    peak_tracks: list[list[np.ndarray]] = [[] for _ in range(int(geometry['packet_count']))]
    raw_peak_counts: list[int] = []
    previous: list[np.ndarray] | None = None
    anchors = [np.asarray(center, dtype=float) for center in geometry['packet_centers']]
    for idx, state in enumerate(states):
        grid = field_grid_2d(positions, state, resolution)
        current, raw_count = ordered_peaks(grid, int(geometry['packet_count']), anchors if idx == 0 else None, previous)
        raw_peak_counts.append(raw_count)
        for track, point in zip(peak_tracks, current):
            track.append(point.copy())
        previous = [point.copy() for point in current]
    center_histories = [np.asarray(track, dtype=float) for track in peak_tracks]

    sample_dt = float(np.median(np.diff(times))) if len(times) >= 2 else dt
    mean_pair_distances, min_pair_distances = mean_pair_separation_series(center_histories)
    close_threshold = 2.0 * bandwidth
    composite_lifetime = float(np.sum(mean_pair_distances <= close_threshold) * sample_dt)
    binding_persistence = float(np.mean(mean_pair_distances <= close_threshold)) if mean_pair_distances.size else 0.0
    encounter_dwell_time = composite_lifetime
    min_separation = float(np.min(min_pair_distances)) if min_pair_distances.size else 0.0
    final_mean_separation = float(mean_pair_distances[-1]) if mean_pair_distances.size else 0.0
    post_collision_trend = final_mean_separation - min_separation
    tail = mean_pair_distances[max(0, int(0.75 * len(mean_pair_distances))):]
    bound_tail_fraction = float(np.mean(tail <= close_threshold)) if tail.size else 0.0
    initial_peak_count = int(raw_peak_counts[0]) if raw_peak_counts else 0

    reflection_votes = 0
    deflection_angles: list[float] = []
    circulation_terms: list[float] = []
    for centers, kick in zip(center_histories, kick_vectors):
        net_disp = periodic_displacement(centers[-1][None, :], centers[0])[0]
        kick_unit = kick / max(np.linalg.norm(kick), 1.0e-12)
        along = float(np.dot(net_disp, kick_unit))
        if along < 0.0:
            reflection_votes += 1
        deflection_angles.append(angle_between(net_disp, kick_unit))
        if len(centers) >= 3:
            velocity_0 = periodic_displacement(centers[1][None, :], centers[0])[0]
            velocity_1 = periodic_displacement(centers[-1][None, :], centers[-2])[0]
            circulation_terms.append(float(velocity_0[0] * velocity_1[1] - velocity_0[1] * velocity_1[0]))
    reflection_fraction = float(reflection_votes / max(len(center_histories), 1))
    deflection_angle_proxy = float(np.mean(deflection_angles)) if deflection_angles else 0.0
    localized_circulation_proxy = float(np.mean(np.abs(circulation_terms))) if circulation_terms else 0.0

    corridor_dwell = 0.0
    if geometry['geometry_id'] == 'counter_propagating_corridor_pair':
        flags = []
        for step in range(len(times)):
            flags.append(float(all(abs(hist[step, 1] - 0.5) <= 0.08 for hist in center_histories)))
        corridor_dwell = float(np.sum(flags) * sample_dt)

    initial_grade = grade_hist[0]
    final_grade = grade_hist[-1]
    grade_transfer_amplitude = float(np.max(np.linalg.norm(grade_hist - initial_grade[None, :], axis=1)))

    collision_label = classify_collision(
        geometry_id=geometry['geometry_id'],
        initial_peak_count=initial_peak_count,
        reflection_fraction=reflection_fraction,
        deflection_angle=deflection_angle_proxy,
        grade_transfer_amplitude=grade_transfer_amplitude,
        bound_tail_fraction=bound_tail_fraction,
        encounter_dwell_time=encounter_dwell_time,
        t_final=float(times[-1]),
    )
    persistence_label = classify_persistence(
        collision_label=collision_label,
        composite_lifetime=composite_lifetime,
        encounter_dwell_time=encounter_dwell_time,
        t_final=float(times[-1]),
    )
    new_collision_class = int(collision_label not in {'pass-through dispersive', 'unresolved / mixed'})
    gate_met = int(persistence_label == 'candidate bounded collision regime' and new_collision_class == 1)

    row = {
        'run_id': run['run_id'],
        'geometry_id': geometry['geometry_id'],
        'phase_id': run['phase_id'],
        'resolution': resolution,
        'graph_type': run['graph_type'],
        'packet_count': int(geometry['packet_count']),
        'initial_peak_count': initial_peak_count,
        'collision_label': collision_label,
        'persistence_label': persistence_label,
        'composite_lifetime': composite_lifetime,
        'binding_persistence': binding_persistence,
        'corridor_dwell': corridor_dwell,
        'coarse_basin_persistence': coarse_basin_persistence,
        'minimum_separation': min_separation,
        'final_mean_separation': final_mean_separation,
        'post_collision_separation_trend': post_collision_trend,
        'encounter_dwell_time': encounter_dwell_time,
        'deflection_angle_proxy': deflection_angle_proxy,
        'reflection_fraction': reflection_fraction,
        'grade_transfer_amplitude': grade_transfer_amplitude,
        'omega0_weight_initial': float(initial_grade[0]),
        'omega1_weight_initial': float(initial_grade[1]),
        'omega2_weight_initial': float(initial_grade[2]),
        'omega0_weight_final': float(final_grade[0]),
        'omega1_weight_final': float(final_grade[1]),
        'omega2_weight_final': float(final_grade[2]),
        'localized_circulation_proxy': localized_circulation_proxy,
        'new_collision_class': new_collision_class,
        'gate_met': gate_met,
        'promoted_followup': gate_met,
        'notes': geometry['notes'],
    }

    detail = {
        'config': {
            'resolution': resolution,
            'epsilon': epsilon,
            'bandwidth': bandwidth,
            'amplitude': amplitude,
            't_final': t_final,
            'dt': dt,
            'steps': steps,
            'kick_cycles': kick_cycles,
        },
        'geometry': geometry,
        'run': run,
        'metrics': row,
        'times': times,
        'mean_pair_distances': mean_pair_distances.tolist(),
        'minimum_pair_distances': min_pair_distances.tolist(),
        'peak_count_series': raw_peak_counts,
        'total_grade_weights': grade_hist.tolist(),
        'total_centers': total_centers.tolist(),
        'peak_tracks': [hist.tolist() for hist in center_histories],
    }

    plot_paths: list[Path] = []
    traj_path = WORK_PLOT_DIR / f"stage23_1_{run['run_id']}_collision_trajectory.png"
    sep_path = WORK_PLOT_DIR / f"stage23_1_{run['run_id']}_separation_trace.png"
    grade_path = WORK_PLOT_DIR / f"stage23_1_{run['run_id']}_grade_weight_trace.png"
    plot_trajectories(traj_path, run['run_id'], center_histories, geometry['geometry_id'])
    plot_trace(sep_path, run['run_id'], times, mean_pair_distances, 'mean peak separation', 'Separation trace')
    plot_grades(grade_path, run['run_id'], times, grade_hist)
    plot_paths.extend([traj_path, sep_path, grade_path])

    return row, detail, plot_paths


def write_note(json_path: Path, csv_path: Path, rows: list[dict[str, Any]], stamped_plots: list[str], note_path: Path | None = None) -> None:
    collision_counts = Counter(row['collision_label'] for row in rows)
    persistence_counts = Counter(row['persistence_label'] for row in rows)
    gated = sum(int(row['gate_met']) for row in rows)
    lines = [
        '# Stage 23.1 DK Collision Geometry Scan v1',
        '',
        f"Timestamped JSON: `{json_path.relative_to(REPO_ROOT)}`",
        f"Timestamped CSV: `{csv_path.relative_to(REPO_ROOT)}`",
        '',
        f'Total runs: {len(rows)}',
        f'Candidate bounded collision regimes: {gated}',
        '',
        '## Collision labels',
    ]
    for label, count in sorted(collision_counts.items()):
        lines.append(f'- {label}: {count}')
    lines.extend([
        '',
        '## Persistence labels',
    ])
    for label, count in sorted(persistence_counts.items()):
        lines.append(f'- {label}: {count}')
    lines.extend([
        '',
        '## Plots',
    ])
    for plot in stamped_plots:
        lines.append(f'- `{plot}`')
    target = note_path or NOTE_PATH
    target.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main() -> None:
    args = parse_args()
    runsheet = load_runsheet(args.runsheet)
    runs = selected_runs(runsheet['runs'], args.run_ids)
    geometry_by_id = lookup_by_id(runsheet['geometries'], 'geometry_id')
    fixed = runsheet['fixed_parameters']

    rows: list[dict[str, Any]] = []
    details: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    for run in runs:
        geometry = geometry_by_id[run['geometry_id']]
        row, detail, run_plots = simulate_run(run, geometry, fixed)
        rows.append(row)
        details.append(detail)
        plot_paths.extend(run_plots)

    rows.sort(key=lambda item: item['run_id'])

    matrix_path = WORK_PLOT_DIR / 'stage23_1_collision_class_matrix.png'
    persistence_path = WORK_PLOT_DIR / 'stage23_1_persistence_comparison_panel.png'
    grade_summary_path = WORK_PLOT_DIR / 'stage23_1_grade_transfer_summary.png'
    plot_collision_matrix(matrix_path, rows)
    plot_persistence_panel(persistence_path, rows)
    plot_grade_transfer_summary(grade_summary_path, rows)
    plot_paths.extend([matrix_path, persistence_path, grade_summary_path])

    result = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'fixed_parameters': fixed,
        'runs': details,
        'summary': {
            'collision_labels': dict(Counter(row['collision_label'] for row in rows)),
            'persistence_labels': dict(Counter(row['persistence_label'] for row in rows)),
            'candidate_bounded_collision_regimes': sum(int(row['gate_met']) for row in rows),
        },
    }
    json_path, csv_path, stamped_plots, _ = save_atlas_payload(
        runsheet['experiment_slug'],
        result,
        rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots, note_path=args.note_path)

    if args.append_log:
        observation = 'the Stage 23.1 DK collision scan was executed with the minimal 10-run geometry-phase matrix'
        conclusion = 'the collision geometry baseline is ready for persistence-first classification in the Dirac-Kaehler sector'
        append_log(
            'Stage 23.1 DK Collision Geometry Scan',
            f"runs={len(rows)}, resolution={fixed['resolution']}, graph={fixed['graph_type']}, bandwidth={fixed['bandwidth']}",
            json_path,
            [REPO_ROOT / rel for rel in stamped_plots],
            observation,
            conclusion,
        )


if __name__ == '__main__':
    main()
