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
import scipy.sparse.linalg as spla

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt, save_atlas_payload
from stage23_1_dk_collision_geometry_scan import (
    CSV_FIELDS as BASE_FIELDS,
    WORK_PLOT_DIR,
    basin_id,
    classify_collision,
    classify_persistence,
    field_grid_2d,
    load_runsheet,
    longest_constant_run,
    mean_pair_separation_series,
    ordered_peaks,
    plot_trace,
    plot_trajectories,
    plot_grades,
    selected_runs,
)
from stage9b_common import build_dk2d_complex, first_order_dt, gaussian_profile, pack_positions, periodic_displacement

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage23_3_dk_clustered_texture_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_23_3_DK_Clustered_Texture_Stabilization_v1.md'

CSV_FIELDS = BASE_FIELDS + [
    'probe_class',
    'strength_label',
    'alpha',
    'beta',
    'gamma',
    'phase_difference_mean',
    'phase_difference_std',
    'bounded_oscillatory_separation',
    'grade_exchange_coherence',
    'promotion_candidate',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage 23.3 DK clustered texture stabilization probe.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    return parser.parse_args()


def packet_state(
    positions: np.ndarray,
    block_sizes: tuple[int, ...],
    center: np.ndarray,
    sigma: float,
    amplitude: float,
    kick_vector: np.ndarray,
    kick_cycles: float,
) -> np.ndarray:
    parts: list[np.ndarray] = []
    start = 0
    for grade_idx, size in enumerate(block_sizes):
        coords = positions[start:start + size]
        profile, delta = gaussian_profile(coords, center, sigma)
        phase_arg = 2.0 * math.pi * kick_cycles * (delta @ kick_vector)
        phase = np.exp(1j * phase_arg)
        if grade_idx in (0, 1):
            part = profile * phase
        else:
            part = np.zeros(size, dtype=complex)
        parts.append(np.asarray(part, dtype=complex))
        start += size
    psi = np.concatenate(parts)
    psi /= max(np.linalg.norm(psi), 1.0e-12)
    return float(amplitude) * psi


def grade_weights(vec: np.ndarray, block_sizes: tuple[int, ...]) -> list[float]:
    weights: list[float] = []
    start = 0
    total = float(np.sum(np.abs(vec) ** 2)) or 1.0
    for size in block_sizes:
        part = vec[start:start + size]
        weights.append(float(np.sum(np.abs(part) ** 2) / total))
        start += size
    return weights


def local_grade_mix(vec: np.ndarray, block_sizes: tuple[int, ...]) -> np.ndarray:
    n0, n1, n2 = block_sizes
    q0 = vec[:n0]
    q1 = vec[n0:n0 + n1]
    out = np.zeros_like(vec)
    m = min(n0, n1)
    out[:m] += q1[:m]
    out[n0:n0 + m] += q0[:m]
    return out


def graded_energy_density(vec: np.ndarray, positions: np.ndarray, resolution: int) -> np.ndarray:
    grid = field_grid_2d(positions, vec, resolution)
    mean = float(np.mean(grid))
    centered = grid - mean
    max_abs = float(np.max(np.abs(centered)))
    if max_abs <= 1.0e-12:
        return np.zeros_like(centered)
    return centered / max_abs


def edge_scale_from_grid(grid: np.ndarray, edge_midpoints: np.ndarray, gamma: float, clamp: float = 0.05) -> np.ndarray:
    n_side = grid.shape[0]
    scale = np.ones(len(edge_midpoints), dtype=float)
    for idx, point in enumerate(edge_midpoints):
        ii = int(math.floor(point[0] * n_side)) % n_side
        jj = int(math.floor(point[1] * n_side)) % n_side
        factor = float(np.clip(gamma * grid[ii, jj], -clamp, clamp))
        scale[idx] = 1.0 + factor
    return scale


def pair_phase_difference(state: np.ndarray, packet_states: list[np.ndarray]) -> float:
    refs = []
    for packet in packet_states:
        overlap = np.vdot(packet, state)
        refs.append(float(np.angle(overlap)))
    if len(refs) < 2:
        return 0.0
    diff = refs[1] - refs[0]
    return float((diff + math.pi) % (2.0 * math.pi) - math.pi)


def bounded_oscillation(series: np.ndarray, close_threshold: float) -> bool:
    if len(series) < 6:
        return False
    tail = series[len(series) // 2:]
    if tail.size == 0:
        return False
    below = tail <= close_threshold
    if np.mean(below) < 0.5:
        return False
    span = float(np.max(tail) - np.min(tail))
    return 0.02 <= span <= 0.18


def evolve_probe(
    base_operator: sp.csr_matrix,
    positions: np.ndarray,
    edge_midpoints: np.ndarray,
    block_sizes: tuple[int, ...],
    packet_states: list[np.ndarray],
    alpha: float,
    beta: float,
    gamma: float,
    dt: float,
    steps: int,
    bandwidth: float,
    resolution: int,
) -> tuple[list[np.ndarray], list[float], list[float]]:
    state = np.sum(packet_states, axis=0)
    states = [state.copy()]
    phase_diffs = [pair_phase_difference(state, packet_states)]
    coherence_flags = [0.0]

    n0, n1, n2 = block_sizes
    edge_slice = slice(n0, n0 + n1)
    base_edge = np.asarray(base_operator[edge_slice, edge_slice].diagonal()).real
    edge_operator = base_operator.copy().tolil()

    for _ in range(steps):
        probe_force = np.zeros_like(state)
        if alpha > 0.0:
            probe_force += alpha * local_grade_mix(state, block_sizes)

        phase_diff = pair_phase_difference(state, packet_states)
        phase_diffs.append(phase_diff)

        density_grid = graded_energy_density(state, positions, resolution)
        if gamma > 0.0:
            edge_scale = edge_scale_from_grid(density_grid, edge_midpoints, gamma)
            edge_operator = base_operator.copy().tolil()
            edge_operator[edge_slice, edge_slice] = sp.diags(base_edge * edge_scale)
            active_operator = edge_operator.tocsr()
        else:
            active_operator = base_operator

        acc = -(active_operator @ state) + probe_force
        if beta > 0.0:
            total_weights = np.abs(state) ** 2
            center_weights = np.sum(total_weights)
            center_close = center_weights > 0.0 and abs(phase_diff) < math.pi
            if center_close:
                probe_phase = -1j * beta * phase_diff * state
                acc += probe_phase
                coherence_flags.append(1.0)
            else:
                coherence_flags.append(0.0)
        else:
            coherence_flags.append(0.0)

        state = state + dt * acc
        norm = np.linalg.norm(state)
        if norm > 1.0e-12:
            state = state / norm * np.linalg.norm(states[0])
        states.append(state.copy())

    return states, phase_diffs[: len(states)], coherence_flags[: len(states)]


def simulate_run(run: dict[str, Any], common: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], list[Path]]:
    resolution = int(common['resolution'])
    epsilon = 0.2
    bandwidth = 0.08
    amplitude = float(common['base_amplitude_scale']) * 0.5
    t_final = 0.6
    kick_cycles = 1.0

    complex_data = build_dk2d_complex(n_side=resolution, epsilon=epsilon)
    positions = pack_positions(complex_data.points, complex_data.edge_midpoints, complex_data.face_centers)
    block_sizes = complex_data.block_sizes

    centers = [np.array([0.46, 0.50], dtype=float), np.array([0.54, 0.50], dtype=float)]
    kicks = [np.array([0.4, 0.0], dtype=float), np.array([-0.4, 0.0], dtype=float)]
    packet_states = [
        packet_state(positions, block_sizes, centers[idx], bandwidth, amplitude, kicks[idx], kick_cycles)
        for idx in range(2)
    ]

    dt = first_order_dt(complex_data.dirac_kahler, 0.35)
    steps = max(2, int(math.ceil(t_final / dt)))
    states, phase_diffs, coherence_flags = evolve_probe(
        base_operator=complex_data.dirac_kahler,
        positions=positions,
        edge_midpoints=complex_data.edge_midpoints,
        block_sizes=block_sizes,
        packet_states=packet_states,
        alpha=float(run['alpha']),
        beta=float(run['beta']),
        gamma=float(run['gamma']),
        dt=dt,
        steps=steps,
        bandwidth=bandwidth,
        resolution=resolution,
    )
    times = [idx * dt for idx in range(len(states))]

    grade_hist = np.asarray([grade_weights(state, block_sizes) for state in states], dtype=float)

    peak_tracks: list[list[np.ndarray]] = [[] for _ in range(2)]
    raw_peak_counts: list[int] = []
    previous: list[np.ndarray] | None = None
    anchors = [center.copy() for center in centers]
    for idx, state in enumerate(states):
        grid = field_grid_2d(positions, state, resolution)
        current, raw_count = ordered_peaks(grid, 2, anchors if idx == 0 else None, previous)
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

    total_centers = []
    basin_ids = []
    for state in states:
        weights = np.abs(state) ** 2
        total = float(np.sum(weights))
        center = np.zeros(2, dtype=float)
        if total > 0.0:
            for axis in range(2):
                angles = 2.0 * math.pi * positions[:, axis]
                s = float(np.sum(weights * np.sin(angles)))
                c = float(np.sum(weights * np.cos(angles)))
                center[axis] = ((math.atan2(s, c) / (2.0 * math.pi)) + 1.0) % 1.0
        total_centers.append(center)
        basin_ids.append(basin_id(center))
    coarse_basin_persistence = longest_constant_run(basin_ids)

    reflection_fraction = 0.0
    deflection_angle_proxy = 0.0
    grade_transfer_amplitude = float(np.max(np.linalg.norm(grade_hist - grade_hist[0][None, :], axis=1)))
    collision_label = classify_collision(
        geometry_id='tight_clustered_pair',
        initial_peak_count=initial_peak_count,
        reflection_fraction=reflection_fraction,
        deflection_angle=deflection_angle_proxy,
        grade_transfer_amplitude=grade_transfer_amplitude,
        bound_tail_fraction=bound_tail_fraction,
        encounter_dwell_time=encounter_dwell_time,
        t_final=float(times[-1]),
    )
    persistence_label = classify_persistence(collision_label, composite_lifetime, encounter_dwell_time, float(times[-1]))

    phase_array = np.asarray(phase_diffs, dtype=float)
    bounded_sep = int(bounded_oscillation(mean_pair_distances, close_threshold))
    grade_coherence = float(np.mean(coherence_flags)) if coherence_flags else 0.0
    promotion_candidate = int(
        composite_lifetime >= 0.2 * float(times[-1])
        and bounded_sep == 1
        and grade_coherence >= 0.25
    )

    row = {
        'run_id': run['run_id'],
        'geometry_id': 'tight_clustered_pair',
        'phase_id': 'in_phase',
        'resolution': resolution,
        'graph_type': common['graph_type'],
        'packet_count': 2,
        'initial_peak_count': initial_peak_count,
        'collision_label': collision_label,
        'persistence_label': persistence_label,
        'composite_lifetime': composite_lifetime,
        'binding_persistence': binding_persistence,
        'corridor_dwell': 0.0,
        'coarse_basin_persistence': coarse_basin_persistence,
        'minimum_separation': min_separation,
        'final_mean_separation': final_mean_separation,
        'post_collision_separation_trend': post_collision_trend,
        'encounter_dwell_time': encounter_dwell_time,
        'deflection_angle_proxy': deflection_angle_proxy,
        'reflection_fraction': reflection_fraction,
        'grade_transfer_amplitude': grade_transfer_amplitude,
        'omega0_weight_initial': float(grade_hist[0, 0]),
        'omega1_weight_initial': float(grade_hist[0, 1]),
        'omega2_weight_initial': float(grade_hist[0, 2]),
        'omega0_weight_final': float(grade_hist[-1, 0]),
        'omega1_weight_final': float(grade_hist[-1, 1]),
        'omega2_weight_final': float(grade_hist[-1, 2]),
        'localized_circulation_proxy': 0.0,
        'new_collision_class': int(collision_label not in {'pass-through dispersive', 'unresolved / mixed'}),
        'gate_met': promotion_candidate,
        'promoted_followup': promotion_candidate,
        'notes': run['notes'],
        'probe_class': run['probe_class'],
        'strength_label': run['strength_label'],
        'alpha': float(run['alpha']),
        'beta': float(run['beta']),
        'gamma': float(run['gamma']),
        'phase_difference_mean': float(np.mean(np.abs(phase_array))) if phase_array.size else 0.0,
        'phase_difference_std': float(np.std(phase_array)) if phase_array.size else 0.0,
        'bounded_oscillatory_separation': bounded_sep,
        'grade_exchange_coherence': grade_coherence,
        'promotion_candidate': promotion_candidate,
    }

    detail = {
        'run': run,
        'metrics': row,
        'times': times,
        'mean_pair_distances': mean_pair_distances.tolist(),
        'minimum_pair_distances': min_pair_distances.tolist(),
        'phase_differences': phase_diffs,
        'grade_weights': grade_hist.tolist(),
        'peak_tracks': [hist.tolist() for hist in center_histories],
    }

    plot_paths: list[Path] = []
    traj_path = WORK_PLOT_DIR / f"stage23_3_{run['run_id']}_collision_trajectory.png"
    sep_path = WORK_PLOT_DIR / f"stage23_3_{run['run_id']}_separation_trace.png"
    grade_path = WORK_PLOT_DIR / f"stage23_3_{run['run_id']}_grade_weight_trace.png"
    phase_path = WORK_PLOT_DIR / f"stage23_3_{run['run_id']}_phase_difference_trace.png"
    plot_trajectories(traj_path, run['run_id'], center_histories, 'tight_clustered_pair')
    plot_trace(sep_path, run['run_id'], times, mean_pair_distances, 'mean peak separation', 'Stage 23.3 separation trace')
    plot_grades(grade_path, run['run_id'], times, grade_hist)
    plot_trace(phase_path, run['run_id'], times, phase_array, 'phase difference', 'Stage 23.3 phase-difference trace')
    plot_paths.extend([traj_path, sep_path, grade_path, phase_path])
    return row, detail, plot_paths


def plot_persistence_panel(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [row['run_id'].replace('S23_3_', '') for row in rows]
    lifetimes = [float(row['composite_lifetime']) for row in rows]
    coherence = [float(row['grade_exchange_coherence']) for row in rows]
    fig, axes = plt.subplots(2, 1, figsize=(11.0, 7.0), sharex=True)
    axes[0].bar(range(len(rows)), lifetimes, color='tab:blue')
    axes[0].set_ylabel('composite lifetime')
    axes[0].set_title('Stage 23.3 persistence-gain comparison')
    axes[0].grid(axis='y', alpha=0.25)
    axes[1].bar(range(len(rows)), coherence, color='tab:green')
    axes[1].set_ylabel('grade coherence')
    axes[1].set_xticks(range(len(rows)), labels, rotation=45, ha='right')
    axes[1].grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_stability_table(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [row['run_id'].replace('S23_3_', '') for row in rows]
    values = [
        2 if int(row['promotion_candidate']) else (1 if int(row['bounded_oscillatory_separation']) else 0)
        for row in rows
    ]
    fig, ax = plt.subplots(figsize=(11.0, 3.8))
    ax.bar(range(len(rows)), values, color='tab:orange')
    ax.set_ylabel('stability score')
    ax.set_xticks(range(len(rows)), labels, rotation=45, ha='right')
    ax.set_title('Stage 23.3 stability classification')
    ax.grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_coherence_map(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [row['run_id'].replace('S23_3_', '') for row in rows]
    values = [float(row['phase_difference_std']) for row in rows]
    fig, ax = plt.subplots(figsize=(11.0, 3.8))
    ax.bar(range(len(rows)), values, color='tab:purple')
    ax.set_ylabel('phase difference std')
    ax.set_xticks(range(len(rows)), labels, rotation=45, ha='right')
    ax.set_title('Stage 23.3 grade-coherence evolution map')
    ax.grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(json_path: Path, csv_path: Path, rows: list[dict[str, Any]], stamped_plots: list[str]) -> None:
    collision_counts = Counter(row['collision_label'] for row in rows)
    persistence_counts = Counter(row['persistence_label'] for row in rows)
    promoted = sum(int(row['promotion_candidate']) for row in rows)
    lines = [
        '# Stage 23.3 DK Clustered Texture Stabilization v1',
        '',
        f"Timestamped JSON: `{json_path.relative_to(REPO_ROOT)}`",
        f"Timestamped CSV: `{csv_path.relative_to(REPO_ROOT)}`",
        '',
        f'Total runs: {len(rows)}',
        f'Promotion candidates: {promoted}',
        '',
        '## Collision labels',
    ]
    for label, count in sorted(collision_counts.items()):
        lines.append(f'- {label}: {count}')
    lines.extend(['', '## Persistence labels'])
    for label, count in sorted(persistence_counts.items()):
        lines.append(f'- {label}: {count}')
    lines.extend(['', '## Plots'])
    for plot in stamped_plots:
        lines.append(f'- `{plot}`')
    NOTE_PATH.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main() -> None:
    args = parse_args()
    runsheet = json.loads(args.runsheet.read_text(encoding='utf-8'))
    runs = selected_runs(runsheet['runs'], args.run_ids)
    common = runsheet['common_fields']

    rows: list[dict[str, Any]] = []
    details: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    for run in runs:
        row, detail, run_plots = simulate_run(run, common)
        rows.append(row)
        details.append(detail)
        plot_paths.extend(run_plots)

    rows.sort(key=lambda item: item['run_id'])
    details.sort(key=lambda item: item['run']['run_id'])

    persistence_path = WORK_PLOT_DIR / 'stage23_3_persistence_gain_panel.png'
    stability_path = WORK_PLOT_DIR / 'stage23_3_stability_classification_table.png'
    coherence_path = WORK_PLOT_DIR / 'stage23_3_grade_coherence_map.png'
    plot_persistence_panel(persistence_path, rows)
    plot_stability_table(stability_path, rows)
    plot_coherence_map(coherence_path, rows)
    plot_paths.extend([persistence_path, stability_path, coherence_path])

    result = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'common_fields': common,
        'runs': details,
        'summary': {
            'collision_labels': dict(Counter(row['collision_label'] for row in rows)),
            'persistence_labels': dict(Counter(row['persistence_label'] for row in rows)),
            'promotion_candidates': sum(int(row['promotion_candidate']) for row in rows),
        },
    }
    json_path, csv_path, stamped_plots, _ = save_atlas_payload(
        'stage23_3_dk_clustered_texture',
        result,
        rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots)


if __name__ == '__main__':
    main()
