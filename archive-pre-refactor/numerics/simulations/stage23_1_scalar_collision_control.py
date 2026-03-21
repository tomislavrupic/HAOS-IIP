#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from DK_stage6_common import build_dk2d_complex
from stage10_common import ATLAS_NOTES, REPO_ROOT, plt, save_atlas_payload
from stage23_1_dk_collision_geometry_scan import (
    CSV_FIELDS,
    WORK_PLOT_DIR,
    angle_between,
    basin_id,
    field_grid_2d,
    load_runsheet,
    longest_constant_run,
    lookup_by_id,
    mean_pair_separation_series,
    ordered_peaks,
    plot_collision_matrix,
    plot_persistence_panel,
    plot_trace,
    plot_trajectories,
    selected_runs,
)
from stage9b_common import first_order_dt, gaussian_profile, pack_positions, periodic_displacement

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage23_1_scalar_collision_control_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_23_1_Scalar_Collision_Control_v1.md'


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the matched scalar bosonic control matrix for Stage 23.1.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    return parser.parse_args()


def scalar_packet_state(
    points: np.ndarray,
    center: np.ndarray,
    sigma: float,
    amplitude: float,
    phase_offset: float,
    kick_vector: np.ndarray,
    kick_cycles: float,
) -> tuple[np.ndarray, np.ndarray]:
    profile, delta = gaussian_profile(points, center, sigma)
    phase_arg = 2.0 * math.pi * kick_cycles * (delta @ kick_vector)
    phase = np.exp(1j * (phase_arg + phase_offset))
    q0 = profile * phase
    q0 /= max(np.linalg.norm(q0), 1.0e-12)
    slope = (delta @ kick_vector) / max(sigma * sigma, 1.0e-12)
    v0 = slope * profile * phase
    v0 /= max(np.linalg.norm(v0), 1.0e-12)
    return float(amplitude) * q0, float(amplitude) * v0


def leapfrog_states(operator, q0: np.ndarray, v0: np.ndarray, dt: float, steps: int) -> list[np.ndarray]:
    q = np.asarray(q0, dtype=complex).copy()
    v = np.asarray(v0, dtype=complex).copy()
    states = [q.copy()]
    for _ in range(steps):
        acc = -(operator @ q)
        v_half = v + 0.5 * dt * acc
        q = q + dt * v_half
        acc_new = -(operator @ q)
        v = v_half + 0.5 * dt * acc_new
        states.append(np.asarray(q, dtype=complex))
    return states


def classify_collision(
    reflection_fraction: float,
    deflection_angle: float,
    bound_tail_fraction: float,
    encounter_dwell_time: float,
    t_final: float,
) -> str:
    if bound_tail_fraction >= 0.6:
        return 'metastable composite'
    if reflection_fraction >= 0.75:
        return 'reflective / exclusion-like'
    if deflection_angle >= 0.35:
        return 'deflective / glancing'
    if encounter_dwell_time >= 0.15 * t_final:
        return 'unresolved / mixed'
    return 'pass-through dispersive'


def classify_persistence(collision_label: str, composite_lifetime: float, encounter_dwell_time: float, t_final: float) -> str:
    if collision_label == 'metastable composite' and composite_lifetime >= 0.2 * t_final:
        return 'candidate bounded collision regime'
    if encounter_dwell_time >= 0.15 * t_final:
        return 'weak persistence gain'
    return 'no persistence gain'


def plot_scalar_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [row['run_id'].replace('S23_1_scalar_', '') for row in rows]
    deflection = [float(row['deflection_angle_proxy']) for row in rows]
    min_sep = [float(row['minimum_separation']) for row in rows]
    fig, axes = plt.subplots(2, 1, figsize=(10.0, 6.8), sharex=True)
    axes[0].bar(range(len(rows)), deflection, color='tab:green')
    axes[0].set_ylabel('deflection angle')
    axes[0].set_title('Scalar control collision summary')
    axes[0].grid(axis='y', alpha=0.25)
    axes[1].bar(range(len(rows)), min_sep, color='tab:purple')
    axes[1].set_ylabel('min separation')
    axes[1].set_xticks(range(len(rows)), labels, rotation=45, ha='right')
    axes[1].grid(axis='y', alpha=0.25)
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
    points = np.asarray(complex_data.points, dtype=float)
    operator = (complex_data.delta1 @ complex_data.d0).tocsr()

    q_parts = []
    v_parts = []
    kick_vectors = [np.asarray(vec, dtype=float) for vec in geometry['kick_vectors']]
    phase_offset = float(run['phase_offset_rad'])
    for idx, center in enumerate(geometry['packet_centers']):
        q0, v0 = scalar_packet_state(
            points=points,
            center=np.asarray(center, dtype=float),
            sigma=bandwidth,
            amplitude=amplitude,
            phase_offset=0.0 if idx == 0 else phase_offset,
            kick_vector=kick_vectors[idx],
            kick_cycles=kick_cycles,
        )
        q_parts.append(q0)
        v_parts.append(v0)

    q_total = np.sum(q_parts, axis=0)
    v_total = np.sum(v_parts, axis=0)
    dt = first_order_dt(operator @ operator, dt_scale)
    steps = max(2, int(math.ceil(t_final / dt)))
    states = leapfrog_states(operator, q_total, v_total, dt, steps)
    times = [idx * dt for idx in range(len(states))]

    peak_tracks: list[list[np.ndarray]] = [[] for _ in range(int(geometry['packet_count']))]
    raw_peak_counts: list[int] = []
    previous: list[np.ndarray] | None = None
    anchors = [np.asarray(center, dtype=float) for center in geometry['packet_centers']]
    for idx, state in enumerate(states):
        grid = field_grid_2d(points, state, resolution)
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

    total_centers = []
    basin_ids = []
    for state in states:
        weights = np.abs(state) ** 2
        total = float(np.sum(weights))
        if total <= 0.0:
            center = np.zeros(2, dtype=float)
        else:
            center = np.zeros(2, dtype=float)
            for axis in range(2):
                angles = 2.0 * math.pi * points[:, axis]
                s = float(np.sum(weights * np.sin(angles)))
                c = float(np.sum(weights * np.cos(angles)))
                center[axis] = ((math.atan2(s, c) / (2.0 * math.pi)) + 1.0) % 1.0
        total_centers.append(center)
        basin_ids.append(basin_id(center))
    coarse_basin_persistence = longest_constant_run(basin_ids)

    reflection_votes = 0
    deflection_angles: list[float] = []
    for centers, kick in zip(center_histories, kick_vectors):
        net_disp = periodic_displacement(centers[-1][None, :], centers[0])[0]
        kick_unit = kick / max(np.linalg.norm(kick), 1.0e-12)
        along = float(np.dot(net_disp, kick_unit))
        if along < 0.0:
            reflection_votes += 1
        deflection_angles.append(angle_between(net_disp, kick_unit))
    reflection_fraction = float(reflection_votes / max(len(center_histories), 1))
    deflection_angle_proxy = float(np.mean(deflection_angles)) if deflection_angles else 0.0

    corridor_dwell = 0.0
    if geometry['geometry_id'] == 'counter_propagating_corridor_pair':
        flags = []
        for step in range(len(times)):
            flags.append(float(all(abs(hist[step, 1] - 0.5) <= 0.08 for hist in center_histories)))
        corridor_dwell = float(np.sum(flags) * sample_dt)

    collision_label = classify_collision(
        reflection_fraction=reflection_fraction,
        deflection_angle=deflection_angle_proxy,
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
    gate_met = int(persistence_label == 'candidate bounded collision regime')

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
        'grade_transfer_amplitude': 0.0,
        'omega0_weight_initial': 1.0,
        'omega1_weight_initial': 0.0,
        'omega2_weight_initial': 0.0,
        'omega0_weight_final': 1.0,
        'omega1_weight_final': 0.0,
        'omega2_weight_final': 0.0,
        'localized_circulation_proxy': 0.0,
        'new_collision_class': int(collision_label not in {'pass-through dispersive', 'unresolved / mixed'}),
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
        'total_centers': [center.tolist() for center in total_centers],
        'peak_tracks': [hist.tolist() for hist in center_histories],
    }

    plot_paths: list[Path] = []
    traj_path = WORK_PLOT_DIR / f"stage23_1_scalar_{run['run_id']}_collision_trajectory.png"
    sep_path = WORK_PLOT_DIR / f"stage23_1_scalar_{run['run_id']}_separation_trace.png"
    plot_trajectories(traj_path, run['run_id'], center_histories, geometry['geometry_id'])
    plot_trace(sep_path, run['run_id'], times, mean_pair_distances, 'mean peak separation', 'Scalar separation trace')
    plot_paths.extend([traj_path, sep_path])
    return row, detail, plot_paths


def write_note(json_path: Path, csv_path: Path, rows: list[dict[str, Any]], stamped_plots: list[str]) -> None:
    collision_counts = Counter(row['collision_label'] for row in rows)
    persistence_counts = Counter(row['persistence_label'] for row in rows)
    gated = sum(int(row['gate_met']) for row in rows)
    lines = [
        '# Stage 23.1 Scalar Collision Control v1',
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
    lines.extend(['', '## Persistence labels'])
    for label, count in sorted(persistence_counts.items()):
        lines.append(f'- {label}: {count}')
    lines.extend(['', '## Plots'])
    for plot in stamped_plots:
        lines.append(f'- `{plot}`')
    NOTE_PATH.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main() -> None:
    args = parse_args()
    runsheet = load_runsheet(args.runsheet)
    runs = selected_runs(runsheet['runs'], args.run_ids)
    geometry_by_id = lookup_by_id(runsheet['geometries'], 'geometry_id')
    fixed = runsheet['fixed_parameters']

    rows = []
    details = []
    plot_paths = []
    for run in runs:
        geometry = geometry_by_id[run['geometry_id']]
        row, detail, run_plots = simulate_run(run, geometry, fixed)
        rows.append(row)
        details.append(detail)
        plot_paths.extend(run_plots)

    rows.sort(key=lambda item: item['run_id'])

    matrix_path = WORK_PLOT_DIR / 'stage23_1_scalar_collision_class_matrix.png'
    persistence_path = WORK_PLOT_DIR / 'stage23_1_scalar_persistence_comparison_panel.png'
    summary_path = WORK_PLOT_DIR / 'stage23_1_scalar_collision_summary.png'
    plot_collision_matrix(matrix_path, rows)
    plot_persistence_panel(persistence_path, rows)
    plot_scalar_summary(summary_path, rows)
    plot_paths.extend([matrix_path, persistence_path, summary_path])

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
    write_note(json_path, csv_path, rows, stamped_plots)


if __name__ == '__main__':
    main()
