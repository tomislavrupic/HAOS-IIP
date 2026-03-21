#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt, save_atlas_payload
from stage23_1_dk_collision_geometry_scan import (
    CSV_FIELDS as BASE_FIELDS,
    WORK_PLOT_DIR,
    field_grid_2d,
    load_runsheet,
    mean_pair_separation_series,
    ordered_peaks,
    plot_trace,
    plot_trajectories,
    plot_grades,
    selected_runs,
    weighted_center,
)
from stage9b_common import build_dk2d_complex, first_order_dt, gaussian_profile, pack_positions, periodic_displacement

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage23_6_dk_phase_sensitivity_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_23_6_DK_Phase_Sensitivity_Braid_Exchange_v1.md'

CSV_FIELDS = BASE_FIELDS + [
    'run_class',
    'variant_label',
    'phase_offset_fraction_of_pi',
    'local_phase_skew_fraction',
    'beta',
    'flow_concentration_index',
    'grade_exchange_coherence',
    'channel_count',
    'loop_count',
    'recirculation_score',
    'phase_alignment_metric',
    'topology_label',
    'braid_robust',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage 23.6 DK phase-sensitivity and braid-exchange stability scan.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    return parser.parse_args()


def packet_state_with_skew(
    positions: np.ndarray,
    block_sizes: tuple[int, ...],
    center: np.ndarray,
    sigma: float,
    amplitude: float,
    phase_offset: float,
    kick_vector: np.ndarray,
    kick_cycles: float,
    phase_skew_fraction: float,
) -> np.ndarray:
    parts: list[np.ndarray] = []
    start = 0
    kick_norm = max(float(np.linalg.norm(kick_vector)), 1.0e-12)
    skew_dir = kick_vector / kick_norm
    for grade_idx, size in enumerate(block_sizes):
        coords = positions[start:start + size]
        profile, delta = gaussian_profile(coords, center, sigma)
        base_phase = 2.0 * math.pi * kick_cycles * (delta @ kick_vector)
        skew_arg = math.pi * phase_skew_fraction * np.clip((delta @ skew_dir) / max(sigma, 1.0e-12), -1.0, 1.0)
        phase = np.exp(1j * (base_phase + phase_offset + skew_arg))
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


def pair_phase_difference(state: np.ndarray, packet_states: list[np.ndarray]) -> float:
    refs = []
    for packet in packet_states:
        overlap = np.vdot(packet, state)
        refs.append(float(np.angle(overlap)))
    if len(refs) < 2:
        return 0.0
    diff = refs[1] - refs[0]
    return float((diff + math.pi) % (2.0 * math.pi) - math.pi)


def local_grade_mix(vec: np.ndarray, block_sizes: tuple[int, ...]) -> np.ndarray:
    n0, n1, _ = block_sizes
    q0 = vec[:n0]
    q1 = vec[n0:n0 + n1]
    out = np.zeros_like(vec)
    m = min(n0, n1)
    out[:m] += q1[:m]
    out[n0:n0 + m] += q0[:m]
    return out


def edge_grid(edge_midpoints: np.ndarray, edge_values: np.ndarray, n_side: int) -> np.ndarray:
    grid = np.zeros((n_side, n_side), dtype=float)
    coords = np.floor(np.asarray(edge_midpoints) * n_side).astype(int) % n_side
    for idx, coord in enumerate(coords):
        grid[coord[0], coord[1]] += float(edge_values[idx])
    return grid


def connected_components(mask: np.ndarray) -> int:
    seen = np.zeros_like(mask, dtype=bool)
    n0, n1 = mask.shape
    count = 0
    for i in range(n0):
        for j in range(n1):
            if not mask[i, j] or seen[i, j]:
                continue
            count += 1
            stack = [(i, j)]
            seen[i, j] = True
            while stack:
                x, y = stack.pop()
                for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    nx = x + dx
                    ny = y + dy
                    if 0 <= nx < n0 and 0 <= ny < n1 and mask[nx, ny] and not seen[nx, ny]:
                        seen[nx, ny] = True
                        stack.append((nx, ny))
    return count


def grade_exchange_signal(grade_hist: np.ndarray) -> tuple[np.ndarray, float]:
    signal = np.linalg.norm(grade_hist - grade_hist[0][None, :], axis=1)
    peak = float(np.max(signal)) if signal.size else 0.0
    if peak <= 1.0e-12:
        return signal, 0.0
    coherence = float(np.mean(signal) / peak)
    return signal, max(0.0, min(1.0, coherence))


def phase_alignment_metric(phase_diffs: list[float]) -> float:
    if not phase_diffs:
        return 0.0
    mean_abs = float(np.mean(np.abs(np.asarray(phase_diffs, dtype=float))))
    return float(np.clip(1.0 - mean_abs / math.pi, 0.0, 1.0))


def recurrence_indicator(series: np.ndarray, close_threshold: float) -> float:
    if len(series) < 6:
        return 0.0
    tail = series[len(series) // 3:]
    if tail.size < 4:
        return 0.0
    deriv = np.diff(tail)
    signs = np.sign(deriv)
    signs = signs[signs != 0.0]
    if signs.size < 2:
        return 0.0
    turns = np.sum(signs[1:] * signs[:-1] < 0.0)
    closeness = float(np.mean(tail <= close_threshold))
    return float((turns / max(1, signs.size - 1)) * closeness)


def estimate_topology(avg_edge_grid: np.ndarray, braid_flag: bool, coherence: float, amplitude: float) -> tuple[str, float, int, int, float]:
    total = float(np.sum(avg_edge_grid))
    if total <= 1.0e-12:
        return 'dispersive_pass', 0.0, 0, 0, 0.0
    flat = np.sort(avg_edge_grid.ravel())[::-1]
    top_k = max(1, int(math.ceil(0.1 * flat.size)))
    concentration = float(np.sum(flat[:top_k]) / total)
    threshold = max(0.15 * float(np.max(avg_edge_grid)), float(np.mean(avg_edge_grid) + np.std(avg_edge_grid)))
    mask = avg_edge_grid >= threshold
    channels = connected_components(mask)
    loops = 0
    recirc = 0.0
    if braid_flag and coherence >= 0.45 and concentration >= 0.86:
        return 'braid_like_exchange', concentration, channels, loops, recirc
    if coherence >= 0.40 and amplitude >= 0.14:
        return 'transfer_smeared', concentration, channels, loops, recirc
    if coherence < 0.35 or concentration < 0.75:
        return 'dispersive_pass', concentration, channels, loops, recirc
    return 'mixed_unresolved', concentration, channels, loops, recirc


def evolve(base_operator, block_sizes: tuple[int, ...], state0: np.ndarray, dt: float, steps: int, beta: float) -> list[np.ndarray]:
    state = state0.copy()
    states = [state.copy()]
    base_norm = np.linalg.norm(state0)
    for _ in range(steps):
        acc = -(base_operator @ state)
        if beta > 0.0:
            acc += beta * local_grade_mix(state, block_sizes)
        state = state + dt * acc
        norm = np.linalg.norm(state)
        if not np.isfinite(norm) or norm <= 1.0e-12:
            raise RuntimeError('numerical instability indicator triggered')
        state = state / norm * base_norm
        states.append(state.copy())
    return states


def simulate_run(run: dict[str, Any], common: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], list[Path]]:
    resolution = int(common['resolution'])
    epsilon = 0.2
    sigma = float(common['mean_width'])
    amplitude = 0.5 * float(common['amplitude_scale'])
    t_final = 0.9
    kick_cycles = 1.0
    separation = 0.08

    complex_data = build_dk2d_complex(n_side=resolution, epsilon=epsilon)
    positions = pack_positions(complex_data.points, complex_data.edge_midpoints, complex_data.face_centers)
    block_sizes = complex_data.block_sizes
    n0, n1, _ = block_sizes

    centers = [
        np.array([0.50 - 0.5 * separation, 0.50], dtype=float),
        np.array([0.50 + 0.5 * separation, 0.50], dtype=float),
    ]
    kicks = [np.array([0.4, 0.0], dtype=float), np.array([-0.4, 0.0], dtype=float)]
    phase_offset = float(run['phase_offset_fraction_of_pi']) * math.pi
    skew = float(run['local_phase_skew_fraction'])
    beta = float(run['beta'])

    packet_states = [
        packet_state_with_skew(
            positions=positions,
            block_sizes=block_sizes,
            center=centers[idx],
            sigma=sigma,
            amplitude=amplitude,
            phase_offset=0.0 if idx == 0 else phase_offset,
            kick_vector=kicks[idx],
            kick_cycles=kick_cycles,
            phase_skew_fraction=0.0 if idx == 0 else skew,
        )
        for idx in range(2)
    ]

    psi0 = np.sum(packet_states, axis=0)
    dt = first_order_dt(complex_data.dirac_kahler, 0.35)
    steps = max(2, int(math.ceil(t_final / dt)))
    states = evolve(complex_data.dirac_kahler, block_sizes, psi0, dt, steps, beta)
    times = [idx * dt for idx in range(len(states))]

    grade_hist = np.asarray([grade_weights(state, block_sizes) for state in states], dtype=float)
    transfer_signal, coherence = grade_exchange_signal(grade_hist)
    grade_amplitude = float(np.max(transfer_signal)) if transfer_signal.size else 0.0

    peak_tracks: list[list[np.ndarray]] = [[] for _ in range(2)]
    previous: list[np.ndarray] | None = None
    anchors = [center.copy() for center in centers]
    braid_votes = []
    edge_grids = []
    phase_diffs = []
    raw_peak_counts: list[int] = []
    for idx, state in enumerate(states):
        grid = field_grid_2d(positions, state, resolution)
        current, raw_count = ordered_peaks(grid, 2, anchors if idx == 0 else None, previous)
        raw_peak_counts.append(raw_count)
        for track, point in zip(peak_tracks, current):
            track.append(point.copy())
        if len(current) == 2:
            braid_votes.append(float(current[0][0] > current[1][0]))
        previous = [point.copy() for point in current]
        edge_values = np.abs(state[n0:n0 + n1])
        edge_grids.append(edge_grid(complex_data.edge_midpoints, edge_values, resolution))
        phase_diffs.append(pair_phase_difference(state, packet_states))
    center_histories = [np.asarray(track, dtype=float) for track in peak_tracks]
    avg_edge_grid = np.mean(np.asarray(edge_grids, dtype=float), axis=0)

    total_centers = np.asarray([weighted_center(positions, np.abs(state) ** 2) for state in states], dtype=float)
    mean_pair_distances, min_pair_distances = mean_pair_separation_series(center_histories)
    close_threshold = 2.0 * sigma
    composite_lifetime = float(np.sum(mean_pair_distances <= close_threshold) * dt)
    binding_persistence = float(np.mean(mean_pair_distances <= close_threshold)) if mean_pair_distances.size else 0.0
    final_mean_separation = float(mean_pair_distances[-1]) if mean_pair_distances.size else 0.0
    min_separation = float(np.min(min_pair_distances)) if min_pair_distances.size else 0.0
    post_collision_trend = final_mean_separation - min_separation

    braid_flag = abs(float(np.mean(braid_votes)) - 0.5) > 0.25 if braid_votes else False
    topo, concentration, channels, loops, recirc = estimate_topology(avg_edge_grid, braid_flag, coherence, grade_amplitude)
    phase_align = phase_alignment_metric(phase_diffs)
    recurrence = recurrence_indicator(mean_pair_distances, close_threshold)

    row = {
        'run_id': run['run_id'],
        'geometry_id': 'tight_clustered_pair',
        'phase_id': run['variant_label'],
        'resolution': resolution,
        'graph_type': common['graph_type'],
        'packet_count': 2,
        'initial_peak_count': int(raw_peak_counts[0]) if raw_peak_counts else 0,
        'collision_label': topo,
        'persistence_label': 'weak persistence gain' if composite_lifetime >= 0.15 * float(times[-1]) else 'no persistence gain',
        'composite_lifetime': composite_lifetime,
        'binding_persistence': binding_persistence,
        'corridor_dwell': 0.0,
        'coarse_basin_persistence': 0.0,
        'minimum_separation': min_separation,
        'final_mean_separation': final_mean_separation,
        'post_collision_separation_trend': post_collision_trend,
        'encounter_dwell_time': composite_lifetime,
        'deflection_angle_proxy': 0.0,
        'reflection_fraction': 0.0,
        'grade_transfer_amplitude': grade_amplitude,
        'omega0_weight_initial': float(grade_hist[0, 0]),
        'omega1_weight_initial': float(grade_hist[0, 1]),
        'omega2_weight_initial': float(grade_hist[0, 2]),
        'omega0_weight_final': float(grade_hist[-1, 0]),
        'omega1_weight_final': float(grade_hist[-1, 1]),
        'omega2_weight_final': float(grade_hist[-1, 2]),
        'localized_circulation_proxy': recirc,
        'new_collision_class': int(topo == 'braid_like_exchange'),
        'gate_met': 0,
        'promoted_followup': 0,
        'notes': run['notes'],
        'run_class': run['run_class'],
        'variant_label': run['variant_label'],
        'phase_offset_fraction_of_pi': float(run['phase_offset_fraction_of_pi']),
        'local_phase_skew_fraction': skew,
        'beta': beta,
        'flow_concentration_index': concentration,
        'grade_exchange_coherence': coherence,
        'channel_count': channels,
        'loop_count': loops,
        'recirculation_score': recirc,
        'phase_alignment_metric': phase_align,
        'topology_label': topo,
        'braid_robust': int(topo == 'braid_like_exchange'),
    }

    detail = {
        'run': run,
        'metrics': row,
        'times': times,
        'mean_pair_distances': mean_pair_distances.tolist(),
        'minimum_pair_distances': min_pair_distances.tolist(),
        'phase_differences': phase_diffs,
        'grade_weights': grade_hist.tolist(),
        'grade_exchange_trace': transfer_signal.tolist(),
        'peak_tracks': [hist.tolist() for hist in center_histories],
        'flow_concentration_trace': [float(np.sum(np.sort(g.ravel())[::-1][:max(1, int(math.ceil(0.1 * g.size)))]) / max(np.sum(g), 1.0e-12)) for g in edge_grids],
        'edge_current_map': avg_edge_grid.tolist(),
        'recurrence_indicator': recurrence,
    }

    plot_paths: list[Path] = []
    traj_path = WORK_PLOT_DIR / f"stage23_6_{run['run_id']}_topology_trajectory.png"
    flow_path = WORK_PLOT_DIR / f"stage23_6_{run['run_id']}_flow_concentration_trace.png"
    grade_path = WORK_PLOT_DIR / f"stage23_6_{run['run_id']}_grade_exchange_trace.png"
    grades_path = WORK_PLOT_DIR / f"stage23_6_{run['run_id']}_grade_weights.png"
    plot_trajectories(traj_path, run['run_id'], center_histories, 'tight_clustered_pair')
    plot_trace(flow_path, run['run_id'], times, np.asarray(detail['flow_concentration_trace'], dtype=float), 'flow concentration', 'Stage 23.6 flow concentration trace')
    plot_trace(grade_path, run['run_id'], times, transfer_signal, 'grade exchange', 'Stage 23.6 grade-exchange trace')
    plot_grades(grades_path, run['run_id'], times, grade_hist)
    plot_paths.extend([traj_path, flow_path, grade_path, grades_path])
    return row, detail, plot_paths


def plot_topology_matrix(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [row['run_id'].replace('S23_6_', '') for row in rows]
    mapping = {
        'braid_like_exchange': 0,
        'transfer_smeared': 1,
        'mixed_unresolved': 2,
        'dispersive_pass': 3,
    }
    values = np.asarray([[mapping[str(row['topology_label'])]] for row in rows], dtype=float)
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    im = ax.imshow(values, aspect='auto', cmap='viridis', vmin=0.0, vmax=3.0)
    ax.set_yticks(range(len(labels)), labels)
    ax.set_xticks([0], ['topology'])
    ax.set_title('Stage 23.6 topology class matrix vs phase')
    for idx, row in enumerate(rows):
        ax.text(0, idx, str(row['topology_label']).replace('_', ' '), ha='center', va='center', color='white', fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_coherence_heatmap(path: Path, details: list[dict[str, Any]]) -> None:
    max_len = max(len(detail['grade_exchange_trace']) for detail in details)
    heat = np.zeros((len(details), max_len), dtype=float)
    for idx, detail in enumerate(details):
        signal = np.asarray(detail['grade_exchange_trace'], dtype=float)
        peak = float(np.max(signal)) if signal.size else 0.0
        if peak > 1.0e-12:
            heat[idx, :signal.size] = signal / peak
    labels = [detail['run']['run_id'].replace('S23_6_', '') for detail in details]
    fig, ax = plt.subplots(figsize=(10.4, 4.8))
    im = ax.imshow(heat, aspect='auto', cmap='magma', vmin=0.0, vmax=1.0)
    ax.set_yticks(range(len(labels)), labels)
    ax.set_xlabel('time index')
    ax.set_title('Stage 23.6 coherence heatmap')
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_braid_robustness(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [row['run_id'].replace('S23_6_', '') for row in rows]
    vals = [float(row['phase_alignment_metric']) for row in rows]
    colors = ['tab:blue' if str(row['topology_label']) == 'braid_like_exchange' else 'tab:orange' for row in rows]
    fig, ax = plt.subplots(figsize=(10.4, 4.8))
    ax.bar(range(len(rows)), vals, color=colors)
    ax.set_ylabel('phase alignment metric')
    ax.set_xticks(range(len(rows)), labels, rotation=45, ha='right')
    ax.set_title('Stage 23.6 braid-family robustness panel')
    ax.grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(json_path: Path, csv_path: Path, rows: list[dict[str, Any]], stamped_plots: list[str], braid_corridor: list[str]) -> None:
    topo_counts = Counter(str(row['topology_label']) for row in rows)
    lines = [
        '# Stage 23.6 DK Phase Sensitivity Braid Exchange v1',
        '',
        f"Timestamped JSON: `{json_path.relative_to(REPO_ROOT)}`",
        f"Timestamped CSV: `{csv_path.relative_to(REPO_ROOT)}`",
        '',
        f'Total runs: {len(rows)}',
        '',
        '## Topology labels',
    ]
    for label, count in sorted(topo_counts.items()):
        lines.append(f'- {label}: {count}')
    lines.extend([
        '',
        '## Braid corridor',
        f"- {braid_corridor if braid_corridor else 'none'}",
        '',
        '## Plots',
    ])
    for plot in stamped_plots:
        lines.append(f'- `{plot}`')
    NOTE_PATH.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main() -> None:
    args = parse_args()
    runsheet = load_runsheet(args.runsheet)
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

    corridor_order = [
        'S23_6_phase_baseline',
        'S23_6_phase_quarter',
        'S23_6_phase_half',
        'S23_6_phase_three_quarter',
        'S23_6_phase_antiphase',
    ]
    braid_corridor = [row['run_id'] for row in rows if row['run_id'] in corridor_order and row['topology_label'] == 'braid_like_exchange']

    matrix_path = WORK_PLOT_DIR / 'stage23_6_topology_class_matrix.png'
    heatmap_path = WORK_PLOT_DIR / 'stage23_6_coherence_heatmap.png'
    robustness_path = WORK_PLOT_DIR / 'stage23_6_braid_robustness_panel.png'
    plot_topology_matrix(matrix_path, rows)
    plot_coherence_heatmap(heatmap_path, details)
    plot_braid_robustness(robustness_path, rows)
    plot_paths.extend([matrix_path, heatmap_path, robustness_path])

    result = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'common_fields': common,
        'runs': details,
        'summary': {
            'topology_labels': dict(Counter(row['topology_label'] for row in rows)),
            'braid_corridor': braid_corridor,
        },
    }
    json_path, csv_path, stamped_plots, _ = save_atlas_payload(
        'stage23_6_dk_phase_sensitivity',
        result,
        rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots, braid_corridor)


if __name__ == '__main__':
    main()
