#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt, save_atlas_payload
from stage23_1_dk_collision_geometry_scan import (
    CSV_FIELDS as BASE_FIELDS,
    WORK_PLOT_DIR,
    basin_id,
    field_grid_2d,
    load_runsheet,
    longest_constant_run,
    mean_pair_separation_series,
    ordered_peaks,
    plot_grades,
    plot_trace,
    plot_trajectories,
    selected_runs,
    weighted_center,
    packet_state,
)
from stage9b_common import build_dk2d_complex, first_order_dt, pack_positions, periodic_displacement

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage23_5_dk_energy_flow_topology_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_23_5_DK_Energy_Flow_Topology_Mapping_v1.md'

CSV_FIELDS = BASE_FIELDS + [
    'family_id',
    'variation_class',
    'initial_phase_offset',
    'width_ratio_a_to_b',
    'separation_scale',
    'edge_current_map_norm',
    'flow_concentration_index',
    'recirculation_score',
    'channel_count',
    'loop_count',
    'grade_asymmetry_index',
    'topology_label',
    'secondary_label',
    'stable_topology_family',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage 23.5 DK energy-flow topology mapping scan.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    return parser.parse_args()


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


def transport_span_series(center_series: np.ndarray) -> np.ndarray:
    if len(center_series) == 0:
        return np.zeros(0, dtype=float)
    cumulative = [0.0]
    for prev, nxt in zip(center_series[:-1], center_series[1:]):
        step = periodic_displacement(nxt[None, :], prev)[0]
        cumulative.append(cumulative[-1] + float(np.linalg.norm(step)))
    return np.asarray(cumulative, dtype=float)


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


def grade_exchange_signal(grade_hist: np.ndarray) -> tuple[np.ndarray, float]:
    signal = np.linalg.norm(grade_hist - grade_hist[0][None, :], axis=1)
    peak = float(np.max(signal)) if signal.size else 0.0
    if peak <= 1.0e-12:
        return signal, 0.0
    coherence = float(np.mean(signal) / peak)
    return signal, max(0.0, min(1.0, coherence))


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


def estimate_flow_topology(
    avg_edge_grid: np.ndarray,
    transport_span: float,
    recurrence: float,
    coherence: float,
    grade_amplitude: float,
    braid_flag: bool,
) -> tuple[str, float, float, int, int]:
    total = float(np.sum(avg_edge_grid))
    if total <= 1.0e-12:
        return 'diffuse_washout', 0.0, 0.0, 0, 0
    flat = np.sort(avg_edge_grid.ravel())[::-1]
    top_k = max(1, int(math.ceil(0.1 * flat.size)))
    concentration = float(np.sum(flat[:top_k]) / total)
    threshold = max(0.15 * float(np.max(avg_edge_grid)), float(np.mean(avg_edge_grid) + np.std(avg_edge_grid)))
    mask = avg_edge_grid >= threshold
    channels = connected_components(mask)
    path_length_proxy = max(transport_span, 1.0e-12)
    net_disp_proxy = max(1.0e-12, transport_span * (1.0 - recurrence))
    recirculation = float(np.clip(1.0 - net_disp_proxy / path_length_proxy, 0.0, 1.0)) if transport_span > 1.0e-12 else 0.0
    loops = 1 if recirculation >= 0.45 and recurrence >= 0.12 else 0

    if loops >= 1:
        label = 'looped_circulation'
    elif braid_flag and coherence >= 0.45:
        label = 'braid_like_exchange'
    elif channels >= 2 and concentration >= 0.28:
        label = 'split_channel_transport'
    elif channels == 1 and concentration >= 0.30 and transport_span >= 0.01:
        label = 'single_channel_transport'
    elif concentration < 0.22 and grade_amplitude < 0.08:
        label = 'diffuse_washout'
    else:
        label = 'unresolved_mixed_topology'
    return label, concentration, recirculation, channels, loops


def secondary_label(grade_hist: np.ndarray, grade_amplitude: float) -> tuple[str, float]:
    final0 = float(grade_hist[-1, 0]) if grade_hist.size else 0.0
    final1 = float(grade_hist[-1, 1]) if grade_hist.shape[1] > 1 else 0.0
    asym = final1 - final0
    asym_index = abs(asym)
    if grade_amplitude >= 0.12:
        return 'transfer_dominant', asym_index
    if asym >= 0.05:
        return 'one_form_dominant', asym_index
    if asym <= -0.05:
        return 'zero_form_dominant', asym_index
    return 'grade_balanced', asym_index


def phase_offset_radians(value: float) -> float:
    return float(2.0 * math.pi * value)


def separation_value(scale: str) -> float:
    mapping = {
        'tight_baseline': 0.08,
        'medium_clustered': 0.12,
        'near_overlap': 0.04,
    }
    return mapping[scale]


def sigma_pair(base_sigma: float, ratio_a_to_b: float) -> tuple[float, float]:
    root = math.sqrt(ratio_a_to_b)
    return base_sigma * root, base_sigma / root


def evolve(base_operator, state0: np.ndarray, dt: float, steps: int) -> list[np.ndarray]:
    state = state0.copy()
    states = [state.copy()]
    for _ in range(steps):
        state = state - dt * (base_operator @ state)
        norm = np.linalg.norm(state)
        if norm > 1.0e-12:
            state = state / norm * np.linalg.norm(states[0])
        states.append(state.copy())
    return states


def simulate_run(run: dict[str, Any], common: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], list[Path]]:
    resolution = int(common['resolution'])
    epsilon = 0.2
    base_sigma = 0.08
    amplitude = 0.5
    t_final = 0.9
    kick_cycles = 1.0

    complex_data = build_dk2d_complex(n_side=resolution, epsilon=epsilon)
    positions = pack_positions(complex_data.points, complex_data.edge_midpoints, complex_data.face_centers)
    block_sizes = complex_data.block_sizes
    n0, n1, _ = block_sizes

    sep = separation_value(str(run['separation_scale']))
    sigma_a, sigma_b = sigma_pair(base_sigma, float(run['width_ratio_a_to_b']))
    centers = [
        np.array([0.50 - 0.5 * sep, 0.50], dtype=float),
        np.array([0.50 + 0.5 * sep, 0.50], dtype=float),
    ]
    kicks = [np.array([0.4, 0.0], dtype=float), np.array([-0.4, 0.0], dtype=float)]
    phase_offset = phase_offset_radians(float(run['initial_phase_offset']))

    packet_states = [
        packet_state(
            positions=positions,
            block_sizes=block_sizes,
            center=centers[idx],
            sigma=sigma_a if idx == 0 else sigma_b,
            amplitude=amplitude,
            phase_offset=0.0 if idx == 0 else phase_offset,
            kick_vector=kicks[idx],
            kick_cycles=kick_cycles,
            active_grades=(0, 1),
        )
        for idx in range(2)
    ]

    psi0 = np.sum(packet_states, axis=0)
    dt = first_order_dt(complex_data.dirac_kahler, 0.35)
    steps = max(2, int(math.ceil(t_final / dt)))
    states = evolve(complex_data.dirac_kahler, psi0, dt, steps)
    times = [idx * dt for idx in range(len(states))]

    grade_hist = np.asarray([grade_weights(state, block_sizes) for state in states], dtype=float)
    transfer_signal, coherence = grade_exchange_signal(grade_hist)
    grade_amplitude = float(np.max(transfer_signal)) if transfer_signal.size else 0.0

    peak_tracks: list[list[np.ndarray]] = [[] for _ in range(2)]
    raw_peak_counts: list[int] = []
    previous: list[np.ndarray] | None = None
    anchors = [center.copy() for center in centers]
    braid_votes = []
    edge_grids = []
    phase_diffs = []
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
    basin_ids = [basin_id(center) for center in total_centers]
    coarse_basin_persistence = longest_constant_run(basin_ids)
    transport_series = transport_span_series(total_centers)
    transport_span = float(transport_series[-1]) if transport_series.size else 0.0

    sample_dt = float(np.median(np.diff(times))) if len(times) >= 2 else dt
    mean_pair_distances, min_pair_distances = mean_pair_separation_series(center_histories)
    close_threshold = sigma_a + sigma_b
    composite_lifetime = float(np.sum(mean_pair_distances <= close_threshold) * sample_dt)
    binding_persistence = float(np.mean(mean_pair_distances <= close_threshold)) if mean_pair_distances.size else 0.0
    encounter_dwell_time = composite_lifetime
    min_separation = float(np.min(min_pair_distances)) if min_pair_distances.size else 0.0
    final_mean_separation = float(mean_pair_distances[-1]) if mean_pair_distances.size else 0.0
    post_collision_trend = final_mean_separation - min_separation
    initial_peak_count = int(raw_peak_counts[0]) if raw_peak_counts else 0

    recurrence = recurrence_indicator(mean_pair_distances, close_threshold)
    braid_flag = abs(float(np.mean(braid_votes)) - 0.5) > 0.25 if braid_votes else False
    topology_label, concentration, recirculation, channels, loops = estimate_flow_topology(
        avg_edge_grid=avg_edge_grid,
        transport_span=transport_span,
        recurrence=recurrence,
        coherence=coherence,
        grade_amplitude=grade_amplitude,
        braid_flag=braid_flag,
    )
    sec_label, asym_index = secondary_label(grade_hist, grade_amplitude)

    row = {
        'run_id': run['run_id'],
        'geometry_id': 'tight_clustered_pair',
        'phase_id': str(run['variation_class']),
        'resolution': resolution,
        'graph_type': common['graph_type'],
        'packet_count': 2,
        'initial_peak_count': initial_peak_count,
        'collision_label': topology_label,
        'persistence_label': 'weak persistence gain' if encounter_dwell_time >= 0.15 * float(times[-1]) else 'no persistence gain',
        'composite_lifetime': composite_lifetime,
        'binding_persistence': binding_persistence,
        'corridor_dwell': 0.0,
        'coarse_basin_persistence': coarse_basin_persistence,
        'minimum_separation': min_separation,
        'final_mean_separation': final_mean_separation,
        'post_collision_separation_trend': post_collision_trend,
        'encounter_dwell_time': encounter_dwell_time,
        'deflection_angle_proxy': 0.0,
        'reflection_fraction': 0.0,
        'grade_transfer_amplitude': grade_amplitude,
        'omega0_weight_initial': float(grade_hist[0, 0]),
        'omega1_weight_initial': float(grade_hist[0, 1]),
        'omega2_weight_initial': float(grade_hist[0, 2]),
        'omega0_weight_final': float(grade_hist[-1, 0]),
        'omega1_weight_final': float(grade_hist[-1, 1]),
        'omega2_weight_final': float(grade_hist[-1, 2]),
        'localized_circulation_proxy': recirculation,
        'new_collision_class': int(topology_label not in {'diffuse_washout', 'unresolved_mixed_topology'}),
        'gate_met': 0,
        'promoted_followup': 0,
        'notes': run['notes'],
        'family_id': run['family_id'],
        'variation_class': run['variation_class'],
        'initial_phase_offset': float(run['initial_phase_offset']),
        'width_ratio_a_to_b': float(run['width_ratio_a_to_b']),
        'separation_scale': str(run['separation_scale']),
        'edge_current_map_norm': float(np.linalg.norm(avg_edge_grid)),
        'flow_concentration_index': concentration,
        'recirculation_score': recirculation,
        'channel_count': channels,
        'loop_count': loops,
        'grade_asymmetry_index': asym_index,
        'topology_label': topology_label,
        'secondary_label': sec_label,
        'stable_topology_family': 0,
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
        'transport_span_series': transport_series.tolist(),
        'grade_flow_trace': transfer_signal.tolist(),
        'edge_current_map': avg_edge_grid.tolist(),
    }

    plot_paths: list[Path] = []
    traj_path = WORK_PLOT_DIR / f"stage23_5_{run['run_id']}_collision_trajectory.png"
    sep_path = WORK_PLOT_DIR / f"stage23_5_{run['run_id']}_separation_trace.png"
    grade_path = WORK_PLOT_DIR / f"stage23_5_{run['run_id']}_grade_weight_trace.png"
    flow_path = WORK_PLOT_DIR / f"stage23_5_{run['run_id']}_grade_flow_trace.png"
    edge_path = WORK_PLOT_DIR / f"stage23_5_{run['run_id']}_edge_current_map.png"
    plot_trajectories(traj_path, run['run_id'], center_histories, 'tight_clustered_pair')
    plot_trace(sep_path, run['run_id'], times, mean_pair_distances, 'mean peak separation', 'Stage 23.5 separation trace')
    plot_grades(grade_path, run['run_id'], times, grade_hist)
    plot_trace(flow_path, run['run_id'], times, transfer_signal, 'grade flow trace', 'Stage 23.5 grade-flow trace')
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    im = ax.imshow(avg_edge_grid, cmap='viridis', origin='lower', aspect='auto')
    ax.set_title(f'Edge current map: {run["run_id"]}')
    ax.set_xlabel('y cell')
    ax.set_ylabel('x cell')
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    fig.tight_layout()
    fig.savefig(edge_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.extend([traj_path, sep_path, grade_path, flow_path, edge_path])
    return row, detail, plot_paths


def plot_topology_matrix(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [row['run_id'].replace('S23_5_', '') for row in rows]
    mapping = {
        'single_channel_transport': 0,
        'split_channel_transport': 1,
        'looped_circulation': 2,
        'braid_like_exchange': 3,
        'diffuse_washout': 4,
        'unresolved_mixed_topology': 5,
    }
    values = np.asarray([[mapping[str(row['topology_label'])]] for row in rows], dtype=float)
    fig, ax = plt.subplots(figsize=(8.8, 5.0))
    im = ax.imshow(values, aspect='auto', cmap='viridis', vmin=0.0, vmax=5.0)
    ax.set_yticks(range(len(labels)), labels)
    ax.set_xticks([0], ['topology'])
    ax.set_title('Stage 23.5 topology class matrix')
    for idx, row in enumerate(rows):
        ax.text(0, idx, str(row['topology_label']).replace('_', ' '), ha='center', va='center', color='white', fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_flow_panel(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [row['run_id'].replace('S23_5_', '') for row in rows]
    concentration = [float(row['flow_concentration_index']) for row in rows]
    span = [float(row['encounter_dwell_time']) for row in rows]
    fig, axes = plt.subplots(2, 1, figsize=(10.5, 6.8), sharex=True)
    axes[0].bar(range(len(rows)), concentration, color='tab:blue')
    axes[0].set_ylabel('flow concentration')
    axes[0].set_title('Stage 23.5 flow concentration panel')
    axes[0].grid(axis='y', alpha=0.25)
    axes[1].bar(range(len(rows)), span, color='tab:orange')
    axes[1].set_ylabel('encounter dwell')
    axes[1].set_xticks(range(len(rows)), labels, rotation=45, ha='right')
    axes[1].grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_grade_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [row['run_id'].replace('S23_5_', '') for row in rows]
    asym = [float(row['grade_asymmetry_index']) for row in rows]
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    ax.bar(range(len(rows)), asym, color='tab:green')
    ax.set_ylabel('grade asymmetry index')
    ax.set_xticks(range(len(rows)), labels, rotation=45, ha='right')
    ax.set_title('Stage 23.5 grade asymmetry summary')
    ax.grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(json_path: Path, csv_path: Path, rows: list[dict[str, Any]], stamped_plots: list[str], stable_families: list[str]) -> None:
    topology_counts = Counter(str(row['topology_label']) for row in rows)
    secondary_counts = Counter(str(row['secondary_label']) for row in rows)
    lines = [
        '# Stage 23.5 DK Energy Flow Topology Mapping v1',
        '',
        f"Timestamped JSON: `{json_path.relative_to(REPO_ROOT)}`",
        f"Timestamped CSV: `{csv_path.relative_to(REPO_ROOT)}`",
        '',
        f'Total runs: {len(rows)}',
        '',
        '## Topology labels',
    ]
    for label, count in sorted(topology_counts.items()):
        lines.append(f'- {label}: {count}')
    lines.extend(['', '## Secondary labels'])
    for label, count in sorted(secondary_counts.items()):
        lines.append(f'- {label}: {count}')
    lines.extend([
        '',
        '## Stable topology families',
        f"- {stable_families if stable_families else 'none'}",
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

    by_variation: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_variation[str(row['variation_class'])].append(row)
    stable_families: list[str] = []
    for variation, group in by_variation.items():
        labels = Counter(str(row['topology_label']) for row in group)
        most_common, count = labels.most_common(1)[0]
        if most_common not in {'diffuse_washout', 'unresolved_mixed_topology'} and count >= 2:
            stable_families.append(variation)
            for row in group:
                if row['topology_label'] == most_common:
                    row['stable_topology_family'] = 1

    matrix_path = WORK_PLOT_DIR / 'stage23_5_topology_class_matrix.png'
    flow_path = WORK_PLOT_DIR / 'stage23_5_flow_concentration_panel.png'
    grade_path = WORK_PLOT_DIR / 'stage23_5_grade_asymmetry_summary.png'
    plot_topology_matrix(matrix_path, rows)
    plot_flow_panel(flow_path, rows)
    plot_grade_summary(grade_path, rows)
    plot_paths.extend([matrix_path, flow_path, grade_path])

    result = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'common_fields': common,
        'runs': details,
        'summary': {
            'topology_labels': dict(Counter(row['topology_label'] for row in rows)),
            'secondary_labels': dict(Counter(row['secondary_label'] for row in rows)),
            'stable_families': stable_families,
        },
    }
    json_path, csv_path, stamped_plots, _ = save_atlas_payload(
        'stage23_5_dk_energy_flow_topology',
        result,
        rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots, stable_families)


if __name__ == '__main__':
    main()
