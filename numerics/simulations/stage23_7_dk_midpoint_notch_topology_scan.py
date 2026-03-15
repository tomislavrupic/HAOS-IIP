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
    plot_grades,
    plot_trace,
    plot_trajectories,
    selected_runs,
    weighted_center,
)
from stage23_6_dk_phase_sensitivity_braid_exchange import (
    connected_components,
    edge_grid,
    evolve,
    grade_exchange_signal,
    grade_weights,
    packet_state_with_skew,
    pair_phase_difference,
    phase_alignment_metric,
)
from stage9b_common import build_dk2d_complex, first_order_dt, pack_positions

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage23_7_dk_midpoint_notch_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_23_7_DK_Midpoint_Notch_Topology_Scan_v1.md'

CSV_FIELDS = BASE_FIELDS + [
    'role',
    'phase_offset_fraction_of_pi',
    'flow_concentration_index',
    'grade_exchange_coherence',
    'channel_count',
    'recirculation_score',
    'grade_asymmetry_index',
    'phase_alignment_metric',
    'topology_label',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage 23.7 DK midpoint-notch topology scan.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    return parser.parse_args()


def estimate_notch_topology(
    avg_edge_grid: np.ndarray,
    braid_flag: bool,
    coherence: float,
    amplitude: float,
    channels: int,
) -> tuple[str, float]:
    total = float(np.sum(avg_edge_grid))
    if total <= 1.0e-12:
        return 'unresolved_mixed_topology', 0.0
    flat = np.sort(avg_edge_grid.ravel())[::-1]
    top_k = max(1, int(math.ceil(0.1 * flat.size)))
    concentration = float(np.sum(flat[:top_k]) / total)
    if braid_flag and coherence >= 0.46 and concentration >= 0.885:
        return 'braid_like_exchange', concentration
    if channels >= 2 and coherence >= 0.44 and amplitude >= 0.17:
        return 'fragmented_exchange', concentration
    if coherence >= 0.44 and amplitude >= 0.17:
        return 'transfer_smeared', concentration
    return 'unresolved_mixed_topology', concentration


def recirculation_score(mean_pair_distances: np.ndarray, close_threshold: float) -> float:
    if len(mean_pair_distances) < 6:
        return 0.0
    tail = mean_pair_distances[len(mean_pair_distances) // 3:]
    if tail.size < 3:
        return 0.0
    centered = tail - float(np.mean(tail))
    amplitude = float(np.max(np.abs(centered)))
    if amplitude <= 1.0e-12:
        return 0.0
    zero_crossings = np.sum(centered[1:] * centered[:-1] < 0.0)
    return float(min(1.0, zero_crossings / max(1, tail.size - 1)) * float(np.mean(tail <= close_threshold)))


def grade_asymmetry_index(grade_hist: np.ndarray) -> float:
    if grade_hist.size == 0:
        return 0.0
    final0 = float(grade_hist[-1, 0])
    final1 = float(grade_hist[-1, 1]) if grade_hist.shape[1] > 1 else 0.0
    return abs(final1 - final0)


def notch_interpretation(rows: list[dict[str, Any]]) -> str:
    labels = [str(row['topology_label']) for row in rows]
    unique = list(dict.fromkeys(labels))
    if len(unique) == 1:
        return 'inconclusive'
    midpoint_idx = next((idx for idx, row in enumerate(rows) if row['run_id'] == 'S23_7_phase_midpoint'), None)
    if midpoint_idx is None:
        return 'inconclusive'
    midpoint_label = labels[midpoint_idx]
    left = labels[max(0, midpoint_idx - 1)]
    right = labels[min(len(labels) - 1, midpoint_idx + 1)]
    if midpoint_label != left and midpoint_label != right and left == right:
        return 'sharp'
    smeared_band = sum(label != labels[0] for label in labels[1:-1])
    if smeared_band >= 3:
        return 'broad'
    if len(unique) >= 3:
        return 'multi-layered'
    return 'inconclusive'


def simulate_run(run: dict[str, Any], common: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], list[Path]]:
    resolution = int(common['resolution'])
    epsilon = 0.2
    sigma = 0.08
    amplitude = 0.5 * float(common['amplitude_scale'])
    t_final = 0.9
    kick_cycles = 1.0
    separation = 0.08
    beta = float(common['beta'])
    skew = float(common['skew'])

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
    braid_votes: list[float] = []
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

    mean_pair_distances, min_pair_distances = mean_pair_separation_series(center_histories)
    close_threshold = 2.0 * sigma
    composite_lifetime = float(np.sum(mean_pair_distances <= close_threshold) * dt)
    binding_persistence = float(np.mean(mean_pair_distances <= close_threshold)) if mean_pair_distances.size else 0.0
    final_mean_separation = float(mean_pair_distances[-1]) if mean_pair_distances.size else 0.0
    min_separation = float(np.min(min_pair_distances)) if min_pair_distances.size else 0.0
    post_collision_trend = final_mean_separation - min_separation

    threshold = max(0.15 * float(np.max(avg_edge_grid)), float(np.mean(avg_edge_grid) + np.std(avg_edge_grid)))
    channels = connected_components(avg_edge_grid >= threshold)
    braid_flag = abs(float(np.mean(braid_votes)) - 0.5) > 0.25 if braid_votes else False
    topo, concentration = estimate_notch_topology(avg_edge_grid, braid_flag, coherence, grade_amplitude, channels)
    phase_align = phase_alignment_metric(phase_diffs)
    recirc = recirculation_score(mean_pair_distances, close_threshold)
    asym = grade_asymmetry_index(grade_hist)

    row = {
        'run_id': run['run_id'],
        'geometry_id': 'tight_clustered_pair',
        'phase_id': str(run['phase_offset_fraction_of_pi']),
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
        'new_collision_class': int(topo in {'braid_like_exchange', 'fragmented_exchange'}),
        'gate_met': 0,
        'promoted_followup': 0,
        'notes': run['notes'],
        'role': run['role'],
        'phase_offset_fraction_of_pi': float(run['phase_offset_fraction_of_pi']),
        'flow_concentration_index': concentration,
        'grade_exchange_coherence': coherence,
        'channel_count': channels,
        'recirculation_score': recirc,
        'grade_asymmetry_index': asym,
        'phase_alignment_metric': phase_align,
        'topology_label': topo,
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
    }

    plot_paths: list[Path] = []
    traj_path = WORK_PLOT_DIR / f"stage23_7_{run['run_id']}_collision_trajectory.png"
    grade_path = WORK_PLOT_DIR / f"stage23_7_{run['run_id']}_grade_weight_trace.png"
    flow_path = WORK_PLOT_DIR / f"stage23_7_{run['run_id']}_flow_field.png"
    plot_trajectories(traj_path, run['run_id'], center_histories, 'tight_clustered_pair')
    plot_grades(grade_path, run['run_id'], times, grade_hist)
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    im = ax.imshow(avg_edge_grid, cmap='viridis', origin='lower', aspect='auto')
    ax.set_title(f'Local flow field: {run["run_id"]}')
    ax.set_xlabel('y cell')
    ax.set_ylabel('x cell')
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    fig.tight_layout()
    fig.savefig(flow_path, dpi=180, bbox_inches='tight')
    plt.close(fig)
    plot_paths.extend([traj_path, grade_path, flow_path])
    return row, detail, plot_paths


def plot_notch_classification(path: Path, rows: list[dict[str, Any]]) -> None:
    mapping = {
        'braid_like_exchange': 0,
        'transfer_smeared': 1,
        'fragmented_exchange': 2,
        'unresolved_mixed_topology': 3,
    }
    xs = [float(row['phase_offset_fraction_of_pi']) for row in rows]
    ys = [mapping[str(row['topology_label'])] for row in rows]
    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    ax.plot(xs, ys, marker='o', color='tab:blue')
    ax.set_xlabel('phase offset fraction of pi')
    ax.set_ylabel('topology class index')
    ax.set_yticks(list(mapping.values()), list(mapping.keys()))
    ax.set_title('Stage 23.7 topology classification vs phase notch')
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_flow_panel(path: Path, rows: list[dict[str, Any]]) -> None:
    xs = [float(row['phase_offset_fraction_of_pi']) for row in rows]
    ys = [float(row['flow_concentration_index']) for row in rows]
    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    ax.plot(xs, ys, marker='o', color='tab:orange')
    ax.set_xlabel('phase offset fraction of pi')
    ax.set_ylabel('flow concentration index')
    ax.set_title('Stage 23.7 flow concentration panel')
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_coherence_curve(path: Path, rows: list[dict[str, Any]]) -> None:
    xs = [float(row['phase_offset_fraction_of_pi']) for row in rows]
    ys = [float(row['grade_exchange_coherence']) for row in rows]
    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    ax.plot(xs, ys, marker='o', color='tab:green')
    ax.set_xlabel('phase offset fraction of pi')
    ax.set_ylabel('grade exchange coherence')
    ax.set_title('Stage 23.7 coherence vs detuning')
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_phase_table(path: Path, rows: list[dict[str, Any]]) -> None:
    mapping = {
        'braid_like_exchange': 0,
        'transfer_smeared': 1,
        'fragmented_exchange': 2,
        'unresolved_mixed_topology': 3,
    }
    vals = np.asarray([[mapping[str(row['topology_label'])]] for row in rows], dtype=float)
    labels = [f"{float(row['phase_offset_fraction_of_pi']):.3f}" for row in rows]
    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    im = ax.imshow(vals, aspect='auto', cmap='viridis', vmin=0.0, vmax=3.0)
    ax.set_yticks(range(len(labels)), labels)
    ax.set_xticks([0], ['topology'])
    ax.set_title('Stage 23.7 phase-corridor topology table')
    for idx, row in enumerate(rows):
        ax.text(0, idx, str(row['topology_label']).replace('_', ' '), ha='center', va='center', color='white', fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(json_path: Path, csv_path: Path, rows: list[dict[str, Any]], stamped_plots: list[str], interpretation: str) -> None:
    counts = Counter(str(row['topology_label']) for row in rows)
    lines = [
        '# Stage 23.7 DK Midpoint Notch Topology Scan v1',
        '',
        f"Timestamped JSON: `{json_path.relative_to(REPO_ROOT)}`",
        f"Timestamped CSV: `{csv_path.relative_to(REPO_ROOT)}`",
        '',
        f'Total runs: {len(rows)}',
        '',
        '## Topology labels',
    ]
    for label, count in sorted(counts.items()):
        lines.append(f'- {label}: {count}')
    lines.extend([
        '',
        '## Midpoint corridor interpretation',
        f'- {interpretation}',
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

    rows.sort(key=lambda item: float(item['phase_offset_fraction_of_pi']))
    details.sort(key=lambda item: float(item['run']['phase_offset_fraction_of_pi']))
    interpretation = notch_interpretation(rows)

    class_path = WORK_PLOT_DIR / 'stage23_7_topology_classification_vs_phase.png'
    flow_path = WORK_PLOT_DIR / 'stage23_7_flow_concentration_panel.png'
    coh_path = WORK_PLOT_DIR / 'stage23_7_coherence_vs_detuning.png'
    table_path = WORK_PLOT_DIR / 'stage23_7_phase_corridor_topology_table.png'
    plot_notch_classification(class_path, rows)
    plot_flow_panel(flow_path, rows)
    plot_coherence_curve(coh_path, rows)
    plot_phase_table(table_path, rows)
    plot_paths.extend([class_path, flow_path, coh_path, table_path])

    result = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'common_fields': common,
        'runs': details,
        'summary': {
            'topology_labels': dict(Counter(row['topology_label'] for row in rows)),
            'midpoint_interpretation': interpretation,
        },
    }
    json_path, csv_path, stamped_plots, _ = save_atlas_payload(
        'stage23_7_dk_midpoint_notch',
        result,
        rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots, interpretation)


if __name__ == '__main__':
    main()
