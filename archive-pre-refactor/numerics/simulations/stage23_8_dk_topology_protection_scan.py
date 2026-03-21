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
    field_grid_2d,
    load_runsheet,
    mean_pair_separation_series,
    ordered_peaks,
    plot_trace,
    plot_trajectories,
    selected_runs,
)
from stage23_6_dk_phase_sensitivity_braid_exchange import (
    connected_components,
    edge_grid,
    evolve,
    grade_exchange_signal,
    grade_weights,
    pair_phase_difference,
    phase_alignment_metric,
)
from stage9b_common import build_dk2d_complex, first_order_dt, pack_positions, periodic_displacement

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage23_8_dk_topology_protection_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_23_8_DK_Topology_Protection_Mechanisms_v1.md'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage23_8_topology_protection'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = BASE_FIELDS + [
    'branch_id',
    'branch_label',
    'role',
    'variant_label',
    'phase_offset_fraction_of_pi',
    'grade0_scale',
    'grade1_scale',
    'width_ratio_a_to_b',
    'anisotropy_x',
    'anisotropy_y',
    'flow_concentration_index',
    'grade_exchange_coherence',
    'channel_count',
    'recirculation_score',
    'grade_asymmetry_index',
    'phase_alignment_metric',
    'topology_label',
]

TOPOLOGY_MAPPING = {
    'braid_like_exchange': 0,
    'transfer_smeared': 1,
    'fragmented_exchange': 2,
    'unresolved_mixed_topology': 3,
}

BRANCH_NAME = {
    'A_grade_balance': 'grade-balance',
    'B_phase_edge': 'phase-edge',
    'C_symmetry_break': 'symmetry',
    'D_anisotropy_break': 'anisotropy',
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage 23.8 DK topology-protection mechanism scan.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    return parser.parse_args()


def sigma_pair(base_sigma: float, ratio_a_to_b: float) -> tuple[float, float]:
    root = math.sqrt(max(ratio_a_to_b, 1.0e-12))
    return base_sigma * root, base_sigma / root


def packet_state_with_profile(
    positions: np.ndarray,
    block_sizes: tuple[int, ...],
    center: np.ndarray,
    sigma: float,
    amplitude: float,
    phase_offset: float,
    kick_vector: np.ndarray,
    kick_cycles: float,
    grade0_scale: float,
    grade1_scale: float,
    anisotropy: tuple[float, float],
) -> np.ndarray:
    axis_x = max(sigma * anisotropy[0], 1.0e-12)
    axis_y = max(sigma * anisotropy[1], 1.0e-12)
    parts: list[np.ndarray] = []
    start = 0
    for grade_idx, size in enumerate(block_sizes):
        coords = positions[start:start + size]
        delta = periodic_displacement(coords, center)
        scaled_r2 = (delta[:, 0] / axis_x) ** 2 + (delta[:, 1] / axis_y) ** 2
        profile = np.exp(-0.5 * scaled_r2)
        base_phase = 2.0 * math.pi * kick_cycles * (delta @ kick_vector)
        phase = np.exp(1j * (base_phase + phase_offset))
        if grade_idx == 0:
            part = float(grade0_scale) * profile * phase
        elif grade_idx == 1:
            part = float(grade1_scale) * profile * phase
        else:
            part = np.zeros(size, dtype=complex)
        parts.append(np.asarray(part, dtype=complex))
        start += size
    psi = np.concatenate(parts)
    psi /= max(np.linalg.norm(psi), 1.0e-12)
    return float(amplitude) * psi


def estimate_topology(
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


def branch_sort_key(row: dict[str, Any]) -> tuple[int, int]:
    return int(row['branch_order']), int(row['variant_order'])


def short_topology(label: str) -> str:
    mapping = {
        'braid_like_exchange': 'braid',
        'transfer_smeared': 'smeared',
        'fragmented_exchange': 'fragmented',
        'unresolved_mixed_topology': 'mixed',
    }
    return mapping.get(label, label)


def humanize_protection(label: str) -> str:
    return label.replace('_', ' ')


def branch_protection_label(rows: list[dict[str, Any]]) -> str:
    rows = sorted(rows, key=branch_sort_key)
    branch_id = str(rows[0]['branch_id'])
    labels = [str(row['topology_label']) for row in rows]
    braid_count = sum(label == 'braid_like_exchange' for label in labels)

    if branch_id == 'B_phase_edge':
        if labels[0] == 'braid_like_exchange' and labels[1] == 'transfer_smeared' and labels[2] == 'braid_like_exchange':
            return 'edge_sensitive'
        if braid_count == len(labels):
            return 'protected'
        return 'unprotected_fragile'

    if braid_count == len(labels):
        return 'protected'

    baseline = labels[0]
    if baseline == 'braid_like_exchange' and any(label != 'braid_like_exchange' for label in labels[1:]):
        if branch_id == 'C_symmetry_break':
            return 'symmetry_sensitive'
        if branch_id == 'D_anisotropy_break':
            return 'anisotropy_sensitive'
        return 'unprotected_fragile'
    return 'unprotected_fragile'


def stage_outcome(branch_labels: dict[str, str]) -> str:
    clear = {'protected', 'edge_sensitive', 'symmetry_sensitive', 'anisotropy_sensitive'}
    return 'positive' if any(label in clear for label in branch_labels.values()) else 'negative_inconclusive'


def interpretation_note(branch_labels: dict[str, str]) -> str:
    clauses: list[str] = []
    if branch_labels.get('B_phase_edge') == 'edge_sensitive':
        clauses.append('phase detuning shows a reproducible braid-to-smeared band edge')
    if branch_labels.get('C_symmetry_break') == 'symmetry_sensitive':
        clauses.append('width asymmetry destabilizes the braid family')
    if branch_labels.get('D_anisotropy_break') == 'anisotropy_sensitive':
        clauses.append('weak directional anisotropy destabilizes the braid family')
    if branch_labels.get('A_grade_balance') == 'protected':
        clauses.append('grade balance perturbations preserve braid topology across the tested range')
    elif branch_labels.get('A_grade_balance') == 'unprotected_fragile':
        clauses.append('grade balance tuning matters and can weaken topology retention')
    if not clauses:
        return 'No branch produced a clean protection rule in this pass.'
    return '; '.join(clauses).capitalize() + '.'


def simulate_run(run: dict[str, Any], common: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], list[Path]]:
    resolution = int(common['resolution'])
    epsilon = 0.2
    sigma_base = float(common['mean_width'])
    amplitude = 0.5 * float(common['amplitude_scale'])
    t_final = 0.9
    kick_cycles = 1.0
    separation = 0.08
    beta = float(common['beta'])

    complex_data = build_dk2d_complex(n_side=resolution, epsilon=epsilon)
    positions = pack_positions(complex_data.points, complex_data.edge_midpoints, complex_data.face_centers)
    block_sizes = complex_data.block_sizes
    n0, n1, _ = block_sizes

    sigma_a, sigma_b = sigma_pair(sigma_base, float(run['width_ratio_a_to_b']))
    anisotropy = (float(run['anisotropy_x']), float(run['anisotropy_y']))
    phase_offset = float(run['phase_offset_fraction_of_pi']) * math.pi
    grade0_scale = float(run['grade0_scale'])
    grade1_scale = float(run['grade1_scale'])

    centers = [
        np.array([0.50 - 0.5 * separation, 0.50], dtype=float),
        np.array([0.50 + 0.5 * separation, 0.50], dtype=float),
    ]
    kicks = [np.array([0.4, 0.0], dtype=float), np.array([-0.4, 0.0], dtype=float)]

    packet_states = [
        packet_state_with_profile(
            positions=positions,
            block_sizes=block_sizes,
            center=centers[idx],
            sigma=sigma_a if idx == 0 else sigma_b,
            amplitude=amplitude,
            phase_offset=0.0 if idx == 0 else phase_offset,
            kick_vector=kicks[idx],
            kick_cycles=kick_cycles,
            grade0_scale=grade0_scale,
            grade1_scale=grade1_scale,
            anisotropy=anisotropy,
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
    peak_transfer = max(grade_amplitude, 1.0e-12)
    coherence_trace = transfer_signal / peak_transfer if transfer_signal.size else np.zeros(0, dtype=float)

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
    close_threshold = 2.0 * sigma_base
    composite_lifetime = float(np.sum(mean_pair_distances <= close_threshold) * dt)
    binding_persistence = float(np.mean(mean_pair_distances <= close_threshold)) if mean_pair_distances.size else 0.0
    final_mean_separation = float(mean_pair_distances[-1]) if mean_pair_distances.size else 0.0
    min_separation = float(np.min(min_pair_distances)) if min_pair_distances.size else 0.0
    post_collision_trend = final_mean_separation - min_separation

    threshold = max(0.15 * float(np.max(avg_edge_grid)), float(np.mean(avg_edge_grid) + np.std(avg_edge_grid)))
    channels = connected_components(avg_edge_grid >= threshold)
    braid_flag = abs(float(np.mean(braid_votes)) - 0.5) > 0.25 if braid_votes else False
    topo, concentration = estimate_topology(avg_edge_grid, braid_flag, coherence, grade_amplitude, channels)
    phase_align = phase_alignment_metric(phase_diffs)
    recirc = recirculation_score(mean_pair_distances, close_threshold)
    asym = grade_asymmetry_index(grade_hist)

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
        'new_collision_class': int(topo in {'braid_like_exchange', 'fragmented_exchange'}),
        'gate_met': 0,
        'promoted_followup': 0,
        'notes': run['notes'],
        'branch_id': run['branch_id'],
        'branch_label': run['branch_label'],
        'role': run['role'],
        'variant_label': run['variant_label'],
        'phase_offset_fraction_of_pi': float(run['phase_offset_fraction_of_pi']),
        'grade0_scale': grade0_scale,
        'grade1_scale': grade1_scale,
        'width_ratio_a_to_b': float(run['width_ratio_a_to_b']),
        'anisotropy_x': anisotropy[0],
        'anisotropy_y': anisotropy[1],
        'flow_concentration_index': concentration,
        'grade_exchange_coherence': coherence,
        'channel_count': channels,
        'recirculation_score': recirc,
        'grade_asymmetry_index': asym,
        'phase_alignment_metric': phase_align,
        'topology_label': topo,
        'branch_order': int(run['branch_order']),
        'variant_order': int(run['variant_order']),
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
        'coherence_trace': coherence_trace.tolist(),
        'peak_tracks': [hist.tolist() for hist in center_histories],
        'flow_concentration_trace': [
            float(np.sum(np.sort(g.ravel())[::-1][:max(1, int(math.ceil(0.1 * g.size)))]) / max(np.sum(g), 1.0e-12))
            for g in edge_grids
        ],
        'edge_current_map': avg_edge_grid.tolist(),
    }

    plot_paths: list[Path] = []
    traj_path = WORK_PLOT_DIR / f"stage23_8_{run['run_id']}_topology_trajectory.png"
    coh_path = WORK_PLOT_DIR / f"stage23_8_{run['run_id']}_coherence_trace.png"
    flow_path = WORK_PLOT_DIR / f"stage23_8_{run['run_id']}_flow_concentration_trace.png"
    plot_trajectories(traj_path, run['run_id'], center_histories, 'tight_clustered_pair')
    plot_trace(coh_path, run['run_id'], times, coherence_trace, 'normalized coherence', 'Stage 23.8 coherence trace')
    plot_trace(
        flow_path,
        run['run_id'],
        times,
        np.asarray(detail['flow_concentration_trace'], dtype=float),
        'flow concentration',
        'Stage 23.8 flow concentration trace',
    )
    plot_paths.extend([traj_path, coh_path, flow_path])
    return row, detail, plot_paths


def plot_protection_matrix(path: Path, rows: list[dict[str, Any]]) -> None:
    ordered = sorted(rows, key=branch_sort_key)
    branches = []
    by_branch: dict[str, list[dict[str, Any]]] = {}
    for row in ordered:
        branch_id = str(row['branch_id'])
        if branch_id not in by_branch:
            branches.append(branch_id)
            by_branch[branch_id] = []
        by_branch[branch_id].append(row)

    values = np.zeros((len(branches), 3), dtype=float)
    annotations: list[list[str]] = []
    for row_idx, branch_id in enumerate(branches):
        branch_rows = by_branch[branch_id]
        branch_annotations: list[str] = []
        for col_idx, row in enumerate(branch_rows):
            values[row_idx, col_idx] = TOPOLOGY_MAPPING[str(row['topology_label'])]
            branch_annotations.append(short_topology(str(row['topology_label'])))
        annotations.append(branch_annotations)

    fig, ax = plt.subplots(figsize=(9.8, 4.8))
    im = ax.imshow(values, aspect='auto', cmap='viridis', vmin=0.0, vmax=3.0)
    ax.set_yticks(range(len(branches)), [BRANCH_NAME[branch] for branch in branches])
    ax.set_xticks([0, 1, 2], ['baseline', 'probe 1', 'probe 2'])
    ax.set_title('Stage 23.8 branch-by-branch topology protection matrix')
    for i, branch_annotations in enumerate(annotations):
        for j, label in enumerate(branch_annotations):
            ax.text(j, i, label, ha='center', va='center', color='white', fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_braid_vs_smeared_panel(path: Path, rows: list[dict[str, Any]]) -> None:
    ordered = sorted(rows, key=branch_sort_key)
    branches = []
    braid_counts = []
    smeared_counts = []
    for branch_id in dict.fromkeys(str(row['branch_id']) for row in ordered):
        branch_rows = [row for row in ordered if str(row['branch_id']) == branch_id]
        branches.append(BRANCH_NAME[branch_id])
        braid_counts.append(sum(str(row['topology_label']) == 'braid_like_exchange' for row in branch_rows))
        smeared_counts.append(sum(str(row['topology_label']) == 'transfer_smeared' for row in branch_rows))

    x = np.arange(len(branches))
    width = 0.34
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    ax.bar(x - 0.5 * width, braid_counts, width=width, label='braid', color='tab:blue')
    ax.bar(x + 0.5 * width, smeared_counts, width=width, label='smeared', color='tab:orange')
    ax.set_xticks(x, branches)
    ax.set_ylabel('run count')
    ax.set_title('Stage 23.8 braid-vs-smeared branch comparison')
    ax.set_ylim(0, 3.4)
    ax.grid(axis='y', alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_coherence_heatmap(path: Path, details: list[dict[str, Any]]) -> None:
    ordered = sorted(details, key=lambda item: branch_sort_key(item['run']))
    max_len = max(len(detail['coherence_trace']) for detail in ordered)
    heat = np.zeros((len(ordered), max_len), dtype=float)
    labels = []
    for idx, detail in enumerate(ordered):
        trace = np.asarray(detail['coherence_trace'], dtype=float)
        heat[idx, :trace.size] = trace
        labels.append(str(detail['run']['run_id']).replace('S23_8', ''))
    fig, ax = plt.subplots(figsize=(10.6, 5.2))
    im = ax.imshow(heat, aspect='auto', cmap='magma', vmin=0.0, vmax=1.0)
    ax.set_yticks(range(len(labels)), labels)
    ax.set_xlabel('time index')
    ax.set_title('Stage 23.8 coherence robustness heatmap')
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(
    json_path: Path,
    csv_path: Path,
    rows: list[dict[str, Any]],
    stamped_plots: list[str],
    branch_labels: dict[str, str],
    outcome: str,
    interpretation: str,
) -> None:
    counts = Counter(str(row['topology_label']) for row in rows)
    ordered_branches = ['A_grade_balance', 'B_phase_edge', 'C_symmetry_break', 'D_anisotropy_break']
    lines = [
        '# Stage 23.8 DK Topology Protection Mechanisms v1',
        '',
        f"Timestamped JSON: `{json_path.relative_to(REPO_ROOT)}`",
        f"Timestamped CSV: `{csv_path.relative_to(REPO_ROOT)}`",
        '',
        f'Total runs: {len(rows)}',
        f'Stage outcome: {outcome}',
        '',
        '## Protection labels',
    ]
    for branch_id in ordered_branches:
        lines.append(f"- {BRANCH_NAME[branch_id]}: {humanize_protection(branch_labels.get(branch_id, 'unprotected_fragile'))}")
    lines.extend([
        '',
        '## Topology labels',
    ])
    for label, count in sorted(counts.items()):
        lines.append(f'- {label}: {count}')
    lines.extend([
        '',
        '## Interpretation note',
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

    rows.sort(key=branch_sort_key)
    details.sort(key=lambda item: branch_sort_key(item['run']))

    grouped: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(str(row['branch_id']), []).append(row)
    branch_labels = {branch_id: branch_protection_label(branch_rows) for branch_id, branch_rows in grouped.items()}
    outcome = stage_outcome(branch_labels)
    interpretation = interpretation_note(branch_labels)

    matrix_path = WORK_PLOT_DIR / 'stage23_8_topology_protection_matrix.png'
    comparison_path = WORK_PLOT_DIR / 'stage23_8_braid_vs_smeared_panel.png'
    heatmap_path = WORK_PLOT_DIR / 'stage23_8_coherence_robustness_heatmap.png'
    plot_protection_matrix(matrix_path, rows)
    plot_braid_vs_smeared_panel(comparison_path, rows)
    plot_coherence_heatmap(heatmap_path, details)
    plot_paths.extend([matrix_path, comparison_path, heatmap_path])

    result = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'common_fields': common,
        'runs': details,
        'summary': {
            'topology_labels': dict(Counter(row['topology_label'] for row in rows)),
            'protection_labels': branch_labels,
            'stage_outcome': outcome,
            'interpretation_note': interpretation,
        },
    }
    clean_rows = [{key: value for key, value in row.items() if key not in {'branch_order', 'variant_order'}} for row in rows]
    json_path, csv_path, stamped_plots, _ = save_atlas_payload(
        'stage23_8_dk_topology_protection',
        result,
        clean_rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots, branch_labels, outcome, interpretation)


if __name__ == '__main__':
    main()
