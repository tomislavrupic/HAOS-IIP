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

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt, save_atlas_payload
from stage23_1_dk_collision_geometry_scan import (
    CSV_FIELDS as BASE_FIELDS,
    WORK_PLOT_DIR,
    field_grid_2d,
    mean_pair_separation_series,
    ordered_peaks,
    plot_grades,
    plot_trace,
    plot_trajectories,
    weighted_center,
)
from stage23_6_dk_phase_sensitivity_braid_exchange import (
    edge_grid,
    estimate_topology,
    evolve,
    grade_exchange_signal,
    grade_weights,
    pair_phase_difference,
    phase_alignment_metric,
    recurrence_indicator,
    packet_state_with_skew,
)
from stage9b_common import build_dk2d_complex, first_order_dt, pack_positions

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_10_harmonic_address_compatibility_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_C0_10_Harmonic_Address_Compatibility_Probe_v1.md'

CSV_FIELDS = BASE_FIELDS + [
    'group_id',
    'group_label',
    'variant_label',
    'field_mode',
    'address_modulus',
    'left_address',
    'right_address',
    'persistence_gain',
    'topology_class',
    'braid_exchange_indicator',
    'transfer_smear_index',
    'flow_concentration_index',
    'encounter_sorting_score',
    'address_selectivity_index',
    'classification_label',
    'grade_exchange_coherence',
    'channel_count',
    'loop_count',
    'recirculation_score',
    'phase_alignment_metric',
    'address_same_activity',
    'address_near_activity',
    'address_mismatch_activity',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.10 harmonic-address compatibility probe.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def selected_runs(runs: list[dict[str, Any]], run_ids: list[str]) -> list[dict[str, Any]]:
    if not run_ids:
        return runs
    wanted = set(run_ids)
    return [run for run in runs if str(run['run_id']) in wanted]


def grade_index_map(block_sizes: tuple[int, ...]) -> np.ndarray:
    labels = np.zeros(sum(block_sizes), dtype=int)
    start = 0
    for grade, size in enumerate(block_sizes):
        labels[start:start + size] = grade
        start += size
    return labels


def minimal_mod_distance(a: np.ndarray, b: np.ndarray, modulus: int) -> np.ndarray:
    diff = np.mod(a - b, modulus)
    return np.minimum(diff, modulus - diff)


def compatibility_weights(distances: np.ndarray, eta_match: float, eta_mismatch: float) -> np.ndarray:
    weights = np.full(distances.shape, float(eta_mismatch), dtype=float)
    weights[distances == 0] = 1.0
    weights[distances == 1] = float(eta_match)
    return weights


def make_rng(run_id: str) -> np.random.Generator:
    seed = sum((idx + 1) * ord(ch) for idx, ch in enumerate(run_id)) + 101
    return np.random.default_rng(seed)


def build_address_labels(
    packet_states: list[np.ndarray],
    block_sizes: tuple[int, ...],
    left_address: int,
    right_address: int,
    modulus: int,
    field_mode: str,
    support_floor_fraction: float,
    dominant_ratio_threshold: float,
    rng: np.random.Generator,
) -> np.ndarray:
    packet_power = np.stack([np.abs(packet) ** 2 for packet in packet_states], axis=1)
    total_power = np.sum(packet_power, axis=1)
    max_power = np.max(total_power) if total_power.size else 0.0
    active = total_power >= float(support_floor_fraction) * max(max_power, 1.0e-12)
    dominant = np.argmax(packet_power, axis=1)
    dominant_power = packet_power[np.arange(packet_power.shape[0]), dominant]
    dominant_ratio = dominant_power / np.maximum(total_power, 1.0e-12)
    strong = active & (dominant_ratio >= float(dominant_ratio_threshold))
    mixed = active & ~strong
    grade_ids = grade_index_map(block_sizes)

    labels = np.full(total_power.shape, int(left_address), dtype=int)
    labels[active & (dominant == 0)] = int(left_address)
    labels[active & (dominant == 1)] = int(right_address)

    if field_mode == 'packet_dominant':
        labels[mixed] = np.where(dominant[mixed] == 0, int(left_address), int(right_address))
        return labels

    if field_mode == 'symmetry_broken':
        labels[mixed] = (labels[mixed] + 1 + grade_ids[mixed]) % modulus
        labels[~active] = (int(right_address) + 2 + grade_ids[~active]) % modulus
        return labels

    if field_mode == 'randomized':
        labels[:] = rng.integers(0, modulus, size=labels.size, dtype=int)
        return labels

    raise ValueError(f'unknown field_mode: {field_mode}')


def address_weighted_operator(
    base_operator: sp.csr_matrix,
    labels: np.ndarray,
    modulus: int,
    eta_match: float,
    eta_mismatch: float,
) -> sp.csr_matrix:
    coo = base_operator.tocoo()
    distances = minimal_mod_distance(labels[coo.row], labels[coo.col], modulus)
    weights = compatibility_weights(distances, eta_match, eta_mismatch)
    weighted = sp.csr_matrix((coo.data * weights, (coo.row, coo.col)), shape=base_operator.shape)
    return weighted.tocsr()


def address_activity_statistics(
    states: list[np.ndarray],
    operator: sp.csr_matrix,
    labels: np.ndarray,
    modulus: int,
    eta_match: float,
    eta_mismatch: float,
) -> tuple[float, float, float, float, float]:
    coo = operator.tocoo()
    mask = coo.row < coo.col
    if not np.any(mask):
        return 0.0, 0.0, 0.0, 0.0, 0.0
    row = coo.row[mask]
    col = coo.col[mask]
    distances = minimal_mod_distance(labels[row], labels[col], modulus)
    same_mask = distances == 0
    near_mask = distances == 1
    mismatch_mask = distances >= 2
    activity = np.zeros(row.shape[0], dtype=float)
    for state in states:
        activity += np.abs(state[row]) * np.abs(state[col])
    total = float(np.sum(activity))
    if total <= 1.0e-12:
        return 0.0, 0.0, 0.0, 0.0, 0.0
    same = float(np.sum(activity[same_mask]) / total)
    near = float(np.sum(activity[near_mask]) / total)
    mismatch = float(np.sum(activity[mismatch_mask]) / total)
    weighted_mean = same + float(eta_match) * near + float(eta_mismatch) * mismatch
    selectivity = float(np.clip((weighted_mean - float(eta_mismatch)) / max(1.0 - float(eta_mismatch), 1.0e-12), 0.0, 1.0))
    sorting = float(np.clip(same - mismatch, -1.0, 1.0))
    return same, near, mismatch, sorting, selectivity


def transfer_smear_index(concentration: float, coherence: float) -> float:
    return float(np.clip(1.0 - 0.5 * (float(concentration) + float(coherence)), 0.0, 1.0))


def classify_row(
    row: dict[str, Any],
    compatible_braid_present: bool,
    compatible_smear_mean: float,
) -> str:
    if str(row['group_id']) == 'fully_compatible' and int(row['braid_exchange_indicator']):
        return 'address_protected_braid'
    if float(row['encounter_sorting_score']) >= 0.12 and float(row['address_selectivity_index']) >= 0.5:
        return 'address_sorted_encounter'
    if compatible_braid_present and not int(row['braid_exchange_indicator']):
        if str(row['topology_class']) == 'transfer_smeared' or float(row['transfer_smear_index']) >= compatible_smear_mean + 0.04:
            return 'address_induced_smear'
        return 'address_fragile_braid'
    return 'no_address_effect'


def simulate_run(run: dict[str, Any], common: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], list[Path]]:
    resolution = int(common['resolution'])
    epsilon = float(common['epsilon'])
    sigma = float(common['mean_width'])
    amplitude = float(common['amplitude'])
    t_final = float(common['t_final'])
    dt_scale = float(common['dt_scale'])
    kick_cycles = float(common['kick_cycles'])
    phase_offset = float(common['phase_offset_fraction_of_pi']) * math.pi
    skew = float(common['local_phase_skew_fraction'])
    beta = float(common['beta'])
    separation = float(common['separation'])
    modulus = int(common['address_modulus'])
    eta_match = float(common['eta_match'])
    eta_mismatch = float(common['eta_mismatch'])
    support_floor_fraction = float(common['address_support_floor_fraction'])
    dominant_ratio_threshold = float(common['dominant_ratio_threshold'])

    complex_data = build_dk2d_complex(n_side=resolution, epsilon=epsilon)
    positions = pack_positions(complex_data.points, complex_data.edge_midpoints, complex_data.face_centers)
    block_sizes = complex_data.block_sizes
    n0, n1, _ = block_sizes

    centers = [
        np.array([0.50 - 0.5 * separation, 0.50], dtype=float),
        np.array([0.50 + 0.5 * separation, 0.50], dtype=float),
    ]
    kicks = [np.array([0.4, 0.0], dtype=float), np.array([-0.4, 0.0], dtype=float)]
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

    rng = make_rng(str(run['run_id']))
    labels = build_address_labels(
        packet_states=packet_states,
        block_sizes=block_sizes,
        left_address=int(run['left_address']),
        right_address=int(run['right_address']),
        modulus=modulus,
        field_mode=str(run['field_mode']),
        support_floor_fraction=support_floor_fraction,
        dominant_ratio_threshold=dominant_ratio_threshold,
        rng=rng,
    )

    base_operator = complex_data.dirac_kahler
    effective_operator = address_weighted_operator(base_operator, labels, modulus, eta_match, eta_mismatch)
    dt = first_order_dt(base_operator, dt_scale)

    psi0 = np.sum(packet_states, axis=0)
    steps = max(2, int(math.ceil(t_final / dt)))
    states = evolve(effective_operator, block_sizes, psi0, dt, steps, beta)
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
    same_activity, near_activity, mismatch_activity, sorting_score, selectivity = address_activity_statistics(
        states=states,
        operator=effective_operator,
        labels=labels,
        modulus=modulus,
        eta_match=eta_match,
        eta_mismatch=eta_mismatch,
    )
    smear_index = transfer_smear_index(concentration, coherence)

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
        'group_id': run['group_id'],
        'group_label': run['group_label'],
        'variant_label': run['variant_label'],
        'field_mode': run['field_mode'],
        'address_modulus': modulus,
        'left_address': int(run['left_address']),
        'right_address': int(run['right_address']),
        'persistence_gain': 0.0,
        'topology_class': topo,
        'braid_exchange_indicator': int(topo == 'braid_like_exchange'),
        'transfer_smear_index': smear_index,
        'flow_concentration_index': concentration,
        'encounter_sorting_score': sorting_score,
        'address_selectivity_index': selectivity,
        'classification_label': 'pending',
        'grade_exchange_coherence': coherence,
        'channel_count': channels,
        'loop_count': loops,
        'recirculation_score': recirc,
        'phase_alignment_metric': phase_align,
        'address_same_activity': same_activity,
        'address_near_activity': near_activity,
        'address_mismatch_activity': mismatch_activity,
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
        'address_label_histogram': {str(key): int(value) for key, value in Counter(labels.tolist()).items()},
        'address_labels': labels.tolist(),
    }

    plot_paths: list[Path] = []
    stem = f"stage_c0_10_{run['run_id']}"
    traj_path = WORK_PLOT_DIR / f'{stem}_topology_trajectory.png'
    flow_path = WORK_PLOT_DIR / f'{stem}_flow_concentration_trace.png'
    grade_path = WORK_PLOT_DIR / f'{stem}_grade_exchange_trace.png'
    grades_path = WORK_PLOT_DIR / f'{stem}_grade_weights.png'
    plot_trajectories(traj_path, run['run_id'], center_histories, 'tight_clustered_pair')
    plot_trace(flow_path, run['run_id'], times, np.asarray(detail['flow_concentration_trace'], dtype=float), 'flow concentration', 'Stage C0.10 flow concentration trace')
    plot_trace(grade_path, run['run_id'], times, transfer_signal, 'grade exchange', 'Stage C0.10 grade-exchange trace')
    plot_grades(grades_path, run['run_id'], times, grade_hist)
    plot_paths.extend([traj_path, flow_path, grade_path, grades_path])
    return row, detail, plot_paths


def plot_topology_heatmap(path: Path, rows: list[dict[str, Any]]) -> None:
    label_to_value = {
        'braid_like_exchange': 0,
        'transfer_smeared': 1,
        'mixed_unresolved': 2,
        'dispersive_pass': 3,
    }
    values = np.full((3, 3), np.nan, dtype=float)
    text: dict[tuple[int, int], str] = {}
    ordered = sorted(rows, key=lambda item: (str(item['group_id']), str(item['run_id'])))
    for row in ordered:
        group_order = {'fully_compatible': 0, 'near_compatible': 1, 'incompatible': 2}[str(row['group_id'])]
        variant_order = next(idx for idx, item in enumerate(sorted([r for r in rows if r['group_id'] == row['group_id']], key=lambda r: r['run_id'])) if item['run_id'] == row['run_id'])
        values[group_order, variant_order] = label_to_value[str(row['topology_class'])]
        text[(group_order, variant_order)] = str(row['topology_class']).replace('_', ' ')
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    im = ax.imshow(values, aspect='auto', cmap='viridis', vmin=0.0, vmax=3.0)
    ax.set_yticks([0, 1, 2], ['full', 'near', 'mismatch'])
    ax.set_xticks([0, 1, 2], ['v0', 'v1', 'v2'])
    ax.set_title('Stage C0.10 topology vs address compatibility')
    for (i, j), label in text.items():
        ax.text(j, i, label, ha='center', va='center', color='white', fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_persistence_panel(path: Path, rows: list[dict[str, Any]]) -> None:
    ordered = sorted(rows, key=lambda item: item['run_id'])
    labels = [str(row['run_id']).replace('C010_', '') for row in ordered]
    composite = [float(row['composite_lifetime']) for row in ordered]
    gains = [float(row['persistence_gain']) for row in ordered]
    fig, axes = plt.subplots(2, 1, figsize=(10.8, 7.0), sharex=True)
    axes[0].bar(range(len(ordered)), composite, color='tab:blue')
    axes[0].set_ylabel('composite lifetime')
    axes[0].set_title('Stage C0.10 persistence distribution')
    axes[0].grid(axis='y', alpha=0.25)
    colors = ['tab:green' if value >= 0.0 else 'tab:red' for value in gains]
    axes[1].bar(range(len(ordered)), gains, color=colors)
    axes[1].set_ylabel('persistence gain')
    axes[1].set_xticks(range(len(ordered)), labels, rotation=45, ha='right')
    axes[1].grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_selectivity_panel(path: Path, rows: list[dict[str, Any]]) -> None:
    ordered = sorted(rows, key=lambda item: item['run_id'])
    labels = [str(row['run_id']).replace('C010_', '') for row in ordered]
    selectivity = [float(row['address_selectivity_index']) for row in ordered]
    sorting = [float(row['encounter_sorting_score']) for row in ordered]
    fig, axes = plt.subplots(2, 1, figsize=(10.8, 7.0), sharex=True)
    axes[0].bar(range(len(ordered)), selectivity, color='tab:purple')
    axes[0].set_ylabel('address selectivity')
    axes[0].set_title('Stage C0.10 address selection and sorting')
    axes[0].grid(axis='y', alpha=0.25)
    axes[1].bar(range(len(ordered)), sorting, color='tab:orange')
    axes[1].set_ylabel('sorting score')
    axes[1].set_xticks(range(len(ordered)), labels, rotation=45, ha='right')
    axes[1].grid(axis='y', alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(
    json_path: Path,
    csv_path: Path,
    rows: list[dict[str, Any]],
    stamped_plots: list[str],
    summary: dict[str, Any],
) -> None:
    lines = [
        '# Stage C0.10 Harmonic Address Compatibility Probe v1',
        '',
        f"Timestamped JSON: `{json_path.relative_to(REPO_ROOT)}`",
        f"Timestamped CSV: `{csv_path.relative_to(REPO_ROOT)}`",
        '',
        'Architecture notice: This branch keeps the clustered Phase III DK seed geometry, the projected-transverse sector, the fixed graph resolution, and the Stage 23 topology measurement stack in place. Only a discrete harmonic-address compatibility factor is added as a static relational multiplier on the combinatorial operator.',
        '',
        f"Stage summary: classification counts {summary['classification_counts']}; topology counts {summary['topology_counts']}.",
        '',
        'Per-run summary:',
        '',
    ]
    ordered = sorted(rows, key=lambda item: item['run_id'])
    for row in ordered:
        lines.extend([
            f"- `{row['run_id']}`",
            f"  - group: `{row['group_label']}`",
            f"  - topology class: `{row['topology_class']}`",
            f"  - classification: `{row['classification_label']}`",
            f"  - persistence gain: `{row['persistence_gain']:.4f}`",
            f"  - flow concentration: `{row['flow_concentration_index']:.4f}`",
            f"  - transfer smear index: `{row['transfer_smear_index']:.4f}`",
            f"  - encounter sorting score: `{row['encounter_sorting_score']:.4f}`",
            f"  - address selectivity index: `{row['address_selectivity_index']:.4f}`",
        ])
    lines.extend([
        '',
        'Interpretive note:',
        f"- {summary['interpretive_note']}",
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

    incompatible_rows = [row for row in rows if str(row['group_id']) == 'incompatible']
    incompatible_mean = float(np.mean([float(row['composite_lifetime']) for row in incompatible_rows])) if incompatible_rows else 0.0
    compatible_rows = [row for row in rows if str(row['group_id']) == 'fully_compatible']
    compatible_braid_present = any(int(row['braid_exchange_indicator']) for row in compatible_rows)
    compatible_smear_mean = float(np.mean([float(row['transfer_smear_index']) for row in compatible_rows])) if compatible_rows else 0.0

    for row in rows:
        row['persistence_gain'] = float(row['composite_lifetime']) - incompatible_mean
    for row in rows:
        row['classification_label'] = classify_row(row, compatible_braid_present, compatible_smear_mean)
    for detail in details:
        run_id = str(detail['run']['run_id'])
        detail['metrics'] = next(row for row in rows if str(row['run_id']) == run_id)

    topology_heatmap = WORK_PLOT_DIR / 'stage_c0_10_topology_vs_address_heatmap.png'
    persistence_panel = WORK_PLOT_DIR / 'stage_c0_10_persistence_distribution_panel.png'
    selectivity_panel = WORK_PLOT_DIR / 'stage_c0_10_address_selectivity_panel.png'
    plot_topology_heatmap(topology_heatmap, rows)
    plot_persistence_panel(persistence_panel, rows)
    plot_selectivity_panel(selectivity_panel, rows)
    plot_paths.extend([topology_heatmap, persistence_panel, selectivity_panel])

    classification_counts = dict(Counter(str(row['classification_label']) for row in rows))
    topology_counts = dict(Counter(str(row['topology_class']) for row in rows))
    sorting_positive = sum(float(row['encounter_sorting_score']) >= 0.12 for row in rows)
    if rows:
        persistence_span = float(max(float(row['persistence_gain']) for row in rows)) - float(min(float(row['persistence_gain']) for row in rows))
    else:
        persistence_span = 0.0
    if classification_counts.get('address_protected_braid', 0) >= 1:
        interpretive_note = 'Harmonic addressing behaves like a weak proto-selection rule on this clustered seed: exact-match labeling preserves the braid class, while near-match and incompatible sectors relax to smeared transfer. The visible signal is topology and sorting bias, not measurable persistence gain.'
    elif classification_counts.get('address_sorted_encounter', 0) >= 2:
        interpretive_note = 'Harmonic addressing does not generate a persistence shift on this clustered seed, but it does bias encounter sorting and selective flow concentration in a compatibility-dependent way.'
    else:
        interpretive_note = 'Harmonic addressing remains dynamically weak on this clustered seed: topology and persistence stay close to the incompatible control envelope, so the address rule is better described as a mild texture bias than a strong proto-selection rule.'

    summary = {
        'classification_counts': classification_counts,
        'topology_counts': topology_counts,
        'compatible_braid_present': compatible_braid_present,
        'incompatible_mean_composite_lifetime': incompatible_mean,
        'sorting_positive_runs': sorting_positive,
        'persistence_gain_span': persistence_span,
        'interpretive_note': interpretive_note,
    }
    result = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'common_fields': common,
        'runs': details,
        'summary': summary,
    }
    json_path, csv_path, stamped_plots, _ = save_atlas_payload(
        'stage_c0_10_harmonic_address_compatibility_probe',
        result,
        rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots, summary)


if __name__ == '__main__':
    main()
