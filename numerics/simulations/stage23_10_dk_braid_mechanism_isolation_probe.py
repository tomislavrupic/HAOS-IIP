#!/usr/bin/env python3

from __future__ import annotations

import argparse
import copy
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt, save_atlas_payload
from stage23_1_dk_collision_geometry_scan import (
    field_grid_2d,
    mean_pair_separation_series,
    ordered_peaks,
    plot_trajectories,
)
from stage23_6_dk_phase_sensitivity_braid_exchange import connected_components, edge_grid, grade_exchange_signal
from stage23_9_dk_collision_phenomenology_consolidation import evolve_signed
from stage9b_common import build_dk2d_complex, first_order_dt, pack_positions, periodic_displacement

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage23_10_dk_braid_mechanism_isolation_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_23_10_DK_Braid_Mechanism_Isolation_Probe_v1.md'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage23_10_mechanism_isolation'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'run_id',
    'run_order',
    'family_id',
    'family_label',
    'family_group',
    'role',
    'packet_count',
    'topology_class',
    'base_topology_class',
    'refined_topology_class',
    'base_braid_survival_time',
    'refined_braid_survival_time',
    'flow_concentration_index',
    'grade_transfer_asymmetry',
    'refinement_topology_stability',
    'refinement_topology_match_flag',
    'topology_repeatability_score',
    'non_clustered_family_flag',
    'notes',
]

TOPOLOGY_CODE = {
    'braid_like_exchange': 0,
    'transfer_smeared': 1,
    'dispersive_pass': 2,
    'localized_capture': 3,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage 23.10 DK braid mechanism isolation probe.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


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


def overlap_phase_alignment(state: np.ndarray, packet_states: list[np.ndarray]) -> float:
    phases = [float(np.angle(np.vdot(packet, state))) for packet in packet_states]
    if len(phases) < 2:
        return 0.0
    anchor = phases[0]
    diffs = []
    for phase in phases[1:]:
        diff = (phase - anchor + math.pi) % (2.0 * math.pi) - math.pi
        diffs.append(abs(diff))
    mean_abs = float(np.mean(diffs)) if diffs else 0.0
    return float(np.clip(1.0 - mean_abs / math.pi, 0.0, 1.0))


def flow_concentration(grid: np.ndarray) -> float:
    total = float(np.sum(grid))
    if total <= 1.0e-12:
        return 0.0
    flat = np.sort(grid.ravel())[::-1]
    top_k = max(1, int(math.ceil(0.1 * flat.size)))
    return float(np.sum(flat[:top_k]) / total)


def count_channels(grid: np.ndarray) -> int:
    threshold = max(0.15 * float(np.max(grid)), float(np.mean(grid) + np.std(grid)))
    mask = grid >= threshold
    return connected_components(mask)


def classify_topology(
    packet_count: int,
    braid_flag: bool,
    concentration: float,
    coherence: float,
    grade_amplitude: float,
    recurrence: float,
    close_fraction: float,
) -> str:
    if packet_count == 2 and braid_flag and concentration >= 0.885 and coherence >= 0.40:
        return 'braid_like_exchange'
    if close_fraction >= 0.58 and recurrence >= 0.08 and concentration >= 0.72:
        return 'localized_capture'
    if coherence >= 0.33 and grade_amplitude >= 0.08:
        return 'transfer_smeared'
    return 'dispersive_pass'


def phase_jittered_specs(packet_specs: list[dict[str, Any]], jitter_fraction: float) -> list[dict[str, Any]]:
    specs = copy.deepcopy(packet_specs)
    for idx, spec in enumerate(specs):
        if idx == 0:
            continue
        direction = 1.0 if idx % 2 else -1.0
        spec['phase_offset_fraction_of_pi'] = float(spec['phase_offset_fraction_of_pi']) + direction * float(jitter_fraction)
    return specs


def analyse_states(
    states: list[np.ndarray],
    packet_states: list[np.ndarray],
    positions: np.ndarray,
    edge_midpoints: np.ndarray,
    block_sizes: tuple[int, ...],
    resolution: int,
    anchors: list[np.ndarray],
    close_threshold: float,
) -> dict[str, Any]:
    n0, n1, _ = block_sizes
    packet_count = len(packet_states)
    peak_tracks: list[list[np.ndarray]] = [[] for _ in range(packet_count)]
    previous: list[np.ndarray] | None = None
    edge_grids: list[np.ndarray] = []
    flow_trace: list[float] = []
    align_trace: list[float] = []
    braid_votes: list[float] = []
    grade_hist: list[list[float]] = []

    for step_idx, state in enumerate(states):
        grid = field_grid_2d(positions, state, resolution)
        current, _ = ordered_peaks(grid, packet_count, anchors if step_idx == 0 else None, previous)
        for track, point in zip(peak_tracks, current):
            track.append(point.copy())
        previous = [point.copy() for point in current]
        if packet_count == 2 and len(current) == 2:
            braid_votes.append(float(current[0][0] > current[1][0]))
        edge_values = np.abs(state[n0:n0 + n1])
        egrid = edge_grid(edge_midpoints, edge_values, resolution)
        edge_grids.append(egrid)
        flow_trace.append(flow_concentration(egrid))
        align_trace.append(overlap_phase_alignment(state, packet_states))
        total = float(np.sum(np.abs(state) ** 2)) or 1.0
        grade_hist.append([
            float(np.sum(np.abs(state[:n0]) ** 2) / total),
            float(np.sum(np.abs(state[n0:n0 + n1]) ** 2) / total),
            float(np.sum(np.abs(state[n0 + n1:]) ** 2) / total),
        ])

    center_histories = [np.asarray(track, dtype=float) for track in peak_tracks]
    mean_pair_distances, _ = mean_pair_separation_series(center_histories)
    close_fraction = float(np.mean(mean_pair_distances <= close_threshold)) if mean_pair_distances.size else 0.0

    grade_hist_arr = np.asarray(grade_hist, dtype=float)
    transfer_signal, coherence = grade_exchange_signal(grade_hist_arr)
    grade_amplitude = float(np.max(transfer_signal)) if transfer_signal.size else 0.0
    recurrence = 0.0
    if mean_pair_distances.size >= 6:
        tail = mean_pair_distances[len(mean_pair_distances) // 3:]
        deriv = np.diff(tail)
        signs = np.sign(deriv)
        signs = signs[signs != 0.0]
        if signs.size >= 2:
            turns = np.sum(signs[1:] * signs[:-1] < 0.0)
            recurrence = float((turns / max(1, signs.size - 1)) * close_fraction)

    braid_flag = abs(float(np.mean(braid_votes)) - 0.5) > 0.25 if braid_votes else False
    avg_edge_grid = np.mean(np.asarray(edge_grids, dtype=float), axis=0)
    concentration = flow_concentration(avg_edge_grid)
    topology_class = classify_topology(
        packet_count=packet_count,
        braid_flag=braid_flag,
        concentration=concentration,
        coherence=coherence,
        grade_amplitude=grade_amplitude,
        recurrence=recurrence,
        close_fraction=close_fraction,
    )

    instant_labels: list[str] = []
    max_grade = max(grade_amplitude, 1.0e-12)
    for idx, egrid in enumerate(edge_grids):
        local_conc = flow_concentration(egrid)
        local_close = 1.0 if idx < len(mean_pair_distances) and float(mean_pair_distances[idx]) <= close_threshold else 0.0
        local_coh = float(np.clip(float(transfer_signal[idx]) / max_grade, 0.0, 1.0)) if idx < len(transfer_signal) else coherence
        local_braid = braid_flag and packet_count == 2 and idx < len(braid_votes) and abs(braid_votes[idx] - 0.5) > 0.25
        instant_labels.append(
            classify_topology(
                packet_count=packet_count,
                braid_flag=local_braid,
                concentration=local_conc,
                coherence=local_coh,
                grade_amplitude=grade_amplitude,
                recurrence=recurrence,
                close_fraction=local_close,
            )
        )

    grade_transfer_asymmetry = abs(float(grade_hist_arr[-1, 1]) - float(grade_hist_arr[-1, 0])) if grade_hist_arr.size else 0.0
    braid_survival_time = float(sum(label == 'braid_like_exchange' for label in instant_labels))
    return {
        'topology_class': topology_class,
        'braid_survival_time': braid_survival_time,
        'flow_concentration_index': concentration,
        'grade_transfer_asymmetry': grade_transfer_asymmetry,
        'grade_exchange_coherence': coherence,
        'topology_trace': instant_labels,
        'flow_trace': flow_trace,
        'alignment_trace': align_trace,
        'grade_exchange_trace': transfer_signal.tolist(),
        'center_histories': [history.tolist() for history in center_histories],
        'channel_count': count_channels(avg_edge_grid),
    }


def simulate_resolution(
    packet_specs: list[dict[str, Any]],
    family_label: str,
    run_id: str,
    resolution: int,
    common: dict[str, Any],
    trajectory_suffix: str | None = None,
) -> tuple[dict[str, Any], dict[str, Any], list[Path]]:
    epsilon = float(common['base_epsilon'])
    kick_cycles = 1.0
    beta = float(common['beta'])
    t_final = float(common['t_final'])

    complex_data = build_dk2d_complex(n_side=resolution, epsilon=epsilon)
    positions = pack_positions(complex_data.points, complex_data.edge_midpoints, complex_data.face_centers)
    block_sizes = complex_data.block_sizes
    anchors = [np.asarray(spec['center'], dtype=float) for spec in packet_specs]
    sigmas = [float(spec['sigma']) for spec in packet_specs]
    close_threshold = 2.0 * float(np.mean(sigmas))

    packet_states = []
    for spec in packet_specs:
        packet_states.append(
            packet_state_with_profile(
                positions=positions,
                block_sizes=block_sizes,
                center=np.asarray(spec['center'], dtype=float),
                sigma=float(spec['sigma']),
                amplitude=float(spec['amplitude']),
                phase_offset=float(spec['phase_offset_fraction_of_pi']) * math.pi,
                kick_vector=np.asarray(spec['kick'], dtype=float),
                kick_cycles=kick_cycles,
                grade0_scale=float(spec['grade0_scale']),
                grade1_scale=float(spec['grade1_scale']),
                anisotropy=(float(spec['anisotropy_x']), float(spec['anisotropy_y'])),
            )
        )

    psi0 = np.sum(packet_states, axis=0)
    dt = first_order_dt(complex_data.dirac_kahler, float(common['dt_scale']))
    steps = max(2, int(math.ceil(t_final / dt)))
    states = evolve_signed(complex_data.dirac_kahler, block_sizes, psi0, dt, steps, beta, direction=1.0)
    metrics = analyse_states(
        states=states,
        packet_states=packet_states,
        positions=positions,
        edge_midpoints=complex_data.edge_midpoints,
        block_sizes=block_sizes,
        resolution=resolution,
        anchors=anchors,
        close_threshold=close_threshold,
    )

    plot_paths: list[Path] = []
    if trajectory_suffix is not None:
        trajectory_path = WORK_PLOT_DIR / f'{run_id}_{trajectory_suffix}_trajectory.png'
        center_histories = [np.asarray(history, dtype=float) for history in metrics['center_histories']]
        plot_trajectories(trajectory_path, f'{run_id}_{trajectory_suffix}', center_histories, family_label.replace(' ', '_').lower())
        plot_paths.append(trajectory_path)

    detail = {
        'resolution': int(resolution),
        'epsilon': epsilon,
        'dt': dt,
        'steps': steps,
        'metrics': metrics,
    }
    return metrics, detail, plot_paths


def refinement_stability(base_metrics: dict[str, Any], refined_metrics: dict[str, Any]) -> tuple[float, int]:
    same_topology = int(str(base_metrics['topology_class']) == str(refined_metrics['topology_class']))
    base_survival = float(base_metrics['braid_survival_time'])
    refined_survival = float(refined_metrics['braid_survival_time'])
    if base_survival <= 1.0e-12 and refined_survival <= 1.0e-12:
        return (1.0 if same_topology else 9.99), same_topology
    if base_survival <= 1.0e-12 or refined_survival <= 1.0e-12:
        return 9.99, same_topology
    ratio = max(base_survival, refined_survival) / max(min(base_survival, refined_survival), 1.0e-12)
    if not same_topology:
        ratio = max(ratio, 1.5)
    return float(ratio), same_topology


def repeatability_score(
    run: dict[str, Any],
    common: dict[str, Any],
    reference_topology: str,
) -> tuple[float, list[dict[str, Any]]]:
    trials: list[dict[str, Any]] = []
    hits = 0
    for jitter in common['repeatability_phase_jitter_fraction_of_pi']:
        jittered = phase_jittered_specs(run['packet_specs'], float(jitter))
        metrics, detail, _ = simulate_resolution(
            packet_specs=jittered,
            family_label=run['family_label'],
            run_id=run['run_id'],
            resolution=int(common['base_resolution']),
            common=common,
            trajectory_suffix=None,
        )
        match = int(str(metrics['topology_class']) == str(reference_topology))
        hits += match
        trials.append(
            {
                'phase_jitter_fraction_of_pi': float(jitter),
                'topology_class': metrics['topology_class'],
                'braid_survival_time': float(metrics['braid_survival_time']) * float(detail['dt']),
                'flow_concentration_index': float(metrics['flow_concentration_index']),
                'match_reference': match,
            }
        )
    score = float(hits / max(1, len(trials)))
    return score, trials


def classify_mechanism(rows: list[dict[str, Any]]) -> tuple[str, dict[str, Any]]:
    non_clustered = [row for row in rows if int(row['non_clustered_family_flag']) == 1]
    robust_non_clustered = [
        row for row in non_clustered
        if str(row['topology_class']) == 'braid_like_exchange'
        and float(row['refinement_topology_stability']) <= 1.1
        and float(row['topology_repeatability_score']) >= 0.7
    ]
    clustered_braid = [
        row for row in rows
        if row['family_group'] == 'clustered' and str(row['topology_class']) == 'braid_like_exchange'
    ]
    clustered_degraded = next(row for row in rows if row['family_id'] == 'clustered_degraded')
    clustered_only = bool(clustered_braid) and all(str(row['topology_class']) != 'braid_like_exchange' for row in non_clustered)
    degraded_collapse = str(clustered_degraded['topology_class']) != 'braid_like_exchange'

    if len(robust_non_clustered) >= 2:
        classification = 'intrinsic_exchange_mechanism'
    elif clustered_only or degraded_collapse:
        classification = 'clustered_texture_artefact'
    else:
        classification = 'mechanism_boundary_inconclusive'

    payload = {
        'robust_non_clustered_braid_hits': int(len(robust_non_clustered)),
        'non_clustered_braid_hits': int(sum(str(row['topology_class']) == 'braid_like_exchange' for row in non_clustered)),
        'clustered_braid_hits': int(len(clustered_braid)),
        'clustered_degraded_collapse': bool(degraded_collapse),
    }
    return classification, payload


def plot_topology_survival_matrix(path: Path, rows: list[dict[str, Any]]) -> None:
    ordered = sorted(rows, key=lambda row: int(row['run_order']))
    grid = np.zeros((len(ordered), 2), dtype=float)
    for idx, row in enumerate(ordered):
        grid[idx, 0] = TOPOLOGY_CODE[str(row['base_topology_class'])]
        grid[idx, 1] = TOPOLOGY_CODE[str(row['refined_topology_class'])]

    fig, ax = plt.subplots(figsize=(7.6, 5.8))
    im = ax.imshow(grid, aspect='auto', cmap='viridis', vmin=0.0, vmax=3.0)
    ax.set_xticks([0, 1], ['n', '2n'])
    ax.set_yticks(range(len(ordered)), [str(row['run_id']).replace('S23_10_', '') for row in ordered])
    ax.set_title('Stage 23.10 topology survival matrix')
    for idx, row in enumerate(ordered):
        ax.text(1.15, idx, f"rep={float(row['topology_repeatability_score']):.2f}", va='center', fontsize=8)
        ax.text(1.75, idx, f"stab={float(row['refinement_topology_stability']):.2f}", va='center', fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_family_comparison_panel(path: Path, rows: list[dict[str, Any]]) -> None:
    ordered = sorted(rows, key=lambda row: int(row['run_order']))
    labels = [str(row['run_id']).replace('S23_10_', '') for row in ordered]
    base_vals = [float(row['base_braid_survival_time']) for row in ordered]
    refined_vals = [float(row['refined_braid_survival_time']) for row in ordered]
    x = np.arange(len(ordered), dtype=float)
    width = 0.38
    fig, ax = plt.subplots(figsize=(10.8, 4.8))
    ax.bar(x - 0.5 * width, base_vals, width=width, color='tab:blue', label='n')
    ax.bar(x + 0.5 * width, refined_vals, width=width, color='tab:orange', label='2n')
    ax.set_ylabel('braid survival time')
    ax.set_xticks(x, labels, rotation=45, ha='right')
    ax.set_title('Stage 23.10 family-comparison persistence panel')
    ax.grid(axis='y', alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(
    json_path: Path,
    csv_path: Path,
    rows: list[dict[str, Any]],
    stamped_plots: list[str],
    classification: str,
    classification_payload: dict[str, Any],
) -> None:
    ordered = sorted(rows, key=lambda row: int(row['run_order']))
    lines = [
        '# Stage 23.10 DK Braid Mechanism Isolation Probe v1',
        '',
        f"Timestamped JSON: `{json_path.relative_to(REPO_ROOT)}`",
        f"Timestamped CSV: `{csv_path.relative_to(REPO_ROOT)}`",
        '',
        f"Mechanism classification: `{classification}`",
        '',
        '## Per-run readout',
    ]
    for row in ordered:
        lines.append(
            f"- `{row['run_id']}`: base=`{row['base_topology_class']}`, refined=`{row['refined_topology_class']}`, "
            f"repeatability=`{float(row['topology_repeatability_score']):.2f}`, "
            f"stability=`{float(row['refinement_topology_stability']):.2f}`"
        )

    lines.extend([
        '',
        '## Decision summary',
        f"- robust non-clustered braid hits: `{classification_payload['robust_non_clustered_braid_hits']}`",
        f"- non-clustered braid hits (raw): `{classification_payload['non_clustered_braid_hits']}`",
        f"- clustered braid hits: `{classification_payload['clustered_braid_hits']}`",
        f"- clustered degraded collapse: `{classification_payload['clustered_degraded_collapse']}`",
        '',
        '## Interpretation',
    ])

    if classification == 'intrinsic_exchange_mechanism':
        lines.append('Braid-like exchange is no longer confined to the clustered seed family: at least two non-clustered families recover the braid label with acceptable refinement stability and phase-jitter repeatability, so the exchange pattern is classified as an intrinsic local mechanism of the frozen DK branch.')
    elif classification == 'clustered_texture_artefact':
        lines.append('The braid window remains clustered-texture dependent in this pass: non-clustered families do not recover a robust braid sector, and the clustered family either stays isolated or collapses under mild degradation. Phase III therefore closes this line as a seed-dependent texture rather than a general exchange law.')
    else:
        lines.append('The present nine-run lattice narrows the mechanism boundary but does not fully close it. The braid sector is neither reproduced broadly enough to classify as intrinsic nor cleanly confined enough to label as a pure clustered artefact.')

    lines.extend([
        '',
        '## Assumption note',
        'Because the frozen DK lattice here is regular, the requested degree-spectrum perturbation was represented by a mild localized graded-support skew surrogate rather than a literal graph-degree modification. The kernel class, grading rule, timestep policy, and normalization were otherwise held fixed.',
        '',
        '## Plots',
    ])
    for plot in stamped_plots:
        lines.append(f'- `{plot}`')

    NOTE_PATH.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main() -> None:
    args = parse_args()
    runsheet = load_runsheet(args.runsheet)
    common = runsheet['common_fields']

    rows: list[dict[str, Any]] = []
    details: list[dict[str, Any]] = []
    plot_paths: list[Path] = []

    for run in sorted(runsheet['runs'], key=lambda item: int(item['run_order'])):
        base_metrics, base_detail, base_plots = simulate_resolution(
            packet_specs=run['packet_specs'],
            family_label=run['family_label'],
            run_id=run['run_id'],
            resolution=int(common['base_resolution']),
            common=common,
            trajectory_suffix='n12',
        )
        refined_metrics, refined_detail, refined_plots = simulate_resolution(
            packet_specs=run['packet_specs'],
            family_label=run['family_label'],
            run_id=run['run_id'],
            resolution=int(common['refined_resolution']),
            common=common,
            trajectory_suffix='n24',
        )
        stability_ratio, match_flag = refinement_stability(base_metrics, refined_metrics)
        repeat_score, repeat_trials = repeatability_score(run, common, str(base_metrics['topology_class']))

        row = {
            'run_id': run['run_id'],
            'run_order': int(run['run_order']),
            'family_id': run['family_id'],
            'family_label': run['family_label'],
            'family_group': run['family_group'],
            'role': run['role'],
            'packet_count': len(run['packet_specs']),
            'topology_class': str(base_metrics['topology_class']),
            'base_topology_class': str(base_metrics['topology_class']),
            'refined_topology_class': str(refined_metrics['topology_class']),
            'base_braid_survival_time': float(base_metrics['braid_survival_time']) * float(base_detail['dt']),
            'refined_braid_survival_time': float(refined_metrics['braid_survival_time']) * float(refined_detail['dt']),
            'flow_concentration_index': float(base_metrics['flow_concentration_index']),
            'grade_transfer_asymmetry': float(base_metrics['grade_transfer_asymmetry']),
            'refinement_topology_stability': float(stability_ratio),
            'refinement_topology_match_flag': int(match_flag),
            'topology_repeatability_score': float(repeat_score),
            'non_clustered_family_flag': int(run['family_group'] != 'clustered'),
            'notes': run['notes'],
        }
        detail = {
            'run_id': run['run_id'],
            'family_id': run['family_id'],
            'family_label': run['family_label'],
            'family_group': run['family_group'],
            'notes': run['notes'],
            'base_resolution': base_detail,
            'refined_resolution': refined_detail,
            'repeatability_trials': repeat_trials,
        }
        rows.append(row)
        details.append(detail)
        plot_paths.extend(base_plots + refined_plots)

    classification, classification_payload = classify_mechanism(rows)
    topology_matrix = WORK_PLOT_DIR / 'stage23_10_topology_survival_matrix.png'
    persistence_panel = WORK_PLOT_DIR / 'stage23_10_family_comparison_persistence_panel.png'
    plot_topology_survival_matrix(topology_matrix, rows)
    plot_family_comparison_panel(persistence_panel, rows)
    plot_paths.extend([topology_matrix, persistence_panel])

    result = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'common_fields': common,
        'runs': details,
        'summary': {
            'classification': classification,
            'classification_payload': classification_payload,
            'topology_counts': dict(Counter(str(row['topology_class']) for row in rows)),
        },
    }
    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        'stage23_10_dk_braid_mechanism_isolation_probe',
        result,
        rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots, classification, classification_payload)


if __name__ == '__main__':
    main()
