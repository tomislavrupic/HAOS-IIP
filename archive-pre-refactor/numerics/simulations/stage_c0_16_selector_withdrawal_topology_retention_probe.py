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
from stage23_1_dk_collision_geometry_scan import CSV_FIELDS as BASE_FIELDS, WORK_PLOT_DIR, plot_trace
from stage23_6_dk_phase_sensitivity_braid_exchange import evolve, packet_state_with_skew
from stage9b_common import build_dk2d_complex, first_order_dt, pack_positions
from stage_c0_12_harmonic_detuning_continuum_scan import address_weighted_operator
from stage_c0_14_topology_locking_mechanism_probe import (
    apply_node_bias,
    build_family_labels_with_mode,
    identify_motif_nodes,
    make_rng,
    measure_window,
)

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_16_selector_withdrawal_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_C0_16_Selector_Withdrawal_Topology_Retention_Probe_v1.md'

CSV_FIELDS = BASE_FIELDS + [
    'run_order',
    'family_id',
    'family_label',
    'family_mode',
    'withdrawal_protocol',
    'withdrawal_label',
    'arming_delta',
    'armed_topology_class',
    'post_withdrawal_topology_class',
    'withdrawal_braid_dwell_time',
    'withdrawal_class_retention',
    'post_withdrawal_selectivity',
    'post_withdrawal_coherence',
    'post_withdrawal_flow_concentration',
    'selector_off_recurrence',
    'locking_after_withdrawal_label',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.16 selector-withdrawal topology retention probe.')
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


def topology_code(label: str) -> float:
    return {
        'braid_like_exchange': 0.0,
        'transfer_smeared': 1.0,
        'unresolved': 2.0,
    }.get(label, 3.0)


def make_initial_packets(positions: np.ndarray, block_sizes: tuple[int, ...], common: dict[str, Any]) -> list[np.ndarray]:
    amplitude = float(common['amplitude'])
    sigma = float(common['mean_width'])
    phase_offset = float(common['phase_offset_fraction_of_pi']) * math.pi
    skew = float(common['local_phase_skew_fraction'])
    kick_cycles = float(common['kick_cycles'])
    separation = float(common['separation'])
    centers = [
        np.array([0.50 - 0.5 * separation, 0.50], dtype=float),
        np.array([0.50 + 0.5 * separation, 0.50], dtype=float),
    ]
    kicks = [np.array([0.4, 0.0], dtype=float), np.array([-0.4, 0.0], dtype=float)]
    return [
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


def majority_label(labels: list[str]) -> str:
    if not labels:
        return 'unresolved'
    return Counter(labels).most_common(1)[0][0]


def classify_withdrawal_label(
    armed_topology: str,
    post_topology: str,
    braid_dwell_fraction: float,
    class_retention: float,
    withdrawal_coherence: float,
    withdrawal_flow: float,
) -> str:
    if (
        armed_topology == 'braid_like_exchange'
        and post_topology == 'braid_like_exchange'
        and braid_dwell_fraction >= 0.67
        and class_retention >= 0.67
        and withdrawal_coherence >= 0.42
        and withdrawal_flow >= 0.88
    ):
        return 'locked_after_withdrawal'
    if (
        armed_topology == 'braid_like_exchange'
        and post_topology == 'braid_like_exchange'
        and braid_dwell_fraction >= 0.34
        and class_retention >= 0.50
    ):
        return 'quasi_locked_retention'
    if armed_topology == 'braid_like_exchange' and braid_dwell_fraction > 0.0:
        return 'topology_afterglow'
    if armed_topology == 'braid_like_exchange' and post_topology != 'braid_like_exchange':
        return 'selector_dependent_collapse'
    return 'no_locking_after_withdrawal'


def simulate_run(run: dict[str, Any], common: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], list[Path]]:
    resolution = int(common['resolution'])
    epsilon = float(common['epsilon'])
    modulus = float(common['address_modulus'])
    eta_match = float(common['eta_match'])
    eta_mismatch = float(common['eta_mismatch'])
    support_floor_fraction = float(common['support_floor_fraction'])
    dominant_ratio_threshold = float(common['dominant_ratio_threshold'])
    beta = float(common['beta'])
    segment_steps = int(common['segment_steps'])
    arming_segments = int(common['arming_segments'])
    withdrawal_segments = int(common['withdrawal_segments'])
    dt_scale = float(common['dt_scale'])

    complex_data = build_dk2d_complex(n_side=resolution, epsilon=epsilon)
    positions = pack_positions(complex_data.points, complex_data.edge_midpoints, complex_data.face_centers)
    block_sizes = complex_data.block_sizes
    packet_states = make_initial_packets(positions, block_sizes, common)
    base_operator = complex_data.dirac_kahler
    dt = first_order_dt(base_operator, dt_scale)
    sigma = float(common['mean_width'])
    rng = make_rng(str(run['run_id']))
    labels = build_family_labels_with_mode(
        packet_states=packet_states,
        block_sizes=block_sizes,
        delta=float(run['arming_delta']),
        family_mode=str(run['family_mode']),
        support_floor_fraction=support_floor_fraction,
        dominant_ratio_threshold=dominant_ratio_threshold,
        modulus=modulus,
        rng=rng,
    )
    selected_operator = address_weighted_operator(base_operator, labels, modulus, eta_match, eta_mismatch)
    motif_nodes = identify_motif_nodes(base_operator, packet_states, common)
    motif_bias = np.zeros(base_operator.shape[0], dtype=float)
    motif_bias[motif_nodes] = 1.0
    degree_bias = np.zeros(base_operator.shape[0], dtype=float)
    degree_bias[: max(4, base_operator.shape[0] // 16)] = 1.0

    psi = np.sum(packet_states, axis=0)
    topology_trace: list[str] = []
    topology_codes: list[float] = []
    selectivity_trace: list[float] = []
    flow_trace: list[float] = []
    coherence_trace: list[float] = []
    segment_times: list[float] = []
    windows: list[dict[str, Any]] = []

    cumulative_steps = 0
    total_segments = arming_segments + withdrawal_segments
    for seg_idx in range(total_segments):
        if seg_idx < arming_segments:
            operator = selected_operator
        else:
            operator = base_operator
            protocol = str(run['withdrawal_protocol'])
            if protocol == 'degree_skew_withdrawal':
                operator = apply_node_bias(operator, degree_bias, float(common['degree_bias_amplitude']))
            elif protocol == 'motif_withdrawal':
                operator = apply_node_bias(operator, motif_bias, float(common['motif_bias_amplitude']))

        states = evolve(operator, block_sizes, psi, dt, segment_steps, beta)
        psi = states[-1]
        window = measure_window(
            states=states,
            packet_states=packet_states,
            positions=positions,
            edge_midpoints=complex_data.edge_midpoints,
            block_sizes=block_sizes,
            resolution=resolution,
            sigma=sigma,
            labels=labels,
            operator=operator,
            modulus=modulus,
            eta_match=eta_match,
            eta_mismatch=eta_mismatch,
        )
        cumulative_steps += segment_steps
        segment_times.append(cumulative_steps * dt)
        windows.append(window)
        label = str(window['topology_class'])
        topology_trace.append(label)
        topology_codes.append(topology_code(label))
        selectivity_trace.append(float(window['address_selectivity_index']))
        flow_trace.append(float(window['flow_concentration_index']))
        coherence_trace.append(float(window['exchange_coherence']))

    armed_windows = windows[:arming_segments]
    withdrawal_windows = windows[arming_segments:]
    armed_topology = majority_label([str(window['topology_class']) for window in armed_windows])
    post_topology = majority_label([str(window['topology_class']) for window in withdrawal_windows])
    withdrawal_braid_dwell_segments = sum(1 for window in withdrawal_windows if str(window['topology_class']) == 'braid_like_exchange')
    withdrawal_braid_dwell_time = float(withdrawal_braid_dwell_segments * segment_steps * dt)
    withdrawal_class_retention = float(
        np.mean([str(window['topology_class']) == armed_topology for window in withdrawal_windows])
    ) if withdrawal_windows else 0.0
    withdrawal_selectivity = float(np.mean([float(window['address_selectivity_index']) for window in withdrawal_windows])) if withdrawal_windows else 0.0
    withdrawal_coherence = float(np.mean([float(window['exchange_coherence']) for window in withdrawal_windows])) if withdrawal_windows else 0.0
    withdrawal_flow = float(np.mean([float(window['flow_concentration_index']) for window in withdrawal_windows])) if withdrawal_windows else 0.0
    withdrawal_recurrence = float(np.mean([float(window['recurrence_indicator']) for window in withdrawal_windows])) if withdrawal_windows else 0.0
    braid_dwell_fraction = float(withdrawal_braid_dwell_segments / max(len(withdrawal_windows), 1))
    withdrawal_label = classify_withdrawal_label(
        armed_topology=armed_topology,
        post_topology=post_topology,
        braid_dwell_fraction=braid_dwell_fraction,
        class_retention=withdrawal_class_retention,
        withdrawal_coherence=withdrawal_coherence,
        withdrawal_flow=withdrawal_flow,
    )

    row = {
        'run_id': str(run['run_id']),
        'geometry_id': 'tight_clustered_pair',
        'phase_id': str(run['withdrawal_protocol']),
        'resolution': resolution,
        'graph_type': common['graph_type'],
        'packet_count': 2,
        'initial_peak_count': int(armed_windows[0]['raw_peak_count']) if armed_windows else 0,
        'collision_label': post_topology,
        'persistence_label': withdrawal_label,
        'composite_lifetime': withdrawal_braid_dwell_time,
        'binding_persistence': braid_dwell_fraction,
        'corridor_dwell': 0.0,
        'coarse_basin_persistence': 0.0,
        'minimum_separation': float(np.min(withdrawal_windows[-1]['mean_pair_distances'])) if withdrawal_windows and withdrawal_windows[-1]['mean_pair_distances'].size else 0.0,
        'final_mean_separation': float(withdrawal_windows[-1]['mean_pair_distances'][-1]) if withdrawal_windows and withdrawal_windows[-1]['mean_pair_distances'].size else 0.0,
        'post_collision_separation_trend': float(withdrawal_windows[-1]['transport_span']) if withdrawal_windows else 0.0,
        'encounter_dwell_time': withdrawal_braid_dwell_time,
        'deflection_angle_proxy': 0.0,
        'reflection_fraction': 0.0,
        'grade_transfer_amplitude': float(np.max(withdrawal_windows[-1]['transfer_signal'])) if withdrawal_windows and withdrawal_windows[-1]['transfer_signal'].size else 0.0,
        'omega0_weight_initial': 0.0,
        'omega1_weight_initial': 0.0,
        'omega2_weight_initial': 0.0,
        'omega0_weight_final': 0.0,
        'omega1_weight_final': 0.0,
        'omega2_weight_final': 0.0,
        'localized_circulation_proxy': 0.0,
        'new_collision_class': int(withdrawal_label in {'locked_after_withdrawal', 'quasi_locked_retention'}),
        'gate_met': int(withdrawal_label in {'locked_after_withdrawal', 'quasi_locked_retention'}),
        'promoted_followup': 0,
        'run_order': int(run['run_order']),
        'family_id': str(run['family_id']),
        'family_label': str(run['family_label']),
        'family_mode': str(run['family_mode']),
        'withdrawal_protocol': str(run['withdrawal_protocol']),
        'withdrawal_label': str(run['withdrawal_label']),
        'arming_delta': float(run['arming_delta']),
        'armed_topology_class': armed_topology,
        'post_withdrawal_topology_class': post_topology,
        'withdrawal_braid_dwell_time': withdrawal_braid_dwell_time,
        'withdrawal_class_retention': withdrawal_class_retention,
        'post_withdrawal_selectivity': withdrawal_selectivity,
        'post_withdrawal_coherence': withdrawal_coherence,
        'post_withdrawal_flow_concentration': withdrawal_flow,
        'selector_off_recurrence': withdrawal_recurrence,
        'locking_after_withdrawal_label': withdrawal_label,
        'notes': str(run['notes']),
    }

    detail = {
        'run': run,
        'metrics': row,
        'segment_times': segment_times,
        'topology_trace': topology_trace,
        'selectivity_trace': selectivity_trace,
        'flow_trace': flow_trace,
        'coherence_trace': coherence_trace,
        'phase_split': {
            'arming_segments': arming_segments,
            'withdrawal_segments': withdrawal_segments,
        },
        'windows': [
            {
                'phase': 'arming' if idx < arming_segments else 'withdrawal',
                'topology_class': str(window['topology_class']),
                'flow_concentration_index': float(window['flow_concentration_index']),
                'address_selectivity_index': float(window['address_selectivity_index']),
                'exchange_coherence': float(window['exchange_coherence']),
            }
            for idx, window in enumerate(windows)
        ],
    }

    plot_paths: list[Path] = []
    topology_path = WORK_PLOT_DIR / f"stage_c0_16_{run['run_id']}_topology_trace.png"
    selectivity_path = WORK_PLOT_DIR / f"stage_c0_16_{run['run_id']}_selector_off_selectivity_trace.png"
    plot_trace(
        topology_path,
        str(run['run_id']),
        segment_times,
        np.asarray(topology_codes, dtype=float),
        'topology code',
        'Stage C0.16 topology trace',
    )
    plot_trace(
        selectivity_path,
        str(run['run_id']),
        segment_times,
        np.asarray(selectivity_trace, dtype=float),
        'address selectivity',
        'Stage C0.16 selector-off selectivity trace',
    )
    plot_paths.extend([topology_path, selectivity_path])
    return row, detail, plot_paths


def plot_withdrawal_dwell_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [str(row['run_id']) for row in rows]
    x = np.arange(len(rows))
    width = 0.35
    dwell = np.asarray([float(row['withdrawal_braid_dwell_time']) for row in rows], dtype=float)
    retention = np.asarray([float(row['withdrawal_class_retention']) for row in rows], dtype=float)
    fig, ax = plt.subplots(figsize=(14, 5.2))
    ax.bar(x - width / 2.0, dwell, width=width, label='withdrawal braid dwell')
    ax.bar(x + width / 2.0, retention, width=width, label='class retention')
    ax.set_xticks(x, labels, rotation=45, ha='right')
    ax.set_ylabel('dwell / retention')
    ax.set_title('Stage C0.16 withdrawal dwell summary')
    ax.grid(axis='y', alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_class_matrix(path: Path, rows: list[dict[str, Any]]) -> None:
    families = sorted({str(row['family_label']) for row in rows})
    protocols = sorted({str(row['withdrawal_label']) for row in rows})
    code_map = {
        'no_locking_after_withdrawal': 0.0,
        'selector_dependent_collapse': 1.0,
        'topology_afterglow': 2.0,
        'quasi_locked_retention': 3.0,
        'locked_after_withdrawal': 4.0,
    }
    text_map = {(family, protocol): '' for family in families for protocol in protocols}
    grid = np.zeros((len(families), len(protocols)), dtype=float)
    for row in rows:
        i = families.index(str(row['family_label']))
        j = protocols.index(str(row['withdrawal_label']))
        label = str(row['locking_after_withdrawal_label'])
        grid[i, j] = code_map[label]
        text_map[(str(row['family_label']), str(row['withdrawal_label']))] = label.replace('_', ' ')
    fig, ax = plt.subplots(figsize=(10.6, 4.4))
    im = ax.imshow(grid, aspect='auto', cmap='Blues', vmin=0.0, vmax=4.0)
    ax.set_xticks(np.arange(len(protocols)), protocols, rotation=20, ha='right')
    ax.set_yticks(np.arange(len(families)), families)
    ax.set_title('Stage C0.16 locking-after-withdrawal matrix')
    for i, family in enumerate(families):
        for j, protocol in enumerate(protocols):
            ax.text(j, i, text_map[(family, protocol)], ha='center', va='center', fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def global_read(rows: list[dict[str, Any]]) -> str:
    locked = [row for row in rows if str(row['locking_after_withdrawal_label']) in {'locked_after_withdrawal', 'quasi_locked_retention'}]
    afterglow = [row for row in rows if str(row['locking_after_withdrawal_label']) == 'topology_afterglow']
    if len(locked) >= 3:
        return 'selector withdrawal leaves a real retained topology signal in several runs, which improves the locking case'
    if afterglow:
        return 'selector withdrawal produces only short-lived afterglow, not stable locking'
    return 'selected topology collapses once the selector is removed, so the locking case does not improve here'


def write_note(json_path: Path, csv_path: Path, rows: list[dict[str, Any]], stamped_plots: list[str], verdict: str) -> None:
    lines = [
        '# Stage C0.16 Selector Withdrawal Topology Retention Probe v1',
        '',
        f"Timestamped JSON: `{json_path.relative_to(REPO_ROOT)}`",
        f"Timestamped CSV: `{csv_path.relative_to(REPO_ROOT)}`",
        '',
        f'Global read: `{verdict}`',
        '',
        '## Per-run summary',
    ]
    for row in rows:
        lines.append(
            f"- `{row['run_id']}`: family=`{row['family_label']}`, withdrawal=`{row['withdrawal_label']}`, armed=`{row['armed_topology_class']}`, post=`{row['post_withdrawal_topology_class']}`, dwell=`{row['withdrawal_braid_dwell_time']:.4f}`, label=`{row['locking_after_withdrawal_label']}`"
        )
    lines.extend([
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

    rows.sort(key=lambda item: int(item['run_order']))
    details.sort(key=lambda item: int(item['metrics']['run_order']))

    dwell_path = WORK_PLOT_DIR / 'stage_c0_16_withdrawal_dwell_summary.png'
    matrix_path = WORK_PLOT_DIR / 'stage_c0_16_locking_after_withdrawal_matrix.png'
    plot_withdrawal_dwell_summary(dwell_path, rows)
    plot_class_matrix(matrix_path, rows)
    plot_paths.extend([dwell_path, matrix_path])

    verdict = global_read(rows)
    summary = {
        'withdrawal_labels': dict(Counter(str(row['locking_after_withdrawal_label']) for row in rows)),
        'post_topology_counts': dict(Counter(str(row['post_withdrawal_topology_class']) for row in rows)),
        'verdict': verdict,
    }
    result = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'common_fields': common,
        'runs': details,
        'summary': summary,
    }
    json_path, csv_path, stamped_plots, _ = save_atlas_payload(
        'stage_c0_16_selector_withdrawal_topology_retention_probe',
        result,
        rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots, verdict)


if __name__ == '__main__':
    main()
