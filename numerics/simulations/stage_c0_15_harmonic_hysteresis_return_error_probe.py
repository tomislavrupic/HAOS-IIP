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
from stage_c0_11_harmonic_address_resonance_ladder import combined_initial_weights
from stage_c0_12_harmonic_detuning_continuum_scan import address_weighted_operator
from stage_c0_14_topology_locking_mechanism_probe import (
    apply_node_bias,
    build_family_labels_with_mode,
    identify_motif_nodes,
    make_rng,
    measure_window,
)

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_15_harmonic_hysteresis_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_C0_15_Harmonic_Hysteresis_Return_Error_Probe_v1.md'

CSV_FIELDS = BASE_FIELDS + [
    'run_order',
    'family_id',
    'family_label',
    'family_mode',
    'protocol',
    'protocol_label',
    'initial_delta',
    'final_delta',
    'initial_topology_class',
    'final_topology_class',
    'topology_return_error',
    'flow_return_error',
    'selectivity_return_error',
    'coherence_return_error',
    'hysteresis_area_flow_concentration',
    'hysteresis_area_selectivity',
    'hysteresis_area_coherence',
    'irreversible_witness',
    'retrace_class',
    'peak_hysteresis_metric',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.15 harmonic hysteresis return-error probe.')
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


def choose_skew_nodes(packet_states: list[np.ndarray], top_count: int = 12) -> np.ndarray:
    weights = combined_initial_weights(packet_states)
    order = np.argsort(weights)[::-1]
    chosen = order[:top_count]
    mask = np.zeros_like(weights, dtype=float)
    mask[chosen] = 1.0
    return mask


def hysteresis_area(deltas: list[float], metric: list[float]) -> float:
    if len(metric) < 3:
        return 0.0
    mid = len(metric) // 2
    forward_y = np.asarray(metric[:mid], dtype=float)
    reverse_y = np.asarray(metric[mid + 1 :], dtype=float)[::-1]
    forward_x = np.asarray(deltas[:mid], dtype=float)
    reverse_x = np.asarray(deltas[mid + 1 :], dtype=float)[::-1]
    n = min(len(forward_y), len(reverse_y), len(forward_x), len(reverse_x))
    if n == 0:
        return 0.0
    x = 0.5 * (forward_x[:n] + reverse_x[:n])
    diff = np.abs(forward_y[:n] - reverse_y[:n])
    if n == 1:
        return float(diff[0])
    return float(np.trapezoid(diff, x))


def retrace_class(
    topology_return_error: float,
    flow_error: float,
    selectivity_error: float,
    coherence_error: float,
    area_flow: float,
    area_selectivity: float,
    area_coherence: float,
    common: dict[str, Any],
) -> tuple[str, int]:
    tol_flow = float(common['return_flow_tolerance'])
    tol_selectivity = float(common['return_selectivity_tolerance'])
    tol_coherence = float(common['return_coherence_tolerance'])
    tol_area = float(common['hysteresis_area_tolerance'])
    threshold_area = float(common['irreversibility_area_threshold'])
    max_area = max(area_flow, area_selectivity, area_coherence)
    if (
        topology_return_error < 0.5
        and flow_error <= tol_flow
        and selectivity_error <= tol_selectivity
        and coherence_error <= tol_coherence
        and max_area <= tol_area
    ):
        return 'reversible_retrace', 0
    if topology_return_error >= 0.5 or max_area >= threshold_area:
        return 'irreversible_witness', 1
    return 'hysteretic_retrace', 0


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
    dt_scale = float(common['dt_scale'])
    protocol_envelope = list(common['protocol_envelope'])

    complex_data = build_dk2d_complex(n_side=resolution, epsilon=epsilon)
    positions = pack_positions(complex_data.points, complex_data.edge_midpoints, complex_data.face_centers)
    block_sizes = complex_data.block_sizes
    packet_states = make_initial_packets(positions, block_sizes, common)
    base_operator = complex_data.dirac_kahler
    dt = first_order_dt(base_operator, dt_scale)
    sigma = float(common['mean_width'])
    rng = make_rng(str(run['run_id']))
    motif_nodes = identify_motif_nodes(base_operator, packet_states, common)
    skew_nodes = choose_skew_nodes(packet_states)

    psi = np.sum(packet_states, axis=0)
    schedule = [float(value) for value in run['delta_schedule']]
    segment_windows: list[dict[str, Any]] = []
    topology_trace: list[str] = []
    topology_codes: list[float] = []
    flow_trace: list[float] = []
    selectivity_trace: list[float] = []
    coherence_trace: list[float] = []
    segment_times: list[float] = []

    cumulative_steps = 0
    for idx, delta in enumerate(schedule):
        labels = build_family_labels_with_mode(
            packet_states=packet_states,
            block_sizes=block_sizes,
            delta=delta,
            family_mode=str(run['family_mode']),
            support_floor_fraction=support_floor_fraction,
            dominant_ratio_threshold=dominant_ratio_threshold,
            modulus=modulus,
            rng=rng,
        )
        operator = address_weighted_operator(base_operator, labels, modulus, eta_match, eta_mismatch)
        envelope = float(protocol_envelope[min(idx, len(protocol_envelope) - 1)])
        if str(run['protocol']) == 'degree_skew_cycle' and envelope > 0.0:
            operator = apply_node_bias(operator, skew_nodes, float(common['degree_bias_amplitude']) * envelope)
        elif str(run['protocol']) == 'motif_cycle' and envelope > 0.0:
            motif_bias = np.zeros(operator.shape[0], dtype=float)
            motif_bias[motif_nodes] = 1.0
            operator = apply_node_bias(operator, motif_bias, float(common['motif_bias_amplitude']) * envelope)

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
        segment_windows.append(window)
        label = str(window['topology_class'])
        topology_trace.append(label)
        topology_codes.append(topology_code(label))
        flow_trace.append(float(window['flow_concentration_index']))
        selectivity_trace.append(float(window['address_selectivity_index']))
        coherence_trace.append(float(window['exchange_coherence']))

    initial = segment_windows[0]
    final = segment_windows[-1]
    topology_error = 0.0 if str(initial['topology_class']) == str(final['topology_class']) else 1.0
    flow_error = abs(float(final['flow_concentration_index']) - float(initial['flow_concentration_index']))
    selectivity_error = abs(float(final['address_selectivity_index']) - float(initial['address_selectivity_index']))
    coherence_error = abs(float(final['exchange_coherence']) - float(initial['exchange_coherence']))
    area_flow = hysteresis_area(schedule, flow_trace)
    area_selectivity = hysteresis_area(schedule, selectivity_trace)
    area_coherence = hysteresis_area(schedule, coherence_trace)
    retrace_label, irreversible_flag = retrace_class(
        topology_return_error=topology_error,
        flow_error=flow_error,
        selectivity_error=selectivity_error,
        coherence_error=coherence_error,
        area_flow=area_flow,
        area_selectivity=area_selectivity,
        area_coherence=area_coherence,
        common=common,
    )

    row = {
        'run_id': str(run['run_id']),
        'geometry_id': 'tight_clustered_pair',
        'phase_id': str(run['protocol']),
        'resolution': resolution,
        'graph_type': common['graph_type'],
        'packet_count': 2,
        'initial_peak_count': int(initial['raw_peak_count']),
        'collision_label': str(final['topology_class']),
        'persistence_label': retrace_label,
        'composite_lifetime': 0.0,
        'binding_persistence': 0.0,
        'corridor_dwell': 0.0,
        'coarse_basin_persistence': 0.0,
        'minimum_separation': float(np.min(final['mean_pair_distances'])) if final['mean_pair_distances'].size else 0.0,
        'final_mean_separation': float(final['mean_pair_distances'][-1]) if final['mean_pair_distances'].size else 0.0,
        'post_collision_separation_trend': float(final['transport_span']),
        'encounter_dwell_time': 0.0,
        'deflection_angle_proxy': 0.0,
        'reflection_fraction': 0.0,
        'grade_transfer_amplitude': float(np.max(final['transfer_signal'])) if final['transfer_signal'].size else 0.0,
        'omega0_weight_initial': 0.0,
        'omega1_weight_initial': 0.0,
        'omega2_weight_initial': 0.0,
        'omega0_weight_final': 0.0,
        'omega1_weight_final': 0.0,
        'omega2_weight_final': 0.0,
        'localized_circulation_proxy': 0.0,
        'new_collision_class': int(retrace_label == 'irreversible_witness'),
        'gate_met': int(irreversible_flag),
        'promoted_followup': 0,
        'run_order': int(run['run_order']),
        'family_id': str(run['family_id']),
        'family_label': str(run['family_label']),
        'family_mode': str(run['family_mode']),
        'protocol': str(run['protocol']),
        'protocol_label': str(run['protocol_label']),
        'initial_delta': float(schedule[0]),
        'final_delta': float(schedule[-1]),
        'initial_topology_class': str(initial['topology_class']),
        'final_topology_class': str(final['topology_class']),
        'topology_return_error': topology_error,
        'flow_return_error': flow_error,
        'selectivity_return_error': selectivity_error,
        'coherence_return_error': coherence_error,
        'hysteresis_area_flow_concentration': area_flow,
        'hysteresis_area_selectivity': area_selectivity,
        'hysteresis_area_coherence': area_coherence,
        'irreversible_witness': irreversible_flag,
        'retrace_class': retrace_label,
        'peak_hysteresis_metric': max(area_flow, area_selectivity, area_coherence),
        'notes': str(run['notes']),
    }

    detail = {
        'run': run,
        'metrics': row,
        'delta_schedule': schedule,
        'segment_times': segment_times,
        'topology_trace': topology_trace,
        'flow_trace': flow_trace,
        'selectivity_trace': selectivity_trace,
        'coherence_trace': coherence_trace,
        'segment_windows': [
            {
                'delta': schedule[idx],
                'topology_class': str(window['topology_class']),
                'flow_concentration_index': float(window['flow_concentration_index']),
                'address_selectivity_index': float(window['address_selectivity_index']),
                'exchange_coherence': float(window['exchange_coherence']),
            }
            for idx, window in enumerate(segment_windows)
        ],
    }

    plot_paths: list[Path] = []
    topology_path = WORK_PLOT_DIR / f"stage_c0_15_{run['run_id']}_topology_trace.png"
    selectivity_path = WORK_PLOT_DIR / f"stage_c0_15_{run['run_id']}_selectivity_trace.png"
    plot_trace(
        topology_path,
        str(run['run_id']),
        segment_times,
        np.asarray(topology_codes, dtype=float),
        'topology code',
        'Stage C0.15 topology trace',
    )
    plot_trace(
        selectivity_path,
        str(run['run_id']),
        segment_times,
        np.asarray(selectivity_trace, dtype=float),
        'address selectivity',
        'Stage C0.15 selectivity trace',
    )
    plot_paths.extend([topology_path, selectivity_path])
    return row, detail, plot_paths


def plot_hysteresis_panel(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [str(row['run_id']) for row in rows]
    x = np.arange(len(rows))
    width = 0.24
    flow = np.asarray([float(row['hysteresis_area_flow_concentration']) for row in rows], dtype=float)
    selectivity = np.asarray([float(row['hysteresis_area_selectivity']) for row in rows], dtype=float)
    coherence = np.asarray([float(row['hysteresis_area_coherence']) for row in rows], dtype=float)
    fig, ax = plt.subplots(figsize=(14, 5.2))
    ax.bar(x - width, flow, width=width, label='flow')
    ax.bar(x, selectivity, width=width, label='selectivity')
    ax.bar(x + width, coherence, width=width, label='coherence')
    ax.set_xticks(x, labels, rotation=45, ha='right')
    ax.set_ylabel('hysteresis area')
    ax.set_title('Stage C0.15 hysteresis panel')
    ax.grid(axis='y', alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_return_error_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    labels = [str(row['run_id']) for row in rows]
    x = np.arange(len(rows))
    width = 0.24
    flow = np.asarray([float(row['flow_return_error']) for row in rows], dtype=float)
    selectivity = np.asarray([float(row['selectivity_return_error']) for row in rows], dtype=float)
    coherence = np.asarray([float(row['coherence_return_error']) for row in rows], dtype=float)
    topology = np.asarray([float(row['topology_return_error']) for row in rows], dtype=float)
    fig, ax = plt.subplots(figsize=(14, 5.2))
    ax.bar(x - width, flow, width=width, label='flow')
    ax.bar(x, selectivity, width=width, label='selectivity')
    ax.bar(x + width, coherence, width=width, label='coherence')
    ax.plot(x, topology, color='tab:red', marker='o', linestyle='none', label='topology error')
    ax.set_xticks(x, labels, rotation=45, ha='right')
    ax.set_ylabel('return error')
    ax.set_title('Stage C0.15 return-error summary')
    ax.grid(axis='y', alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def global_read(rows: list[dict[str, Any]]) -> str:
    irreversible = sum(int(row['irreversible_witness']) for row in rows)
    if irreversible >= 6:
        return 'broad return failure appears across the matrix, consistent with a genuine hysteresis-based irreversibility witness'
    if irreversible >= 3:
        return 'path dependence appears in part of the matrix, but broad irreversibility is not yet established'
    hysteretic = sum(1 for row in rows if str(row['retrace_class']) == 'hysteretic_retrace')
    if hysteretic:
        return 'the branch mostly retraces, but nonzero hysteresis remains in several runs without opening broad irreversibility'
    return 'the branch retraces cleanly enough that this scan does not improve the irreversibility case'


def write_note(json_path: Path, csv_path: Path, rows: list[dict[str, Any]], stamped_plots: list[str], verdict: str) -> None:
    lines = [
        '# Stage C0.15 Harmonic Hysteresis Return-Error Probe v1',
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
            f"- `{row['run_id']}`: family=`{row['family_label']}`, protocol=`{row['protocol_label']}`, retrace=`{row['retrace_class']}`, topology_return_error=`{row['topology_return_error']:.0f}`, peak_hysteresis=`{row['peak_hysteresis_metric']:.4f}`"
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

    hysteresis_path = WORK_PLOT_DIR / 'stage_c0_15_hysteresis_panel.png'
    return_path = WORK_PLOT_DIR / 'stage_c0_15_return_error_summary.png'
    plot_hysteresis_panel(hysteresis_path, rows)
    plot_return_error_summary(return_path, rows)
    plot_paths.extend([hysteresis_path, return_path])

    verdict = global_read(rows)
    summary = {
        'retrace_counts': dict(Counter(str(row['retrace_class']) for row in rows)),
        'irreversible_count': int(sum(int(row['irreversible_witness']) for row in rows)),
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
        'stage_c0_15_harmonic_hysteresis_return_error_probe',
        result,
        rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots, verdict)


if __name__ == '__main__':
    main()
