#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import scipy.sparse as sp

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt, save_atlas_payload
from stage_c0_17_bridge_cycle_latch_withdrawal_probe import (
    BRIDGE_RUNSHEET_PATH,
    apply_basis_latch,
    bridge_setup,
    dwell_time_from_trace,
    estimate_lambda_max,
    one_step_cn,
    recurrence_indicator,
    retention_label,
    window_metrics,
)
from stage_c0_dk_bridge_validation import (
    address_weighted_operator,
    block_dirac,
    build_address_labels,
    build_d0,
    build_d1,
    mismatch_scores,
    orient_edges,
    representative_graph,
    simulate_family_variant,
    triangle_list,
)


RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_19_peak_braid_withdrawal_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_C0_19_Peak_Braid_Withdrawal_Probe_v1.md'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage_c0_19_peak_braid'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'run_id',
    'family_id',
    'family_label',
    'protocol_id',
    'protocol_label',
    'armed_topology_class',
    'armed_window_start_step',
    'armed_window_end_step',
    'armed_window_length',
    'armed_peak_loop_score',
    'post_withdrawal_topology_class',
    'withdrawal_braid_dwell_time',
    'loop_retention_score',
    'path_reuse_retention_score',
    'selector_off_recurrence_indicator',
    'locking_label',
    'trace_plot',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.19 peak-braid withdrawal probe.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def selected_pairs(
    families: list[dict[str, Any]],
    protocols: list[dict[str, Any]],
    run_ids: list[str],
) -> list[tuple[str, dict[str, Any], dict[str, Any]]]:
    pairs: list[tuple[str, dict[str, Any], dict[str, Any]]] = []
    for family in families:
        for protocol in protocols:
            run_id = f"C019_{family['family_id']}_{protocol['protocol_id']}"
            pairs.append((run_id, family, protocol))
    if not run_ids:
        return pairs
    wanted = set(run_ids)
    return [item for item in pairs if item[0] in wanted]


def baseline_result_for_family(family: dict[str, Any], bridge_common: dict[str, Any]) -> dict[str, float]:
    result = simulate_family_variant(
        family=family,
        common=bridge_common,
        delta=0.0,
        left_sigma_scale=1.0,
        right_sigma_scale=1.0,
        grade0_scale=float(bridge_common['grade0_scale']),
        grade1_scale=float(bridge_common['grade1_scale']),
        grade2_scale=float(bridge_common['grade2_scale']),
        degree_skew_strength=0.0,
    )
    return {
        'flow_concentration_index': float(result['flow_concentration_index']),
        'grade_exchange_coherence': float(result['grade_exchange_coherence']),
        'address_selectivity_index': float(result['address_selectivity_index']),
        'loop_score': float(result['loop_score']),
        'grade_asymmetry_index': float(result['grade_asymmetry_index']),
    }


def base_operator_for_family(family: dict[str, Any]) -> sp.csr_matrix:
    family_id = str(family['family_id'])
    G = representative_graph(family_id)
    left_seed = int(family['left_seed'])
    right_seed = int(family['right_seed'])
    nodes = sorted(int(node) for node in G.nodes())
    scores = mismatch_scores(G, left_seed=left_seed, right_seed=right_seed)
    oriented_edges = orient_edges(G, scores)
    triangles = triangle_list(G)
    d0, edge_lookup, _node_to_idx = build_d0(nodes, oriented_edges)
    d1 = build_d1(triangles, edge_lookup)
    D, _Delta = block_dirac(d0, d1)
    return D


def operator_for_delta(
    bridge: dict[str, Any],
    family: dict[str, Any],
    bridge_common: dict[str, Any],
    delta: float,
    frozen_weights: np.ndarray | None = None,
    latch_strength: float = 0.0,
) -> tuple[sp.csr_matrix, np.ndarray]:
    base_operator = base_operator_for_family(family)
    labels = build_address_labels(
        packet_states=bridge['packet_states'],
        block_sizes=bridge['block_sizes'],
        left_address=0.0,
        right_address=2.0 * float(delta),
        support_floor_fraction=float(bridge_common['support_floor_fraction']),
        dominant_ratio_threshold=float(bridge_common['dominant_ratio_threshold']),
    )
    operator = address_weighted_operator(
        base_operator=base_operator,
        labels=labels,
        modulus=float(bridge_common['address_modulus']),
        eta_match=float(bridge_common['eta_match']),
        eta_mismatch=float(bridge_common['eta_mismatch']),
    )
    if frozen_weights is not None and latch_strength > 0.0:
        operator = apply_basis_latch(operator, frozen_weights, latch_strength)
    return operator, labels


def arm_peak_window(
    bridge: dict[str, Any],
    family: dict[str, Any],
    stage_common: dict[str, Any],
    bridge_common: dict[str, Any],
    baseline_result: dict[str, float],
) -> dict[str, Any]:
    exact_operator, _labels = operator_for_delta(bridge, family, bridge_common, delta=0.0)
    dt = float(stage_common['dt_scale']) / max(estimate_lambda_max(exact_operator), 1.0e-12)
    state = bridge['psi0'].copy()
    states = [state.copy()]
    records: list[dict[str, Any]] = []

    for step in range(1, int(stage_common['arming_steps']) + 1):
        state = one_step_cn(exact_operator, state, dt)
        states.append(state.copy())
        metrics = window_metrics(
            states=states[max(0, len(states) - 16):],
            bridge=bridge,
            labels=bridge['exact_labels'],
            operator=exact_operator,
            baseline_result=baseline_result,
            common=bridge_common,
        )
        records.append({'step': step, 'state': state.copy(), 'metrics': metrics})

    windows: list[list[dict[str, Any]]] = []
    current: list[dict[str, Any]] = []
    for record in records:
        if str(record['metrics']['topology_class']) == 'braid_like_exchange':
            current.append(record)
        elif current:
            windows.append(current)
            current = []
    if current:
        windows.append(current)

    if windows:
        chosen = max(
            windows,
            key=lambda window: (
                len(window),
                float(window[-1]['metrics']['loop_score']),
                float(window[-1]['metrics']['flow_concentration_index']),
            ),
        )
        armed_topology = 'braid_like_exchange'
    else:
        chosen = [max(records, key=lambda item: (float(item['metrics']['loop_score']), float(item['metrics']['flow_concentration_index'])))]
        armed_topology = str(chosen[-1]['metrics']['topology_class'])

    n0, n1, n2 = bridge['block_sizes']
    accum_edge = np.zeros(len(bridge['oriented_edges']), dtype=float)
    accum_face = np.zeros(len(bridge['triangles']), dtype=float)
    for record in chosen:
        state_i = record['state']
        accum_edge += np.abs(state_i[n0:n0 + n1]) ** 2
        if n2 > 0:
            accum_face += np.abs(state_i[n0 + n1:n0 + n1 + n2]) ** 2
    frozen_weights = np.zeros(sum(bridge['block_sizes']), dtype=float)
    frozen_weights[n0:n0 + n1] = accum_edge
    if n2 > 0:
        frozen_weights[n0 + n1:n0 + n1 + n2] = accum_face
    if np.max(frozen_weights) > 0.0:
        frozen_weights = frozen_weights / float(np.max(frozen_weights))

    return {
        'dt': dt,
        'armed_state': chosen[-1]['state'].copy(),
        'armed_metrics': chosen[-1]['metrics'],
        'armed_topology_class': armed_topology,
        'window_start_step': int(chosen[0]['step']),
        'window_end_step': int(chosen[-1]['step']),
        'window_length': int(len(chosen)),
        'frozen_weights': frozen_weights,
    }


def run_protocol(
    protocol_id: str,
    family: dict[str, Any],
    bridge: dict[str, Any],
    bridge_common: dict[str, Any],
    stage_common: dict[str, Any],
    baseline_result: dict[str, float],
    armed: dict[str, Any],
) -> tuple[dict[str, Any], list[str], list[float]]:
    latch_strength = float(stage_common['combined_latch_strength'])
    protocol_trace: list[str] = []
    latch_trace: list[float] = [0.0]
    state = armed['armed_state'].copy()
    states = [state.copy()]

    if protocol_id == 'peak_abrupt_release':
        schedule = [float(stage_common['withdrawal_delta'])] * int(stage_common['withdrawal_steps'])
        use_latch = False
    elif protocol_id == 'peak_frozen_combined_release':
        schedule = [float(stage_common['withdrawal_delta'])] * int(stage_common['withdrawal_steps'])
        use_latch = True
    else:
        ramp_steps = int(stage_common['ramp_steps'])
        tail = max(int(stage_common['withdrawal_steps']) - ramp_steps, 0)
        schedule = [float(stage_common['withdrawal_delta']) * (step + 1) / ramp_steps for step in range(ramp_steps)]
        schedule.extend([float(stage_common['withdrawal_delta'])] * tail)
        use_latch = True

    for delta in schedule:
        operator, labels = operator_for_delta(
            bridge=bridge,
            family=family,
            bridge_common=bridge_common,
            delta=delta,
            frozen_weights=armed['frozen_weights'] if use_latch else None,
            latch_strength=latch_strength if use_latch else 0.0,
        )
        state = one_step_cn(operator, state, armed['dt'])
        states.append(state.copy())
        metrics = window_metrics(
            states=states[max(0, len(states) - 12):],
            bridge=bridge,
            labels=labels,
            operator=operator,
            baseline_result=baseline_result,
            common=bridge_common,
        )
        protocol_trace.append(str(metrics['topology_class']))
        duty = 0.0
        if use_latch and np.sum(armed['frozen_weights']) > 0.0:
            duty = float(np.max(armed['frozen_weights']) / max(float(np.sum(armed['frozen_weights'])), 1.0e-12))
        latch_trace.append(duty)

    post_operator, post_labels = operator_for_delta(
        bridge=bridge,
        family=family,
        bridge_common=bridge_common,
        delta=float(stage_common['withdrawal_delta']),
        frozen_weights=armed['frozen_weights'] if use_latch else None,
        latch_strength=latch_strength if use_latch else 0.0,
    )
    post_metrics = window_metrics(
        states=states[max(0, len(states) - 16):],
        bridge=bridge,
        labels=post_labels,
        operator=post_operator,
        baseline_result=baseline_result,
        common=bridge_common,
    )
    row = {
        'armed_topology_class': str(armed['armed_topology_class']),
        'armed_window_start_step': int(armed['window_start_step']),
        'armed_window_end_step': int(armed['window_end_step']),
        'armed_window_length': int(armed['window_length']),
        'armed_peak_loop_score': float(armed['armed_metrics']['loop_score']),
        'post_withdrawal_topology_class': str(post_metrics['topology_class']),
        'withdrawal_braid_dwell_time': dwell_time_from_trace(protocol_trace, 'braid_like_exchange'),
        'loop_retention_score': float(post_metrics['loop_score']) / max(float(armed['armed_metrics']['loop_score']), 1.0e-12),
        'path_reuse_retention_score': float(post_metrics['flow_concentration_index']) / max(float(armed['armed_metrics']['flow_concentration_index']), 1.0e-12),
        'selector_off_recurrence_indicator': recurrence_indicator(
            np.asarray([1.0 if label == 'braid_like_exchange' else 0.0 for label in protocol_trace], dtype=float),
            close_threshold=0.5,
        ),
    }
    if row['armed_topology_class'] != 'braid_like_exchange':
        row['locking_label'] = 'no_armed_braid_window'
    else:
        row['locking_label'] = retention_label(
            armed_topology=str(row['armed_topology_class']),
            post_topology=str(row['post_withdrawal_topology_class']),
            dwell=float(row['withdrawal_braid_dwell_time']),
            loop_retention_score=float(row['loop_retention_score']),
            path_retention_score=float(row['path_reuse_retention_score']),
        )
    return row, protocol_trace, latch_trace


def trace_plot(path: Path, run_id: str, topology_trace: list[str], latch_trace: list[float]) -> None:
    mapping = {'braid_like_exchange': 0, 'transfer_smeared': 1, 'unresolved_mixed': 2}
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    axes[0].plot(range(len(topology_trace)), [mapping.get(label, 2) for label in topology_trace], marker='o')
    axes[0].set_yticks([0, 1, 2])
    axes[0].set_yticklabels(['braid', 'smeared', 'mixed'])
    axes[0].set_title('Withdrawal topology trace')
    axes[0].grid(alpha=0.25)
    axes[1].plot(range(len(latch_trace)), latch_trace, marker='o', color='tab:red')
    axes[1].set_title('Frozen-latch duty trace')
    axes[1].grid(alpha=0.25)
    fig.suptitle(run_id)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def summary_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    x = np.arange(len(rows))
    dwell = [float(row['withdrawal_braid_dwell_time']) for row in rows]
    window = [float(row['armed_window_length']) for row in rows]
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.5))
    axes[0].bar(x, window, color='tab:blue')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([str(row['run_id']) for row in rows], rotation=30, ha='right')
    axes[0].set_title('Recovered braid window length')
    axes[0].grid(alpha=0.25, axis='y')
    axes[1].bar(x, dwell, color='tab:orange')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([str(row['run_id']) for row in rows], rotation=30, ha='right')
    axes[1].set_title('Post-withdrawal braid dwell')
    axes[1].grid(alpha=0.25, axis='y')
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def note_text(json_rel: str, csv_rel: str, rows: list[dict[str, Any]], plots: list[str]) -> str:
    counts = Counter(str(row['locking_label']) for row in rows)
    global_read = 'fair bridge withdrawal still does not retain topology after selector removal'
    if counts.get('locked_after_withdrawal', 0) > 0 or counts.get('quasi_locked_retention', 0) > 0:
        global_read = 'fair bridge withdrawal opens partial topology retention after selector removal'
    lines = [
        '# Stage C0.19 Peak-Braid Withdrawal Probe v1',
        '',
        f'Timestamped JSON: `{json_rel}`',
        f'Timestamped CSV: `{csv_rel}`',
        '',
        f'Global read: `{global_read}`',
        '',
        '## Per-run summary',
    ]
    for row in rows:
        lines.append(
            f"- `{row['run_id']}`: family=`{row['family_label']}`, protocol=`{row['protocol_label']}`, "
            f"armed=`{row['armed_topology_class']}`, armed_window=`{row['armed_window_start_step']}-{row['armed_window_end_step']}` "
            f"(len `{int(row['armed_window_length'])}`), post=`{row['post_withdrawal_topology_class']}`, "
            f"dwell=`{float(row['withdrawal_braid_dwell_time']):.4f}`, label=`{row['locking_label']}`"
        )
    lines.extend(['', '## Plots'])
    for rel in plots:
        lines.append(f'- `{rel}`')
    return '\n'.join(lines) + '\n'


def main() -> None:
    args = parse_args()
    runsheet = load_runsheet(args.runsheet)
    bridge_common = json.loads(BRIDGE_RUNSHEET_PATH.read_text(encoding='utf-8'))['common_fields']
    selected = selected_pairs(runsheet['families'], runsheet['protocols'], args.run_ids)

    rows: list[dict[str, Any]] = []
    plot_paths: list[Path] = []

    for run_id, family, protocol in selected:
        bridge = bridge_setup(family, bridge_common, runsheet['common_fields'])
        baseline_result = baseline_result_for_family(family, bridge_common)
        armed = arm_peak_window(bridge, family, runsheet['common_fields'], bridge_common, baseline_result)
        result, topology_trace, latch_trace = run_protocol(
            protocol_id=str(protocol['protocol_id']),
            family=family,
            bridge=bridge,
            bridge_common=bridge_common,
            stage_common=runsheet['common_fields'],
            baseline_result=baseline_result,
            armed=armed,
        )
        result.update(
            {
                'run_id': run_id,
                'family_id': str(family['family_id']),
                'family_label': str(family['family_label']),
                'protocol_id': str(protocol['protocol_id']),
                'protocol_label': str(protocol['protocol_label']),
                'notes': '',
            }
        )
        trace_path = WORK_PLOT_DIR / f'{run_id}_trace.png'
        trace_plot(trace_path, run_id, topology_trace, latch_trace)
        result['trace_plot'] = str(trace_path)
        rows.append(result)
        plot_paths.append(trace_path)

    summary_path = WORK_PLOT_DIR / 'stage_c0_19_withdrawal_summary.png'
    summary_plot(summary_path, rows)
    plot_paths.append(summary_path)

    result_payload = {
        'stage': 'c0_19_peak_braid_withdrawal_probe',
        'runsheet': str(Path(args.runsheet).relative_to(REPO_ROOT)),
        'rows': rows,
    }
    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug='stage_c0_19_peak_braid_withdrawal_probe',
        result=result_payload,
        csv_rows=rows,
        csv_fieldnames=CSV_FIELDS,
        plot_paths=plot_paths,
    )
    note = note_text(
        str(json_path.relative_to(REPO_ROOT)),
        str(csv_path.relative_to(REPO_ROOT)),
        rows,
        stamped_plots,
    )
    NOTE_PATH.write_text(note, encoding='utf-8')


if __name__ == '__main__':
    main()
