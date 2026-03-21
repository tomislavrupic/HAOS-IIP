#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import ATLAS_NOTES, REPO_ROOT, save_atlas_payload
from stage_c0_17_bridge_cycle_latch_withdrawal_probe import (
    BRIDGE_RUNSHEET_PATH,
    CSV_FIELDS,
    WORK_PLOT_DIR,
    address_activity_statistics,
    apply_basis_latch,
    bridge_setup,
    build_latch_weights,
    classify_relative_topology,
    dwell_time_from_trace,
    edge_component_count,
    estimate_lambda_max,
    grade_exchange_signal,
    grade_histories,
    latch_strength_for_rule,
    note_text as _unused_note_text,
    one_step_cn,
    recurrence_indicator,
    retention_label,
    summary_plot,
    trace_plot,
    triangle_loop_metrics,
    window_metrics,
)


RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_18_deferred_bridge_latch_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_C0_18_Deferred_Bridge_Latch_Withdrawal_Probe_v1.md'


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.18 deferred bridge-latch withdrawal probe.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def selected_pairs(
    families: list[dict[str, Any]],
    latches: list[dict[str, Any]],
    run_ids: list[str],
) -> list[tuple[str, dict[str, Any], dict[str, Any]]]:
    pairs: list[tuple[str, dict[str, Any], dict[str, Any]]] = []
    for family in families:
        for latch in latches:
            run_id = f"C018_{family['family_id']}_{latch['latch_id']}"
            pairs.append((run_id, family, latch))
    if not run_ids:
        return pairs
    wanted = set(run_ids)
    return [item for item in pairs if item[0] in wanted]


def note_text(json_rel: str, csv_rel: str, rows: list[dict[str, Any]], plots: list[str]) -> str:
    counts = Counter(str(row['locking_label']) for row in rows)
    global_read = 'deferred bridge-local latching still does not retain topology after selector withdrawal'
    if counts.get('locked_after_withdrawal', 0) > 0 or counts.get('quasi_locked_retention', 0) > 0:
        global_read = 'deferred bridge-local latching opens partial retention after selector withdrawal'
    lines = [
        '# Stage C0.18 Deferred Bridge Latch Withdrawal Probe v1',
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
            f"- `{row['run_id']}`: family=`{row['family_label']}`, latch=`{row['latch_label']}`, "
            f"armed=`{row['armed_topology_class']}`, post=`{row['post_withdrawal_topology_class']}`, "
            f"dwell=`{float(row['withdrawal_braid_dwell_time']):.4f}`, label=`{row['locking_label']}`"
        )
    lines.extend(['', '## Plots'])
    for rel in plots:
        lines.append(f'- `{rel}`')
    return '\n'.join(lines) + '\n'


def simulate_run(run_id: str, family: dict[str, Any], latch: dict[str, Any], common: dict[str, Any], bridge_common: dict[str, Any]) -> tuple[dict[str, Any], list[Path]]:
    bridge = bridge_setup(family, bridge_common, common)
    baseline_result = {
        'flow_concentration_index': 0.2907950679287607 if bridge['family_id'] == 'clustered_composite_anchor' else (0.4823 if bridge['family_id'] == 'counter_propagating_corridor' else 0.2820),
        'grade_exchange_coherence': 0.49694013187909375 if bridge['family_id'] == 'clustered_composite_anchor' else (0.5734 if bridge['family_id'] == 'counter_propagating_corridor' else 0.5822),
        'address_selectivity_index': 1.0,
        'loop_score': 0.21852819979502663 if bridge['family_id'] == 'clustered_composite_anchor' else (0.1142 if bridge['family_id'] == 'counter_propagating_corridor' else 0.1034),
        'grade_asymmetry_index': 0.18,
    }

    latch_id = str(latch['latch_id'])
    latch_strength = latch_strength_for_rule(latch_id, common)
    upper_arm = bridge['exact_operator']
    upper_withdraw = apply_basis_latch(bridge['withdrawal_operator'], np.ones_like(bridge['psi0'], dtype=float), latch_strength)
    lam_upper = max(estimate_lambda_max(upper_arm), estimate_lambda_max(upper_withdraw), 1.0e-12)
    dt = float(common['dt_scale']) / lam_upper

    accum_edge = np.zeros(len(bridge['oriented_edges']), dtype=float)
    accum_face = np.zeros(len(bridge['triangles']), dtype=float)
    state = bridge['psi0'].copy()
    arming_states = [state.copy()]
    latch_trace = [0.0]
    topology_trace: list[str] = []
    frozen_weights = np.zeros_like(bridge['psi0'], dtype=float)

    for _ in range(int(common['arming_steps'])):
        frozen_weights = build_latch_weights(latch_id, state, accum_edge, accum_face, bridge['block_sizes'])
        state = one_step_cn(bridge['exact_operator'], state, dt)
        arming_states.append(state.copy())
        duty = float(np.max(frozen_weights) / max(float(np.sum(frozen_weights)), 1.0e-12)) if np.sum(frozen_weights) > 0 else 0.0
        latch_trace.append(duty)

    armed_metrics = window_metrics(
        states=arming_states[max(0, len(arming_states) - 16):],
        bridge=bridge,
        labels=bridge['exact_labels'],
        operator=bridge['exact_operator'],
        baseline_result=baseline_result,
        common=common,
    )
    armed_topology = str(armed_metrics['topology_class'])
    topology_trace.extend([armed_topology] * len(arming_states))

    frozen_weights = frozen_weights / max(float(np.max(frozen_weights)), 1.0e-12) if np.max(frozen_weights) > 0 else frozen_weights
    withdrawal_operator = apply_basis_latch(bridge['withdrawal_operator'], frozen_weights, latch_strength)
    withdrawal_states = [state.copy()]
    withdrawal_labels_trace: list[str] = []
    for _ in range(int(common['withdrawal_steps'])):
        state = one_step_cn(withdrawal_operator, state, dt)
        withdrawal_states.append(state.copy())
        metrics = window_metrics(
            states=withdrawal_states[max(0, len(withdrawal_states) - 12):],
            bridge=bridge,
            labels=bridge['withdrawal_labels'],
            operator=withdrawal_operator,
            baseline_result=baseline_result,
            common=common,
        )
        withdrawal_labels_trace.append(str(metrics['topology_class']))
        topology_trace.append(str(metrics['topology_class']))
        duty = float(np.max(frozen_weights) / max(float(np.sum(frozen_weights)), 1.0e-12)) if np.sum(frozen_weights) > 0 else 0.0
        latch_trace.append(duty)

    post_metrics = window_metrics(
        states=withdrawal_states[max(0, len(withdrawal_states) - 16):],
        bridge=bridge,
        labels=bridge['withdrawal_labels'],
        operator=withdrawal_operator,
        baseline_result=baseline_result,
        common=common,
    )
    post_topology = str(post_metrics['topology_class'])
    withdrawal_braid_dwell = dwell_time_from_trace(withdrawal_labels_trace, 'braid_like_exchange')
    armed_loop = float(armed_metrics['loop_score'])
    post_loop = float(post_metrics['loop_score'])
    loop_retention_score = float(post_loop / max(armed_loop, 1.0e-12)) if armed_loop > 1.0e-12 else 0.0
    armed_path = float(armed_metrics['flow_concentration_index'])
    post_path = float(post_metrics['flow_concentration_index'])
    path_reuse_retention_score = float(post_path / max(armed_path, 1.0e-12)) if armed_path > 1.0e-12 else 0.0
    latch_duty_cycle = float(np.mean(latch_trace))
    series = np.asarray([float(metrics == 'braid_like_exchange') for metrics in withdrawal_labels_trace], dtype=float)
    selector_off_recurrence_indicator = recurrence_indicator(series, close_threshold=0.5)
    label = retention_label(
        armed_topology=armed_topology,
        post_topology=post_topology,
        dwell=withdrawal_braid_dwell,
        loop_retention_score=loop_retention_score,
        path_retention_score=path_reuse_retention_score,
    )

    trace_path = WORK_PLOT_DIR / f'{run_id}_trace.png'
    trace_plot(trace_path, run_id, topology_trace, latch_trace)

    row = {
        'run_id': run_id,
        'family_id': str(family['family_id']),
        'family_label': str(family['family_label']),
        'latch_id': latch_id,
        'latch_label': str(latch['latch_label']),
        'armed_topology_class': armed_topology,
        'post_withdrawal_topology_class': post_topology,
        'withdrawal_braid_dwell_time': withdrawal_braid_dwell,
        'loop_retention_score': loop_retention_score,
        'path_reuse_retention_score': path_reuse_retention_score,
        'latch_duty_cycle': latch_duty_cycle,
        'flow_concentration_index_post': float(post_metrics['flow_concentration_index']),
        'grade_exchange_coherence_post': float(post_metrics['grade_exchange_coherence']),
        'selector_off_recurrence_indicator': selector_off_recurrence_indicator,
        'locking_label': label,
        'notes': '',
        'trace_plot': trace_path,
    }
    return row, [trace_path]


def main() -> None:
    args = parse_args()
    runsheet = load_runsheet(args.runsheet)
    bridge_common = json.loads(BRIDGE_RUNSHEET_PATH.read_text(encoding='utf-8'))['common_fields']
    selected = selected_pairs(runsheet['families'], runsheet['latch_rules'], args.run_ids)
    rows: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    for run_id, family, latch in selected:
        row, run_plots = simulate_run(run_id, family, latch, runsheet['common_fields'], bridge_common)
        rows.append(row)
        plot_paths.extend(run_plots)

    summary_path = WORK_PLOT_DIR / 'stage_c0_18_withdrawal_summary.png'
    summary_plot(summary_path, rows)
    plot_paths.append(summary_path)

    payload = {
        'stage': runsheet['stage'],
        'description': runsheet['description'],
        'rows': [{key: value for key, value in row.items() if key != 'trace_plot'} for row in rows],
    }
    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug='stage_c0_18_deferred_bridge_latch_withdrawal_probe',
        result=payload,
        csv_rows=[{key: value for key, value in row.items() if key in CSV_FIELDS} for row in rows],
        csv_fieldnames=CSV_FIELDS,
        plot_paths=plot_paths,
    )
    NOTE_PATH.write_text(
        note_text(
            json_rel=str(json_path.relative_to(REPO_ROOT)),
            csv_rel=str(csv_path.relative_to(REPO_ROOT)),
            rows=rows,
            plots=stamped_plots,
        ),
        encoding='utf-8',
    )
    print(f'JSON: {json_path}')
    print(f'CSV: {csv_path}')
    print(f'NOTE: {NOTE_PATH}')
    for rel in stamped_plots:
        print(f'PLOT: {REPO_ROOT / rel}')


if __name__ == '__main__':
    main()
