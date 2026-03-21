#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt, save_atlas_payload
from stage_c0_17_bridge_cycle_latch_withdrawal_probe import (
    apply_basis_latch,
    bridge_setup,
    one_step_cn,
    window_metrics,
)
from stage_c0_19_peak_braid_withdrawal_probe import (
    baseline_result_for_family,
    operator_for_delta,
    arm_peak_window,
)


RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_20_armed_bridge_hysteresis_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_C0_20_Armed_Bridge_Hysteresis_Probe_v1.md'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage_c0_20_bridge_hysteresis'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'run_id',
    'family_id',
    'family_label',
    'protocol_id',
    'protocol_label',
    'armed_topology_class',
    'armed_window_length',
    'final_return_topology_class',
    'topology_return_error',
    'return_state_overlap',
    'hysteresis_area_flow',
    'hysteresis_area_loop',
    'braid_fraction_during_cycle',
    'irreversibility_label',
    'trace_plot',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.20 armed bridge hysteresis probe.')
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
            run_id = f"C020_{family['family_id']}_{protocol['protocol_id']}"
            pairs.append((run_id, family, protocol))
    if not run_ids:
        return pairs
    wanted = set(run_ids)
    return [item for item in pairs if item[0] in wanted]


def irreversibility_label(
    armed_topology: str,
    final_topology: str,
    return_error: int,
    area_flow: float,
    area_loop: float,
) -> str:
    if armed_topology != 'braid_like_exchange':
        return 'no_armed_braid_window'
    if final_topology != 'braid_like_exchange' or return_error > 0:
        return 'irreversible_witness'
    if area_flow >= 0.03 or area_loop >= 0.01:
        return 'hysteretic_retrace'
    return 'reversible_retrace'


def area_between_curves(forward: list[float], reverse: list[float]) -> float:
    if not forward or not reverse:
        return 0.0
    n = min(len(forward), len(reverse))
    if n == 0:
        return 0.0
    f = np.asarray(forward[:n], dtype=float)
    r = np.asarray(list(reversed(reverse[-n:])), dtype=float)
    return float(np.mean(np.abs(f - r)))


def trace_plot(path: Path, run_id: str, deltas: list[float], topology_codes: list[int], flow: list[float], loop: list[float]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.1))
    axes[0].plot(range(len(topology_codes)), topology_codes, marker='o')
    axes[0].set_yticks([0, 1, 2])
    axes[0].set_yticklabels(['braid', 'smeared', 'mixed'])
    axes[0].set_title('Topology during cycle')
    axes[0].grid(alpha=0.25)
    axes[1].plot(range(len(deltas)), deltas, marker='o', color='tab:orange')
    axes[1].set_title('Detuning schedule')
    axes[1].grid(alpha=0.25)
    axes[2].plot(range(len(flow)), flow, marker='o', label='flow')
    axes[2].plot(range(len(loop)), loop, marker='s', label='loop')
    axes[2].set_title('Flow / loop response')
    axes[2].legend(fontsize=8)
    axes[2].grid(alpha=0.25)
    fig.suptitle(run_id)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def summary_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    x = np.arange(len(rows))
    flow_area = [float(row['hysteresis_area_flow']) for row in rows]
    loop_area = [float(row['hysteresis_area_loop']) for row in rows]
    fig, ax = plt.subplots(figsize=(12.0, 4.5))
    ax.plot(x, flow_area, marker='o', label='flow hysteresis area')
    ax.plot(x, loop_area, marker='s', label='loop hysteresis area')
    ax.set_xticks(x)
    ax.set_xticklabels([str(row['run_id']) for row in rows], rotation=30, ha='right')
    ax.set_title('Bridge hysteresis areas')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def note_text(json_rel: str, csv_rel: str, rows: list[dict[str, Any]], plots: list[str]) -> str:
    counts = Counter(str(row['irreversibility_label']) for row in rows)
    global_read = 'armed bridge hysteresis remains mostly retraceable'
    if counts.get('irreversible_witness', 0) > 0:
        global_read = 'armed bridge hysteresis opens explicit irreversibility witnesses even without locking'
    lines = [
        '# Stage C0.20 Armed Bridge Hysteresis Probe v1',
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
            f"armed=`{row['armed_topology_class']}`, final=`{row['final_return_topology_class']}`, "
            f"return_error=`{row['topology_return_error']}`, flow_area=`{float(row['hysteresis_area_flow']):.4f}`, "
            f"loop_area=`{float(row['hysteresis_area_loop']):.4f}`, label=`{row['irreversibility_label']}`"
        )
    lines.extend(['', '## Plots'])
    for rel in plots:
        lines.append(f'- `{rel}`')
    return '\n'.join(lines) + '\n'


def main() -> None:
    args = parse_args()
    runsheet = load_runsheet(args.runsheet)
    bridge_common = json.loads((REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_dk_bridge_validation_runs.json').read_text(encoding='utf-8'))['common_fields']
    selected = selected_pairs(runsheet['families'], runsheet['protocols'], args.run_ids)

    rows: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    mapping = {'braid_like_exchange': 0, 'transfer_smeared': 1, 'unresolved_mixed': 2}

    for run_id, family, protocol in selected:
        bridge = bridge_setup(family, bridge_common, {'withdrawal_delta': 1.0})
        baseline_result = baseline_result_for_family(family, bridge_common)
        armed = arm_peak_window(bridge, family, runsheet['common_fields'], bridge_common, baseline_result)

        state = armed['armed_state'].copy()
        states = [state.copy()]
        flow_trace: list[float] = []
        loop_trace: list[float] = []
        topology_trace: list[str] = []
        delta_trace: list[float] = []

        if str(protocol['protocol_id']) == 'plain_forward_reverse':
            schedule = [step / int(runsheet['common_fields']['forward_steps']) for step in range(int(runsheet['common_fields']['forward_steps']) + 1)]
            schedule.extend([step / int(runsheet['common_fields']['forward_steps']) for step in range(int(runsheet['common_fields']['forward_steps']) - 1, -1, -1)])
            use_latch = False
        elif str(protocol['protocol_id']) == 'frozen_combined_forward_reverse':
            schedule = [step / int(runsheet['common_fields']['forward_steps']) for step in range(int(runsheet['common_fields']['forward_steps']) + 1)]
            schedule.extend([step / int(runsheet['common_fields']['forward_steps']) for step in range(int(runsheet['common_fields']['forward_steps']) - 1, -1, -1)])
            use_latch = True
        else:
            schedule = [step / int(runsheet['common_fields']['forward_steps']) for step in range(int(runsheet['common_fields']['forward_steps']) + 1)]
            schedule.extend([1.0] * int(runsheet['common_fields']['hold_steps']))
            schedule.extend([step / int(runsheet['common_fields']['forward_steps']) for step in range(int(runsheet['common_fields']['forward_steps']) - 1, -1, -1)])
            use_latch = True

        for delta in schedule:
            operator, labels = operator_for_delta(
                bridge=bridge,
                family=family,
                bridge_common=bridge_common,
                delta=float(delta),
                frozen_weights=armed['frozen_weights'] if use_latch else None,
                latch_strength=float(runsheet['common_fields']['combined_latch_strength']) if use_latch else 0.0,
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
            delta_trace.append(float(delta))
            topology_trace.append(str(metrics['topology_class']))
            flow_trace.append(float(metrics['flow_concentration_index']))
            loop_trace.append(float(metrics['loop_score']))

        operator0, labels0 = operator_for_delta(
            bridge=bridge,
            family=family,
            bridge_common=bridge_common,
            delta=0.0,
            frozen_weights=armed['frozen_weights'] if use_latch else None,
            latch_strength=float(runsheet['common_fields']['combined_latch_strength']) if use_latch else 0.0,
        )
        final_metrics = None
        for _ in range(int(runsheet['common_fields']['return_settle_steps'])):
            state = one_step_cn(operator0, state, armed['dt'])
            states.append(state.copy())
            final_metrics = window_metrics(
                states=states[max(0, len(states) - 12):],
                bridge=bridge,
                labels=labels0,
                operator=operator0,
                baseline_result=baseline_result,
                common=bridge_common,
            )
            topology_trace.append(str(final_metrics['topology_class']))
            flow_trace.append(float(final_metrics['flow_concentration_index']))
            loop_trace.append(float(final_metrics['loop_score']))
            delta_trace.append(0.0)

        assert final_metrics is not None
        forward_count = int(runsheet['common_fields']['forward_steps']) + 1
        reverse_count = int(runsheet['common_fields']['forward_steps']) + 1
        if str(protocol['protocol_id']) == 'frozen_hold_forward_reverse':
            reverse_start = forward_count + int(runsheet['common_fields']['hold_steps'])
        else:
            reverse_start = forward_count
        forward_flow = flow_trace[:forward_count]
        reverse_flow = flow_trace[reverse_start:reverse_start + reverse_count]
        forward_loop = loop_trace[:forward_count]
        reverse_loop = loop_trace[reverse_start:reverse_start + reverse_count]
        return_overlap = float(np.abs(np.vdot(armed['armed_state'], state)))
        final_label = str(final_metrics['topology_class'])
        return_error = int(final_label != str(armed['armed_topology_class']))
        row = {
            'run_id': run_id,
            'family_id': str(family['family_id']),
            'family_label': str(family['family_label']),
            'protocol_id': str(protocol['protocol_id']),
            'protocol_label': str(protocol['protocol_label']),
            'armed_topology_class': str(armed['armed_topology_class']),
            'armed_window_length': int(armed['window_length']),
            'final_return_topology_class': final_label,
            'topology_return_error': return_error,
            'return_state_overlap': return_overlap,
            'hysteresis_area_flow': area_between_curves(forward_flow, reverse_flow),
            'hysteresis_area_loop': area_between_curves(forward_loop, reverse_loop),
            'braid_fraction_during_cycle': float(np.mean([1.0 if label == 'braid_like_exchange' else 0.0 for label in topology_trace])) if topology_trace else 0.0,
            'irreversibility_label': irreversibility_label(
                armed_topology=str(armed['armed_topology_class']),
                final_topology=final_label,
                return_error=return_error,
                area_flow=area_between_curves(forward_flow, reverse_flow),
                area_loop=area_between_curves(forward_loop, reverse_loop),
            ),
            'notes': '',
        }
        trace_path = WORK_PLOT_DIR / f'{run_id}_trace.png'
        trace_plot(
            trace_path,
            run_id,
            deltas=delta_trace,
            topology_codes=[mapping.get(label, 2) for label in topology_trace],
            flow=flow_trace,
            loop=loop_trace,
        )
        row['trace_plot'] = str(trace_path)
        rows.append(row)
        plot_paths.append(trace_path)

    summary_path = WORK_PLOT_DIR / 'stage_c0_20_hysteresis_summary.png'
    summary_plot(summary_path, rows)
    plot_paths.append(summary_path)

    result_payload = {
        'stage': 'c0_20_armed_bridge_hysteresis_probe',
        'runsheet': str(Path(args.runsheet).relative_to(REPO_ROOT)),
        'rows': rows,
    }
    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug='stage_c0_20_armed_bridge_hysteresis_probe',
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
