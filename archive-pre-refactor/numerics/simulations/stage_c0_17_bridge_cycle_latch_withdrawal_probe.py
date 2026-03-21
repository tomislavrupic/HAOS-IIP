#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt, save_atlas_payload
from stage_c0_dk_bridge_validation import (
    RUNSHEET_PATH as BRIDGE_RUNSHEET_PATH,
    address_activity_statistics,
    address_weighted_operator,
    block_dirac,
    build_address_labels,
    build_d0,
    build_d1,
    classify_relative_topology,
    edge_component_count,
    estimate_lambda_max,
    grade_exchange_signal,
    grade_histories,
    initial_packet,
    mismatch_scores,
    orient_edges,
    representative_graph,
    triangle_list,
    triangle_loop_metrics,
)


RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage_c0_17_bridge_cycle_latch_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_C0_17_Bridge_Cycle_Latch_Withdrawal_Probe_v1.md'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage_c0_17_bridge_cycle_latch'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'run_id',
    'family_id',
    'family_label',
    'latch_id',
    'latch_label',
    'armed_topology_class',
    'post_withdrawal_topology_class',
    'withdrawal_braid_dwell_time',
    'loop_retention_score',
    'path_reuse_retention_score',
    'latch_duty_cycle',
    'flow_concentration_index_post',
    'grade_exchange_coherence_post',
    'selector_off_recurrence_indicator',
    'locking_label',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run the Stage C0.17 bridge cycle-latch withdrawal probe.')
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
            run_id = f"C017_{family['family_id']}_{latch['latch_id']}"
            pairs.append((run_id, family, latch))
    if not run_ids:
        return pairs
    wanted = set(run_ids)
    return [item for item in pairs if item[0] in wanted]


def bridge_setup(family: dict[str, Any], common: dict[str, Any], stage_common: dict[str, Any]) -> dict[str, Any]:
    family_id = str(family['family_id'])
    left_seed = int(family['left_seed'])
    right_seed = int(family['right_seed'])
    G = representative_graph(family_id)
    nodes = sorted(int(node) for node in G.nodes())
    scores = mismatch_scores(G, left_seed=left_seed, right_seed=right_seed)
    oriented_edges = orient_edges(G, scores)
    triangles = triangle_list(G)
    d0, edge_lookup, _node_to_idx = build_d0(nodes, oriented_edges)
    d1 = build_d1(triangles, edge_lookup) if triangles else sp.csr_matrix((0, len(oriented_edges)), dtype=float)
    D, _Delta = block_dirac(d0, d1)
    sigma_shell = float(common['sigma_shell'])
    kick_cycles = float(common['kick_cycles'])
    packet_a = initial_packet(
        G=G,
        nodes=nodes,
        oriented_edges=oriented_edges,
        triangles=triangles,
        seed=left_seed,
        sigma_shell=sigma_shell,
        kick_cycles=kick_cycles,
        direction_sign=+1.0,
        grade0_scale=float(common['grade0_scale']),
        grade1_scale=float(common['grade1_scale']),
        grade2_scale=float(common['grade2_scale']),
    )
    packet_b = initial_packet(
        G=G,
        nodes=nodes,
        oriented_edges=oriented_edges,
        triangles=triangles,
        seed=right_seed,
        sigma_shell=sigma_shell,
        kick_cycles=kick_cycles,
        direction_sign=-1.0,
        grade0_scale=float(common['grade0_scale']),
        grade1_scale=float(common['grade1_scale']),
        grade2_scale=float(common['grade2_scale']),
    )
    block_sizes = (len(nodes), len(oriented_edges), len(triangles))
    exact_labels = build_address_labels(
        packet_states=[packet_a, packet_b],
        block_sizes=block_sizes,
        left_address=0.0,
        right_address=0.0,
        support_floor_fraction=float(common['support_floor_fraction']),
        dominant_ratio_threshold=float(common['dominant_ratio_threshold']),
    )
    withdrawal_delta = float(stage_common['withdrawal_delta'])
    withdrawal_labels = build_address_labels(
        packet_states=[packet_a, packet_b],
        block_sizes=block_sizes,
        left_address=0.0,
        right_address=2.0 * withdrawal_delta,
        support_floor_fraction=float(common['support_floor_fraction']),
        dominant_ratio_threshold=float(common['dominant_ratio_threshold']),
    )
    exact_operator = address_weighted_operator(
        base_operator=D,
        labels=exact_labels,
        modulus=float(common['address_modulus']),
        eta_match=float(common['eta_match']),
        eta_mismatch=float(common['eta_mismatch']),
    )
    withdrawal_operator = address_weighted_operator(
        base_operator=D,
        labels=withdrawal_labels,
        modulus=float(common['address_modulus']),
        eta_match=float(common['eta_match']),
        eta_mismatch=float(common['eta_mismatch']),
    )
    psi0 = packet_a + packet_b
    psi0 /= max(float(np.linalg.norm(psi0)), 1.0e-12)
    return {
        'family_id': family_id,
        'family_label': str(family['family_label']),
        'nodes': nodes,
        'oriented_edges': oriented_edges,
        'triangles': triangles,
        'd1': d1,
        'block_sizes': block_sizes,
        'exact_operator': exact_operator,
        'withdrawal_operator': withdrawal_operator,
        'packet_states': [packet_a, packet_b],
        'exact_labels': exact_labels,
        'withdrawal_labels': withdrawal_labels,
        'psi0': psi0,
    }


def apply_basis_latch(
    operator: sp.csr_matrix,
    basis_weights: np.ndarray,
    amplitude: float,
) -> sp.csr_matrix:
    if amplitude <= 0.0 or basis_weights.size == 0:
        return operator
    weights = np.asarray(basis_weights, dtype=float)
    peak = float(np.max(weights))
    if peak <= 1.0e-12:
        return operator
    weights = weights / peak
    coo = operator.tocoo()
    mult = 1.0 + float(amplitude) * 0.5 * (weights[coo.row] + weights[coo.col])
    return sp.csr_matrix((coo.data * mult, (coo.row, coo.col)), shape=operator.shape).tocsr()


def one_step_cn(operator: sp.csr_matrix, state: np.ndarray, dt: float) -> np.ndarray:
    ident = sp.identity(operator.shape[0], dtype=complex, format='csr')
    lhs = (ident + 0.5j * dt * operator).tocsc()
    rhs = (ident - 0.5j * dt * operator) @ state
    next_state = spla.spsolve(lhs, rhs)
    next_state = np.asarray(next_state, dtype=complex)
    next_state /= max(float(np.linalg.norm(next_state)), 1.0e-12)
    return next_state


def latch_strength_for_rule(latch_id: str, common: dict[str, Any]) -> float:
    if latch_id == 'cycle_occupancy_latch':
        return float(common['cycle_latch_strength'])
    if latch_id == 'edge_path_reuse_latch':
        return float(common['path_latch_strength'])
    return float(common['combined_latch_strength'])


def build_latch_weights(
    latch_id: str,
    state: np.ndarray,
    accum_edge: np.ndarray,
    accum_face: np.ndarray,
    block_sizes: tuple[int, int, int],
) -> np.ndarray:
    n0, n1, n2 = block_sizes
    edge_power = np.abs(state[n0:n0 + n1]) ** 2
    face_power = np.abs(state[n0 + n1:n0 + n1 + n2]) ** 2 if n2 > 0 else np.zeros(0, dtype=float)
    accum_edge += edge_power
    if n2 > 0:
        accum_face += face_power

    weights = np.zeros(sum(block_sizes), dtype=float)
    if latch_id == 'cycle_occupancy_latch':
        if n2 > 0:
            weights[n0 + n1:n0 + n1 + n2] = accum_face
    elif latch_id == 'edge_path_reuse_latch':
        weights[n0:n0 + n1] = accum_edge
    else:
        weights[n0:n0 + n1] = accum_edge
        if n2 > 0:
            weights[n0 + n1:n0 + n1 + n2] = accum_face
    return weights


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


def window_metrics(
    states: list[np.ndarray],
    bridge: dict[str, Any],
    labels: np.ndarray,
    operator: sp.csr_matrix,
    baseline_result: dict[str, Any],
    common: dict[str, Any],
) -> dict[str, Any]:
    n0, n1, _n2 = bridge['block_sizes']
    edge_series = np.asarray([state[n0:n0 + n1] for state in states], dtype=complex)
    grade_hist = grade_histories(states, bridge['block_sizes'])
    _signal, coherence = grade_exchange_signal(grade_hist)
    edge_power = np.mean(np.abs(edge_series) ** 2, axis=0) if edge_series.size else np.zeros(n1, dtype=float)
    total_edge_power = float(np.sum(edge_power))
    if total_edge_power > 1.0e-12:
        flat = np.sort(edge_power)[::-1]
        top_k = max(1, int(math.ceil(0.1 * flat.size)))
        concentration = float(np.sum(flat[:top_k]) / total_edge_power)
        threshold = max(0.12 * float(np.max(edge_power)), float(np.mean(edge_power) + 0.50 * np.std(edge_power)))
        active_mask = edge_power >= threshold
    else:
        concentration = 0.0
        active_mask = np.zeros(n1, dtype=bool)
    channel_count = int(edge_component_count(bridge['oriented_edges'], active_mask))
    loop_count, loop_score = triangle_loop_metrics(edge_series, bridge['triangles'], bridge['d1'])
    sorting, selectivity = address_activity_statistics(
        states=states,
        operator=operator,
        labels=labels,
        modulus=float(common['address_modulus']),
        eta_match=float(common['eta_match']),
        eta_mismatch=float(common['eta_mismatch']),
    )
    row = {
        'flow_concentration_index': concentration,
        'grade_exchange_coherence': coherence,
        'address_selectivity_index': selectivity,
        'address_sorting_score': sorting,
        'loop_count': loop_count,
        'loop_score': loop_score,
        'channel_count': channel_count,
        'grade_asymmetry_index': float(abs(grade_hist[-1, 0] - grade_hist[-1, 1])) if grade_hist.shape[1] >= 2 else 0.0,
    }
    row['topology_class'] = classify_relative_topology(row, baseline_result)
    return row


def dwell_time_from_trace(labels: list[str], target: str) -> float:
    count = 0
    for label in labels:
        if label == target:
            count += 1
        else:
            break
    return float(count)


def retention_label(
    armed_topology: str,
    post_topology: str,
    dwell: float,
    loop_retention_score: float,
    path_retention_score: float,
) -> str:
    if armed_topology == 'braid_like_exchange' and post_topology == 'braid_like_exchange' and dwell >= 8.0:
        return 'locked_after_withdrawal'
    if armed_topology == 'braid_like_exchange' and dwell >= 3.0:
        return 'quasi_locked_retention'
    if armed_topology == 'braid_like_exchange' and (loop_retention_score >= 0.75 or path_retention_score >= 0.75):
        return 'topology_afterglow'
    if armed_topology == 'braid_like_exchange':
        return 'selector_dependent_collapse'
    return 'no_locking_after_withdrawal'


def trace_plot(path: Path, run_id: str, topology_trace: list[str], latch_trace: list[float]) -> None:
    mapping = {'braid_like_exchange': 0, 'transfer_smeared': 1, 'unresolved_mixed': 2}
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    axes[0].plot(range(len(topology_trace)), [mapping[label] for label in topology_trace], marker='o')
    axes[0].set_yticks([0, 1, 2])
    axes[0].set_yticklabels(['braid', 'smeared', 'mixed'])
    axes[0].set_title('Topology trace')
    axes[0].grid(alpha=0.25)
    axes[1].plot(range(len(latch_trace)), latch_trace, marker='o', color='tab:red')
    axes[1].set_title('Latch duty trace')
    axes[1].grid(alpha=0.25)
    fig.suptitle(run_id)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def summary_plot(path: Path, rows: list[dict[str, Any]]) -> None:
    x = np.arange(len(rows))
    dwell = [float(row['withdrawal_braid_dwell_time']) for row in rows]
    loop_retention = [float(row['loop_retention_score']) for row in rows]
    path_retention = [float(row['path_reuse_retention_score']) for row in rows]
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.5))
    axes[0].bar(x, dwell, color='tab:blue')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([str(row['run_id']) for row in rows], rotation=30, ha='right')
    axes[0].set_title('Withdrawal braid dwell')
    axes[0].grid(alpha=0.25, axis='y')
    axes[1].plot(x, loop_retention, marker='o', label='loop retention')
    axes[1].plot(x, path_retention, marker='s', label='path retention')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([str(row['run_id']) for row in rows], rotation=30, ha='right')
    axes[1].set_title('Retention scores')
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.25)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


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
    upper_operator = apply_basis_latch(bridge['exact_operator'], np.ones_like(bridge['psi0'], dtype=float), latch_strength)
    upper_operator = apply_basis_latch(bridge['withdrawal_operator'], np.ones_like(bridge['psi0'], dtype=float), latch_strength)
    lam_upper = max(estimate_lambda_max(upper_operator), 1.0e-12)
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
        active_operator = apply_basis_latch(bridge['exact_operator'], frozen_weights, latch_strength)
        state = one_step_cn(active_operator, state, dt)
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


def note_text(json_rel: str, csv_rel: str, rows: list[dict[str, Any]], plots: list[str]) -> str:
    counts = Counter(str(row['locking_label']) for row in rows)
    global_read = 'no bridge-local latch retains topology after selector withdrawal' if counts.get('locked_after_withdrawal', 0) == 0 and counts.get('quasi_locked_retention', 0) == 0 else 'bridge-local latching opens partial retention after selector withdrawal'
    lines = [
        '# Stage C0.17 Bridge Cycle-Latch Withdrawal Probe v1',
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

    summary_path = WORK_PLOT_DIR / 'stage_c0_17_withdrawal_summary.png'
    summary_plot(summary_path, rows)
    plot_paths.append(summary_path)

    payload = {
        'stage': runsheet['stage'],
        'description': runsheet['description'],
        'rows': [{key: value for key, value in row.items() if key != 'trace_plot'} for row in rows],
    }
    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug='stage_c0_17_bridge_cycle_latch_withdrawal_probe',
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
