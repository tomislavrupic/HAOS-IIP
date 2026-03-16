#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt, save_atlas_payload
from stage23_1_dk_collision_geometry_scan import mean_pair_separation_series
from stage23_10_dk_braid_mechanism_isolation_probe import analyse_states, packet_state_with_profile
from stage23_9_dk_collision_phenomenology_consolidation import evolve_signed, local_grade_mix
from stage9b_common import build_dk2d_complex, first_order_dt, pack_positions

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage23_13_dk_mechanism_isolation_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_III_C4_Mechanism_Isolation_Probe_v1.md'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage23_13_mechanism_isolation'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'run_id',
    'run_order',
    'branch_id',
    'branch_label',
    'role',
    'ablation_strength',
    'phase_offset_fraction_of_pi',
    'separation',
    'width_ratio_a_to_b',
    'grade_transfer_lock_alpha',
    'operator_01_scale',
    'topology_class',
    'braid_like_exchange',
    'transfer_smeared',
    'localized_encounter',
    'unresolved_mixed',
    'braid_survival_time',
    'grade_exchange_coherence',
    'flow_concentration_index',
    'separation_oscillation_indicator',
    'topology_return_error',
    'mechanism_sensitivity_score',
    'branch_mechanism_label',
    'notes',
]

TOPOLOGY_CODE = {
    'braid_like_exchange': 0,
    'transfer_smeared': 1,
    'localized_capture': 2,
    'dispersive_pass': 3,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run Stage 23.13 DK mechanism isolation probe.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def sigma_pair(base_sigma: float, ratio_a_to_b: float) -> tuple[float, float]:
    root = math.sqrt(max(ratio_a_to_b, 1.0e-12))
    return base_sigma * root, base_sigma / root


def separation_oscillation_indicator(center_histories: list[list[list[float]]]) -> float:
    histories = [np.asarray(track, dtype=float) for track in center_histories]
    mean_pair_distances, _ = mean_pair_separation_series(histories)
    if mean_pair_distances.size < 6:
        return 0.0
    tail = mean_pair_distances[len(mean_pair_distances) // 3:]
    centered = tail - float(np.mean(tail))
    amplitude = float(np.max(np.abs(centered))) if centered.size else 0.0
    if amplitude <= 1.0e-12:
        return 0.0
    zero_crossings = int(np.sum(centered[1:] * centered[:-1] < 0.0))
    return float(zero_crossings / max(1, tail.size - 1))


def topology_survival_time(metrics: dict[str, Any], dt: float) -> float:
    topology = str(metrics['topology_class'])
    trace = [str(label) for label in metrics['topology_trace']]
    return float(sum(label == topology for label in trace) * dt)


def build_operator_scaled_01(base_operator, block_sizes: tuple[int, ...], scale: float):
    if abs(float(scale) - 1.0) <= 1.0e-12:
        return base_operator
    n0, n1, _ = block_sizes
    lil = base_operator.tocsr().astype(complex).tolil(copy=True)
    for row_start, row_stop, col_start, col_stop in (
        (0, n0, n0, n0 + n1),
        (n0, n0 + n1, 0, n0),
    ):
        for row in range(row_start, row_stop):
            cols = lil.rows[row]
            data = lil.data[row]
            for idx, col in enumerate(cols):
                if col_start <= col < col_stop:
                    data[idx] *= float(scale)
    return lil.tocsr()


def evolve_grade_locked(
    base_operator,
    block_sizes: tuple[int, ...],
    state0: np.ndarray,
    dt: float,
    steps: int,
    beta: float,
    lock_alpha: float,
) -> list[np.ndarray]:
    state = state0.copy()
    states = [state.copy()]
    base_norm = float(np.linalg.norm(state0))
    n0, n1, _ = block_sizes
    init_energy = np.asarray(
        [
            float(np.sum(np.abs(state0[:n0]) ** 2)),
            float(np.sum(np.abs(state0[n0:n0 + n1]) ** 2)),
            float(np.sum(np.abs(state0[n0 + n1:]) ** 2)),
        ],
        dtype=float,
    )
    init_fractions = init_energy / max(float(np.sum(init_energy)), 1.0e-12)

    for _ in range(steps):
        acc = -(base_operator @ state)
        if beta > 0.0:
            acc += float(beta) * local_grade_mix(state, block_sizes)
        state = state + dt * acc
        if lock_alpha > 0.0:
            parts = [state[:n0].copy(), state[n0:n0 + n1].copy(), state[n0 + n1:].copy()]
            total = float(np.sum(np.abs(state) ** 2))
            for idx, part in enumerate(parts):
                current = float(np.sum(np.abs(part) ** 2))
                target = (1.0 - float(lock_alpha)) * current + float(lock_alpha) * init_fractions[idx] * total
                if current > 1.0e-12 and target > 0.0:
                    part *= math.sqrt(target / current)
            state = np.concatenate(parts)
        norm = float(np.linalg.norm(state))
        if not np.isfinite(norm) or norm <= 1.0e-12:
            raise RuntimeError('numerical instability indicator triggered')
        state = state / norm * base_norm
        states.append(state.copy())
    return states


def simulate_run(run: dict[str, Any], common: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any], list[Path]]:
    resolution = int(common['base_resolution'])
    epsilon = float(common['base_epsilon'])
    sigma_base = float(common['mean_width'])
    amplitude = 0.5 * float(common['amplitude_scale'])
    t_final = float(common['t_final'])
    beta = float(common['beta'])

    complex_data = build_dk2d_complex(n_side=resolution, epsilon=epsilon)
    positions = pack_positions(complex_data.points, complex_data.edge_midpoints, complex_data.face_centers)
    block_sizes = complex_data.block_sizes
    dt = first_order_dt(complex_data.dirac_kahler, float(common['dt_scale']))
    steps = max(2, int(math.ceil(t_final / dt)))

    sigma_a, sigma_b = sigma_pair(sigma_base, float(run['width_ratio_a_to_b']))
    separation = float(run['separation'])
    centers = [
        np.array([0.50 - 0.5 * separation, 0.50], dtype=float),
        np.array([0.50 + 0.5 * separation, 0.50], dtype=float),
    ]
    kicks = [np.array([0.4, 0.0], dtype=float), np.array([-0.4, 0.0], dtype=float)]

    packet_states = [
        packet_state_with_profile(
            positions=positions,
            block_sizes=block_sizes,
            center=centers[0],
            sigma=sigma_a,
            amplitude=amplitude,
            phase_offset=0.0,
            kick_vector=kicks[0],
            kick_cycles=1.0,
            grade0_scale=1.0,
            grade1_scale=1.0,
            anisotropy=(1.0, 1.0),
        ),
        packet_state_with_profile(
            positions=positions,
            block_sizes=block_sizes,
            center=centers[1],
            sigma=sigma_b,
            amplitude=amplitude,
            phase_offset=float(run['phase_offset_fraction_of_pi']) * math.pi,
            kick_vector=kicks[1],
            kick_cycles=1.0,
            grade0_scale=1.0,
            grade1_scale=1.0,
            anisotropy=(1.0, 1.0),
        ),
    ]

    psi0 = np.sum(packet_states, axis=0)
    close_threshold = 2.0 * float(np.mean([sigma_a, sigma_b]))
    operator = build_operator_scaled_01(
        complex_data.dirac_kahler,
        block_sizes,
        float(run['operator_01_scale']),
    )

    if float(run['grade_transfer_lock_alpha']) > 0.0:
        forward_states = evolve_grade_locked(
            operator,
            block_sizes,
            psi0,
            dt,
            steps,
            beta,
            float(run['grade_transfer_lock_alpha']),
        )
    else:
        forward_states = evolve_signed(operator, block_sizes, psi0, dt, steps, beta, direction=1.0)

    reverse_states = evolve_signed(operator, block_sizes, forward_states[-1], dt, steps, beta, direction=-1.0)

    forward_metrics = analyse_states(
        states=forward_states,
        packet_states=packet_states,
        positions=positions,
        edge_midpoints=complex_data.edge_midpoints,
        block_sizes=block_sizes,
        resolution=resolution,
        anchors=centers,
        close_threshold=close_threshold,
    )
    reverse_metrics = analyse_states(
        states=reverse_states,
        packet_states=packet_states,
        positions=positions,
        edge_midpoints=complex_data.edge_midpoints,
        block_sizes=block_sizes,
        resolution=resolution,
        anchors=centers,
        close_threshold=close_threshold,
    )

    topology_trace = [str(label) for label in forward_metrics['topology_trace']]
    braid_indicator_trace = [1 if label == 'braid_like_exchange' else 0 for label in topology_trace]
    coherence_trace = [float(value) for value in forward_metrics['grade_exchange_trace']]
    times = (np.arange(len(topology_trace), dtype=float) * dt).tolist()
    braid_survival = topology_survival_time(forward_metrics, dt)
    separation_oscillation = separation_oscillation_indicator(forward_metrics['center_histories'])
    topology_return_error = int(str(reverse_metrics['topology_class']) != str(forward_metrics['topology_class']))

    row = {
        'run_id': str(run['run_id']),
        'run_order': int(run['run_order']),
        'branch_id': str(run['branch_id']),
        'branch_label': str(run['branch_label']),
        'role': str(run['role']),
        'ablation_strength': float(run['ablation_strength']),
        'phase_offset_fraction_of_pi': float(run['phase_offset_fraction_of_pi']),
        'separation': float(run['separation']),
        'width_ratio_a_to_b': float(run['width_ratio_a_to_b']),
        'grade_transfer_lock_alpha': float(run['grade_transfer_lock_alpha']),
        'operator_01_scale': float(run['operator_01_scale']),
        'topology_class': str(forward_metrics['topology_class']),
        'braid_like_exchange': int(str(forward_metrics['topology_class']) == 'braid_like_exchange'),
        'transfer_smeared': int(str(forward_metrics['topology_class']) == 'transfer_smeared'),
        'localized_encounter': int(str(forward_metrics['topology_class']) == 'localized_capture'),
        'unresolved_mixed': int(str(forward_metrics['topology_class']) not in {'braid_like_exchange', 'transfer_smeared', 'localized_capture'}),
        'braid_survival_time': float(braid_survival),
        'grade_exchange_coherence': float(forward_metrics['grade_exchange_coherence']),
        'flow_concentration_index': float(forward_metrics['flow_concentration_index']),
        'separation_oscillation_indicator': float(separation_oscillation),
        'topology_return_error': int(topology_return_error),
        'mechanism_sensitivity_score': 0.0,
        'branch_mechanism_label': '',
        'notes': str(run['notes']),
    }
    detail = {
        'run_id': str(run['run_id']),
        'branch_id': str(run['branch_id']),
        'role': str(run['role']),
        'times': times,
        'topology_trace': topology_trace,
        'braid_indicator_trace': braid_indicator_trace,
        'coherence_trace': coherence_trace,
        'ablation_strength_log': {
            'ablation_strength': float(run['ablation_strength']),
            'phase_offset_fraction_of_pi': float(run['phase_offset_fraction_of_pi']),
            'phase_shift_from_control': float(run['phase_offset_fraction_of_pi']) - float(common['reference_phase_offset_fraction_of_pi']),
            'separation': float(run['separation']),
            'separation_delta': float(run['separation']) - float(common['separation']),
            'width_ratio_a_to_b': float(run['width_ratio_a_to_b']),
            'grade_transfer_lock_alpha': float(run['grade_transfer_lock_alpha']),
            'operator_01_scale': float(run['operator_01_scale']),
        },
        'forward_metrics': forward_metrics,
        'reverse_metrics': reverse_metrics,
        'dt': float(dt),
        'steps': int(steps),
    }

    trace_path = WORK_PLOT_DIR / f"{run['run_id']}_trace_panel.png"
    fig, axes = plt.subplots(3, 1, figsize=(6.8, 7.8), sharex=True)
    axes[0].plot(times, [TOPOLOGY_CODE[label] for label in topology_trace], color='tab:blue')
    axes[0].set_ylabel('topology')
    axes[0].set_yticks([0, 1, 2, 3], ['braid', 'smeared', 'capture', 'mixed'])
    axes[0].set_title(f"{run['run_id']} | {run['branch_label']} | {run['role']}")
    axes[0].grid(alpha=0.25)
    axes[1].plot(times, braid_indicator_trace, color='tab:orange')
    axes[1].set_ylabel('braid')
    axes[1].set_yticks([0, 1], ['off', 'on'])
    axes[1].grid(alpha=0.25)
    axes[2].plot(times, coherence_trace, color='tab:green')
    axes[2].set_ylabel('coherence')
    axes[2].set_xlabel('time')
    axes[2].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(trace_path, dpi=180, bbox_inches='tight')
    plt.close(fig)

    return row, detail, [trace_path]


def mechanism_sensitivity_score(baseline: dict[str, Any], row: dict[str, Any]) -> float:
    if str(row['role']) == 'control':
        return 0.0
    topology_term = 0.0
    if str(row['topology_class']) != str(baseline['topology_class']):
        if str(row['topology_class']) == 'transfer_smeared':
            topology_term = 1.0
        elif str(row['topology_class']) == 'localized_capture':
            topology_term = 0.8
        else:
            topology_term = 0.6
    baseline_survival = max(float(baseline['braid_survival_time']), 1.0e-12)
    survival_drop = float(np.clip((float(baseline['braid_survival_time']) - float(row['braid_survival_time'])) / baseline_survival, 0.0, 1.0))
    baseline_coherence = max(float(baseline['grade_exchange_coherence']), 1.0e-12)
    coherence_drop = float(np.clip((float(baseline['grade_exchange_coherence']) - float(row['grade_exchange_coherence'])) / baseline_coherence, 0.0, 1.0))
    return float(np.clip(0.55 * topology_term + 0.30 * survival_drop + 0.15 * coherence_drop, 0.0, 1.0))


def meaningful_collapse(row: dict[str, Any]) -> bool:
    return (
        str(row['topology_class']) != 'braid_like_exchange'
        and float(row['flow_concentration_index']) >= 0.82
        and float(row['grade_exchange_coherence']) >= 0.35
    )


def primary_label_for_branch(branch_id: str) -> str:
    mapping = {
        'A_grade_transfer_suppression': 'grade_transfer_primary',
        'B_phase_corridor_flattening': 'phase_corridor_primary',
        'C_geometry_overlap_weakening': 'geometry_overlap_primary',
        'D_operator_cross_coupling_suppression': 'operator_coupling_primary',
    }
    return mapping[str(branch_id)]


def classify_branch(rows: list[dict[str, Any]]) -> tuple[str, dict[str, Any]]:
    ordered = sorted(
        rows,
        key=lambda row: {'control': 0, 'weak_suppression': 1, 'moderate_suppression': 2}[str(row['role'])],
    )
    baseline, weak, moderate = ordered
    weak_collapse = meaningful_collapse(weak)
    moderate_collapse = meaningful_collapse(moderate)
    weak_retains = str(weak['topology_class']) == 'braid_like_exchange'
    moderate_retains = str(moderate['topology_class']) == 'braid_like_exchange'
    if weak_collapse and moderate_collapse:
        label = primary_label_for_branch(str(baseline['branch_id']))
    elif weak_retains and moderate_collapse:
        label = 'coupled_mechanism_no_single_primary'
    elif weak_retains and moderate_retains:
        label = 'no_clear_mechanism_signal'
    else:
        label = 'no_clear_mechanism_signal'
    summary = {
        'weak_collapse': bool(weak_collapse),
        'moderate_collapse': bool(moderate_collapse),
        'mean_sensitivity_score': float(np.mean([float(weak['mechanism_sensitivity_score']), float(moderate['mechanism_sensitivity_score'])])),
        'max_sensitivity_score': float(max(float(weak['mechanism_sensitivity_score']), float(moderate['mechanism_sensitivity_score']))),
        'control_topology': str(baseline['topology_class']),
        'weak_topology': str(weak['topology_class']),
        'moderate_topology': str(moderate['topology_class']),
    }
    return label, summary


def classify_overall(branch_payloads: dict[str, dict[str, Any]]) -> tuple[str, str]:
    def join_names(names: list[str]) -> str:
        if len(names) == 1:
            return names[0]
        if len(names) == 2:
            return f'{names[0]} and {names[1]}'
        return ', '.join(names[:-1]) + f', and {names[-1]}'

    primary_labels = {
        'grade_transfer_primary',
        'phase_corridor_primary',
        'geometry_overlap_primary',
        'operator_coupling_primary',
    }
    primary_ids = [
        branch_id
        for branch_id, payload in branch_payloads.items()
        if str(payload['label']) in primary_labels
    ]
    coupled_ids = [
        branch_id
        for branch_id, payload in branch_payloads.items()
        if str(payload['label']) == 'coupled_mechanism_no_single_primary'
    ]
    no_clear_ids = [
        branch_id
        for branch_id, payload in branch_payloads.items()
        if str(payload['label']) == 'no_clear_mechanism_signal'
    ]

    if len(primary_ids) == 1:
        primary = primary_ids[0]
        primary_name = {
            'A_grade_transfer_suppression': 'grade transfer',
            'B_phase_corridor_flattening': 'phase corridor',
            'C_geometry_overlap_weakening': 'geometry overlap',
            'D_operator_cross_coupling_suppression': 'operator cross-coupling',
        }[primary]
        secondary_names = [
            {
                'A_grade_transfer_suppression': 'grade transfer',
                'B_phase_corridor_flattening': 'phase corridor',
                'C_geometry_overlap_weakening': 'geometry overlap',
                'D_operator_cross_coupling_suppression': 'operator cross-coupling',
            }[branch_id]
            for branch_id in coupled_ids
        ]
        if secondary_names:
            statement = (
                f"braid-like exchange in the DK sector is primarily controlled by {primary_name}, "
                f"with {join_names(secondary_names)} acting as secondary support"
            )
        else:
            statement = f"braid-like exchange in the DK sector is primarily controlled by {primary_name}"
        return str(branch_payloads[primary]['label']), statement

    if len(primary_ids) >= 2 or coupled_ids:
        return 'coupled_mechanism_no_single_primary', 'the DK braid / smear split remains coupled at the present reduction level'

    if no_clear_ids:
        return 'no_clear_mechanism_signal', 'no single mechanism branch produced a discriminative collapse signal'

    return 'no_clear_mechanism_signal', 'no single mechanism branch produced a discriminative collapse signal'


def plot_branch_comparison_matrix(path: Path, rows: list[dict[str, Any]]) -> None:
    branch_order = [
        'A_grade_transfer_suppression',
        'B_phase_corridor_flattening',
        'C_geometry_overlap_weakening',
        'D_operator_cross_coupling_suppression',
    ]
    role_order = ['control', 'weak_suppression', 'moderate_suppression']
    label_map = {
        'A_grade_transfer_suppression': 'grade transfer',
        'B_phase_corridor_flattening': 'phase corridor',
        'C_geometry_overlap_weakening': 'geometry overlap',
        'D_operator_cross_coupling_suppression': 'operator coupling',
    }
    grid = np.zeros((len(branch_order), len(role_order)), dtype=float)
    annotations: list[list[str]] = [['' for _ in role_order] for _ in branch_order]
    short = {
        'braid_like_exchange': 'braid',
        'transfer_smeared': 'smear',
        'localized_capture': 'capture',
        'dispersive_pass': 'mixed',
    }
    for i, branch_id in enumerate(branch_order):
        for j, role in enumerate(role_order):
            row = next(item for item in rows if str(item['branch_id']) == branch_id and str(item['role']) == role)
            grid[i, j] = float(row['mechanism_sensitivity_score'])
            annotations[i][j] = short[str(row['topology_class'])]

    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    im = ax.imshow(grid, aspect='auto', cmap='viridis', vmin=0.0, vmax=1.0)
    ax.set_xticks(range(len(role_order)), ['control', 'weak', 'moderate'])
    ax.set_yticks(range(len(branch_order)), [label_map[item] for item in branch_order])
    ax.set_title('Stage III-C4 mechanism branch comparison matrix')
    for i in range(len(branch_order)):
        for j in range(len(role_order)):
            ax.text(j, i, annotations[i][j], ha='center', va='center', color='white', fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_braid_survival_panel(path: Path, rows: list[dict[str, Any]]) -> None:
    branch_order = [
        'A_grade_transfer_suppression',
        'B_phase_corridor_flattening',
        'C_geometry_overlap_weakening',
        'D_operator_cross_coupling_suppression',
    ]
    label_map = {
        'A_grade_transfer_suppression': 'grade transfer',
        'B_phase_corridor_flattening': 'phase corridor',
        'C_geometry_overlap_weakening': 'geometry overlap',
        'D_operator_cross_coupling_suppression': 'operator coupling',
    }
    role_order = ['control', 'weak_suppression', 'moderate_suppression']
    x = np.arange(len(role_order), dtype=float)
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    for branch_id in branch_order:
        vals = [
            float(next(item for item in rows if str(item['branch_id']) == branch_id and str(item['role']) == role)['braid_survival_time'])
            for role in role_order
        ]
        ax.plot(x, vals, marker='o', linewidth=2.0, label=label_map[branch_id])
    ax.set_xticks(x, ['control', 'weak', 'moderate'])
    ax.set_ylabel('braid survival time')
    ax.set_title('Stage III-C4 braid survival vs ablation')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_coherence_collapse_panel(path: Path, rows: list[dict[str, Any]]) -> None:
    branch_order = [
        'A_grade_transfer_suppression',
        'B_phase_corridor_flattening',
        'C_geometry_overlap_weakening',
        'D_operator_cross_coupling_suppression',
    ]
    label_map = {
        'A_grade_transfer_suppression': 'grade transfer',
        'B_phase_corridor_flattening': 'phase corridor',
        'C_geometry_overlap_weakening': 'geometry overlap',
        'D_operator_cross_coupling_suppression': 'operator coupling',
    }
    role_order = ['control', 'weak_suppression', 'moderate_suppression']
    x = np.arange(len(role_order), dtype=float)
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    for branch_id in branch_order:
        vals = [
            float(next(item for item in rows if str(item['branch_id']) == branch_id and str(item['role']) == role)['grade_exchange_coherence'])
            for role in role_order
        ]
        ax.plot(x, vals, marker='o', linewidth=2.0, label=label_map[branch_id])
    ax.set_xticks(x, ['control', 'weak', 'moderate'])
    ax.set_ylabel('grade-exchange coherence')
    ax.set_title('Stage III-C4 coherence collapse panel')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(
    json_path: Path,
    csv_path: Path,
    rows: list[dict[str, Any]],
    stamped_plots: list[str],
    branch_payloads: dict[str, dict[str, Any]],
    overall_label: str,
    overall_statement: str,
) -> None:
    lines = [
        '# Stage III-C4 Mechanism Isolation Probe v1',
        '',
        'This mechanism probe stays inside the frozen HAOS (Harmonic Address Operating System) / IIP (Interaction Invariance Physics) Dirac-Kahler clustered collision sector.',
        '',
        f"Timestamped JSON: `{json_path.relative_to(REPO_ROOT)}`",
        f"Timestamped CSV: `{csv_path.relative_to(REPO_ROOT)}`",
        '',
        f"overall_mechanism_label = `{overall_label}`",
        '',
        '## Assumption note',
        'The grade-transfer branch is implemented as bounded blockwise grade-fraction anchoring, because the frozen DK branch does not expose a literal transfer-off switch without changing operator class. The operator-coupling branch scales only the `0 <-> 1` DK blocks and leaves the rest of the propagator intact.',
        '',
        '## Branch outcomes',
    ]
    for branch_id in [
        'A_grade_transfer_suppression',
        'B_phase_corridor_flattening',
        'C_geometry_overlap_weakening',
        'D_operator_cross_coupling_suppression',
    ]:
        payload = branch_payloads[branch_id]
        summary = payload['summary']
        lines.append(
            f"- `{payload['branch_label']}`: label=`{payload['label']}`, "
            f"control=`{summary['control_topology']}`, weak=`{summary['weak_topology']}`, moderate=`{summary['moderate_topology']}`, "
            f"mean_sensitivity=`{summary['mean_sensitivity_score']:.3f}`"
        )

    lines.extend([
        '',
        '## Per-run summary',
    ])
    for row in sorted(rows, key=lambda item: int(item['run_order'])):
        lines.append(
            f"- `{row['run_id']}`: topology=`{row['topology_class']}`, braid_survival=`{row['braid_survival_time']:.3f}`, "
            f"coherence=`{row['grade_exchange_coherence']:.3f}`, flow=`{row['flow_concentration_index']:.3f}`, "
            f"sensitivity=`{row['mechanism_sensitivity_score']:.3f}`"
        )

    lines.extend([
        '',
        '## Mechanism conclusion',
    ])
    if overall_label == 'phase_corridor_primary':
        lines.append(
            'The strongest discriminative collapse comes from phase-corridor flattening: both weak and moderate phase flattening move the clustered control out of the braid class while the control remains braid-like and the run stays bounded rather than dissolving into noise.'
        )
    elif overall_label == 'geometry_overlap_primary':
        lines.append(
            'The strongest discriminative collapse comes from geometry-overlap weakening: braid is lost under controlled overlap reduction while the other branches retain the braid control more effectively.'
        )
    elif overall_label == 'operator_coupling_primary':
        lines.append(
            'The strongest discriminative collapse comes from suppressing the `0 <-> 1` operator blocks: braid depends primarily on the mixed DK operator structure rather than only on seed geometry or phase.'
        )
    elif overall_label == 'grade_transfer_primary':
        lines.append(
            'The strongest discriminative collapse comes from suppressing net cross-grade redistribution: braid depends primarily on active grade transfer rather than only on phase or geometry.'
        )
    elif overall_label == 'coupled_mechanism_no_single_primary':
        lines.append(
            'No single ablation dominates cleanly. The braid / smear split remains coupled at the present reduction level, so the correct read is a multi-factor mechanism rather than a single isolated driver.'
        )
    else:
        lines.append(
            'No branch produced a clean discriminative mechanism signal beyond generic weakening. At this level the Phase III result remains phenomenological rather than mechanistically isolated.'
        )
    lines.append(f'- Conclusion statement: `{overall_statement}`')

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
    common = runsheet['common_fields']

    rows: list[dict[str, Any]] = []
    details: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    by_branch: dict[str, list[dict[str, Any]]] = defaultdict(list)

    for run in sorted(runsheet['runs'], key=lambda item: int(item['run_order'])):
        row, detail, run_plots = simulate_run(run, common)
        rows.append(row)
        details.append(detail)
        plot_paths.extend(run_plots)
        by_branch[str(row['branch_id'])].append(row)

    branch_payloads: dict[str, dict[str, Any]] = {}
    for branch_id, branch_rows in by_branch.items():
        control = next(row for row in branch_rows if str(row['role']) == 'control')
        for row in branch_rows:
            row['mechanism_sensitivity_score'] = mechanism_sensitivity_score(control, row)
        label, summary = classify_branch(branch_rows)
        for row in branch_rows:
            row['branch_mechanism_label'] = label
        branch_payloads[branch_id] = {
            'branch_id': branch_id,
            'branch_label': str(control['branch_label']),
            'label': label,
            'summary': summary,
        }

    overall_label, overall_statement = classify_overall(branch_payloads)

    matrix_path = WORK_PLOT_DIR / 'stage23_13_mechanism_branch_comparison_matrix.png'
    braid_panel_path = WORK_PLOT_DIR / 'stage23_13_braid_survival_vs_ablation_panel.png'
    coherence_panel_path = WORK_PLOT_DIR / 'stage23_13_coherence_collapse_panel.png'
    plot_branch_comparison_matrix(matrix_path, rows)
    plot_braid_survival_panel(braid_panel_path, rows)
    plot_coherence_collapse_panel(coherence_panel_path, rows)
    plot_paths.extend([matrix_path, braid_panel_path, coherence_panel_path])

    result = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'common_fields': common,
        'runs': details,
        'branch_summary': branch_payloads,
        'overall_summary': {
            'overall_mechanism_label': overall_label,
            'overall_statement': overall_statement,
        },
    }
    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        'stage23_13_dk_mechanism_isolation_probe',
        result,
        rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots, branch_payloads, overall_label, overall_statement)


if __name__ == '__main__':
    main()
