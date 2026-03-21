#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import scipy.sparse as sp

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt, save_atlas_payload
from stage_11_integrated_austerity_resolution_probe import (
    basis_supports,
    field_update,
    initial_packet,
    matched_family_graph,
    node_activity_from_state,
    node_kernel_matrix,
    node_triangle_counts,
    motif_nodes_for_family,
    orient_edges_from_mismatch,
    support_distance_vector,
    triangle_list,
    radius_operator,
)
from stage_c0_dk_bridge_validation import (
    block_dirac,
    build_d0,
    build_d1,
    edge_component_count,
    grade_exchange_signal,
    grade_histories,
    triangle_loop_metrics,
)


RUNSHEET_PATH = (
    REPO_ROOT
    / 'numerics'
    / 'simulations'
    / 'stage_12_intrinsic_grading_and_structural_irreversibility_runs.json'
)
NOTE_PATH = ATLAS_NOTES / 'Stage_12_Intrinsic_Grading_and_Structural_Irreversibility_v1.md'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage_12_intrinsic_grading'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'run_id',
    'family_id',
    'family_label',
    'protocol_id',
    'protocol_label',
    'resolution_id',
    'resolution_label',
    'dominant_forward_topology_class',
    'return_checkpoint_topology_class',
    'final_topology_class',
    'survival_time',
    'protected_fraction',
    'topology_return_error',
    'mismatch_ledger_alignment',
    'signed_field_flux_projection',
    'family_ordering_score',
    'structural_irreversibility_score',
    'recovery_entropy',
    'sector_separation',
    'detuning_monotonicity',
    'refinement_drift',
    'trace_plot',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run Stage 12 intrinsic grading and structural irreversibility probe.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    parser.add_argument('--run-ids', nargs='*', default=[])
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def classify_topology_label_free(
    flow_concentration: float,
    coherence: float,
    loop_score: float,
    channel_count: int,
    capture_localization: float,
) -> str:
    if loop_score >= 0.16 and flow_concentration >= 0.22 and coherence >= 0.38 and channel_count <= 2:
        return 'braid_like_exchange'
    if capture_localization >= 0.36 and channel_count <= 1:
        return 'localized_capture'
    if flow_concentration >= 0.17 and coherence >= 0.24:
        return 'transfer_smeared'
    if loop_score >= 0.08 or channel_count >= 2:
        return 'unresolved_mixed'
    return 'no_protected_structure'


def metrics_for_state_label_free(
    states_window: list[np.ndarray],
    nodes: list[int],
    oriented_edges: list[tuple[int, int]],
    triangles: list[tuple[int, int, int]],
    d1: sp.csr_matrix,
    block_sizes: tuple[int, int, int],
) -> dict[str, float | int | str]:
    n0, n1, _n2 = block_sizes
    edge_series = np.asarray([state[n0:n0 + n1] for state in states_window], dtype=complex)
    grade_hist = grade_histories(states_window, block_sizes)
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
    channel_count = int(edge_component_count(oriented_edges, active_mask))
    loop_count, loop_score = triangle_loop_metrics(edge_series, triangles, d1)
    node_power = np.abs(states_window[-1][:n0]) ** 2
    capture_localization = float(np.max(node_power) / max(float(np.sum(node_power)), 1.0e-12))
    topology = classify_topology_label_free(
        flow_concentration=concentration,
        coherence=float(coherence),
        loop_score=float(loop_score),
        channel_count=channel_count,
        capture_localization=capture_localization,
    )
    return {
        'topology_class': topology,
        'flow_concentration_index': concentration,
        'grade_exchange_coherence': float(coherence),
        'loop_count': int(loop_count),
        'loop_score': float(loop_score),
        'channel_count': int(channel_count),
        'capture_localization_index': capture_localization,
    }


def correlation_or_zero(a: list[float], b: list[float]) -> float:
    if len(a) < 3 or len(b) < 3:
        return 0.0
    arr_a = np.asarray(a, dtype=float)
    arr_b = np.asarray(b, dtype=float)
    if float(np.std(arr_a)) <= 1.0e-12 or float(np.std(arr_b)) <= 1.0e-12:
        return 0.0
    return float(np.corrcoef(arr_a, arr_b)[0, 1])


def conditional_entropy(labels: list[str], bands: list[str]) -> float:
    if not labels or not bands or len(labels) != len(bands):
        return 0.0
    total = len(labels)
    grouped: dict[str, list[str]] = defaultdict(list)
    for label, band in zip(labels, bands):
        grouped[str(band)].append(str(label))
    entropy = 0.0
    for values in grouped.values():
        weight = len(values) / total
        counts = Counter(values)
        probs = np.asarray([count / len(values) for count in counts.values()], dtype=float)
        band_entropy = -float(np.sum(probs * np.log(np.maximum(probs, 1.0e-12)))) / max(math.log(5.0), 1.0)
        entropy += weight * band_entropy
    return float(entropy)


def mode_or_default(values: list[str], default: str) -> str:
    if not values:
        return default
    return Counter(values).most_common(1)[0][0]


def area_between_curves(forward: list[float], reverse: list[float]) -> float:
    if not forward or not reverse:
        return 0.0
    n = min(len(forward), len(reverse))
    if n == 0:
        return 0.0
    f = np.asarray(forward[:n], dtype=float)
    r = np.asarray(list(reversed(reverse[-n:])), dtype=float)
    return float(np.mean(np.abs(f - r)))


def recovered_grading_band(
    states_window: list[np.ndarray],
    nodes: list[int],
    oriented_edges: list[tuple[int, int]],
    triangles: list[tuple[int, int, int]],
) -> tuple[str, float]:
    if len(states_window) < 2:
        return 'mixed', 0.0
    activity = np.asarray(
        [
            node_activity_from_state(
                psi=state,
                nodes=nodes,
                oriented_edges=oriented_edges,
                triangles=triangles,
            )
            for state in states_window
        ],
        dtype=float,
    )
    if activity.shape[1] < 2:
        return 'mixed', 0.0
    centered = activity - np.mean(activity, axis=0, keepdims=True)
    cov = centered.T @ centered
    evals, evecs = np.linalg.eigh(cov)
    lead = np.asarray(evecs[:, int(np.argmax(evals))], dtype=float)
    pivot = float(np.median(lead))
    pos_mask = lead >= pivot
    neg_mask = ~pos_mask
    if not np.any(pos_mask) or not np.any(neg_mask):
        return 'mixed', 0.0
    latest = activity[-1]
    total = max(float(np.sum(latest)), 1.0e-12)
    pos_power = float(np.sum(latest[pos_mask])) / total
    neg_power = float(np.sum(latest[neg_mask])) / total
    separation = abs(pos_power - neg_power)
    if separation < 0.12:
        return 'mixed', separation
    return ('sector_pos' if pos_power >= neg_power else 'sector_neg'), separation


def inject_latent_seed_kick(
    psi: np.ndarray,
    G: nx.Graph,
    nodes: list[int],
    oriented_edges: list[tuple[int, int]],
    triangles: list[tuple[int, int, int]],
    family: dict[str, Any],
    common: dict[str, Any],
    detuning: float,
    amplitude_scale: float,
) -> np.ndarray:
    left_packet = initial_packet(
        G=G,
        nodes=nodes,
        oriented_edges=oriented_edges,
        triangles=triangles,
        seed=int(family['left_seed']),
        kernel_width=float(common['kernel_width']),
        direction_sign=+1.0,
        grade0_scale=float(common['grade0_scale']),
        grade1_scale=float(common['grade1_scale']),
        grade2_scale=float(common['grade2_scale']),
    )
    right_packet = initial_packet(
        G=G,
        nodes=nodes,
        oriented_edges=oriented_edges,
        triangles=triangles,
        seed=int(family['right_seed']),
        kernel_width=float(common['kernel_width']),
        direction_sign=-1.0,
        grade0_scale=float(common['grade0_scale']),
        grade1_scale=float(common['grade1_scale']),
        grade2_scale=float(common['grade2_scale']),
    )
    phase = np.exp(1j * math.pi * float(detuning))
    kicked = psi + float(amplitude_scale) * (left_packet + phase * right_packet)
    norm = max(float(np.linalg.norm(kicked)), 1.0e-12)
    return kicked / norm


def build_protocol_schedule(common: dict[str, Any], protocol_id: str) -> list[dict[str, Any]]:
    detunings = [float(value) for value in common['latent_detunings']]
    protections = [float(value) for value in common['protection_strengths']]
    center = int(common['kernel_radius_center'])
    radii = [max(1, center + int(offset)) for offset in common['radius_offsets']]
    phase_steps = int(common['phase_steps'])
    restore_steps = int(common['restore_steps'])

    forward = [
        {
            'phase': 'forward',
            'detuning': detuning,
            'protection': protection,
            'radius': radius,
            'direction': +1.0,
            'steps': phase_steps,
            'apply_kick': True,
        }
        for detuning, protection, radius in zip(detunings, protections, radii)
    ]
    reverse = [
        {
            'phase': 'reverse',
            'detuning': detuning,
            'protection': protection,
            'radius': radius,
            'direction': -1.0,
            'steps': phase_steps,
            'apply_kick': False,
        }
        for detuning, protection, radius in zip(reversed(detunings), reversed(protections), reversed(radii))
    ]
    restore = [
        {
            'phase': 'restore',
            'detuning': 0.0,
            'protection': 1.0,
            'radius': center,
            'direction': +1.0,
            'steps': restore_steps,
            'apply_kick': False,
        }
    ]
    reperturb = [
        {
            'phase': 'reperturb',
            'detuning': detuning,
            'protection': protection,
            'radius': radius,
            'direction': +1.0,
            'steps': phase_steps,
            'apply_kick': True,
        }
        for detuning, protection, radius in zip(detunings, protections, radii)
    ]

    if protocol_id == 'forward_only':
        return forward
    if protocol_id == 'forward_reverse_restore':
        return forward + reverse + restore
    if protocol_id == 'forward_reverse_reperturb':
        return forward + reverse + reperturb
    raise ValueError(f'unknown protocol_id: {protocol_id}')


def trace_plot(
    path: Path,
    run_id: str,
    topology_codes: list[int],
    mismatch_gradient: list[float],
    flux_projection: list[float],
    grading_trace: list[float],
    return_error_trace: list[int],
) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(16.0, 4.2))
    axes[0].plot(range(len(topology_codes)), topology_codes, marker='o')
    axes[0].set_yticks([0, 1, 2, 3, 4])
    axes[0].set_yticklabels(['braid', 'smear', 'mixed', 'capture', 'none'])
    axes[0].set_title('Topology trace')
    axes[0].grid(alpha=0.25)
    axes[1].plot(range(len(mismatch_gradient)), mismatch_gradient, marker='o', color='tab:orange')
    axes[1].set_title('Mismatch ledger')
    axes[1].grid(alpha=0.25)
    axes[2].plot(range(len(flux_projection)), flux_projection, marker='o', color='tab:red')
    axes[2].set_title('Signed-field flux mismatch')
    axes[2].grid(alpha=0.25)
    axes[3].plot(range(len(grading_trace)), grading_trace, marker='o', label='grading separation')
    axes[3].plot(range(len(return_error_trace)), return_error_trace, marker='s', label='return error')
    axes[3].set_title('Recovered grading / return')
    axes[3].legend(fontsize=8)
    axes[3].grid(alpha=0.25)
    fig.suptitle(run_id)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def summary_plot(path: Path, rows: list[dict[str, Any]], field: str, title: str) -> None:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row['family_label'])].append(row)
    fig, ax = plt.subplots(figsize=(11.5, 4.5))
    for family, family_rows in grouped.items():
        family_rows = sorted(family_rows, key=lambda item: (str(item['protocol_id']), str(item['resolution_id'])))
        x = np.arange(len(family_rows))
        y = [float(item[field]) for item in family_rows]
        ax.plot(x, y, marker='o', label=family)
    ax.set_title(title)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def closure_mode(
    orient_success: bool,
    evolution_success: bool,
    grading_success: bool,
    recovery_success: bool,
    irreversibility_success: bool,
) -> str:
    if orient_success and evolution_success and grading_success and recovery_success and irreversibility_success:
        return 'positive_structural_closure'
    if orient_success and evolution_success and irreversibility_success:
        return 'split_closure'
    return 'negative_closure'


def note_text(json_rel: str, csv_rel: str, rows: list[dict[str, Any]], plots: list[str], summary: dict[str, Any]) -> str:
    lines = [
        '# Stage 12 Intrinsic Grading and Structural Irreversibility v1',
        '',
        f'Timestamped JSON: `{json_rel}`',
        f'Timestamped CSV: `{csv_rel}`',
        '',
        f"Global read: `{summary['closure_mode']}`",
        '',
        '## Obligation diagnostics',
        f"- `O_C0_ORIENT`: positive-family count = `{summary['positive_alignment_families']}`; success = `{summary['orient_success']}`",
        f"- `O_C0_EVOLUTION`: median refinement ratio = `{summary['median_refinement_ratio']:.4f}`; success = `{summary['evolution_success']}`",
        f"- `Intrinsic grading`: entropy = `{summary['mean_recovery_entropy']:.4f}`, separation = `{summary['mean_sector_separation']:.4f}`, monotonicity = `{summary['mean_detuning_monotonicity']:.4f}`; success = `{summary['grading_success']}`",
        f"- `Family-wide recovery`: ordered settings fraction = `{summary['ordered_fraction']:.4f}`; success = `{summary['recovery_success']}`",
        f"- `Structural irreversibility`: bidirectional positive fraction = `{summary['irreversible_fraction']:.4f}`; success = `{summary['irreversibility_success']}`",
        '',
        '## Family medians',
    ]
    for family_id, values in summary['family_medians'].items():
        lines.append(
            f"- `{family_id}`: alignment=`{values['mismatch_ledger_alignment']:.4f}`, "
            f"protected_fraction=`{values['protected_fraction']:.4f}`, "
            f"irreversibility=`{values['structural_irreversibility_score']:.4f}`, "
            f"entropy=`{values['recovery_entropy']:.4f}`"
        )
    lines.extend(['', '## Per-run summary'])
    for row in rows:
        lines.append(
            f"- `{row['run_id']}`: family=`{row['family_label']}`, protocol=`{row['protocol_label']}`, "
            f"resolution=`{row['resolution_label']}`, forward=`{row['dominant_forward_topology_class']}`, "
            f"return=`{row['return_checkpoint_topology_class']}`, final=`{row['final_topology_class']}`, "
            f"protected=`{float(row['protected_fraction']):.4f}`, irreversibility=`{float(row['structural_irreversibility_score']):.4f}`, "
            f"entropy=`{float(row['recovery_entropy']):.4f}`"
        )
    lines.extend(['', '## Plots'])
    for rel in plots:
        lines.append(f'- `{rel}`')
    return '\n'.join(lines) + '\n'


def main() -> None:
    args = parse_args()
    runsheet = load_runsheet(args.runsheet)
    common = runsheet['common_fields']

    selected: list[tuple[str, dict[str, Any], dict[str, Any], dict[str, Any]]] = []
    for family in runsheet['families']:
        for protocol in runsheet['protocols']:
            for resolution in runsheet['resolutions']:
                run_id = f"S12_{family['family_id']}_{protocol['protocol_id']}_{resolution['resolution_id']}"
                selected.append((run_id, family, protocol, resolution))
    if args.run_ids:
        wanted = set(args.run_ids)
        selected = [item for item in selected if item[0] in wanted]

    rows: list[dict[str, Any]] = []
    plot_paths: list[Path] = []
    topology_map = {
        'braid_like_exchange': 0,
        'transfer_smeared': 1,
        'unresolved_mixed': 2,
        'localized_capture': 3,
        'no_protected_structure': 4,
    }

    for run_id, family, protocol, resolution in selected:
        G = matched_family_graph(str(family['family_id']))
        nodes = sorted(int(node) for node in G.nodes())
        triangles = triangle_list(G)
        tri_counts = node_triangle_counts(G, triangles)
        motif_nodes = motif_nodes_for_family(G, triangles, int(family['left_seed']), int(family['right_seed']))
        path_nodes = nx.shortest_path(G, source=int(family['left_seed']), target=int(family['right_seed']))
        path_mask = np.asarray([1.0 if node in path_nodes else 0.0 for node in nodes], dtype=float)
        motif_mask = np.asarray([1.0 if node in motif_nodes else 0.0 for node in nodes], dtype=float)
        substeps = int(resolution['substeps'])
        schedule = build_protocol_schedule(common, str(protocol['protocol_id']))

        mismatch = np.zeros(len(nodes), dtype=float)
        initial_edges = orient_edges_from_mismatch(G, mismatch, nodes, tri_counts)
        psi = inject_latent_seed_kick(
            psi=np.zeros(len(nodes) + len(initial_edges) + len(triangles), dtype=complex),
            G=G,
            nodes=nodes,
            oriented_edges=initial_edges,
            triangles=triangles,
            family=family,
            common=common,
            detuning=float(common['latent_detunings'][0]),
            amplitude_scale=1.0,
        )
        states: list[np.ndarray] = [psi.copy()]
        topology_trace: list[str] = []
        topology_codes: list[int] = []
        mismatch_gradient_trace: list[float] = []
        protected_index_trace: list[float] = []
        flux_projection_trace: list[float] = []
        grading_band_trace: list[str] = []
        grading_separation_trace: list[float] = []
        return_error_trace: list[int] = []
        segment_separation_means: list[float] = []
        segment_detunings: list[float] = []
        segment_forward_topologies: list[str] = []
        phase_lookup: list[str] = []
        flow_trace: list[float] = []
        reverse_flow_trace: list[float] = []
        forward_flow_trace: list[float] = []

        for segment in schedule:
            oriented_edges = orient_edges_from_mismatch(G, mismatch, nodes, tri_counts)
            if bool(segment['apply_kick']):
                psi = inject_latent_seed_kick(
                    psi=psi,
                    G=G,
                    nodes=nodes,
                    oriented_edges=oriented_edges,
                    triangles=triangles,
                    family=family,
                    common=common,
                    detuning=float(segment['detuning']),
                    amplitude_scale=float(common['segment_kick']),
                )
                states.append(psi.copy())

            segment_separations: list[float] = []
            segment_topologies: list[str] = []
            for _ in range(int(segment['steps'])):
                for _substep in range(substeps):
                    direction = float(segment['direction'])
                    radius = int(segment['radius'])
                    node_activity = node_activity_from_state(
                        psi=psi,
                        nodes=nodes,
                        oriented_edges=oriented_edges,
                        triangles=triangles,
                    )
                    node_kernel = node_kernel_matrix(G, radius=radius)
                    forcing = node_kernel @ node_activity
                    forcing += float(segment['protection']) * (0.65 * path_mask + 0.35 * motif_mask)
                    delta_mismatch = direction * (float(common['epsilon_mismatch']) / substeps) * forcing
                    mismatch = mismatch + delta_mismatch

                    oriented_edges = orient_edges_from_mismatch(G, mismatch, nodes, tri_counts)
                    d0, edge_lookup, _node_to_idx = build_d0(nodes, oriented_edges)
                    d1 = build_d1(triangles, edge_lookup)
                    D_signed, _Delta = block_dirac(d0, d1)
                    supports = basis_supports(nodes, oriented_edges, triangles)
                    coo = D_signed.tocoo()
                    support_distances = support_distance_vector(G, supports, coo.row, coo.col)
                    operator = radius_operator(D_signed, support_distances, radius)

                    psi_next = field_update(operator, psi, direction * (float(common['epsilon_field']) / substeps))
                    n0 = len(nodes)
                    n1 = len(oriented_edges)
                    mismatch_projection = np.asarray(d0 @ delta_mismatch, dtype=float)
                    flux = np.abs(
                        (psi_next[n0:n0 + n1] - psi[n0:n0 + n1])
                        / max(float(common['epsilon_field']) / substeps, 1.0e-12)
                    )
                    flux_projection_trace.append(float(np.linalg.norm(np.abs(mismatch_projection) - flux)))

                    states.append(psi_next.copy())
                    window = states[max(0, len(states) - 12):]
                    metrics = metrics_for_state_label_free(
                        states_window=window,
                        nodes=nodes,
                        oriented_edges=oriented_edges,
                        triangles=triangles,
                        d1=d1,
                        block_sizes=(len(nodes), len(oriented_edges), len(triangles)),
                    )
                    topology = str(metrics['topology_class'])
                    topology_trace.append(topology)
                    topology_codes.append(topology_map[topology])
                    current_protected = 1.0 if topology in {'braid_like_exchange', 'transfer_smeared'} else 0.0
                    protected_index_trace.append(
                        (protected_index_trace[-1] if protected_index_trace else 0.0) + current_protected
                    )
                    edge_grad = np.asarray(
                        [abs(mismatch[nodes.index(u)] - mismatch[nodes.index(v)]) for u, v in oriented_edges],
                        dtype=float,
                    )
                    mismatch_gradient_trace.append(float(np.mean(edge_grad)) if edge_grad.size else 0.0)
                    band, separation = recovered_grading_band(window, nodes, oriented_edges, triangles)
                    grading_band_trace.append(band)
                    grading_separation_trace.append(float(separation))
                    segment_separations.append(float(separation))
                    segment_topologies.append(topology)
                    phase_lookup.append(str(segment['phase']))

                    flow_value = float(metrics['flow_concentration_index'])
                    flow_trace.append(flow_value)
                    if str(segment['phase']) == 'forward':
                        forward_flow_trace.append(flow_value)
                    elif str(segment['phase']) == 'reverse':
                        reverse_flow_trace.append(flow_value)

                    psi = psi_next

            segment_separation_means.append(float(np.mean(segment_separations)) if segment_separations else 0.0)
            segment_detunings.append(float(segment['detuning']))
            if str(segment['phase']) == 'forward':
                segment_forward_topologies.extend(segment_topologies)

        dominant_forward_topology = mode_or_default(segment_forward_topologies, 'no_protected_structure')
        return_indices = [idx for idx, phase in enumerate(phase_lookup) if phase in {'reverse', 'restore'}]
        return_topologies = [topology_trace[idx] for idx in return_indices]
        return_checkpoint_topology = mode_or_default(return_topologies, topology_trace[-1] if topology_trace else 'no_protected_structure')
        final_topology = topology_trace[-1] if topology_trace else 'no_protected_structure'
        topology_return_error = 0 if return_checkpoint_topology == dominant_forward_topology else 1
        structural_irreversibility = (
            float(topology_return_error)
            + area_between_curves(forward_flow_trace, reverse_flow_trace)
        )
        protected_fraction = float(
            sum(1 for topo in topology_trace if topo in {'braid_like_exchange', 'transfer_smeared'})
            / max(len(topology_trace), 1)
        )
        detuning_monotonicity = correlation_or_zero(segment_detunings, segment_separation_means)
        row = {
            'run_id': run_id,
            'family_id': str(family['family_id']),
            'family_label': str(family['family_label']),
            'protocol_id': str(protocol['protocol_id']),
            'protocol_label': str(protocol['protocol_label']),
            'resolution_id': str(resolution['resolution_id']),
            'resolution_label': str(resolution['resolution_label']),
            'dominant_forward_topology_class': dominant_forward_topology,
            'return_checkpoint_topology_class': return_checkpoint_topology,
            'final_topology_class': final_topology,
            'survival_time': float(
                sum(1 for topo in topology_trace if topo in {'braid_like_exchange', 'transfer_smeared'})
            ),
            'protected_fraction': protected_fraction,
            'topology_return_error': int(topology_return_error),
            'mismatch_ledger_alignment': correlation_or_zero(protected_index_trace, mismatch_gradient_trace),
            'signed_field_flux_projection': float(np.mean(flux_projection_trace)) if flux_projection_trace else 0.0,
            'family_ordering_score': 0.0,
            'structural_irreversibility_score': float(structural_irreversibility),
            'recovery_entropy': conditional_entropy(topology_trace, grading_band_trace),
            'sector_separation': float(np.mean(grading_separation_trace)) if grading_separation_trace else 0.0,
            'detuning_monotonicity': detuning_monotonicity,
            'refinement_drift': 0.0,
            'notes': '',
        }
        trace_path = WORK_PLOT_DIR / f'{run_id}_trace.png'
        trace_plot(
            trace_path,
            run_id,
            topology_codes,
            mismatch_gradient_trace,
            flux_projection_trace,
            grading_separation_trace,
            [topology_return_error] * len(topology_codes),
        )
        row['trace_plot'] = str(trace_path)
        rows.append(row)
        plot_paths.append(trace_path)

    row_lookup = {
        (str(row['family_id']), str(row['protocol_id']), str(row['resolution_id'])): row
        for row in rows
    }
    for family in runsheet['families']:
        for protocol in runsheet['protocols']:
            base_key = (str(family['family_id']), str(protocol['protocol_id']), 'base')
            refined_key = (str(family['family_id']), str(protocol['protocol_id']), 'refined')
            if base_key in row_lookup and refined_key in row_lookup:
                base_row = row_lookup[base_key]
                refined_row = row_lookup[refined_key]
                drift = abs(float(base_row['protected_fraction']) - float(refined_row['protected_fraction']))
                drift += abs(
                    float(base_row['structural_irreversibility_score'])
                    - float(refined_row['structural_irreversibility_score'])
                )
                base_row['refinement_drift'] = drift
                refined_row['refinement_drift'] = drift

    grouped_family: dict[str, list[dict[str, Any]]] = defaultdict(list)
    grouped_setting: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    grouped_pair: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped_family[str(row['family_id'])].append(row)
        grouped_setting[(str(row['protocol_id']), str(row['resolution_id']))].append(row)
        grouped_pair[(str(row['family_id']), str(row['protocol_id']))].append(row)

    family_medians: dict[str, dict[str, float]] = {}
    for family_id, family_rows in grouped_family.items():
        family_medians[family_id] = {
            'mismatch_ledger_alignment': float(np.median([float(item['mismatch_ledger_alignment']) for item in family_rows])),
            'protected_fraction': float(np.median([float(item['protected_fraction']) for item in family_rows])),
            'structural_irreversibility_score': float(
                np.median([float(item['structural_irreversibility_score']) for item in family_rows])
            ),
            'recovery_entropy': float(np.median([float(item['recovery_entropy']) for item in family_rows])),
        }

    positive_alignment_families = sum(
        1 for values in family_medians.values() if values['mismatch_ledger_alignment'] > 0.15
    )
    orient_success = positive_alignment_families >= 2

    refinement_ratios: list[float] = []
    for key, pair_rows in grouped_pair.items():
        lookup = {str(item['resolution_id']): item for item in pair_rows}
        if 'base' in lookup and 'refined' in lookup:
            base_flux = max(float(lookup['base']['signed_field_flux_projection']), 1.0e-12)
            refined_flux = float(lookup['refined']['signed_field_flux_projection'])
            refinement_ratios.append(refined_flux / base_flux)
    median_refinement_ratio = float(np.median(refinement_ratios)) if refinement_ratios else 0.0
    evolution_success = 0.5 <= median_refinement_ratio <= 1.5

    mean_recovery_entropy = float(np.mean([float(row['recovery_entropy']) for row in rows])) if rows else 0.0
    mean_sector_separation = float(np.mean([float(row['sector_separation']) for row in rows])) if rows else 0.0
    mean_detuning_monotonicity = float(np.mean([float(row['detuning_monotonicity']) for row in rows])) if rows else 0.0
    mean_refinement_drift = float(np.mean([float(row['refinement_drift']) for row in rows])) if rows else 0.0
    grading_success = (
        mean_recovery_entropy <= 0.32
        and mean_sector_separation >= 0.18
        and mean_detuning_monotonicity >= 0.15
        and mean_refinement_drift <= 0.35
    )

    ordered_settings = 0
    for setting, setting_rows in grouped_setting.items():
        lookup = {str(item['family_id']): float(item['protected_fraction']) for item in setting_rows}
        ordered = (
            lookup.get('clustered_braid_seed', -1.0) >= lookup.get('corridor_transport_seed', -1.0)
            and lookup.get('corridor_transport_seed', -1.0) >= lookup.get('triad_competition_seed', -1.0)
        )
        if ordered:
            ordered_settings += 1
        for item in setting_rows:
            item['family_ordering_score'] = 1.0 if ordered else 0.0
    ordered_fraction = float(ordered_settings / max(len(grouped_setting), 1))
    recovery_success = ordered_fraction >= 0.70

    bidirectional_rows = [
        row for row in rows
        if str(row['protocol_id']) in {'forward_reverse_restore', 'forward_reverse_reperturb'}
    ]
    irreversible_fraction = float(
        sum(1 for row in bidirectional_rows if float(row['structural_irreversibility_score']) >= 0.18)
        / max(len(bidirectional_rows), 1)
    )
    irreversibility_success = irreversible_fraction >= 0.30

    mode = closure_mode(
        orient_success=orient_success,
        evolution_success=evolution_success,
        grading_success=grading_success,
        recovery_success=recovery_success,
        irreversibility_success=irreversibility_success,
    )
    for row in rows:
        row['notes'] = mode

    entropy_plot = WORK_PLOT_DIR / 'stage_12_entropy_panel.png'
    irrev_plot = WORK_PLOT_DIR / 'stage_12_irreversibility_panel.png'
    drift_plot = WORK_PLOT_DIR / 'stage_12_refinement_panel.png'
    summary_plot(entropy_plot, rows, 'recovery_entropy', 'Stage 12 recovered grading entropy')
    summary_plot(irrev_plot, rows, 'structural_irreversibility_score', 'Stage 12 structural irreversibility')
    summary_plot(drift_plot, rows, 'refinement_drift', 'Stage 12 refinement drift')
    plot_paths.extend([entropy_plot, irrev_plot, drift_plot])

    summary = {
        'closure_mode': mode,
        'positive_alignment_families': positive_alignment_families,
        'orient_success': orient_success,
        'median_refinement_ratio': median_refinement_ratio,
        'evolution_success': evolution_success,
        'mean_recovery_entropy': mean_recovery_entropy,
        'mean_sector_separation': mean_sector_separation,
        'mean_detuning_monotonicity': mean_detuning_monotonicity,
        'mean_refinement_drift': mean_refinement_drift,
        'grading_success': grading_success,
        'ordered_fraction': ordered_fraction,
        'recovery_success': recovery_success,
        'irreversible_fraction': irreversible_fraction,
        'irreversibility_success': irreversibility_success,
        'family_medians': family_medians,
    }

    result_payload = {
        'stage': 'stage_12_intrinsic_grading_and_structural_irreversibility_probe',
        'runsheet': str(Path(args.runsheet).relative_to(REPO_ROOT)),
        'rows': rows,
        'summary': summary,
    }
    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug='stage_12_intrinsic_grading_and_structural_irreversibility_probe',
        result=result_payload,
        csv_rows=rows,
        csv_fieldnames=CSV_FIELDS,
        plot_paths=plot_paths,
    )
    NOTE_PATH.write_text(
        note_text(
            str(json_path.relative_to(REPO_ROOT)),
            str(csv_path.relative_to(REPO_ROOT)),
            rows,
            stamped_plots,
            summary,
        ),
        encoding='utf-8',
    )


if __name__ == '__main__':
    main()
