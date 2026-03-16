#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, deque
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt, save_atlas_payload
from stage23_1_dk_collision_geometry_scan import mean_pair_separation_series
from stage23_10_dk_braid_mechanism_isolation_probe import analyse_states, packet_state_with_profile
from stage23_9_dk_collision_phenomenology_consolidation import evolve_signed
from stage9b_common import build_dk2d_complex, first_order_dt, pack_positions

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage23_12_dk_effective_phase_diagram_closure_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_III_C3_Effective_Phase_Diagram_Closure_v1.md'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage23_12_effective_phase_diagram_closure'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'run_id',
    'run_order',
    'phase_label',
    'phase_offset_fraction_of_pi',
    'width_regime',
    'width_ratio_a_to_b',
    'topology_class',
    'dominant_reverse_topology',
    'topology_survival_time',
    'refinement_stability_flag',
    'weak_coupling_stability_flag',
    'bidirectional_stability_flag',
    'flow_concentration_index',
    'grade_exchange_coherence',
    'separation_oscillation_indicator',
    'derived_phase_label',
    'low_beta_base_topology',
    'low_beta_refined_topology',
    'high_beta_base_topology',
    'high_beta_refined_topology',
    'notes',
]

PHASE_CODE = {
    'stable_braid_phase': 0,
    'smeared_transfer_phase': 1,
    'localized_encounter_phase': 2,
    'transient_mixed_phase': 3,
}

PHASE_COLORBAR_LABELS = ['stable braid', 'smeared', 'localized', 'transient']


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run Stage 23.12 DK effective phase diagram closure.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    return parser.parse_args()


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def sigma_pair(base_sigma: float, ratio_a_to_b: float) -> tuple[float, float]:
    root = math.sqrt(max(ratio_a_to_b, 1.0e-12))
    return base_sigma * root, base_sigma / root


def dominant_label(labels: list[str], anchor_label: str) -> str:
    counts = Counter(labels)
    max_count = max(counts.values())
    winners = [label for label, count in counts.items() if count == max_count]
    if anchor_label in winners:
        return anchor_label
    priority = ['braid_like_exchange', 'transfer_smeared', 'localized_capture', 'dispersive_pass']
    return min(winners, key=lambda label: priority.index(label) if label in priority else len(priority))


def topology_survival_time(metrics: dict[str, Any], dt: float) -> float:
    topology = str(metrics['topology_class'])
    trace = [str(label) for label in metrics['topology_trace']]
    return float(sum(label == topology for label in trace) * dt)


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


def simulate_variant(
    phase_fraction: float,
    width_ratio: float,
    resolution: int,
    beta: float,
    common: dict[str, Any],
) -> dict[str, Any]:
    epsilon = float(common['base_epsilon'])
    sigma_base = float(common['mean_width'])
    amplitude = 0.5 * float(common['amplitude_scale'])
    t_final = float(common['t_final'])
    separation = float(common['separation'])

    complex_data = build_dk2d_complex(n_side=int(resolution), epsilon=epsilon)
    positions = pack_positions(complex_data.points, complex_data.edge_midpoints, complex_data.face_centers)
    block_sizes = complex_data.block_sizes

    sigma_a, sigma_b = sigma_pair(sigma_base, width_ratio)
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
            phase_offset=float(phase_fraction) * math.pi,
            kick_vector=kicks[1],
            kick_cycles=1.0,
            grade0_scale=1.0,
            grade1_scale=1.0,
            anisotropy=(1.0, 1.0),
        ),
    ]

    psi0 = np.sum(packet_states, axis=0)
    dt = first_order_dt(complex_data.dirac_kahler, float(common['dt_scale']))
    steps = max(2, int(math.ceil(t_final / dt)))
    close_threshold = 2.0 * float(np.mean([sigma_a, sigma_b]))

    forward_states = evolve_signed(complex_data.dirac_kahler, block_sizes, psi0, dt, steps, float(beta), direction=1.0)
    reverse_states = evolve_signed(
        complex_data.dirac_kahler,
        block_sizes,
        forward_states[-1],
        dt,
        steps,
        float(beta),
        direction=-1.0,
    )

    forward_metrics = analyse_states(
        states=forward_states,
        packet_states=packet_states,
        positions=positions,
        edge_midpoints=complex_data.edge_midpoints,
        block_sizes=block_sizes,
        resolution=int(resolution),
        anchors=centers,
        close_threshold=close_threshold,
    )
    reverse_metrics = analyse_states(
        states=reverse_states,
        packet_states=packet_states,
        positions=positions,
        edge_midpoints=complex_data.edge_midpoints,
        block_sizes=block_sizes,
        resolution=int(resolution),
        anchors=centers,
        close_threshold=close_threshold,
    )

    return {
        'resolution': int(resolution),
        'beta': float(beta),
        'dt': float(dt),
        'steps': int(steps),
        'forward': forward_metrics,
        'reverse': reverse_metrics,
        'forward_survival_time': topology_survival_time(forward_metrics, dt),
        'reverse_survival_time': topology_survival_time(reverse_metrics, dt),
        'separation_oscillation_indicator': separation_oscillation_indicator(forward_metrics['center_histories']),
    }


def derive_phase_label(
    topology_class: str,
    refinement_flag: bool,
    weak_coupling_flag: bool,
    bidirectional_flag: bool,
) -> str:
    stable = bool(refinement_flag and weak_coupling_flag and bidirectional_flag)
    if topology_class == 'braid_like_exchange' and stable:
        return 'stable_braid_phase'
    if topology_class == 'transfer_smeared' and stable:
        return 'smeared_transfer_phase'
    if topology_class == 'localized_capture' and stable:
        return 'localized_encounter_phase'
    return 'transient_mixed_phase'


def build_phase_grid(rows: list[dict[str, Any]], phase_values: list[float], width_values: list[float]) -> np.ndarray:
    grid = np.zeros((len(width_values), len(phase_values)), dtype=float)
    for row in rows:
        i = width_values.index(float(row['width_ratio_a_to_b']))
        j = phase_values.index(float(row['phase_offset_fraction_of_pi']))
        grid[i, j] = PHASE_CODE[str(row['derived_phase_label'])]
    return grid


def build_survival_grid(rows: list[dict[str, Any]], phase_values: list[float], width_values: list[float]) -> np.ndarray:
    grid = np.zeros((len(width_values), len(phase_values)), dtype=float)
    for row in rows:
        i = width_values.index(float(row['width_ratio_a_to_b']))
        j = phase_values.index(float(row['phase_offset_fraction_of_pi']))
        grid[i, j] = float(row['topology_survival_time'])
    return grid


def stable_components(rows: list[dict[str, Any]], phase_values: list[float], width_values: list[float]) -> list[dict[str, Any]]:
    cell_map = {
        (width_values.index(float(row['width_ratio_a_to_b'])), phase_values.index(float(row['phase_offset_fraction_of_pi']))): row
        for row in rows
    }
    visited: set[tuple[int, int]] = set()
    components: list[dict[str, Any]] = []

    for key, row in cell_map.items():
        if key in visited or str(row['derived_phase_label']) == 'transient_mixed_phase':
            continue
        label = str(row['derived_phase_label'])
        queue: deque[tuple[int, int]] = deque([key])
        visited.add(key)
        cells: list[tuple[int, int]] = []
        while queue:
            ci, cj = queue.popleft()
            cells.append((ci, cj))
            for ni, nj in ((ci - 1, cj), (ci + 1, cj), (ci, cj - 1), (ci, cj + 1)):
                neighbor = (ni, nj)
                if neighbor in visited or neighbor not in cell_map:
                    continue
                neighbor_row = cell_map[neighbor]
                if str(neighbor_row['derived_phase_label']) != label:
                    continue
                visited.add(neighbor)
                queue.append(neighbor)
        width_span = [float(width_values[idx]) for idx, _ in cells]
        phase_span = [float(phase_values[idx]) for _, idx in cells]
        components.append(
            {
                'label': label,
                'size': len(cells),
                'cells': [(int(i), int(j)) for i, j in cells],
                'width_min': min(width_span),
                'width_max': max(width_span),
                'phase_min': min(phase_span),
                'phase_max': max(phase_span),
            }
        )
    return components


def count_boundary_segments(rows: list[dict[str, Any]], phase_values: list[float], width_values: list[float]) -> int:
    label_grid = [
        [
            next(
                str(row['derived_phase_label'])
                for row in rows
                if float(row['width_ratio_a_to_b']) == width
                and float(row['phase_offset_fraction_of_pi']) == phase
            )
            for phase in phase_values
        ]
        for width in width_values
    ]
    segments = 0
    for i in range(len(width_values)):
        for j in range(len(phase_values) - 1):
            if label_grid[i][j] != label_grid[i][j + 1]:
                segments += 1
    for i in range(len(width_values) - 1):
        for j in range(len(phase_values)):
            if label_grid[i][j] != label_grid[i + 1][j]:
                segments += 1
    return segments


def minimal_boundary_summary(components: list[dict[str, Any]], rows: list[dict[str, Any]]) -> list[str]:
    lines: list[str] = []
    stable_rows = [row for row in rows if str(row['derived_phase_label']) != 'transient_mixed_phase']
    if not stable_rows:
        return ['No stable non-transient region closed on the tested manifold.']
    for component in sorted(components, key=lambda item: (-int(item['size']), str(item['label']))):
        label = str(component['label'])
        lines.append(
            f"{label}: size={int(component['size'])}, phase in [{component['phase_min']:.3f}, {component['phase_max']:.3f}], "
            f"width ratio in [{component['width_min']:.2f}, {component['width_max']:.2f}]"
        )
    return lines


def plot_phase_map(path: Path, rows: list[dict[str, Any]], phase_values: list[float], width_labels: list[str], width_values: list[float]) -> None:
    grid = build_phase_grid(rows, phase_values, width_values)
    fig, ax = plt.subplots(figsize=(9.0, 4.6))
    im = ax.imshow(grid, aspect='auto', cmap='viridis', vmin=0.0, vmax=3.0)
    ax.set_xticks(range(len(phase_values)), [f'{phase:.3f}' for phase in phase_values])
    ax.set_yticks(range(len(width_labels)), width_labels)
    ax.set_xlabel('phase offset (fraction of pi)')
    ax.set_ylabel('width regime')
    ax.set_title('Stage III-C3 parameter-space phase map')
    for i, width in enumerate(width_values):
        for j, phase in enumerate(phase_values):
            row = next(
                item for item in rows
                if float(item['width_ratio_a_to_b']) == width
                and float(item['phase_offset_fraction_of_pi']) == phase
            )
            ax.text(j, i, str(row['topology_class']).replace('_', '\n'), ha='center', va='center', fontsize=7, color='white')
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cbar.set_ticks([0, 1, 2, 3])
    cbar.set_ticklabels(PHASE_COLORBAR_LABELS)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_survival_heatmap(path: Path, rows: list[dict[str, Any]], phase_values: list[float], width_labels: list[str], width_values: list[float]) -> None:
    grid = build_survival_grid(rows, phase_values, width_values)
    fig, ax = plt.subplots(figsize=(9.0, 4.6))
    im = ax.imshow(grid, aspect='auto', cmap='magma')
    ax.set_xticks(range(len(phase_values)), [f'{phase:.3f}' for phase in phase_values])
    ax.set_yticks(range(len(width_labels)), width_labels)
    ax.set_xlabel('phase offset (fraction of pi)')
    ax.set_ylabel('width regime')
    ax.set_title('Stage III-C3 topology survival heatmap')
    for i in range(len(width_values)):
        for j in range(len(phase_values)):
            ax.text(j, i, f'{grid[i, j]:.2f}', ha='center', va='center', fontsize=8, color='white')
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def plot_stability_table(path: Path, rows: list[dict[str, Any]], phase_values: list[float], width_labels: list[str], width_values: list[float]) -> None:
    score_grid = np.zeros((len(width_values), len(phase_values)), dtype=float)
    fig, ax = plt.subplots(figsize=(9.6, 4.8))
    for i, width in enumerate(width_values):
        for j, phase in enumerate(phase_values):
            row = next(
                item for item in rows
                if float(item['width_ratio_a_to_b']) == width
                and float(item['phase_offset_fraction_of_pi']) == phase
            )
            score = (
                int(row['refinement_stability_flag'])
                + int(row['weak_coupling_stability_flag'])
                + int(row['bidirectional_stability_flag'])
            )
            score_grid[i, j] = float(score)
    im = ax.imshow(score_grid, aspect='auto', cmap='cividis', vmin=0.0, vmax=3.0)
    ax.set_xticks(range(len(phase_values)), [f'{phase:.3f}' for phase in phase_values])
    ax.set_yticks(range(len(width_labels)), width_labels)
    ax.set_xlabel('phase offset (fraction of pi)')
    ax.set_ylabel('width regime')
    ax.set_title('Stage III-C3 refinement / beta / reverse stability table')
    for i, width in enumerate(width_values):
        for j, phase in enumerate(phase_values):
            row = next(
                item for item in rows
                if float(item['width_ratio_a_to_b']) == width
                and float(item['phase_offset_fraction_of_pi']) == phase
            )
            ax.text(
                j,
                i,
                f"R{int(row['refinement_stability_flag'])}\nB{int(row['weak_coupling_stability_flag'])}\nD{int(row['bidirectional_stability_flag'])}",
                ha='center',
                va='center',
                fontsize=7,
                color='white',
            )
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cbar.set_ticks([0, 1, 2, 3])
    cbar.set_ticklabels(['0', '1', '2', '3'])
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(
    json_path: Path,
    csv_path: Path,
    rows: list[dict[str, Any]],
    stamped_plots: list[str],
    summary: dict[str, Any],
    boundary_lines: list[str],
) -> None:
    lines = [
        '# Stage III-C3 Effective Phase Diagram Closure v1',
        '',
        'This closure pass stays inside the frozen HAOS (Harmonic Address Operating System) / IIP (Interaction Invariance Physics) clustered Dirac-Kahler collision sector.',
        '',
        f"Timestamped JSON: `{json_path.relative_to(REPO_ROOT)}`",
        f"Timestamped CSV: `{csv_path.relative_to(REPO_ROOT)}`",
        '',
        f"phase_diagram_closed = `{'TRUE' if bool(summary['phase_diagram_closed']) else 'FALSE'}`",
        f"closure_classification = `{summary['closure_classification']}`",
        '',
        '## Closure note',
        'The source brief specified a three-axis manifold but a 15-run primary grid. This implementation therefore uses a `5 x 3` phase-width lattice and treats `beta in {0.01, 0.02}` as an internal weak-coupling stability probe inside each primary cell rather than as an external 45-point expansion.',
        '',
        '## Summary metrics',
        f"- stable connected region present: `{summary['stable_connected_region_present']}`",
        f"- finite boundary flag: `{summary['finite_boundary_flag']}`",
        f"- bidirectional preservation flag: `{summary['bidirectional_preservation_flag']}`",
        f"- boundary segment count: `{summary['boundary_segment_count']}`",
        f"- stable phase counts: `{summary['stable_phase_counts']}`",
        f"- dominant topology counts: `{summary['dominant_topology_counts']}`",
        '',
        '## Per-cell summary',
    ]
    for row in sorted(rows, key=lambda item: int(item['run_order'])):
        lines.append(
            f"- `{row['run_id']}`: phase=`{row['phase_offset_fraction_of_pi']:.3f}`, width=`{row['width_ratio_a_to_b']:.2f}`, "
            f"topology=`{row['topology_class']}`, phase_label=`{row['derived_phase_label']}`, "
            f"R=`{int(row['refinement_stability_flag'])}`, B=`{int(row['weak_coupling_stability_flag'])}`, D=`{int(row['bidirectional_stability_flag'])}`"
        )

    lines.extend([
        '',
        '## Minimal regime boundary summary',
    ])
    for line in boundary_lines:
        lines.append(f'- {line}')

    lines.extend([
        '',
        '## Interpretation',
    ])
    if bool(summary['phase_diagram_closed']):
        lines.append(
            'The clustered DK sector admits a closed effective phase description on this manifold, but the closure is one-sided: the stable connected region is a smeared-transfer phase, while the nominal braid corridor collapses into transient mixed cells once refinement, weak-coupling consistency, and reverse evolution are enforced together. Phase III therefore closes this branch as a sector-local phenomenology map rather than a general intrinsic braid mechanism.'
        )
    else:
        lines.append(
            'The tested manifold does not close into a stable effective phase diagram under the joint refinement and reverse-evolution checks. The DK collision sector therefore remains descriptive at this level, and a later phase should pivot to a new operator or constraint sector rather than continue local closure tuning.'
        )

    if int(summary['stable_phase_counts'].get('localized_encounter_phase', 0)) == 0:
        lines.append('No localized encounter phase closed on this lattice.')
    if int(summary['stable_phase_counts'].get('stable_braid_phase', 0)) == 0:
        lines.append('No stable braid phase survived the full closure criteria on this lattice.')

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
    phase_values = [float(value) for value in common['phase_values_fraction_of_pi']]
    width_specs = list(common['width_regimes'])
    width_values = [float(item['width_ratio_a_to_b']) for item in width_specs]
    width_labels = [str(item['width_regime']).replace('_', ' ') for item in width_specs]
    resolutions = [int(common['base_resolution']), int(common['refined_resolution'])]
    beta_values = [float(value) for value in common['weak_beta_values']]
    low_beta = min(beta_values)
    high_beta = max(beta_values)
    base_resolution = int(common['base_resolution'])
    refined_resolution = int(common['refined_resolution'])

    rows: list[dict[str, Any]] = []
    details: list[dict[str, Any]] = []

    for run in sorted(runsheet['runs'], key=lambda item: int(item['run_order'])):
        variant_lookup: dict[tuple[float, int], dict[str, Any]] = {}
        for beta in beta_values:
            for resolution in resolutions:
                variant_lookup[(float(beta), int(resolution))] = simulate_variant(
                    phase_fraction=float(run['phase_offset_fraction_of_pi']),
                    width_ratio=float(run['width_ratio_a_to_b']),
                    resolution=int(resolution),
                    beta=float(beta),
                    common=common,
                )

        forward_lookup = {
            key: str(value['forward']['topology_class'])
            for key, value in variant_lookup.items()
        }
        reverse_lookup = {
            key: str(value['reverse']['topology_class'])
            for key, value in variant_lookup.items()
        }
        anchor_key = (high_beta, base_resolution)
        dominant_forward = dominant_label(list(forward_lookup.values()), forward_lookup[anchor_key])
        dominant_reverse = dominant_label(list(reverse_lookup.values()), reverse_lookup[anchor_key])

        refinement_stability_flag = int(
            all(
                forward_lookup[(beta, base_resolution)] == forward_lookup[(beta, refined_resolution)]
                for beta in beta_values
            )
        )
        weak_coupling_stability_flag = int(
            all(
                forward_lookup[(low_beta, resolution)] == forward_lookup[(high_beta, resolution)]
                for resolution in resolutions
            )
        )
        bidirectional_stability_flag = int(
            all(forward_lookup[key] == reverse_lookup[key] for key in forward_lookup)
        )

        dominant_variants = [value for value in variant_lookup.values() if str(value['forward']['topology_class']) == dominant_forward]
        topology_survival = float(np.mean([float(value['forward_survival_time']) for value in dominant_variants])) if dominant_variants else float(variant_lookup[anchor_key]['forward_survival_time'])
        flow_concentration = float(np.mean([float(value['forward']['flow_concentration_index']) for value in variant_lookup.values()]))
        coherence = float(np.mean([float(value['forward']['grade_exchange_coherence']) for value in variant_lookup.values()]))
        separation_oscillation = float(np.mean([float(value['separation_oscillation_indicator']) for value in variant_lookup.values()]))
        derived_label = derive_phase_label(
            topology_class=dominant_forward,
            refinement_flag=bool(refinement_stability_flag),
            weak_coupling_flag=bool(weak_coupling_stability_flag),
            bidirectional_flag=bool(bidirectional_stability_flag),
        )

        row = {
            'run_id': str(run['run_id']),
            'run_order': int(run['run_order']),
            'phase_label': str(run['phase_label']),
            'phase_offset_fraction_of_pi': float(run['phase_offset_fraction_of_pi']),
            'width_regime': str(run['width_regime']),
            'width_ratio_a_to_b': float(run['width_ratio_a_to_b']),
            'topology_class': dominant_forward,
            'dominant_reverse_topology': dominant_reverse,
            'topology_survival_time': topology_survival,
            'refinement_stability_flag': refinement_stability_flag,
            'weak_coupling_stability_flag': weak_coupling_stability_flag,
            'bidirectional_stability_flag': bidirectional_stability_flag,
            'flow_concentration_index': flow_concentration,
            'grade_exchange_coherence': coherence,
            'separation_oscillation_indicator': separation_oscillation,
            'derived_phase_label': derived_label,
            'low_beta_base_topology': forward_lookup[(low_beta, base_resolution)],
            'low_beta_refined_topology': forward_lookup[(low_beta, refined_resolution)],
            'high_beta_base_topology': forward_lookup[(high_beta, base_resolution)],
            'high_beta_refined_topology': forward_lookup[(high_beta, refined_resolution)],
            'notes': str(run['notes']),
        }
        detail = {
            'run_id': str(run['run_id']),
            'phase_label': str(run['phase_label']),
            'phase_offset_fraction_of_pi': float(run['phase_offset_fraction_of_pi']),
            'width_regime': str(run['width_regime']),
            'width_ratio_a_to_b': float(run['width_ratio_a_to_b']),
            'dominant_forward_topology': dominant_forward,
            'dominant_reverse_topology': dominant_reverse,
            'derived_phase_label': derived_label,
            'variant_checks': [
                {
                    'beta': key[0],
                    'resolution': key[1],
                    'dt': float(value['dt']),
                    'steps': int(value['steps']),
                    'forward_topology_class': str(value['forward']['topology_class']),
                    'reverse_topology_class': str(value['reverse']['topology_class']),
                    'forward_survival_time': float(value['forward_survival_time']),
                    'reverse_survival_time': float(value['reverse_survival_time']),
                    'flow_concentration_index': float(value['forward']['flow_concentration_index']),
                    'grade_exchange_coherence': float(value['forward']['grade_exchange_coherence']),
                    'separation_oscillation_indicator': float(value['separation_oscillation_indicator']),
                }
                for key, value in sorted(variant_lookup.items(), key=lambda item: (item[0][0], item[0][1]))
            ],
        }
        rows.append(row)
        details.append(detail)

    components = stable_components(rows, phase_values, width_values)
    stable_connected_region_present = any(int(component['size']) >= 2 for component in components)
    stable_singletons = sum(1 for component in components if int(component['size']) == 1)
    boundary_segment_count = count_boundary_segments(rows, phase_values, width_values)
    finite_boundary_flag = bool(stable_connected_region_present and stable_singletons == 0 and boundary_segment_count > 0)
    bidirectional_preservation_flag = bool(
        all(int(row['bidirectional_stability_flag']) == 1 for row in rows if str(row['derived_phase_label']) != 'transient_mixed_phase')
    )
    phase_diagram_closed = bool(stable_connected_region_present and finite_boundary_flag and bidirectional_preservation_flag)

    stable_phase_counts = dict(Counter(str(row['derived_phase_label']) for row in rows))
    dominant_topology_counts = dict(Counter(str(row['topology_class']) for row in rows))
    if phase_diagram_closed and int(stable_phase_counts.get('smeared_transfer_phase', 0)) > 0 and int(stable_phase_counts.get('stable_braid_phase', 0)) == 0:
        closure_classification = 'smeared_dominant_closed_phase_diagram'
    elif phase_diagram_closed:
        closure_classification = 'closed_effective_phase_diagram'
    else:
        closure_classification = 'exploratory_only'

    phase_map_path = WORK_PLOT_DIR / 'stage23_12_effective_phase_map.png'
    survival_heatmap_path = WORK_PLOT_DIR / 'stage23_12_topology_survival_heatmap.png'
    stability_table_path = WORK_PLOT_DIR / 'stage23_12_refinement_stability_table.png'
    plot_phase_map(phase_map_path, rows, phase_values, width_labels, width_values)
    plot_survival_heatmap(survival_heatmap_path, rows, phase_values, width_labels, width_values)
    plot_stability_table(stability_table_path, rows, phase_values, width_labels, width_values)
    plot_paths = [phase_map_path, survival_heatmap_path, stability_table_path]

    summary = {
        'phase_diagram_closed': phase_diagram_closed,
        'closure_classification': closure_classification,
        'stable_connected_region_present': stable_connected_region_present,
        'finite_boundary_flag': finite_boundary_flag,
        'bidirectional_preservation_flag': bidirectional_preservation_flag,
        'boundary_segment_count': int(boundary_segment_count),
        'stable_phase_counts': stable_phase_counts,
        'dominant_topology_counts': dominant_topology_counts,
        'stable_components': components,
    }
    boundary_lines = minimal_boundary_summary(components, rows)

    result = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'common_fields': common,
        'cells': details,
        'summary': summary,
    }
    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        'stage23_12_dk_effective_phase_diagram_closure',
        result,
        rows,
        CSV_FIELDS,
        plot_paths,
    )
    write_note(json_path, csv_path, rows, stamped_plots, summary, boundary_lines)


if __name__ == '__main__':
    main()
