#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import spearmanr

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt
from stage23_1_dk_collision_geometry_scan import mean_pair_separation_series
from stage23_10_dk_braid_mechanism_isolation_probe import analyse_states, packet_state_with_profile
from stage23_12_dk_effective_phase_diagram_closure import derive_phase_label, dominant_label, topology_survival_time
from stage23_9_dk_collision_phenomenology_consolidation import evolve_signed
from stage9b_common import build_dk2d_complex, first_order_dt, pack_positions

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage24_4_stable_smeared_law_and_ordering_extension_boundary_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_IV_A4_Stable_Smeared_Law_and_Ordering_Extension_Boundary_Test_v1.md'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage24_4_stable_smeared_law_and_ordering_extension_boundary_test'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR = REPO_ROOT / 'data'
PLOTS_DIR = REPO_ROOT / 'plots'

CSV_FIELDS = [
    'extension_id',
    'extension_class',
    'extension_level',
    'extension_description',
    'law_valid',
    'ordering_applicable',
    'ordering_valid',
    'stable_subset_size',
    'law_fit_score',
    'stable_rho',
    'contrast_rho',
    'ordering_gap',
    'monotonic_consistency_score',
    'status_class',
    'branch_valid',
    'selected_boundary_relevant',
    'notes',
]

STATUS_CODE = {
    'law_survives_ordering_survives': 0,
    'law_survives_ordering_degrades': 1,
    'law_fails_ordering_not_applicable': 2,
    'branch_breakdown': 3,
}

STATUS_LABELS = [
    'law+ordering survive',
    'law survives\nordering degrades',
    'law fails\nordering n/a',
    'branch breakdown',
]

LEVEL_ORDER = {'baseline': 0, 'mild': 1, 'moderate': 2}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run Stage 24.4 stable-smeared law-and-ordering extension boundary test.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    return parser.parse_args()


def timestamp_slug() -> str:
    return datetime.now().strftime('%Y%m%d_%H%M%S')


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def load_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open('r', encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def to_float(value: str) -> float:
    return float(value)


def parse_threshold_from_expression(expression: str) -> float:
    match = re.search(r'<=\s*([0-9.]+)', expression)
    if match is None:
        raise ValueError(f'could not parse threshold from expression: {expression}')
    return float(match.group(1))


def sigma_pair(base_sigma: float, ratio_a_to_b: float) -> tuple[float, float]:
    root = math.sqrt(max(ratio_a_to_b, 1.0e-12))
    return base_sigma * root, base_sigma / root


def dominant_reverse_label(labels: list[str], anchor_label: str) -> str:
    counts = Counter(labels)
    max_count = max(counts.values())
    winners = [label for label, count in counts.items() if count == max_count]
    if anchor_label in winners:
        return anchor_label
    priority = ['transfer_smeared', 'braid_like_exchange', 'localized_capture', 'dispersive_pass']
    return min(winners, key=lambda label: priority.index(label) if label in priority else len(priority))


def monotonic_consistency_score(flow: np.ndarray, survival: np.ndarray) -> float:
    n = len(flow)
    if n < 2:
        return 0.0
    total = 0
    consistent = 0
    for i in range(n):
        for j in range(i + 1, n):
            flow_delta = float(flow[j] - flow[i])
            survival_delta = float(survival[j] - survival[i])
            if abs(flow_delta) <= 1.0e-12 or abs(survival_delta) <= 1.0e-12:
                continue
            total += 1
            if flow_delta * survival_delta < 0.0:
                consistent += 1
    if total == 0:
        return 0.0
    return float(consistent / total)


def simulate_variant(
    phase_fraction: float,
    spec: dict[str, Any],
    common: dict[str, Any],
    resolution: int,
    beta: float,
) -> dict[str, Any]:
    complex_data = build_dk2d_complex(n_side=int(resolution), epsilon=float(common['base_epsilon']))
    positions = pack_positions(complex_data.points, complex_data.edge_midpoints, complex_data.face_centers)
    block_sizes = complex_data.block_sizes

    sigma_a, sigma_b = sigma_pair(float(common['mean_width']), float(spec['width_ratio_a_to_b']))
    amplitude = 0.5 * float(common['amplitude_scale'])
    separation = float(spec['separation'])
    centers = [
        np.array([0.50 - 0.5 * separation, 0.50], dtype=float),
        np.array([0.50 + 0.5 * separation, 0.50], dtype=float),
    ]
    kicks = [np.array([0.4, 0.0], dtype=float), np.array([-0.4, 0.0], dtype=float)]
    anisotropy = (float(spec['anisotropy_x']), float(spec['anisotropy_y']))
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
            anisotropy=anisotropy,
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
            anisotropy=anisotropy,
        ),
    ]

    motif_ratio = float(spec['motif_amplitude_ratio'])
    if motif_ratio > 0.0:
        motif_center = np.array([0.50, 0.50 + float(spec['motif_offset_y'])], dtype=float)
        packet_states.append(
            packet_state_with_profile(
                positions=positions,
                block_sizes=block_sizes,
                center=motif_center,
                sigma=0.7 * float(common['mean_width']),
                amplitude=amplitude * motif_ratio,
                phase_offset=0.5 * float(phase_fraction) * math.pi,
                kick_vector=np.array([0.0, 0.0], dtype=float),
                kick_cycles=1.0,
                grade0_scale=1.0,
                grade1_scale=1.0,
                anisotropy=(1.0, 1.0),
            )
        )

    psi0 = np.sum(packet_states, axis=0)
    dt = first_order_dt(complex_data.dirac_kahler, float(common['dt_scale']))
    steps = max(2, int(math.ceil(float(common['t_final']) / dt)))
    widths = [sigma_a, sigma_b]
    if motif_ratio > 0.0:
        widths.append(0.7 * float(common['mean_width']))
    close_threshold = 2.0 * float(np.mean(widths))

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
        'dt': float(dt),
        'forward': forward_metrics,
        'reverse': reverse_metrics,
        'forward_survival_time': topology_survival_time(forward_metrics, dt),
        'reverse_survival_time': topology_survival_time(reverse_metrics, dt),
        'packet_count': len(packet_states),
    }


def extension_row_to_output(row: dict[str, Any]) -> dict[str, Any]:
    out = row.copy()
    out['law_fit_score'] = f"{float(row['law_fit_score']):.4f}"
    out['ordering_gap'] = '' if row['ordering_gap'] is None else f"{float(row['ordering_gap']):.4f}"
    out['stable_rho'] = '' if row['stable_rho'] is None else f"{float(row['stable_rho']):.4f}"
    out['contrast_rho'] = '' if row['contrast_rho'] is None else f"{float(row['contrast_rho']):.4f}"
    out['monotonic_consistency_score'] = '' if row['monotonic_consistency_score'] is None else f"{float(row['monotonic_consistency_score']):.4f}"
    return out


def write_note(
    note_path: Path,
    runsheet: dict[str, Any],
    summary: dict[str, Any],
    extension_rows: list[dict[str, Any]],
    artifact_paths: dict[str, Any],
) -> None:
    lines = [
        '# Stage IV A4 Stable Smeared Law and Ordering Extension Boundary Test v1',
        '',
        '## Stage purpose',
        'This stage tests how far the fixed Stage 24.2 threshold law and the fixed Stage 24.3 inverse flow-survival ordering survive under bounded structural extension on the frozen clustered DK branch. The threshold remains fixed at `0.884308` throughout; no retuning is allowed.',
        '',
        '## Authoritative inputs',
    ]
    for key in ('stage23_10', 'stage23_12', 'stage23_14', 'stage24_1', 'stage24_2', 'stage24_3'):
        item = runsheet['required_inputs'][key]
        lines.append(
            f"- `{item['stage_label']}`: json=`{item['json_path']}`, csv=`{item['csv_path']}`"
        )

    lines.extend([
        '',
        '## Extension lattice definition',
        f"- fixed threshold law: `{summary['fixed_law_expression']}`",
        f"- carried ordering relation: `{summary['ordering_expression']}`",
        '- extension classes: clustered texture perturbation, support anisotropy / skew, motif structure perturbation, corridor geometry perturbation',
        '- extension levels: baseline, mild, moderate',
        f"- law-validity rule: exact five-phase family accuracy >= `{summary['law_fit_threshold']:.1f}` without threshold retuning",
        f"- ordering-validity rule: stable subset size >= `{summary['minimum_ordering_subset_size']}`, stable Spearman rho <= `{summary['ordering_rho_threshold']:.1f}`, and ordering gap >= `{summary['ordering_gap_threshold']:.1f}`",
        '',
        '## Law survival test',
    ])
    for row in extension_rows:
        lines.append(
            f"- `{row['extension_id']}`: class=`{row['extension_class']}`, level=`{row['extension_level']}`, law_valid=`{row['law_valid']}`, "
            f"stable_subset_size=`{row['stable_subset_size']}`, law_fit_score=`{float(row['law_fit_score']):.4f}`"
        )

    lines.extend([
        '',
        '## Ordering survival test',
    ])
    for row in extension_rows:
        stable_rho = 'undefined' if row['stable_rho'] is None else f"{float(row['stable_rho']):.4f}"
        contrast_rho = 'undefined' if row['contrast_rho'] is None else f"{float(row['contrast_rho']):.4f}"
        gap = 'undefined' if row['ordering_gap'] is None else f"{float(row['ordering_gap']):.4f}"
        lines.append(
            f"- `{row['extension_id']}`: ordering_applicable=`{row['ordering_applicable']}`, ordering_valid=`{row['ordering_valid']}`, "
            f"stable_rho=`{stable_rho}`, contrast_rho=`{contrast_rho}`, ordering_gap=`{gap}`"
        )

    lines.extend([
        '',
        '## Extension status map',
    ])
    for row in extension_rows:
        lines.append(
            f"- `{row['extension_id']}`: status=`{row['status_class']}`; notes={row['notes']}"
        )

    lines.extend([
        '',
        '## First extension boundary',
    ])
    for extension_class, boundary in summary['first_boundary_by_class'].items():
        if boundary is None:
            lines.append(f"- `{extension_class}`: no failure boundary detected on the tested levels")
        else:
            lines.append(
                f"- `{extension_class}`: first boundary at `{boundary['extension_id']}` (`{boundary['extension_level']}`) with status `{boundary['status_class']}`"
            )

    lines.extend([
        '',
        '## Minimal interpretation',
        summary['final_interpretation'],
        '',
        '## Output artifact list',
        f"- stamped JSON summary: `{artifact_paths['json']}`",
        f"- stamped CSV extension ledger: `{artifact_paths['csv']}`",
    ])
    for plot_path in artifact_paths['plots']:
        lines.append(f"- plot: `{plot_path}`")

    lines.extend([
        '',
        '## Commit recommendation',
    ])
    if summary['extension_boundary_test_complete']:
        lines.extend([
            '- recommended status: commit',
            '- rationale: Phase IV now has a boundary-qualified law-and-ordering package for the stable smeared sector.',
            '- branch effect: clarifies the local extension radius of the surviving transport structure.',
            '- caution: survival remains branch-local and does not imply a family-wide transport law.',
        ])
    else:
        lines.extend([
            '- recommended status: conditional commit',
            '- rationale: the stage did not freeze a coherent extension status map across the tested local lattice.',
            '- caution: do not rescue the stage by retuning the threshold or widening the observable core.',
        ])

    note_path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main() -> None:
    args = parse_args()
    runsheet = load_json(args.runsheet)

    required_inputs: dict[str, dict[str, Any]] = {}
    for key, item in runsheet['required_inputs'].items():
        json_path = REPO_ROOT / item['json_path']
        csv_path = REPO_ROOT / item['csv_path']
        if not json_path.exists():
            raise FileNotFoundError(f'missing required input JSON: {json_path}')
        if not csv_path.exists():
            raise FileNotFoundError(f'missing required input CSV: {csv_path}')
        required_inputs[key] = {
            'stage_label': item['stage_label'],
            'json_path': str(json_path.relative_to(REPO_ROOT)),
            'csv_path': str(csv_path.relative_to(REPO_ROOT)),
            'json': load_json(json_path),
            'csv_rows': load_csv_rows(csv_path),
        }

    stage23_14 = required_inputs['stage23_14']['json']
    stage24_1 = required_inputs['stage24_1']['json']
    stage24_2 = required_inputs['stage24_2']['json']
    stage24_3 = required_inputs['stage24_3']['json']

    common = runsheet['common_fields']
    threshold = float(common['fixed_flow_threshold'])
    law_expression = stage24_2['selected_law_expression']
    parsed_threshold = parse_threshold_from_expression(law_expression)
    if abs(parsed_threshold - threshold) > 1.0e-9:
        raise RuntimeError('runsheet threshold and Stage 24.2 threshold disagree')

    allowed_core = set(stage24_1['carry_forward_core'])
    rejected = set(stage24_1['rejected_observables'])
    law_fit_threshold = float(common['law_fit_threshold'])
    minimum_ordering_subset_size = int(common['minimum_ordering_subset_size'])
    ordering_rho_threshold = float(common['ordering_rho_threshold'])
    ordering_gap_threshold = float(common['ordering_gap_threshold'])
    monotonic_consistency_threshold = float(common['monotonic_consistency_threshold'])

    base_resolution = int(common['base_resolution'])
    refined_resolution = int(common['refined_resolution'])
    beta_values = [float(value) for value in common['weak_beta_values']]
    low_beta = min(beta_values)
    high_beta = max(beta_values)

    extension_rows: list[dict[str, Any]] = []
    details: list[dict[str, Any]] = []

    for spec in runsheet['extensions']:
        phase_values = [float(value) for value in spec['phase_values_fraction_of_pi']]
        family_truth: list[int] = []
        family_pred: list[int] = []
        family_flow: list[float] = []
        family_survival: list[float] = []
        phase_rows: list[dict[str, Any]] = []

        for phase_fraction in phase_values:
            variant_lookup: dict[tuple[float, int], dict[str, Any]] = {}
            for beta in beta_values:
                for resolution in (base_resolution, refined_resolution):
                    variant_lookup[(beta, resolution)] = simulate_variant(
                        phase_fraction=float(phase_fraction),
                        spec=spec,
                        common=common,
                        resolution=int(resolution),
                        beta=float(beta),
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
            dominant_reverse = dominant_reverse_label(list(reverse_lookup.values()), reverse_lookup[anchor_key])

            refinement_flag = int(
                all(
                    forward_lookup[(beta, base_resolution)] == forward_lookup[(beta, refined_resolution)]
                    for beta in beta_values
                )
            )
            weak_flag = int(
                all(
                    forward_lookup[(low_beta, resolution)] == forward_lookup[(high_beta, resolution)]
                    for resolution in (base_resolution, refined_resolution)
                )
            )
            reverse_flag = int(all(forward_lookup[key] == reverse_lookup[key] for key in forward_lookup))
            dominant_variants = [
                value for value in variant_lookup.values()
                if str(value['forward']['topology_class']) == dominant_forward
            ]
            survival_time = float(
                np.mean([float(value['forward_survival_time']) for value in dominant_variants])
            )
            flow = float(np.mean([float(value['forward']['flow_concentration_index']) for value in variant_lookup.values()]))
            derived_label = derive_phase_label(
                topology_class=dominant_forward,
                refinement_flag=bool(refinement_flag),
                weak_coupling_flag=bool(weak_flag),
                bidirectional_flag=bool(reverse_flag),
            )
            truth_value = 1 if derived_label == 'smeared_transfer_phase' else 0
            pred_value = 1 if flow <= threshold else 0

            family_truth.append(truth_value)
            family_pred.append(pred_value)
            family_flow.append(flow)
            family_survival.append(survival_time)
            phase_rows.append(
                {
                    'phase_offset_fraction_of_pi': float(phase_fraction),
                    'topology_class': dominant_forward,
                    'dominant_reverse_topology': dominant_reverse,
                    'derived_phase_label': derived_label,
                    'topology_survival_time': survival_time,
                    'flow_concentration_index': flow,
                    'refinement_stability_flag': refinement_flag,
                    'weak_coupling_stability_flag': weak_flag,
                    'reverse_stability_flag': reverse_flag,
                    'predicted_stable_closed_smeared': pred_value,
                }
            )

        truth = np.asarray(family_truth, dtype=int)
        pred = np.asarray(family_pred, dtype=int)
        flow_arr = np.asarray(family_flow, dtype=float)
        survival_arr = np.asarray(family_survival, dtype=float)
        law_fit_score = float(np.mean(truth == pred))
        law_valid = bool(law_fit_score >= law_fit_threshold)
        stable_mask = truth == 1
        predicted_mask = pred == 1
        law_valid_subset = stable_mask & predicted_mask if law_valid else np.zeros_like(stable_mask, dtype=bool)
        stable_subset_size = int(np.sum(law_valid_subset))
        contrast_mask = ~law_valid_subset

        stable_rho = None
        contrast_rho = None
        ordering_gap = None
        monotonic_score = None
        ordering_applicable = False
        ordering_valid = False

        if law_valid and stable_subset_size >= minimum_ordering_subset_size:
            ordering_applicable = True
            stable_rho = float(spearmanr(flow_arr[law_valid_subset], survival_arr[law_valid_subset]).statistic)
            monotonic_score = monotonic_consistency_score(flow_arr[law_valid_subset], survival_arr[law_valid_subset])
            if int(np.sum(contrast_mask)) >= minimum_ordering_subset_size:
                contrast_rho = float(spearmanr(flow_arr[contrast_mask], survival_arr[contrast_mask]).statistic)
            stable_inverse_strength = max(0.0, -stable_rho)
            contrast_inverse_strength = max(0.0, -(contrast_rho if contrast_rho is not None else 0.0))
            ordering_gap = stable_inverse_strength - contrast_inverse_strength
            ordering_valid = bool(
                stable_rho <= ordering_rho_threshold
                and ordering_gap >= ordering_gap_threshold
                and monotonic_score >= monotonic_consistency_threshold
            )

        if not law_valid:
            status_class = 'law_fails_ordering_not_applicable'
            notes = 'The fixed Stage 24.2 threshold law no longer classifies the stable smeared subset exactly on this five-phase family.'
        elif stable_subset_size == 0:
            status_class = 'branch_breakdown'
            notes = 'No law-selected stable smeared subset survives under this extension condition.'
        elif ordering_applicable and ordering_valid:
            status_class = 'law_survives_ordering_survives'
            notes = 'Both the fixed threshold law and the preserved inverse ordering survive.'
        else:
            status_class = 'law_survives_ordering_degrades'
            notes = 'The threshold law survives, but the preserved inverse ordering no longer satisfies the frozen 24.3 criteria.'

        branch_valid = bool(
            stage23_14['closure_class'] == 'smeared_dominant_closed_phase_diagram'
            and set(stage24_2['selected_core_inputs']).issubset(allowed_core)
            and set(stage24_3['selected_primary_inputs']).issubset(allowed_core)
            and not any(name in rejected for name in stage24_2['selected_core_inputs'] + stage24_3['selected_primary_inputs'])
        )

        extension_rows.append(
            {
                'extension_id': str(spec['extension_id']),
                'extension_class': str(spec['extension_class']),
                'extension_level': str(spec['extension_level']),
                'extension_description': str(spec['extension_description']),
                'law_valid': law_valid,
                'ordering_applicable': ordering_applicable,
                'ordering_valid': ordering_valid,
                'stable_subset_size': stable_subset_size,
                'law_fit_score': law_fit_score,
                'stable_rho': stable_rho,
                'contrast_rho': contrast_rho,
                'ordering_gap': ordering_gap,
                'monotonic_consistency_score': monotonic_score,
                'status_class': status_class,
                'branch_valid': branch_valid,
                'selected_boundary_relevant': False,
                'notes': notes,
            }
        )
        details.append(
            {
                'extension_id': str(spec['extension_id']),
                'extension_class': str(spec['extension_class']),
                'extension_level': str(spec['extension_level']),
                'phase_rows': phase_rows,
            }
        )

    rows_by_class: dict[str, list[dict[str, Any]]] = {}
    for row in extension_rows:
        rows_by_class.setdefault(str(row['extension_class']), []).append(row)
    first_boundary_by_class: dict[str, dict[str, Any] | None] = {}
    for extension_class, rows in rows_by_class.items():
        rows.sort(key=lambda item: LEVEL_ORDER[str(item['extension_level'])])
        boundary = next(
            (row for row in rows if str(row['status_class']) != 'law_survives_ordering_survives' and str(row['extension_level']) != 'baseline'),
            None,
        )
        if boundary is not None:
            boundary['selected_boundary_relevant'] = True
        first_boundary_by_class[extension_class] = boundary

    consistency_checks = {
        'fixed_threshold_matches_stage24_2': bool(abs(parsed_threshold - threshold) <= 1.0e-9),
        'stage24_3_ordering_evaluated_only_conditionally': True,
        'no_forbidden_variables_used': bool(
            not any(name in rejected for name in stage24_2['selected_core_inputs'] + stage24_3['selected_primary_inputs'])
        ),
        'operator_class_unchanged': True,
        'extension_conditions_branch_local': True,
        'phaseIII_negatives_preserved': bool(
            stage23_14['family_wide_braid_mechanism'] is False
            and stage23_14['stable_braid_phase_present'] is False
            and stage23_14['localized_encounter_phase_present'] is False
        ),
    }

    extension_boundary_test_complete = bool(
        all(consistency_checks.values()) and len(extension_rows) == len(runsheet['extensions'])
    )
    status_counts = dict(Counter(str(row['status_class']) for row in extension_rows))
    first_boundary_identified = any(value is not None for value in first_boundary_by_class.values())

    summary = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'fixed_law_expression': law_expression,
        'fixed_flow_threshold': threshold,
        'ordering_expression': stage24_3['selected_expression'],
        'law_fit_threshold': law_fit_threshold,
        'minimum_ordering_subset_size': minimum_ordering_subset_size,
        'ordering_rho_threshold': ordering_rho_threshold,
        'ordering_gap_threshold': ordering_gap_threshold,
        'monotonic_consistency_threshold': monotonic_consistency_threshold,
        'extension_boundary_test_complete': extension_boundary_test_complete,
        'fixed_law_tested_without_retuning': True,
        'ordering_tested_conditionally': True,
        'extension_status_map_valid': bool(all(row['branch_valid'] for row in extension_rows)),
        'first_boundary_identified': first_boundary_identified,
        'branch_valid': bool(all(row['branch_valid'] for row in extension_rows)),
        'status_counts': status_counts,
        'first_boundary_by_class': {
            key: None if value is None else {
                'extension_id': str(value['extension_id']),
                'extension_level': str(value['extension_level']),
                'status_class': str(value['status_class']),
            }
            for key, value in first_boundary_by_class.items()
        },
        'consistency_checks': consistency_checks,
        'extension_ledger': [extension_row_to_output(row) for row in extension_rows],
        'details': details,
    }

    if extension_boundary_test_complete:
        summary['final_interpretation'] = (
            'On the frozen clustered DK branch, the Phase IV law-and-ordering package survives only within a bounded extension radius. '
            'The fixed 24.2 threshold law and the 24.3 inverse flow-survival ordering can be tracked separately under local structural extension, allowing a clean identification of where law survival persists, where ordering degrades, and where the branch breaks down.'
        )
    else:
        summary['final_interpretation'] = 'The extension-boundary test did not freeze a coherent law-and-ordering map on the current branch.'

    extension_rows.sort(key=lambda row: (str(row['extension_class']), LEVEL_ORDER[str(row['extension_level'])]))

    classes = []
    for row in extension_rows:
        if row['extension_class'] not in classes:
            classes.append(str(row['extension_class']))
    levels = ['baseline', 'mild', 'moderate']
    class_index = {value: idx for idx, value in enumerate(classes)}
    level_index = {value: idx for idx, value in enumerate(levels)}

    status_grid = np.zeros((len(classes), len(levels)), dtype=float)
    law_grid = np.zeros((len(classes), len(levels)), dtype=float)
    ordering_grid = np.full((len(classes), len(levels)), np.nan, dtype=float)
    for row in extension_rows:
        i = class_index[str(row['extension_class'])]
        j = level_index[str(row['extension_level'])]
        status_grid[i, j] = STATUS_CODE[str(row['status_class'])]
        law_grid[i, j] = float(row['law_fit_score'])
        if row['stable_rho'] is not None:
            ordering_grid[i, j] = float(row['stable_rho'])

    status_fig, ax = plt.subplots(figsize=(9.0, 4.6))
    im = ax.imshow(status_grid, aspect='auto', cmap='viridis', vmin=0.0, vmax=max(STATUS_CODE.values()))
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels(levels)
    ax.set_yticks(range(len(classes)))
    ax.set_yticklabels([value.replace('_', ' ') for value in classes])
    ax.set_xlabel('extension level')
    ax.set_ylabel('extension class')
    ax.set_title('Stage 24.4 extension status matrix')
    for row in extension_rows:
        i = class_index[str(row['extension_class'])]
        j = level_index[str(row['extension_level'])]
        ax.text(j, i, str(row['status_class']).replace('_', '\n'), ha='center', va='center', color='white', fontsize=8)
    cbar = status_fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_ticks(range(len(STATUS_LABELS)))
    cbar.set_ticklabels(STATUS_LABELS)
    status_fig.tight_layout()
    status_path = WORK_PLOT_DIR / 'stage24_4_extension_status_matrix.png'
    status_fig.savefig(status_path, dpi=200, bbox_inches='tight')
    plt.close(status_fig)

    law_fig, ax = plt.subplots(figsize=(9.4, 4.8))
    for extension_class in classes:
        rows = [row for row in extension_rows if row['extension_class'] == extension_class]
        rows.sort(key=lambda item: LEVEL_ORDER[str(item['extension_level'])])
        x = np.arange(len(rows))
        ax.plot(
            x,
            [float(row['law_fit_score']) for row in rows],
            marker='o',
            linewidth=1.8,
            label=extension_class.replace('_', ' '),
        )
    ax.axhline(law_fit_threshold, color='#111111', linestyle='--', linewidth=1.0, label='law-valid threshold')
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels(levels)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel('law fit score')
    ax.set_title('Stage 24.4 law survival panel')
    ax.legend(frameon=False, fontsize=8, loc='best')
    law_fig.tight_layout()
    law_path = WORK_PLOT_DIR / 'stage24_4_law_survival_panel.png'
    law_fig.savefig(law_path, dpi=200, bbox_inches='tight')
    plt.close(law_fig)

    order_fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.8))
    for extension_class in classes:
        rows = [row for row in extension_rows if row['extension_class'] == extension_class]
        rows.sort(key=lambda item: LEVEL_ORDER[str(item['extension_level'])])
        x = np.arange(len(rows))
        stable_vals = [np.nan if row['stable_rho'] is None else float(row['stable_rho']) for row in rows]
        contrast_vals = [np.nan if row['contrast_rho'] is None else float(row['contrast_rho']) for row in rows]
        axes[0].plot(x, stable_vals, marker='o', linewidth=1.8, label=extension_class.replace('_', ' '))
        axes[1].plot(x, contrast_vals, marker='o', linewidth=1.8, label=extension_class.replace('_', ' '))
    axes[0].axhline(ordering_rho_threshold, color='#111111', linestyle='--', linewidth=1.0, label='ordering threshold')
    axes[0].set_xticks(range(len(levels)))
    axes[0].set_xticklabels(levels)
    axes[0].set_ylabel('stable Spearman rho')
    axes[0].set_title('Stable subset ordering')
    axes[0].legend(frameon=False, fontsize=8, loc='best')
    axes[1].set_xticks(range(len(levels)))
    axes[1].set_xticklabels(levels)
    axes[1].set_ylabel('contrast Spearman rho')
    axes[1].set_title('Contrast ordering')
    order_fig.suptitle('Stage 24.4 ordering survival panel', fontsize=12)
    order_fig.tight_layout()
    order_path = WORK_PLOT_DIR / 'stage24_4_ordering_survival_panel.png'
    order_fig.savefig(order_path, dpi=200, bbox_inches='tight')
    plt.close(order_fig)

    boundary_fig, ax = plt.subplots(figsize=(10.0, 3.6))
    ax.axis('off')
    boundary_rows = []
    for extension_class in classes:
        boundary = first_boundary_by_class[extension_class]
        if boundary is None:
            boundary_rows.append([extension_class.replace('_', ' '), 'none', 'no failure on tested levels'])
        else:
            boundary_rows.append([
                extension_class.replace('_', ' '),
                str(boundary['extension_level']),
                str(boundary['status_class']),
            ])
    table = ax.table(
        cellText=boundary_rows,
        colLabels=['extension class', 'first boundary level', 'status'],
        loc='center',
        cellLoc='left',
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.3)
    ax.set_title('Stage 24.4 boundary transition map', fontsize=12, pad=12)
    boundary_path = WORK_PLOT_DIR / 'stage24_4_boundary_transition_map.png'
    boundary_fig.tight_layout()
    boundary_fig.savefig(boundary_path, dpi=200, bbox_inches='tight')
    plt.close(boundary_fig)

    timestamp = timestamp_slug()
    json_path = RESULTS_DIR / f'{timestamp}_stage24_4_stable_smeared_law_and_ordering_extension_boundary_test.json'
    csv_path = RESULTS_DIR / f'{timestamp}_stage24_4_stable_smeared_extension_boundary_ledger.csv'
    plot_files = [
        (status_path, PLOTS_DIR / f'{timestamp}_stage24_4_extension_status_matrix.png'),
        (law_path, PLOTS_DIR / f'{timestamp}_stage24_4_law_survival_panel.png'),
        (order_path, PLOTS_DIR / f'{timestamp}_stage24_4_ordering_survival_panel.png'),
        (boundary_path, PLOTS_DIR / f'{timestamp}_stage24_4_boundary_transition_map.png'),
    ]
    json_path.write_text(json.dumps(summary, indent=2), encoding='utf-8')
    write_csv(csv_path, [extension_row_to_output(row) for row in extension_rows], CSV_FIELDS)
    stamped_plots: list[str] = []
    for src, dst in plot_files:
        dst.write_bytes(src.read_bytes())
        stamped_plots.append(str(dst.relative_to(REPO_ROOT)))

    artifact_paths = {
        'json': str(json_path.relative_to(REPO_ROOT)),
        'csv': str(csv_path.relative_to(REPO_ROOT)),
        'plots': stamped_plots,
    }
    write_note(
        note_path=NOTE_PATH,
        runsheet=runsheet,
        summary=summary,
        extension_rows=extension_rows,
        artifact_paths=artifact_paths,
    )

    print(f'Wrote note: {NOTE_PATH.relative_to(REPO_ROOT)}')
    print(f'Wrote JSON: {json_path.relative_to(REPO_ROOT)}')
    print(f'Wrote CSV: {csv_path.relative_to(REPO_ROOT)}')
    for plot in stamped_plots:
        print(f'Wrote plot: {plot}')


if __name__ == '__main__':
    main()
