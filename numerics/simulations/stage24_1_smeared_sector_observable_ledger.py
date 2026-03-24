#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt, save_atlas_payload

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage24_1_smeared_sector_observable_ledger_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_IV_A1_Smeared_Sector_Observable_Ledger_v1.md'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage24_1_smeared_sector_observable_ledger'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)

CSV_FIELDS = [
    'observable_id',
    'observable_label',
    'status',
    'allowed_into_phaseIV_core',
    'source_stage',
    'scope',
    'available_on_frozen_branch',
    'sensitivity_score',
    'stable_range',
    'contrast_range',
    'rationale',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run Stage 24.1 smeared-sector observable ledger.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    return parser.parse_args()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def load_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open('r', encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))


def to_bool(value: str) -> bool:
    return str(value).strip() in {'1', 'True', 'true'}


def to_float(value: str) -> float:
    return float(value)


def format_range(values: list[float]) -> str:
    if not values:
        return 'n/a'
    return f"{min(values):.4f} -> {max(values):.4f}"


def interval_sensitivity(stable_values: list[float], contrast_values: list[float]) -> float:
    if not stable_values or not contrast_values:
        return 0.0
    stable_min = min(stable_values)
    stable_max = max(stable_values)
    contrast_min = min(contrast_values)
    contrast_max = max(contrast_values)
    total_min = min(stable_min, contrast_min)
    total_max = max(stable_max, contrast_max)
    total_span = max(total_max - total_min, 1.0e-12)
    overlap = max(0.0, min(stable_max, contrast_max) - max(stable_min, contrast_min))
    return float(max(0.0, 1.0 - overlap / total_span))


def binary_sensitivity(stable_values: list[int], contrast_values: list[int]) -> float:
    if not stable_values or not contrast_values:
        return 0.0
    stable_mean = float(np.mean(stable_values))
    contrast_mean = float(np.mean(contrast_values))
    return float(abs(stable_mean - contrast_mean))


def build_stability_grids(rows: list[dict[str, str]]) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[float], list[float]]:
    width_values = sorted({to_float(row['width_ratio_a_to_b']) for row in rows})
    phase_values = sorted({to_float(row['phase_offset_fraction_of_pi']) for row in rows})
    width_index = {value: idx for idx, value in enumerate(width_values)}
    phase_index = {value: idx for idx, value in enumerate(phase_values)}

    refinement = np.zeros((len(width_values), len(phase_values)), dtype=float)
    weak = np.zeros((len(width_values), len(phase_values)), dtype=float)
    reverse = np.zeros((len(width_values), len(phase_values)), dtype=float)

    for row in rows:
        i = width_index[to_float(row['width_ratio_a_to_b'])]
        j = phase_index[to_float(row['phase_offset_fraction_of_pi'])]
        refinement[i, j] = 1.0 if to_bool(row['refinement_stability_flag']) else 0.0
        weak[i, j] = 1.0 if to_bool(row['weak_coupling_stability_flag']) else 0.0
        reverse[i, j] = 1.0 if to_bool(row['bidirectional_stability_flag']) else 0.0

    return refinement, weak, reverse, width_values, phase_values


def write_note(
    note_path: Path,
    runsheet: dict[str, Any],
    summary: dict[str, Any],
    ledger_rows: list[dict[str, Any]],
    artifact_paths: dict[str, str],
) -> None:
    primary_rows = [row for row in ledger_rows if row['status'] == 'primary_observable']
    secondary_rows = [row for row in ledger_rows if row['status'] == 'secondary_observable']
    derived_rows = [row for row in ledger_rows if row['status'] == 'derived']
    reject_rows = [row for row in ledger_rows if row['status'] == 'reject']

    lines = [
        '# Stage IV A1 Smeared Sector Observable Ledger v1',
        '',
        '## Stage purpose',
        'This note defines the minimal operational observable basis of the surviving smeared-transfer sector on the frozen clustered Dirac-Kahler branch. It is a reduction-only ledger extraction from the frozen Phase III closure outputs.',
        '',
        '## Frozen inputs',
    ]
    for key in ('stage23_11', 'stage23_12', 'stage23_13', 'stage23_14'):
        item = runsheet['required_inputs'][key]
        lines.append(
            f"- `{item['stage_label']}`: json=`{item['json_path']}`, csv=`{item['csv_path']}`"
        )

    lines.extend([
        '',
        '## Stable object and contrast policy',
        f"- stable smeared-sector cells = `{summary['stable_smeared_cell_count']}`",
        f"- contrast transient-mixed cells = `{summary['contrast_control_cell_count']}`",
        f"- stable phase span = `{summary['stable_phase_min']:.3f} -> {summary['stable_phase_max']:.3f}`",
        f"- derived phase-corridor width = `{summary['phase_corridor_width']:.3f}`",
        '- contrast controls are used only to test separation; they are not promoted into the Phase IV core.',
        '',
        '## Validation checks',
    ])
    for key, value in summary['validation_checks'].items():
        lines.append(f"- `{key}` = `{value}`")

    lines.extend([
        '',
        '## Primary observable basis',
    ])
    for row in primary_rows:
        lines.append(
            f"- `{row['observable_id']}`: sensitivity=`{row['sensitivity_score']}`; stable=`{row['stable_range']}`; contrast=`{row['contrast_range']}`; rationale={row['rationale']}"
        )

    lines.extend([
        '',
        '## Secondary and derived observables',
    ])
    for row in secondary_rows + derived_rows:
        lines.append(
            f"- `{row['observable_id']}`: status=`{row['status']}`; stable=`{row['stable_range']}`; contrast=`{row['contrast_range']}`; rationale={row['rationale']}"
        )

    lines.extend([
        '',
        '## Explicit reject table',
        '',
        '| Observable | Reason rejected |',
        '| --- | --- |',
    ])
    for row in reject_rows:
        lines.append(f"| `{row['observable_id']}` | {row['rationale']} |")

    lines.extend([
        '',
        '## Carry-forward set for 24.2',
        f"- core carry-forward observables: `{', '.join(summary['carry_forward_core'])}`",
        f"- context-only carry-forward observables: `{', '.join(summary['carry_forward_context'])}`",
        '- excluded braid-centered quantities remain outside the Phase IV core and may appear only as contrast descriptors if explicitly needed.',
        '',
        '## Final read',
        summary['final_read'],
        '',
        '## Output artifact list',
        f"- stamped JSON summary: `{artifact_paths['json']}`",
        f"- stamped CSV ledger: `{artifact_paths['csv']}`",
    ])
    for plot_path in artifact_paths['plots']:
        lines.append(f"- plot: `{plot_path}`")

    note_path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main() -> None:
    args = parse_args()
    runsheet = load_json(args.runsheet)

    required_inputs = {}
    for key, item in runsheet['required_inputs'].items():
        json_path = REPO_ROOT / item['json_path']
        csv_path = REPO_ROOT / item['csv_path']
        if not json_path.exists():
            raise FileNotFoundError(f"missing required input JSON: {json_path}")
        if not csv_path.exists():
            raise FileNotFoundError(f"missing required input CSV: {csv_path}")
        required_inputs[key] = {
            'stage_label': item['stage_label'],
            'json': load_json(json_path),
            'csv_rows': load_csv_rows(csv_path),
            'json_path': str(json_path.relative_to(REPO_ROOT)),
            'csv_path': str(csv_path.relative_to(REPO_ROOT)),
        }

    s11 = required_inputs['stage23_11']['json']
    s12 = required_inputs['stage23_12']['json']
    s12_rows = required_inputs['stage23_12']['csv_rows']
    s13_rows = required_inputs['stage23_13']['csv_rows']
    s14 = required_inputs['stage23_14']['json']

    stable_rows = [row for row in s12_rows if row['derived_phase_label'] == 'smeared_transfer_phase']
    contrast_rows = [row for row in s12_rows if row['derived_phase_label'] == 'transient_mixed_phase']
    if not stable_rows:
        raise RuntimeError('no stable smeared-transfer cells found in Stage 23.12 closure output')
    if not contrast_rows:
        raise RuntimeError('no transient mixed contrast cells found in Stage 23.12 closure output')

    stable_phase_values = [to_float(row['phase_offset_fraction_of_pi']) for row in stable_rows]
    phase_corridor_width = float(max(stable_phase_values) - min(stable_phase_values))

    stable_topology_survival = [to_float(row['topology_survival_time']) for row in stable_rows]
    contrast_topology_survival = [to_float(row['topology_survival_time']) for row in contrast_rows]
    stable_flow = [to_float(row['flow_concentration_index']) for row in stable_rows]
    contrast_flow = [to_float(row['flow_concentration_index']) for row in contrast_rows]
    stable_coherence = [to_float(row['grade_exchange_coherence']) for row in stable_rows]
    contrast_coherence = [to_float(row['grade_exchange_coherence']) for row in contrast_rows]
    stable_refinement = [1 if to_bool(row['refinement_stability_flag']) else 0 for row in stable_rows]
    contrast_refinement = [1 if to_bool(row['refinement_stability_flag']) else 0 for row in contrast_rows]
    stable_weak = [1 if to_bool(row['weak_coupling_stability_flag']) else 0 for row in stable_rows]
    contrast_weak = [1 if to_bool(row['weak_coupling_stability_flag']) else 0 for row in contrast_rows]
    stable_reverse = [1 if to_bool(row['bidirectional_stability_flag']) else 0 for row in stable_rows]
    contrast_reverse = [1 if to_bool(row['bidirectional_stability_flag']) else 0 for row in contrast_rows]

    phase_sensitivity = interval_sensitivity(stable_phase_values, [to_float(row['phase_offset_fraction_of_pi']) for row in contrast_rows])
    survival_sensitivity = interval_sensitivity(stable_topology_survival, contrast_topology_survival)
    flow_sensitivity = interval_sensitivity(stable_flow, contrast_flow)
    coherence_sensitivity = interval_sensitivity(stable_coherence, contrast_coherence)
    refinement_sensitivity = binary_sensitivity(stable_refinement, contrast_refinement)
    weak_sensitivity = binary_sensitivity(stable_weak, contrast_weak)
    reverse_sensitivity = binary_sensitivity(stable_reverse, contrast_reverse)

    braid_rows = [
        row for row in s13_rows
        if row['topology_class'] == 'braid_like_exchange'
    ]
    braid_indicator_sensitivity = binary_sensitivity(
        [0 for _ in stable_rows],
        [1 for _ in braid_rows] if braid_rows else [1],
    )

    ledger_rows: list[dict[str, Any]] = [
        {
            'observable_id': 'phase_corridor_position',
            'observable_label': 'Phase corridor position',
            'status': 'secondary_observable',
            'allowed_into_phaseIV_core': False,
            'source_stage': '23.12 / III-C3',
            'scope': 'context_coordinate',
            'available_on_frozen_branch': True,
            'sensitivity_score': f'{phase_sensitivity:.3f}',
            'stable_range': format_range(stable_phase_values),
            'contrast_range': format_range([to_float(row['phase_offset_fraction_of_pi']) for row in contrast_rows]),
            'rationale': 'Phase position organizes where the stable smeared row sits, but the same sampled phases also appear in transient mixed controls, so it is contextual rather than a primary separator.',
        },
        {
            'observable_id': 'phase_corridor_width',
            'observable_label': 'Phase corridor width',
            'status': 'derived',
            'allowed_into_phaseIV_core': False,
            'source_stage': '23.12 / III-C3',
            'scope': 'derived_boundary_summary',
            'available_on_frozen_branch': True,
            'sensitivity_score': 'n/a',
            'stable_range': f'{phase_corridor_width:.3f}',
            'contrast_range': 'n/a',
            'rationale': 'Corridor width is a boundary summary of the stable smeared component, not a per-cell operational observable.',
        },
        {
            'observable_id': 'topology_survival_time',
            'observable_label': 'Topology survival time',
            'status': 'primary_observable',
            'allowed_into_phaseIV_core': True,
            'source_stage': '23.12 / III-C3',
            'scope': 'stable_smeared_region',
            'available_on_frozen_branch': True,
            'sensitivity_score': f'{survival_sensitivity:.3f}',
            'stable_range': format_range(stable_topology_survival),
            'contrast_range': format_range(contrast_topology_survival),
            'rationale': 'This is the minimal persistence-carrying observable on the closed branch. It is not sufficient alone, but it is required to carry smeared-sector lifetime structure into 24.2.',
        },
        {
            'observable_id': 'transfer_asymmetry',
            'observable_label': 'Transfer asymmetry',
            'status': 'reject',
            'allowed_into_phaseIV_core': False,
            'source_stage': 'none',
            'scope': 'unavailable_on_frozen_branch',
            'available_on_frozen_branch': False,
            'sensitivity_score': 'n/a',
            'stable_range': 'n/a',
            'contrast_range': 'n/a',
            'rationale': 'No branch-stable transfer-asymmetry quantity is present in the frozen Phase III closure outputs without introducing a new bookkeeping transport model.',
        },
        {
            'observable_id': 'refinement_stability_flag',
            'observable_label': 'Refinement stability flag',
            'status': 'primary_observable',
            'allowed_into_phaseIV_core': True,
            'source_stage': '23.12 / III-C3',
            'scope': 'closure_support',
            'available_on_frozen_branch': True,
            'sensitivity_score': f'{refinement_sensitivity:.3f}',
            'stable_range': format_range([float(value) for value in stable_refinement]),
            'contrast_range': format_range([float(value) for value in contrast_refinement]),
            'rationale': 'Refinement stability is a necessary closure-support observable and perfectly separates the stable smeared cells from the transient mixed controls on the frozen lattice.',
        },
        {
            'observable_id': 'weak_coupling_stability_flag',
            'observable_label': 'Weak-coupling stability flag',
            'status': 'primary_observable',
            'allowed_into_phaseIV_core': True,
            'source_stage': '23.12 / III-C3',
            'scope': 'closure_support',
            'available_on_frozen_branch': True,
            'sensitivity_score': f'{weak_sensitivity:.3f}',
            'stable_range': format_range([float(value) for value in stable_weak]),
            'contrast_range': format_range([float(value) for value in contrast_weak]),
            'rationale': 'The stable smeared component is beta-stable by construction. This flag is not sufficient on its own, but it is indispensable because the closed regime is defined only after the weak-coupling consistency check is passed.',
        },
        {
            'observable_id': 'reverse_stability_flag',
            'observable_label': 'Reverse stability flag',
            'status': 'primary_observable',
            'allowed_into_phaseIV_core': True,
            'source_stage': '23.12 / III-C3',
            'scope': 'closure_support',
            'available_on_frozen_branch': True,
            'sensitivity_score': f'{reverse_sensitivity:.3f}',
            'stable_range': format_range([float(value) for value in stable_reverse]),
            'contrast_range': format_range([float(value) for value in contrast_reverse]),
            'rationale': 'Reverse stability is required by the Phase III closure rule and cleanly separates the stable smeared cells from the transient mixed controls.',
        },
        {
            'observable_id': 'coherence_decay_rate',
            'observable_label': 'Coherence decay rate',
            'status': 'reject',
            'allowed_into_phaseIV_core': False,
            'source_stage': 'none',
            'scope': 'unavailable_on_frozen_branch',
            'available_on_frozen_branch': False,
            'sensitivity_score': 'n/a',
            'stable_range': 'n/a',
            'contrast_range': 'n/a',
            'rationale': 'The frozen closure outputs retain coherence levels, not a regime-stable coherence-decay-rate observable. Promoting a decay metric here would require a new analysis layer not licensed by Phase III.',
        },
        {
            'observable_id': 'grade_exchange_coherence_level',
            'observable_label': 'Grade-exchange coherence level',
            'status': 'secondary_observable',
            'allowed_into_phaseIV_core': False,
            'source_stage': '23.12 / III-C3',
            'scope': 'descriptive_support',
            'available_on_frozen_branch': True,
            'sensitivity_score': f'{coherence_sensitivity:.3f}',
            'stable_range': format_range(stable_coherence),
            'contrast_range': format_range(contrast_coherence),
            'rationale': 'Coherence remains well-defined across the smeared sector, but its stable and contrast ranges overlap too strongly to justify promotion into the minimal Phase IV core.',
        },
        {
            'observable_id': 'flow_concentration_index',
            'observable_label': 'Flow concentration index',
            'status': 'primary_observable',
            'allowed_into_phaseIV_core': True,
            'source_stage': '23.12 / III-C3',
            'scope': 'stable_smeared_region',
            'available_on_frozen_branch': True,
            'sensitivity_score': f'{flow_sensitivity:.3f}',
            'stable_range': format_range(stable_flow),
            'contrast_range': format_range(contrast_flow),
            'rationale': 'This is the strongest continuous transport-like separator on the frozen branch: the stable smeared cells occupy a clean low-flow-concentration band disjoint from the transient mixed controls.',
        },
        {
            'observable_id': 'topology_return_error',
            'observable_label': 'Topology return error',
            'status': 'reject',
            'allowed_into_phaseIV_core': False,
            'source_stage': '23.13 / III-C4',
            'scope': 'ablation_specific_metric',
            'available_on_frozen_branch': True,
            'sensitivity_score': 'n/a',
            'stable_range': 'n/a',
            'contrast_range': 'n/a',
            'rationale': 'Return error exists only as a branch-specific ablation diagnostic in Stage 23.13 and does not define the stable smeared sector across the frozen closure lattice.',
        },
        {
            'observable_id': 'smeared_regime_membership_flag',
            'observable_label': 'Smeared regime membership flag',
            'status': 'derived',
            'allowed_into_phaseIV_core': False,
            'source_stage': '23.12 / III-C3',
            'scope': 'derived_target_label',
            'available_on_frozen_branch': True,
            'sensitivity_score': '1.000',
            'stable_range': '1',
            'contrast_range': '0',
            'rationale': 'This is the target label produced by the Phase III closure pass. It is necessary as an outcome variable, but it is not itself a measurable Phase IV observable.',
        },
        {
            'observable_id': 'braid_like_exchange_indicator',
            'observable_label': 'Braid-like exchange indicator',
            'status': 'reject',
            'allowed_into_phaseIV_core': False,
            'source_stage': '23.12 / III-C3; 23.14 / III-C5',
            'scope': 'phaseIII_excluded_claim_space',
            'available_on_frozen_branch': True,
            'sensitivity_score': f'{braid_indicator_sensitivity:.3f}',
            'stable_range': '0',
            'contrast_range': 'contrast-only',
            'rationale': 'Braid-centered quantities are explicitly excluded from promotion into the Phase IV core. They may appear only as contrast descriptors, because the family-wide braid claim and stable braid-phase claim failed in Phase III.',
        },
        {
            'observable_id': 'mechanism_sensitivity_score',
            'observable_label': 'Mechanism sensitivity score',
            'status': 'reject',
            'allowed_into_phaseIV_core': False,
            'source_stage': '23.13 / III-C4',
            'scope': 'ablation_specific_metric',
            'available_on_frozen_branch': True,
            'sensitivity_score': 'n/a',
            'stable_range': 'n/a',
            'contrast_range': 'n/a',
            'rationale': 'Mechanism sensitivity is tied to the ablation design of Stage 23.13 and does not characterize the closed smeared sector itself.',
        },
    ]

    primary_ids = [row['observable_id'] for row in ledger_rows if row['status'] == 'primary_observable']
    context_ids = ['phase_corridor_position', 'phase_corridor_width']

    validation_checks = {
        'stage23_11_effective_model_valid_true': bool(s11['summary']['effective_model_valid']),
        'stage23_12_phase_diagram_closed_true': bool(s12['summary']['phase_diagram_closed']),
        'stage23_12_stable_region_is_smeared_only': bool(
            s12['summary']['stable_phase_counts'].get('smeared_transfer_phase', 0) > 0
            and 'stable_braid_phase' not in s12['summary']['stable_phase_counts']
            and 'localized_encounter_phase' not in s12['summary']['stable_phase_counts']
        ),
        'stage23_13_phase_corridor_primary': bool(
            s13_rows[3]['branch_mechanism_label'] == 'phase_corridor_primary'
            or required_inputs['stage23_13']['json']['overall_summary']['overall_mechanism_label'] == 'phase_corridor_primary'
        ),
        'stage23_14_freeze_valid_true': bool(s14['phaseIII_freeze_valid']),
        'primary_basis_smaller_than_candidate_list': len(primary_ids) < len(runsheet['candidate_observables']),
        'primary_basis_distinguishes_stable_region': bool(
            all(
                to_bool(row['refinement_stability_flag'])
                and to_bool(row['weak_coupling_stability_flag'])
                and to_bool(row['bidirectional_stability_flag'])
                for row in stable_rows
            )
            and all(
                not (
                    to_bool(row['refinement_stability_flag'])
                    and to_bool(row['weak_coupling_stability_flag'])
                    and to_bool(row['bidirectional_stability_flag'])
                )
                for row in contrast_rows
            )
        ),
        'no_braid_quantity_promoted_to_primary': all(
            'braid' not in row['observable_id'] for row in ledger_rows if row['status'] == 'primary_observable'
        ),
    }

    observable_ledger_valid = all(validation_checks.values())

    summary = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'observable_ledger_valid': observable_ledger_valid,
        'stable_smeared_cell_count': len(stable_rows),
        'contrast_control_cell_count': len(contrast_rows),
        'stable_phase_min': min(stable_phase_values),
        'stable_phase_max': max(stable_phase_values),
        'phase_corridor_width': phase_corridor_width,
        'minimal_observable_basis': primary_ids,
        'carry_forward_core': primary_ids,
        'carry_forward_context': context_ids,
        'secondary_observables': [row['observable_id'] for row in ledger_rows if row['status'] == 'secondary_observable'],
        'derived_observables': [row['observable_id'] for row in ledger_rows if row['status'] == 'derived'],
        'rejected_observables': [row['observable_id'] for row in ledger_rows if row['status'] == 'reject'],
        'contrast_policy': runsheet['contrast_policy'],
        'validation_checks': validation_checks,
        'final_read': (
            'On the frozen clustered DK branch, the closed smeared-transfer regime is minimally characterized by '
            'refinement stability, weak-coupling stability, reverse stability, topology survival time, and flow concentration index, '
            'while phase-corridor position and grade-exchange coherence remain secondary or derived descriptors.'
        ),
        'candidate_ledger': ledger_rows,
    }

    table_fig, ax = plt.subplots(figsize=(12.5, 5.8))
    ax.axis('off')
    table_rows = [
        [
            row['observable_id'],
            row['status'].replace('_', ' '),
            'yes' if row['allowed_into_phaseIV_core'] else 'no',
            row['sensitivity_score'],
        ]
        for row in ledger_rows
    ]
    table = ax.table(
        cellText=table_rows,
        colLabels=['observable', 'status', 'phase IV core', 'sensitivity'],
        loc='center',
        cellLoc='left',
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.25)
    ax.set_title('Stage 24.1 smeared-sector observable ledger', fontsize=12, pad=12)
    table_path = WORK_PLOT_DIR / 'stage24_1_smeared_sector_observable_table.png'
    table_fig.tight_layout()
    table_fig.savefig(table_path, dpi=200, bbox_inches='tight')
    plt.close(table_fig)

    corr_fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    axes[0].scatter(contrast_topology_survival, contrast_flow, color='#8f8f8f', label='transient mixed controls', s=45)
    axes[0].scatter(stable_topology_survival, stable_flow, color='#c23b22', label='stable smeared cells', s=55)
    axes[0].set_xlabel('topology survival time')
    axes[0].set_ylabel('flow concentration index')
    axes[0].set_title('Primary transport separation')
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].scatter(
        [to_float(row['phase_offset_fraction_of_pi']) for row in contrast_rows],
        contrast_coherence,
        color='#8f8f8f',
        label='transient mixed controls',
        s=45,
    )
    axes[1].scatter(stable_phase_values, stable_coherence, color='#c23b22', label='stable smeared cells', s=55)
    axes[1].set_xlabel('phase corridor position')
    axes[1].set_ylabel('grade-exchange coherence')
    axes[1].set_title('Secondary descriptor overlap')
    axes[1].legend(frameon=False, fontsize=8)
    corr_fig.tight_layout()
    corr_path = WORK_PLOT_DIR / 'stage24_1_smeared_sector_correlation_panel.png'
    corr_fig.savefig(corr_path, dpi=200, bbox_inches='tight')
    plt.close(corr_fig)

    refinement_grid, weak_grid, reverse_grid, width_values, phase_values = build_stability_grids(s12_rows)
    stab_fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.2), sharey=True)
    for ax, grid, title in zip(
        axes,
        [refinement_grid, weak_grid, reverse_grid],
        ['Refinement stability', 'Weak-coupling stability', 'Reverse stability'],
    ):
        image = ax.imshow(grid, vmin=0.0, vmax=1.0, cmap='Greys', aspect='auto')
        ax.set_title(title)
        ax.set_xticks(range(len(phase_values)))
        ax.set_xticklabels([f'{value:.3f}' for value in phase_values], rotation=45, ha='right', fontsize=8)
        ax.set_xlabel('phase / pi')
        ax.set_yticks(range(len(width_values)))
        ax.set_yticklabels([f'{value:.2f}' for value in width_values], fontsize=8)
        ax.set_ylabel('width ratio')
        for i in range(len(width_values)):
            for j in range(len(phase_values)):
                ax.text(j, i, f'{int(grid[i, j])}', ha='center', va='center', color='#c23b22' if grid[i, j] == 1 else '#111111', fontsize=8)
    stab_fig.suptitle('Stage 24.1 closure-support stability panel', fontsize=12)
    stab_fig.tight_layout()
    stab_path = WORK_PLOT_DIR / 'stage24_1_smeared_sector_stability_panel.png'
    stab_fig.savefig(stab_path, dpi=200, bbox_inches='tight')
    plt.close(stab_fig)

    json_path, csv_path, stamped_plots, _timestamp = save_atlas_payload(
        experiment_slug='stage24_1_smeared_sector_observable_ledger',
        result=summary,
        csv_rows=ledger_rows,
        csv_fieldnames=CSV_FIELDS,
        plot_paths=[table_path, corr_path, stab_path],
    )

    artifact_paths = {
        'json': str(json_path.relative_to(REPO_ROOT)),
        'csv': str(csv_path.relative_to(REPO_ROOT)),
        'plots': stamped_plots,
    }
    write_note(
        note_path=NOTE_PATH,
        runsheet=runsheet,
        summary=summary,
        ledger_rows=ledger_rows,
        artifact_paths=artifact_paths,
    )

    print(f'Wrote note: {NOTE_PATH.relative_to(REPO_ROOT)}')
    print(f'Wrote JSON: {json_path.relative_to(REPO_ROOT)}')
    print(f'Wrote CSV: {csv_path.relative_to(REPO_ROOT)}')
    for plot in stamped_plots:
        print(f'Wrote plot: {plot}')


if __name__ == '__main__':
    main()
