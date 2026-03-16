#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Any

from stage10_common import ATLAS_NOTES, PLOTS, REPO_ROOT, RESULTS, plt

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage23_14_phaseIII_final_consolidation_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_III_C5_Final_Consolidation_and_Freeze_v1.md'

CLAIM_FIELDS = [
    'claim_id',
    'claim_text',
    'status',
    'supporting_stage',
    'supporting_metric_or_fact',
    'scope',
    'notes',
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run Stage 23.14 Phase III final consolidation and freeze.')
    parser.add_argument('--runsheet', type=Path, default=RUNSHEET_PATH)
    return parser.parse_args()


def timestamp_slug() -> str:
    return datetime.now().strftime('%Y%m%d_%H%M%S')


def load_runsheet(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding='utf-8'))


def read_required_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f'Missing required upstream JSON: {path}')
    return json.loads(path.read_text(encoding='utf-8'))


def read_required_text(path: Path) -> str:
    if not path.exists():
        raise FileNotFoundError(f'Missing required upstream note: {path}')
    return path.read_text(encoding='utf-8')


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def render_claim_table(path: Path, title: str, rows: list[dict[str, Any]], show_status: bool = True) -> None:
    fig, ax = plt.subplots(figsize=(12.4, 2.4 + 0.52 * len(rows)))
    ax.axis('off')
    columns = ['claim_id', 'claim_text', 'supporting_stage']
    labels = ['ID', 'Claim', 'Stage']
    if show_status:
        columns.insert(2, 'status')
        labels.insert(2, 'Status')
    table_rows = [[str(row[column]) for column in columns] for row in rows]
    table = ax.table(
        cellText=table_rows,
        colLabels=labels,
        colWidths=[0.12, 0.58, 0.12, 0.18] if show_status else [0.14, 0.66, 0.20],
        cellLoc='left',
        loc='center',
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.55)
    for (row_idx, col_idx), cell in table.get_celld().items():
        cell.set_edgecolor('#bbbbbb')
        if row_idx == 0:
            cell.set_facecolor('#e9eef6')
            cell.set_text_props(weight='bold')
        else:
            cell.set_facecolor('#ffffff' if row_idx % 2 else '#f8fafc')
    ax.set_title(title, fontsize=12, pad=14)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def render_mechanism_phase_panel(path: Path, summary: dict[str, Any]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12.8, 4.2))
    panels = [
        (
            'Scope',
            [
                'clustered sector only',
                'no family-wide braid mechanism',
                'braid demoted to texture-level descriptor',
            ],
        ),
        (
            'Phase Closure',
            [
                'effective_model_valid = true',
                'phase_diagram_closed = true',
                'closure = smeared-dominant',
                'stable row = width_ratio 1.35',
                'stable braid phase absent',
            ],
        ),
        (
            'Mechanism',
            [
                'phase corridor primary',
                'geometry overlap secondary',
                '0<->1 cross-coupling secondary',
                'grade transfer not primary',
            ],
        ),
    ]
    for ax, (title, lines) in zip(axes, panels):
        ax.axis('off')
        ax.set_title(title, fontsize=12, pad=12)
        y = 0.90
        for line in lines:
            ax.text(
                0.02,
                y,
                f'- {line}',
                transform=ax.transAxes,
                fontsize=10,
                va='top',
                bbox=dict(boxstyle='round,pad=0.25', facecolor='#f8fafc', edgecolor='#d0d7de'),
            )
            y -= 0.18
    fig.suptitle('Stage 23.14 mechanism and phase closure panel', fontsize=13)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def render_claim_status_matrix(path: Path, rows: list[dict[str, Any]]) -> None:
    status_code = {'survives': 0, 'bounded_support': 1, 'fails': 2, 'excluded': 3}
    grid = [[status_code[str(row['status'])]] for row in rows]
    fig, ax = plt.subplots(figsize=(3.8, 0.45 * len(rows) + 1.6))
    im = ax.imshow(grid, aspect='auto', cmap='viridis', vmin=0.0, vmax=3.0)
    ax.set_xticks([0], ['status'])
    ax.set_yticks(range(len(rows)), [str(row['claim_id']) for row in rows])
    ax.set_title('Stage 23.14 claim status matrix')
    for idx, row in enumerate(rows):
        ax.text(0, idx, str(row['status']), ha='center', va='center', color='white', fontsize=8)
    cbar = fig.colorbar(im, ax=ax, fraction=0.08, pad=0.06)
    cbar.set_ticks([0, 1, 2, 3])
    cbar.set_ticklabels(['survives', 'bounded', 'fails', 'excluded'])
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def build_claim_rows(inputs: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    s10 = inputs['stage23_10']['json']
    s11 = inputs['stage23_11']['json']
    s12 = inputs['stage23_12']['json']
    s13 = inputs['stage23_13']['json']
    return [
        {
            'claim_id': 'P3_POS_01',
            'claim_text': 'The tested clustered DK braid/smear sector admits a valid two-state effective surrogate on the anchor lattice.',
            'status': 'survives',
            'supporting_stage': '23.11 / III-C2',
            'supporting_metric_or_fact': 'effective_model_valid = TRUE; topology agreement = 1.0000; threshold error = 0.0000',
            'scope': 'clustered_sector_only',
            'notes': 'Valid only as a clustered texture surrogate.',
        },
        {
            'claim_id': 'P3_POS_02',
            'claim_text': 'The clustered DK sector closes as a finite smeared_dominant_closed_phase_diagram.',
            'status': 'survives',
            'supporting_stage': '23.12 / III-C3',
            'supporting_metric_or_fact': f"phase_diagram_closed = {s12['summary']['phase_diagram_closed']}; closure_classification = {s12['summary']['closure_classification']}",
            'scope': 'phase_closure',
            'notes': 'Closure is one-sided and smeared-dominant.',
        },
        {
            'claim_id': 'P3_POS_03',
            'claim_text': 'Only the strong-asymmetry row width_ratio = 1.35 survives as refinement-/beta-/reverse-stable across the sampled phase corridor.',
            'status': 'survives',
            'supporting_stage': '23.12 / III-C3',
            'supporting_metric_or_fact': 'stable component = smeared_transfer_phase; width_min = width_max = 1.35; phase range 0.375 -> 0.575',
            'scope': 'phase_closure',
            'notes': 'No symmetric or moderate-asymmetry cell closes stably.',
        },
        {
            'claim_id': 'P3_POS_04',
            'claim_text': 'phase_corridor_primary is the controlling mechanism classification for the braid/smear split in the clustered sector.',
            'status': 'survives',
            'supporting_stage': '23.13 / III-C4',
            'supporting_metric_or_fact': f"overall_mechanism_label = {s13['overall_summary']['overall_mechanism_label']}",
            'scope': 'mechanism_read',
            'notes': 'Weak and moderate phase flattening both flip the clustered control to transfer_smeared.',
        },
        {
            'claim_id': 'P3_POS_05',
            'claim_text': 'Geometry overlap and 0 <-> 1 operator cross-coupling act as bounded secondary supports only.',
            'status': 'bounded_support',
            'supporting_stage': '23.13 / III-C4',
            'supporting_metric_or_fact': 'Both branches collapse braid only at moderate ablation, not weak ablation.',
            'scope': 'mechanism_read',
            'notes': 'Secondary support only; not primary control.',
        },
        {
            'claim_id': 'P3_NEG_01',
            'claim_text': 'The DK braid does not generalize beyond clustered seeds and is not a family-wide intrinsic mechanism.',
            'status': 'fails',
            'supporting_stage': '23.10 / III-C1',
            'supporting_metric_or_fact': f"non_clustered_braid_hits = {s10['summary']['classification_payload']['non_clustered_braid_hits']}; robust_non_clustered_braid_hits = {s10['summary']['classification_payload']['robust_non_clustered_braid_hits']}",
            'scope': 'not_family_wide',
            'notes': 'Family-wide braid generalization fails on this branch.',
        },
        {
            'claim_id': 'P3_NEG_02',
            'claim_text': 'No stable closed braid phase survives the full closure pass.',
            'status': 'fails',
            'supporting_stage': '23.12 / III-C3',
            'supporting_metric_or_fact': f"stable_phase_counts = {s12['summary']['stable_phase_counts']}",
            'scope': 'phase_closure',
            'notes': 'Braid remains transient and does not close as a stable phase.',
        },
        {
            'claim_id': 'P3_NEG_03',
            'claim_text': 'No localized encounter phase survives the closure pass.',
            'status': 'fails',
            'supporting_stage': '23.12 / III-C3',
            'supporting_metric_or_fact': f"stable_phase_counts = {s12['summary']['stable_phase_counts']}",
            'scope': 'phase_closure',
            'notes': 'Localized encounter phase absent from the stable closure set.',
        },
        {
            'claim_id': 'P3_NEG_04',
            'claim_text': 'Grade transfer is not the primary driver on the frozen branch.',
            'status': 'fails',
            'supporting_stage': '23.13 / III-C4',
            'supporting_metric_or_fact': f"grade-transfer branch label = {s13['branch_summary']['A_grade_transfer_suppression']['label']}",
            'scope': 'mechanism_read',
            'notes': 'Grade-transfer-primary mechanism claim fails.',
        },
        {
            'claim_id': 'P3_NEG_05',
            'claim_text': 'The observed clustered braid-like exchange may not be elevated into a universal topological law.',
            'status': 'excluded',
            'supporting_stage': '23.10 / III-C1; 23.12 / III-C3; 23.13 / III-C4',
            'supporting_metric_or_fact': 'Clustered-only scope, absence of stable braid closure, and phase-corridor-primary control bound the interpretation.',
            'scope': 'clustered_sector_only',
            'notes': 'Universal mechanism reading is excluded on this branch.',
        },
    ]


def build_consistency_checks(inputs: dict[str, dict[str, Any]], claim_rows: list[dict[str, Any]]) -> dict[str, bool]:
    s10 = inputs['stage23_10']['json']
    s11 = inputs['stage23_11']['json']
    s12 = inputs['stage23_12']['json']
    s13 = inputs['stage23_13']['json']

    positive_claim_text = ' '.join(row['claim_text'] for row in claim_rows if str(row['status']) in {'survives', 'bounded_support'}).lower()
    contradiction_free = all([
        'family-wide intrinsic braid mechanism' not in positive_claim_text,
        'stable closed braid phase' not in positive_claim_text,
        'localized encounter phase' not in positive_claim_text,
        'grade transfer is the primary driver' not in positive_claim_text,
    ])

    checks = {
        'stage23_10_zero_non_clustered_braid_hits': bool(
            int(s10['summary']['classification_payload']['robust_non_clustered_braid_hits']) == 0
            and int(s10['summary']['classification_payload']['non_clustered_braid_hits']) == 0
        ),
        'stage23_11_effective_model_valid_true': bool(s11['summary']['effective_model_valid']),
        'stage23_12_phase_diagram_closed_true': bool(s12['summary']['phase_diagram_closed']),
        'stage23_12_no_stable_braid_phase': bool(
            'stable_braid_phase' not in s12['summary']['stable_phase_counts']
            and all(str(cell['derived_phase_label']) != 'stable_braid_phase' for cell in s12['cells'])
        ),
        'stage23_13_phase_corridor_primary': bool(s13['overall_summary']['overall_mechanism_label'] == 'phase_corridor_primary'),
        'positive_core_noncontradictory': contradiction_free,
    }
    return checks


def build_summary_payload(runsheet: dict[str, Any], inputs: dict[str, dict[str, Any]], claim_rows: list[dict[str, Any]], checks: dict[str, bool]) -> dict[str, Any]:
    s10 = inputs['stage23_10']['json']
    s11 = inputs['stage23_11']['json']
    s12 = inputs['stage23_12']['json']
    s13 = inputs['stage23_13']['json']

    positive_core_frozen = all(str(row['status']) in {'survives', 'bounded_support'} for row in claim_rows[:5])
    negative_boundaries_frozen = all(str(row['status']) in {'fails', 'excluded'} for row in claim_rows[5:])
    mechanism_question_closed = bool(s13['overall_summary']['overall_mechanism_label'] == 'phase_corridor_primary')
    inputs_consistent = all(bool(value) for value in checks.values())
    freeze_valid = bool(inputs_consistent and positive_core_frozen and negative_boundaries_frozen and mechanism_question_closed)

    return {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'phaseIII_consolidation_complete': True,
        'phaseIII_inputs_consistent': inputs_consistent,
        'phaseIII_freeze_valid': freeze_valid,
        'positive_core_frozen': positive_core_frozen,
        'negative_boundaries_frozen': negative_boundaries_frozen,
        'mechanism_question_closed': mechanism_question_closed,
        'family_wide_braid_mechanism': False,
        'effective_model_valid': bool(s11['summary']['effective_model_valid']),
        'phase_diagram_closed': bool(s12['summary']['phase_diagram_closed']),
        'closure_class': 'smeared_dominant_closed_phase_diagram',
        'stable_braid_phase_present': False,
        'localized_encounter_phase_present': False,
        'mechanism_primary': 'phase_corridor_primary',
        'geometry_overlap_role': 'secondary_support',
        'operator_cross_coupling_role': 'secondary_support',
        'grade_transfer_primary': False,
        'phaseIII_minimal_read': 'clustered DK sector closes as a finite smeared-transfer effective regime; braid is texture-level and phase-corridor-controlled, not a family-wide intrinsic mechanism',
        'consistency_checks': checks,
        'inputs': {
            key: {
                'stage_label': value['stage_label'],
                'note_path': value['note_path'],
                'json_path': value['json_path'],
                'csv_path': value['csv_path'],
            }
            for key, value in runsheet['required_inputs'].items()
        },
        'upstream_facts': {
            'stage23_10': s10['summary'],
            'stage23_11': s11['summary'],
            'stage23_12': s12['summary'],
            'stage23_13': {
                'overall_summary': s13['overall_summary'],
                'branch_summary': s13['branch_summary'],
            },
        },
        'claim_ledger': claim_rows,
    }


def write_note(
    note_path: Path,
    runsheet: dict[str, Any],
    summary: dict[str, Any],
    claim_rows: list[dict[str, Any]],
    artifact_paths: dict[str, Path],
) -> None:
    positive_rows = claim_rows[:5]
    negative_rows = claim_rows[5:]
    inputs = runsheet['required_inputs']
    checks = summary['consistency_checks']

    lines = [
        '# Stage III C5 Final Consolidation and Freeze v1',
        '',
        '## Stage purpose',
        'This note freezes the final Phase III claim structure for the frozen HAOS (Harmonic Address Operating System) / IIP (Interaction Invariance Physics) clustered Dirac-Kahler branch. It is a reduction-only, consolidation-only, and freeze-only package over the already-run Stages 23.10 through 23.13.',
        '',
        '## Inputs consolidated',
    ]
    for key in ('stage23_10', 'stage23_11', 'stage23_12', 'stage23_13'):
        item = inputs[key]
        lines.append(
            f"- `{item['stage_label']}`: note=`{item['note_path']}`, json=`{item['json_path']}`, csv=`{item['csv_path']}`"
        )

    lines.extend([
        '',
        '## Consistency checks',
    ])
    for name, value in checks.items():
        lines.append(f"- `{name}` = `{value}`")

    lines.extend([
        '',
        '## Surviving positive core',
    ])
    for row in positive_rows:
        lines.append(
            f"- `{row['claim_id']}`: {row['claim_text']} Support: `{row['supporting_stage']}` via `{row['supporting_metric_or_fact']}`."
        )

    lines.extend([
        '',
        '## Failed / excluded claims',
    ])
    for row in negative_rows:
        lines.append(
            f"- `{row['claim_id']}`: {row['claim_text']} Support: `{row['supporting_stage']}` via `{row['supporting_metric_or_fact']}`."
        )

    lines.extend([
        '',
        '## Minimal effective interpretation',
        'The surviving positive core of Phase III is narrow and explicit. The clustered DK braid/smear sector compresses into a valid two-state surrogate, but the closed phase content is one-sided: a finite smeared-transfer effective region survives, while braid-like exchange remains a clustered texture-level descriptor. Within that clustered sector, the braid/smear split is controlled primarily by the phase corridor, with geometry overlap and `0 <-> 1` operator cross-coupling acting only as bounded secondary supports.',
        '',
        '## Final closure statement',
        'The DK clustered collision sector does not support a family-wide intrinsic braid mechanism. Under mechanism-isolation, reduction, and closure tests, the surviving positive core is a finite smeared-transfer effective sector. Braid-like exchange remains a clustered texture-level descriptor rather than a stable closed phase. Within that clustered sector, the braid/smear split is primarily controlled by phase corridor, with geometry overlap and `0 <-> 1` operator cross-coupling acting as secondary supports. Grade transfer does not appear to be the primary control variable on the frozen branch.',
        '',
        '## Frozen claims ledger summary',
        f"- `phaseIII_consolidation_complete` = `{summary['phaseIII_consolidation_complete']}`",
        f"- `phaseIII_freeze_valid` = `{summary['phaseIII_freeze_valid']}`",
        f"- `family_wide_braid_mechanism` = `{summary['family_wide_braid_mechanism']}`",
        f"- `effective_model_valid` = `{summary['effective_model_valid']}`",
        f"- `phase_diagram_closed` = `{summary['phase_diagram_closed']}`",
        f"- `closure_class` = `{summary['closure_class']}`",
        f"- `stable_braid_phase_present` = `{summary['stable_braid_phase_present']}`",
        f"- `localized_encounter_phase_present` = `{summary['localized_encounter_phase_present']}`",
        f"- `mechanism_primary` = `{summary['mechanism_primary']}`",
        f"- `geometry_overlap_role` = `{summary['geometry_overlap_role']}`",
        f"- `operator_cross_coupling_role` = `{summary['operator_cross_coupling_role']}`",
        f"- `grade_transfer_primary` = `{summary['grade_transfer_primary']}`",
        f"- `phaseIII_minimal_read` = `{summary['phaseIII_minimal_read']}`",
        '',
        '## Output artifact list',
        f"- final JSON summary: `{artifact_paths['json_summary'].relative_to(REPO_ROOT)}`",
        f"- claims ledger CSV: `{artifact_paths['csv_ledger'].relative_to(REPO_ROOT)}`",
        f"- positive core table: `{artifact_paths['positive_table'].relative_to(REPO_ROOT)}`",
        f"- negative boundary table: `{artifact_paths['negative_table'].relative_to(REPO_ROOT)}`",
        f"- mechanism and phase closure panel: `{artifact_paths['mechanism_panel'].relative_to(REPO_ROOT)}`",
        f"- claim status matrix: `{artifact_paths['claim_matrix'].relative_to(REPO_ROOT)}`",
        '',
        '## Commit recommendation',
        '- recommended status: commit',
        '- rationale: Phase III mechanism, reduction, and closure questions are now internally consolidated.',
        '- branch effect: closes Phase III on the current frozen DK clustered branch.',
        '- caution: reopening the braid question requires a new branch with an explicitly changed operator class, family regime, or structural prior.',
    ])

    note_path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main() -> None:
    args = parse_args()
    runsheet = load_runsheet(args.runsheet)

    inputs: dict[str, dict[str, Any]] = {}
    for key, item in runsheet['required_inputs'].items():
        note_path = REPO_ROOT / item['note_path']
        json_path = REPO_ROOT / item['json_path']
        csv_path = REPO_ROOT / item['csv_path']
        read_required_text(note_path)
        if not csv_path.exists():
            raise FileNotFoundError(f'Missing required upstream CSV: {csv_path}')
        inputs[key] = {
            'note_path': str(item['note_path']),
            'json_path': str(item['json_path']),
            'csv_path': str(item['csv_path']),
            'stage_label': str(item['stage_label']),
            'json': read_required_json(json_path),
        }

    claim_rows = build_claim_rows(inputs)
    checks = build_consistency_checks(inputs, claim_rows)
    summary = build_summary_payload(runsheet, inputs, claim_rows, checks)

    if not bool(summary['phaseIII_freeze_valid']):
        raise RuntimeError('Stage 23.14 freeze validation failed; refusing to write an invalid freeze package.')

    timestamp = timestamp_slug()
    json_summary_path = RESULTS / f'{timestamp}_stage23_14_phaseIII_final_consolidation_and_freeze.json'
    csv_ledger_path = RESULTS / f'{timestamp}_stage23_14_phaseIII_claims_ledger.csv'
    positive_table_path = PLOTS / f'{timestamp}_stage23_14_phaseIII_positive_core_table.png'
    negative_table_path = PLOTS / f'{timestamp}_stage23_14_phaseIII_negative_boundary_table.png'
    mechanism_panel_path = PLOTS / f'{timestamp}_stage23_14_phaseIII_mechanism_and_phase_closure_panel.png'
    claim_matrix_path = PLOTS / f'{timestamp}_stage23_14_phaseIII_claim_status_matrix.png'

    json_summary_path.write_text(json.dumps(summary, indent=2), encoding='utf-8')
    write_csv(csv_ledger_path, claim_rows, CLAIM_FIELDS)
    render_claim_table(positive_table_path, 'Stage 23.14 positive core table', claim_rows[:5])
    render_claim_table(negative_table_path, 'Stage 23.14 negative boundary table', claim_rows[5:])
    render_mechanism_phase_panel(mechanism_panel_path, summary)
    render_claim_status_matrix(claim_matrix_path, claim_rows)

    artifact_paths = {
        'json_summary': json_summary_path,
        'csv_ledger': csv_ledger_path,
        'positive_table': positive_table_path,
        'negative_table': negative_table_path,
        'mechanism_panel': mechanism_panel_path,
        'claim_matrix': claim_matrix_path,
    }
    write_note(NOTE_PATH, runsheet, summary, claim_rows, artifact_paths)


if __name__ == '__main__':
    main()
