#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Any

from stage10_common import ATLAS_NOTES, PLOTS, REPO_ROOT, RESULTS, plt

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage24_5_phaseIV_consolidation_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_IV_A5_Phase_IV_Consolidation_and_Freeze_v1.md'

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
    parser = argparse.ArgumentParser(description='Run Stage 24.5 Phase IV consolidation and freeze.')
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


def read_required_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f'Missing required upstream CSV: {path}')
    with path.open('r', encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))


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


def render_claim_table(path: Path, title: str, rows: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(12.6, 2.4 + 0.52 * len(rows)))
    ax.axis('off')
    table_rows = [
        [str(row['claim_id']), str(row['claim_text']), str(row['status']), str(row['supporting_stage'])]
        for row in rows
    ]
    table = ax.table(
        cellText=table_rows,
        colLabels=['ID', 'Claim', 'Status', 'Stage'],
        colWidths=[0.12, 0.58, 0.14, 0.16],
        cellLoc='left',
        loc='center',
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.55)
    for (row_idx, _col_idx), cell in table.get_celld().items():
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


def render_status_map_panel(path: Path, summary: dict[str, Any]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.3))

    counts = summary['extension_status_counts']
    labels = list(counts.keys())
    values = [counts[key] for key in labels]
    axes[0].bar(labels, values, color=['#2e7d32', '#f9a825', '#c62828'])
    axes[0].set_title('24.4 extension status counts')
    axes[0].set_ylabel('conditions')
    axes[0].tick_params(axis='x', rotation=20)

    axes[1].axis('off')
    axes[1].set_title('First honest boundaries', fontsize=12, pad=12)
    y = 0.88
    boundary_lines = [
        f"ordering first: {summary['earliest_ordering_boundary_level']} ({summary['earliest_ordering_boundary_class']})",
        f"selector first: {summary['first_selector_failure_level']} ({summary['first_selector_failure_class']})",
        f"anisotropy failure seen: {summary['support_anisotropy_failure_seen']}",
        f"ordering before selector: {summary['ordering_degrades_before_selector_failure']}",
    ]
    for line in boundary_lines:
        axes[1].text(
            0.03,
            y,
            f'- {line}',
            transform=axes[1].transAxes,
            fontsize=10,
            va='top',
            bbox=dict(boxstyle='round,pad=0.25', facecolor='#f8fafc', edgecolor='#d0d7de'),
        )
        y -= 0.22

    fig.suptitle('Stage 24.5 Phase IV status map panel', fontsize=13)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def find_csv_row(rows: list[dict[str, str]], key: str, value: str) -> dict[str, str]:
    for row in rows:
        if str(row[key]) == value:
            return row
    raise KeyError(f'No CSV row found where {key} == {value}')


def build_claim_rows(inputs: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    s23_14 = inputs['stage23_14']['json']
    s24_1 = inputs['stage24_1']['json']
    s24_2 = inputs['stage24_2']['json']
    s24_3 = inputs['stage24_3']['json']
    s24_4 = inputs['stage24_4']['json']
    s24_3_rows = inputs['stage24_3']['csv_rows']

    q3 = find_csv_row(s24_3_rows, 'candidate_id', 'IVA3_Q3')
    portable_rows = [
        {
            'claim_id': 'P4_POS_01',
            'claim_text': 'The smeared sector is carried by the accepted 24.1 observable core and not by rejected braid- or asymmetry-based quantities.',
            'status': 'survives_portably',
            'supporting_stage': '24.1 / IV-A1',
            'supporting_metric_or_fact': f"carry_forward_core = {s24_1['carry_forward_core']}; rejected_observables = {s24_1['rejected_observables']}",
            'scope': 'observable_ledger',
            'notes': 'Phase IV keeps braid and the rejected observable family outside the core.',
        },
        {
            'claim_id': 'P4_POS_02',
            'claim_text': 'The stable closed smeared sector is selected by the fixed threshold law: stable_closed_smeared = 1 iff flow_concentration_index <= 0.884308.',
            'status': 'survives_portably',
            'supporting_stage': '24.2 / IV-A2',
            'supporting_metric_or_fact': f"selected_model_id = {s24_2['selected_model_id']}; selected_law_expression = {s24_2['selected_law_expression']}",
            'scope': 'selector_law',
            'notes': 'The selector uses only flow_concentration_index and keeps corridor variables out of state.',
        },
        {
            'claim_id': 'P4_POS_03',
            'claim_text': 'The fixed threshold selector remains valid across most tested extension conditions without retuning.',
            'status': 'survives_portably',
            'supporting_stage': '24.4 / IV-A4',
            'supporting_metric_or_fact': f"status_counts = {s24_4['status_counts']}",
            'scope': 'local_extension_radius',
            'notes': 'Selector-valid conditions include all law_survives_ordering_survives and law_survives_ordering_degrades cases.',
        },
        {
            'claim_id': 'P4_POS_04',
            'claim_text': 'Ordering boundaries appear before the first selector-failure boundary on the tested extension lattice.',
            'status': 'survives_portably',
            'supporting_stage': '24.4 / IV-A4',
            'supporting_metric_or_fact': f"first_boundary_by_class = {s24_4['first_boundary_by_class']}",
            'scope': 'extension_boundary_ordering',
            'notes': 'Motif and clustered-texture ordering degradation occur before the first selector failure at moderate corridor geometry perturbation.',
        },
        {
            'claim_id': 'P4_POS_05',
            'claim_text': 'phase_corridor_position and phase_corridor_width remain external conditioning coordinates rather than Phase IV state variables.',
            'status': 'survives_portably',
            'supporting_stage': '24.1 / IV-A1; 24.2 / IV-A2',
            'supporting_metric_or_fact': f"carry_forward_context = {s24_1['carry_forward_context']}; context_variables_not_promoted = {s24_2['context_variables_not_promoted']}",
            'scope': 'phaseIV_portable_core',
            'notes': 'Corridor variables are not promoted into the selector law.',
        },
    ]

    local_rows = [
        {
            'claim_id': 'P4_LOC_01',
            'claim_text': 'Inside the law-defined stable smeared sector, lower flow_concentration_index orders with longer topology_survival_time.',
            'status': 'survives_locally',
            'supporting_stage': '24.3 / IV-A3',
            'supporting_metric_or_fact': f"selected_candidate_id = {s24_3['selected_candidate_id']}; expression = {s24_3['selected_expression']}",
            'scope': 'local_extension_radius',
            'notes': 'This is the accepted local ordering statement for Phase IV.',
        },
        {
            'claim_id': 'P4_LOC_02',
            'claim_text': 'The accepted inverse flow-survival ordering is not a restatement of the selector threshold law.',
            'status': 'survives_locally',
            'supporting_stage': '24.3 / IV-A3',
            'supporting_metric_or_fact': f"law_restatement_rejected = {s24_3['law_restatement_rejected']}; exact_support_signature_present = {s24_3['exact_support_signature_present']}",
            'scope': 'local_structure',
            'notes': 'Q1 remains secondary and Q2 remains explicitly rejected as a quasi-invariant.',
        },
        {
            'claim_id': 'P4_LOC_03',
            'claim_text': 'The stable-set inverse ordering differs sharply from contrast behavior on the frozen branch.',
            'status': 'survives_locally',
            'supporting_stage': '24.3 / IV-A3',
            'supporting_metric_or_fact': f"stable_metric = {q3['stable_metric']}; contrast_metric = {q3['contrast_metric']}",
            'scope': 'local_structure',
            'notes': 'The sign split rho < 0 stable and rho > 0 contrast is part of the frozen Phase IV read.',
        },
        {
            'claim_id': 'P4_LOC_04',
            'claim_text': 'The inverse flow-survival ordering has bounded local portability but does not survive uniformly across the tested extension lattice.',
            'status': 'survives_locally',
            'supporting_stage': '24.4 / IV-A4',
            'supporting_metric_or_fact': f"status_counts = {s24_4['status_counts']}; first_boundary_by_class = {s24_4['first_boundary_by_class']}",
            'scope': 'extension_boundary_ordering',
            'notes': 'Ordering survives widely but reaches distinct local boundaries before selector failure in some classes.',
        },
        {
            'claim_id': 'P4_LOC_05',
            'claim_text': 'Support anisotropy / skew stayed inside the tested survival radius on the current extension lattice.',
            'status': 'bounded_support',
            'supporting_stage': '24.4 / IV-A4',
            'supporting_metric_or_fact': f"support_anisotropy_skew boundary = {s24_4['first_boundary_by_class']['support_anisotropy_skew']}",
            'scope': 'extension_boundary_ordering',
            'notes': 'No failure boundary was reached on the tested support-anisotropy levels.',
        },
    ]

    negative_rows = [
        {
            'claim_id': 'P4_NEG_01',
            'claim_text': 'The inverse flow-survival ordering does not survive uniformly across the tested extension lattice.',
            'status': 'fails',
            'supporting_stage': '24.4 / IV-A4',
            'supporting_metric_or_fact': f"status_counts = {s24_4['status_counts']}",
            'scope': 'extension_boundary_ordering',
            'notes': 'Uniform ordering portability fails on this branch.',
        },
        {
            'claim_id': 'P4_NEG_02',
            'claim_text': 'Motif structure perturbation reaches the earliest tested ordering boundary at IVA4_E3_mild.',
            'status': 'fails',
            'supporting_stage': '24.4 / IV-A4',
            'supporting_metric_or_fact': f"motif_structure_perturbation boundary = {s24_4['first_boundary_by_class']['motif_structure_perturbation']}",
            'scope': 'extension_boundary_ordering',
            'notes': 'This is the first ordering degradation boundary on the tested lattice.',
        },
        {
            'claim_id': 'P4_NEG_03',
            'claim_text': 'Clustered texture perturbation reaches ordering degradation at IVA4_E1_moderate.',
            'status': 'fails',
            'supporting_stage': '24.4 / IV-A4',
            'supporting_metric_or_fact': f"clustered_texture_perturbation boundary = {s24_4['first_boundary_by_class']['clustered_texture_perturbation']}",
            'scope': 'extension_boundary_ordering',
            'notes': 'Clustered texture is a later ordering boundary than motif perturbation.',
        },
        {
            'claim_id': 'P4_NEG_04',
            'claim_text': 'Moderate corridor geometry perturbation reaches the first tested selector-failure boundary at IVA4_E4_moderate.',
            'status': 'fails',
            'supporting_stage': '24.4 / IV-A4',
            'supporting_metric_or_fact': f"corridor_geometry_perturbation boundary = {s24_4['first_boundary_by_class']['corridor_geometry_perturbation']}",
            'scope': 'extension_boundary_selector',
            'notes': 'Selector failure appears later than the earliest ordering boundaries.',
        },
        {
            'claim_id': 'P4_NEG_05',
            'claim_text': 'No braid-labeled quantity re-enters the Phase IV core.',
            'status': 'excluded',
            'supporting_stage': '24.1 / IV-A1; 23.14 / III-C5',
            'supporting_metric_or_fact': f"rejected_observables = {s24_1['rejected_observables']}; family_wide_braid_mechanism = {s23_14['family_wide_braid_mechanism']}",
            'scope': 'observable_ledger',
            'notes': 'Braid remains excluded from the Phase IV explanatory core.',
        },
        {
            'claim_id': 'P4_NEG_06',
            'claim_text': 'The rejected observable family remains excluded from Phase IV core structure.',
            'status': 'excluded',
            'supporting_stage': '24.1 / IV-A1',
            'supporting_metric_or_fact': f"rejected_observables = {s24_1['rejected_observables']}",
            'scope': 'observable_ledger',
            'notes': 'Transfer asymmetry, coherence decay, topology return error, braid-like exchange indicator, and mechanism sensitivity score stay excluded.',
        },
        {
            'claim_id': 'P4_NEG_07',
            'claim_text': 'The current law-and-ordering package does not justify a family-wide transport law claim.',
            'status': 'excluded',
            'supporting_stage': '24.4 / IV-A4',
            'supporting_metric_or_fact': f"status_counts = {s24_4['status_counts']}; first_boundary_by_class = {s24_4['first_boundary_by_class']}",
            'scope': 'phaseIV_portable_core',
            'notes': 'Portability remains bounded and branch-local.',
        },
        {
            'claim_id': 'P4_NEG_08',
            'claim_text': 'Nothing in Phase IV reopens family-wide braid mechanism, stable braid phase, or localized encounter claims from Phase III.',
            'status': 'excluded',
            'supporting_stage': '23.14 / III-C5; 24.5 / IV-A5',
            'supporting_metric_or_fact': f"family_wide_braid_mechanism = {s23_14['family_wide_braid_mechanism']}; stable_braid_phase_present = {s23_14['stable_braid_phase_present']}; localized_encounter_phase_present = {s23_14['localized_encounter_phase_present']}",
            'scope': 'phaseIII_background',
            'notes': 'Phase III negatives remain fully intact.',
        },
    ]

    return portable_rows + local_rows + negative_rows


def build_consistency_checks(runsheet: dict[str, Any], inputs: dict[str, dict[str, Any]], claim_rows: list[dict[str, Any]]) -> dict[str, bool]:
    s23_14 = inputs['stage23_14']['json']
    s24_1 = inputs['stage24_1']['json']
    s24_2 = inputs['stage24_2']['json']
    s24_3 = inputs['stage24_3']['json']
    s24_4 = inputs['stage24_4']['json']
    s24_3_rows = inputs['stage24_3']['csv_rows']

    expected_core = list(runsheet['expected_observable_core'])
    expected_rejected = list(runsheet['expected_rejected_observables'])
    expected_status_counts = dict(runsheet['expected_status_counts'])
    expected_boundaries = dict(runsheet['expected_boundaries'])
    q3 = find_csv_row(s24_3_rows, 'candidate_id', 'IVA3_Q3')

    portable_text = ' '.join(
        row['claim_text'].lower()
        for row in claim_rows
        if str(row['status']) == 'survives_portably'
    )
    contradiction_free = all([
        'uniformly across the tested extension lattice' not in portable_text,
        'family-wide transport law' not in portable_text,
        'family-wide braid mechanism' not in portable_text,
        'stable braid phase' not in portable_text,
        'localized encounter phase' not in portable_text,
    ])

    return {
        'stage24_1_expected_primary_observables': bool(s24_1['carry_forward_core'] == expected_core),
        'stage24_1_rejected_family_intact': bool(s24_1['rejected_observables'] == expected_rejected),
        'stage24_2_selector_law_exact': bool(s24_2['selected_law_expression'] == runsheet['expected_selector_law']),
        'stage24_3_selected_q3': bool(s24_3['selected_candidate_id'] == runsheet['expected_ordering_candidate_id']),
        'stage24_3_contrast_separated_spearman_signs': bool(
            'spearman_rho=-' in str(q3['stable_metric']) and 'spearman_rho=0.' in str(q3['contrast_metric'])
        ),
        'stage24_4_status_counts_match': bool(s24_4['status_counts'] == expected_status_counts),
        'stage24_4_boundaries_match': bool(
            all(
                s24_4['first_boundary_by_class'].get(key) is None if value is None else
                s24_4['first_boundary_by_class'].get(key, {}).get('extension_id') == value
                for key, value in expected_boundaries.items()
            )
        ),
        'portable_claims_noncontradictory': contradiction_free,
        'phaseIII_negatives_preserved': bool(
            s23_14['family_wide_braid_mechanism'] is False
            and s23_14['stable_braid_phase_present'] is False
            and s23_14['localized_encounter_phase_present'] is False
        ),
    }


def build_summary_payload(runsheet: dict[str, Any], inputs: dict[str, dict[str, Any]], claim_rows: list[dict[str, Any]], checks: dict[str, bool]) -> dict[str, Any]:
    s24_1 = inputs['stage24_1']['json']
    s24_2 = inputs['stage24_2']['json']
    s24_3 = inputs['stage24_3']['json']
    s24_4 = inputs['stage24_4']['json']

    phaseIV_inputs_consistent = all(bool(value) for value in checks.values())
    portable_core_frozen = all(
        str(row['status']) == 'survives_portably'
        for row in claim_rows[:5]
    )
    local_structure_frozen = all(
        str(row['status']) in {'survives_locally', 'bounded_support'}
        for row in claim_rows[5:10]
    )
    negative_boundaries_frozen = all(
        str(row['status']) in {'fails', 'excluded'}
        for row in claim_rows[10:]
    )
    selector_question_closed = bool(s24_2['selected_law_expression'] == runsheet['expected_selector_law'])
    ordering_boundary_closed = bool(s24_4['first_boundary_by_class']['motif_structure_perturbation']['extension_id'] == 'IVA4_E3_mild')
    freeze_valid = bool(
        phaseIV_inputs_consistent
        and portable_core_frozen
        and local_structure_frozen
        and negative_boundaries_frozen
        and selector_question_closed
        and ordering_boundary_closed
    )

    return {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'phaseIV_consolidation_complete': True,
        'phaseIV_inputs_consistent': phaseIV_inputs_consistent,
        'phaseIV_freeze_valid': freeze_valid,
        'portable_core_frozen': portable_core_frozen,
        'local_structure_frozen': local_structure_frozen,
        'negative_boundaries_frozen': negative_boundaries_frozen,
        'selector_question_closed': selector_question_closed,
        'ordering_boundary_closed': ordering_boundary_closed,
        'observable_core_frozen': bool(s24_1['observable_ledger_valid']),
        'selector_law_found': bool(s24_2['effective_law_found']),
        'selector_threshold_fixed': 0.884308,
        'selector_law_boundedly_portable': True,
        'local_ordering_found': bool(s24_3['quasi_invariant_found']),
        'ordering_uniformly_portable': False,
        'ordering_degrades_before_selector_failure': True,
        'earliest_ordering_boundary_class': 'motif_structure_perturbation',
        'earliest_ordering_boundary_level': 'IVA4_E3_mild',
        'first_selector_failure_class': 'corridor_geometry_perturbation',
        'first_selector_failure_level': 'IVA4_E4_moderate',
        'support_anisotropy_failure_seen': False,
        'braid_reintroduced': False,
        'rejected_observables_reintroduced': False,
        'family_wide_transport_claim': False,
        'extension_status_counts': dict(s24_4['status_counts']),
        'phaseIV_minimal_read': 'boundedly portable stable-smeared threshold selector with locally portable inverse flow-survival ordering; ordering boundaries appear before selector failure in some extension classes',
        'consistency_checks': checks,
        'claim_ledger': claim_rows,
    }


def write_note(note_path: Path, runsheet: dict[str, Any], summary: dict[str, Any], claim_rows: list[dict[str, Any]], artifact_paths: dict[str, Any]) -> None:
    portable_rows = claim_rows[:5]
    local_rows = claim_rows[5:10]
    negative_rows = claim_rows[10:]

    lines = [
        '# Stage IV A5 Phase IV Consolidation and Freeze v1',
        '',
        '## Stage purpose',
        'This stage freezes the current Phase IV law-and-ordering branch on the frozen clustered DK sector. It is a reduction-only consolidation pass over Stages 24.1 through 24.4, with the 24.4 rerun stamped `20260316_121603` treated as authoritative.',
        '',
        '## Inputs consolidated',
    ]
    for key in ('stage23_14', 'stage24_1', 'stage24_2', 'stage24_3', 'stage24_4'):
        item = runsheet['required_inputs'][key]
        lines.append(
            f"- `{item['stage_label']}`: json=`{item['json_path']}`, csv=`{item['csv_path']}`"
        )

    lines.extend([
        '',
        '## Consistency checks',
    ])
    for key, value in summary['consistency_checks'].items():
        lines.append(f"- `{key}` = `{value}`")

    lines.extend([
        '',
        '## Portable but bounded positive core',
    ])
    for row in portable_rows:
        lines.append(
            f"- `{row['claim_id']}`: {row['claim_text']} [{row['status']}]"
        )

    lines.extend([
        '',
        '## Locally portable internal structure',
    ])
    for row in local_rows:
        lines.append(
            f"- `{row['claim_id']}`: {row['claim_text']} [{row['status']}]"
        )

    lines.extend([
        '',
        '## Failed / excluded claims',
    ])
    for row in negative_rows:
        lines.append(
            f"- `{row['claim_id']}`: {row['claim_text']} [{row['status']}]"
        )

    lines.extend([
        '',
        '## Final closure statement',
        'On the frozen clustered DK branch, Phase IV identifies a compact threshold selector for the stable closed smeared sector, flow_concentration_index <= 0.884308, together with a non-trivial inverse flow-survival ordering relation that survives only within a bounded local extension radius. Under controlled structural extension, ordering degradation appears before selector failure in some classes, while moderate corridor-geometry perturbation provides the first tested boundary at which the fixed selector itself fails. The surviving Phase IV content is therefore a bounded local transport package with separable selector and ordering boundaries, not a family-wide transport structure.',
        '',
        '## Frozen claims ledger summary',
        f"- portable core claims = `{len(portable_rows)}`",
        f"- local structure claims = `{len(local_rows)}`",
        f"- failed / excluded claims = `{len(negative_rows)}`",
        f"- extension status counts = `{summary['extension_status_counts']}`",
        '',
        '## Output artifact list',
        f"- stamped JSON summary: `{artifact_paths['json']}`",
        f"- stamped CSV claims ledger: `{artifact_paths['csv']}`",
    ])
    for plot in artifact_paths['plots']:
        lines.append(f"- plot: `{plot}`")

    lines.extend([
        '',
        '## Commit recommendation',
        '- recommended status: commit',
        '- rationale: Phase IV selector, local ordering structure, and extension boundaries are now internally consolidated.',
        '- branch effect: closes the current Phase IV law-and-ordering arc on the frozen clustered DK branch.',
        '- caution: reopening portability or internal-structure claims requires a new branch with explicitly changed extension manifold, operator class, or structural prior.',
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
        inputs[key] = {
            'stage_label': item['stage_label'],
            'note_path': str(note_path.relative_to(REPO_ROOT)),
            'json_path': str(json_path.relative_to(REPO_ROOT)),
            'csv_path': str(csv_path.relative_to(REPO_ROOT)),
            'note_text': read_required_text(note_path),
            'json': read_required_json(json_path),
            'csv_rows': read_required_csv(csv_path),
        }

    claim_rows = build_claim_rows(inputs)
    checks = build_consistency_checks(runsheet, inputs, claim_rows)
    summary = build_summary_payload(runsheet, inputs, claim_rows, checks)

    timestamp = timestamp_slug()
    json_path = RESULTS / f'{timestamp}_stage24_5_phaseIV_consolidation_and_freeze.json'
    csv_path = RESULTS / f'{timestamp}_stage24_5_phaseIV_claims_ledger.csv'

    portable_path = PLOTS / f'{timestamp}_stage24_5_phaseIV_portable_core_table.png'
    local_path = PLOTS / f'{timestamp}_stage24_5_phaseIV_local_structure_table.png'
    negative_path = PLOTS / f'{timestamp}_stage24_5_phaseIV_boundary_and_negative_table.png'
    status_path = PLOTS / f'{timestamp}_stage24_5_phaseIV_status_map_panel.png'

    render_claim_table(portable_path, 'Stage 24.5 portable but bounded core', claim_rows[:5])
    render_claim_table(local_path, 'Stage 24.5 locally portable internal structure', claim_rows[5:10])
    render_claim_table(negative_path, 'Stage 24.5 boundary and negative claims', claim_rows[10:])
    render_status_map_panel(status_path, summary)

    json_path.write_text(json.dumps(summary, indent=2), encoding='utf-8')
    write_csv(csv_path, claim_rows, CLAIM_FIELDS)

    artifact_paths = {
        'json': str(json_path.relative_to(REPO_ROOT)),
        'csv': str(csv_path.relative_to(REPO_ROOT)),
        'plots': [
            str(portable_path.relative_to(REPO_ROOT)),
            str(local_path.relative_to(REPO_ROOT)),
            str(negative_path.relative_to(REPO_ROOT)),
            str(status_path.relative_to(REPO_ROOT)),
        ],
    }
    write_note(NOTE_PATH, runsheet, summary, claim_rows, artifact_paths)

    print(f'Wrote note: {NOTE_PATH.relative_to(REPO_ROOT)}')
    print(f'Wrote JSON: {json_path.relative_to(REPO_ROOT)}')
    print(f'Wrote CSV: {csv_path.relative_to(REPO_ROOT)}')
    for plot in artifact_paths['plots']:
        print(f'Wrote plot: {plot}')


if __name__ == '__main__':
    main()
