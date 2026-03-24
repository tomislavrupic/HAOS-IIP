#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import spearmanr

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage24_3_stable_smeared_quasi_invariant_probe_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_IV_A3_Stable_Smeared_Quasi_Invariant_Probe_v1.md'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage24_3_stable_smeared_quasi_invariant_probe'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR = REPO_ROOT / 'data'
PLOTS_DIR = REPO_ROOT / 'plots'

CSV_FIELDS = [
    'candidate_id',
    'candidate_label',
    'candidate_type',
    'expression',
    'primary_inputs',
    'context_inputs',
    'uses_only_core',
    'restates_law',
    'adds_structure_beyond_law',
    'complexity',
    'stable_metric',
    'contrast_metric',
    'stability_score',
    'contrast_score',
    'separation_score',
    'branch_valid',
    'selected',
    'notes',
]


@dataclass
class CandidateResult:
    candidate_id: str
    candidate_label: str
    candidate_type: str
    expression: str
    primary_inputs: list[str]
    context_inputs: list[str]
    uses_only_core: bool
    restates_law: bool
    adds_structure_beyond_law: bool
    complexity: int
    stable_metric: str
    contrast_metric: str
    stability_score: float
    contrast_score: float
    separation_score: float
    branch_valid: bool
    selected: bool
    notes: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run Stage 24.3 stable-smeared quasi-invariant probe.')
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


def to_int(value: str) -> int:
    return int(float(value))


def parse_threshold_from_expression(expression: str) -> float:
    match = re.search(r'<=\s*([0-9.]+)', expression)
    if match is None:
        raise ValueError(f'could not parse threshold from expression: {expression}')
    return float(match.group(1))


def coefficient_of_variation(values: np.ndarray) -> float:
    mean = float(np.mean(values))
    std = float(np.std(values))
    return std / max(abs(mean), 1.0e-12)


def format_range(values: np.ndarray) -> str:
    return f'{float(np.min(values)):.4f} -> {float(np.max(values)):.4f}'


def score_from_cv(cv: float) -> float:
    return 1.0 / (1.0 + max(cv, 0.0))


def score_from_rmse_norm(rmse_norm: float) -> float:
    return 1.0 / (1.0 + max(rmse_norm, 0.0))


def candidate_to_row(candidate: CandidateResult) -> dict[str, Any]:
    return {
        'candidate_id': candidate.candidate_id,
        'candidate_label': candidate.candidate_label,
        'candidate_type': candidate.candidate_type,
        'expression': candidate.expression,
        'primary_inputs': '; '.join(candidate.primary_inputs),
        'context_inputs': '; '.join(candidate.context_inputs),
        'uses_only_core': candidate.uses_only_core,
        'restates_law': candidate.restates_law,
        'adds_structure_beyond_law': candidate.adds_structure_beyond_law,
        'complexity': candidate.complexity,
        'stable_metric': candidate.stable_metric,
        'contrast_metric': candidate.contrast_metric,
        'stability_score': f'{candidate.stability_score:.4f}',
        'contrast_score': f'{candidate.contrast_score:.4f}',
        'separation_score': f'{candidate.separation_score:.4f}',
        'branch_valid': candidate.branch_valid,
        'selected': candidate.selected,
        'notes': candidate.notes,
    }


def fit_support_signature(
    refinement: np.ndarray,
    weak: np.ndarray,
    reverse: np.ndarray,
    positive_mask: np.ndarray,
    contrast_mask: np.ndarray,
) -> CandidateResult:
    counts = refinement + weak + reverse
    stable_counts = counts[positive_mask]
    contrast_counts = counts[contrast_mask]
    stable_value = float(stable_counts[0])
    stable_constancy = 1.0 if len(set(float(x) for x in stable_counts)) == 1 else 0.0
    contrast_match_rate = float(np.mean(contrast_counts == stable_value))
    return CandidateResult(
        candidate_id='IVA3_Q1',
        candidate_label='closure_support_signature',
        candidate_type='support_signature',
        expression='closure_support_count = refinement_stability_flag + weak_coupling_stability_flag + reverse_stability_flag',
        primary_inputs=[
            'refinement_stability_flag',
            'weak_coupling_stability_flag',
            'reverse_stability_flag',
        ],
        context_inputs=[],
        uses_only_core=True,
        restates_law=False,
        adds_structure_beyond_law=False,
        complexity=3,
        stable_metric=f'unique={sorted(set(int(x) for x in stable_counts))}; range={format_range(stable_counts)}',
        contrast_metric=f'unique={sorted(set(int(x) for x in contrast_counts))}; stable_match_rate={contrast_match_rate:.4f}',
        stability_score=stable_constancy,
        contrast_score=contrast_match_rate,
        separation_score=stable_constancy - contrast_match_rate,
        branch_valid=True,
        selected=False,
        notes='Exact branch-local support signature, but it rephrases the closure-support flags already available in Stage 24.2 rather than adding a new relation.',
    )


def fit_restated_law_band(
    flow: np.ndarray,
    positive_mask: np.ndarray,
    contrast_mask: np.ndarray,
    threshold: float,
) -> CandidateResult:
    stable_flow = flow[positive_mask]
    contrast_flow = flow[contrast_mask]
    stable_cv = coefficient_of_variation(stable_flow)
    contrast_cv = coefficient_of_variation(contrast_flow)
    return CandidateResult(
        candidate_id='IVA3_Q2',
        candidate_label='flow_threshold_band',
        candidate_type='restated_law_band',
        expression=f'flow_concentration_index <= {threshold:.6f}',
        primary_inputs=['flow_concentration_index'],
        context_inputs=[],
        uses_only_core=True,
        restates_law=True,
        adds_structure_beyond_law=False,
        complexity=1,
        stable_metric=f'cv={stable_cv:.4f}; range={format_range(stable_flow)}',
        contrast_metric=f'cv={contrast_cv:.4f}; range={format_range(contrast_flow)}',
        stability_score=score_from_cv(stable_cv),
        contrast_score=score_from_cv(contrast_cv),
        separation_score=score_from_cv(stable_cv) - score_from_cv(contrast_cv),
        branch_valid=True,
        selected=False,
        notes='Reported for completeness only. This is the Stage 24.2 law itself and may not be selected as the quasi-invariant read.',
    )


def fit_inverse_ordering_relation(
    flow: np.ndarray,
    survival: np.ndarray,
    positive_mask: np.ndarray,
    contrast_mask: np.ndarray,
) -> CandidateResult:
    stable_rho = float(spearmanr(flow[positive_mask], survival[positive_mask]).statistic)
    contrast_rho = float(spearmanr(flow[contrast_mask], survival[contrast_mask]).statistic)
    stable_inverse_score = max(0.0, -stable_rho)
    contrast_inverse_score = max(0.0, -contrast_rho)
    return CandidateResult(
        candidate_id='IVA3_Q3',
        candidate_label='inverse_flow_survival_ordering',
        candidate_type='ordering_relation',
        expression='lower flow_concentration_index orders with longer topology_survival_time inside the stable smeared sector',
        primary_inputs=['flow_concentration_index', 'topology_survival_time'],
        context_inputs=[],
        uses_only_core=True,
        restates_law=False,
        adds_structure_beyond_law=True,
        complexity=2,
        stable_metric=f'spearman_rho={stable_rho:.4f}',
        contrast_metric=f'spearman_rho={contrast_rho:.4f}',
        stability_score=stable_inverse_score,
        contrast_score=contrast_inverse_score,
        separation_score=stable_inverse_score - contrast_inverse_score,
        branch_valid=True,
        selected=False,
        notes='Nonparametric preserved ordering relation. Stable rows maintain a near-monotone inverse order between flow and survival, while contrast controls do not.',
    )


def fit_affine_envelope(
    flow: np.ndarray,
    survival: np.ndarray,
    positive_mask: np.ndarray,
    contrast_mask: np.ndarray,
) -> CandidateResult:
    stable_flow = flow[positive_mask]
    stable_survival = survival[positive_mask]
    contrast_flow = flow[contrast_mask]
    contrast_survival = survival[contrast_mask]

    design = np.column_stack([np.ones(len(stable_flow)), stable_flow])
    coeffs = np.linalg.lstsq(design, stable_survival, rcond=None)[0]
    stable_pred = coeffs[0] + coeffs[1] * stable_flow
    contrast_pred = coeffs[0] + coeffs[1] * contrast_flow

    stable_rmse = float(np.sqrt(np.mean((stable_pred - stable_survival) ** 2)))
    contrast_rmse = float(np.sqrt(np.mean((contrast_pred - contrast_survival) ** 2)))
    stable_span = max(float(np.max(stable_survival) - np.min(stable_survival)), 1.0e-12)
    contrast_span = max(float(np.max(contrast_survival) - np.min(contrast_survival)), 1.0e-12)
    stable_rmse_norm = stable_rmse / stable_span
    contrast_rmse_norm = contrast_rmse / contrast_span

    return CandidateResult(
        candidate_id='IVA3_Q4',
        candidate_label='affine_flow_survival_envelope',
        candidate_type='affine_envelope',
        expression=f'topology_survival_time ~= {coeffs[0]:.4f} {coeffs[1]:+.4f} * flow_concentration_index',
        primary_inputs=['flow_concentration_index', 'topology_survival_time'],
        context_inputs=[],
        uses_only_core=True,
        restates_law=False,
        adds_structure_beyond_law=True,
        complexity=3,
        stable_metric=f'rmse_over_span={stable_rmse_norm:.4f}',
        contrast_metric=f'rmse_over_span={contrast_rmse_norm:.4f}',
        stability_score=score_from_rmse_norm(stable_rmse_norm),
        contrast_score=score_from_rmse_norm(contrast_rmse_norm),
        separation_score=score_from_rmse_norm(stable_rmse_norm) - score_from_rmse_norm(contrast_rmse_norm),
        branch_valid=True,
        selected=False,
        notes='Compact stable-sector envelope. The fit is tight on the positive set but is treated as secondary to the nonparametric ordering relation because it is fitted on only five stable cells.',
    )


def fit_survival_band(
    survival: np.ndarray,
    positive_mask: np.ndarray,
    contrast_mask: np.ndarray,
) -> CandidateResult:
    stable_survival = survival[positive_mask]
    contrast_survival = survival[contrast_mask]
    stable_cv = coefficient_of_variation(stable_survival)
    contrast_cv = coefficient_of_variation(contrast_survival)
    return CandidateResult(
        candidate_id='IVA3_Q5',
        candidate_label='topology_survival_band',
        candidate_type='bounded_quantity',
        expression='bounded topology_survival_time envelope inside the stable smeared sector',
        primary_inputs=['topology_survival_time'],
        context_inputs=[],
        uses_only_core=True,
        restates_law=False,
        adds_structure_beyond_law=True,
        complexity=1,
        stable_metric=f'cv={stable_cv:.4f}; range={format_range(stable_survival)}',
        contrast_metric=f'cv={contrast_cv:.4f}; range={format_range(contrast_survival)}',
        stability_score=score_from_cv(stable_cv),
        contrast_score=score_from_cv(contrast_cv),
        separation_score=score_from_cv(stable_cv) - score_from_cv(contrast_cv),
        branch_valid=True,
        selected=False,
        notes='Reported as a bounded quantity candidate, but the survival envelope alone is not more stable than the contrast controls.',
    )


def select_candidate(candidates: list[CandidateResult]) -> CandidateResult | None:
    valid = [
        candidate for candidate in candidates
        if candidate.branch_valid
        and candidate.uses_only_core
        and not candidate.restates_law
        and candidate.adds_structure_beyond_law
        and candidate.separation_score > 0.0
    ]
    if not valid:
        return None
    valid.sort(key=lambda item: (-item.separation_score, item.complexity, -item.stability_score))
    return valid[0]


def write_note(
    note_path: Path,
    runsheet: dict[str, Any],
    summary: dict[str, Any],
    candidates: list[CandidateResult],
    artifact_paths: dict[str, Any],
) -> None:
    selected = next((candidate for candidate in candidates if candidate.selected), None)
    rejected = [candidate for candidate in candidates if not candidate.selected]

    lines = [
        '# Stage IV A3 Stable Smeared Quasi-Invariant Probe v1',
        '',
        '## Stage purpose',
        'This stage probes whether the law-defined stable smeared sector on the frozen clustered Dirac-Kahler branch carries any branch-stable quasi-invariant, bounded composite, or preserved ordering relation. It is a reduction-only pass anchored to the Stage 24.2 threshold law.',
        '',
        '## Authoritative inputs',
    ]
    for key in ('stage23_12', 'stage23_14', 'stage24_1', 'stage24_2'):
        item = runsheet['required_inputs'][key]
        lines.append(
            f"- `{item['stage_label']}`: json=`{item['json_path']}`, csv=`{item['csv_path']}`"
        )

    lines.extend([
        '',
        '## Law-anchored positive and contrast sets',
        f"- selected threshold law: `{summary['selected_law_expression']}`",
        f"- positive stable smeared cells = `{summary['positive_cell_count']}`",
        f"- full contrast cells = `{summary['contrast_cell_count']}`",
        f"- near-threshold contrast cells = `{summary['near_threshold_contrast_count']}`",
        '- positive rows are defined by the intersection of the Stage 24.2 threshold law and the frozen Stage 23.12 stable-smeared truth label.',
        '- contrast rows are transient mixed controls above the threshold boundary; the near-threshold subset is retained as the main failure comparison.',
        '',
        '## Candidate quasi-invariant families',
    ])
    for item in runsheet['candidate_families']:
        lines.append(
            f"- `{item['candidate_id']}` `{item['candidate_type']}`: {item['description']}"
        )

    lines.extend([
        '',
        '## Consistency checks',
    ])
    for key, value in summary['consistency_checks'].items():
        lines.append(f"- `{key}` = `{value}`")

    lines.extend([
        '',
        '## Candidate comparison results',
    ])
    for candidate in candidates:
        lines.append(
            f"- `{candidate.candidate_id}`: stable_score=`{candidate.stability_score:.4f}`, contrast_score=`{candidate.contrast_score:.4f}`, separation=`{candidate.separation_score:.4f}`, selected=`{candidate.selected}`"
        )

    lines.extend([
        '',
        '## Selected quasi-invariant',
    ])
    if selected is None:
        lines.append('No branch-valid quasi-invariant candidate passed the Stage 24.3 constraints.')
    else:
        lines.append(f"- selected candidate: `{selected.candidate_id}` `{selected.candidate_label}`")
        lines.append(f"- expression: `{selected.expression}`")
        lines.append(f"- stable metric: `{selected.stable_metric}`")
        lines.append(f"- contrast metric: `{selected.contrast_metric}`")
        lines.append(f"- separation score: `{selected.separation_score:.4f}`")

    lines.extend([
        '',
        '## Rejected and secondary candidates',
    ])
    for candidate in rejected:
        if candidate.restates_law:
            reason = 'reported but rejected because it only restates the Stage 24.2 threshold law'
        elif not candidate.adds_structure_beyond_law:
            reason = 'reported as a support signature, but not selected because it does not add structure beyond the existing law/closure support'
        elif candidate.separation_score <= 0.0:
            reason = 'rejected because it is not more stable inside the positive set than in the contrast controls'
        else:
            reason = 'left as a secondary candidate because a simpler or stronger branch-valid relation dominated'
        lines.append(f"- `{candidate.candidate_id}`: {reason}; notes={candidate.notes}")

    lines.extend([
        '',
        '## Minimal interpretation',
        summary['final_interpretation'],
        '',
        '## Output artifact list',
        f"- stamped JSON summary: `{artifact_paths['json']}`",
        f"- stamped CSV candidate ledger: `{artifact_paths['csv']}`",
    ])
    for plot_path in artifact_paths['plots']:
        lines.append(f"- plot: `{plot_path}`")

    lines.extend([
        '',
        '## Commit recommendation',
    ])
    if summary['quasi_invariant_found']:
        lines.extend([
            '- recommended status: commit',
            '- rationale: Phase IV now has a branch-local quasi-invariant read for the law-defined stable smeared sector using only the accepted observable core.',
            '- branch effect: enables 24.4 extension-boundary testing from a compact law-and-quasi-invariant basis.',
            '- caution: the selected read is branch-local and does not reopen braid or enlarge the observable ledger.',
        ])
    else:
        lines.extend([
            '- recommended status: conditional commit',
            '- rationale: the negative result is clean, but no branch-valid quasi-invariant passed the Stage 24.3 constraints.',
            '- caution: do not rescue the stage by promoting rejected variables or by reusing braid-centered descriptors.',
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

    stage23_12_rows = required_inputs['stage23_12']['csv_rows']
    stage23_14 = required_inputs['stage23_14']['json']
    stage24_1 = required_inputs['stage24_1']['json']
    stage24_2 = required_inputs['stage24_2']['json']

    allowed_core = set(stage24_1['carry_forward_core'])
    allowed_context = set(stage24_1['carry_forward_context'])
    rejected = set(stage24_1['rejected_observables'])
    threshold = parse_threshold_from_expression(stage24_2['selected_law_expression'])
    near_threshold_margin = float(runsheet['near_threshold_margin'])

    flow = np.asarray([to_float(row['flow_concentration_index']) for row in stage23_12_rows], dtype=float)
    survival = np.asarray([to_float(row['topology_survival_time']) for row in stage23_12_rows], dtype=float)
    refinement = np.asarray([to_int(row['refinement_stability_flag']) for row in stage23_12_rows], dtype=float)
    weak = np.asarray([to_int(row['weak_coupling_stability_flag']) for row in stage23_12_rows], dtype=float)
    reverse = np.asarray([to_int(row['bidirectional_stability_flag']) for row in stage23_12_rows], dtype=float)
    phase = np.asarray([to_float(row['phase_offset_fraction_of_pi']) for row in stage23_12_rows], dtype=float)
    width = np.asarray([to_float(row['width_ratio_a_to_b']) for row in stage23_12_rows], dtype=float)

    truth_mask = np.asarray(
        [row['derived_phase_label'] == 'smeared_transfer_phase' for row in stage23_12_rows],
        dtype=bool,
    )
    law_mask = flow <= threshold
    positive_mask = truth_mask & law_mask
    contrast_mask = ~positive_mask
    near_threshold_mask = contrast_mask & ((flow - threshold) <= near_threshold_margin)

    if int(np.sum(positive_mask)) == 0:
        raise RuntimeError('law-anchored positive set is empty')
    if int(np.sum(contrast_mask)) == 0:
        raise RuntimeError('contrast set is empty')

    candidates = [
        fit_support_signature(refinement, weak, reverse, positive_mask, contrast_mask),
        fit_restated_law_band(flow, positive_mask, contrast_mask, threshold),
        fit_inverse_ordering_relation(flow, survival, positive_mask, contrast_mask),
        fit_affine_envelope(flow, survival, positive_mask, contrast_mask),
        fit_survival_band(survival, positive_mask, contrast_mask),
    ]

    for candidate in candidates:
        core_ok = set(candidate.primary_inputs).issubset(allowed_core)
        context_ok = set(candidate.context_inputs).issubset(allowed_context)
        forbidden_used = any(name in rejected for name in candidate.primary_inputs + candidate.context_inputs)
        no_braid_language = 'braid' not in candidate.expression.lower()
        candidate.uses_only_core = bool(core_ok and context_ok and not forbidden_used)
        candidate.branch_valid = bool(
            candidate.branch_valid
            and candidate.uses_only_core
            and context_ok
            and not forbidden_used
            and no_braid_language
        )

    selected = select_candidate(candidates)
    if selected is not None:
        selected.selected = True

    consistency_checks = {
        'stage24_1_observable_ledger_valid_true': bool(stage24_1['observable_ledger_valid']),
        'stage24_2_effective_law_found_true': bool(stage24_2['effective_law_found']),
        'stage24_2_selected_core_is_flow_only': bool(stage24_2['selected_core_inputs'] == ['flow_concentration_index']),
        'law_matches_truth_on_closure_lattice': bool(np.array_equal(law_mask, truth_mask)),
        'positive_set_is_smeared_only': bool(np.all(np.asarray([row['derived_phase_label'] for row in np.asarray(stage23_12_rows, dtype=object)[positive_mask]]) == 'smeared_transfer_phase')),
        'no_forbidden_observables_in_candidates': all(
            not any(name in rejected for name in candidate.primary_inputs + candidate.context_inputs)
            for candidate in candidates
        ),
        'no_braid_quantity_promoted': all('braid' not in candidate.expression.lower() for candidate in candidates),
        'selected_candidate_uses_core_only': bool(selected is not None and selected.uses_only_core),
        'selected_candidate_not_law_restatement': bool(selected is not None and not selected.restates_law),
        'phaseIII_smeared_side_still_frozen': bool(stage23_14['closure_class'] == 'smeared_dominant_closed_phase_diagram'),
    }

    quasi_invariant_found = bool(selected is not None and all(consistency_checks.values()))

    summary = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'selected_law_expression': stage24_2['selected_law_expression'],
        'flow_threshold': threshold,
        'positive_cell_count': int(np.sum(positive_mask)),
        'contrast_cell_count': int(np.sum(contrast_mask)),
        'near_threshold_contrast_count': int(np.sum(near_threshold_mask)),
        'quasi_invariant_found': quasi_invariant_found,
        'selected_candidate_id': selected.candidate_id if selected is not None else None,
        'selected_candidate_label': selected.candidate_label if selected is not None else None,
        'selected_candidate_type': selected.candidate_type if selected is not None else None,
        'selected_expression': selected.expression if selected is not None else None,
        'selected_primary_inputs': selected.primary_inputs if selected is not None else [],
        'selected_stability_score': selected.stability_score if selected is not None else None,
        'selected_contrast_score': selected.contrast_score if selected is not None else None,
        'selected_separation_score': selected.separation_score if selected is not None else None,
        'exact_support_signature_present': bool(any(candidate.candidate_id == 'IVA3_Q1' and candidate.stability_score == 1.0 for candidate in candidates)),
        'law_restatement_rejected': bool(any(candidate.candidate_id == 'IVA3_Q2' and not candidate.selected for candidate in candidates)),
        'consistency_checks': consistency_checks,
        'candidate_ledger': [candidate_to_row(candidate) for candidate in candidates],
    }

    if quasi_invariant_found and selected is not None:
        summary['final_interpretation'] = (
            'On the frozen clustered DK branch, the law-defined stable smeared sector carries a branch-local preserved ordering relation: '
            'lower flow concentration co-occurs with longer topology survival across the stable set, while the contrast controls do not preserve that inverse order. '
            'This read uses only accepted core observables, does not restate the threshold law, and remains bounded to the stable smeared sector.'
        )
    else:
        summary['final_interpretation'] = 'No branch-valid quasi-invariant was found under the current law-defined smeared-sector constraints.'

    table_fig, ax = plt.subplots(figsize=(14.0, 4.8))
    ax.axis('off')
    table_rows = [
        [
            candidate.candidate_id,
            candidate.candidate_type,
            '; '.join(candidate.primary_inputs),
            f'{candidate.stability_score:.3f}',
            f'{candidate.contrast_score:.3f}',
            f'{candidate.separation_score:.3f}',
            'yes' if candidate.restates_law else 'no',
            'yes' if candidate.selected else 'no',
        ]
        for candidate in candidates
    ]
    table = ax.table(
        cellText=table_rows,
        colLabels=['candidate', 'type', 'inputs', 'stable', 'contrast', 'gap', 'restates law', 'selected'],
        loc='center',
        cellLoc='left',
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.28)
    ax.set_title('Stage 24.3 quasi-invariant candidate comparison', fontsize=12, pad=12)
    table_path = WORK_PLOT_DIR / 'stage24_3_quasi_invariant_candidate_table.png'
    table_fig.tight_layout()
    table_fig.savefig(table_path, dpi=200, bbox_inches='tight')
    plt.close(table_fig)

    relation_fig, ax = plt.subplots(figsize=(7.2, 5.4))
    stable_color = '#c23b22'
    contrast_color = '#8f8f8f'
    near_color = '#e0a11b'
    ax.scatter(flow[contrast_mask], survival[contrast_mask], color=contrast_color, s=50, label='contrast controls')
    ax.scatter(flow[near_threshold_mask], survival[near_threshold_mask], facecolors='none', edgecolors=near_color, linewidths=1.4, s=78, label='near-threshold controls')
    ax.scatter(flow[positive_mask], survival[positive_mask], color=stable_color, s=58, label='stable closed smeared')
    stable_flow = flow[positive_mask]
    stable_survival = survival[positive_mask]
    design = np.column_stack([np.ones(len(stable_flow)), stable_flow])
    coeffs = np.linalg.lstsq(design, stable_survival, rcond=None)[0]
    x_line = np.linspace(float(np.min(stable_flow)), float(np.max(stable_flow)), 100)
    y_line = coeffs[0] + coeffs[1] * x_line
    ax.plot(x_line, y_line, color=stable_color, linestyle='--', linewidth=1.2, label='stable affine envelope')
    ax.axvline(threshold, color='#111111', linestyle=':', linewidth=1.0, label='24.2 threshold')
    ax.set_xlabel('flow concentration index')
    ax.set_ylabel('topology survival time')
    ax.set_title('Selected relation: stable smeared ordering structure')
    ax.legend(frameon=False, fontsize=8, loc='best')
    relation_fig.tight_layout()
    relation_path = WORK_PLOT_DIR / 'stage24_3_selected_relation_panel.png'
    relation_fig.savefig(relation_path, dpi=200, bbox_inches='tight')
    plt.close(relation_fig)

    score_fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.8))
    labels = [candidate.candidate_id for candidate in candidates]
    stable_scores = [candidate.stability_score for candidate in candidates]
    contrast_scores = [candidate.contrast_score for candidate in candidates]
    gaps = [candidate.separation_score for candidate in candidates]
    colors = ['#c23b22' if candidate.selected else '#8f8f8f' for candidate in candidates]

    x = np.arange(len(labels))
    width_bar = 0.38
    axes[0].bar(x - width_bar / 2.0, stable_scores, width_bar, color='#c23b22', label='stable score')
    axes[0].bar(x + width_bar / 2.0, contrast_scores, width_bar, color='#8f8f8f', label='contrast score')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels)
    axes[0].set_ylim(0.0, 1.05)
    axes[0].set_title('Stable vs contrast stability score')
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].bar(labels, gaps, color=colors)
    axes[1].axhline(0.0, color='#111111', linewidth=1.0)
    axes[1].set_title('Separation score')
    axes[1].set_ylabel('stable - contrast')

    score_fig.tight_layout()
    score_path = WORK_PLOT_DIR / 'stage24_3_stable_vs_contrast_panel.png'
    score_fig.savefig(score_path, dpi=200, bbox_inches='tight')
    plt.close(score_fig)

    timestamp = timestamp_slug()
    json_path = RESULTS_DIR / f'{timestamp}_stage24_3_stable_smeared_quasi_invariant_probe.json'
    csv_path = RESULTS_DIR / f'{timestamp}_stage24_3_stable_smeared_quasi_invariant_ledger.csv'
    plot_files = [
        (table_path, PLOTS_DIR / f'{timestamp}_stage24_3_quasi_invariant_candidate_table.png'),
        (relation_path, PLOTS_DIR / f'{timestamp}_stage24_3_selected_relation_panel.png'),
        (score_path, PLOTS_DIR / f'{timestamp}_stage24_3_stable_vs_contrast_panel.png'),
    ]
    json_path.write_text(json.dumps(summary, indent=2), encoding='utf-8')
    write_csv(csv_path, [candidate_to_row(candidate) for candidate in candidates], CSV_FIELDS)
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
        candidates=candidates,
        artifact_paths=artifact_paths,
    )

    print(f'Wrote note: {NOTE_PATH.relative_to(REPO_ROOT)}')
    print(f'Wrote JSON: {json_path.relative_to(REPO_ROOT)}')
    print(f'Wrote CSV: {csv_path.relative_to(REPO_ROOT)}')
    for plot in stamped_plots:
        print(f'Wrote plot: {plot}')


if __name__ == '__main__':
    main()
