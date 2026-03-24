#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

from stage10_common import ATLAS_NOTES, REPO_ROOT, plt

RUNSHEET_PATH = REPO_ROOT / 'numerics' / 'simulations' / 'stage24_2_smeared_sector_effective_law_runs.json'
NOTE_PATH = ATLAS_NOTES / 'Stage_IV_A2_Smeared_Sector_Effective_Law_Extraction_v1.md'
WORK_PLOT_DIR = Path('/tmp') / 'haos_iip_stage24_2_smeared_sector_effective_law'
WORK_PLOT_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR = REPO_ROOT / 'data'
PLOTS_DIR = REPO_ROOT / 'plots'

CSV_FIELDS = [
    'model_id',
    'model_family',
    'target_definition',
    'core_inputs',
    'context_inputs',
    'forbidden_inputs_used',
    'parameter_count',
    'rule_count',
    'accuracy',
    'balanced_accuracy',
    'baseline_comparison',
    'compression_score',
    'interpretability_score',
    'branch_valid',
    'selected',
    'law_expression',
    'notes',
]


@dataclass
class ModelResult:
    model_id: str
    model_family: str
    target_definition: str
    core_inputs: list[str]
    context_inputs: list[str]
    forbidden_inputs_used: bool
    parameter_count: int
    rule_count: int
    accuracy: float
    balanced_accuracy: float
    baseline_comparison: float
    compression_score: float
    interpretability_score: float
    branch_valid: bool
    selected: bool
    law_expression: str
    notes: str
    predictions: np.ndarray
    score_values: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Run Stage 24.2 smeared-sector effective law extraction.')
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


def accuracy(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    return float(np.mean(y_true == y_pred))


def balanced_accuracy(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    positives = y_true == 1
    negatives = y_true == 0
    tpr = float(np.mean(y_pred[positives] == 1)) if np.any(positives) else 0.0
    tnr = float(np.mean(y_pred[negatives] == 0)) if np.any(negatives) else 0.0
    return 0.5 * (tpr + tnr)


def majority_baseline(y_true: np.ndarray) -> tuple[np.ndarray, float, float]:
    majority = 1 if int(np.sum(y_true)) > (len(y_true) / 2.0) else 0
    pred = np.full_like(y_true, fill_value=majority)
    return pred, accuracy(y_true, pred), balanced_accuracy(y_true, pred)


def threshold_candidates(values: np.ndarray) -> list[float]:
    unique = sorted(set(float(value) for value in values))
    if len(unique) == 1:
        return [unique[0]]
    return [0.5 * (left + right) for left, right in zip(unique[:-1], unique[1:])]


def fit_single_threshold_feature(
    feature_arrays: dict[str, np.ndarray],
    y_true: np.ndarray,
    target_definition: str,
    baseline_bal_acc: float,
) -> ModelResult:
    best_pred = None
    best_expr = ''
    best_bal_acc = -1.0
    best_acc = -1.0
    best_threshold = 0.0
    best_feature_name = ''
    best_values = None
    best_is_continuous = False

    for feature_name, values in feature_arrays.items():
        is_continuous = len(set(float(value) for value in values)) > 2
        for threshold in threshold_candidates(values):
            for direction in ('<=', '>='):
                pred = (values <= threshold).astype(int) if direction == '<=' else (values >= threshold).astype(int)
                bal_acc = balanced_accuracy(y_true, pred)
                acc = accuracy(y_true, pred)
                if (bal_acc, acc, is_continuous) > (best_bal_acc, best_acc, best_is_continuous):
                    best_pred = pred
                    best_bal_acc = bal_acc
                    best_acc = acc
                    best_threshold = float(threshold)
                    best_feature_name = feature_name
                    best_values = values.copy()
                    best_is_continuous = is_continuous
                    best_expr = f'stable_closed_smeared = 1 iff {feature_name} {direction} {threshold:.6f}'

    assert best_pred is not None and best_values is not None
    compression = max(0.0, best_bal_acc - baseline_bal_acc) / (1 + 1 + 1)
    interpretability = 1.0 / (1.0 + 1 + 1)
    return ModelResult(
        model_id='IVA2_M1',
        model_family='threshold_rule',
        target_definition=target_definition,
        core_inputs=[best_feature_name],
        context_inputs=[],
        forbidden_inputs_used=False,
        parameter_count=1,
        rule_count=1,
        accuracy=best_acc,
        balanced_accuracy=best_bal_acc,
        baseline_comparison=best_bal_acc - baseline_bal_acc,
        compression_score=compression,
        interpretability_score=interpretability,
        branch_valid=True,
        selected=False,
        law_expression=best_expr,
        notes=f'Best single-threshold rule found on `{best_feature_name}` with threshold `{best_threshold:.6f}`.',
        predictions=best_pred,
        score_values=best_values,
    )


def fit_linear_surrogate(
    X: np.ndarray,
    feature_names: list[str],
    y_true: np.ndarray,
    target_definition: str,
    baseline_bal_acc: float,
) -> ModelResult:
    means = np.mean(X, axis=0)
    stds = np.std(X, axis=0)
    safe_stds = np.where(stds <= 1.0e-12, 1.0, stds)
    Xz = (X - means) / safe_stds
    design = np.column_stack([np.ones(len(Xz)), Xz])
    coeffs = np.linalg.lstsq(design, y_true.astype(float), rcond=None)[0]
    scores = design @ coeffs

    best_pred = None
    best_threshold = 0.5
    best_bal_acc = -1.0
    best_acc = -1.0
    for threshold in threshold_candidates(scores):
        pred = (scores >= threshold).astype(int)
        bal_acc = balanced_accuracy(y_true, pred)
        acc = accuracy(y_true, pred)
        if (bal_acc, acc) > (best_bal_acc, best_acc):
            best_pred = pred
            best_bal_acc = bal_acc
            best_acc = acc
            best_threshold = float(threshold)

    assert best_pred is not None
    coeff_terms = [f'{coeffs[0]:+.4f}']
    for name, coeff in zip(feature_names, coeffs[1:]):
        coeff_terms.append(f'{coeff:+.4f}*{name}')
    expr = 'score = ' + ' '.join(coeff_terms) + f'; stable_closed_smeared = 1 iff score >= {best_threshold:.6f}'
    compression = max(0.0, best_bal_acc - baseline_bal_acc) / (1 + len(feature_names) + 1)
    interpretability = 1.0 / (1.0 + len(feature_names) + 1)
    return ModelResult(
        model_id='IVA2_M2',
        model_family='simple_reduced_surrogate',
        target_definition=target_definition,
        core_inputs=feature_names[:],
        context_inputs=[],
        forbidden_inputs_used=False,
        parameter_count=len(feature_names) + 1,
        rule_count=1,
        accuracy=best_acc,
        balanced_accuracy=best_bal_acc,
        baseline_comparison=best_bal_acc - baseline_bal_acc,
        compression_score=compression,
        interpretability_score=interpretability,
        branch_valid=True,
        selected=False,
        law_expression=expr,
        notes='Least-squares linear score model over the accepted core observables.',
        predictions=best_pred,
        score_values=scores,
    )


def fit_context_threshold(
    phase_values: np.ndarray,
    feature_name: str,
    values: np.ndarray,
    y_true: np.ndarray,
    target_definition: str,
    baseline_bal_acc: float,
) -> ModelResult:
    split = float(np.median(sorted(set(float(value) for value in phase_values))))
    lower_mask = phase_values <= split
    upper_mask = phase_values > split

    def best_threshold(sub_values: np.ndarray, sub_targets: np.ndarray) -> tuple[float, np.ndarray]:
        best_t = float(np.mean(sub_values))
        best_pred = np.zeros_like(sub_targets)
        best_bal = -1.0
        best_acc = -1.0
        for threshold in threshold_candidates(sub_values):
            pred = (sub_values <= threshold).astype(int)
            bal_acc = balanced_accuracy(sub_targets, pred)
            acc = accuracy(sub_targets, pred)
            if (bal_acc, acc) > (best_bal, best_acc):
                best_bal = bal_acc
                best_acc = acc
                best_t = float(threshold)
                best_pred = pred
        return best_t, best_pred

    lower_t, lower_pred = best_threshold(values[lower_mask], y_true[lower_mask])
    upper_t, upper_pred = best_threshold(values[upper_mask], y_true[upper_mask])

    pred = np.zeros_like(y_true)
    pred[lower_mask] = lower_pred
    pred[upper_mask] = upper_pred
    bal_acc = balanced_accuracy(y_true, pred)
    acc = accuracy(y_true, pred)
    expr = (
        f'if phase_corridor_position <= {split:.6f}: stable_closed_smeared = 1 iff {feature_name} <= {lower_t:.6f}; '
        f'else stable_closed_smeared = 1 iff {feature_name} <= {upper_t:.6f}'
    )
    compression = max(0.0, bal_acc - baseline_bal_acc) / (1 + 2 + 3)
    interpretability = 1.0 / (1.0 + 2 + 3)
    return ModelResult(
        model_id='IVA2_M3',
        model_family='context_conditioned_threshold_rule',
        target_definition=target_definition,
        core_inputs=[feature_name],
        context_inputs=['phase_corridor_position', 'phase_corridor_width'],
        forbidden_inputs_used=False,
        parameter_count=2,
        rule_count=3,
        accuracy=acc,
        balanced_accuracy=bal_acc,
        baseline_comparison=bal_acc - baseline_bal_acc,
        compression_score=compression,
        interpretability_score=interpretability,
        branch_valid=True,
        selected=False,
        law_expression=expr,
        notes='Two-bin context-conditioned threshold rule. Phase enters only as an external selector, not as a state variable.',
        predictions=pred,
        score_values=values.copy(),
    )


def fit_closure_conjunction(
    refinement: np.ndarray,
    weak: np.ndarray,
    reverse: np.ndarray,
    y_true: np.ndarray,
    target_definition: str,
    baseline_bal_acc: float,
) -> ModelResult:
    pred = ((refinement == 1) & (weak == 1) & (reverse == 1)).astype(int)
    bal_acc = balanced_accuracy(y_true, pred)
    acc = accuracy(y_true, pred)
    expr = (
        'stable_closed_smeared = 1 iff refinement_stability_flag = 1 '
        'and weak_coupling_stability_flag = 1 and reverse_stability_flag = 1'
    )
    compression = max(0.0, bal_acc - baseline_bal_acc) / (1 + 0 + 3)
    interpretability = 1.0 / (1.0 + 0 + 3)
    return ModelResult(
        model_id='IVA2_M4',
        model_family='closure_support_conjunction',
        target_definition=target_definition,
        core_inputs=[
            'refinement_stability_flag',
            'weak_coupling_stability_flag',
            'reverse_stability_flag',
        ],
        context_inputs=[],
        forbidden_inputs_used=False,
        parameter_count=0,
        rule_count=3,
        accuracy=acc,
        balanced_accuracy=bal_acc,
        baseline_comparison=bal_acc - baseline_bal_acc,
        compression_score=compression,
        interpretability_score=interpretability,
        branch_valid=True,
        selected=False,
        law_expression=expr,
        notes='Minimal conjunction over the binary closure-support observables.',
        predictions=pred,
        score_values=pred.astype(float),
    )


def select_model(models: list[ModelResult]) -> ModelResult | None:
    valid = [
        model for model in models
        if model.branch_valid and not model.forbidden_inputs_used and model.baseline_comparison > 0.0
    ]
    if not valid:
        return None
    valid.sort(
        key=lambda model: (
            model.parameter_count + model.rule_count,
            -model.interpretability_score,
            -model.balanced_accuracy,
            -model.compression_score,
        )
    )
    return valid[0]


def model_to_row(model: ModelResult) -> dict[str, Any]:
    return {
        'model_id': model.model_id,
        'model_family': model.model_family,
        'target_definition': model.target_definition,
        'core_inputs': '; '.join(model.core_inputs),
        'context_inputs': '; '.join(model.context_inputs),
        'forbidden_inputs_used': model.forbidden_inputs_used,
        'parameter_count': model.parameter_count,
        'rule_count': model.rule_count,
        'accuracy': f'{model.accuracy:.4f}',
        'balanced_accuracy': f'{model.balanced_accuracy:.4f}',
        'baseline_comparison': f'{model.baseline_comparison:.4f}',
        'compression_score': f'{model.compression_score:.4f}',
        'interpretability_score': f'{model.interpretability_score:.4f}',
        'branch_valid': model.branch_valid,
        'selected': model.selected,
        'law_expression': model.law_expression,
        'notes': model.notes,
    }


def write_note(
    note_path: Path,
    runsheet: dict[str, Any],
    summary: dict[str, Any],
    models: list[ModelResult],
    artifact_paths: dict[str, Any],
) -> None:
    selected = next((model for model in models if model.selected), None)
    failed = [model for model in models if not model.selected]

    lines = [
        '# Stage IV A2 Smeared Sector Effective Law Extraction v1',
        '',
        '## Stage purpose',
        'This stage runs the first true effective-law extraction pass for the surviving smeared-transfer sector on the frozen clustered DK branch. It compares a minimal family of reduced laws under the Stage 24.1 observable ledger and keeps only the smallest branch-valid law that achieves genuine predictive compression.',
        '',
        '## Authoritative inputs',
    ]
    for key in ('stage23_11', 'stage23_12', 'stage23_13', 'stage23_14', 'stage23_15', 'stage24_1'):
        item = runsheet['required_inputs'][key]
        if 'paper_path' in item:
            lines.append(f"- `{item['stage_label']}`: paper=`{item['paper_path']}`")
        else:
            lines.append(
                f"- `{item['stage_label']}`: json=`{item['json_path']}`, csv=`{item['csv_path']}`"
            )

    lines.extend([
        '',
        '## Law target definition',
        summary['target_definition'],
        '',
        '## Candidate model families',
    ])
    for model in models:
        lines.append(
            f"- `{model.model_id}` `{model.model_family}`: core=`{', '.join(model.core_inputs)}` context=`{', '.join(model.context_inputs) if model.context_inputs else 'none'}`"
        )

    lines.extend([
        '',
        '## Consistency checks',
    ])
    for key, value in summary['consistency_checks'].items():
        lines.append(f"- `{key}` = `{value}`")

    lines.extend([
        '',
        '## Model comparison results',
    ])
    for model in models:
        lines.append(
            f"- `{model.model_id}`: balanced_accuracy=`{model.balanced_accuracy:.4f}`, baseline_delta=`{model.baseline_comparison:.4f}`, complexity=`{model.parameter_count + model.rule_count}`, branch_valid=`{model.branch_valid}`, selected=`{model.selected}`"
        )

    lines.extend([
        '',
        '## Selected effective law',
    ])
    if selected is None:
        lines.append('No branch-valid compact law passed the Stage 24.2 constraints.')
    else:
        lines.append(f"- selected model: `{selected.model_id}`")
        lines.append(f"- law: `{selected.law_expression}`")
        lines.append(f"- core inputs: `{', '.join(selected.core_inputs)}`")
        lines.append(
            f"- fit: accuracy=`{selected.accuracy:.4f}`, balanced_accuracy=`{selected.balanced_accuracy:.4f}`, baseline_delta=`{selected.baseline_comparison:.4f}`"
        )

    lines.extend([
        '',
        '## Failed candidates and why they failed',
    ])
    for model in failed:
        reason = 'unselected because a simpler branch-valid model achieved equal or better compression'
        if not model.branch_valid:
            reason = 'failed branch-validity checks'
        elif model.baseline_comparison <= 0.0:
            reason = 'failed to beat the naive baseline'
        lines.append(f"- `{model.model_id}`: {reason}; notes={model.notes}")

    lines.extend([
        '',
        '## Minimal interpretation',
    ])
    if selected is None:
        lines.append('No compact law was found under the current observable ledger.')
    else:
        lines.append(summary['final_interpretation'])

    lines.extend([
        '',
        '## Output artifact list',
        f"- stamped JSON summary: `{artifact_paths['json']}`",
        f"- stamped CSV model ledger: `{artifact_paths['csv']}`",
    ])
    for plot_path in artifact_paths['plots']:
        lines.append(f"- plot: `{plot_path}`")

    lines.extend([
        '',
        '## Commit recommendation',
    ])
    if summary['effective_law_found']:
        lines.extend([
            '- recommended status: commit',
            '- rationale: Phase IV now has a valid reduced law candidate for the smeared sector using only the frozen observable core.',
            '- branch effect: enables 24.3 quasi-invariant probing from a compact law basis.',
            '- caution: rejected variables remain excluded unless a later explicit branch changes the measurement regime.',
        ])
    else:
        lines.extend([
            '- recommended status: conditional commit',
            '- rationale: the negative result is clean, but no compact law passed the current ledger constraints.',
            '- caution: do not rescue the law by enlarging the observable basis on this branch.',
        ])

    note_path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main() -> None:
    args = parse_args()
    runsheet = load_json(args.runsheet)

    required_inputs: dict[str, dict[str, Any]] = {}
    for key, item in runsheet['required_inputs'].items():
        if 'paper_path' in item:
            paper_path = REPO_ROOT / item['paper_path']
            if not paper_path.exists():
                raise FileNotFoundError(f'missing required paper input: {paper_path}')
            required_inputs[key] = {
                'stage_label': item['stage_label'],
                'paper_path': str(paper_path.relative_to(REPO_ROOT)),
            }
            continue

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

    stage24_1 = required_inputs['stage24_1']['json']
    stage23_12_rows = required_inputs['stage23_12']['csv_rows']
    stage23_14 = required_inputs['stage23_14']['json']

    allowed_core = list(stage24_1['carry_forward_core'])
    allowed_context = list(stage24_1['carry_forward_context'])
    rejected = set(stage24_1['rejected_observables'])

    y_true = np.asarray(
        [1 if row['derived_phase_label'] == 'smeared_transfer_phase' else 0 for row in stage23_12_rows],
        dtype=int,
    )
    phase_position = np.asarray([to_float(row['phase_offset_fraction_of_pi']) for row in stage23_12_rows], dtype=float)
    phase_corridor_width = np.full_like(phase_position, fill_value=float(stage24_1['phase_corridor_width']), dtype=float)
    topology_survival_time = np.asarray([to_float(row['topology_survival_time']) for row in stage23_12_rows], dtype=float)
    refinement = np.asarray([to_int(row['refinement_stability_flag']) for row in stage23_12_rows], dtype=float)
    weak = np.asarray([to_int(row['weak_coupling_stability_flag']) for row in stage23_12_rows], dtype=float)
    reverse = np.asarray([to_int(row['bidirectional_stability_flag']) for row in stage23_12_rows], dtype=float)
    flow = np.asarray([to_float(row['flow_concentration_index']) for row in stage23_12_rows], dtype=float)
    width_ratio = np.asarray([to_float(row['width_ratio_a_to_b']) for row in stage23_12_rows], dtype=float)

    feature_map = {
        'topology_survival_time': topology_survival_time,
        'refinement_stability_flag': refinement,
        'weak_coupling_stability_flag': weak,
        'reverse_stability_flag': reverse,
        'flow_concentration_index': flow,
    }

    target_definition = runsheet['target_definition']
    baseline_pred, baseline_acc, baseline_bal_acc = majority_baseline(y_true)

    model_m1 = fit_single_threshold_feature(
        feature_arrays=feature_map,
        y_true=y_true,
        target_definition=target_definition,
        baseline_bal_acc=baseline_bal_acc,
    )
    model_m2 = fit_linear_surrogate(
        X=np.column_stack([feature_map[name] for name in allowed_core]),
        feature_names=allowed_core,
        y_true=y_true,
        target_definition=target_definition,
        baseline_bal_acc=baseline_bal_acc,
    )
    model_m3 = fit_context_threshold(
        phase_values=phase_position,
        feature_name=model_m1.core_inputs[0],
        values=feature_map[model_m1.core_inputs[0]],
        y_true=y_true,
        target_definition=target_definition,
        baseline_bal_acc=baseline_bal_acc,
    )
    model_m4 = fit_closure_conjunction(
        refinement=refinement,
        weak=weak,
        reverse=reverse,
        y_true=y_true,
        target_definition=target_definition,
        baseline_bal_acc=baseline_bal_acc,
    )

    models = [model_m1, model_m2, model_m3, model_m4]

    for model in models:
        core_ok = set(model.core_inputs).issubset(set(allowed_core))
        context_ok = set(model.context_inputs).issubset(set(allowed_context))
        forbidden_used = any(name in rejected for name in model.core_inputs + model.context_inputs)
        simpler_than_lookup = (model.parameter_count + model.rule_count) < len(y_true)
        beats_baseline = model.balanced_accuracy > baseline_bal_acc
        no_braid_language = 'braid' not in model.law_expression.lower()
        smeared_side_consistent = stage23_14['closure_class'] == 'smeared_dominant_closed_phase_diagram'
        model.forbidden_inputs_used = forbidden_used
        model.branch_valid = bool(
            core_ok
            and context_ok
            and not forbidden_used
            and simpler_than_lookup
            and beats_baseline
            and no_braid_language
            and smeared_side_consistent
        )

    selected = select_model(models)
    if selected is not None:
        selected.selected = True

    consistency_checks = {
        'stage23_11_effective_model_valid_true': bool(required_inputs['stage23_11']['json']['summary']['effective_model_valid']),
        'stage23_13_phase_corridor_primary': bool(required_inputs['stage23_13']['json']['overall_summary']['overall_mechanism_label'] == 'phase_corridor_primary'),
        'stage23_14_freeze_valid_true': bool(stage23_14['phaseIII_freeze_valid']),
        'only_accepted_primary_observables_used': all(
            set(model.core_inputs).issubset(set(allowed_core)) for model in models
        ),
        'context_variables_not_promoted': all(
            set(model.context_inputs).issubset(set(allowed_context)) for model in models
        ),
        'rejected_variables_unused': all(
            not any(name in rejected for name in model.core_inputs + model.context_inputs) for model in models
        ),
        'winning_law_simpler_than_lookup_table': bool(selected is not None and (selected.parameter_count + selected.rule_count) < len(y_true)),
        'winning_law_beats_naive_baseline': bool(selected is not None and selected.balanced_accuracy > baseline_bal_acc),
        'winning_law_smeared_side_consistent': bool(selected is not None and stage23_14['closure_class'] == 'smeared_dominant_closed_phase_diagram'),
    }

    effective_law_found = bool(selected is not None and all(consistency_checks.values()))
    summary = {
        'experiment': runsheet['stage'],
        'description': runsheet['description'],
        'target_definition': target_definition,
        'effective_law_found': effective_law_found,
        'minimal_model_selected': bool(selected is not None),
        'uses_only_ledger_core': consistency_checks['only_accepted_primary_observables_used'],
        'context_variables_not_promoted': consistency_checks['context_variables_not_promoted'],
        'rejected_variables_unused': consistency_checks['rejected_variables_unused'],
        'predictive_compression_valid': bool(selected is not None and selected.balanced_accuracy > baseline_bal_acc),
        'baseline_accuracy': baseline_acc,
        'baseline_balanced_accuracy': baseline_bal_acc,
        'selected_model_id': selected.model_id if selected is not None else None,
        'selected_model_family': selected.model_family if selected is not None else None,
        'selected_law_expression': selected.law_expression if selected is not None else None,
        'selected_core_inputs': selected.core_inputs if selected is not None else [],
        'selected_context_inputs': selected.context_inputs if selected is not None else [],
        'selected_accuracy': selected.accuracy if selected is not None else None,
        'selected_balanced_accuracy': selected.balanced_accuracy if selected is not None else None,
        'selected_baseline_delta': selected.baseline_comparison if selected is not None else None,
        'consistency_checks': consistency_checks,
        'model_ledger': [model_to_row(model) for model in models],
    }

    if effective_law_found and selected is not None:
        summary['final_interpretation'] = (
            'On the frozen clustered DK branch, the smeared sector admits a compact effective law built from the accepted Phase IV observable core. '
            'The selected law predicts branch-valid smeared-sector persistence/stability without reintroducing rejected braid- or asymmetry-based quantities. '
            'Phase-corridor variables remain external conditioning coordinates rather than promoted state variables.'
        )
    else:
        summary['final_interpretation'] = 'No compact law was found under the current observable ledger.'

    table_fig, ax = plt.subplots(figsize=(13.5, 4.8))
    ax.axis('off')
    table_rows = [
        [
            model.model_id,
            model.model_family,
            '; '.join(model.core_inputs),
            '; '.join(model.context_inputs) if model.context_inputs else 'none',
            f'{model.balanced_accuracy:.3f}',
            f'{model.baseline_comparison:.3f}',
            f'{model.parameter_count + model.rule_count}',
            'yes' if model.branch_valid else 'no',
            'yes' if model.selected else 'no',
        ]
        for model in models
    ]
    table = ax.table(
        cellText=table_rows,
        colLabels=['model', 'family', 'core inputs', 'context', 'bal. acc.', 'delta', 'complexity', 'valid', 'selected'],
        loc='center',
        cellLoc='left',
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.0, 1.25)
    ax.set_title('Stage 24.2 law family comparison', fontsize=12, pad=12)
    table_path = WORK_PLOT_DIR / 'stage24_2_law_family_comparison_table.png'
    table_fig.tight_layout()
    table_fig.savefig(table_path, dpi=200, bbox_inches='tight')
    plt.close(table_fig)

    decision_fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    labels = [model.model_id for model in models]
    bal_scores = [model.balanced_accuracy for model in models]
    complexities = [model.parameter_count + model.rule_count for model in models]
    colors = ['#c23b22' if model.selected else '#8f8f8f' for model in models]
    axes[0].bar(labels, bal_scores, color=colors)
    axes[0].axhline(baseline_bal_acc, color='#1f1f1f', linestyle='--', linewidth=1, label='naive baseline')
    axes[0].set_ylabel('balanced accuracy')
    axes[0].set_title('Predictive adequacy')
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].bar(labels, complexities, color=colors)
    axes[1].set_ylabel('parameter + rule count')
    axes[1].set_title('Model complexity')

    if selected is not None:
        decision_fig.suptitle(f'Selected law: {selected.model_id} | {selected.law_expression}', fontsize=11)
    else:
        decision_fig.suptitle('No branch-valid minimal law selected', fontsize=11)
    decision_fig.tight_layout()
    decision_path = WORK_PLOT_DIR / 'stage24_2_minimal_law_decision_panel.png'
    decision_fig.savefig(decision_path, dpi=200, bbox_inches='tight')
    plt.close(decision_fig)

    phase_values = sorted(set(float(value) for value in phase_position))
    width_values = sorted(set(float(value) for value in width_ratio))
    phase_index = {value: idx for idx, value in enumerate(phase_values)}
    width_index = {value: idx for idx, value in enumerate(width_values)}
    truth_grid = np.zeros((len(width_values), len(phase_values)), dtype=float)
    pred_grid = np.zeros((len(width_values), len(phase_values)), dtype=float)
    selected_pred = selected.predictions if selected is not None else baseline_pred
    for row, truth, pred in zip(stage23_12_rows, y_true, selected_pred):
        i = width_index[to_float(row['width_ratio_a_to_b'])]
        j = phase_index[to_float(row['phase_offset_fraction_of_pi'])]
        truth_grid[i, j] = float(truth)
        pred_grid[i, j] = float(pred)

    pred_fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.6), sharey=True)
    for ax, grid, title in zip(axes, [truth_grid, pred_grid], ['truth', 'selected law prediction']):
        image = ax.imshow(grid, vmin=0.0, vmax=1.0, cmap='OrRd', aspect='auto')
        ax.set_title(title)
        ax.set_xticks(range(len(phase_values)))
        ax.set_xticklabels([f'{value:.3f}' for value in phase_values], rotation=45, ha='right', fontsize=8)
        ax.set_xlabel('phase corridor position')
        ax.set_yticks(range(len(width_values)))
        ax.set_yticklabels([f'{value:.2f}' for value in width_values], fontsize=8)
        ax.set_ylabel('width ratio')
        for i in range(len(width_values)):
            for j in range(len(phase_values)):
                ax.text(j, i, str(int(grid[i, j])), ha='center', va='center', color='#111111', fontsize=9)
    pred_fig.suptitle('Stage 24.2 prediction vs truth', fontsize=12)
    pred_fig.tight_layout()
    pred_path = WORK_PLOT_DIR / 'stage24_2_prediction_vs_truth_panel.png'
    pred_fig.savefig(pred_path, dpi=200, bbox_inches='tight')
    plt.close(pred_fig)

    context_fig, ax = plt.subplots(figsize=(7.2, 4.8))
    stable_mask = y_true == 1
    transient_mask = y_true == 0
    ax.scatter(phase_position[transient_mask], flow[transient_mask], color='#8f8f8f', s=50, label='transient mixed')
    ax.scatter(phase_position[stable_mask], flow[stable_mask], color='#c23b22', s=55, label='stable closed smeared')
    if selected is not None and selected.model_id == 'IVA2_M1':
        if '<=' in selected.law_expression:
            threshold = float(selected.law_expression.split('<=')[-1].strip())
            ax.axhline(threshold, color='#c23b22', linestyle='--', linewidth=1.2, label='selected threshold')
        elif '>=' in selected.law_expression:
            threshold = float(selected.law_expression.split('>=')[-1].strip())
            ax.axhline(threshold, color='#c23b22', linestyle='--', linewidth=1.2, label='selected threshold')
    ax.set_xlabel('phase corridor position')
    ax.set_ylabel('flow concentration index')
    ax.set_title('Context conditioning panel')
    ax.legend(frameon=False, fontsize=8)
    context_fig.tight_layout()
    context_path = WORK_PLOT_DIR / 'stage24_2_context_conditioning_panel.png'
    context_fig.savefig(context_path, dpi=200, bbox_inches='tight')
    plt.close(context_fig)

    timestamp = timestamp_slug()
    json_path = RESULTS_DIR / f'{timestamp}_stage24_2_smeared_sector_effective_law_extraction.json'
    csv_path = RESULTS_DIR / f'{timestamp}_stage24_2_smeared_sector_effective_law_model_ledger.csv'
    plot_files = [
        (table_path, PLOTS_DIR / f'{timestamp}_stage24_2_law_family_comparison_table.png'),
        (decision_path, PLOTS_DIR / f'{timestamp}_stage24_2_minimal_law_decision_panel.png'),
        (pred_path, PLOTS_DIR / f'{timestamp}_stage24_2_prediction_vs_truth_panel.png'),
        (context_path, PLOTS_DIR / f'{timestamp}_stage24_2_context_conditioning_panel.png'),
    ]
    json_path.write_text(json.dumps(summary, indent=2), encoding='utf-8')
    write_csv(csv_path, [model_to_row(model) for model in models], CSV_FIELDS)
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
        models=models,
        artifact_paths=artifact_paths,
    )

    print(f'Wrote note: {NOTE_PATH.relative_to(REPO_ROOT)}')
    print(f'Wrote JSON: {json_path.relative_to(REPO_ROOT)}')
    print(f'Wrote CSV: {csv_path.relative_to(REPO_ROOT)}')
    for plot in stamped_plots:
        print(f'Wrote plot: {plot}')


if __name__ == '__main__':
    main()
