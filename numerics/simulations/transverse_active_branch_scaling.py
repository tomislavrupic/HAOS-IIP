#!/usr/bin/env python3

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linear_sum_assignment

from L1_stage3_common import PLOTS, REPO_ROOT, append_log, save_result_payload, plt
from L1_transverse_band_test import analyze_case
from transverse_active_branch_identity import evaluate_pair
from transverse_active_sector_transport import normalized_matrix, pair_transport_metrics, transported_modes

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'transverse_active_branch_scaling': {
        'sizes': [12, 16, 20],
        'variants': ['baseline', 'line_defect', 'mild_disorder'],
        'restricted_modes': 6,
        'transport_modes': 4,
        'probe_n_side': 6,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
        'flux_tube_phase': math.pi / 2.0,
        'disorder_strength': 0.12,
        'candidate_exponents': [1.5, 1.75, 2.0, 2.25, 2.5],
        'thresholds': {
            'max_global_max_drift': 0.05,
            'max_best_exponent_deviation': 0.25,
            'min_identity_pairs': 2,
        },
    },
}

V1_THRESHOLDS = {
    'min_mean_overlap': 0.75,
    'min_min_overlap': 0.50,
    'min_mean_principal_cosine': 0.80,
    'max_max_scaled_eigen_drift': 0.05,
    'min_mean_assignment_margin': 0.10,
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['transverse_active_branch_scaling'] = dict(DEFAULT_CONFIG['transverse_active_branch_scaling'])
    merged['transverse_active_branch_scaling']['thresholds'] = dict(DEFAULT_CONFIG['transverse_active_branch_scaling']['thresholds'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_active_branch_scaling'})
        if isinstance(on_disk.get('transverse_active_branch_scaling'), dict):
            block = dict(on_disk['transverse_active_branch_scaling'])
            thresholds = block.pop('thresholds', None)
            merged['transverse_active_branch_scaling'].update(block)
            if isinstance(thresholds, dict):
                merged['transverse_active_branch_scaling']['thresholds'].update(thresholds)
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_active_branch_scaling'})
        if isinstance(config.get('transverse_active_branch_scaling'), dict):
            block = dict(config['transverse_active_branch_scaling'])
            thresholds = block.pop('thresholds', None)
            merged['transverse_active_branch_scaling'].update(block)
            if isinstance(thresholds, dict):
                merged['transverse_active_branch_scaling']['thresholds'].update(thresholds)
    return merged


def scaling_metrics(
    case_from: dict[str, Any],
    case_to: dict[str, Any],
    n_from: int,
    n_to: int,
    probe_n_side: int,
    transport_modes_count: int,
    exponent: float,
) -> dict[str, Any]:
    modes_from = transported_modes(case_from, probe_n_side, transport_modes_count)
    modes_to = transported_modes(case_to, probe_n_side, transport_modes_count)
    if not modes_from or not modes_to:
        return {
            'transport_modes': 0,
            'matched_pairs': [],
            'mean_scaled_eigen_drift': math.nan,
            'max_scaled_eigen_drift': math.nan,
        }

    matrix_from = normalized_matrix(modes_from)
    matrix_to = normalized_matrix(modes_to)
    overlap = np.abs(matrix_from.conj().T @ matrix_to)
    row_ind, col_ind = linear_sum_assignment(-overlap)
    order = np.argsort(row_ind)
    row_ind = row_ind[order]
    col_ind = col_ind[order]

    scaled_from = [(n_from**exponent) * float(value) for value in case_from['restricted_transverse_spectrum'][: len(modes_from)]]
    scaled_to = [(n_to**exponent) * float(value) for value in case_to['restricted_transverse_spectrum'][: len(modes_to)]]
    drifts: list[float] = []
    matched_pairs: list[dict[str, Any]] = []
    for left, right in zip(row_ind, col_ind):
        scaled_drift = abs(scaled_from[left] - scaled_to[right]) / max(abs(scaled_from[left]), abs(scaled_to[right]), 1.0e-12)
        matched_pairs.append(
            {
                'from_mode': int(left),
                'to_mode': int(right),
                'scaled_eigen_from': float(scaled_from[left]),
                'scaled_eigen_to': float(scaled_to[right]),
                'scaled_eigen_drift': float(scaled_drift),
            }
        )
        drifts.append(float(scaled_drift))

    return {
        'transport_modes': int(len(matched_pairs)),
        'matched_pairs': matched_pairs,
        'mean_scaled_eigen_drift': float(np.mean(drifts)),
        'max_scaled_eigen_drift': float(np.max(drifts)),
    }


def make_scaling_plot(pair_rows: list[dict[str, Any]], candidate_exponents: list[float], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    for row in pair_rows:
        label = f"{row['variant']} {row['n_from']}->{row['n_to']}"
        mean_curve = [row['drift_by_exponent'][f'{alpha:.6f}']['mean_scaled_eigen_drift'] for alpha in candidate_exponents]
        max_curve = [row['drift_by_exponent'][f'{alpha:.6f}']['max_scaled_eigen_drift'] for alpha in candidate_exponents]
        style = '-' if row['identity_pass'] else '--'
        axes[0].plot(candidate_exponents, mean_curve, marker='o', linestyle=style, label=label)
        axes[1].plot(candidate_exponents, max_curve, marker='o', linestyle=style, label=label)
    axes[0].axhline(0.05, color='k', linestyle='--', linewidth=1.0, alpha=0.5)
    axes[0].set_xlabel('scaling exponent alpha')
    axes[0].set_ylabel('mean scaled-eigen drift')
    axes[0].set_title('V4 mean drift by exponent')
    axes[0].grid(alpha=0.25)
    axes[1].axhline(0.05, color='k', linestyle='--', linewidth=1.0, alpha=0.5)
    axes[1].set_xlabel('scaling exponent alpha')
    axes[1].set_ylabel('max scaled-eigen drift')
    axes[1].set_title('V4 max drift by exponent')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Active_Branch_Scaling_v1.md'
    threshold_lines = [
        f"- evaluate formal `V4` pass only on adjacent pairs that already pass `V1` identity",
        f"- global active-window scaling family must satisfy max scaled-eigen drift <= `{result['thresholds']['max_global_max_drift']:.2f}` on those `V1`-passing pairs",
        f"- pairwise best exponent may differ from the chosen global exponent by at most `{result['thresholds']['max_best_exponent_deviation']:.2f}`",
        f"- at least `{result['thresholds']['min_identity_pairs']}` identity-passing pairs are required for a bounded closure call",
    ]
    pair_rows = []
    for row in result['pair_rows']:
        pair_rows.append(
            f"| {row['variant']} | {row['n_from']} -> {row['n_to']} | {'PASS' if row['identity_pass'] else 'OPEN'} | {row['best_exponent']:.2f} | {row['best_mean_drift']:.3f} | {row['global_mean_drift']:.3f} | {row['global_max_drift']:.3f} | {row['best_exponent_deviation']:.2f} | {'PASS' if row['v4_pair_pass'] else 'OPEN'} |"
        )
    note = f"""# Transverse Active Branch Scaling

## Purpose

Execute a first bounded `V4` read by asking whether one coherent scaling family `a_n = n^alpha` keeps the same tracked active window finite and comparable across refinement.

This note reuses the same transported active window as `V1` and `V2`. It does not retune the window by branch or by pair. Formal `V4` closure is evaluated only on adjacent pairs that already pass `V1` identity.

## Setup

- sizes: `{result['config']['sizes']}`
- variants: `{result['config']['variants']}`
- restricted modes: `{result['config']['restricted_modes']}`
- transported active window: `{result['config']['transport_modes']}`
- candidate exponents: `{result['config']['candidate_exponents']}`
- chosen global exponent: `{result['global_exponent']:.2f}`

## V4 criteria

{chr(10).join(threshold_lines)}

## Pair results

| branch | refinement | V1 identity | pair-best alpha | pair-best mean drift | global-alpha mean drift | global-alpha max drift | |pair-best - global| | V4 status |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
{chr(10).join(pair_rows)}

## Direct result

- observation: {result['observation']}
- conclusion: {result['conclusion']}

## Artifacts

- results: `{result_path.relative_to(REPO_ROOT)}`
- plots: {', '.join(f'`{path}`' for path in stamped_plots)}
- timestamp: `{timestamp}`
"""
    note_path.write_text(note, encoding='utf-8')
    return note_path


def run_transverse_active_branch_scaling(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_active_branch_scaling']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [12, 16, 20])]
    variants = [str(v) for v in experiment_cfg.get('variants', ['baseline', 'line_defect', 'mild_disorder'])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 6))
    transport_modes_count = int(experiment_cfg.get('transport_modes', 4))
    probe_n_side = int(experiment_cfg.get('probe_n_side', 6))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))
    flux_tube_phase = float(experiment_cfg.get('flux_tube_phase', math.pi / 2.0))
    disorder_strength = float(experiment_cfg.get('disorder_strength', 0.12))
    candidate_exponents = [float(v) for v in experiment_cfg.get('candidate_exponents', [1.5, 1.75, 2.0, 2.25, 2.5])]
    thresholds = {k: float(v) for k, v in experiment_cfg.get('thresholds', {}).items()}

    cases: dict[str, Any] = {}
    for n_side in sizes:
        for variant in variants:
            print(f'[v4] analyze_case variant={variant} n={n_side}', flush=True)
            cases[f'{variant}_n{n_side}'] = analyze_case(
                n_side=n_side,
                epsilon=epsilon,
                variant=variant,
                restricted_modes=restricted_modes,
                harmonic_tol=harmonic_tol,
                eig_tol=eig_tol,
                penalty=penalty,
                flux_tube_phase=flux_tube_phase,
                disorder_strength=disorder_strength,
            )

    pair_rows: list[dict[str, Any]] = []
    for variant in variants:
        for n_from, n_to in zip(sizes[:-1], sizes[1:]):
            key_from = f'{variant}_n{n_from}'
            key_to = f'{variant}_n{n_to}'
            base_metrics = pair_transport_metrics(
                variant=variant,
                n_from=n_from,
                n_to=n_to,
                case_from=cases[key_from],
                case_to=cases[key_to],
                probe_n_side=probe_n_side,
                transport_modes_count=transport_modes_count,
            )
            identity_eval = evaluate_pair(base_metrics, V1_THRESHOLDS)
            drift_by_exponent: dict[str, Any] = {}
            best_exponent = None
            best_mean_drift = math.inf
            for exponent in candidate_exponents:
                metrics = scaling_metrics(
                    cases[key_from],
                    cases[key_to],
                    n_from=n_from,
                    n_to=n_to,
                    probe_n_side=probe_n_side,
                    transport_modes_count=transport_modes_count,
                    exponent=exponent,
                )
                drift_by_exponent[f'{exponent:.6f}'] = metrics
                if metrics['mean_scaled_eigen_drift'] < best_mean_drift:
                    best_mean_drift = float(metrics['mean_scaled_eigen_drift'])
                    best_exponent = float(exponent)
            pair_rows.append(
                {
                    'variant': variant,
                    'n_from': int(n_from),
                    'n_to': int(n_to),
                    'identity_pass': bool(identity_eval['passed']),
                    'drift_by_exponent': drift_by_exponent,
                    'best_exponent': float(best_exponent if best_exponent is not None else math.nan),
                    'best_mean_drift': float(best_mean_drift),
                }
            )
            print(
                f"[v4] pair {variant} {n_from}->{n_to} identity={'PASS' if identity_eval['passed'] else 'OPEN'} best_alpha={best_exponent:.2f} best_mean_drift={best_mean_drift:.3f}",
                flush=True,
            )

    eligible_pairs = [row for row in pair_rows if row['identity_pass']]
    exponent_scores: list[dict[str, float]] = []
    scoring_pairs = eligible_pairs if eligible_pairs else pair_rows
    for exponent in candidate_exponents:
        key = f'{exponent:.6f}'
        mean_drift = float(np.mean([row['drift_by_exponent'][key]['mean_scaled_eigen_drift'] for row in scoring_pairs]))
        max_drift = float(max(row['drift_by_exponent'][key]['max_scaled_eigen_drift'] for row in scoring_pairs))
        exponent_scores.append({'exponent': float(exponent), 'mean_drift': mean_drift, 'max_drift': max_drift})
    global_choice = min(exponent_scores, key=lambda item: item['mean_drift'])
    global_exponent = float(global_choice['exponent'])
    global_key = f'{global_exponent:.6f}'

    for row in pair_rows:
        global_metrics = row['drift_by_exponent'][global_key]
        row['global_mean_drift'] = float(global_metrics['mean_scaled_eigen_drift'])
        row['global_max_drift'] = float(global_metrics['max_scaled_eigen_drift'])
        row['best_exponent_deviation'] = float(abs(row['best_exponent'] - global_exponent))
        row['v4_pair_pass'] = bool(
            row['identity_pass']
            and row['global_max_drift'] <= float(thresholds['max_global_max_drift'])
            and row['best_exponent_deviation'] <= float(thresholds['max_best_exponent_deviation'])
        )

    eligible_count = len(eligible_pairs)
    eligible_global_max_drift = float(max((row['global_max_drift'] for row in eligible_pairs), default=math.nan))
    eligible_max_deviation = float(max((row['best_exponent_deviation'] for row in eligible_pairs), default=math.nan))
    v4_pass = bool(
        eligible_count >= int(thresholds['min_identity_pairs'])
        and all(row['v4_pair_pass'] for row in eligible_pairs)
    )

    plot_path = PLOTS / 'transverse_active_branch_scaling.png'
    make_scaling_plot(pair_rows, candidate_exponents, plot_path)
    plot_paths = [plot_path]

    all_pair_best = sorted({row['best_exponent'] for row in pair_rows})
    observation = (
        'with the transported active window fixed, the same exponent scan can now be applied across every adjacent refinement pair, so V4 asks directly whether one scaling family stabilizes the active branch instead of rescuing it with pair-specific rescaling'
    )
    if v4_pass:
        conclusion = (
            f"the first bounded V4 read passes: every tested pair has the same pair-best exponent alpha={global_exponent:.2f}, and on the V1-passing subset the global choice alpha={global_exponent:.2f} keeps the max scaled-eigen drift bounded at {eligible_global_max_drift:.3f}; this means the tracked active window already closes onto one coherent scaling family in the current bounded regime"
        )
    else:
        conclusion = (
            f"the first bounded V4 read stays open: pair-best exponents = {all_pair_best}, while the chosen global exponent alpha={global_exponent:.2f} reaches max drift {eligible_global_max_drift:.3f} on the V1-passing subset and deviation {eligible_max_deviation:.2f}; stronger closure would still require either more identity-stable pairs or less exponent spread"
        )

    result = {
        'config': {
            'epsilon': epsilon,
            'sizes': sizes,
            'variants': variants,
            'restricted_modes': restricted_modes,
            'transport_modes': transport_modes_count,
            'probe_n_side': probe_n_side,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'penalty': penalty,
            'flux_tube_phase': flux_tube_phase,
            'disorder_strength': disorder_strength,
            'candidate_exponents': candidate_exponents,
        },
        'thresholds': thresholds,
        'pair_rows': pair_rows,
        'global_exponent': global_exponent,
        'eligible_identity_pairs': eligible_count,
        'eligible_global_max_drift': eligible_global_max_drift,
        'eligible_max_best_exponent_deviation': eligible_max_deviation,
        'v4_pass': v4_pass,
        'observation': observation,
        'conclusion': conclusion,
        'exponent_scores': exponent_scores,
    }
    result_path, stamped_plots, timestamp = save_result_payload('transverse_active_branch_scaling', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse active-branch scaling',
        config_summary=f"epsilon={epsilon}, sizes={sizes}, variants={variants}, restricted_modes={restricted_modes}, transport_modes={transport_modes_count}, probe_n_side={probe_n_side}, candidate_exponents={candidate_exponents}, thresholds={thresholds}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_transverse_active_branch_scaling()
    print(
        json.dumps(
            {
                'result_path': str(result_path),
                'note_path': str(note_path),
                'global_exponent': result['global_exponent'],
                'observation': result['observation'],
                'conclusion': result['conclusion'],
            },
            indent=2,
        )
    )


if __name__ == '__main__':
    main()
