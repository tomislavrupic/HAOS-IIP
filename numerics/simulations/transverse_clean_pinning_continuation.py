#!/usr/bin/env python3

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from L1_stage3_common import PLOTS, REPO_ROOT, append_log, save_result_payload, plt
from L1_transverse_band_test import analyze_case
from transverse_clean_subspace_transport import best_window_metrics, evaluate_pair, window_metrics

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'transverse_clean_pinning_continuation': {
        'sizes': [12, 16, 20],
        'pinning_strengths': [0.0, 0.002, 0.005, 0.01, 0.015],
        'restricted_modes': 8,
        'window_size': 4,
        'scan_depth': 8,
        'probe_n_side': 6,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
        'flux_tube_phase': math.pi / 2.0,
        'thresholds': {
            'min_best_projector_affinity': 0.80,
            'min_best_mean_principal_cosine': 0.80,
            'max_best_scaled_window_drift': 0.05,
            'min_affinity_gain': 0.20,
        },
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['transverse_clean_pinning_continuation'] = dict(DEFAULT_CONFIG['transverse_clean_pinning_continuation'])
    merged['transverse_clean_pinning_continuation']['thresholds'] = dict(DEFAULT_CONFIG['transverse_clean_pinning_continuation']['thresholds'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_clean_pinning_continuation'})
        if isinstance(on_disk.get('transverse_clean_pinning_continuation'), dict):
            block = dict(on_disk['transverse_clean_pinning_continuation'])
            thresholds = block.pop('thresholds', None)
            merged['transverse_clean_pinning_continuation'].update(block)
            if isinstance(thresholds, dict):
                merged['transverse_clean_pinning_continuation']['thresholds'].update(thresholds)
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_clean_pinning_continuation'})
        if isinstance(config.get('transverse_clean_pinning_continuation'), dict):
            block = dict(config['transverse_clean_pinning_continuation'])
            thresholds = block.pop('thresholds', None)
            merged['transverse_clean_pinning_continuation'].update(block)
            if isinstance(thresholds, dict):
                merged['transverse_clean_pinning_continuation']['thresholds'].update(thresholds)
    return merged


def variant_for_strength(strength: float) -> str:
    return 'baseline' if abs(strength) <= 1.0e-15 else 'mild_disorder'


def summarize_strength(strength: float, pair_rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        'pinning_strength': float(strength),
        'resolved_pairs': int(sum(1 for row in pair_rows if row['evaluation']['resolved'])),
        'total_pairs': int(len(pair_rows)),
        'all_pairs_resolved': bool(all(row['evaluation']['resolved'] for row in pair_rows)),
        'all_pairs_stable_fixed': bool(all(row['evaluation']['status'] == 'stable_fixed' for row in pair_rows)),
        'worst_fixed_affinity': float(min(row['fixed_window']['projector_affinity'] for row in pair_rows)),
        'worst_best_affinity': float(min(row['best_window']['projector_affinity'] for row in pair_rows)),
        'worst_best_mean_principal_cosine': float(min(row['best_window']['mean_principal_cosine'] for row in pair_rows)),
        'worst_best_scaled_window_drift': float(max(row['best_window']['scaled_window_drift'] for row in pair_rows)),
        'largest_affinity_gain': float(max(row['evaluation']['affinity_gain'] for row in pair_rows)),
    }


def make_affinity_plot(summary_rows: list[dict[str, Any]], path: Path) -> None:
    strengths = [row['pinning_strength'] for row in summary_rows]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    axes[0].plot(strengths, [row['worst_fixed_affinity'] for row in summary_rows], marker='o', label='worst fixed affinity')
    axes[0].plot(strengths, [row['worst_best_affinity'] for row in summary_rows], marker='s', label='worst best affinity')
    axes[0].axhline(0.80, color='k', linestyle='--', linewidth=1.0, alpha=0.5)
    axes[0].set_xlabel('pinning strength')
    axes[0].set_ylabel('subspace affinity')
    axes[0].set_title('Worst clean-window affinity under pinning')
    axes[0].set_ylim(0.0, 1.02)
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=8)

    axes[1].plot(strengths, [row['worst_best_mean_principal_cosine'] for row in summary_rows], marker='o', label='worst best cosine')
    axes[1].plot(strengths, [row['worst_best_scaled_window_drift'] for row in summary_rows], marker='s', label='worst best drift')
    axes[1].axhline(0.80, color='k', linestyle='--', linewidth=1.0, alpha=0.5)
    axes[1].axhline(0.05, color='k', linestyle=':', linewidth=1.0, alpha=0.5)
    axes[1].set_xlabel('pinning strength')
    axes[1].set_ylabel('cosine / drift')
    axes[1].set_title('Best-window alignment and drift')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_resolution_plot(summary_rows: list[dict[str, Any]], path: Path) -> None:
    strengths = [row['pinning_strength'] for row in summary_rows]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(strengths, [row['resolved_pairs'] for row in summary_rows], marker='o', label='resolved pairs')
    ax.plot(strengths, [row['largest_affinity_gain'] for row in summary_rows], marker='s', label='largest affinity gain')
    ax.set_xlabel('pinning strength')
    ax.set_ylabel('resolution count / gain')
    ax.set_title('Clean-branch resolution under infinitesimal pinning')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Clean_Pinning_Continuation_v1.md'
    rows = []
    for item in result['summary_rows']:
        rows.append(
            f"| {item['pinning_strength']:.3f} | {item['resolved_pairs']}/{item['total_pairs']} | {item['worst_fixed_affinity']:.3f} | {item['worst_best_affinity']:.3f} | {item['worst_best_mean_principal_cosine']:.3f} | {item['worst_best_scaled_window_drift']:.3f} | {item['largest_affinity_gain']:.3f} | {'RESOLVED' if item['all_pairs_resolved'] else 'OPEN'} |"
        )
    note = f"""# Transverse Clean Pinning Continuation

## Purpose

Execute `CB2` in bounded form by testing whether infinitesimal symmetry pinning resolves the clean-baseline identity question when measured with the same windowed subspace metric introduced in `CB1-lite`.

This scan keeps the same active coexact sector, the same probe-space transport rule, and the same low-window scan. Only the smooth pinning strength changes.

## Setup

- sizes: `{result['config']['sizes']}`
- pinning strengths: `{result['config']['pinning_strengths']}`
- restricted modes analyzed: `{result['config']['restricted_modes']}`
- window size: `{result['config']['window_size']}`
- scan depth: `{result['config']['scan_depth']}`
- probe lattice side: `{result['config']['probe_n_side']}`

## Strength summary

| pinning strength | resolved pairs | worst fixed affinity | worst best affinity | worst best cosine | worst best drift | largest affinity gain | CB2 status |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
{chr(10).join(rows)}

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


def run_transverse_clean_pinning_continuation(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_clean_pinning_continuation']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [12, 16, 20])]
    pinning_strengths = [float(v) for v in experiment_cfg.get('pinning_strengths', [0.0, 0.002, 0.005, 0.01, 0.015])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 8))
    window_size = int(experiment_cfg.get('window_size', 4))
    scan_depth = int(experiment_cfg.get('scan_depth', 8))
    probe_n_side = int(experiment_cfg.get('probe_n_side', 6))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))
    flux_tube_phase = float(experiment_cfg.get('flux_tube_phase', math.pi / 2.0))
    thresholds = {k: float(v) for k, v in experiment_cfg.get('thresholds', {}).items()}

    continuation_rows: list[dict[str, Any]] = []
    detailed_rows: list[dict[str, Any]] = []
    for strength in pinning_strengths:
        variant = variant_for_strength(strength)
        print(f'[cb2] starting strength={strength:.3f} variant={variant}', flush=True)
        cases: dict[int, dict[str, Any]] = {}
        for n_side in sizes:
            print(f'[cb2] analyze_case strength={strength:.3f} n={n_side}', flush=True)
            kwargs: dict[str, Any] = {
                'n_side': n_side,
                'epsilon': epsilon,
                'variant': variant,
                'restricted_modes': restricted_modes,
                'harmonic_tol': harmonic_tol,
                'eig_tol': eig_tol,
                'penalty': penalty,
                'flux_tube_phase': flux_tube_phase,
            }
            if variant == 'mild_disorder':
                kwargs['disorder_strength'] = strength
            cases[n_side] = analyze_case(**kwargs)

        pair_rows: list[dict[str, Any]] = []
        for n_from, n_to in zip(sizes[:-1], sizes[1:]):
            fixed_metrics = window_metrics(cases[n_from], cases[n_to], n_from, n_to, 0, 0, window_size, probe_n_side)
            best_metrics = best_window_metrics(cases[n_from], cases[n_to], n_from, n_to, scan_depth, window_size, probe_n_side)
            evaluation = evaluate_pair(best_metrics, fixed_metrics, thresholds)
            pair_row = {
                'pinning_strength': strength,
                'variant': variant,
                'n_from': int(n_from),
                'n_to': int(n_to),
                'fixed_window': fixed_metrics,
                'best_window': best_metrics,
                'evaluation': evaluation,
            }
            pair_rows.append(pair_row)
            detailed_rows.append(pair_row)
            print(
                f"[cb2] strength={strength:.3f} {n_from}->{n_to} "
                f"fixed_affinity={fixed_metrics['projector_affinity']:.3f} "
                f"best_affinity={best_metrics['projector_affinity']:.3f} "
                f"status={evaluation['status']}",
                flush=True,
            )

        continuation_rows.append(summarize_strength(strength, pair_rows))

    affinity_plot = PLOTS / 'transverse_clean_pinning_continuation_affinity.png'
    resolution_plot = PLOTS / 'transverse_clean_pinning_continuation_resolution.png'
    make_affinity_plot(continuation_rows, affinity_plot)
    make_resolution_plot(continuation_rows, resolution_plot)
    plot_paths = [affinity_plot, resolution_plot]

    baseline_row = next(row for row in continuation_rows if abs(row['pinning_strength']) <= 1.0e-15)
    resolved_row = next((row for row in continuation_rows if row['all_pairs_resolved']), None)
    observation = (
        'the clean-baseline identity question can be sharpened by introducing an infinitesimal symmetry pinning and measuring the same low-window subspace metric as the pinning is removed'
    )
    if resolved_row is None:
        conclusion = (
            f"even with infinitesimal pinning up to strength {continuation_rows[-1]['pinning_strength']:.3f}, the clean family never fully resolves under the present windowed subspace metric: "
            f"the baseline starts at worst best affinity {baseline_row['worst_best_affinity']:.3f} and the strongest tested pinning only reaches {max(row['worst_best_affinity'] for row in continuation_rows):.3f}; "
            f"this points to a deeper clean-baseline issue than simple infinitesimal symmetry lifting"
        )
    else:
        conclusion = (
            f"the clean family resolves under infinitesimal pinning: the baseline starts at worst best affinity {baseline_row['worst_best_affinity']:.3f}, and all adjacent refinement pairs first resolve at pinning strength {resolved_row['pinning_strength']:.3f}; "
            f"this strongly supports a degenerate clean-family interpretation rather than branch failure"
        )

    result = {
        'config': {
            'sizes': sizes,
            'pinning_strengths': pinning_strengths,
            'restricted_modes': restricted_modes,
            'window_size': window_size,
            'scan_depth': scan_depth,
            'probe_n_side': probe_n_side,
            'penalty': penalty,
        },
        'thresholds': thresholds,
        'summary_rows': continuation_rows,
        'pair_rows': detailed_rows,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('transverse_clean_pinning_continuation', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse clean pinning continuation',
        config_summary=(
            f"epsilon={epsilon}, sizes={sizes}, pinning_strengths={pinning_strengths}, restricted_modes={restricted_modes}, "
            f"window_size={window_size}, scan_depth={scan_depth}, probe_n_side={probe_n_side}, thresholds={thresholds}"
        ),
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


if __name__ == '__main__':
    result, result_path, _, _, note_path = run_transverse_clean_pinning_continuation()
    print(json.dumps(result, indent=2))
    print(f'results: {result_path.relative_to(REPO_ROOT)}')
    print(f'note: {note_path.relative_to(REPO_ROOT)}')
