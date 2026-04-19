#!/usr/bin/env python3

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from L1_stage3_common import (
    PLOTS,
    REPO_ROOT,
    VARIANT_LABELS,
    VARIANT_MARKERS,
    append_log,
    save_result_payload,
    plt,
)
from L1_transverse_band_test import analyze_case
from transverse_active_sector_transport import pair_transport_metrics

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'transverse_active_branch_identity': {
        'sizes': [12, 16, 20, 24, 28],
        'variants': ['baseline', 'puncture', 'line_defect', 'flux_tube', 'mild_disorder'],
        'restricted_modes': 6,
        'transport_modes': 6,
        'probe_n_side': 6,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
        'flux_tube_phase': math.pi / 2.0,
        'disorder_strength': 0.12,
        'thresholds': {
            'min_mean_overlap': 0.75,
            'min_min_overlap': 0.50,
            'min_mean_principal_cosine': 0.80,
            'max_max_scaled_eigen_drift': 0.05,
            'min_mean_assignment_margin': 0.10,
        },
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['transverse_active_branch_identity'] = dict(DEFAULT_CONFIG['transverse_active_branch_identity'])
    merged['transverse_active_branch_identity']['thresholds'] = dict(DEFAULT_CONFIG['transverse_active_branch_identity']['thresholds'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_active_branch_identity'})
        if isinstance(on_disk.get('transverse_active_branch_identity'), dict):
            on_disk_identity = dict(on_disk['transverse_active_branch_identity'])
            thresholds = on_disk_identity.pop('thresholds', None)
            merged['transverse_active_branch_identity'].update(on_disk_identity)
            if isinstance(thresholds, dict):
                merged['transverse_active_branch_identity']['thresholds'].update(thresholds)
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_active_branch_identity'})
        if isinstance(config.get('transverse_active_branch_identity'), dict):
            runtime_identity = dict(config['transverse_active_branch_identity'])
            thresholds = runtime_identity.pop('thresholds', None)
            merged['transverse_active_branch_identity'].update(runtime_identity)
            if isinstance(thresholds, dict):
                merged['transverse_active_branch_identity']['thresholds'].update(thresholds)
    return merged


def evaluate_pair(metrics: dict[str, Any], thresholds: dict[str, float]) -> dict[str, Any]:
    checks = {
        'mean_overlap': float(metrics['mean_overlap']) >= float(thresholds['min_mean_overlap']),
        'min_overlap': float(metrics['min_overlap']) >= float(thresholds['min_min_overlap']),
        'mean_principal_cosine': float(metrics['mean_principal_cosine']) >= float(thresholds['min_mean_principal_cosine']),
        'max_scaled_eigen_drift': float(metrics['max_scaled_eigen_drift']) <= float(thresholds['max_max_scaled_eigen_drift']),
        'mean_assignment_margin': float(metrics['mean_assignment_margin']) >= float(thresholds['min_mean_assignment_margin']),
    }
    return {
        'checks': checks,
        'passed': bool(all(checks.values())),
    }


def summarize_variant(variant: str, pair_results: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        'variant': variant,
        'variant_label': VARIANT_LABELS[variant],
        'passing_pairs': int(sum(1 for item in pair_results if item['evaluation']['passed'])),
        'total_pairs': int(len(pair_results)),
        'all_pairs_pass': bool(all(item['evaluation']['passed'] for item in pair_results)),
        'worst_mean_overlap': float(min(item['metrics']['mean_overlap'] for item in pair_results)),
        'worst_min_overlap': float(min(item['metrics']['min_overlap'] for item in pair_results)),
        'worst_mean_principal_cosine': float(min(item['metrics']['mean_principal_cosine'] for item in pair_results)),
        'worst_max_scaled_eigen_drift': float(max(item['metrics']['max_scaled_eigen_drift'] for item in pair_results)),
        'worst_mean_assignment_margin': float(min(item['metrics']['mean_assignment_margin'] for item in pair_results)),
        'mean_exact_match_fraction': float(sum(item['metrics']['exact_match_fraction'] for item in pair_results) / len(pair_results)),
    }


def make_overlap_plot(pair_results_by_variant: dict[str, list[dict[str, Any]]], sizes: list[int], path: Path) -> None:
    x_values = sizes[:-1]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    for variant, rows in pair_results_by_variant.items():
        axes[0].plot(
            x_values,
            [row['metrics']['mean_overlap'] for row in rows],
            marker=VARIANT_MARKERS[variant],
            label=VARIANT_LABELS[variant],
        )
        axes[1].plot(
            x_values,
            [row['metrics']['mean_principal_cosine'] for row in rows],
            marker=VARIANT_MARKERS[variant],
            label=VARIANT_LABELS[variant],
        )
    axes[0].axhline(0.75, color='k', linestyle='--', linewidth=1.0, alpha=0.5)
    axes[0].set_xlabel('coarse size n')
    axes[0].set_ylabel('mean matched overlap')
    axes[0].set_title('V1 overlap criterion')
    axes[0].set_ylim(0.0, 1.02)
    axes[0].grid(alpha=0.25)
    axes[1].axhline(0.80, color='k', linestyle='--', linewidth=1.0, alpha=0.5)
    axes[1].set_xlabel('coarse size n')
    axes[1].set_ylabel('mean principal cosine')
    axes[1].set_title('V1 subspace-alignment criterion')
    axes[1].set_ylim(0.0, 1.02)
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_identity_plot(pair_results_by_variant: dict[str, list[dict[str, Any]]], sizes: list[int], path: Path) -> None:
    x_values = sizes[:-1]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    for variant, rows in pair_results_by_variant.items():
        axes[0].plot(
            x_values,
            [row['metrics']['mean_assignment_margin'] for row in rows],
            marker=VARIANT_MARKERS[variant],
            label=VARIANT_LABELS[variant],
        )
        axes[1].plot(
            x_values,
            [row['metrics']['max_scaled_eigen_drift'] for row in rows],
            marker=VARIANT_MARKERS[variant],
            label=VARIANT_LABELS[variant],
        )
    axes[0].axhline(0.10, color='k', linestyle='--', linewidth=1.0, alpha=0.5)
    axes[0].set_xlabel('coarse size n')
    axes[0].set_ylabel('mean assignment margin')
    axes[0].set_title('V1 relabeling-coherence criterion')
    axes[0].set_ylim(0.0, 1.02)
    axes[0].grid(alpha=0.25)
    axes[1].axhline(0.05, color='k', linestyle='--', linewidth=1.0, alpha=0.5)
    axes[1].set_xlabel('coarse size n')
    axes[1].set_ylabel('max scaled-eigen drift')
    axes[1].set_title('V1 eigen-drift criterion')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Active_Branch_Identity_v1.md'
    threshold_lines = [
        f"- mean matched overlap >= `{result['thresholds']['min_mean_overlap']:.2f}`",
        f"- minimum matched overlap >= `{result['thresholds']['min_min_overlap']:.2f}`",
        f"- mean principal cosine >= `{result['thresholds']['min_mean_principal_cosine']:.2f}`",
        f"- max scaled-eigen drift <= `{result['thresholds']['max_max_scaled_eigen_drift']:.2f}`",
        f"- mean assignment margin >= `{result['thresholds']['min_mean_assignment_margin']:.2f}`",
    ]
    pair_rows = []
    for item in result['pair_results']:
        metrics = item['metrics']
        pair_rows.append(
            f"| {item['variant']} | {item['n_from']} -> {item['n_to']} | {metrics['mean_overlap']:.3f} | {metrics['min_overlap']:.3f} | {metrics['mean_principal_cosine']:.3f} | {metrics['max_scaled_eigen_drift']:.3f} | {metrics['mean_assignment_margin']:.3f} | {metrics['exact_match_fraction']:.3f} | {'PASS' if item['evaluation']['passed'] else 'FAIL'} |"
        )
    variant_rows = []
    for summary in result['variant_summary']:
        variant_rows.append(
            f"| {summary['variant']} | {summary['passing_pairs']}/{summary['total_pairs']} | {summary['worst_mean_overlap']:.3f} | {summary['worst_mean_principal_cosine']:.3f} | {summary['worst_max_scaled_eigen_drift']:.3f} | {summary['worst_mean_assignment_margin']:.3f} | {summary['mean_exact_match_fraction']:.3f} | {'PASS' if summary['all_pairs_pass'] else 'OPEN'} |"
        )
    note = f"""# Transverse Active Branch Identity

## Purpose

Execute `V1` of the post-foundations validation line by asking whether the same low active coexact family survives as the same physical object across refinement.

The active-sector transport map `J_n` is kept fixed. This note does not test threshold stabilization yet; it tests identity on the transported branch itself.

## Setup

- sizes: `{result['config']['sizes']}`
- variants: `{result['config']['variants']}`
- restricted modes: `{result['config']['restricted_modes']}`
- transported low-window modes: `{result['config']['transport_modes']}`
- probe lattice side: `{result['config']['probe_n_side']}`
- mild-disorder strength: `{result['config']['disorder_strength']}`

## V1 criteria

The present bounded `V1` diagnostic declares an adjacent refinement pair `PASS` only if all of the following hold:

{chr(10).join(threshold_lines)}

Assignment margin is the chosen-match overlap minus the strongest competing overlap on the same row / column pair, clipped below at zero. It is used here as the bounded operational proxy for "no arbitrary relabeling."

## Pair results

| branch | refinement | mean overlap | min overlap | mean principal cosine | max scaled-eigen drift | mean assignment margin | exact-match fraction | status |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
{chr(10).join(pair_rows)}

## Variant summary

| branch | passing pairs | worst mean overlap | worst mean principal cosine | worst max drift | worst mean margin | mean exact-match fraction | V1 status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
{chr(10).join(variant_rows)}

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


def run_transverse_active_branch_identity(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_active_branch_identity']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [12, 16, 20, 24, 28])]
    variants = [str(v) for v in experiment_cfg.get('variants', ['baseline', 'puncture', 'line_defect', 'flux_tube', 'mild_disorder'])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 6))
    transport_modes_count = int(experiment_cfg.get('transport_modes', 6))
    probe_n_side = int(experiment_cfg.get('probe_n_side', 6))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))
    flux_tube_phase = float(experiment_cfg.get('flux_tube_phase', math.pi / 2.0))
    disorder_strength = float(experiment_cfg.get('disorder_strength', 0.12))
    thresholds = {k: float(v) for k, v in experiment_cfg.get('thresholds', {}).items()}

    cases: dict[str, Any] = {}
    for n_side in sizes:
        for variant in variants:
            print(f'[v1] analyze_case variant={variant} n={n_side}', flush=True)
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

    pair_results: list[dict[str, Any]] = []
    pair_results_by_variant: dict[str, list[dict[str, Any]]] = {variant: [] for variant in variants}
    for variant in variants:
        for n_from, n_to in zip(sizes[:-1], sizes[1:]):
            print(f'[v1] transport_metrics variant={variant} {n_from}->{n_to}', flush=True)
            metrics = pair_transport_metrics(
                variant=variant,
                n_from=n_from,
                n_to=n_to,
                case_from=cases[f'{variant}_n{n_from}'],
                case_to=cases[f'{variant}_n{n_to}'],
                probe_n_side=probe_n_side,
                transport_modes_count=transport_modes_count,
            )
            evaluation = evaluate_pair(metrics, thresholds)
            item = {
                'variant': variant,
                'n_from': int(n_from),
                'n_to': int(n_to),
                'metrics': metrics,
                'evaluation': evaluation,
            }
            pair_results.append(item)
            pair_results_by_variant[variant].append(item)
            print(
                '[v1] evaluated '
                f"variant={variant} {n_from}->{n_to} "
                f"mean_overlap={metrics['mean_overlap']:.3f} "
                f"mean_principal_cosine={metrics['mean_principal_cosine']:.3f} "
                f"mean_assignment_margin={metrics['mean_assignment_margin']:.3f} "
                f"status={'PASS' if evaluation['passed'] else 'FAIL'}",
                flush=True,
            )

    variant_summary = [summarize_variant(variant, pair_results_by_variant[variant]) for variant in variants]

    overlap_plot = PLOTS / 'transverse_active_branch_identity_overlap.png'
    identity_plot = PLOTS / 'transverse_active_branch_identity_margin.png'
    make_overlap_plot(pair_results_by_variant, sizes, overlap_plot)
    make_identity_plot(pair_results_by_variant, sizes, identity_plot)
    plot_paths = [overlap_plot, identity_plot]

    passing_variants = [summary['variant'] for summary in variant_summary if summary['all_pairs_pass']]
    failing_variants = [summary['variant'] for summary in variant_summary if not summary['all_pairs_pass']]
    weakest_variant = min(variant_summary, key=lambda item: item['worst_mean_overlap'])
    strongest_variant = max(variant_summary, key=lambda item: item['passing_pairs'] / max(item['total_pairs'], 1))
    observation = (
        'with the active coexact restriction and probe-space transport map held fixed, branch identity can now be tested directly against explicit V1 thresholds instead of being inferred from raw qualitative stability language'
    )
    if len(failing_variants) == 0:
        conclusion = (
            f"all tested variants pass every adjacent refinement pair under the present V1 thresholds, with {strongest_variant['variant']} giving the strongest aggregate identity retention; in this bounded window, V1 would be closed"
        )
    else:
        passing_text = ', '.join(passing_variants) if passing_variants else 'none'
        failing_text = ', '.join(failing_variants)
        conclusion = (
            f"V1 is not yet closed on the full tested family: passing variants = {passing_text}, open variants = {failing_text}; the weakest case is {weakest_variant['variant']} with worst mean overlap {weakest_variant['worst_mean_overlap']:.3f} and worst mean assignment margin {weakest_variant['worst_mean_assignment_margin']:.3f}, so the active branch is now measurable but still not uniformly stable as the same transported physical object across refinement"
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
        },
        'thresholds': thresholds,
        'pair_results': pair_results,
        'variant_summary': variant_summary,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('transverse_active_branch_identity', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse active-branch identity',
        config_summary=(
            f"epsilon={epsilon}, sizes={sizes}, variants={variants}, restricted_modes={restricted_modes}, "
            f"transport_modes={transport_modes_count}, probe_n_side={probe_n_side}, disorder_strength={disorder_strength}, thresholds={thresholds}"
        ),
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_transverse_active_branch_identity()
    print(json.dumps({'result_path': str(result_path), 'note_path': str(note_path), 'observation': result['observation'], 'conclusion': result['conclusion']}, indent=2))


if __name__ == '__main__':
    main()
