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

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'transverse_active_branch_purity': {
        'sizes': [12, 16, 20],
        'variants': ['baseline', 'line_defect', 'mild_disorder'],
        'restricted_modes': 6,
        'purity_modes': 4,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
        'flux_tube_phase': math.pi / 2.0,
        'disorder_strength': 0.12,
        'thresholds': {
            'max_exact_fraction': 1e-8,
            'max_harmonic_fraction': 1e-8,
            'min_coexact_fraction': 0.999999,
        },
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['transverse_active_branch_purity'] = dict(DEFAULT_CONFIG['transverse_active_branch_purity'])
    merged['transverse_active_branch_purity']['thresholds'] = dict(DEFAULT_CONFIG['transverse_active_branch_purity']['thresholds'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_active_branch_purity'})
        if isinstance(on_disk.get('transverse_active_branch_purity'), dict):
            block = dict(on_disk['transverse_active_branch_purity'])
            thresholds = block.pop('thresholds', None)
            merged['transverse_active_branch_purity'].update(block)
            if isinstance(thresholds, dict):
                merged['transverse_active_branch_purity']['thresholds'].update(thresholds)
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_active_branch_purity'})
        if isinstance(config.get('transverse_active_branch_purity'), dict):
            block = dict(config['transverse_active_branch_purity'])
            thresholds = block.pop('thresholds', None)
            merged['transverse_active_branch_purity'].update(block)
            if isinstance(thresholds, dict):
                merged['transverse_active_branch_purity']['thresholds'].update(thresholds)
    return merged


def evaluate_case(records: list[dict[str, Any]], thresholds: dict[str, float], purity_modes: int) -> dict[str, Any]:
    window = records[:purity_modes]
    worst_exact = max(float(record['exact_fraction']) for record in window)
    worst_harmonic = max(float(record['harmonic_fraction']) for record in window)
    worst_coexact = min(float(record['coexact_fraction']) for record in window)
    checks = {
        'exact_fraction': worst_exact <= float(thresholds['max_exact_fraction']),
        'harmonic_fraction': worst_harmonic <= float(thresholds['max_harmonic_fraction']),
        'coexact_fraction': worst_coexact >= float(thresholds['min_coexact_fraction']),
    }
    return {
        'worst_exact_fraction': worst_exact,
        'worst_harmonic_fraction': worst_harmonic,
        'worst_coexact_fraction': worst_coexact,
        'checks': checks,
        'passed': bool(all(checks.values())),
    }


def summarize_variant(variant: str, rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        'variant': variant,
        'variant_label': VARIANT_LABELS[variant],
        'passing_sizes': int(sum(1 for row in rows if row['evaluation']['passed'])),
        'total_sizes': int(len(rows)),
        'all_sizes_pass': bool(all(row['evaluation']['passed'] for row in rows)),
        'worst_exact_fraction': float(max(row['evaluation']['worst_exact_fraction'] for row in rows)),
        'worst_harmonic_fraction': float(max(row['evaluation']['worst_harmonic_fraction'] for row in rows)),
        'worst_coexact_fraction': float(min(row['evaluation']['worst_coexact_fraction'] for row in rows)),
    }


def make_purity_plot(case_rows_by_variant: dict[str, list[dict[str, Any]]], sizes: list[int], thresholds: dict[str, float], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    for variant, rows in case_rows_by_variant.items():
        axes[0].plot(
            sizes,
            [row['evaluation']['worst_harmonic_fraction'] for row in rows],
            marker=VARIANT_MARKERS[variant],
            label=VARIANT_LABELS[variant],
        )
        axes[1].plot(
            sizes,
            [row['evaluation']['worst_coexact_fraction'] for row in rows],
            marker=VARIANT_MARKERS[variant],
            label=VARIANT_LABELS[variant],
        )
    axes[0].axhline(float(thresholds['max_harmonic_fraction']), color='k', linestyle='--', linewidth=1.0, alpha=0.5)
    axes[0].set_xlabel('size n')
    axes[0].set_ylabel('worst harmonic fraction')
    axes[0].set_title('V3 harmonic-leakage bound')
    axes[0].grid(alpha=0.25)
    axes[1].axhline(float(thresholds['min_coexact_fraction']), color='k', linestyle='--', linewidth=1.0, alpha=0.5)
    axes[1].set_xlabel('size n')
    axes[1].set_ylabel('worst coexact fraction')
    axes[1].set_title('V3 coexact-purity bound')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Active_Branch_Purity_v1.md'
    threshold_lines = [
        f"- worst exact fraction <= `{result['thresholds']['max_exact_fraction']:.1e}`",
        f"- worst harmonic fraction <= `{result['thresholds']['max_harmonic_fraction']:.1e}`",
        f"- worst coexact fraction >= `{result['thresholds']['min_coexact_fraction']:.6f}`",
    ]
    case_rows = []
    for row in result['case_rows']:
        ev = row['evaluation']
        case_rows.append(
            f"| {row['variant']} | {row['n_side']} | {ev['worst_exact_fraction']:.3e} | {ev['worst_harmonic_fraction']:.3e} | {ev['worst_coexact_fraction']:.9f} | {'PASS' if ev['passed'] else 'FAIL'} |"
        )
    variant_rows = []
    for row in result['variant_summary']:
        variant_rows.append(
            f"| {row['variant']} | {row['passing_sizes']}/{row['total_sizes']} | {row['worst_exact_fraction']:.3e} | {row['worst_harmonic_fraction']:.3e} | {row['worst_coexact_fraction']:.9f} | {'PASS' if row['all_sizes_pass'] else 'OPEN'} |"
        )
    note = f"""# Transverse Active Branch Purity

## Purpose

Execute a first bounded `V3` read by testing whether the active restricted branch remains physically coexact under refinement rather than drifting back into harmonic or mixed sectors.

This first freeze uses the restricted active window itself:

- sizes: `{result['config']['sizes']}`
- variants: `{result['config']['variants']}`
- restricted modes: `{result['config']['restricted_modes']}`
- purity window: first `{result['config']['purity_modes']}` restricted modes

## V3 criteria

For this first bounded `V3` freeze, a size slice counts as `PASS` only if all of the following hold on the tested low restricted window:

{chr(10).join(threshold_lines)}

This is an internal active-branch purity test. It does not yet test transported purity against a separately transported harmonic subspace.

## Case results

| branch | size | worst exact fraction | worst harmonic fraction | worst coexact fraction | status |
| --- | ---: | ---: | ---: | ---: | --- |
{chr(10).join(case_rows)}

## Variant summary

| branch | passing sizes | worst exact fraction | worst harmonic fraction | worst coexact fraction | V3 status |
| --- | ---: | ---: | ---: | ---: | --- |
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


def run_transverse_active_branch_purity(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_active_branch_purity']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [12, 16, 20])]
    variants = [str(v) for v in experiment_cfg.get('variants', ['baseline', 'line_defect', 'mild_disorder'])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 6))
    purity_modes = int(experiment_cfg.get('purity_modes', 4))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))
    flux_tube_phase = float(experiment_cfg.get('flux_tube_phase', math.pi / 2.0))
    disorder_strength = float(experiment_cfg.get('disorder_strength', 0.12))
    thresholds = {k: float(v) for k, v in experiment_cfg.get('thresholds', {}).items()}

    case_rows: list[dict[str, Any]] = []
    case_rows_by_variant: dict[str, list[dict[str, Any]]] = {variant: [] for variant in variants}
    for n_side in sizes:
        for variant in variants:
            print(f'[v3] analyze_case variant={variant} n={n_side}', flush=True)
            case = analyze_case(
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
            evaluation = evaluate_case(case['restricted_transverse_modes'], thresholds, purity_modes)
            row = {
                'variant': variant,
                'n_side': int(n_side),
                'evaluation': evaluation,
            }
            case_rows.append(row)
            case_rows_by_variant[variant].append(row)
            print(
                '[v3] evaluated '
                f"variant={variant} n={n_side} "
                f"harmonic={evaluation['worst_harmonic_fraction']:.3e} "
                f"coexact={evaluation['worst_coexact_fraction']:.9f} "
                f"status={'PASS' if evaluation['passed'] else 'FAIL'}",
                flush=True,
            )

    variant_summary = [summarize_variant(variant, case_rows_by_variant[variant]) for variant in variants]
    purity_plot = PLOTS / 'transverse_active_branch_purity.png'
    make_purity_plot(case_rows_by_variant, sizes, thresholds, purity_plot)
    plot_paths = [purity_plot]

    failing_variants = [item['variant'] for item in variant_summary if not item['all_sizes_pass']]
    worst_variant = max(variant_summary, key=lambda item: item['worst_harmonic_fraction'])
    observation = (
        'the restricted active branch now carries explicit exact / harmonic / coexact fractions, so harmonic recontamination can be tested directly on the low active window instead of being assumed away'
    )
    if failing_variants:
        conclusion = (
            f"V3 is still open on the bounded trio: failing variants = {', '.join(failing_variants)}; the worst harmonic leakage occurs on {worst_variant['variant']} at {worst_variant['worst_harmonic_fraction']:.3e}, so the restricted branch is not yet uniformly clean under the present purity thresholds"
        )
    else:
        conclusion = (
            f"all tested variants pass the first bounded V3 purity thresholds, and the worst harmonic leakage remains {worst_variant['worst_harmonic_fraction']:.3e}; in this first bounded read, the low restricted active window remains physically coexact to leading order under refinement"
        )

    result = {
        'config': {
            'epsilon': epsilon,
            'sizes': sizes,
            'variants': variants,
            'restricted_modes': restricted_modes,
            'purity_modes': purity_modes,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'penalty': penalty,
            'flux_tube_phase': flux_tube_phase,
            'disorder_strength': disorder_strength,
        },
        'thresholds': thresholds,
        'case_rows': case_rows,
        'variant_summary': variant_summary,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('transverse_active_branch_purity', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse active-branch purity',
        config_summary=(
            f"epsilon={epsilon}, sizes={sizes}, variants={variants}, restricted_modes={restricted_modes}, "
            f"purity_modes={purity_modes}, disorder_strength={disorder_strength}, thresholds={thresholds}"
        ),
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_transverse_active_branch_purity()
    print(json.dumps({'result_path': str(result_path), 'note_path': str(note_path), 'observation': result['observation'], 'conclusion': result['conclusion']}, indent=2))


if __name__ == '__main__':
    main()
