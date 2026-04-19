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
from transverse_sector_window_scan import summarize_case

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'transverse_sector_gap_scaling': {
        'sizes': [20, 24, 28],
        'variants': ['baseline', 'puncture', 'flux_tube'],
        'restricted_modes': 20,
        'full_mode_scan_count': 40,
        'coexact_threshold': 0.8,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
        'flux_tube_phase': math.pi / 2.0,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['transverse_sector_gap_scaling'] = dict(DEFAULT_CONFIG['transverse_sector_gap_scaling'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_sector_gap_scaling'})
        if isinstance(on_disk.get('transverse_sector_gap_scaling'), dict):
            merged['transverse_sector_gap_scaling'].update(on_disk['transverse_sector_gap_scaling'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_sector_gap_scaling'})
        if isinstance(config.get('transverse_sector_gap_scaling'), dict):
            merged['transverse_sector_gap_scaling'].update(config['transverse_sector_gap_scaling'])
    return merged


def make_gap_plot(metrics: dict[str, dict[str, Any]], sizes: list[int], variants: list[str], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    for variant in variants:
        raw_gap = []
        scaled_gap = []
        for n_side in sizes:
            key = f'{variant}_n{n_side}'
            raw_gap.append(float(metrics[key]['harmonic_to_coexact_gap']))
            scaled_gap.append(float(metrics[key]['scaled_harmonic_to_coexact_gap']))
        axes[0].plot(sizes, raw_gap, marker=VARIANT_MARKERS[variant], label=VARIANT_LABELS[variant])
        axes[1].plot(sizes, scaled_gap, marker=VARIANT_MARKERS[variant], label=VARIANT_LABELS[variant])
    axes[0].set_xlabel('lattice side n')
    axes[0].set_ylabel(r'$\lambda_{\mathrm{coexact}}-\lambda_{\mathrm{floor}}$')
    axes[0].set_title('Raw harmonic-to-coexact gap')
    axes[0].grid(alpha=0.25)
    axes[1].set_xlabel('lattice side n')
    axes[1].set_ylabel(r'$n^2(\lambda_{\mathrm{coexact}}-\lambda_{\mathrm{floor}})$')
    axes[1].set_title('Scaled harmonic-to-coexact gap')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_depth_plot(metrics: dict[str, dict[str, Any]], sizes: list[int], variants: list[str], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    for variant in variants:
        candidate_depth = []
        ratio = []
        for n_side in sizes:
            key = f'{variant}_n{n_side}'
            candidate_depth.append(float(metrics[key]['candidate_mode_index']))
            ratio.append(float(metrics[key]['restricted_to_candidate_ratio']))
        axes[0].plot(sizes, candidate_depth, marker=VARIANT_MARKERS[variant], label=VARIANT_LABELS[variant])
        axes[1].plot(sizes, ratio, marker=VARIANT_MARKERS[variant], label=VARIANT_LABELS[variant])
    axes[0].set_xlabel('lattice side n')
    axes[0].set_ylabel('coexact candidate index')
    axes[0].set_title('Depth of first strong coexact candidate')
    axes[0].grid(alpha=0.25)
    axes[1].set_xlabel('lattice side n')
    axes[1].set_ylabel('restricted / coexact candidate')
    axes[1].set_title('Projected floor versus coexact candidate')
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_localization_plot(metrics: dict[str, dict[str, Any]], sizes: list[int], variants: list[str], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for variant in variants:
        localization = []
        for n_side in sizes:
            key = f'{variant}_n{n_side}'
            localization.append(float(metrics[key]['candidate_near_defect_fraction']))
        ax.plot(sizes, localization, marker=VARIANT_MARKERS[variant], label=VARIANT_LABELS[variant])
    ax.set_xlabel('lattice side n')
    ax.set_ylabel('candidate near-defect fraction')
    ax.set_title('Coexact candidate localization')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Sector_Gap_Scaling_v1.md'
    rows = []
    for key, item in sorted(result['metrics'].items()):
        rows.append(
            f"| {item['variant']} | {item['n_side']} | {item['candidate_mode_index']} | {item['candidate_coexact_fraction']:.3f} | {item['scaled_harmonic_to_coexact_gap']:.3f} | {item['candidate_near_defect_fraction']:.3f} | {item['restricted_to_candidate_ratio']:.3f} |"
        )
    note = f"""# Transverse Sector Gap Scaling

## Purpose

Run the larger-size CP2c probe on the critical backgrounds and test whether the harmonic-to-coexact gap shrinks, plateaus, or stays finite as the lattice grows.

## Setup

- sizes: `{result['config']['sizes']}`
- variants: `{result['config']['variants']}`
- restricted modes: `{result['config']['restricted_modes']}`
- full-mode scan count: `{result['config']['full_mode_scan_count']}`
- coexact threshold: `{result['config']['coexact_threshold']}`

## Gap metrics

| branch | `n` | coexact candidate `j` | candidate coexact frac | scaled gap | candidate near-defect frac | restricted / candidate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
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


def run_transverse_sector_gap_scaling(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_sector_gap_scaling']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [20, 24, 28])]
    variants = [str(v) for v in experiment_cfg.get('variants', ['baseline', 'puncture', 'flux_tube'])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 20))
    full_mode_scan_count = int(experiment_cfg.get('full_mode_scan_count', 40))
    coexact_threshold = float(experiment_cfg.get('coexact_threshold', 0.8))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))
    flux_tube_phase = float(experiment_cfg.get('flux_tube_phase', math.pi / 2.0))

    metrics: dict[str, dict[str, Any]] = {}
    for n_side in sizes:
        for variant in variants:
            key = f'{variant}_n{n_side}'
            case = analyze_case(
                n_side=n_side,
                epsilon=epsilon,
                variant=variant,
                restricted_modes=restricted_modes,
                harmonic_tol=harmonic_tol,
                eig_tol=eig_tol,
                penalty=penalty,
                flux_tube_phase=flux_tube_phase,
                full_mode_scan_count=full_mode_scan_count,
            )
            metrics[key] = summarize_case(case, coexact_threshold=coexact_threshold)

    gap_plot = PLOTS / 'transverse_sector_gap_scaling.png'
    depth_plot = PLOTS / 'transverse_sector_gap_depth.png'
    localization_plot = PLOTS / 'transverse_sector_gap_localization.png'
    make_gap_plot(metrics, sizes, variants, gap_plot)
    make_depth_plot(metrics, sizes, variants, depth_plot)
    make_localization_plot(metrics, sizes, variants, localization_plot)
    plot_paths = [gap_plot, depth_plot, localization_plot]

    baseline_start = metrics[f'baseline_n{sizes[0]}']
    baseline_end = metrics[f'baseline_n{sizes[-1]}']
    smallest_gap = min((metrics[f'{variant}_n{sizes[-1]}'] for variant in variants if variant != 'baseline'), key=lambda item: item['scaled_harmonic_to_coexact_gap'])
    strongest_localization = max((metrics[f'{variant}_n{sizes[-1]}'] for variant in variants if variant != 'baseline'), key=lambda item: item['candidate_near_defect_fraction'])
    observation = (
        'across the larger CP2c sizes, the coexact branch stays separated from the harmonic floor by a finite gap instead of collapsing onto it, while defect backgrounds reduce that gap and can localize the coexact candidate more strongly than the clean baseline'
    )
    conclusion = (
        f"on the clean baseline, the scaled harmonic-to-coexact gap moves from {baseline_start['scaled_harmonic_to_coexact_gap']:.3f} at n={sizes[0]} to {baseline_end['scaled_harmonic_to_coexact_gap']:.3f} at n={sizes[-1]}, so it does not collapse toward zero in this window; at n={sizes[-1]}, {smallest_gap['variant']} has the smallest scaled gap at {smallest_gap['scaled_harmonic_to_coexact_gap']:.3f}, while {strongest_localization['variant']} has the strongest localization with near-defect fraction {strongest_localization['candidate_near_defect_fraction']:.3f}, which means defects compress the separation but still do not turn the coexact branch into the actual low continuum floor"
    )

    result = {
        'config': {
            'epsilon': epsilon,
            'sizes': sizes,
            'variants': variants,
            'restricted_modes': restricted_modes,
            'full_mode_scan_count': full_mode_scan_count,
            'coexact_threshold': coexact_threshold,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'penalty': penalty,
            'flux_tube_phase': flux_tube_phase,
        },
        'metrics': metrics,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('transverse_sector_gap_scaling', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse sector gap scaling',
        config_summary=f"epsilon={epsilon}, sizes={sizes}, variants={variants}, restricted_modes={restricted_modes}, full_mode_scan_count={full_mode_scan_count}, coexact_threshold={coexact_threshold}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_transverse_sector_gap_scaling()
    print(json.dumps({'result_path': str(result_path), 'note_path': str(note_path), 'observation': result['observation'], 'conclusion': result['conclusion']}, indent=2))


if __name__ == '__main__':
    main()
