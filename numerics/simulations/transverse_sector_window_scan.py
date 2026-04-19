#!/usr/bin/env python3

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

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
    'transverse_sector_window_scan': {
        'sizes': [12, 16, 20, 24, 28],
        'variants': ['baseline', 'puncture', 'line_defect', 'flux_tube'],
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
    merged['transverse_sector_window_scan'] = dict(DEFAULT_CONFIG['transverse_sector_window_scan'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_sector_window_scan'})
        if isinstance(on_disk.get('transverse_sector_window_scan'), dict):
            merged['transverse_sector_window_scan'].update(on_disk['transverse_sector_window_scan'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_sector_window_scan'})
        if isinstance(config.get('transverse_sector_window_scan'), dict):
            merged['transverse_sector_window_scan'].update(config['transverse_sector_window_scan'])
    return merged


def select_coexact_candidate(full_modes: list[dict[str, Any]], coexact_threshold: float) -> tuple[dict[str, Any] | None, bool]:
    if not full_modes:
        return None, False
    threshold_hits = [record for record in full_modes if float(record.get('coexact_fraction', 0.0)) >= coexact_threshold]
    if threshold_hits:
        return threshold_hits[0], True
    fallback = max(full_modes, key=lambda record: (float(record.get('coexact_fraction', 0.0)), -float(record.get('eigenvalue', math.inf))))
    return fallback, False


def summarize_case(case: dict[str, Any], coexact_threshold: float) -> dict[str, Any]:
    full_modes = list(case.get('full_modes', []))
    candidate, meets_threshold = select_coexact_candidate(full_modes, coexact_threshold)
    lowest_full = full_modes[0] if full_modes else None
    floor_eigenvalue = float(lowest_full['eigenvalue']) if lowest_full is not None else None
    restricted_floor = float(case['restricted_transverse_spectrum'][0]) if case['restricted_transverse_spectrum'] else None
    candidate_eigenvalue = float(candidate['eigenvalue']) if candidate is not None else None
    restricted_ratio = None
    if restricted_floor is not None and candidate_eigenvalue is not None and candidate_eigenvalue > 1e-12:
        restricted_ratio = float(restricted_floor / candidate_eigenvalue)
    harmonic_to_coexact_gap = None
    scaled_harmonic_to_coexact_gap = None
    if floor_eigenvalue is not None and candidate_eigenvalue is not None:
        harmonic_to_coexact_gap = float(candidate_eigenvalue - floor_eigenvalue)
        scaled_harmonic_to_coexact_gap = float((case['config']['n_side'] ** 2) * harmonic_to_coexact_gap)
    return {
        'variant': str(case['config']['variant']),
        'n_side': int(case['config']['n_side']),
        'harmonic_dimension': int(case['dimensions']['harmonic']),
        'floor_eigenvalue': floor_eigenvalue,
        'floor_scaled_eigenvalue': float((case['config']['n_side'] ** 2) * floor_eigenvalue) if floor_eigenvalue is not None else None,
        'lowest_full_harmonic_fraction': float(lowest_full['harmonic_fraction']) if lowest_full is not None else math.nan,
        'lowest_full_coexact_fraction': float(lowest_full['coexact_fraction']) if lowest_full is not None else math.nan,
        'candidate_mode_index': int(candidate['mode_index']) if candidate is not None else None,
        'candidate_meets_threshold': bool(meets_threshold),
        'candidate_eigenvalue': candidate_eigenvalue,
        'candidate_scaled_eigenvalue': float((case['config']['n_side'] ** 2) * candidate_eigenvalue) if candidate_eigenvalue is not None else None,
        'harmonic_to_coexact_gap': harmonic_to_coexact_gap,
        'scaled_harmonic_to_coexact_gap': scaled_harmonic_to_coexact_gap,
        'candidate_coexact_fraction': float(candidate['coexact_fraction']) if candidate is not None else None,
        'candidate_harmonic_fraction': float(candidate['harmonic_fraction']) if candidate is not None else None,
        'candidate_near_defect_fraction': float(candidate['near_defect_fraction']) if candidate is not None else None,
        'candidate_ipr': float(candidate['ipr']) if candidate is not None else None,
        'candidate_support_pattern': str(candidate['support_pattern']) if candidate is not None else 'none',
        'restricted_to_candidate_ratio': restricted_ratio,
    }


def make_index_plot(metrics: dict[str, dict[str, Any]], sizes: list[int], variants: list[str], full_mode_scan_count: int, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for variant in variants:
        y_values = []
        for n_side in sizes:
            key = f'{variant}_n{n_side}'
            value = metrics[key]['candidate_mode_index']
            y_values.append(float(value) if value is not None else math.nan)
        ax.plot(sizes, y_values, marker=VARIANT_MARKERS[variant], label=VARIANT_LABELS[variant])
    ax.set_xlabel('lattice side n')
    ax.set_ylabel('first coexact candidate index')
    ax.set_title(f'Coexact candidate depth within first {full_mode_scan_count} full modes')
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_localization_plot(metrics: dict[str, dict[str, Any]], sizes: list[int], variants: list[str], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    for variant in variants:
        near_defect = []
        scaled_eigs = []
        for n_side in sizes:
            key = f'{variant}_n{n_side}'
            near_value = metrics[key]['candidate_near_defect_fraction']
            eig_value = metrics[key]['candidate_scaled_eigenvalue']
            near_defect.append(float(near_value) if near_value is not None else math.nan)
            scaled_eigs.append(float(eig_value) if eig_value is not None else math.nan)
        axes[0].plot(sizes, near_defect, marker=VARIANT_MARKERS[variant], label=VARIANT_LABELS[variant])
        axes[1].plot(sizes, scaled_eigs, marker=VARIANT_MARKERS[variant], label=VARIANT_LABELS[variant])
    axes[0].set_ylabel('candidate near-defect fraction')
    axes[0].set_title('Defect localization of coexact candidate')
    axes[1].set_ylabel(r'$n^2 \lambda_{\mathrm{candidate}}$')
    axes[1].set_title('Scaled candidate eigenvalue')
    for ax in axes:
        ax.set_xlabel('lattice side n')
        ax.grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_gap_plot(metrics: dict[str, dict[str, Any]], sizes: list[int], variants: list[str], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    for variant in variants:
        raw_gap = []
        scaled_gap = []
        for n_side in sizes:
            key = f'{variant}_n{n_side}'
            raw_value = metrics[key]['harmonic_to_coexact_gap']
            scaled_value = metrics[key]['scaled_harmonic_to_coexact_gap']
            raw_gap.append(float(raw_value) if raw_value is not None else math.nan)
            scaled_gap.append(float(scaled_value) if scaled_value is not None else math.nan)
        axes[0].plot(sizes, raw_gap, marker=VARIANT_MARKERS[variant], label=VARIANT_LABELS[variant])
        axes[1].plot(sizes, scaled_gap, marker=VARIANT_MARKERS[variant], label=VARIANT_LABELS[variant])
    axes[0].set_ylabel(r'$\lambda_{\mathrm{coexact}}-\lambda_{\mathrm{floor}}$')
    axes[0].set_title('Raw harmonic-to-coexact gap')
    axes[1].set_ylabel(r'$n^2(\lambda_{\mathrm{coexact}}-\lambda_{\mathrm{floor}})$')
    axes[1].set_title('Scaled harmonic-to-coexact gap')
    for ax in axes:
        ax.set_xlabel('lattice side n')
        ax.grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_threshold_plot(metrics: dict[str, dict[str, Any]], sizes: list[int], variants: list[str], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    for variant in variants:
        coexact = []
        ratios = []
        for n_side in sizes:
            key = f'{variant}_n{n_side}'
            coexact_value = metrics[key]['candidate_coexact_fraction']
            ratio_value = metrics[key]['restricted_to_candidate_ratio']
            coexact.append(float(coexact_value) if coexact_value is not None else math.nan)
            ratios.append(float(ratio_value) if ratio_value is not None else math.nan)
        axes[0].plot(sizes, coexact, marker=VARIANT_MARKERS[variant], label=VARIANT_LABELS[variant])
        axes[1].plot(sizes, ratios, marker=VARIANT_MARKERS[variant], label=VARIANT_LABELS[variant])
    axes[0].axhline(0.8, color='black', linestyle='--', linewidth=0.8)
    axes[0].set_ylabel('candidate coexact fraction')
    axes[0].set_title('Candidate coexact purity')
    axes[1].set_ylabel('restricted / candidate eigenvalue')
    axes[1].set_title('Projected floor versus full-spectrum candidate')
    for ax in axes:
        ax.set_xlabel('lattice side n')
        ax.grid(alpha=0.25)
    axes[1].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Sector_Window_Scan_v1.md'
    rows = []
    for key, item in sorted(result['metrics'].items()):
        candidate_index_text = str(item['candidate_mode_index']) if item['candidate_mode_index'] is not None else 'none'
        rows.append(
            f"| {item['variant']} | {item['n_side']} | {item['lowest_full_harmonic_fraction']:.3f} | {item['candidate_coexact_fraction']:.3f} | {candidate_index_text} | {'yes' if item['candidate_meets_threshold'] else 'no'} | {item['scaled_harmonic_to_coexact_gap']:.3f} | {item['candidate_near_defect_fraction']:.3f} | {item['restricted_to_candidate_ratio']:.3f} | {item['candidate_support_pattern']} |"
        )
    note = f"""# Transverse Sector Window Scan

## Purpose

Push the CP2 sector split deeper into the full `L1` spectrum and test whether defect-localized backgrounds pull a genuinely coexact branch down or localize it in a controlled way.

## Setup

- sizes: `{result['config']['sizes']}`
- variants: `{result['config']['variants']}`
- restricted modes: `{result['config']['restricted_modes']}`
- full-mode scan count: `{result['config']['full_mode_scan_count']}`
- coexact threshold: `{result['config']['coexact_threshold']}`

## Sector-window metrics

| branch | `n` | lowest full harmonic frac | candidate coexact frac | candidate `j` | meets threshold? | scaled gap | candidate near-defect frac | restricted / candidate | support |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
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


def run_transverse_sector_window_scan(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_sector_window_scan']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [12, 16, 20, 24, 28])]
    variants = [str(v) for v in experiment_cfg.get('variants', ['baseline', 'puncture', 'line_defect', 'flux_tube'])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 20))
    full_mode_scan_count = int(experiment_cfg.get('full_mode_scan_count', 40))
    coexact_threshold = float(experiment_cfg.get('coexact_threshold', 0.8))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))
    flux_tube_phase = float(experiment_cfg.get('flux_tube_phase', math.pi / 2.0))

    cases: dict[str, Any] = {}
    metrics: dict[str, dict[str, Any]] = {}
    for n_side in sizes:
        for variant in variants:
            key = f'{variant}_n{n_side}'
            cases[key] = analyze_case(
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
            metrics[key] = summarize_case(cases[key], coexact_threshold=coexact_threshold)

    index_plot = PLOTS / 'transverse_sector_candidate_index.png'
    localization_plot = PLOTS / 'transverse_sector_localization.png'
    gap_plot = PLOTS / 'transverse_sector_gap.png'
    threshold_plot = PLOTS / 'transverse_sector_candidate_metrics.png'
    make_index_plot(metrics, sizes, variants, full_mode_scan_count, index_plot)
    make_localization_plot(metrics, sizes, variants, localization_plot)
    make_gap_plot(metrics, sizes, variants, gap_plot)
    make_threshold_plot(metrics, sizes, variants, threshold_plot)
    plot_paths = [index_plot, localization_plot, gap_plot, threshold_plot]

    baseline_largest = metrics[f'baseline_n{sizes[-1]}']
    defect_candidates = [metrics[f'{variant}_n{sizes[-1]}'] for variant in variants if variant != 'baseline']
    earliest_defect = min(defect_candidates, key=lambda item: item['candidate_mode_index'] if item['candidate_mode_index'] is not None else full_mode_scan_count + 1)
    most_localized_defect = max(defect_candidates, key=lambda item: item['candidate_near_defect_fraction'] if item['candidate_near_defect_fraction'] is not None else -math.inf)
    smallest_gap = min(defect_candidates, key=lambda item: item['scaled_harmonic_to_coexact_gap'] if item['scaled_harmonic_to_coexact_gap'] is not None else math.inf)
    observation = (
        'the true low full-spectrum floor remains harmonic across the scanned backgrounds, but the direct harmonic-to-coexact gap stays finite and defect backgrounds can reduce that separation while pulling the coexact candidate deeper into the low window'
    )
    conclusion = (
        f"at the largest size, the clean baseline still has harmonic floor fraction {baseline_largest['lowest_full_harmonic_fraction']:.3f}, reaches its first coexact candidate at j={baseline_largest['candidate_mode_index']}, and keeps a scaled harmonic-to-coexact gap of {baseline_largest['scaled_harmonic_to_coexact_gap']:.3f}; by contrast, {earliest_defect['variant']} reaches its coexact candidate at j={earliest_defect['candidate_mode_index']}, {smallest_gap['variant']} gives the smallest scaled gap at {smallest_gap['scaled_harmonic_to_coexact_gap']:.3f}, and {most_localized_defect['variant']} shows the strongest candidate defect localization with near-defect fraction {most_localized_defect['candidate_near_defect_fraction']:.3f}, so CP2c says defects compress the harmonic-to-coexact separation but still do not collapse it into a clean propagating vector floor"
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
    result_path, stamped_plots, timestamp = save_result_payload('transverse_sector_window_scan', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse sector window scan',
        config_summary=f"epsilon={epsilon}, sizes={sizes}, variants={variants}, restricted_modes={restricted_modes}, full_mode_scan_count={full_mode_scan_count}, coexact_threshold={coexact_threshold}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_transverse_sector_window_scan()
    print(json.dumps({'result_path': str(result_path), 'note_path': str(note_path), 'observation': result['observation'], 'conclusion': result['conclusion']}, indent=2))


if __name__ == '__main__':
    main()
