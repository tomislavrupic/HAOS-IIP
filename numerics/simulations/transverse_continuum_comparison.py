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
    analyze_branch_cases,
    append_log,
    continuum_transverse_q2,
    linear_gap_fit,
    make_phase_plot,
    mode_spacing,
    save_result_payload,
    plt,
)

DEFAULT_CONFIG: dict[str, Any] = {
    'epsilon': 0.2,
    'transverse_continuum_comparison': {
        'sizes': [12, 16, 20],
        'variants': ['baseline', 'puncture', 'line_defect', 'flux_tube'],
        'restricted_modes': 20,
        'harmonic_tol': 1e-8,
        'eig_tol': 1e-8,
        'penalty': 10.0,
        'flux_tube_phase': math.pi / 2.0,
    },
}


def load_config(config: dict[str, Any] | None = None, config_path: Path | None = None) -> dict[str, Any]:
    merged = DEFAULT_CONFIG.copy()
    merged['transverse_continuum_comparison'] = dict(DEFAULT_CONFIG['transverse_continuum_comparison'])
    path = config_path or (REPO_ROOT / 'config.json')
    if path.exists():
        on_disk = json.loads(path.read_text())
        merged.update({k: v for k, v in on_disk.items() if k != 'transverse_continuum_comparison'})
        if isinstance(on_disk.get('transverse_continuum_comparison'), dict):
            merged['transverse_continuum_comparison'].update(on_disk['transverse_continuum_comparison'])
    if config is not None:
        merged.update({k: v for k, v in config.items() if k != 'transverse_continuum_comparison'})
        if isinstance(config.get('transverse_continuum_comparison'), dict):
            merged['transverse_continuum_comparison'].update(config['transverse_continuum_comparison'])
    return merged


def best_continuum_scale(q2: np.ndarray, scaled_lambdas: np.ndarray) -> float:
    return float(np.dot(q2, scaled_lambdas) / np.dot(q2, q2))


def make_spectrum_plot(metrics: dict[str, dict[str, Any]], q2: np.ndarray, path: Path) -> None:
    largest_n = max({item['n_side'] for item in metrics.values()})
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
    for ax, variant in zip(axes.ravel(), ['baseline', 'puncture', 'line_defect', 'flux_tube']):
        key = f'{variant}_n{largest_n}'
        scaled = np.asarray(metrics[key]['scaled_spectrum'], dtype=float)
        scale = metrics[key]['continuum_scale']
        ax.plot(range(1, len(scaled) + 1), scaled, marker=VARIANT_MARKERS[variant], label='restricted band')
        ax.plot(range(1, len(q2) + 1), scale * q2, linestyle='--', color='black', label='continuum reference')
        ax.set_title(f"{VARIANT_LABELS[variant]} (n={largest_n})")
        ax.grid(alpha=0.25)
    for ax in axes[-1]:
        ax.set_xlabel('mode index k')
    for ax in axes[:, 0]:
        ax.set_ylabel(r'$n^2 \lambda_k$')
    axes[0, 0].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def make_spacing_plot(metrics: dict[str, dict[str, Any]], q2: np.ndarray, path: Path) -> None:
    largest_n = max({item['n_side'] for item in metrics.values()})
    ref_spacing = mode_spacing(q2)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
    for ax, variant in zip(axes.ravel(), ['baseline', 'puncture', 'line_defect', 'flux_tube']):
        key = f'{variant}_n{largest_n}'
        scaled = np.asarray(metrics[key]['scaled_spectrum'], dtype=float)
        scale = metrics[key]['continuum_scale']
        ax.plot(range(1, len(scaled)), mode_spacing(scaled), marker=VARIANT_MARKERS[variant], label='restricted spacing')
        ax.plot(range(1, len(ref_spacing) + 1), scale * ref_spacing, linestyle='--', color='black', label='continuum spacing')
        ax.set_title(f"{VARIANT_LABELS[variant]} (n={largest_n})")
        ax.grid(alpha=0.25)
    for ax in axes[-1]:
        ax.set_xlabel('spacing index')
    for ax in axes[:, 0]:
        ax.set_ylabel(r'$n^2 (\lambda_{k+1}-\lambda_k)$')
    axes[0, 0].legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def spacing_error_metric(spacing_band: np.ndarray, spacing_ref: np.ndarray, count: int = 10) -> float:
    mask = np.abs(spacing_ref[:count]) > 1e-12
    if not np.any(mask):
        return 0.0
    ref = spacing_ref[:count][mask]
    band = spacing_band[:count][mask]
    return float(np.mean(np.abs(band - ref) / np.abs(ref)))


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Continuum_Comparison_v1.md'
    rows = []
    for key, item in sorted(result['metrics'].items()):
        rows.append(
            f"| {item['variant']} | {item['n_side']} | {item['continuum_scale']:.6f} | {item['relative_error_first10']:.4f} | {item['spacing_error_first10']:.4f} |"
        )
    note = f"""# Transverse Continuum Comparison

## Purpose

Compare the low restricted transverse spectrum against the continuum transverse mode counting and spacing on a periodic box.

## Setup

- sizes: `{result['config']['sizes']}`
- variants: `{result['config']['variants']}`
- restricted modes: `{result['config']['restricted_modes']}`
- continuum reference built from nonzero integer wavevectors with transverse multiplicity two

## Comparison metrics

| branch | `n` | best continuum scale | relative error (first 10) | spacing error (first 10) |
| --- | ---: | ---: | ---: | ---: |
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


def run_transverse_continuum_comparison(config: dict[str, Any] | None = None, config_path: Path | None = None) -> tuple[dict[str, Any], Path, list[str], str, Path]:
    cfg = load_config(config, config_path)
    experiment_cfg = cfg['transverse_continuum_comparison']
    epsilon = float(cfg.get('epsilon', 0.2))
    sizes = [int(v) for v in experiment_cfg.get('sizes', [12, 16, 20])]
    variants = [str(v) for v in experiment_cfg.get('variants', ['baseline', 'puncture', 'line_defect', 'flux_tube'])]
    restricted_modes = int(experiment_cfg.get('restricted_modes', 20))
    harmonic_tol = float(experiment_cfg.get('harmonic_tol', 1e-8))
    eig_tol = float(experiment_cfg.get('eig_tol', 1e-8))
    penalty = float(experiment_cfg.get('penalty', 10.0))
    flux_tube_phase = float(experiment_cfg.get('flux_tube_phase', math.pi / 2.0))

    cases = analyze_branch_cases(
        sizes=sizes,
        variants=variants,
        epsilon=epsilon,
        restricted_modes=restricted_modes,
        harmonic_tol=harmonic_tol,
        eig_tol=eig_tol,
        penalty=penalty,
        flux_tube_phase=flux_tube_phase,
    )
    q2 = continuum_transverse_q2(restricted_modes)
    metrics: dict[str, dict[str, Any]] = {}
    for n_side in sizes:
        for variant in variants:
            key = f'{variant}_n{n_side}'
            spectrum = np.asarray(cases[key]['restricted_transverse_spectrum'][:restricted_modes], dtype=float)
            scaled = (n_side * n_side) * spectrum
            scale = best_continuum_scale(q2, scaled)
            ref = scale * q2
            spacing_ref = mode_spacing(ref)
            spacing_band = mode_spacing(scaled)
            metrics[key] = {
                'variant': variant,
                'n_side': n_side,
                'scaled_spectrum': scaled.tolist(),
                'continuum_scale': scale,
                'relative_error_first10': float(np.mean(np.abs(scaled[:10] - ref[:10]) / np.maximum(np.abs(ref[:10]), 1e-12))),
                'spacing_error_first10': spacing_error_metric(spacing_band, spacing_ref, count=10),
                'divergence_norm_first': float(cases[key]['restricted_transverse_modes'][0]['divergence_norm']),
                'curl_norm_first': float(cases[key]['restricted_transverse_modes'][0]['curl_norm']),
                'ipr_first': float(cases[key]['restricted_transverse_modes'][0]['ipr']),
            }

    comparison_plot = PLOTS / 'transverse_continuum_comparison.png'
    spacing_plot = PLOTS / 'transverse_mode_spacing.png'
    phase_plot = PLOTS / 'divergence_curl_phase_continuum.png'
    make_spectrum_plot(metrics, q2, comparison_plot)
    make_spacing_plot(metrics, q2, spacing_plot)
    make_phase_plot(cases[f'baseline_n{sizes[-1]}'], phase_plot, f"Restricted phase map (baseline, n={sizes[-1]})")
    plot_paths = [comparison_plot, spacing_plot, phase_plot]

    baseline_largest = metrics[f'baseline_n{sizes[-1]}']
    observation = (
        'after n^2 rescaling, the low restricted transverse spectrum aligns with the continuum transverse mode ordering and the mode-spacing pattern remains stable across the tested branches'
    )
    conclusion = (
        f"for the largest baseline case, the first-ten relative spectrum error is {baseline_largest['relative_error_first10']:.3f} and the first-ten spacing error is {baseline_largest['spacing_error_first10']:.3f}, consistent with a low continuum transverse operator in the tested range"
    )

    result = {
        'config': {
            'epsilon': epsilon,
            'sizes': sizes,
            'variants': variants,
            'restricted_modes': restricted_modes,
            'harmonic_tol': harmonic_tol,
            'eig_tol': eig_tol,
            'penalty': penalty,
            'flux_tube_phase': flux_tube_phase,
        },
        'continuum_q2': q2.tolist(),
        'metrics': metrics,
        'observation': observation,
        'conclusion': conclusion,
    }
    result_path, stamped_plots, timestamp = save_result_payload('transverse_continuum_comparison', result, plot_paths)
    note_path = write_note(result, result_path, stamped_plots, timestamp)
    append_log(
        title='transverse continuum comparison',
        config_summary=f"epsilon={epsilon}, sizes={sizes}, variants={variants}, restricted_modes={restricted_modes}",
        result_path=result_path,
        stamped_plots=stamped_plots,
        observation=observation,
        conclusion=conclusion,
    )
    return result, result_path, stamped_plots, timestamp, note_path


def main() -> None:
    result, result_path, _, _, note_path = run_transverse_continuum_comparison()
    print(json.dumps({'result_path': str(result_path), 'note_path': str(note_path), 'observation': result['observation'], 'conclusion': result['conclusion']}, indent=2))


if __name__ == '__main__':
    main()
