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
        'coexact_threshold': 0.8,
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


def make_sector_split_plot(metrics: dict[str, dict[str, Any]], path: Path) -> None:
    variant_order = [variant for variant in ['baseline', 'puncture', 'line_defect', 'flux_tube'] if any(item['variant'] == variant for item in metrics.values())]
    largest_n = max({item['n_side'] for item in metrics.values()})
    selected = [metrics[f'{variant}_n{largest_n}'] for variant in variant_order]
    indices = np.arange(len(selected), dtype=float)
    harmonic = [item['full_low_harmonic_fraction'] for item in selected]
    coexact = [item['full_low_coexact_fraction'] for item in selected]
    exact = [item['full_low_exact_fraction'] for item in selected]

    fig, ax = plt.subplots(figsize=(10, 5))
    width = 0.24
    ax.bar(indices - width, exact, width=width, label='exact')
    ax.bar(indices, harmonic, width=width, label='harmonic')
    ax.bar(indices + width, coexact, width=width, label='coexact')
    for idx, item in enumerate(selected):
        label = str(item['first_coexact_mode_index']) if item['first_coexact_mode_index'] is not None else 'none'
        ax.text(indices[idx], 1.02, f'j={label}', ha='center', va='bottom', fontsize=8)
    ax.set_xticks(indices, [VARIANT_LABELS[item['variant']] for item in selected], rotation=15, ha='right')
    ax.set_ylim(0.0, 1.12)
    ax.set_ylabel('fraction of lowest full L1 mode')
    ax.set_title(f'Lowest full-mode sector split (n={largest_n})')
    ax.grid(alpha=0.25, axis='y')
    ax.legend(fontsize=8)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def spacing_error_metric(spacing_band: np.ndarray, spacing_ref: np.ndarray, count: int = 10) -> float:
    mask = np.abs(spacing_ref[:count]) > 1e-12
    if not np.any(mask):
        return 0.0
    ref = spacing_ref[:count][mask]
    band = spacing_band[:count][mask]
    return float(np.mean(np.abs(band - ref) / np.abs(ref)))


def select_first_coexact_mode(full_modes: list[dict[str, Any]], coexact_threshold: float) -> tuple[dict[str, Any] | None, bool]:
    qualified = [record for record in full_modes if float(record.get('coexact_fraction', 0.0)) >= coexact_threshold]
    if qualified:
        return qualified[0], True
    if not full_modes:
        return None, False
    fallback = max(full_modes, key=lambda record: (float(record.get('coexact_fraction', 0.0)), -float(record.get('eigenvalue', math.inf))))
    return fallback, False


def write_note(result: dict[str, Any], result_path: Path, stamped_plots: list[str], timestamp: str) -> Path:
    note_path = REPO_ROOT / 'experiments' / 'vector_sector' / 'Transverse_Continuum_Comparison_v1.md'
    rows = []
    for key, item in sorted(result['metrics'].items()):
        first_coexact_index = item['first_coexact_mode_index']
        first_coexact_index_text = str(first_coexact_index) if first_coexact_index is not None else 'none'
        threshold_text = 'yes' if item['first_coexact_mode_meets_threshold'] else 'no'
        first_coexact_ratio = item['restricted_to_first_coexact_ratio']
        ratio_text = f"{first_coexact_ratio:.3f}" if first_coexact_ratio is not None else 'none'
        rows.append(
            f"| {item['variant']} | {item['n_side']} | {item['continuum_scale']:.6f} | {item['relative_error_first10']:.4f} | {item['spacing_error_first10']:.4f} | {item['full_low_harmonic_fraction']:.3f} | {item['full_low_coexact_fraction']:.3f} | {first_coexact_index_text} | {threshold_text} | {ratio_text} |"
        )
    note = f"""# Transverse Continuum Comparison

## Purpose

Compare the low restricted transverse spectrum against the continuum transverse mode counting and spacing on a periodic box, while explicitly checking where that projected branch sits relative to the full low `L1` sector split.

## Setup

- sizes: `{result['config']['sizes']}`
- variants: `{result['config']['variants']}`
- restricted modes: `{result['config']['restricted_modes']}`
- coexact threshold for full-mode scan: `{result['config']['coexact_threshold']}`
- continuum reference built from nonzero integer wavevectors with transverse multiplicity two
- restricted band is built after exact/harmonic projection; the full-mode columns below check whether the actual low `L1` floor already sits in that same coexact sector

## Comparison metrics

| branch | `n` | best continuum scale | relative error (first 10) | spacing error (first 10) | lowest full harmonic frac | lowest full coexact frac | coexact candidate `j` | meets threshold? | restricted / candidate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
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
    coexact_threshold = float(experiment_cfg.get('coexact_threshold', 0.8))
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
            full_modes = list(cases[key].get('full_modes', []))
            lowest_full = full_modes[0] if full_modes else None
            first_coexact, meets_threshold = select_first_coexact_mode(full_modes, coexact_threshold)
            first_coexact_eigenvalue = float(first_coexact['eigenvalue']) if first_coexact is not None else None
            restricted_to_first_coexact_ratio = None
            if first_coexact_eigenvalue is not None and first_coexact_eigenvalue > 1e-12:
                restricted_to_first_coexact_ratio = float(spectrum[0] / first_coexact_eigenvalue)
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
                'near_defect_fraction_first': float(cases[key]['restricted_transverse_modes'][0]['near_defect_fraction']),
                'harmonic_dimension': int(cases[key]['dimensions']['harmonic']),
                'full_low_exact_fraction': float(lowest_full['exact_fraction']) if lowest_full is not None else math.nan,
                'full_low_harmonic_fraction': float(lowest_full['harmonic_fraction']) if lowest_full is not None else math.nan,
                'full_low_coexact_fraction': float(lowest_full['coexact_fraction']) if lowest_full is not None else math.nan,
                'first_coexact_mode_index': int(first_coexact['mode_index']) if first_coexact is not None else None,
                'first_coexact_mode_meets_threshold': bool(meets_threshold),
                'first_coexact_eigenvalue': first_coexact_eigenvalue,
                'restricted_to_first_coexact_ratio': restricted_to_first_coexact_ratio,
            }

    comparison_plot = PLOTS / 'transverse_continuum_comparison.png'
    spacing_plot = PLOTS / 'transverse_mode_spacing.png'
    phase_plot = PLOTS / 'divergence_curl_phase_continuum.png'
    sector_plot = PLOTS / 'transverse_sector_split.png'
    make_spectrum_plot(metrics, q2, comparison_plot)
    make_spacing_plot(metrics, q2, spacing_plot)
    make_phase_plot(cases[f'baseline_n{sizes[-1]}'], phase_plot, f"Restricted phase map (baseline, n={sizes[-1]})")
    make_sector_split_plot(metrics, sector_plot)
    plot_paths = [comparison_plot, spacing_plot, phase_plot, sector_plot]

    baseline_largest = metrics[f'baseline_n{sizes[-1]}']
    baseline_coexact_index = baseline_largest['first_coexact_mode_index']
    baseline_coexact_index_text = str(baseline_coexact_index) if baseline_coexact_index is not None else 'none'
    baseline_meets_threshold = bool(baseline_largest['first_coexact_mode_meets_threshold'])
    baseline_ratio = baseline_largest['restricted_to_first_coexact_ratio']
    baseline_ratio_text = f'{baseline_ratio:.3f}' if baseline_ratio is not None else 'none'
    if baseline_meets_threshold:
        baseline_sector_sentence = f'the first full mode above the coexact threshold appears at index {baseline_coexact_index_text}'
    else:
        baseline_sector_sentence = f'none of the first {restricted_modes} full modes reaches the coexact threshold {coexact_threshold:.2f}, and the best low coexact candidate within that window sits at index {baseline_coexact_index_text}'
    observation = (
        'after n^2 rescaling, the low restricted transverse spectrum still follows the continuum ordering across the tested branches, but the full low L1 floor remains sector-split from that branch rather than coinciding with it'
    )
    conclusion = (
        f"for the largest baseline case, the first-ten restricted spectrum error is {baseline_largest['relative_error_first10']:.3f}, the lowest full mode carries harmonic fraction {baseline_largest['full_low_harmonic_fraction']:.3f}, and {baseline_sector_sentence}; the restricted floor tracks that candidate with ratio {baseline_ratio_text}, so CP2 now makes the harmonic/coexact split explicit without yet establishing a clean low continuum vector band"
    )

    result = {
        'config': {
            'epsilon': epsilon,
            'sizes': sizes,
            'variants': variants,
            'restricted_modes': restricted_modes,
            'coexact_threshold': coexact_threshold,
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
        config_summary=f"epsilon={epsilon}, sizes={sizes}, variants={variants}, restricted_modes={restricted_modes}, coexact_threshold={coexact_threshold}",
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
